//
// Created by Serban Porumbescu on 5/4/23.
//

#include <set>
#include <type_traits>
#include <igl/Timer.h>
#include <wmtk/ExecutionScheduler.hpp>
#include <wmtk/utils/TupleUtils.hpp>

namespace app::def {

template <typename Scalar>
DelaunayEdgeFlip<Scalar>::DelaunayEdgeFlip(
    std::vector< Eigen::Vector3<Scalar> >_m_vertex_positions,
    const std::vector<std::array<size_t, 3>>& tris,
    int num_threads)
{
    NUM_THREADS = num_threads;
    p_vertex_attrs = &vertex_attrs;

    // Create the mesh connectivity
    wmtk::TriMesh::create_mesh(_m_vertex_positions.size(), tris);
    assert(vertex_attrs.size() == _m_vertex_positions.size());

    // Partition the mesh
    auto partition_ids = wmtk::partition_morton(_m_vertex_positions, NUM_THREADS);
    assert(partition_ids.size() == _m_vertex_positions.size());
    wmtk::logger().info("partition_ids: {}", partition_ids.size());

    std::set<unsigned long> partition_set;
    for (auto i = 0; i < _m_vertex_positions.size(); i++) {
        partition_set.insert(partition_ids[i]);
        vertex_attrs[i] = {_m_vertex_positions[i], partition_ids[i]};
    }

    for (auto pid : partition_set) {
        wmtk::logger().info("pid: {}\n", pid);
    }
}

template <typename Scalar>
bool DelaunayEdgeFlip<Scalar>::invariants(const std::vector<Tuple>& new_tris)
{
    return TriMesh::invariants(new_tris);
}

auto edge_locker = [](auto& m, const auto& e, int task_id) {
    // TODO: this should not be here
    return m.try_set_edge_mutex_two_ring(e, task_id);
};

template <typename Scalar>
void DelaunayEdgeFlip<Scalar>::swap_all_edges()
{
    wmtk::logger().info("***** DelaunayEdgeFlip Collecting edges to swap *****");

    igl::Timer timer;
    timer.start();

    // vector of pair("edge_swap", Tuple)
    auto collect_all_ops = std::vector<std::pair<std::string, Tuple>>();

    // concurrent vector representing our Tuples (edges to swap)
    auto collect_tuples = tbb::concurrent_vector<Tuple>();

    // grab all our edges in parallel
    for_each_edge([&](auto& tup) { collect_tuples.emplace_back(tup); });

    // prepare the operations (e.g., "edge_swap") to perform when we're ready
    collect_all_ops.reserve(collect_tuples.size());
    for (auto& t : collect_tuples) {
        collect_all_ops.emplace_back("edge_swap", t);
    }

    wmtk::logger().info("***** swap get edges time *****: {} ms", timer.getElapsedTimeInMilliSec());
    wmtk::logger().info("DelaunayEdgeFlip: size for edges to swap is {}", collect_all_ops.size());

    auto setup_and_execute = [&](auto executor) {
        executor.num_threads = NUM_THREADS;

        executor.stopping_criterion_checking_frequency = collect_tuples.size() / 2;

        // new operations should NOT be added to priority queue
        executor.should_renew = [](auto) { return false; };

        executor.priority = [](auto& m, auto op, const Tuple& e) {
            return 1.0; //m.compute_vertex_valence(e);
        };
        executor(*this, collect_all_ops);
    };

    if(NUM_THREADS > 1) {
        wmtk::logger().info("DelaunayEdgeFlip NUM_THREADS kPartition {}", NUM_THREADS);
        auto executor = wmtk::ExecutePass<DelaunayEdgeFlip, wmtk::ExecutionPolicy::kPartition>();
        executor.lock_vertices = edge_locker;
        setup_and_execute(executor);
    } else {
        wmtk::logger().info("DelaunayEdgeFlip NUM_THREADS kSeq {}", NUM_THREADS);
        auto executor = wmtk::ExecutePass<DelaunayEdgeFlip, wmtk::ExecutionPolicy::kSeq>();
        setup_and_execute(executor);
    }
}

template <typename Scalar>
bool DelaunayEdgeFlip<Scalar>::write_triangle_mesh(std::string path)
{
    // write the collapsed mesh into a obj and assert the mesh is manifold
    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(vert_capacity(), 3);
    for (auto& t : get_vertices()) {
        auto i = t.vid(*this);
        V.row(i) = vertex_attrs[i].pos;
    }

    Eigen::MatrixXi F = Eigen::MatrixXi::Constant(tri_capacity(), 3, -1);
    for (auto& t : get_faces()) {
        auto i = t.fid(*this);
        auto vs = oriented_tri_vertices(t);
        for (int j = 0; j < 3; j++) {
            F(i, j) = vs[j].vid(*this);
        }
    }
    igl::write_triangle_mesh(path, V, F);
    bool manifold = check_edge_manifold();
    assert(manifold);

    return manifold;
}

template <typename Scalar>
bool DelaunayEdgeFlip<Scalar>::swap_edge_before(const Tuple& t)
{
    if (!TriMesh::swap_edge_before(t))
    {
        return false;
    }

    if (vertex_attrs[t.vid(*this)].freeze &&
        vertex_attrs[t.switch_vertex(*this).vid(*this)].freeze)
    {
        return false;
    }

    // save the current edge's positions
    edge_position_cache.local().v0p = vertex_attrs[t.vid(*this)].pos;
    edge_position_cache.local().v1p = vertex_attrs[t.switch_vertex(*this).vid(*this)].pos;

    return true;
}


template <typename Scalar>
bool DelaunayEdgeFlip<Scalar>::swap_edge_after(const TriMesh::Tuple& t)
{
    // Follow from:
    // https://github.com/Ahdhn/GPUDelaunayEdgeFlip/blob/master/GDEF/exclusion_processing.cuh

    // first check if the edge formed by v0-v1 is a delaunay edge
    // where v2 and v3 are the opposite vertices to the edge
    //    0
    //  / | \
    // 3  |  2
    // \  |  /
    //    1
    // if not delaunay, then we check if flipping it won't create a foldover
    // The case below would create a fold over
    //      0
    //    / | \
    //   /  1  \
    //  / /  \  \
    //  2       3
    // http://www.dma.fi.upm.es/personal/mabellanas/tfcs/flips/Intercambios/html/teoria/teoria_del_ing.htm

    // Grab the positions saved before the edge swap
    auto v0 = edge_position_cache.local().v0p;
    auto v1 = edge_position_cache.local().v1p;

    // Grab the current, post edge swap, positions
    auto v2 = vertex_attrs[t.vid(*this)].pos;
    auto v3 = vertex_attrs[t.switch_vertex(*this).vid(*this)].pos;

    // Compute the Delaunay condition
    // Return false if it fails

    // find the angle between S, M, Q vertices (i.e., angle at M)
    auto angle_between_three_vertices =
        [](const auto& S, const auto& M, const auto& Q) {
//            using Scalar = typename decltype(S)::Scalar;
            using Vector3 = Eigen::Vector3<Scalar>;

            Vector3 p1 = S - M;
            Vector3 p2 = Q - M;
            Scalar dot_pro = p1.dot(p2);

            auto acos_helper = [](auto x) {
                if constexpr (std::is_same_v<std::decay_t<decltype(x)>, float>) {
                    return std::acosf(x);
                } else {
                    return std::acos(x);
                }
            };

            return acos_helper(dot_pro / (p1.norm() * p2.norm()));
        };

    auto lambda = angle_between_three_vertices(v0, v2, v1);
    auto gamma  = angle_between_three_vertices(v0, v3, v1);

    if (lambda + gamma <= static_cast<Scalar>(M_PI) + std::numeric_limits<Scalar>::epsilon())
    {
        // wmtk::logger().info("returning false 1");
        return false;
    } else
    {
        // check if flipping won't create foldover
        double alpha, beta;

        alpha = angle_between_three_vertices(v3, v0, v1);
        beta  = angle_between_three_vertices(v2, v0, v1);

        if (alpha + beta >= static_cast<Scalar>(M_PI) - std::numeric_limits<Scalar>::epsilon())
        {
            // wmtk::logger().info("returning false 2");
            return false;
        }

        alpha = angle_between_three_vertices(v3, v1, v0);
        beta  = angle_between_three_vertices(v2, v1, v0);

        if (alpha + beta >= static_cast<Scalar>(M_PI) - std::numeric_limits<Scalar>::epsilon())
        {
            // wmtk::logger().info("returning false 3");
            return false;
        }
    }

    // wmtk::logger().info("returning true 1");
    return true;
}

} // namespace app::def