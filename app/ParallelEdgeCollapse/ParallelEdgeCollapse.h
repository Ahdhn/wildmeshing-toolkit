#pragma once

#include <igl/write_triangle_mesh.h>
#include <tbb/concurrent_priority_queue.h>
#include <tbb/concurrent_vector.h>
#include <tbb/parallel_for.h>
#include <tbb/task_group.h>
#include <wmtk/ConcurrentTriMesh.h>
#include <wmtk/utils/VectorUtils.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <queue>


namespace Edge2d {

class ElementInQueue
{
public:
    wmtk::ConcurrentTriMesh::Tuple edge;
    double weight;
    int times_skipped = 0;

    ElementInQueue() {}
    ElementInQueue(const wmtk::ConcurrentTriMesh::Tuple& e, double w)
        : edge(e)
        , weight(w)
    {}
};

struct cmp_l
{
    bool operator()(const ElementInQueue& e1, const ElementInQueue& e2)
    {
        if (e1.weight == e2.weight) return e1.edge.vid() > e2.edge.vid();
        return e1.weight < e2.weight;
    }
};

struct cmp_s
{
    bool operator()(const ElementInQueue& e1, const ElementInQueue& e2)
    {
        if (e1.weight == e2.weight) return e1.edge.vid() < e2.edge.vid();
        return e1.weight > e2.weight;
    }
};

class ParallelEdgeCollapse : public wmtk::ConcurrentTriMesh
{
public:
    tbb::concurrent_vector<Eigen::Vector3d> m_vertex_positions;
    tbb::concurrent_vector<int> m_vertex_partition_id;
    tbb::spin_mutex rw_lock;
    int NUM_THREADS;

    ParallelEdgeCollapse(
        tbb::concurrent_vector<Eigen::Vector3d> _m_vertex_positions,
        tbb::concurrent_vector<int> _m_vertex_partition_id,
        int num_t = 1)
        : m_vertex_positions(_m_vertex_positions)
        , m_vertex_partition_id(_m_vertex_partition_id)
        , NUM_THREADS(num_t)
    {}

    ~ParallelEdgeCollapse() {}

    void print_num_attributes()
    {
        std::cout << m_vertex_positions.size() << std::endl;
        std::cout << m_vertex_partition_id.size() << std::endl;
    }

    // TODO cannot be used for parallel
    struct CollapseInfoCache
    {
        Eigen::Vector3d v1p;
        Eigen::Vector3d v2p;
    } collapse_cache;

    bool collapse_before(const Tuple& t) override;
    bool collapse_after(const Tuple& t) override;

    // write the collapsed mesh into a obj
    bool write_triangle_mesh(std::string path)
    {
        Eigen::MatrixXd V = Eigen::MatrixXd::Zero(m_vertex_positions.size(), 3);
        for (auto& t : get_vertices()) {
            auto i = t.vid();
            V.row(i) = m_vertex_positions[i];
        }

        Eigen::MatrixXi F = Eigen::MatrixXi::Constant(tri_capacity(), 3, -1);
        for (auto& t : get_faces()) {
            auto i = t.fid();
            auto vs = oriented_tri_vertices(t);
            for (int j = 0; j < 3; j++) {
                F(i, j) = vs[j].vid();
            }
        }

        return igl::write_triangle_mesh(path, V, F);
    }

    // old implementation
    bool collapse_shortest();
    void collapse_shortest_stuff(
        tbb::concurrent_priority_queue<ElementInQueue, cmp_s>& ec_queue,
        int& target_vertex_count,
        int task_id);

    bool collapse_shortest(int target_vertex_count);

    bool collapse_qec();
    // get the quadrix in form of an array of 10 floating point numbers
    std::array<double, 10> compute_Q_f(wmtk::ConcurrentTriMesh::Tuple& t);

    std::array<double, 10> compute_Q_v(wmtk::ConcurrentTriMesh::Tuple& t);

    double compute_cost_for_v(wmtk::TriMesh::Tuple& v_tuple);

    void update_position(size_t v1, size_t v2, Tuple& new_vert);

    void move_vertex_attribute(size_t from, size_t to) override
    {
        m_vertex_positions[to] = m_vertex_positions[from];
    }

    void resize_attributes(size_t v, size_t e, size_t t) override
    {
        m_vertex_positions.grow_to_at_least(v);
        m_vertex_partition_id.grow_to_at_least(v);
    }
};


} // namespace Edge2d