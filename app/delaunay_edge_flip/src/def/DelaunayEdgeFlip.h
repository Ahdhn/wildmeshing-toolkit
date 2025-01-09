#pragma once
#include <wmtk/utils/PartitionMesh.h>
#include <wmtk/utils/VectorUtils.h>
#include "wmtk/AttributeCollection.hpp"

// clang-format off
#include <wmtk/utils/DisableWarnings.hpp>
#include <igl/write_triangle_mesh.h>
#include <tbb/concurrent_priority_queue.h>
#include <tbb/concurrent_vector.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>
//#include <fastenvelope/FastEnvelope.h>
#include <wmtk/utils/EnableWarnings.hpp>
// clang-format on

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <atomic>
#include <memory>
#include <queue>
namespace app::def {

struct VertexAttributes
{
    Eigen::Vector3d pos;
    size_t partition_id;
    bool freeze = false;
};

class DelaunayEdgeFlip : public wmtk::TriMesh
{
public:

    // Flag to determine whether or not to continue
    // edge flips in subsequent iterations. We want to stop
    // when no more edges can be flipped.
    bool any_edges_flipped;

    using VertAttCol = wmtk::AttributeCollection<VertexAttributes>;
    VertAttCol vertex_attrs;

    explicit DelaunayEdgeFlip(
        std::vector<Eigen::Vector3d> _m_vertex_positions,
        const std::vector<std::array<size_t, 3>>& tris,
        int num_threads = 1);

    ~DelaunayEdgeFlip() override = default;

    bool invariants(const std::vector<Tuple>& new_tris) override;

    bool swap_edge_before(const Tuple& t) override;
    bool swap_edge_after(const Tuple& t) override;

    // Swap all edges in the mesh
    void swap_all_edges();

    bool write_triangle_mesh(std::string path);

private:
    
    std::vector<std::pair<std::string, Tuple>> collect_all_ops;

    // Use to save positions of an edge before
    // the edge is flipped.
    struct PositionInfoCache
    {
        Eigen::Vector3d v0p;
        Eigen::Vector3d v1p;
    };
    tbb::enumerable_thread_specific<PositionInfoCache> edge_position_cache;
};

} // namespace app::def
