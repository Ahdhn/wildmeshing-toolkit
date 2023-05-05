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

template <typename Scalar>
struct VertexAttributes
{
    Eigen::Vector3<Scalar> pos;
    size_t partition_id;
    bool freeze = false;
};

template <typename Scalar>
class DelaunayEdgeFlip : public wmtk::TriMesh
{
public:

    using VertAttCol = wmtk::AttributeCollection<VertexAttributes<Scalar>>;
    VertAttCol vertex_attrs;

    explicit DelaunayEdgeFlip(
        std::vector<Eigen::Vector3<Scalar>> _m_vertex_positions,
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
    // Use to save positions of an edge before
    // the edge is flipped.
    struct PositionInfoCache
    {
        Eigen::Vector3<Scalar> v0p;
        Eigen::Vector3<Scalar> v1p;
    };
    tbb::enumerable_thread_specific<PositionInfoCache> edge_position_cache;
};

} // namespace app::def

#include "DelaunayEdgeFlip_impl.h"
