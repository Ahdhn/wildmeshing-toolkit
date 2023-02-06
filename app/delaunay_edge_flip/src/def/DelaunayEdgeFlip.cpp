#include "DelaunayEdgeFlip.h"
#include <set>

using namespace app::def;
using namespace wmtk;

DelaunayEdgeFlip::DelaunayEdgeFlip(
    std::vector<Eigen::Vector3d> _m_vertex_positions,
    int num_threads)
{
    NUM_THREADS = num_threads;
    p_vertex_attrs = &vertex_attrs;

    vertex_attrs.resize(_m_vertex_positions.size());

    auto partition_ids = wmtk::partition_morton(_m_vertex_positions, NUM_THREADS);
    assert(partition_ids.size() == _m_vertex_positions.size());
    wmtk::logger().info("partition_ids: {}", partition_ids.size());

    std::set<unsigned long> partition_set;
    for (auto i = 0; i < _m_vertex_positions.size(); i++)
    {
        partition_set.insert(partition_ids[i]);
        //vertex_attrs[i] = {_m_vertex_positions[i], 0};
        vertex_attrs[i] = {_m_vertex_positions[i], partition_ids[i]};
    }

    for(auto pid : partition_set)
    {
        wmtk::logger().info("pid: {}\n", pid);
    }

    //SDP, let's go ahead and setup the partitioning here when ready?
    //wmtk::partition_TriMesh(*this, NUM_THREADS);
}
