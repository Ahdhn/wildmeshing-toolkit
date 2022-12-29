#include "DelaunayEdgeFlip.h"

using namespace app::def;
using namespace wmtk;

DelaunayEdgeFlip::DelaunayEdgeFlip(
    std::vector<Eigen::Vector3d> _m_vertex_positions,
    int num_threads)
{
    NUM_THREADS = num_threads;
    p_vertex_attrs = &vertex_attrs;

    vertex_attrs.resize(_m_vertex_positions.size());

    for (auto i = 0; i < _m_vertex_positions.size(); i++)
        vertex_attrs[i] = {_m_vertex_positions[i], 0};
}
