#pragma once

#include <wmtk/TetMesh.h>

#include <Eigen/Core>
#include <memory>

namespace harmonic_tet {

struct VertexAttributes
{
    Eigen::Vector3d m_pos;
};

struct TetAttributes
{
    double m_qualities = 0.;
};

class HarmonicTet : public wmtk::TetMesh
{
public:
    HarmonicTet(
        const std::vector<VertexAttributes>& _vertex_attribute,
        const std::vector<TetAttributes>& _tet_attribute)
    {
        m_vertex_attribute = _vertex_attribute;
        m_tet_attribute = _tet_attribute;
    }

    ////// Attributes related
    // Stores the attributes attached to simplices
    std::vector<VertexAttributes> m_vertex_attribute;
    std::vector<TetAttributes> m_tet_attribute;

    void resize_attributes(size_t v, size_t e, size_t f, size_t t) override
    {
        m_vertex_attribute.resize(v);
        m_tet_attribute.resize(t);
    }

    void move_tet_attribute(size_t from, size_t to) override
    {
        m_tet_attribute[to] = std::move(m_tet_attribute[from]);
    }
    void move_vertex_attribute(size_t from, size_t to) override
    {
        m_vertex_attribute[to] = std::move(m_vertex_attribute[from]);
    }

    void output_mesh(std::string file);

    ////// Operations

    struct SwapInfoCache
    {
        double max_energy = 0.;
    } edgeswap_cache, faceswap_cache; // todo: change for parallel

    void smooth_all_vertices();
    bool smooth_before(const Tuple& t) override;
    bool smooth_after(const Tuple& t) override;

    void swap_all_edges();
    bool swap_edge_before(const Tuple& t) override;
    bool swap_edge_after(const Tuple& t) override;

    void swap_all_faces();
    bool swap_face_before(const Tuple& t) override;
    bool swap_face_after(const Tuple& t) override;

    bool is_inverted(const Tuple& loc);
    double get_quality(const Tuple& loc);
};

class ElementInQueue
{
public:
    wmtk::TetMesh::Tuple edge;
    double weight;

    ElementInQueue() {}
    ElementInQueue(const wmtk::TetMesh::Tuple& e, double w)
        : edge(e)
        , weight(w)
    {}
};

} // namespace harmonic_tet
