#pragma once

#include <tetwild/TetWild.h>

#include <Eigen/Core>
#include <stdexcept>
#include <tetwild/Rational.hpp>

#include "TetOperations.hpp"

using RowMatd = Eigen::Matrix<double, -1, -1, Eigen::RowMajor>;
using RowMati = Eigen::Matrix<int, -1, -1, Eigen::RowMajor>;

namespace wmtk {


using Vector3r = Eigen::Matrix<apps::Rational, 3, 1>;
using Vector2r = Eigen::Matrix<apps::Rational, 2, 1>;
using Eigen::Vector3d;
struct AdaMesh : public tetwild::TetWild
{
public:
    AdaMesh(tetwild::Parameters& params, fastEnvelope::FastEnvelope& envelope)
        : tetwild::TetWild(params, envelope, 1)
    {}
    void insert_all_points(
        const std::vector<Vector3d>& points,
        std::vector<int>& new_vid); // insert points

    void insert_triangles_to_mesh(
        const std::vector<Eigen::Vector3d>& vertices,
        const std::vector<std::array<size_t, 3>>& faces);

public:
    void serialize(std::string filename);
    void deserialize(std::string filename);
};
} // namespace wmtk


namespace adamesh {
struct SplitFace : public wmtk::SplitFace
{
    struct
    {
        Eigen::Vector3d pos;
        tetwild::FaceAttributes face_tag;
        std::map<std::array<size_t, 3>, tetwild::FaceAttributes> f_attrs;
    } cache;
    wmtk::AdaMesh& m_app;
    SplitFace(wmtk::AdaMesh& m_)
        : wmtk::SplitFace(m_)
        , m_app(m_)
    {}

    bool before(const wmtk::TetMesh::Tuple&);
    bool after(const std::vector<wmtk::TetMesh::Tuple>&);
};

struct DivideTet : public wmtk::DivideTet
{
    struct
    {
        Eigen::Vector3d pos;
        std::vector<std::pair<std::array<size_t, 3>, tetwild::FaceAttributes>> f_attrs;
        std::vector<int> old_fids;
    } cache;

    wmtk::AdaMesh& m_app;
    DivideTet(wmtk::AdaMesh& m_)
        : wmtk::DivideTet(m_)
        , m_app(m_)
    {}
    bool before(const wmtk::TetMesh::Tuple&);
    bool after(const std::vector<wmtk::TetMesh::Tuple>&);
};
} // namespace adamesh