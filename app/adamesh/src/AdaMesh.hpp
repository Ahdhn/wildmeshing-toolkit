#pragma once

#include <tetwild/TetWild.h>

#include <Eigen/Core>
#include <stdexcept>
#include <tetwild/Rational.hpp>

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
    void insert_all_triangles(const std::vector<std::array<size_t, 3>>& tris);

public: // callbacks
};
} // namespace wmtk