#include <wmtk/utils/BSP2D.hpp>

#include <catch2/catch.hpp>

#include <limits>

namespace {
inline wmtk::Rational compute_tri_area(
    const wmtk::bsp2d::Point& a,
    const wmtk::bsp2d::Point& b,
    const wmtk::bsp2d::Point& c)
{
    auto a0=a[0];
    auto a1=a[1];
    auto b0=b[0];
    auto b1=b[1];
    auto c0=c[0];
    auto c1=c[1];

    return (-(a1*b0) + a0*b1 + a1*c0 - b1*c0 - a0*c1 + b0*c1)/2.0;
}

} // namespace


TEST_CASE("BSP2D", "[BSP][2d]")
{
    using namespace wmtk;

    auto validate = [](const auto& vertices, const auto& triangles) {
        constexpr double EPS = std::numeric_limits<double>::epsilon();
        const size_t num_vertices = vertices.size();
        const size_t num_triangles = triangles.size();

        for (size_t i = 0; i < num_triangles; i++) {
            const auto& tri = triangles[i];
            REQUIRE(tri[0] < num_vertices);
            REQUIRE(tri[1] < num_vertices);
            REQUIRE(tri[2] < num_vertices);

            // tri must be positive oriented and non-degenerate.
            REQUIRE(
                compute_tri_area( 
                    vertices[tri[0]],
                    vertices[tri[1]],
                    vertices[tri[2]]) > EPS);
        }
    }

    SECTION("Intersection")
    {
        wmtk::bsp2d::Point a = {Rational(0),Rational(0)};
        wmtk::bsp2d::Point b = {Rational(0),Rational(1)};
        wmtk::bsp2d::Point c = {Rational(1),Rational(1)};
        wmtk::bsp2d::Point d = {Rational(1),Rational(0)};

        bool is_cross_c;
        bool is_cross_d;
        wmtk::bsp2d::Point intersection;
        bool res;

        res = wmtk::bsp2d::segment_intersection(a,b,c,d,is_cross_c,is_cross_d,intersection);
        REQUIRE(((res == false) && (is_cross_c == false) && (is_cross_d == false)));
        res = wmtk::bsp2d::segment_intersection(a,b,b,c,is_cross_c,is_cross_d,intersection);
        REQUIRE(((res == true) && (is_cross_c == true) && (is_cross_d == false) && (intersection == b)));
    }

    SECTION("Simple")
    {
        //std::vector<Point2D> points{{{0, 0}}, {{1, 0}}, {{1, 1}}, {{0, 1}}};

        // auto [vertices, triangles] = delaunay2D(points);
        // REQUIRE(vertices.size() == 4);
        // REQUIRE(triangles.size() == 2);
        // validate(vertices, triangles);
        REQUIRE(true);
    }

    // SECTION("Insufficient points should not fail")
    // {
    //     std::vector<Point2D> points{{{1, 0}}, {{0, 1}}};

    //     auto [vertices, triangles] = delaunay2D(points);
    //     REQUIRE(triangles.size() == 0);
    //     validate(vertices, triangles);
    // }
    // SECTION("Coplanar pts")
    // {
    //     std::vector<Point2D> points{{{0, 0}}, {{1, 0}}, {{2, 0}}, {{3, 0}}};

    //     auto [vertices, triangles] = delaunay2D(points);
    //     REQUIRE(vertices.size() == 4);
    //     REQUIRE(triangles.size() == 0);
    //     validate(vertices, triangles);
    // }
    // SECTION("Duplicate pts")
    // {
    //     std::vector<Point2D>
    //         points{{{0, 0}}, {{1, 0}}, {{0, 1}}, {{1, 1}}, {{0, 1}}, {{0, 1}}};

    //     auto [vertices, triangles] = delaunay2D(points);
    //     REQUIRE(vertices.size() == 6); // duplicate pts are kept in the output.
    //     REQUIRE(triangles.size() == 2);
    //     validate(vertices, triangles);
    // }
    // SECTION("Square with centroid")
    // {
    //     std::vector<Point2D> points{
    //         {{0, 0}},
    //         {{1, 0}},
    //         {{1, 1}},
    //         {{0, 1}},
    //         {{0.5, 0.5}}
    //     };

    //     auto [vertices, triangles] = delaunay2D(points);
    //     REQUIRE(vertices.size() == 5);
    //     REQUIRE(triangles.size() == 4);
    //     validate(vertices, triangles);
    // }
    // SECTION("Regular grid")
    // {
    //     constexpr size_t N = 10;
    //     std::vector<Point2D> points;
    //     points.reserve(N * N);
    //     for (size_t i = 0; i < N; i++) {
    //         for (size_t j = 0; j < N; j++) {
    //             double x = i;
    //             double y = j;
    //             points.push_back({{x, y}});
    //         }
    //     }
    //     auto [vertices, triangles] = delaunay2D(points);
    //     REQUIRE(vertices.size() == N * N);
    //     REQUIRE(triangles.size() == (N - 1) * (N - 1) * 2);
    //     validate(vertices, triangles);
    // }
}
