#pragma once

#include <array>
#include <cstddef>
#include <vector>
#include <wmtk/utils/Rational.hpp>

namespace wmtk {
namespace bsp2d {

typedef size_t Index;
typedef std::array<Rational, 2> Point;
typedef std::array<Index, 3> Triangle;
typedef std::array<Index, 2> Edge;

/**
 * Compute the exact intersection of a segment ab and a segment cd
 * @note Only handles non-degenerate configurations, it return false if the segments are parallel
 * @param[a] first point of segment ab
 * @param[b] second point of segment ab
 * @param[c] first point of segment cd
 * @param[d] second point of segment cd
 * @param[intersection] the unique intersection if it exists
 *
 * @returns true if the two segments intersect and are not parallel
 */
bool segment_intersection(
    const Point& a,
    const Point& b,
    const Point& c,
    const Point& d,
    bool& is_cross_c,
    bool& is_cross_d,
    Point& intersection);

/**
 * Checks if the point p is contained in the closed segment ab
 * @param[a] first point of segment ab
 * @param[b] second point of segment ab
 * @param[p] the point to check containment for
 *
 * @returns true if the p is in the segment ab
 */
bool in_segment(
    const Point& a,
    const Point& b,
    const Point& p);

/**
 * Inserts the edges of E inside the triangulation with vertice V and triangles F
 * @param[V] rational vertex coordinates
 * @param[F] triangles
 * @param[E] edges to insert
 *
 * @returns A new set of triangles, and an unsorted map from every edge to insert to all its vertices in the triangulation
 */
std::pair<std::vector<Triangle>&,std::vector<std::vector<Index>>> bsp_subdivision(
    const std::vector<Point>& V, 
    const std::vector<Triangle>& F, 
    const std::vector<Edge>& E);

} // namespace bsp2d
} // namespace wmtk
