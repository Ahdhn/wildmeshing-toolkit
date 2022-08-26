#include "BSP2D.hpp"


namespace wmtk {
namespace bsp2d {

bool segment_intersection(
    const Point& a,
    const Point& b,
    const Point& c,
    const Point& d,
    bool& is_cross_c,
    bool& is_cross_d,
    Point& intersection)
{
    
    // Check the bounding box first
    for (int j = 0; j < 2; j++) {
        if (a[j] < b[j]) {
            if (c[j] < d[j]) {
                if (c[j] > b[j]) return false;
            } else {
                if (d[j] > b[j]) return false;
            }
        } else {
            if (c[j] < d[j]) {
                if (c[j] > a[j]) return false;
            } else {
                if (d[j] > a[j]) return false;
            }
        }
    }

    // If the bounding box intersects, compute the intersection
    auto& x1 = a[0];
    auto& y1 = a[1];
    auto& x2 = b[0];
    auto& y2 = b[1];
    auto& x3 = c[0];
    auto& y3 = c[1];
    auto& x4 = d[0];
    auto& y4 = d[1];

    // Mathematica derivation
    // a = {x1, y1}
    // b = {x2, y2}
    // c = {x3, y3}
    // d = {x4, y4}
    // Solve[a + t1 (b - a) == c + t2 (d - c), {t1, t2}]

    Rational td = (x3*y1 - x4*y1 - x3*y2 + x4*y2 - x1*y3 + x2*y3 + x1*y4 - x2*y4);
    
    if (td.get_sign() == 0) // The two segments are parallel 
        return false;

    Rational tn = -(-(x3*y1) + x4*y1 + x1*y3 - x4*y3 - x1*y4 + x3*y4);
    Rational t = tn / td;
    if (t < 0 || t>1)
        return false;

    tn = -(-(x2*y1) + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3);
    t = tn / td;
    if (t < 0 || t>1)
        return false;

    intersection[0] = x3 + t * (x4 - x3);
    intersection[1] = y3 + t * (y4 - y3);

    intersection[0].canonicalize();
    intersection[1].canonicalize();

    is_cross_c = intersection == c;
    is_cross_d = intersection == d;

    return true;
}

} // namespace bsp2d
} // namespace wmtk
