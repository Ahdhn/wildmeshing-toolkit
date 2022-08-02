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
    is_cross_c = false;
    is_cross_d = false;

    // Check the bounding box first
    for (int j = 0; j < 2; j++) {
        if (a.pos[j] < b.pos[j]) {
            if (c.pos[j] < d) {
                if (c.pos[j] > b.pos[j]) return false;
            } else {
                if (d > b.pos[j]) return false;
            }
        } else {
            if (c.pos[j] < d) {
                if (c.pos[j] > a.pos[j]) return false;
            } else {
                if (d > a.pos[j]) return false;
            }
        }
    }

    // If the bounding box intersects, compute the intersection
    auto& x1 = a[0];
    auto& y1 = a[1];
    auto& x2 = b.pos[0];
    auto& y2 = b.pos[1];
    auto& x3 = c.pos[0];
    auto& y3 = c.pos[1];
    auto& x4 = d.pos[0];
    auto& y4 = d.pos[1];

    Rational td = (x4 - x3) * (y1 - y2) - (x1 - x2) * (y4 - y3);
    int td_sign = td.get_sign();
    if (td_sign == 0) return false;
    Rational tn = (y3 - y4) * (x1 - x3) + (x4 - x3) * (y1 - y3);
    if (tn.get_sign() * td_sign < 0 || tn > td) return false;
    tn = (y1 - y2) * (x1 - x3) + (x2 - x1) * (y1 - y3);
    if (tn.get_sign() * td_sign < 0 || tn > td) return false;

    if (tn.get_sign() == 0) {
        is_cross_p1 = true;
        return true;
    } else if (tn == td) {
        is_cross_p2 = true;
        return true;
    }

    auto& t = tn;
    t = tn / td;

    intersection.pos[0] = x3 + t * (x4 - x3);
    intersection.pos[1] = y3 + t * (y4 - y3);

    intersection.pos[0].canonicalize();
    intersection.pos[1].canonicalize();

    return true;
}

} // namespace bsp2d
} // namespace wmtk
