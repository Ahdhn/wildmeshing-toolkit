#include "QuadricIntegral.h"

#include <wmtk/quadrature/ClippedQuadrature.h>
#include <wmtk/utils/autodiff.h>
#include "helpers.h"

#include <lagrange/utils/assert.h>
#include <lagrange/utils/triangle_area.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

#include <fstream>

namespace wmtk {

namespace {

struct QuadratureCache
{
    wmtk::Quadrature quad;
    wmtk::Quadrature tmp;
};

// Column 0 is the center pixel, and other points are neighbors in clockwise order.
using Stencil = Eigen::Matrix<float, 3, 9>;

Stencil get_stencil(const std::array<wmtk::Image, 3>& images, int x, int y)
{
    auto w = images[0].width();
    auto h = images[0].height();

    auto get_point = [&](int sx, int sy) -> Eigen::Vector3f {
        sx = std::clamp(sx, 0, w - 1);
        sy = std::clamp(sy, 0, h - 1);
        return Eigen::Vector3f(
            images[0].get_raw_image()(sy, sx),
            images[1].get_raw_image()(sy, sx),
            images[2].get_raw_image()(sy, sx));
    };

    Stencil stencil;
    stencil.col(0) = get_point(x, y);
    stencil.col(1) = get_point(x - 1, y - 1);
    stencil.col(2) = get_point(x, y - 1);
    stencil.col(3) = get_point(x + 1, y - 1);
    stencil.col(4) = get_point(x + 1, y);
    stencil.col(5) = get_point(x + 1, y + 1);
    stencil.col(6) = get_point(x, y + 1);
    stencil.col(7) = get_point(x - 1, y + 1);
    stencil.col(8) = get_point(x - 1, y);
    return stencil;
}

Stencil adjusted_stencil(const Stencil& p)
{
    Stencil q = p;
    q.col(1) = 0.25 * (p.col(0) + p.col(1) + p.col(8) + p.col(2));
    q.col(3) = 0.25 * (p.col(0) + p.col(2) + p.col(3) + p.col(4));
    q.col(5) = 0.25 * (p.col(0) + p.col(4) + p.col(5) + p.col(6));
    q.col(7) = 0.25 * (p.col(0) + p.col(6) + p.col(7) + p.col(8));
    q.col(2) = 0.5 * (p.col(1) + p.col(3));
    q.col(4) = 0.5 * (p.col(3) + p.col(5));
    q.col(6) = 0.5 * (p.col(5) + p.col(7));
    q.col(8) = 0.5 * (p.col(7) + p.col(1));
    return q;
}

// Return unnormalized vector to weight the resulting quadric by the area of the patch
Eigen::Vector3d normal_from_patch(const Stencil& patch)
{
    Eigen::Vector3d n = Eigen::Vector3d::Zero();
    for (int i = 1; i <= 8; ++i) {
        const Eigen::Vector3d p0 = patch.col(0).cast<double>();
        const Eigen::Vector3d p1 = patch.col(i).cast<double>();
        const Eigen::Vector3d p2 = patch.col((i % 8) + 1).cast<double>();
        n += (p1 - p0).cross(p2 - p0);
    }
    return n;
}

Eigen::Vector3d stdev_from_patch(const Stencil& patch)
{
    // Mean
    Eigen::Vector3d mean = patch.cast<double>().rowwise().mean();

    // Standard Deviation
    Eigen::Vector3d std = (patch.cast<double>().colwise() - mean).array().pow(2).rowwise().mean();

    return std;
}

Quadric<double> compute_pixel_point_quadric(const Stencil& patch)
{
    const Eigen::Vector3d p0 = patch.col(0).cast<double>();
    return Quadric<double>::point_quadric(p0);
}

Quadric<double> compute_pixel_normal_quadric(const Eigen::Vector3f& normal)
{
    Eigen::Vector3d n = normal.cast<double>();
    n = std::sqrt(n.norm()) * n.normalized();
    return Quadric<double>::point_quadric(n);
}

Quadric<double> compute_pixel_plane_quadric(
    const Stencil& patch,
    const Eigen::Vector3f& normal,
    double sigma_q,
    double sigma_n)
{
    const auto mean_q = patch.col(0).cast<double>();
    // const auto mean_n = normal_from_patch(patch);
    // const double stdev_q = stdev_from_patch(patch).norm();
    Eigen::Vector3d mean_n = normal.cast<double>();
    mean_n = std::sqrt(mean_n.norm()) * mean_n.normalized();

    // We want more uncertainty when the area is small, so we divide by the 3d area
    // sigma_n /= std::max(1e-6, mean_n.norm());
    sigma_n *= std::sqrt(normal.norm());

    // We want more uncertainty when the patch is small, so we divide by the 3d pos standard
    // deviation
    // sigma_q /= std::max(1e-6, stdev_q);
    sigma_q = 0;

    return Quadric<double>::probabilistic_plane_quadric(mean_q, mean_n, sigma_q, sigma_n);
}

Quadric<double> compute_pixel_triangle_quadric(const Stencil& patch, double sigma_q)
{
    Quadric<double> q;
    const double stdev_q = stdev_from_patch(patch).norm();
    // We want more uncertainty when the patch is small, so we divide by the standard deviation
    sigma_q /= std::max(1e-6, stdev_q);

    for (int i = 1; i <= 8; ++i) {
        const Eigen::Vector3d p0 = patch.col(0).cast<double>();
        const Eigen::Vector3d p1 = patch.col(i).cast<double>();
        const Eigen::Vector3d p2 = patch.col((i % 8) + 1).cast<double>();
        q += Quadric<double>::probabilistic_triangle_quadric(p0, p1, p2, sigma_q);
    }
    return q;
}

void set_coefficients_from_quadric(
    std::array<wmtk::Image, 10>& coeffs,
    int x,
    int y,
    const Quadric<double>& q)
{
    coeffs[0].set(y, x, q.A00);
    coeffs[1].set(y, x, q.A01);
    coeffs[2].set(y, x, q.A02);
    coeffs[3].set(y, x, q.A11);
    coeffs[4].set(y, x, q.A12);
    coeffs[5].set(y, x, q.A22);
    coeffs[6].set(y, x, q.b0);
    coeffs[7].set(y, x, q.b1);
    coeffs[8].set(y, x, q.b2);
    coeffs[9].set(y, x, q.c);
}

Quadric<double> quadric_from_coefficients(const std::array<wmtk::Image, 10>& coeffs, int x, int y)
{
    return Quadric<double>::from_coefficients(
        coeffs[0].get_raw_image()(y, x),
        coeffs[1].get_raw_image()(y, x),
        coeffs[2].get_raw_image()(y, x),
        coeffs[3].get_raw_image()(y, x),
        coeffs[4].get_raw_image()(y, x),
        coeffs[5].get_raw_image()(y, x),
        coeffs[6].get_raw_image()(y, x),
        coeffs[7].get_raw_image()(y, x),
        coeffs[8].get_raw_image()(y, x),
        coeffs[9].get_raw_image()(y, x));
}

Quadric<double> quadric_from_coefficients(const Eigen::Vector<float, 10>& coeffs)
{
    return Quadric<double>::from_coefficients(
        coeffs[0],
        coeffs[1],
        coeffs[2],
        coeffs[3],
        coeffs[4],
        coeffs[5],
        coeffs[6],
        coeffs[7],
        coeffs[8],
        coeffs[9]);
}

// Optimized implementation that switches between nearest and bilinear interpolation
template <bool Exact, typename QuadricsFunc, typename DisplacementFunc>
Quadric<double> get_quadric_per_triangle_adaptive(
    const std::array<wmtk::Image, 10>& images,
    const Eigen::Matrix<double, 3, 2, Eigen::RowMajor>& triangle_uv,
    tbb::enumerable_thread_specific<QuadratureCache>& m_cache,
    const int order,
    QuadricsFunc get_quadrics,
    DisplacementFunc get_displacement)
{
    Eigen::AlignedBox2d bbox_uv;
    for (const auto& p : triangle_uv.rowwise()) {
        bbox_uv.extend(p.transpose());
    }

    const int w = images[0].width();
    const int h = images[0].height();
    auto get_coordinate = [&](double u, double v) -> Eigen::Vector2i {
        auto x = u * static_cast<float>(w);
        auto y = v * static_cast<float>(h);
        const auto sx = std::clamp(static_cast<int>(x), 0, w - 1);
        const auto sy = std::clamp(static_cast<int>(y), 0, h - 1);
        return {sx, sy};
    };

    const auto min_pixel = get_coordinate(bbox_uv.min()(0), bbox_uv.min()(1));
    const auto max_pixel = get_coordinate(bbox_uv.max()(0), bbox_uv.max()(1));
    const Eigen::Vector2d pixel_size(1.0 / w, 1.0 / h);

    const double u1 = triangle_uv(0, 0);
    const double v1 = triangle_uv(0, 1);
    const double u2 = triangle_uv(1, 0);
    const double v2 = triangle_uv(1, 1);
    const double u3 = triangle_uv(2, 0);
    const double v3 = triangle_uv(2, 1);
    const double denom = ((v2 - v3) * (u1 - u3) + (u3 - u2) * (v1 - v3));
    if (denom < std::numeric_limits<double>::denorm_min()) {
        // Degenerate triangle
        // logger().warn("Degenerate 2d triangle detected");
        return {};
    }

    const double area_2d = lagrange::triangle_area_2d<double>(
        triangle_uv.row(0),
        triangle_uv.row(1),
        triangle_uv.row(2));
    // if (area_2d < std::numeric_limits<double>::epsilon()) {
    //     logger().warn("Small-area 2d triangle detected");
    //     return {};
    // }

    const std::array<Eigen::Hyperplane<double, 2>, 3> edges = {
        Eigen::Hyperplane<double, 2>::Through(triangle_uv.row(0), triangle_uv.row(1)),
        Eigen::Hyperplane<double, 2>::Through(triangle_uv.row(1), triangle_uv.row(2)),
        Eigen::Hyperplane<double, 2>::Through(triangle_uv.row(2), triangle_uv.row(0)),
    };

    Quadric<double> q;
    const double pixel_radius = pixel_size.norm() * 0.5;
    for (int y = min_pixel.y(); y <= max_pixel.y(); ++y) {
        for (int x = min_pixel.x(); x <= max_pixel.x(); ++x) {
            Eigen::Vector2i pixel_coord(x, y);
            Eigen::AlignedBox2d box;
            box.extend(pixel_coord.cast<double>().cwiseProduct(pixel_size));
            box.extend(
                (pixel_coord + Eigen::Vector2i::Ones()).cast<double>().cwiseProduct(pixel_size));
            internal::Classification sign = internal::Classification::Unknown;
            if constexpr (!Exact) {
                auto sign = internal::point_in_triangle_quick(edges, box.center(), pixel_radius);
                if (sign == internal::Classification::Unknown) {
                    sign = internal::pixel_inside_triangle(triangle_uv, box);
                }
                if (sign == internal::Classification::Outside) {
                    continue;
                }
            }
            if (sign == internal::Classification::Inside) {
                q += quadric_from_coefficients(images, x, y) * pixel_size(0) * pixel_size(1);
            } else {
                wmtk::ClippedQuadrature rules;
                auto& quadr = m_cache.local().quad;
                rules.clipped_triangle_box_quadrature(
                    order,
                    triangle_uv,
                    box,
                    quadr,
                    &m_cache.local().tmp);
                for (size_t i = 0; i < quadr.size(); ++i) {
                    float u = quadr.points()(i, 0);
                    float v = quadr.points()(i, 1);
                    q += quadric_from_coefficients(get_quadrics(u, v)) * quadr.weights()[i];
                }
            }
        }
    }

    // Compute area ratio
    Eigen::Vector3d p1 = get_displacement(u1, v1).template cast<double>();
    Eigen::Vector3d p2 = get_displacement(u2, v2).template cast<double>();
    Eigen::Vector3d p3 = get_displacement(u3, v3).template cast<double>();
    const double area_3d = lagrange::triangle_area_3d<double>(p1, p2, p3);
    if (q.A().isZero()) {
        logger().info(
            "Zero quadric detected for triangle {} / {} / {}",
            triangle_uv.row(0),
            triangle_uv.row(1),
            triangle_uv.row(2));
    }
    return q; // * area_3d / area_2d;
}

std::array<wmtk::Image, 3> compute_adjusted_normals(
    const std::array<wmtk::Image, 3>& displaced,
    const std::array<wmtk::Image, 3>& normals)
{
    // Compute surface normal a ru x rv where ru and rv are the partial derivatives obtained via
    // bicubic interpolation. Since we have non-regular points where the derivatives vanish, we
    // weight in the baked normals from the input coarse mesh. We use an exponential decay formula
    // to favor the normal calculated via derivatives.

    using DScalar = DScalar2<double, Eigen::Vector2d, Eigen::Matrix2d>;

    std::array<wmtk::Image, 3> adjusted_normals = normals;

    auto w = displaced[0].width();
    auto h = displaced[0].height();

    // First pass, compute ru x rv
    tbb::parallel_for(0, h, [&](int y) {
        DScalar::setVariableCount(2);
        for (int x = 0; x < w; ++x) {
            DScalar u(0, (static_cast<float>(x) + 0.5f) / w);
            DScalar v(1, (static_cast<float>(y) + 0.5f) / h);
            auto p = internal::sample_bicubic(displaced, u, v);
            Eigen::Vector3d ru;
            Eigen::Vector3d rv;
            for (size_t k = 0; k < 3; ++k) {
                ru[k] = p[k].getGradient()[0];
                rv[k] = p[k].getGradient()[1];
            }
            Eigen::Vector3d n = ru.cross(rv);
            for (size_t k = 0; k < 3; ++k) {
                adjusted_normals[k].set(y, x, n[k]);
            }
        }
    });

    // Second pass, compute |ru x rv| and determine average area
    wmtk::Image area(w, h);
    tbb::parallel_for(0, h, [&](int y) {
        for (int x = 0; x < w; ++x) {
            area.set(y, x, internal::fetch_texels(adjusted_normals, x, y).norm());
        }
    });
    const float lambda = std::log(2.0) / 1e-2 * area.get_raw_image().mean();

    // Third pass, adjust and normalize
    // tbb::parallel_for(0, h, [&](int y) {
    //     for (int x = 0; x < w; ++x) {
    //         const float w = std::exp(-lambda * internal::fetch_texel(area, x, y));
    //         Eigen::Vector3f nadj = internal::fetch_texels(adjusted_normals, x, y);
    //         Eigen::Vector3f nref =
    //             2.f * internal::fetch_texels(normals, x, y) - Eigen::Vector3f::Ones();
    //         Eigen::Vector3f n = ((1.f - w) * nadj + w * nref).stableNormalized();
    //         for (size_t k = 0; k < 3; ++k) {
    //             adjusted_normals[k].set(y, x, adj[k]);
    //         }
    //     }
    // });

    return adjusted_normals;
}

} // namespace

struct QuadricIntegral::Cache
{
    mutable tbb::enumerable_thread_specific<QuadratureCache> quadrature_cache;
};

QuadricIntegral::QuadricIntegral() = default;
QuadricIntegral::~QuadricIntegral() = default;
QuadricIntegral::QuadricIntegral(QuadricIntegral&&) = default;
QuadricIntegral& QuadricIntegral::operator=(QuadricIntegral&&) = default;

QuadricIntegral::QuadricIntegral(
    std::array<wmtk::Image, 3> displaced_positions,
    const std::array<wmtk::Image, 3>& normals,
    QuadricType quadric_type)
    : m_cache(lagrange::make_value_ptr<Cache>())
    , m_displaced(std::move(displaced_positions))
{
    // Default option for quadric integrator
    m_sampling_method = SamplingMethod::Bilinear;
    m_integration_method = IntegrationMethod::Adaptive;

    // Relative uncertainty on point positions
    const double sigma_q = 1e-4;

    // Relative uncertainty on normal directions
    const double sigma_n = 1e-3;

    auto w = m_displaced[0].width();
    auto h = m_displaced[0].height();
    for (auto& img : m_quadrics) {
        img = wmtk::Image(w, h);
    }

    auto adjusted_normals = compute_adjusted_normals(m_displaced, normals);

    // std::ofstream out("points.csv");
    // out << "x,y,z,nx,ny,nz" << std::endl;
    // for (int y = 0; y < h; ++y) {
    //     for (int x = 0; x < w; ++x) {
    //         auto p = internal::fetch_texels(m_displaced, x, y);
    //         auto n = internal::fetch_texels(adjusted_normals, x, y);
    //         out << p.x() << ", " << p.y() << ", " << p.z() << ", " << n.x() << ", " << n.y() <<
    //         ", "
    //             << n.z() << std::endl;
    //     }
    // }

    tbb::parallel_for(0, h, [&](int y) {
        for (int x = 0; x < w; ++x) {
            auto stencil = adjusted_stencil(get_stencil(m_displaced, x, y));
            auto normal = internal::fetch_texels(adjusted_normals, x, y);
            switch (quadric_type) {
            case QuadricType::Point: {
                auto quadric = normal.norm() * compute_pixel_point_quadric(stencil);
                set_coefficients_from_quadric(m_quadrics, x, y, quadric);
                break;
            }
            case QuadricType::Plane: {
                auto quadric = compute_pixel_plane_quadric(stencil, normal, sigma_q, sigma_n);
                set_coefficients_from_quadric(m_quadrics, x, y, quadric);
                break;
            }
            case QuadricType::Triangle: {
                auto quadric = compute_pixel_triangle_quadric(stencil, sigma_q);
                set_coefficients_from_quadric(m_quadrics, x, y, quadric);
                break;
            }
            case QuadricType::Normal: {
                auto quadric = compute_pixel_normal_quadric(normal);
                set_coefficients_from_quadric(m_quadrics, x, y, quadric);
                break;
            }
            default: break;
            }
        }
    });
}

template <QuadricIntegral::SamplingMethod Sampling, QuadricIntegral::IntegrationMethod Integration>
void QuadricIntegral::get_quadric_per_triangle_internal(
    int num_triangles,
    lagrange::function_ref<std::array<float, 6>(int)> get_triangle,
    lagrange::span<wmtk::Quadric<double>> output_quadrics,
    int order) const
{
    assert(num_triangles == output_quadrics.size());
    tbb::parallel_for(0, num_triangles, [&](int i) {
        auto input_triangle = get_triangle(i);
        Eigen::Matrix<double, 3, 2, Eigen::RowMajor> triangle;
        triangle.row(0) << input_triangle[0], input_triangle[1];
        triangle.row(1) << input_triangle[2], input_triangle[3];
        triangle.row(2) << input_triangle[4], input_triangle[5];
        auto sampling_quadrics = [&](float u, float v) -> Eigen::Matrix<float, 10, 1> {
            if constexpr (Sampling == SamplingMethod::Bicubic) {
                return internal::sample_bicubic(m_quadrics, u, v);
            } else if constexpr (Sampling == SamplingMethod::Nearest) {
                return internal::sample_nearest(m_quadrics, u, v);
            } else if constexpr (Sampling == SamplingMethod::Bilinear) {
                return internal::sample_bilinear(m_quadrics, u, v);
            }
        };
        auto sampling_displaced = [&](float u, float v) -> Eigen::Matrix<float, 3, 1> {
            if constexpr (Sampling == SamplingMethod::Bicubic) {
                return internal::sample_bicubic(m_displaced, u, v);
            } else if constexpr (Sampling == SamplingMethod::Nearest) {
                return internal::sample_nearest(m_displaced, u, v);
            } else if constexpr (Sampling == SamplingMethod::Bilinear) {
                return internal::sample_bilinear(m_displaced, u, v);
            }
        };

        if constexpr (Integration == IntegrationMethod::Exact) {
            output_quadrics[i] = get_quadric_per_triangle_adaptive<true>(
                m_quadrics,
                triangle,
                m_cache->quadrature_cache,
                order,
                sampling_quadrics,
                sampling_displaced);
        } else {
            output_quadrics[i] = get_quadric_per_triangle_adaptive<false>(
                m_quadrics,
                triangle,
                m_cache->quadrature_cache,
                order,
                sampling_quadrics,
                sampling_displaced);
        }
    });
}

void QuadricIntegral::get_quadric_per_triangle(
    int num_triangles,
    lagrange::function_ref<std::array<float, 6>(int)> get_triangle,
    lagrange::span<wmtk::Quadric<double>> output_quadrics) const
{
    int order = 1;
    if (m_quadrature_order == QuadratureOrder::Full) {
        constexpr int Degree = 2;
        order = 2 * (Degree - 1);
    }

    assert(num_triangles == output_quadrics.size());
    switch (m_sampling_method) {
    case SamplingMethod::Bicubic:
        if (m_integration_method == IntegrationMethod::Exact) {
            get_quadric_per_triangle_internal<SamplingMethod::Bicubic, IntegrationMethod::Exact>(
                num_triangles,
                get_triangle,
                output_quadrics,
                order);
        } else {
            get_quadric_per_triangle_internal<SamplingMethod::Bicubic, IntegrationMethod::Adaptive>(
                num_triangles,
                get_triangle,
                output_quadrics,
                order);
        }
        break;
    case SamplingMethod::Nearest:
        if (m_integration_method == IntegrationMethod::Exact) {
            get_quadric_per_triangle_internal<SamplingMethod::Nearest, IntegrationMethod::Exact>(
                num_triangles,
                get_triangle,
                output_quadrics,
                order);
        } else {
            get_quadric_per_triangle_internal<SamplingMethod::Nearest, IntegrationMethod::Adaptive>(
                num_triangles,
                get_triangle,
                output_quadrics,
                order);
        }
        break;
    case SamplingMethod::Bilinear:
        if (m_integration_method == IntegrationMethod::Exact) {
            get_quadric_per_triangle_internal<SamplingMethod::Bilinear, IntegrationMethod::Exact>(
                num_triangles,
                get_triangle,
                output_quadrics,
                order);
        } else {
            get_quadric_per_triangle_internal<
                SamplingMethod::Bilinear,
                IntegrationMethod::Adaptive>(num_triangles, get_triangle, output_quadrics, order);
        }
        break;
    }
}

} // namespace wmtk
