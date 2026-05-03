///|/ Copyright (c) Prusa Research 2018 - 2021 Vojtěch Bubník @bubnikv, Lukáš Matěna @lukasmatena, Enrico Turri @enricoturri1966
///|/ Copyright (c) SuperSlicer 2018 - 2019 Remi Durand @supermerill
///|/
///|/ PrusaSlicer is released under the terms of the AGPLv3 or higher
///|/
#include "../ClipperUtils.hpp"
#include "../MarchingSquares.hpp"
#include "../ShortestPath.hpp"
#include "../Surface.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

#include "FillGyroid.hpp"

// ---------------------------------------------------------------------------
// Marching-squares scalar field for the Z-buckling-bias optimization.
// F(x,y,z) = sin(fx*x)cos(fy*y) + sin(fy*y)cos(fz*z) + sin(fz*z)cos(fx*x)
// fz = omega * baseline tightens the wave along the layer-stacking axis,
// shortening the effective vertical strand length and improving column-
// buckling resistance under Z-axis compression.
// ---------------------------------------------------------------------------
namespace marchsq {
using namespace Slic3r;

using coordr_t = long;
using Pointf   = Vec2d;

struct GyroidField
{
    static constexpr float gsizef = 0.40f;
    static constexpr float rsizef = 0.004f;
    const coord_t          rsize  = scaled(rsizef);
    const coordr_t         gsize  = std::round(gsizef / rsizef);
    Point                  size;
    Point                  offs;
    coordf_t               z;
    float                  fx;
    float                  fy;
    float                  fz;
    float                  isoval = 0.0f;

    explicit GyroidField(const BoundingBox bb, const coordf_t z, const float period, const float omega = 1.0f)
        : size{bb.size()}, offs{bb.min}, z{z}
    {
        const float baseline = float(2.0 * PI) / std::max(period, 1e-3f);
        fx = baseline;
        fy = baseline;
        fz = omega * baseline;
    }

    float get_scalar(coordf_t x, coordf_t y, coordf_t z_arg) const
    {
        const float a = fx * float(x);
        const float b = fy * float(y);
        const float c = fz * float(z_arg);
        return std::sin(a) * std::cos(b) + std::sin(b) * std::cos(c) + std::sin(c) * std::cos(a);
    }

    float get_scalar(Coord p) const
    {
        Pointf pf = to_Pointf(p);
        return get_scalar(pf.x(), pf.y(), z);
    }

    inline coord_t  to_coord (const coordr_t& x) const { return x * rsize; }
    inline coordr_t to_coordr(const coord_t& x)  const { return x / rsize; }
    inline Point  to_Point (const Coord& p) const { return Point(to_coord(p.c) + offs.x(), to_coord(p.r) + offs.y()); }
    inline Coord  to_Coord (const Point& p) const { return Coord(to_coordr(p.y() - offs.y()), to_coordr(p.x() - offs.x())); }
    inline Pointf to_Pointf(const Point& p) const { return Pointf(unscaled(p.x()), unscaled(p.y())); }
    inline Pointf to_Pointf(const Coord& p) const { return to_Pointf(to_Point(p)); }
};

template<> struct _RasterTraits<GyroidField>
{
    using ValueType = float;
    static float  get (const GyroidField& sf, size_t row, size_t col) { return sf.get_scalar(Coord(row, col)); }
    static size_t rows(const GyroidField& sf) { return sf.to_coordr(sf.size.y()); }
    static size_t cols(const GyroidField& sf) { return sf.to_coordr(sf.size.x()); }
};

inline Polylines get_gyroid_polylines(const GyroidField& sf, const double tolerance = SCALED_EPSILON)
{
    std::vector<Ring> rings = execute_with_policy(ex_tbb, sf, sf.isoval, {sf.gsize, sf.gsize});
    Polylines polys;
    polys.reserve(rings.size());
    for (const Ring& ring : rings) {
        Polyline poly;
        Points&  pts = poly.points;
        pts.reserve(ring.size() + 1);
        for (const Coord& crd : ring)
            pts.emplace_back(sf.to_Point(crd));
        pts.push_back(pts.front());
        if (tolerance >= 0.0)
            poly.simplify(tolerance);
        polys.emplace_back(poly);
    }
    return polys;
}

} // namespace marchsq

namespace Slic3r {

// Z-buckling bias: omega = sqrt(1 / density_adj) clamped [1, 2].
static inline double compute_omega_factor(double density_adjusted, double line_spacing, double layer_height)
{
    double lh_ratio   = (line_spacing > 0.) ? layer_height / line_spacing : 0.5;
    double correction = 1.0 / std::sqrt(1.0 + lh_ratio);
    double raw        = std::sqrt(1.0 / std::max(density_adjusted, 0.1)) * correction;
    return std::clamp(raw, 1.0, 2.0);
}

static inline double f(double x, double z_sin, double z_cos, bool vertical, bool flip)
{
    if (vertical) {
        double phase_offset = (z_cos < 0 ? M_PI : 0) + M_PI;
        double a   = sin(x + phase_offset);
        double b   = - z_cos;
        double res = z_sin * cos(x + phase_offset + (flip ? M_PI : 0.));
        double r   = sqrt(sqr(a) + sqr(b));
        return asin(a/r) + asin(res/r) + M_PI;
    } else {
        double phase_offset = z_sin < 0 ? M_PI : 0.;
        double a   = cos(x + phase_offset);
        double b   = - z_sin;
        double res = z_cos * sin(x + phase_offset + (flip ? 0 : M_PI));
        double r   = sqrt(sqr(a) + sqr(b));
        return (asin(a/r) + asin(res/r) + 0.5 * M_PI);
    }
}

static inline Polyline make_wave(
    const std::vector<Vec2d>& one_period, double width, double height, double offset, double scaleFactor,
    double z_cos, double z_sin, bool vertical, bool flip)
{
    std::vector<Vec2d> points = one_period;
    double period = points.back()(0);
    if (width != period) // do not extend if already truncated
    {
        points.reserve(one_period.size() * size_t(floor(width / period)));
        points.pop_back();

        size_t n = points.size();
        do {
            points.emplace_back(points[points.size()-n].x() + period, points[points.size()-n].y());
        } while (points.back()(0) < width - EPSILON);

        points.emplace_back(Vec2d(width, f(width, z_sin, z_cos, vertical, flip)));
    }

    // and construct the final polyline to return:
    Polyline polyline;
    polyline.points.reserve(points.size());
    for (auto& point : points) {
        point(1) += offset;
        point(1) = std::clamp(double(point.y()), 0., height);
        if (vertical)
            std::swap(point(0), point(1));
        polyline.points.push_back(Point::round(point * scaleFactor));
    }

    return polyline;
}

static std::vector<Vec2d> make_one_period(double width, double scaleFactor, double z_cos, double z_sin, bool vertical, bool flip, double tolerance)
{
    std::vector<Vec2d> points;
    double dx = M_PI_2; // exact coordinates on main inflexion lobes
    double limit = std::min(2*M_PI, width);
    points.reserve(size_t(ceil(limit / tolerance / 3)));

    for (double x = 0.; x < limit - EPSILON; x += dx) {
        points.emplace_back(Vec2d(x, f(x, z_sin, z_cos, vertical, flip)));
    }
    points.emplace_back(Vec2d(limit, f(limit, z_sin, z_cos, vertical, flip)));

    // piecewise increase in resolution up to requested tolerance
    for(;;)
    {
        size_t size = points.size();
        for (unsigned int i = 1;i < size; ++i) {
            auto& lp = points[i-1]; // left point
            auto& rp = points[i];   // right point
            double x = lp(0) + (rp(0) - lp(0)) / 2;
            double y = f(x, z_sin, z_cos, vertical, flip);
            Vec2d ip = {x, y};
            if (std::abs(cross2(Vec2d(ip - lp), Vec2d(ip - rp))) > sqr(tolerance)) {
                points.emplace_back(std::move(ip));
            }
        }

        if (size == points.size())
            break;
        else
        {
            // insert new points in order
            std::sort(points.begin(), points.end(),
                      [](const Vec2d &lhs, const Vec2d &rhs) { return lhs(0) < rhs(0); });
        }
    }

    return points;
}

static Polylines make_gyroid_waves(coordf_t gridZ, coordf_t scaleFactor, double width, double height, double tolerance)
{

    //scale factor for 5% : 8 712 388
    // 1z = 10^-6 mm ?
    const double z     = gridZ / scaleFactor;
    const double z_sin = sin(z);
    const double z_cos = cos(z);

    bool vertical = (std::abs(z_sin) <= std::abs(z_cos));
    double lower_bound = 0.;
    double upper_bound = height;
    bool flip = true;
    if (vertical) {
        flip = false;
        lower_bound = -M_PI;
        upper_bound = width - M_PI_2;
        std::swap(width,height);
    }

    std::vector<Vec2d> one_period_odd = make_one_period(width, scaleFactor, z_cos, z_sin, vertical, flip, tolerance); // creates one period of the waves, so it doesn't have to be recalculated all the time
    flip = !flip;                                                                   // even polylines are a bit shifted
    std::vector<Vec2d> one_period_even = make_one_period(width, scaleFactor, z_cos, z_sin, vertical, flip, tolerance);
    Polylines result;

    for (double y0 = lower_bound; y0 < upper_bound + EPSILON; y0 += M_PI) {
        // creates odd polylines
        result.emplace_back(make_wave(one_period_odd, width, height, y0, scaleFactor, z_cos, z_sin, vertical, flip));
        // creates even polylines
        y0 += M_PI;
        if (y0 < upper_bound + EPSILON) {
            result.emplace_back(make_wave(one_period_even, width, height, y0, scaleFactor, z_cos, z_sin, vertical, flip));
        }
    }

    return result;
}

void FillGyroid::_fill_surface_single(
    const FillParams                &params, 
    unsigned int                     thickness_layers,
    const std::pair<float, Point>   &direction, 
    ExPolygon                        expolygon, 
    Polylines                       &polylines_out) const
{
    float infill_angle = direction.first;
    if(std::abs(infill_angle) >= EPSILON)
        expolygon.rotate(-infill_angle);

    BoundingBox bb = expolygon.contour.bounding_box();
    // Distance between the gyroid waves in scaled coordinates.
    coord_t     line_spacing = _line_spacing_for_density(params);
    // Density adjusted to have a good %of weight.
    line_spacing /= FillGyroid::DENSITY_ADJUST;

    // align bounding box to a multiple of our grid module
    bb.merge(align_to_grid(bb.min, Point(2*M_PI*line_spacing, 2*M_PI*line_spacing)));

    // tolerance in scaled units. clamp the maximum tolerance as there's
    // no processing-speed benefit to do so beyond a certain point
    const double tolerance = params.config->get_computed_value("resolution_internal") / unscaled(line_spacing);

    // generate pattern
    Polylines polylines;
    if (params.gyroid_optimized) {
        // Marching-squares iso-extraction on the gyroid implicit field with anisotropic Z.
        const double lh = (params.layer_height > 0.) ? double(params.layer_height) : double(this->get_spacing());
        const double density_adj = std::max(0.001, double(params.density) * FillGyroid::DENSITY_ADJUST);
        const double omega       = compute_omega_factor(density_adj, this->get_spacing(), lh);
        const float  period      = float(2.0 * M_PI) * float(this->get_spacing()) / float(density_adj);

        marchsq::GyroidField sf(bb, this->z, period, float(omega));
        polylines = marchsq::get_gyroid_polylines(sf, SCALED_EPSILON);
    } else {
        polylines = make_gyroid_waves(
            scale_d(this->z),
            coordf_t(line_spacing),
            ceil(bb.size()(0) / line_spacing) + 1.,
            ceil(bb.size()(1) / line_spacing) + 1.,
            tolerance);

        // shift the parametric output to the grid origin; marching squares already
        // emits absolute coords via GyroidField::to_Point so it skips this.
        for (Polyline &pl : polylines)
            pl.translate(bb.min);
    }

    polylines = intersection_pl(polylines, expolygon);

    if (! polylines.empty()) {
        // Remove very small bits, but be careful to not remove infill lines connecting thin walls!
        // The infill perimeter lines should be separated by around a single infill line width.
        const coordf_t minlength = scale_d(0.8 * this->get_spacing());
        polylines.erase(
            std::remove_if(polylines.begin(), polylines.end(), [minlength](const Polyline &pl) { return pl.length() < minlength; }),
            polylines.end());
    }

    if (! polylines.empty()) {
        ensure_valid(polylines, params.fill_resolution);
        // connect lines
        size_t polylines_out_first_idx = polylines_out.size();
        if (params.connection == icNotConnected){
            append(polylines_out, chain_polylines(polylines));
        } else {
            this->connect_infill(chain_polylines(polylines), expolygon, polylines_out, scale_t(this->get_spacing()), params);
        }
        // new paths must be rotated back
        if (std::abs(infill_angle) >= EPSILON) {
            for (auto it = polylines_out.begin() + polylines_out_first_idx; it != polylines_out.end(); ++ it)
                it->rotate(infill_angle);
        }
    }
}

} // namespace Slic3r
