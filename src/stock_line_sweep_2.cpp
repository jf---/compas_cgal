#include "stock_line_sweep_2.h"

#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/Fraction_traits.h>

#include <cmath>
#include <iterator>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

using FT = Epeck::FT;
using Arr = Gps::Arrangement_2;

// Exact rational square root, if one exists: value = p/q in lowest terms is the
// square of a rational iff p and q are both perfect squares. Decided by CGAL's
// own number-type traits on the exact numerator and denominator.
bool rational_square_root(const FT& value, FT& root)
{
    using Traits = CGAL::Fraction_traits<FT>;
    static_assert(std::is_same_v<typename Traits::Is_fraction, CGAL::Tag_true>,
                  "the stock kernel's field type must decompose into a fraction");
    using Numerator = typename Traits::Numerator_type;
    using Denominator = typename Traits::Denominator_type;
    static_assert(std::is_same_v<Numerator, Denominator>,
                  "numerator and denominator share one integer type");
    typename CGAL::Algebraic_structure_traits<Numerator>::Is_square is_square;

    if (CGAL::sign(value) == CGAL::NEGATIVE) {
        return false;
    }
    Numerator numerator;
    Denominator denominator;
    typename Traits::Decompose()(value, numerator, denominator);
    Numerator numerator_root;
    Denominator denominator_root;
    if (!is_square(numerator, numerator_root) || !is_square(denominator, denominator_root)) {
        return false;
    }
    root = typename Traits::Compose()(numerator_root, denominator_root);
    return true;
}

GpsPolygon linear_polygon(const std::vector<EPoint>& vertices)
{
    GpsPolygon polygon;
    for (std::size_t index = 0; index < vertices.size(); ++index) {
        polygon.push_back(GpsXCurve(vertices[index], vertices[(index + 1) % vertices.size()]));
    }
    return polygon;
}

// disk(A, r) U rect(A -/+ h, B -/+ h) U disk(B, r): the two EXACT end disks of
// the tool's endpoint positions joined with the rectangle along the segment.
// With |h| = r (Pythagorean direction) the rectangle's sides are tangent to the
// disks and the union IS the capsule; otherwise |h| is a certified band about r
// and only the side band is inexact. One aggregated join, built in place.
void line_sweep_into(Gps& target, const EPoint& start, const EPoint& end,
                     const FT& tool_radius, const EVector& half_width)
{
    const std::vector<EPoint> corners = {
        start - half_width,
        end - half_width,
        end + half_width,
        start + half_width,
    };
    // Counterclockwise by construction (h is the CCW perpendicular of the
    // segment direction, scaled by a positive rational); checked exactly, as a
    // degenerate or clockwise general polygon corrupts the boolean set silently.
    if (CGAL::orientation(corners[0], corners[1], corners[2]) != CGAL::LEFT_TURN) {
        throw DegenerateLineSweepRectangleError(
            "line sweep rectangle corners are not counterclockwise; the rectangle is degenerate.");
    }
    std::vector<GpsPolygon> parts;
    parts.reserve(3);
    parts.push_back(disk_polygon(start, tool_radius));
    parts.push_back(disk_polygon(end, tool_radius));
    parts.push_back(linear_polygon(corners));
    target.join(parts.begin(), parts.end());
}

// A certified rectangle half-width vector along `unit_normal`: the double scale is
// a CONSTRUCTION parameter aimed at `aim * radius`, and what is trusted is only
// the exact band check lower <= |h| <= upper, compared as squares.
EVector certified_half_width(const EVector& unit_normal,
                             double unit_length,
                             double radius,
                             double aim,
                             const FT& lower,
                             const FT& upper)
{
    const double scale = aim * radius / unit_length;
    if (!std::isfinite(scale) || scale <= 0.0) {
        throw LineSweepBandCertificateError(
            "line sweep rectangle scale is not a positive finite number (radius="
            + std::to_string(radius) + ").");
    }
    const EVector half_width = FT(scale) * unit_normal;
    const FT squared = half_width.squared_length();
    if (CGAL::compare(squared, upper * upper) == CGAL::LARGER
        || CGAL::compare(squared, lower * lower) == CGAL::SMALLER) {
        throw LineSweepBandCertificateError(
            "line sweep rectangle half-width falls outside its certified band.");
    }
    return half_width;
}

// The curve of `halfedge`, directed along the halfedge: the stored edge curve
// when its direction already agrees, its exact opposite otherwise.
GpsXCurve curve_along(const Arr::Halfedge_const_handle& halfedge, const GpsTraits& traits)
{
    const bool halfedge_right = halfedge->direction() == CGAL::ARR_LEFT_TO_RIGHT;
    const GpsXCurve& stored = halfedge->curve();
    if (stored.is_directed_right() == halfedge_right) {
        return stored;
    }
    return traits.construct_opposite_2_object()(stored);
}

} // namespace

bool is_pythagorean_direction(const EVector& along)
{
    FT length;
    return rational_square_root(along.squared_length(), length);
}

LineSweepRegion2::LineSweepRegion2(EPoint start, EPoint end, Epeck::FT tool_radius,
                                   LineSweepModel2 model, EVector under_half_width,
                                   EVector over_half_width)
    : start_(std::move(start)),
      end_(std::move(end)),
      tool_radius_(std::move(tool_radius)),
      model_(model),
      under_half_width_(std::move(under_half_width)),
      over_half_width_(std::move(over_half_width))
{
}

LineSweepRegion2 LineSweepRegion2::build(double x0, double y0, double x1, double y1, double radius)
{
    // Finiteness is a DOUBLE concept, checked at the boundary before anything is
    // injected; the radius sign is then decided exactly on the injected value.
    if (!std::isfinite(x0) || !std::isfinite(y0) || !std::isfinite(x1)
        || !std::isfinite(y1) || !std::isfinite(radius)) {
        throw NonFiniteCapsuleInputError(
            "line sweep requires finite endpoints and radius; got x0=" + std::to_string(x0)
            + ", y0=" + std::to_string(y0) + ", x1=" + std::to_string(x1)
            + ", y1=" + std::to_string(y1) + ", radius=" + std::to_string(radius) + ".");
    }
    const FT tool_radius(radius);
    if (CGAL::sign(tool_radius) != CGAL::POSITIVE) {
        throw InvalidLineSweepRadiusError(
            "line sweep requires a positive tool radius; got " + std::to_string(radius) + ".");
    }

    // Each double is injected as itself -- no snapping, no tolerance at the seam.
    const EPoint start(x0, y0);
    const EPoint end(x1, y1);
    if (start == end) {
        const EVector none(0, 0);
        return LineSweepRegion2(start, end, tool_radius, LineSweepModel2::DISK, none, none);
    }

    const EVector along = end - start;
    FT length;
    if (rational_square_root(along.squared_length(), length)) {
        // |h| = r exactly: the region IS the swept capsule.
        const EVector half_width = (tool_radius / length) * along.perpendicular(CGAL::COUNTERCLOCKWISE);
        return LineSweepRegion2(start, end, tool_radius, LineSweepModel2::EXACT_CAPSULE,
                                half_width, half_width);
    }

    // Canonical direction of the class: divide by the larger |component|, an
    // exact rational with one coordinate +/-1. Every collinear sub-segment of
    // either travel direction maps to +/- the same u, hence the same h.
    const FT abs_x = CGAL::abs(along.x());
    const FT abs_y = CGAL::abs(along.y());
    const FT norm_infinity = (CGAL::compare(abs_x, abs_y) == CGAL::SMALLER) ? abs_y : abs_x;
    const EVector unit = along / norm_infinity;
    const EVector unit_normal = unit.perpendicular(CGAL::COUNTERCLOCKWISE);
    // Construction parameter: |u| in doubles, rounded from the EXACT value, so
    // the double is a deterministic function of the direction class and not of
    // the lazy kernel's interval history. Its error is answered by the exact
    // band check, never trusted.
    const double unit_length = std::hypot(
        CGAL::to_double(CGAL::exact(unit.x())),
        CGAL::to_double(CGAL::exact(unit.y())));
    const FT band(LINE_SWEEP_BAND_FRACTION);
    const EVector under = certified_half_width(
        unit_normal, unit_length, radius, 1.0 - 0.5 * LINE_SWEEP_BAND_FRACTION,
        tool_radius * (FT(1) - band), tool_radius);
    const EVector over = certified_half_width(
        unit_normal, unit_length, radius, 1.0 + 0.5 * LINE_SWEEP_BAND_FRACTION,
        tool_radius, tool_radius * (FT(1) + band));
    return LineSweepRegion2(start, end, tool_radius, LineSweepModel2::CANONICAL_QUAD,
                            under, over);
}

void LineSweepRegion2::under_region_into(Gps& target) const
{
    if (model_ == LineSweepModel2::DISK) {
        target.insert(disk_polygon(start_, tool_radius_));
        return;
    }
    line_sweep_into(target, start_, end_, tool_radius_, under_half_width_);
}

void LineSweepRegion2::over_region_into(Gps& target) const
{
    if (model_ == LineSweepModel2::DISK) {
        target.insert(disk_polygon(start_, tool_radius_));
        return;
    }
    line_sweep_into(target, start_, end_, tool_radius_, over_half_width_);
}

LocalRemovalRegion LineSweepRegion2::local_under_region() const
{
    // Build R exactly in a scratch set -- constant size, at most six input
    // curves -- and read its TRUE boundary off the arrangement, already split
    // at the exact (possibly one-root) junctions of caps and sides.
    Gps scratch;
    under_region_into(scratch);
    if (scratch.number_of_polygons_with_holes() != 1) {
        throw LineSweepRegionTopologyError(
            "line sweep region is not one polygon.");
    }
    std::vector<GpsPolygonWithHoles> polygons;
    scratch.polygons_with_holes(std::back_inserter(polygons));
    if (polygons.front().number_of_holes() != 0) {
        throw LineSweepRegionTopologyError(
            "line sweep region has a hole.");
    }

    LocalRemovalRegion region;
    const Arr& arrangement = scratch.arrangement();
    for (auto edge = arrangement.edges_begin(); edge != arrangement.edges_end(); ++edge) {
        // The removal keeps the material OUTSIDE R, so each boundary curve is
        // traversed with R's exterior on its left (the MATERIAL_OUTSIDE rule,
        // docs/local_depletion.md): take the halfedge whose face is the
        // exterior and orient the curve along it.
        Arr::Halfedge_const_handle outward = edge;
        if (outward->face()->contained()) {
            outward = outward->twin();
        }
        if (outward->face()->contained() || !outward->twin()->face()->contained()) {
            throw LineSweepRegionTopologyError(
                "line sweep region edge does not separate its interior from its exterior.");
        }
        region.boundary_curves.push_back(curve_along(outward, scratch.traits()));
    }
    // The exact midpoint of A and B is at distance 0 < |h| from the axis, hence
    // strictly interior, and the region's interior is connected: one seed.
    const EPoint middle = CGAL::midpoint(start_, end_);
    region.seeds.emplace_back(middle.x(), middle.y());
    return region;
}

// --- Stock2 entry points -----------------------------------------------------

void Stock2::subtract_line_sweep(double x0, double y0, double x1, double y1, double radius)
{
    const LineSweepRegion2 sweep = LineSweepRegion2::build(x0, y0, x1, y1, radius);
    Gps region;
    sweep.under_region_into(region);
    set().difference(region);
}

void Stock2::subtract_line_sweep_over(double x0, double y0, double x1, double y1, double radius)
{
    const LineSweepRegion2 sweep = LineSweepRegion2::build(x0, y0, x1, y1, radius);
    Gps region;
    sweep.over_region_into(region);
    set().difference(region);
}

void Stock2::subtract_line_sweep_local(double x0, double y0, double x1, double y1, double radius)
{
    const LineSweepRegion2 sweep = LineSweepRegion2::build(x0, y0, x1, y1, radius);
    subtract_region_local(set(), sweep.local_under_region());
}

SweepRemovalVerdict Stock2::line_sweep_removal(double x0, double y0, double x1, double y1,
                                               double radius) const
{
    const LineSweepRegion2 sweep = LineSweepRegion2::build(x0, y0, x1, y1, radius);
    Gps under;
    sweep.under_region_into(under);
    // Regularized interior intersection (Gps_do_intersect_functor): what the
    // replay's own region removes, decided exactly.
    if (set_->do_intersect(under)) {
        return SweepRemovalVerdict::REMOVES;
    }
    if (sweep.is_exact()) {
        return SweepRemovalVerdict::NOTHING;
    }
    Gps over;
    sweep.over_region_into(over);
    return set_->do_intersect(over) ? SweepRemovalVerdict::WITHIN_BAND
                                    : SweepRemovalVerdict::NOTHING;
}
