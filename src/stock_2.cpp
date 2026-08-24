#include "stock_2.h"
#include "audit_classification_2.h"
#include "engagement_2.h"
#include "stock_local_2.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>

namespace {

// Convert an Nx3 double matrix (rationals by construction) to a linear
// circle-segment general polygon, validating simplicity via Polygon_2.
GpsPolygon data_to_gps_polygon(Eigen::Ref<const compas::RowMatrixXd> vertices)
{
    if (vertices.rows() < 3) {
        throw std::invalid_argument("Expected at least three polygon vertices.");
    }
    CGAL::Polygon_2<Epeck> simple_check;
    for (int i = 0; i < vertices.rows(); ++i) {
        simple_check.push_back(EPoint(vertices(i, 0), vertices(i, 1)));
    }
    if (!simple_check.is_simple()) {
        throw std::invalid_argument("Polygon boundary must be simple (no self-intersections).");
    }
    if (simple_check.is_clockwise_oriented()) {
        simple_check.reverse_orientation();
    }

    GpsPolygon polygon;
    const std::size_t n = simple_check.size();
    for (std::size_t i = 0; i < n; ++i) {
        const EPoint& a = simple_check[i];
        const EPoint& b = simple_check[(i + 1) % n];
        // Vendored CGAL 6.0.1: X_monotone_curve_2 (_X_monotone_circle_segment_2)
        // has a two-kernel-point constructor but no Segment_2 constructor, so
        // build the linear arc directly from the exact endpoints.
        polygon.push_back(GpsXCurve(a, b));
    }
    return polygon;
}

GpsPolygon exact_linear_polygon(const std::vector<EPoint>& vertices)
{
    GpsPolygon polygon;
    for (std::size_t index = 0; index < vertices.size(); ++index) {
        polygon.push_back(GpsXCurve(
            vertices[index],
            vertices[(index + 1) % vertices.size()]));
    }
    return polygon;
}

Gps exact_disk_union(
    const std::vector<EPoint>& centers,
    const Epeck::FT& radius)
{
    std::vector<GpsPolygon> disks;
    disks.reserve(centers.size());
    for (const EPoint& center : centers) {
        disks.push_back(disk_polygon(center, radius));
    }
    Gps region;
    region.join(disks.begin(), disks.end());
    return region;
}

// The exact annulus between `inner_radius` and `outer_radius` about `center`.
// One builder, shared by the depletion path (Stock2::subtract_annulus_exact) and
// by the full-circle sweep ORACLE the depletion certificates are proved against,
// so "the region we remove" and "the region we prove under-coverage against" are
// the same construction by identity rather than by assertion.
//
// `inner_radius == 0` yields the plain disk: a guide no wider than the tool
// sweeps a filled disk, with no hole to punch. That branch is structural -- an
// exact sign test on an exact quantity -- not an epsilon.
Gps exact_annulus_region(
    const EPoint& center,
    const Epeck::FT& inner_radius,
    const Epeck::FT& outer_radius)
{
    Gps region;
    region.insert(disk_polygon(center, outer_radius));
    if (CGAL::sign(inner_radius) == CGAL::POSITIVE) {
        Gps hole;
        hole.insert(disk_polygon(center, inner_radius));
        region.difference(hole);
    }
    return region;
}

void validate_depletion_trace(const DepletionTrace& trace)
{
    if (trace.center_count != trace.center_parameters.size()
        || !trace.exact_incidence
        || !trace.exact_parameters_in_range
        || !trace.exact_anchors_present
        || !trace.exact_removal_radius_valid
        || !trace.exact_chord_bound_holds
        || (trace.cyclic && !trace.exact_seam_chord_bound_holds)
        || trace.strategy_version.empty()) {
        throw ExactDepletionConstructionError(
            "exact depletion structural trace failed validation.");
    }
}

std::vector<ExactCenterParameter2> to_exact_center_parameters(
    const std::vector<std::tuple<int, std::size_t, std::size_t>>& parameters)
{
    std::vector<ExactCenterParameter2> result;
    result.reserve(parameters.size());
    for (const auto& [chart, numerator, denominator] : parameters) {
        result.push_back({chart, numerator, denominator});
    }
    return result;
}

bool exact_set_is_subset(const Gps& subset, const Gps& superset)
{
    Gps difference(subset);
    difference.difference(superset);
    return difference.is_empty();
}

// Exact bounds of the annulus a tool of radius `tool_radius` sweeps over a FULL
// turn about a guide of radius `guide_r`: [max(rho - r, 0), rho + r]. Both
// arrive as doubles, hence as exact rationals, and the bounds are formed
// EXACTLY -- forming them in double arithmetic first would round rho +/- r to
// the nearest double and put the swept region about an ulp off the sweep oracle
// the depletion certificates are proved against. The rho <= r case (a guide no
// wider than the tool sweeps a filled disk) is an exact comparison, not a
// tolerance. Shared verbatim by the global and local arc-sweep paths, so the two
// remove the SAME region by construction rather than by assertion.
std::pair<Epeck::FT, Epeck::FT> full_turn_annulus_bounds(double guide_r, double tool_radius)
{
    const Epeck::FT guide(guide_r);
    const Epeck::FT tool(tool_radius);
    const Epeck::FT inner = (CGAL::compare(guide, tool) == CGAL::LARGER)
        ? Epeck::FT(guide - tool)
        : Epeck::FT(0);
    return {inner, guide + tool};
}

} // namespace

// Full disk as a two-arc general polygon (counterclockwise). Split at the two
// x-extreme (vertical-tangency) points; radius comes from a double so it is an
// exact rational, hence center.x() +/- radius are plain rationals that GpsPoint
// accepts directly (Task 1 pattern). This mirrors CGAL's own full-circle
// subdivision in Arr_circle_segment_traits_2::Make_x_monotone_2.
// Declared in stock_2.h with external linkage so engagement_2.cpp can build the
// same exact cutter disk it intersects the stock against.
GpsPolygon disk_polygon(const EPoint& center, const Epeck::FT& radius)
{
    const Epeck::FT r_sq = radius * radius;
    const ECircle circle(center, r_sq, CGAL::COUNTERCLOCKWISE);
    const GpsPoint p_min(center.x() - radius, center.y());  // leftmost point
    const GpsPoint p_max(center.x() + radius, center.y());  // rightmost point
    GpsPolygon polygon;
    // Vendored CGAL 6.0.1: X_monotone_curve_2(Circle_2, source, target,
    // Orientation) — verified against Arr_geometry_traits/Circle_segment_2.h.
    polygon.push_back(GpsXCurve(circle, p_min, p_max, CGAL::COUNTERCLOCKWISE));  // lower half
    polygon.push_back(GpsXCurve(circle, p_max, p_min, CGAL::COUNTERCLOCKWISE));  // upper half
    return polygon;
}

Stock2::Stock2(Eigen::Ref<const compas::RowMatrixXd> boundary,
               const std::vector<compas::RowMatrixXd>& holes)
    : set_(std::make_unique<Gps>())
{
    set_->insert(data_to_gps_polygon(boundary));
    for (const auto& hole : holes) {
        Gps hole_set;
        hole_set.insert(data_to_gps_polygon(hole));
        set_->difference(hole_set);
    }
}

Stock2::Stock2(std::unique_ptr<Gps> set) noexcept
    : set_(std::move(set))
{
}

bool Stock2::contains(double x, double y) const
{
    return set_->oriented_side(GpsPoint(Epeck::FT(x), Epeck::FT(y))) == CGAL::ON_POSITIVE_SIDE;
}

bool Stock2::is_empty() const
{
    return set_->is_empty();
}

Stock2 Stock2::clone() const
{
    return Stock2(std::make_unique<Gps>(*set_));
}

bool Stock2::is_subset_of(const Stock2& other) const
{
    return exact_set_is_subset(*set_, *other.set_);
}

bool Stock2::exactly_equals(const Stock2& other) const
{
    return exact_set_is_subset(*set_, *other.set_)
        && exact_set_is_subset(*other.set_, *set_);
}

// Fraction of the tool radius allowed as chain under-coverage slack; the
// generator's radial_clearance margin (1e-3 * D = 2e-3 * r) dominates it.
constexpr double CHAIN_SLACK_FRACTION = 1e-4;

void Stock2::subtract_disk(double cx, double cy, double radius)
{
    if (radius <= 0.0) throw std::invalid_argument("radius should be positive.");
    Gps region;
    region.insert(disk_polygon(EPoint(cx, cy), Epeck::FT(radius)));
    set_->difference(region);
}

// Exact swept region of a disk of radius r carried about a circular guide of
// radius rho: the annulus between rho - r and rho + r. This is the region
// EXACTLY -- an equality, not a bound -- so unlike the disk chain it carries no
// approximation to compensate for.
//
// It is representable because both radii arrive as doubles and every double IS a
// rational: disk_polygon squares them into the rational squared radii
// Gps_circle_segment_traits_2 requires. Two boundary circles (four x-monotone
// arcs) replace the hundreds of disks a chain would need over the same turn.
void Stock2::subtract_annulus(double cx, double cy, double inner_radius, double outer_radius)
{
    // Finiteness is a DOUBLE concept, so it is checked here, at the boundary,
    // before anything is injected: Epeck::FT(NaN) has no meaning to build on.
    if (!std::isfinite(cx) || !std::isfinite(cy)
        || !std::isfinite(inner_radius) || !std::isfinite(outer_radius)) {
        throw NonFiniteAnnulusInputError(
            "subtract_annulus requires finite center coordinates and radii; got "
            "cx=" + std::to_string(cx) + ", cy=" + std::to_string(cy)
            + ", inner_radius=" + std::to_string(inner_radius)
            + ", outer_radius=" + std::to_string(outer_radius) + ".");
    }
    // Each double is injected as itself -- no snapping, no tolerance at the seam.
    // The radius ORDERING is then decided exactly, downstream, by the exact core.
    subtract_annulus_exact(
        EPoint(cx, cy),
        Epeck::FT(inner_radius),
        Epeck::FT(outer_radius));
}

void Stock2::subtract_annulus_exact(
    const EPoint& center,
    const Epeck::FT& inner_radius,
    const Epeck::FT& outer_radius)
{
    if (CGAL::sign(inner_radius) == CGAL::NEGATIVE) {
        throw InvalidAnnulusRadiiError(
            "subtract_annulus requires inner_radius >= 0.");
    }
    if (CGAL::compare(outer_radius, inner_radius) != CGAL::LARGER) {
        throw InvalidAnnulusRadiiError(
            "subtract_annulus requires outer_radius > inner_radius.");
    }
    Gps region = exact_annulus_region(center, inner_radius, outer_radius);
    set_->difference(region);
}

// --- Local depletion ---------------------------------------------------------
// Same argument validation, same exact region, different removal mechanism: the
// arrangement is edited around the region instead of being rebuilt by an overlay
// against the whole stock. See stock_local_2.cpp for the invariants.

void Stock2::subtract_disk_local(double cx, double cy, double radius)
{
    if (radius <= 0.0) throw std::invalid_argument("radius should be positive.");
    subtract_region_local(*set_, local_disk_region(EPoint(cx, cy), Epeck::FT(radius)));
}

void Stock2::subtract_annulus_local(double cx, double cy,
                                    double inner_radius, double outer_radius)
{
    // Finiteness is a DOUBLE concept, checked at the boundary before anything is
    // injected -- the same seam contract subtract_annulus states.
    if (!std::isfinite(cx) || !std::isfinite(cy)
        || !std::isfinite(inner_radius) || !std::isfinite(outer_radius)) {
        throw NonFiniteAnnulusInputError(
            "subtract_annulus requires finite center coordinates and radii; got "
            "cx=" + std::to_string(cx) + ", cy=" + std::to_string(cy)
            + ", inner_radius=" + std::to_string(inner_radius)
            + ", outer_radius=" + std::to_string(outer_radius) + ".");
    }
    subtract_annulus_exact_local(
        EPoint(cx, cy),
        Epeck::FT(inner_radius),
        Epeck::FT(outer_radius));
}

void Stock2::subtract_annulus_exact_local(
    const EPoint& center,
    const Epeck::FT& inner_radius,
    const Epeck::FT& outer_radius)
{
    if (CGAL::sign(inner_radius) == CGAL::NEGATIVE) {
        throw InvalidAnnulusRadiiError(
            "subtract_annulus requires inner_radius >= 0.");
    }
    if (CGAL::compare(outer_radius, inner_radius) != CGAL::LARGER) {
        throw InvalidAnnulusRadiiError(
            "subtract_annulus requires outer_radius > inner_radius.");
    }
    subtract_region_local(*set_, local_annulus_region(center, inner_radius, outer_radius));
}

void Stock2::subtract_arc_sweep_local(double cx, double cy, double sx, double sy,
                                      double ex, double ey, bool cw, double tool_radius)
{
    if (tool_radius <= 0.0) throw std::invalid_argument("tool_radius should be positive.");
    const double guide_r = std::hypot(sx - cx, sy - cy);
    if (guide_r == 0.0) { subtract_disk_local(cx, cy, tool_radius); return; }
    if (sx == ex && sy == ey) {
        const auto [inner, outer] = full_turn_annulus_bounds(guide_r, tool_radius);
        subtract_annulus_exact_local(EPoint(cx, cy), inner, outer);
        return;
    }
    // Partial arc: the removed region is a disk-chain union, whose boundary the
    // local update cannot identify exactly (a chain circle carries both boundary
    // and interior arcs), so this defers to the global path rather than guessing.
    subtract_arc_sweep(cx, cy, sx, sy, ex, ey, cw, tool_radius);
}

// Subtract the union of exact tool disks centered at the given points. One
// chain implementation shared by the capsule (Task 2) and arc (Task 3) paths:
// callers only choose where the centers sit; exact predicates still decide
// point-in on the constructed disks, so the union under-covers the true sweep.
void Stock2::subtract_point_chain(const std::vector<std::pair<double, double>>& centers,
                                  double radius)
{
    // AGGREGATED union: build every disk, union them in ONE sweep-based pass, then a
    // single difference. The prior per-disk `region.join(disk_set)` was the O(N^2)
    // Boolean_set_operations anti-pattern -- each incremental join re-swept the whole
    // growing union. General_polygon_set_2::join(first, last) unions the range in one
    // arrangement. The result is identical (union is associative/commutative) -- this
    // is purely a cost fix, not a semantic change.
    std::vector<GpsPolygon> disks;
    disks.reserve(centers.size());
    for (const auto& [x, y] : centers) {
        disks.push_back(disk_polygon(EPoint(x, y), Epeck::FT(radius)));
    }
    Gps region;
    region.join(disks.begin(), disks.end());
    set_->difference(region);
}

void Stock2::subtract_capsule(double x0, double y0, double x1, double y1, double radius)
{
    if (radius <= 0.0) throw std::invalid_argument("radius should be positive.");
    const double dx = x1 - x0;
    const double dy = y1 - y0;
    const double len_sq = dx * dx + dy * dy;
    if (len_sq == 0.0) {
        subtract_disk(x0, y0, radius);
        return;
    }
    // Disk chain spacing s with per-disk radius r: the chain covers the true
    // capsule shrunk by delta = r - sqrt(r^2 - (s/2)^2) <= s^2/(4r). Choose s
    // so delta <= CHAIN_SLACK_FRACTION * r  =>  s = 2r*sqrt(fraction). s only
    // sets the (certified, under-covering) construction density.
    const double len = std::sqrt(len_sq);
    const double spacing = 2.0 * radius * std::sqrt(CHAIN_SLACK_FRACTION);
    const int n = std::max(1, static_cast<int>(std::ceil(len / spacing)));
    std::vector<std::pair<double, double>> centers;
    centers.reserve(n + 1);
    for (int i = 0; i <= n; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(n);
        centers.emplace_back(x0 + t * dx, y0 + t * dy);
    }
    subtract_point_chain(centers, radius);
}

// The same swept region as the chain above, under-covered the same way, but in
// six curves rather than in hundreds of disks.
//
// WHY A CHAIN WAS USED AT ALL. The capsule's two side lines stand off from the
// segment by r/sqrt(dx^2+dy^2) * (dy, -dx) -- an irrational offset that
// Gps_circle_segment_traits_2 cannot hold, since a line there is a rational
// triple (a, b, c). The chain sidesteps that by never building the sides.
//
// WHAT REPLACES IT. The side lines do not have to be MET, only UNDER-CUT: the
// chain is already an under-approximation, budgeted at CHAIN_SLACK_FRACTION * r.
// So take the EXACT perpendicular (-dy, dx) -- exact because each double is a
// rational and the arrangement holds rationals -- and scale it by a rational
// `scale` whose length h = scale * ||(-dy, dx)|| lands inside the admissible band
//
//     (1 - CHAIN_SLACK_FRACTION) * r  <=  h  <=  r.
//
// Then region = disk(A, r) U disk(B, r) U rect(A, B, h), where the rectangle has
// corners A -+ h_vec and B -+ h_vec, all four EXACT rational points.
//
// Both halves of the contract are exact consequences, not estimates:
//   SUBSET   a point of the rectangle is A + t*d + s*h_vec with t in [0,1],
//            |s| <= 1, and h_vec is perpendicular to d, so its distance to the
//            segment is exactly |s| * h <= h <= r. The end disks are the exact
//            caps. So the removed region is contained in the true capsule --
//            the safety direction: this never over-cuts.
//   COVERAGE anything within (1 - CHAIN_SLACK_FRACTION)*r of the segment is
//            either within that distance of an endpoint (inside an end disk,
//            which has the FULL radius r) or has a perpendicular foot in the
//            interior of the segment at offset <= (1 - f)*r <= h, hence inside
//            the rectangle.
//
// `scale` is a CONSTRUCTION parameter -- the quad twin of the chain's spacing --
// so computing it in doubles is legitimate; what is NOT legitimate is trusting
// it. It is aimed at the MIDDLE of the band (half-width (1 - f/2)*r), which
// leaves f/2 = 5e-5 of relative headroom on each side against a double rounding
// of ~2e-16: eleven orders of magnitude of margin. The band is then CHECKED, as
// an exact rational comparison of squares, and a failure throws rather than
// silently removing a region the certificate does not cover.
void Stock2::subtract_capsule_quad(double x0, double y0, double x1, double y1, double radius)
{
    if (radius <= 0.0) throw std::invalid_argument("radius should be positive.");
    // Finiteness is a DOUBLE concept, checked here at the boundary before
    // anything is injected: Epeck::FT(NaN) has no meaning to build on.
    if (!std::isfinite(x0) || !std::isfinite(y0)
        || !std::isfinite(x1) || !std::isfinite(y1) || !std::isfinite(radius)) {
        throw NonFiniteCapsuleInputError(
            "subtract_capsule_quad requires finite endpoints and radius; got "
            "x0=" + std::to_string(x0) + ", y0=" + std::to_string(y0)
            + ", x1=" + std::to_string(x1) + ", y1=" + std::to_string(y1)
            + ", radius=" + std::to_string(radius) + ".");
    }

    // Each double is injected as itself -- no snapping, no tolerance at the seam.
    const EPoint start(x0, y0);
    const EPoint end(x1, y1);
    // A zero-length motion sweeps the plain disk. Exact point equality decides
    // it, on the injected points; there is no short-segment threshold anywhere
    // below, because the rectangle stays an exact rectangle at every length.
    if (start == end) {
        subtract_disk(x0, y0, radius);
        return;
    }

    const Epeck::FT tool_radius(radius);
    const EVector along = end - start;
    const EVector normal = along.perpendicular(CGAL::COUNTERCLOCKWISE);

    // Construction parameter, not a decision: the rational scale that puts the
    // half-width in the middle of the band. Endpoint differences in double are
    // exact when the endpoints are close (Sterbenz) and relatively accurate
    // otherwise, and any inaccuracy at all is answered by the exact band check.
    const double length = std::hypot(x1 - x0, y1 - y0);
    const double scale_value = (1.0 - 0.5 * CHAIN_SLACK_FRACTION) * radius / length;
    if (!std::isfinite(scale_value) || scale_value <= 0.0) {
        throw CapsuleQuadCertificateError(
            "capsule quad half-width scale is not a positive finite number "
            "(radius=" + std::to_string(radius) + ", length=" + std::to_string(length)
            + "); the radius-to-length ratio is outside the representable range.");
    }

    const EVector half_width_vector = Epeck::FT(scale_value) * normal;
    const Epeck::FT half_width_squared = half_width_vector.squared_length();
    const Epeck::FT slack_half_width =
        tool_radius * (Epeck::FT(1) - Epeck::FT(CHAIN_SLACK_FRACTION));

    // The certificate. Two exact rational comparisons of SQUARED lengths -- no
    // square root is taken, so nothing here is inexact.
    if (CGAL::compare(half_width_squared, tool_radius * tool_radius) == CGAL::LARGER) {
        throw CapsuleQuadCertificateError(
            "capsule quad half-width exceeds the tool radius; the removed region "
            "would not be contained in the swept capsule.");
    }
    if (CGAL::compare(half_width_squared, slack_half_width * slack_half_width) == CGAL::SMALLER) {
        throw CapsuleQuadCertificateError(
            "capsule quad half-width falls below the documented slack budget "
            "(1 - CHAIN_SLACK_FRACTION) * radius.");
    }

    // Counterclockwise by construction: the first edge is `along` and the second
    // is 2 * half_width_vector, whose cross product is scale * ||along||^2 > 0.
    // Checked exactly all the same -- a silently clockwise or degenerate general
    // polygon would corrupt the boolean set with no error of its own.
    const std::vector<EPoint> corners = {
        start - half_width_vector,
        end - half_width_vector,
        end + half_width_vector,
        start + half_width_vector,
    };
    if (CGAL::orientation(corners[0], corners[1], corners[2]) != CGAL::LEFT_TURN) {
        throw CapsuleQuadCertificateError(
            "capsule quad corners are not counterclockwise; the rectangle is "
            "degenerate.");
    }

    // AGGREGATED union of the three parts, then ONE difference -- the same
    // single-sweep discipline subtract_point_chain follows.
    std::vector<GpsPolygon> parts;
    parts.reserve(3);
    parts.push_back(disk_polygon(start, tool_radius));
    parts.push_back(disk_polygon(end, tool_radius));
    parts.push_back(exact_linear_polygon(corners));
    Gps region;
    region.join(parts.begin(), parts.end());
    set_->difference(region);
}

void Stock2::subtract_arc_sweep(double cx, double cy, double sx, double sy,
                                double ex, double ey, bool cw, double tool_radius)
{
    if (tool_radius <= 0.0) throw std::invalid_argument("tool_radius should be positive.");
    const double rx = sx - cx, ry = sy - cy;
    const double guide_r = std::hypot(rx, ry);
    if (guide_r == 0.0) { subtract_disk(cx, cy, tool_radius); return; }

    // A FULL turn sweeps the exact annulus between guide_r - tool_radius and
    // guide_r + tool_radius, so it needs no chain at all. This is a CORRECTNESS
    // improvement, not a loosening of a tolerance: the chain deliberately
    // UNDER-covers the true swept region, leaving the model with sagitta slivers
    // of material the tool had in fact removed -- material that raises every
    // later engagement reading. The annulus is the swept region exactly, so
    // nothing is left over and nothing has to be compensated for.
    //
    // Direction of travel does not enter: a full turn sweeps the same set either
    // way, so `cw` is irrelevant here (it still selects the arc for a partial
    // sweep below).
    //
    // guide_r is a double surrogate for the guide radius, which is irrational in
    // general (sqrt of a rational) and therefore NOT representable as a rational
    // squared radius. That surrogate is a pre-existing property of this
    // double-valued API -- the chain samples the very same approximate circle --
    // and it is injected exactly, with no snapping and no correction constant.
    if (sx == ex && sy == ey) {
        const auto [inner, outer] = full_turn_annulus_bounds(guide_r, tool_radius);
        subtract_annulus_exact(EPoint(cx, cy), inner, outer);
        return;
    }

    double a0 = std::atan2(ry, rx);
    double a1 = std::atan2(ey - cy, ex - cx);
    double sweep = cw ? a0 - a1 : a1 - a0;
    if (sweep <= 0.0) sweep += 2.0 * std::numbers::pi;

    // Chain spacing along the guide: disks of tool_radius at arc-length step s
    // under-cover the true sweep by delta <= s^2/(4*tool_radius) (chord
    // sagitta bound; the straight-chord bound is conservative for the arc chain
    // because consecutive centers are closer along the chord than the arc).
    const double spacing = 2.0 * tool_radius * std::sqrt(CHAIN_SLACK_FRACTION);
    const double arc_len = guide_r * sweep;
    const int n = std::max(4, static_cast<int>(std::ceil(arc_len / spacing)));

    std::vector<std::pair<double, double>> centers;
    centers.reserve(n + 1);
    for (int i = 0; i <= n; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(n);
        const double a = cw ? a0 - sweep * t : a0 + sweep * t;
        centers.emplace_back(cx + guide_r * std::cos(a), cy + guide_r * std::sin(a));
    }
    subtract_point_chain(centers, tool_radius);
}

DepletionTrace Stock2::subtract_exact_segment(
    const ExactSegmentMotion2& motion,
    const Epeck::FT& tool_radius,
    const Epeck::FT& max_chord,
    std::size_t center_count_limit)
{
    ExactDepletionConstruction2 construction = construct_exact_segment_depletion(
        motion,
        tool_radius,
        max_chord,
        center_count_limit);
    std::unique_ptr<Gps> trial = std::make_unique<Gps>(*set_);
    Gps removal = exact_disk_union(construction.centers, tool_radius);
    trial->difference(removal);
    validate_depletion_trace(construction.trace);
    set_.swap(trial);
    return std::move(construction.trace);
}

DepletionTrace Stock2::subtract_exact_full_circle(
    const ExactCircleMotion2& motion,
    const Epeck::FT& tool_radius,
    const Epeck::FT& max_chord,
    std::size_t center_count_limit)
{
    ExactDepletionConstruction2 construction = construct_exact_full_circle_depletion(
        motion,
        tool_radius,
        max_chord,
        center_count_limit);
    std::unique_ptr<Gps> trial = std::make_unique<Gps>(*set_);
    Gps removal = exact_disk_union(construction.centers, tool_radius);
    trial->difference(removal);
    validate_depletion_trace(construction.trace);
    set_.swap(trial);
    return std::move(construction.trace);
}

ExactArcDepletionTrace2 Stock2::subtract_exact_arc(
    const AuditArcMotion2& motion,
    const Epeck::FT& tool_radius,
    const Epeck::FT& max_chord,
    std::size_t center_count_limit)
{
    ExactArcDepletionConstruction2 construction = construct_exact_arc_depletion(
        motion,
        tool_radius,
        max_chord,
        center_count_limit);
    std::unique_ptr<Gps> trial = std::make_unique<Gps>(*set_);
    Gps removal = exact_disk_union(construction.centers, tool_radius);
    trial->difference(removal);
    if (!exact_arc_structural_density_holds(
            motion,
            max_chord,
            construction.trace.parameters())
        || !construction.trace.matches_exact_inputs(
            tool_radius,
            max_chord,
            center_count_limit)
        || !construction.trace.matches_motion(motion)) {
        throw ExactArcForgedTraceError(
            "exact arc depletion failed atomic trace validation");
    }
    set_.swap(trial);
    return std::move(construction.trace);
}

namespace {

struct ExactSweepOracle {
    Gps removal;
    Gps sweep;
};

ExactSweepOracle exact_segment_sweep_oracle(
    const ExactSegmentMotion2& motion,
    const Epeck::FT& exact_length,
    const Epeck::FT& tool_radius,
    const Epeck::FT& max_chord,
    std::size_t center_count_limit)
{
    if (CGAL::compare(exact_length, Epeck::FT(0)) != CGAL::LARGER
        || CGAL::compare(
               (motion.end - motion.start).squared_length(),
               exact_length * exact_length)
            != CGAL::EQUAL) {
        throw ExactDepletionConstructionError(
            "segment sweep oracle requires an exact positive Pythagorean length.");
    }
    ExactDepletionConstruction2 construction = construct_exact_segment_depletion(
        motion,
        tool_radius,
        max_chord,
        center_count_limit);
    Gps removal = exact_disk_union(construction.centers, tool_radius);

    const EVector direction = motion.end - motion.start;
    const EVector normal(
        -direction.y() * tool_radius / exact_length,
        direction.x() * tool_radius / exact_length);
    const std::vector<EPoint> rectangle{
        motion.start - normal,
        motion.end - normal,
        motion.end + normal,
        motion.start + normal,
    };
    std::vector<GpsPolygon> sweep_parts{
        exact_linear_polygon(rectangle),
        disk_polygon(motion.start, tool_radius),
        disk_polygon(motion.end, tool_radius),
    };
    Gps sweep;
    sweep.join(sweep_parts.begin(), sweep_parts.end());
    return {std::move(removal), std::move(sweep)};
}

ExactSweepOracle exact_full_circle_sweep_oracle(
    const ExactCircleMotion2& motion,
    const Epeck::FT& guide_radius,
    const Epeck::FT& tool_radius,
    const Epeck::FT& max_chord,
    std::size_t center_count_limit)
{
    if (CGAL::compare(guide_radius, Epeck::FT(0)) != CGAL::LARGER
        || CGAL::compare(
               motion.phase_vector.squared_length(),
               guide_radius * guide_radius)
            != CGAL::EQUAL) {
        throw ExactDepletionConstructionError(
            "circle sweep oracle requires an exact positive rational guide radius.");
    }
    ExactDepletionConstruction2 construction = construct_exact_full_circle_depletion(
        motion,
        tool_radius,
        max_chord,
        center_count_limit);
    Gps removal = exact_disk_union(construction.centers, tool_radius);

    // The true swept region of a full turn IS the annulus, built by the same
    // shared exact_annulus_region that Stock2::subtract_annulus_exact removes --
    // so the fast path removes precisely the region this oracle certifies.
    const Epeck::FT inner = (CGAL::compare(guide_radius, tool_radius) == CGAL::LARGER)
        ? Epeck::FT(guide_radius - tool_radius)
        : Epeck::FT(0);
    Gps sweep = exact_annulus_region(
        motion.center,
        inner,
        guide_radius + tool_radius);
    return {std::move(removal), std::move(sweep)};
}

bool exact_induction_holds(
    const Stock2& initial,
    const ExactSweepOracle& oracle)
{
    Gps true_remaining(initial.set());
    true_remaining.difference(oracle.sweep);
    Gps modeled_remaining(initial.set());
    modeled_remaining.difference(oracle.removal);
    return exact_set_is_subset(true_remaining, modeled_remaining);
}

} // namespace

bool exact_segment_undercover_holds(
    const ExactSegmentMotion2& motion,
    const Epeck::FT& exact_length,
    const Epeck::FT& tool_radius,
    const Epeck::FT& max_chord,
    std::size_t center_count_limit)
{
    const ExactSweepOracle oracle = exact_segment_sweep_oracle(
        motion,
        exact_length,
        tool_radius,
        max_chord,
        center_count_limit);
    return exact_set_is_subset(oracle.removal, oracle.sweep);
}

bool exact_full_circle_undercover_holds(
    const ExactCircleMotion2& motion,
    const Epeck::FT& guide_radius,
    const Epeck::FT& tool_radius,
    const Epeck::FT& max_chord,
    std::size_t center_count_limit)
{
    const ExactSweepOracle oracle = exact_full_circle_sweep_oracle(
        motion,
        guide_radius,
        tool_radius,
        max_chord,
        center_count_limit);
    return exact_set_is_subset(oracle.removal, oracle.sweep);
}

bool exact_segment_induction_holds(
    const Stock2& initial,
    const ExactSegmentMotion2& motion,
    const Epeck::FT& exact_length,
    const Epeck::FT& tool_radius,
    const Epeck::FT& max_chord,
    std::size_t center_count_limit)
{
    const ExactSweepOracle oracle = exact_segment_sweep_oracle(
        motion,
        exact_length,
        tool_radius,
        max_chord,
        center_count_limit);
    return exact_induction_holds(initial, oracle);
}

bool exact_full_circle_induction_holds(
    const Stock2& initial,
    const ExactCircleMotion2& motion,
    const Epeck::FT& guide_radius,
    const Epeck::FT& tool_radius,
    const Epeck::FT& max_chord,
    std::size_t center_count_limit)
{
    const ExactSweepOracle oracle = exact_full_circle_sweep_oracle(
        motion,
        guide_radius,
        tool_radius,
        max_chord,
        center_count_limit);
    return exact_induction_holds(initial, oracle);
}

// --- Instrumentation --------------------------------------------------------
// Read-only probes on the arrangement the boolean engine maintains. They exist to
// measure how the exact representation GROWS across a run of subtractions; nothing
// they return feeds a geometric decision, so exact-kernel discipline is untouched.

Stock2::ArrangementStats Stock2::arrangement_stats() const
{
    // General_polygon_set_2 exposes a CONST arrangement accessor (vendored CGAL
    // 6.0.1, General_polygon_set_2.h:105 declares `const Arrangement_2&
    // arrangement() const` alongside the mutable overload), so unlike
    // engagement_2.cpp::engaged_arcs_zone -- which needs a non-const handle for
    // Arrangement_zone_2 and therefore documents a read-only const_cast -- this
    // probe needs no cast at all. These are plain counters: no exact evaluation is
    // triggered, so reading them does not perturb a timed run.
    const Gps::Arrangement_2& arr = set_->arrangement();
    return { arr.number_of_vertices(), arr.number_of_halfedges(), arr.number_of_faces() };
}

bool Stock2::representation_is_valid() const
{
    // Gps_on_surface_base_2::is_valid is non-const because it hands the
    // arrangement to CGAL::is_valid by mutable reference; the set itself is only
    // read. unique_ptr does not propagate constness to its pointee, so this
    // const probe needs no cast.
    return set_->is_valid();
}

namespace {

// Printed decimal length of one exact rational, measured by streaming it.
// Backend-agnostic on purpose: the repo rule is to use kernel/number-type
// abstractions rather than naming a concrete rational backend (this build is
// CGAL_DISABLE_GMP + CGAL_USE_BOOST_MP), and operator<< is guaranteed by the
// number type's concept. The printed length -- numerator, '/', denominator and
// any sign -- is a proxy for bit length (bits ~ digits * log2(10)); only its
// GROWTH across successive operations is ever interpreted, never its absolute
// magnitude.
std::size_t decimal_digits(const Epeck::FT& value)
{
    std::ostringstream stream;
    stream << value.exact();
    return stream.str().size();
}

// Running max/mean over the printed lengths of the exact parts sampled so far.
struct DigitAccumulator {
    std::size_t max_digits = 0;
    double sum = 0.0;
    std::size_t count = 0;

    void add_rational(const Epeck::FT& part)
    {
        const std::size_t digits = decimal_digits(part);
        max_digits = std::max(max_digits, digits);
        sum += static_cast<double>(digits);
        ++count;
    }

    // a1() and root() are ONLY defined on an EXTENDED Sqrt_extension -- a rational
    // coordinate carries a0() alone, and reading the extension parts
    // unconditionally is undefined behaviour. Same guard as
    // engagement_2.cpp::as_radpoint and exact_stock_region_2.cpp::lift_coordinate.
    void add_coordinate(const GpsPoint::CoordNT& coordinate)
    {
        add_rational(coordinate.a0());
        if (!coordinate.is_extended()) {
            return;
        }
        add_rational(coordinate.a1());
        add_rational(coordinate.root());
    }
};

} // namespace

Stock2::CoordinateDigits Stock2::coordinate_digits() const
{
    const Gps::Arrangement_2& arr = set_->arrangement();
    DigitAccumulator accumulator;
    for (auto vertex = arr.vertices_begin(); vertex != arr.vertices_end(); ++vertex) {
        accumulator.add_coordinate(vertex->point().x());
        accumulator.add_coordinate(vertex->point().y());
    }
    const double mean = (accumulator.count == 0)
        ? 0.0
        : accumulator.sum / static_cast<double>(accumulator.count);
    return { accumulator.max_digits, mean, accumulator.count };
}

NB_MODULE(_stock_2, m)
{
    register_audit_classification_2(m);
    nb::exception<ExactDepletionConstructionError> construction_error(
        m,
        "ExactDepletionConstructionError");
    nb::exception<ExactDepletionCenterLimitError>(
        m,
        "ExactDepletionCenterLimitError",
        construction_error.ptr());
    nb::exception<ExactArcDepletionPolicyError>(
        m,
        "ExactArcDepletionPolicyError",
        construction_error.ptr());
    nb::exception<ExactArcForgedTraceError>(
        m,
        "ExactArcForgedTraceError",
        construction_error.ptr());
    nb::exception<NonFiniteExactArcDepletionInputError>(
        m,
        "NonFiniteExactArcDepletionInputError",
        PyExc_ValueError);

    // Argument faults at the double boundary: ValueError-derived so they read the
    // same way as every other malformed-input rejection in this module.
    nb::exception<InvalidAnnulusRadiiError>(m, "InvalidAnnulusRadiiError", PyExc_ValueError);
    nb::exception<NonFiniteAnnulusInputError>(m, "NonFiniteAnnulusInputError", PyExc_ValueError);
    nb::exception<NonFiniteCapsuleInputError>(m, "NonFiniteCapsuleInputError", PyExc_ValueError);

    // Not an argument fault: the quad capsule's exact half-width certificate did
    // not hold. Registered explicitly under RuntimeError -- nanobind's default
    // base is Exception, which would leave "broken invariant" and "malformed
    // input" with no common ancestor a caller could tell apart.
    nb::exception<CapsuleQuadCertificateError>(m, "CapsuleQuadCertificateError", PyExc_RuntimeError);

    // Not an argument fault: a broken internal invariant of the local depletion
    // path, surfaced as a RuntimeError so it can never be confused with -- or
    // caught alongside -- a malformed-input rejection.
    nb::exception<LocalDepletionEscapedError>(m, "LocalDepletionEscapedError");

    nb::class_<DepletionTrace>(m, "DepletionTrace")
        .def_prop_ro("center_count", [](const DepletionTrace& trace) {
            return trace.center_count;
        })
        .def_prop_ro("center_parameters", [](const DepletionTrace& trace) {
            std::vector<std::tuple<int, std::size_t, std::size_t>> result;
            result.reserve(trace.center_parameters.size());
            for (const ExactCenterParameter2& parameter : trace.center_parameters) {
                result.emplace_back(
                    parameter.chart,
                    parameter.numerator,
                    parameter.denominator);
            }
            return result;
        })
        .def(
            "matches_exact_inputs",
            [](const DepletionTrace& trace,
               double expected_tool_radius,
               double expected_max_chord,
               std::size_t expected_center_count_limit) {
                return trace.matches_exact_inputs(
                    Epeck::FT(expected_tool_radius),
                    Epeck::FT(expected_max_chord),
                    expected_center_count_limit);
            },
            "expected_tool_radius"_a,
            "expected_max_chord"_a,
            "expected_center_count_limit"_a)
        .def_prop_ro("strategy_version", [](const DepletionTrace& trace) {
            return nb::bytes(
                trace.strategy_version.data(),
                trace.strategy_version.size());
        })
        .def_ro("cyclic", &DepletionTrace::cyclic)
        .def_ro("exact_incidence", &DepletionTrace::exact_incidence)
        .def_ro("exact_parameters_in_range", &DepletionTrace::exact_parameters_in_range)
        .def_ro("exact_anchors_present", &DepletionTrace::exact_anchors_present)
        .def_ro("exact_removal_radius_valid", &DepletionTrace::exact_removal_radius_valid)
        .def_ro("exact_chord_bound_holds", &DepletionTrace::exact_chord_bound_holds)
        .def_ro("exact_seam_chord_bound_holds", &DepletionTrace::exact_seam_chord_bound_holds);

    nb::class_<ExactArcDepletionTrace2>(m, "ExactArcDepletionTrace2")
        .def_prop_ro("center_count", [](const ExactArcDepletionTrace2& trace) {
            return trace.parameters().size();
        })
        .def(
            "matches_exact_inputs",
            [](const ExactArcDepletionTrace2& trace,
               double tool_radius,
               double max_chord,
               std::int64_t center_count_limit) {
                if (!std::isfinite(tool_radius)
                    || !std::isfinite(max_chord)) {
                    throw NonFiniteExactArcDepletionInputError(
                        "exact arc matcher inputs must be finite");
                }
                if (center_count_limit <= 0) {
                    throw ExactDepletionCenterLimitError(
                        "exact arc center-count limit must be positive");
                }
                return trace.matches_exact_inputs(
                    Epeck::FT(tool_radius),
                    Epeck::FT(max_chord),
                    static_cast<std::size_t>(center_count_limit));
            },
            "tool_radius"_a,
            "max_chord"_a,
            "center_count_limit"_a)
        .def(
            "matches_motion",
            &ExactArcDepletionTrace2::matches_motion,
            "motion"_a)
        .def_prop_ro("canonical_bytes", [](const ExactArcDepletionTrace2& trace) {
            return nb::bytes(
                trace.canonical_bytes().data(),
                trace.canonical_bytes().size());
        })
        .def_prop_ro("digest", [](const ExactArcDepletionTrace2& trace) {
            return nb::bytes(
                trace.digest().bytes().data(),
                trace.digest().bytes().size());
        })
        .def_prop_ro("strategy_version", [](const ExactArcDepletionTrace2& trace) {
            return nb::bytes(
                trace.strategy_version().data(),
                trace.strategy_version().size());
        })
        .def_prop_ro("cyclic", &ExactArcDepletionTrace2::cyclic);

    nb::class_<Stock2>(m, "Stock2")
        .def(nb::init<Eigen::Ref<const compas::RowMatrixXd>,
                      const std::vector<compas::RowMatrixXd>&>(),
             "boundary"_a, "holes"_a)
        .def("contains", &Stock2::contains, "x"_a, "y"_a)
        .def("is_empty", &Stock2::is_empty)
        .def("clone", &Stock2::clone)
        .def("is_subset_of", &Stock2::is_subset_of, "other"_a)
        .def("exactly_equals", &Stock2::exactly_equals, "other"_a)
        .def("subtract_capsule", &Stock2::subtract_capsule,
             "x0"_a, "y0"_a, "x1"_a, "y1"_a, "radius"_a)
        .def("subtract_capsule_quad", &Stock2::subtract_capsule_quad,
             "x0"_a, "y0"_a, "x1"_a, "y1"_a, "radius"_a)
        .def("subtract_arc_sweep", &Stock2::subtract_arc_sweep,
             "cx"_a, "cy"_a, "sx"_a, "sy"_a, "ex"_a, "ey"_a, "cw"_a, "tool_radius"_a)
        .def(
            "subtract_exact_segment",
            [](Stock2& stock,
               double x0,
               double y0,
               double x1,
               double y1,
               double tool_radius,
               double max_chord,
               std::size_t center_count_limit) {
                return stock.subtract_exact_segment(
                    {EPoint(x0, y0), EPoint(x1, y1)},
                    Epeck::FT(tool_radius),
                    Epeck::FT(max_chord),
                    center_count_limit);
            },
            "x0"_a,
            "y0"_a,
            "x1"_a,
            "y1"_a,
            "tool_radius"_a,
            "max_chord"_a,
            "center_count_limit"_a)
        .def(
            "subtract_exact_full_circle",
            [](Stock2& stock,
               double cx,
               double cy,
               double phase_x,
               double phase_y,
               bool clockwise,
               double tool_radius,
               double max_chord,
               std::size_t center_count_limit) {
                return stock.subtract_exact_full_circle(
                    {EPoint(cx, cy), EVector(phase_x, phase_y), clockwise},
                    Epeck::FT(tool_radius),
                    Epeck::FT(max_chord),
                    center_count_limit);
            },
            "cx"_a,
            "cy"_a,
            "phase_x"_a,
            "phase_y"_a,
            "clockwise"_a,
            "tool_radius"_a,
            "max_chord"_a,
            "center_count_limit"_a)
        .def(
            "subtract_exact_arc",
            [](Stock2& stock,
               const AuditArcMotion2& motion,
               double tool_radius,
               double max_chord,
               std::int64_t center_count_limit) {
                if (!std::isfinite(tool_radius)
                    || !std::isfinite(max_chord)) {
                    throw NonFiniteExactArcDepletionInputError(
                        "exact arc tool radius and chord bound must be finite");
                }
                if (center_count_limit <= 0) {
                    throw ExactDepletionCenterLimitError(
                        "exact arc center-count limit must be positive");
                }
                return stock.subtract_exact_arc(
                    motion,
                    Epeck::FT(tool_radius),
                    Epeck::FT(max_chord),
                    static_cast<std::size_t>(center_count_limit));
            },
            "motion"_a,
            "tool_radius"_a,
            "max_chord"_a,
            "center_count_limit"_a)
        .def("subtract_disk", &Stock2::subtract_disk, "cx"_a, "cy"_a, "radius"_a)
        .def("subtract_annulus", &Stock2::subtract_annulus,
             "cx"_a, "cy"_a, "inner_radius"_a, "outer_radius"_a)
        .def("subtract_disk_local", &Stock2::subtract_disk_local,
             "cx"_a, "cy"_a, "radius"_a)
        .def("subtract_annulus_local", &Stock2::subtract_annulus_local,
             "cx"_a, "cy"_a, "inner_radius"_a, "outer_radius"_a)
        .def("subtract_arc_sweep_local", &Stock2::subtract_arc_sweep_local,
             "cx"_a, "cy"_a, "sx"_a, "sy"_a, "ex"_a, "ey"_a, "cw"_a, "tool_radius"_a)
        .def("representation_is_valid", &Stock2::representation_is_valid)
        .def(
            "arrangement_stats",
            [](const Stock2& stock) {
                const Stock2::ArrangementStats stats = stock.arrangement_stats();
                return std::make_tuple(stats.vertices, stats.halfedges, stats.faces);
            })
        .def(
            "coordinate_digits",
            [](const Stock2& stock) {
                const Stock2::CoordinateDigits digits = stock.coordinate_digits();
                return std::make_tuple(digits.max_digits, digits.mean_digits, digits.sampled);
            });

    m.def(
        "exact_segment_point_is_incident",
        [](double x0, double y0, double x1, double y1, double px, double py) {
            return exact_segment_point_is_incident(
                {EPoint(x0, y0), EPoint(x1, y1)},
                EPoint(px, py));
        },
        "x0"_a,
        "y0"_a,
        "x1"_a,
        "y1"_a,
        "px"_a,
        "py"_a);
    m.def(
        "exact_circle_point_is_incident",
        [](double cx,
           double cy,
           double phase_x,
           double phase_y,
           double px,
           double py) {
            return exact_circle_point_is_incident(
                {EPoint(cx, cy), EVector(phase_x, phase_y), false},
                EPoint(px, py));
        },
        "cx"_a,
        "cy"_a,
        "phase_x"_a,
        "phase_y"_a,
        "px"_a,
        "py"_a);
    m.def(
        "exact_segment_structural_density_holds",
        [](double x0,
           double y0,
           double x1,
           double y1,
           double max_chord,
           const std::vector<std::tuple<int, std::size_t, std::size_t>>& parameters) {
            return exact_segment_structural_density_holds(
                {EPoint(x0, y0), EPoint(x1, y1)},
                Epeck::FT(max_chord),
                to_exact_center_parameters(parameters));
        },
        "x0"_a,
        "y0"_a,
        "x1"_a,
        "y1"_a,
        "max_chord"_a,
        "parameters"_a);
    m.def(
        "exact_full_circle_structural_density_holds",
        [](double cx,
           double cy,
           double phase_x,
           double phase_y,
           bool clockwise,
           double max_chord,
           const std::vector<std::tuple<int, std::size_t, std::size_t>>& parameters) {
            return exact_full_circle_structural_density_holds(
                {EPoint(cx, cy), EVector(phase_x, phase_y), clockwise},
                Epeck::FT(max_chord),
                to_exact_center_parameters(parameters));
        },
        "cx"_a,
        "cy"_a,
        "phase_x"_a,
        "phase_y"_a,
        "clockwise"_a,
        "max_chord"_a,
        "parameters"_a);
    m.def(
        "exact_segment_undercover_holds",
        [](double x0,
           double y0,
           double x1,
           double y1,
           double exact_length,
           double tool_radius,
           double max_chord,
           std::size_t center_count_limit) {
            return exact_segment_undercover_holds(
                {EPoint(x0, y0), EPoint(x1, y1)},
                Epeck::FT(exact_length),
                Epeck::FT(tool_radius),
                Epeck::FT(max_chord),
                center_count_limit);
        },
        "x0"_a,
        "y0"_a,
        "x1"_a,
        "y1"_a,
        "exact_length"_a,
        "tool_radius"_a,
        "max_chord"_a,
        "center_count_limit"_a);
    m.def(
        "exact_full_circle_undercover_holds",
        [](double cx,
           double cy,
           double phase_x,
           double phase_y,
           double guide_radius,
           double tool_radius,
           double max_chord,
           std::size_t center_count_limit) {
            return exact_full_circle_undercover_holds(
                {EPoint(cx, cy), EVector(phase_x, phase_y), false},
                Epeck::FT(guide_radius),
                Epeck::FT(tool_radius),
                Epeck::FT(max_chord),
                center_count_limit);
        },
        "cx"_a,
        "cy"_a,
        "phase_x"_a,
        "phase_y"_a,
        "guide_radius"_a,
        "tool_radius"_a,
        "max_chord"_a,
        "center_count_limit"_a);
    m.def(
        "exact_segment_induction_holds",
        [](const Stock2& initial,
           double x0,
           double y0,
           double x1,
           double y1,
           double exact_length,
           double tool_radius,
           double max_chord,
           std::size_t center_count_limit) {
            return exact_segment_induction_holds(
                initial,
                {EPoint(x0, y0), EPoint(x1, y1)},
                Epeck::FT(exact_length),
                Epeck::FT(tool_radius),
                Epeck::FT(max_chord),
                center_count_limit);
        },
        "initial"_a,
        "x0"_a,
        "y0"_a,
        "x1"_a,
        "y1"_a,
        "exact_length"_a,
        "tool_radius"_a,
        "max_chord"_a,
        "center_count_limit"_a);
    m.def(
        "exact_full_circle_induction_holds",
        [](const Stock2& initial,
           double cx,
           double cy,
           double phase_x,
           double phase_y,
           double guide_radius,
           double tool_radius,
           double max_chord,
           std::size_t center_count_limit) {
            return exact_full_circle_induction_holds(
                initial,
                {EPoint(cx, cy), EVector(phase_x, phase_y), false},
                Epeck::FT(guide_radius),
                Epeck::FT(tool_radius),
                Epeck::FT(max_chord),
                center_count_limit);
        },
        "initial"_a,
        "cx"_a,
        "cy"_a,
        "phase_x"_a,
        "phase_y"_a,
        "guide_radius"_a,
        "tool_radius"_a,
        "max_chord"_a,
        "center_count_limit"_a);
    m.def("exact_depletion_strategy_version", []() {
        const std::string& version = exact_depletion_strategy_version();
        return nb::bytes(version.data(), version.size());
    });
    m.def("exact_entry_depletion_strategy_version", []() {
        static const std::string version = "exact-precleared-entry-disk-v1";
        return nb::bytes(version.data(), version.size());
    });

    register_engagement(m);
}
