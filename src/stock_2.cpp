#include "stock_2.h"
#include "engagement_2.h"
#include "exact_boundary.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

// The seam guards, shared verbatim with engagement_2.cpp (rationale in
// exact_boundary.h). Every public entry point below runs them BEFORE any exact
// construction and before any mutation of set_.
using exact_boundary::format_double;
using exact_boundary::require_finite;
using exact_boundary::require_positive_radius;

namespace {

// Guard both coordinates of one polygon vertex. The parameter NAME identifies the
// offending coordinate ("boundary vertex 3 x", "hole 0 vertex 1 y"), and is built
// ONLY when a guard fires: the success path pays two isfinite tests and allocates
// nothing. The `prefix + axis` temporaries live to the end of the full-expression,
// so the pointer handed to require_finite stays valid for that call.
void require_finite_vertex(double x, double y, const char* role, int index)
{
    if (std::isfinite(x) && std::isfinite(y)) return;
    const std::string prefix = std::string(role) + " vertex " + std::to_string(index) + " ";
    require_finite(x, (prefix + "x").c_str());
    require_finite(y, (prefix + "y").c_str());
}

// Convert an Nx3 double matrix (rationals by construction) to a linear
// circle-segment general polygon, validating simplicity via Polygon_2. `role`
// names the ring in guard messages ("boundary", "hole 0"), since both the outer
// boundary and every hole arrive here. A non-finite vertex must be refused BEFORE
// is_simple(): +/-Inf coordinates passed that test and produced a usable object
// whose exact state was meaningless, while NaN ones were misdiagnosed as a
// self-intersection ("Polygon boundary must be simple").
GpsPolygon data_to_gps_polygon(Eigen::Ref<const compas::RowMatrixXd> vertices, const char* role)
{
    if (vertices.rows() < 3) {
        throw std::invalid_argument("Expected at least three polygon vertices.");
    }
    CGAL::Polygon_2<Epeck> simple_check;
    for (int i = 0; i < vertices.rows(); ++i) {
        const double x = vertices(i, 0);
        const double y = vertices(i, 1);
        require_finite_vertex(x, y, role, i);
        simple_check.push_back(EPoint(x, y));
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
{
    set_.insert(data_to_gps_polygon(boundary, "boundary"));
    for (std::size_t i = 0; i < holes.size(); ++i) {
        const std::string role = "hole " + std::to_string(i);
        Gps hole_set;
        hole_set.insert(data_to_gps_polygon(holes[i], role.c_str()));
        set_.difference(hole_set);
    }
}

bool Stock2::contains(double x, double y) const
{
    // The query point is injected exactly (FT(x), FT(y)) on the next line, so it
    // must be a real point first. This entry point had no validation of any kind.
    require_finite(x, "x");
    require_finite(y, "y");
    return set_.oriented_side(GpsPoint(Epeck::FT(x), Epeck::FT(y))) == CGAL::ON_POSITIVE_SIDE;
}

bool Stock2::is_empty() const
{
    return set_.is_empty();
}

// Fraction of the tool radius allowed as chain under-coverage slack; the
// generator's radial_clearance margin (1e-3 * D = 2e-3 * r) dominates it.
constexpr double CHAIN_SLACK_FRACTION = 1e-4;

namespace {

// Upper bound on the disk-chain interval count.
//
// WHY A COUNT NEEDS A BOUND AT ALL. `n` is an allocation size in machine terms,
// not a quality knob: the chain is materialised as n+1 centre pairs (16 B each)
// and then as n+1 exact two-arc disk polygons, all live at once for the
// aggregated join. `static_cast<int>` of a quotient above INT_MAX does not raise,
// it SATURATES to 2147483647 (measured on this toolchain, ARM64 fcvtzs), so
// `centers.reserve(n + 1)` then asks for 34.4 GB and the fill loop runs 2^31
// times until the OS kills the process.
//
// FINITENESS DOES NOT CLOSE THIS DOOR, which is why the bound exists on top of
// the seam guards: measured, a fully FINITE 1e6-long move with r = 0.02 gives
// ceil(len/spacing) = 2.5e9 and saturates exactly as +/-Inf does.
//
// THE VALUE. spacing = 2*r*sqrt(CHAIN_SLACK_FRACTION) = 0.02*r, so n = 50*len/r
// and this bound admits 2e5 tool radii of travel. For scale: the heaviest chain
// in this repo's suite is n = 600 (1571 for a full arc sweep), and a demanding
// real move -- 1 m of travel with a 1 mm tool radius -- is n = 50000, still 200x
// under the bound. At the bound itself the centres alone are 160 MB and the disk
// polygons above them are orders of magnitude more, so this is a RUNAWAY
// DETECTOR, not a machining policy: it cannot refuse a toolpath anyone meant.
//
// Being well below INT_MAX also makes the cast TOTAL rather than merely survivable:
// every value that passes the bound converts exactly (no saturation, no undefined
// conversion) and `n + 1` cannot overflow.
constexpr double MAX_CHAIN_INTERVALS = 1e7;

// Chain interval count for a sweep of `path_length` at `spacing`, validated
// BEFORE the otherwise-saturating cast to int. `what` names the quantity in the
// refusal so the caller can act on it.
//
// REFUSED, NEVER CLAMPED. Clamping to MAX_CHAIN_INTERVALS would widen the
// spacing s, and the chain's certified under-coverage delta <= s^2/(4r) =
// CHAIN_SLACK_FRACTION * r is exactly the property the generator's
// radial_clearance margin is sized against. Silently exceeding it would leave the
// stock model believing material remains where the tool has already cut: a wrong
// answer in place of a refusal, which is the one outcome this module may not
// produce.
//
// Spelled `!(intervals <= MAX)` so a NaN quotient is refused too, keeping the
// helper total independently of the finiteness guards its callers already ran.
int chain_intervals(double path_length, double spacing, int minimum, const char* what)
{
    const double intervals = std::ceil(path_length / spacing);
    if (!(intervals <= MAX_CHAIN_INTERVALS)) {
        throw std::invalid_argument(
            std::string(what) + " / chain spacing = " + format_double(intervals)
            + " disk-chain intervals, above the limit of " + format_double(MAX_CHAIN_INTERVALS)
            + " (shorten the motion or use a larger tool radius).");
    }
    return std::max(minimum, static_cast<int>(intervals));
}

} // namespace

void Stock2::subtract_disk(double cx, double cy, double radius)
{
    // Guarded BEFORE the difference below: a guard that fires after the model has
    // been mutated is not a guard. The old `radius <= 0.0` test let +Inf through
    // (Inf <= 0 is false), and the resulting cutter disk silently emptied the
    // ENTIRE stock -- while subtract_capsule and subtract_arc_sweep raised on the
    // same value, so three siblings disagreed about one input.
    require_finite(cx, "cx");
    require_finite(cy, "cy");
    require_positive_radius(radius, "radius");
    Gps region;
    region.insert(disk_polygon(EPoint(cx, cy), Epeck::FT(radius)));
    set_.difference(region);
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
    set_.difference(region);
}

void Stock2::subtract_capsule(double x0, double y0, double x1, double y1, double radius)
{
    // Finiteness here is not only about exact injection: the chain SIZING below
    // divides the segment length by the spacing, and the cast of that quotient to
    // int saturates rather than raising. Measured, the split is by VALUE CLASS and
    // is symmetric in the four coordinates: a NaN endpoint makes the quotient NaN
    // and the cast 0, so n = 1 and the FT injection raises promptly; a +/-Inf
    // endpoint makes the quotient infinite and the cast INT_MAX, so the chain fill
    // asks for 34.4 GB and the process is SIGKILLed (the pytest process died with
    // exit 137). Finiteness alone does NOT close the second door -- a large FINITE
    // ratio saturates identically -- which is what chain_intervals is for.
    require_finite(x0, "x0");
    require_finite(y0, "y0");
    require_finite(x1, "x1");
    require_finite(y1, "y1");
    require_positive_radius(radius, "radius");
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
    const int n = chain_intervals(len, spacing, 1, "capsule length");
    std::vector<std::pair<double, double>> centers;
    centers.reserve(n + 1);
    for (int i = 0; i <= n; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(n);
        centers.emplace_back(x0 + t * dx, y0 + t * dy);
    }
    subtract_point_chain(centers, radius);
}

void Stock2::subtract_arc_sweep(double cx, double cy, double sx, double sy,
                                double ex, double ey, bool cw, double tool_radius)
{
    // Same contract as subtract_capsule, over the six arc coordinates: exact
    // injection downstream, and a real guide circle. Finiteness does NOT bound the
    // chain count -- chain_intervals does that, separately, below.
    require_finite(cx, "cx");
    require_finite(cy, "cy");
    require_finite(sx, "sx");
    require_finite(sy, "sy");
    require_finite(ex, "ex");
    require_finite(ey, "ey");
    require_positive_radius(tool_radius, "tool_radius");
    const double rx = sx - cx, ry = sy - cy;
    const double guide_r = std::hypot(rx, ry);
    if (guide_r == 0.0) { subtract_disk(cx, cy, tool_radius); return; }

    double a0 = std::atan2(ry, rx);
    double a1 = std::atan2(ey - cy, ex - cx);
    double sweep = cw ? a0 - a1 : a1 - a0;
    if (sweep <= 0.0) sweep += 2.0 * std::numbers::pi;
    const bool full = (sx == ex && sy == ey);
    if (full) sweep = 2.0 * std::numbers::pi;

    // Chain spacing along the guide: disks of tool_radius at arc-length step s
    // under-cover the true sweep by delta <= s^2/(4*tool_radius) (chord
    // sagitta bound; the straight-chord bound is conservative for the arc chain
    // because consecutive centers are closer along the chord than the arc).
    const double spacing = 2.0 * tool_radius * std::sqrt(CHAIN_SLACK_FRACTION);
    const double arc_len = guide_r * sweep;
    const int n = chain_intervals(arc_len, spacing, 4, "arc length");

    std::vector<std::pair<double, double>> centers;
    centers.reserve(n + 1);
    for (int i = 0; i <= n; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(n);
        const double a = cw ? a0 - sweep * t : a0 + sweep * t;
        centers.emplace_back(cx + guide_r * std::cos(a), cy + guide_r * std::sin(a));
    }
    subtract_point_chain(centers, tool_radius);
}

Stock2::ArrangementStats Stock2::arrangement_stats() const
{
    // General_polygon_set_2 exposes a CONST arrangement accessor (verified in the
    // vendored CGAL 6.0.1: General_polygon_set_2.h declares both
    // `const Arrangement_2& arrangement() const` and the mutable overload), so
    // unlike engagement_2.cpp::engaged_arcs_zone -- which needs a non-const handle
    // for Arrangement_zone_2 and therefore documents a read-only const_cast --
    // this probe needs no cast at all. These are plain counters: no exact
    // evaluation is triggered, so reading them does not perturb a timed run.
    const Gps::Arrangement_2& arr = set_.arrangement();
    return { arr.number_of_vertices(), arr.number_of_halfedges(), arr.number_of_faces() };
}

namespace {

// Decimal length of an exact rational, measured by streaming it. Backend-agnostic
// on purpose: the repo rule is to use kernel/number-type abstractions rather than
// naming a concrete rational backend (this build is CGAL_DISABLE_GMP +
// CGAL_USE_BOOST_MP), and operator<< is guaranteed by the number type's concept.
// The printed length is a proxy for bit length (bits ~ digits * log2(10)); only
// its GROWTH is interpreted, never its absolute magnitude.
std::size_t decimal_digits(const Epeck::FT& v)
{
    std::ostringstream os;
    os << v.exact();
    return os.str().size();
}

// a1() and root() are only defined on an EXTENDED Sqrt_extension -- a rational
// coordinate carries a0() alone. engagement_2.cpp::as_radpoint guards the same
// way; reading the extension parts unconditionally is undefined behaviour.
void accumulate(const GpsPoint::CoordNT& c, std::size_t& max_digits, double& sum, std::size_t& count)
{
    auto take = [&](const Epeck::FT& part) {
        const std::size_t d = decimal_digits(part);
        max_digits = std::max(max_digits, d);
        sum += static_cast<double>(d);
        ++count;
    };
    take(c.a0());
    if (c.is_extended()) {
        take(c.a1());
        take(c.root());
    }
}

} // namespace

Stock2::CoordinateDigits Stock2::coordinate_digits() const
{
    const Gps::Arrangement_2& arr = set_.arrangement();
    std::size_t max_digits = 0;
    double sum = 0.0;
    std::size_t count = 0;
    for (auto v = arr.vertices_begin(); v != arr.vertices_end(); ++v) {
        accumulate(v->point().x(), max_digits, sum, count);
        accumulate(v->point().y(), max_digits, sum, count);
    }
    const double mean = (count == 0) ? 0.0 : sum / static_cast<double>(count);
    return { max_digits, mean, count };
}

NB_MODULE(_stock_2, m)
{
    nb::class_<Stock2>(m, "Stock2")
        .def(nb::init<Eigen::Ref<const compas::RowMatrixXd>,
                      const std::vector<compas::RowMatrixXd>&>(),
             "boundary"_a, "holes"_a)
        .def("contains", &Stock2::contains, "x"_a, "y"_a)
        .def("is_empty", &Stock2::is_empty)
        .def("subtract_capsule", &Stock2::subtract_capsule,
             "x0"_a, "y0"_a, "x1"_a, "y1"_a, "radius"_a)
        .def("subtract_arc_sweep", &Stock2::subtract_arc_sweep,
             "cx"_a, "cy"_a, "sx"_a, "sy"_a, "ex"_a, "ey"_a, "cw"_a, "tool_radius"_a)
        .def("subtract_disk", &Stock2::subtract_disk, "cx"_a, "cy"_a, "radius"_a)
        .def("arrangement_stats",
             [](const Stock2& s) {
                 const Stock2::ArrangementStats a = s.arrangement_stats();
                 return std::make_tuple(a.vertices, a.halfedges, a.faces);
             })
        .def("coordinate_digits",
             [](const Stock2& s) {
                 const Stock2::CoordinateDigits d = s.coordinate_digits();
                 return std::make_tuple(d.max_digits, d.mean_digits, d.sampled);
             });

    register_engagement(m);
}
