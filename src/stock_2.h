#pragma once

#include "compas.h"
#include "exact_depletion_2.h"
#include "exact_motion_2.h"

// compas.h supplies nanobind/stl/bind_vector.h (opaque nb::bind_vector) but NOT
// the automatic std::vector<T> <-> Python-list type caster. The Stock2
// constructor takes `holes` as a plain Python list, so this module needs the
// caster in its own (NB_STATIC) translation unit.
#include <nanobind/stl/vector.h>

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <utility>

#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/Gps_circle_segment_traits_2.h>
#include <CGAL/Polygon_2.h>

// Exact constructions kernel: boolean set operations on the 2D stock model
// build new geometry (sweep boundaries, intersection vertices) whose exact
// coordinates must be representable, so Epeck (not the repo-default Epick) is
// mandatory here. See CLAUDE.md: exact predicates, no epsilon decisions.
typedef CGAL::Gps_circle_segment_traits_2<Epeck> GpsTraits;
typedef CGAL::General_polygon_set_2<GpsTraits> Gps;
typedef GpsTraits::Polygon_2 GpsPolygon;             // circle-segment general polygon
typedef GpsTraits::Polygon_with_holes_2 GpsPolygonWithHoles;
typedef GpsTraits::Point_2 GpsPoint;                 // one-root coordinates
typedef GpsTraits::X_monotone_curve_2 GpsXCurve;

// One named exception per failure mode of the annulus sweep, so a caller can
// tell a malformed radius pair from a non-finite coordinate without parsing
// message text. Both are argument faults at the double boundary and are exposed
// to Python as ValueError subclasses.
class InvalidAnnulusRadiiError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class NonFiniteAnnulusInputError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// Full disk of the given radius centred at `center`, as a two-arc CCW general
// polygon (split at its x-extreme vertical-tangency points). Shared by the
// stock-subtraction paths and the engagement query (which intersects the stock
// with this exact cutter disk), so it carries external linkage here rather than
// hiding in stock_2.cpp's anonymous namespace.
GpsPolygon disk_polygon(const EPoint& center, const Epeck::FT& radius);

// Exact 2D stock model: the remaining material as a general polygon set of
// linear (this task) and circular (later tasks) boundary arcs. Later tasks
// subtract tool sweeps and query engagement; this task owns init + point-in.
class Stock2 {
public:
    Stock2(Eigen::Ref<const compas::RowMatrixXd> boundary,
           const std::vector<compas::RowMatrixXd>& holes);
    Stock2(Stock2&&) noexcept = default;
    Stock2& operator=(Stock2&&) noexcept = default;
    Stock2(const Stock2&) = delete;
    Stock2& operator=(const Stock2&) = delete;

    bool contains(double x, double y) const;
    bool is_empty() const;
    Stock2 clone() const;
    bool is_subset_of(const Stock2& other) const;
    bool exactly_equals(const Stock2& other) const;

    // Remove the exact disk of the given radius centred at (cx, cy).
    void subtract_disk(double cx, double cy, double radius);

    // Remove the exact annulus between `inner_radius` and `outer_radius` about
    // (cx, cy) -- the swept region of a disk of radius r carried about a circular
    // guide of radius rho, with inner = rho - r and outer = rho + r. Both radii
    // are doubles, hence exact rationals, so their squares are exactly what
    // Gps_circle_segment_traits_2 needs: the region is TWO curves, not a chain.
    // `inner_radius == 0` is the degenerate case (a full disk) and is handled
    // structurally, not by a tolerance.
    void subtract_annulus(double cx, double cy, double inner_radius, double outer_radius);

    // Exact core of the above. Callers that already hold exact radii use this
    // rather than round-tripping them through doubles: rho + r and rho - r are
    // exact rationals, and rounding them to the nearest double would be a snap at
    // a seam where exactness is free -- and would put the swept region roughly an
    // ulp off the sweep oracle the depletion certificates are proved against.
    void subtract_annulus_exact(const EPoint& center,
                                const Epeck::FT& inner_radius,
                                const Epeck::FT& outer_radius);

    // --- Local depletion ------------------------------------------------------
    // Same point set, same canonical Gps representation, computed by editing the
    // arrangement around the removed region instead of overlaying the whole
    // stock (stock_local_2.h, docs/local_depletion.md). Added ALONGSIDE the
    // global path above, which stays the reference the equivalence test decides
    // against; neither implementation shares a line with the other.

    void subtract_disk_local(double cx, double cy, double radius);
    void subtract_annulus_local(double cx, double cy, double inner_radius, double outer_radius);
    void subtract_annulus_exact_local(const EPoint& center,
                                      const Epeck::FT& inner_radius,
                                      const Epeck::FT& outer_radius);

    // Full turn -> the exact annulus, removed locally. A PARTIAL arc is still a
    // disk chain and still goes through the global path: its removed region is
    // a many-arc union that the local update has not been proved out for.
    void subtract_arc_sweep_local(double cx, double cy, double sx, double sy,
                                  double ex, double ey, bool cw, double tool_radius);

    // Remove the tool sweep along segment (x0,y0)->(x1,y1) as a certified
    // under-covering disk chain (the exact oriented capsule has irrational
    // side lines, so it is not representable in the circle-segment traits).
    void subtract_capsule(double x0, double y0, double x1, double y1, double radius);

    // Remove the tool sweep along the circular guide arc from (sx,sy) to
    // (ex,ey) about (cx,cy) — cw selects the sweep direction, start == end
    // means the full circle — as a certified under-covering disk chain.
    void subtract_arc_sweep(double cx, double cy, double sx, double sy,
                            double ex, double ey, bool cw, double tool_radius);

    DepletionTrace subtract_exact_segment(
        const ExactSegmentMotion2& motion,
        const Epeck::FT& tool_radius,
        const Epeck::FT& max_chord,
        std::size_t center_count_limit);

    DepletionTrace subtract_exact_full_circle(
        const ExactCircleMotion2& motion,
        const Epeck::FT& tool_radius,
        const Epeck::FT& max_chord,
        std::size_t center_count_limit);

    const Gps& set() const { return *set_; }          // engagement kernel reads this
    Gps& set() { return *set_; }

    // --- Instrumentation (diagnostics only) ----------------------------------
    // Neither accessor participates in any decision: they report the SIZE and the
    // exact-rational COMPLEXITY of the arrangement the boolean engine has built,
    // so bit growth over a run of subtractions can be measured. No epsilon, no
    // tolerance, no feedback into geometry.

    // Feature counts of the underlying arrangement. Cheap: pure counters, no
    // exact evaluation, so this MAY be read inside a timed run.
    struct ArrangementStats {
        std::size_t vertices;
        std::size_t halfedges;
        std::size_t faces;
    };
    ArrangementStats arrangement_stats() const;

    // Does the underlying set still satisfy General_polygon_set_2's own
    // representation invariant -- every edge separating faces of DIFFERENT
    // containment, oriented with the contained side on its left? A structural
    // gate for the local depletion path, orthogonal to point-set equality:
    // exactly_equals can pass on a set whose arrangement is no longer canonical.
    bool representation_is_valid() const;

    // Printed decimal length of the exact coordinates carried by the
    // arrangement's vertices. WARNING: this calls .exact() on every sampled
    // coordinate, collapsing the lazy filter and changing subsequent timings --
    // it is a DIAGNOSTIC and must never be read inside a timed measurement.
    struct CoordinateDigits {
        std::size_t max_digits;
        double mean_digits;
        std::size_t sampled;
    };
    CoordinateDigits coordinate_digits() const;

private:
    explicit Stock2(std::unique_ptr<Gps> set) noexcept;

    // Subtract the union of exact tool disks of the given radius centred at the
    // listed points — the one chain implementation shared by capsule and arc.
    void subtract_point_chain(const std::vector<std::pair<double, double>>& centers,
                              double radius);

    std::unique_ptr<Gps> set_;
};
