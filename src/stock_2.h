#pragma once

#include "compas_matrix.h"
#include "exact_depletion_2.h"
#include "exact_motion_2.h"

#include <array>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include <CGAL/Arr_trapezoid_ric_point_location.h>
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
typedef CGAL::Arr_trapezoid_ric_point_location<Gps::Arrangement_2>
    GpsPointLocation;

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

// Same split for the quad capsule: a non-finite coordinate or radius is an
// argument fault at the double boundary (ValueError), whereas a quad whose exact
// half-width falls outside its certified band is a broken INTERNAL invariant of
// the construction (RuntimeError) and must never be caught alongside one.
class NonFiniteCapsuleInputError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class CapsuleQuadCertificateError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class InvalidCircleRemovalInputError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// A Gps was offered as a stock's storage whose arrangement reads geometry traits
// that no live Gps in the offered set owns. There is no repair at this point:
// the caller has to build the set on traits something keeps alive. This is a
// broken INTERNAL invariant of the construction, never an argument fault.
class StockTraitsUnownedError : public std::runtime_error {
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
    Stock2(Stock2&& other) noexcept;
    Stock2& operator=(Stock2&& other) noexcept;
    Stock2(const Stock2&) = delete;
    Stock2& operator=(const Stock2&) = delete;

    bool contains(double x, double y) const;
    bool is_empty() const;
    Stock2 clone() const;
    void swap(Stock2& other) noexcept;
    bool is_subset_of(const Stock2& other) const;
    bool exactly_equals(const Stock2& other) const;

    // True when the geometry traits this stock's ARRANGEMENT reads belong to an
    // object this stock keeps alive.
    //
    // CGAL's Gps copy constructor gives the copy a fresh Traits_2 of its own but
    // builds the copy's arrangement as Aos_2(*(ps.m_arr)), and
    // Arrangement_on_surface_2::assign propagates a BORROWED traits pointer
    // verbatim (m_geom_traits = arr.m_own_traits ? new Traits_adaptor_2
    // : arr.m_geom_traits). Every Gps builds its arrangement in borrow mode, so
    // a copied set reads the traits of whatever set it was copied from -- and
    // its own freshly allocated traits is never used. The first two-operand
    // boolean operation then rebuilds the arrangement on the copy's OWN traits,
    // which is why no single traits object can stand for a whole clone family:
    // the owner has to be resolved from the object graph at every point where
    // the arrangement changes.
    //
    // This accessor states the lifetime the object graph assumes; it decides
    // nothing and no geometry depends on it.
    bool arrangement_traits_are_owned_for_audit() const;

    // Sufficient deletion test on this remaining stock, which may conservatively
    // include material already cleared by an emitted path. Circle triples are
    // world-XY (center x, center y, guide radius) in mm. Connector samples must
    // be the preserved old route: this test authorizes removing a circle only.
    bool can_remove_circle(
        const std::array<double, 3>& previous,
        const std::array<double, 3>& current,
        const std::array<double, 3>& next,
        Eigen::Ref<const compas::RowMatrixXd> connector_samples,
        double tool_radius) const;

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

    // Full machining-circle sweep: add/subtract the separately injected guide
    // and cutter radii in the exact kernel, retaining an uncut centre if rho > r.
    void subtract_circle_sweep(double cx, double cy, double guide_radius, double tool_radius);

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

    // Same swept region, same UNDER-covering contract, in SIX curves instead of
    // a chain of hundreds: the two end disks (exact) unioned with the rectangle
    // along the segment at a slightly reduced half-width.
    //
    // The capsule's true side lines stand off at r*sqrt(dx^2+dy^2), which is
    // irrational and therefore not representable here -- but they do not have to
    // be MET, only under-cut. The rectangle is built from the EXACT perpendicular
    // (dx, dy) rotated a quarter turn, scaled by a rational chosen so its exact
    // length h satisfies (1 - CHAIN_SLACK_FRACTION)*r <= h <= r, which is checked
    // as an exact rational comparison, not assumed. So:
    //   * every removed point is within r of the segment (SUBSET of the true
    //     capsule -- the safety direction, never over-cutting), and
    //   * every point within (1 - CHAIN_SLACK_FRACTION)*r of the segment IS
    //     removed -- the same slack budget the disk chain documents.
    // Added ALONGSIDE subtract_capsule, which stays the reference.
    void subtract_capsule_quad(double x0, double y0, double x1, double y1, double radius);

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

    ExactArcDepletionTrace2 subtract_exact_arc(
        const AuditArcMotion2& motion,
        const Epeck::FT& tool_radius,
        const Epeck::FT& max_chord,
        std::size_t center_count_limit);

    const Gps& set() const { return *set_; }          // engagement kernel reads this
    // Any mutable access starts a new arrangement epoch. Detach the CGAL
    // observer before returning the set because the caller may rebuild it.
    Gps& set()
    {
        point_location_.reset();
        return *set_;
    }
    // Built once per read-only arrangement epoch, then reused by every station.
    GpsPointLocation& point_location() const;

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
    // Inputs validated by can_remove_circle before any local geometry changes.
    void intersect_circle_sweep(double cx, double cy,
                                double guide_radius, double tool_radius);
    Stock2(std::shared_ptr<const Gps> traits_owner, std::shared_ptr<Gps> set);

    // Subtract the union of exact tool disks of the given radius centred at the
    // listed points — the one chain implementation shared by capsule and arc.
    void subtract_point_chain(const std::vector<std::pair<double, double>>& centers,
                              double radius);
    void replace_set(std::shared_ptr<Gps> replacement);

    // Which live Gps owns the geometry traits `set`'s arrangement reads: `set`
    // itself when it is a root, otherwise the offered candidate whose traits
    // object matches. Offer every set `set` may have been seeded by copy from.
    // No match throws StockTraitsUnownedError rather than adopting a pointer
    // into memory nothing keeps alive.
    static std::shared_ptr<const Gps> resolve_traits_owner(
        const std::shared_ptr<Gps>& set,
        const std::vector<std::shared_ptr<const Gps>>& candidates);

    // The Gps that owns the geometry traits this stock's ARRANGEMENT reads --
    // this stock's own set when that set is a root.
    //
    // A shared traits OBJECT per clone family cannot express this, and that is
    // the shape this replaces. CGAL's Gps copy constructor gives the copy a
    // fresh Traits_2 of its own but builds the copy's arrangement as
    // Aos_2(*(ps.m_arr)), and Arrangement_on_surface_2::assign propagates the
    // BORROWED traits pointer verbatim -- so a fresh clone reads the traits of
    // the set it was copied from. Then the first two-operand boolean operation
    // rebuilds the arrangement on the copy's OWN traits
    // (_difference(const Aos_2&) does `new Aos_2(m_traits)`), and from that
    // point the family object is no longer the object being read. The owner
    // therefore has to be RESOLVED from the object graph wherever set_ changes:
    // the constructor, clone() and replace_set().
    //
    // Declared BEFORE set_ so it is destroyed AFTER the arrangement that reads
    // it: members are destroyed in reverse declaration order. Aliases set_ when
    // set_ is a root.
    std::shared_ptr<const Gps> traits_owner_;
    std::shared_ptr<Gps> set_;
    mutable std::unique_ptr<GpsPointLocation> point_location_;
};

bool exact_segment_undercover_holds(
    const ExactSegmentMotion2& motion,
    const Epeck::FT& exact_length,
    const Epeck::FT& tool_radius,
    const Epeck::FT& max_chord,
    std::size_t center_count_limit);

bool exact_full_circle_undercover_holds(
    const ExactCircleMotion2& motion,
    const Epeck::FT& guide_radius,
    const Epeck::FT& tool_radius,
    const Epeck::FT& max_chord,
    std::size_t center_count_limit);

bool exact_segment_induction_holds(
    const Stock2& initial,
    const ExactSegmentMotion2& motion,
    const Epeck::FT& exact_length,
    const Epeck::FT& tool_radius,
    const Epeck::FT& max_chord,
    std::size_t center_count_limit);

bool exact_full_circle_induction_holds(
    const Stock2& initial,
    const ExactCircleMotion2& motion,
    const Epeck::FT& guide_radius,
    const Epeck::FT& tool_radius,
    const Epeck::FT& max_chord,
    std::size_t center_count_limit);

// True when BOTH sets of the sweep oracle the four certificate entry points
// above are built on read geometry traits they own. The oracle is built exactly
// as those entry points build it and then observed; the answer decides nothing.
//
// The oracle carries the modeled removal and the true swept region side by
// side, so a member whose arrangement reads freed traits is a certificate
// computed against freed memory. See arrangement_traits_are_owned_for_audit.
bool exact_segment_sweep_oracle_traits_are_owned_for_audit(
    const ExactSegmentMotion2& motion,
    const Epeck::FT& exact_length,
    const Epeck::FT& tool_radius,
    const Epeck::FT& max_chord,
    std::size_t center_count_limit);

bool exact_full_circle_sweep_oracle_traits_are_owned_for_audit(
    const ExactCircleMotion2& motion,
    const Epeck::FT& guide_radius,
    const Epeck::FT& tool_radius,
    const Epeck::FT& max_chord,
    std::size_t center_count_limit);
