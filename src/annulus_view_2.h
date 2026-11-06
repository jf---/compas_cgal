#pragma once

// The arrangement-free rim engine: one station's engaged rim arcs in
//
//     int P \ union_{j in view} A_j
//
// where P is the PART polygon (the material domain, never the tool-centre offset)
// and A_j = {x : max(rho_j - r, 0) <= |x - c_j| <= rho_j + r} is the closed swept
// annulus of emitted circle j. The arrangement engine (Stock2 plus the exact zone
// walk) decides the same set; AnnulusView2 builds no arrangement.
//
// One station:
//   events  the rim's intersections with the candidate annuli's two circles and
//           with P's edges, as one-root points, plus the rim's two x-extreme
//           points. Each is compared, never combined, so no cross-root arithmetic.
//   order   CCW along the rim (exact_rim_order_2.h), duplicates merged by ==.
//   cells   each open cell is decided at a rational point strictly inside it: in
//           int P and outside every closed A_j, a rational predicate. A cell never
//           contains an event, so the witness is on no wall and on no annulus
//           boundary it crosses.
//   arcs    the material cells go to the unchanged clip, run and cap code.
//
// Every double enters once, at the binding, by exact injection. Doubles appear
// after that only as outward-rounded candidate boxes (a box miss proves a miss)
// and as a construction hint for a witness (the witness is then checked exactly).
//
// The model is exact (ruling M3, 2026-11-04): a sub-ulp crossing of the rim with a
// wall is real stock. Measure-zero contacts match the regularized Stock2: a contact
// splits nothing (runs with equal endpoints merge), and a rim on an annulus
// boundary is in the closed annulus.

#include "annulus_filter_2.h"
#include "annulus_store_2.h"
#include "audit_certification_2.h"
#include "compas_matrix.h"
#include "engagement_2.h"
#include "exact_motion_2.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

class Stock2;

class InvalidAnnulusViewPolygonError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

// append or truncate on a view that AnnulusView2::build fixed at construction.
class AnnulusViewNotAppendableError : public std::logic_error {
public:
    using std::logic_error::logic_error;
};

// A rim-cell witness could not be constructed or lies on P's boundary. Both are
// breaches of the event construction, never input conditions.
class AnnulusViewWitnessError : public AuditCertificationError {
public:
    using AuditCertificationError::AuditCertificationError;
};

// Whether one station's two classifications agree, decided exactly: the
// disposition, and each run's endpoints by exact point equality.
struct StationAgreement2 {
    bool disposition;
    bool runs;
};

// What the filtered decision did since the view opened or the counters were reset.
struct AnnulusFilterCounters2 {
    std::int64_t decisions = 0;
    std::int64_t exact_decisions = 0;      // inputs the filter does not take: wholly exact
    std::int64_t stations = 0;
    std::int64_t certified_last_one = 0;
    std::int64_t certified_last_sixteen = 0;
    std::int64_t certified_range = 0;
    std::int64_t exact_stations = 0;       // fell through to the exact path
    std::int64_t exact_station_zero = 0;   // of those, the entry station
    std::int64_t refused = 0;
};

// One station's certifying view and chord bound, off the decision path.
struct AnnulusFilterReport2 {
    int tier;
    double chord_squared_bound;
};

class AnnulusView2 {
public:
    // Raises:
    //   InvalidAnnulusViewPolygonError: fewer than 3 vertices or not simple.
    //   NonFiniteAnnulusViewInputError: a non-finite coordinate or radius.
    //   InvalidAnnulusViewRadiiError: a negative guide or non-positive tool radius.
    static AnnulusView2 build(
        Eigen::Ref<const compas::RowMatrixXd> polygon,
        Eigen::Ref<const compas::RowMatrixXd> annuli,
        double tool_radius);

    // An empty view for the emission stream: `append` adds the next circle's annulus,
    // `truncate` drops a tail, and candidates come from the store's grid. Raises as
    // `build` does on the polygon and the tool radius.
    static AnnulusView2 open(
        Eigen::Ref<const compas::RowMatrixXd> polygon,
        double tool_radius);

    AnnulusView2(AnnulusView2&&) noexcept;
    AnnulusView2& operator=(AnnulusView2&&) noexcept;
    ~AnnulusView2();

    std::size_t size() const;

    // Raises AnnulusViewNotAppendableError on a view from `build`, and as
    // AnnulusStore2::append and AnnulusStore2::truncate do.
    std::size_t append(double cx, double cy, double guide_radius);
    void truncate(std::size_t size);

    // The store's selection, off the decision path: annuli of `range` whose box meets
    // `box`, and the total over many boxes (a query's cost without a Python list).
    std::vector<std::uint32_t> candidates(
        AnnulusViewRange2 range, const AnnulusSearchBox2& box) const;
    std::size_t count_candidates(
        AnnulusViewRange2 range, Eigen::Ref<const compas::RowMatrixXd> boxes) const;
    std::size_t stale_references() const;
    std::size_t store_rebuilds() const;

    // One static station: the rim at (cx, cy) on the view, with no own prefix.
    EngagementSample station(
        AnnulusViewRange2 range, double cx, double cy,
        double cap_chord_ratio, double gap_close_ratio) const;

    // The same as engagement_circle_prefixes on a Stock2 of the same view.
    std::vector<EngagementSample> circle_prefixes(
        AnnulusViewRange2 range, double guide_cx, double guide_cy,
        double guide_radius, double entry_x, double entry_y, bool clockwise,
        const std::vector<double>& parameters, double cap_chord_ratio,
        double gap_close_ratio) const;

    // The same as engagement_circle_exceeded on a Stock2 of the same view.
    bool circle_exceeded(
        AnnulusViewRange2 range, double guide_cx, double guide_cy,
        double guide_radius, double entry_x, double entry_y, bool clockwise,
        const std::vector<double>& parameters, double cap_chord_ratio,
        double gap_close_ratio) const;

    // The exact decision on the nested views of `range`: last-1, last-16, then all of
    // it. The views are nested and a pass on a larger stock is a pass on a smaller, so
    // a pass at any tier is the verdict, and only a refusal on `range` refuses.
    bool circle_exceeded_tiered(
        AnnulusViewRange2 range, double guide_cx, double guide_cy,
        double guide_radius, double entry_x, double entry_y, bool clockwise,
        const std::vector<double>& parameters, double cap_chord_ratio,
        double gap_close_ratio) const;

    // The production decision: AnnulusFilter2 certifies each station on the nested views,
    // and the exact path decides only what no view certified. Equal to circle_exceeded
    // on `range`. Inputs the filter does not take (clockwise, a gap closure, a parameter
    // of 1, invalid values) go wholly to the exact path, which raises as it always does.
    bool circle_exceeded_filtered(
        AnnulusViewRange2 range, double guide_cx, double guide_cy,
        double guide_radius, double entry_x, double entry_y, bool clockwise,
        const std::vector<double>& parameters, double cap_chord_ratio,
        double gap_close_ratio) const;

    // The filter alone, off the decision path: per station its certifying tier, or 3.
    std::vector<AnnulusFilterReport2> circle_filter_report(
        AnnulusViewRange2 range, double guide_cx, double guide_cy,
        double guide_radius, double entry_x, double entry_y,
        const std::vector<double>& parameters, double cap_chord_ratio) const;

    const AnnulusFilterCounters2& filter_counters() const;
    void reset_filter_counters();

    // Test hook: the filter certifies nothing, so every station takes the exact path.
    void force_filter_fallthrough(bool forced);

    // Differential oracle hooks, off the production path: the engine's exact
    // classification against the arrangement engine's on `stock`.
    StationAgreement2 station_agrees_with_stock(
        const Stock2& stock, AnnulusViewRange2 range, double cx, double cy,
        double cap_chord_ratio, double gap_close_ratio) const;

    std::vector<StationAgreement2> circle_agrees_with_stock(
        const Stock2& stock, AnnulusViewRange2 range, double guide_cx,
        double guide_cy, double guide_radius, double entry_x, double entry_y,
        bool clockwise, const std::vector<double>& parameters,
        double cap_chord_ratio, double gap_close_ratio) const;

private:
    struct Index;
    static std::unique_ptr<Index> make_index(
        Eigen::Ref<const compas::RowMatrixXd> polygon, double tool_radius);
    explicit AnnulusView2(std::unique_ptr<Index> index);
    std::unique_ptr<Index> index_;
};
