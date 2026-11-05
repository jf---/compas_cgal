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

#include "audit_certification_2.h"
#include "compas_matrix.h"
#include "engagement_2.h"
#include "exact_motion_2.h"

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

class Stock2;

class InvalidAnnulusViewPolygonError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

class NonFiniteAnnulusViewInputError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

class InvalidAnnulusViewRadiiError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

class InvalidAnnulusViewRangeError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

// A rim-cell witness could not be constructed or lies on P's boundary. Both are
// breaches of the event construction, never input conditions.
class AnnulusViewWitnessError : public AuditCertificationError {
public:
    using AuditCertificationError::AuditCertificationError;
};

// The half-open emission range [first, end) of annuli a decision reads.
struct AnnulusViewRange2 {
    std::size_t first;
    std::size_t end;
};

// Whether one station's two classifications agree, decided exactly: the
// disposition, and each run's endpoints by exact point equality.
struct StationAgreement2 {
    bool disposition;
    bool runs;
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

    AnnulusView2(AnnulusView2&&) noexcept;
    AnnulusView2& operator=(AnnulusView2&&) noexcept;
    ~AnnulusView2();

    std::size_t size() const;

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
    explicit AnnulusView2(std::unique_ptr<const Index> index);
    std::unique_ptr<const Index> index_;
};
