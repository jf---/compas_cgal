#pragma once

// Target sides for the sweep-edge rule (docs/held_coverage_certificate.md,
// "Theorem"): an edge e on a swept boundary is exposed only if the target
// members cover its checked (outer) side, i.e. the face beside e is in T.
//
// T = Omega u (U_e R_e) u (U_v D_r(v)) with the edge rectangles R_e and vertex
// disks D_r(v) of the capsule decomposition. dT lies in the union of the flank
// lines of the R_e and the circles of the D_r(v): Omega's own edges lie inside
// their rectangles, and the end segments of R_e inside the closed vertex disks.
// Every swept piece is cut at those supports (target cuts), so along one of its
// edges T-membership of the checked side is constant, and is decided at the
// edge's exact witness w:
//   1. w in the interior of a member: the checked side is in T;
//   2. w in a member's closed set but no interior: w lies on that member's
//      boundary along a support it shares with the piece (a cut would otherwise
//      have split the edge -- checked, a breach raises), and the member lies on
//      the checked side or not (P6, rational normal signs);
//   3. otherwise the checked side is outside T.
// A VERTEX_DISK member (one-root centre, the exact Omega*) is decided where its
// sign is decidable at w and reported UNDECIDED otherwise; its boundary is not
// cut yet, so a certificate with vertex disks never passes.

#include "coverage_certificate_2/counters_2.h"
#include "coverage_certificate_2/lines_2.h"
#include "coverage_certificate_2/regions_2.h"
#include "coverage_certificate_2/target_2.h"

#include <cstdint>
#include <vector>

namespace compas_cgal::coverage_certificate {

/// The target's pieces and members, as the certificate built them.
struct TargetModel {
    const MachinableTarget* target;        ///< domain, tool radius, vertex disks
    const std::vector<RectRegion>* rects;  ///< target edge rectangles
    const LineTable* lines;
    const CircleRegionSet* set;            ///< holds the target vertex circles
};

/// Target members whose boxes overlap a piece's box.
struct TargetCandidates {
    std::vector<std::uint32_t> rects;    ///< indices into TargetModel::rects
    std::vector<std::uint32_t> circles;  ///< target vertex circle ids in the CircleRegionSet
};

/// The support of the piece being decided, for coincidence analysis.
struct CheckedSupport {
    bool is_line;
    std::uint32_t id;        ///< line id, or boundary-circle id
    CGAL::Sign check_side;   ///< line pieces: side of the line's normal that is checked
};

enum class TargetSide : std::uint8_t { Inside = 0, Outside = 1, Undecided = 2 };

/// Decide whether the face on the checked side of a swept piece at its exact
/// witness w lies in T.
///
/// Raises:
///     CoverageCertificateNativeError: w lies on a target boundary along a
///         support the piece does not share (an uncut target support).
[[nodiscard]] TargetSide decide_checked_side(const TargetModel& model, const TargetCandidates& candidates,
                                             const CheckedSupport& support, const ExactPoint& w, FilterMode mode,
                                             DecisionCounters& counters);

/// Exact closed membership of the rational point u in T. Complete for u inside
/// the box the candidates were gathered for.
[[nodiscard]] bool target_contains(const TargetModel& model, const TargetCandidates& candidates, const Point_2& u,
                                   FilterMode mode, DecisionCounters& counters);

}  // namespace compas_cgal::coverage_certificate
