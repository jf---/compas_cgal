#pragma once

// Exact witnesses on an exposed flank sub-edge (P3, P4, P7 on a line chart).
//
// P3: a rational chart parameter t strictly inside the sub-edge (proposal from
//     the EXACT endpoint towers' enclosures, accepted by exact tower signs, exact
//     bisection otherwise; line_sweep_2.h). The on-edge witness w = x(t) is a
//     one-root point in Q(sqrt(L)).
// P4: w is re-classified exactly against every counted candidate; it must be in
//     the exterior of all of them.
// P7: the uncovered point is the RATIONAL point u = b0 + t d + lambda n, with
//     lambda a rational within h = 2^-k of the flank's normal offset
//     alpha = (k1 / |n|^2) sqrt(L), on the checked side (exactly:
//     sign(lambda |n|^2 - k1 sqrt(L)) is the checked side). lambda comes from an
//     exact rational bracket of alpha refined by bisection; for a rational line
//     alpha = 0 and lambda = +/- h. u is accepted when it lies inside the piece's
//     (widened) box, outside every candidate sweep (owner and coincident ones
//     included) and satisfies the caller's predicate (u in T).

#include "coverage_certificate_2/counters_2.h"
#include "coverage_certificate_2/line_sweep_2.h"
#include "coverage_certificate_2/lines_2.h"
#include "coverage_certificate_2/regions_2.h"

#include <cstdint>
#include <functional>
#include <vector>

namespace compas_cgal::coverage_certificate {

/// What P4 and P7 test a flank witness against.
struct LineWitnessContext {
    const CircleRegionSet& set;
    const LineTable& lines;
    const std::vector<RectRegion>& rects;  ///< swept rectangles
    FilterMode mode;
};

/// P3: the rational chart parameter of an exposed sub-edge's witness.
[[nodiscard]] FT line_edge_parameter(LineSweep& sweep, std::size_t edge, DecisionCounters& counters);

/// P4: the witness at chart parameter t lies in the exterior of every counted
/// candidate (filtered, exact on doubt).
///
/// Raises:
///     CoverageCertificateNativeError: the witness lies in, or on, a counted
///         candidate (the counters disagreed with P4).
void confirm_line_witness(LineSweep& sweep, const FT& t, LineWitnessContext& context, DecisionCounters& counters);

/// Audit oracle for a sub-edge the counters call covered (pure exact P4).
///
/// Raises:
///     CoverageCertificateNativeError: the witness is in no candidate's interior
///         and on no covering coincident support, or on an unpredicted boundary.
void audit_covered_line_edge(LineSweep& sweep, std::size_t edge, LineWitnessContext& context, DecisionCounters& counters);

/// P7 beside the witness at chart parameter t.
///
/// Raises:
///     CoverageCertificateNativeError: P7 exceeded its contract bound.
[[nodiscard]] Point_2 uncovered_beside_line(LineSweep& sweep, const FT& t, LineWitnessContext& context,
                                            const std::function<bool(const Point_2&)>& accept, DecisionCounters& counters);

}  // namespace compas_cgal::coverage_certificate
