#pragma once

// Exact side decisions at a witness that lies ON a support (P5/P6 for regions
// whose boundary is not a single circle or line: exact connector regions and the
// exact target T).
//
// An edge e of a checked piece is cut at every support of a region Y (or of T),
// so on each side of e, close to e, Y-membership is constant: it is the
// membership of the face of the support arrangement that borders e there. When
// the witness w of e lies on one of Y's supports (a coincident support), the
// closed membership at w does not tell the sides apart. The probe point
// q = w + h dir (dir pointing to the checked side, h > 0 rational) lies in that
// bordering face as soon as the half-open segment (w, q] meets no support of Y:
// then Y-membership at q (never Boundary) is the side's membership.
//
// h is decided exactly: every support's crossing parameters tau with the ray
// w + tau dir are exact (rational, one-root or a tower m + n sqrt(d) over the
// support's radicand), and h halves from 1 until h < tau for every positive tau.
// A tau == 0 is w itself (the coincident support or one through w) and imposes
// nothing. The halving is bounded by the root-separation bound of the smallest
// positive tau (separation_2.h). LAZY filters first: a tau whose certified
// enclosure lies above 0 is cleared by any h below the enclosure's lower end,
// and one whose enclosure lies below 0 imposes nothing; only the rest are
// decided exactly.

#include "coverage_certificate_2/counters_2.h"
#include "coverage_certificate_2/kernel_2.h"
#include "coverage_certificate_2/lines_2.h"
#include "coverage_certificate_2/vertex_events_2.h"

#include <utility>
#include <vector>

namespace compas_cgal::coverage_certificate {

/// The supports a probe segment must not meet.
struct ProbeSupports {
    std::vector<const LineSupport*> lines;
    std::vector<std::pair<Point_2, FT>> circles;  ///< rational centre, squared radius
    std::vector<const VertexCircle*> vertex_circles;
};

/// A rational h in (0, 1] such that the half-open segment (w, w + h dir] meets
/// none of `supports`. `mode` selects the filtered (LAZY) or exact-only path.
///
/// Raises:
///     SeparationBoundExceededError: the halving ran past the separation bound
///         of the smallest positive crossing parameter (a wrong sign upstream).
[[nodiscard]] FT probe_step(const Point_2& w, const FT& dx, const FT& dy, const ProbeSupports& supports,
                            FilterMode mode, DecisionCounters& counters);

/// The probe point w + h dir of probe_step.
[[nodiscard]] Point_2 probe_point(const Point_2& w, const FT& dx, const FT& dy, const ProbeSupports& supports,
                                  FilterMode mode, DecisionCounters& counters);

/// dir scaled by a power of two so that its larger |component| lies in
/// [2^PROBE_DIRECTION_LOG2 / 2, 2^PROBE_DIRECTION_LOG2]: a construction
/// parameter only (any positive multiple of dir points to the same side), which
/// keeps the probe segment, and so its box-filtered support candidates, short.
[[nodiscard]] std::pair<FT, FT> probe_direction(const FT& dx, const FT& dy);

/// log2 of the probe direction's larger |component|, mm (2^-10 mm).
inline constexpr int PROBE_DIRECTION_LOG2 = -10;

}  // namespace compas_cgal::coverage_certificate
