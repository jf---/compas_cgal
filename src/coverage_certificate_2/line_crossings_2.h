#pragma once

// P1 for the capsule pieces: crossings of a rational circle with a line support,
// and of two line supports, behind the same certified interval front as the
// circle-circle core (lazy_crossings_2.h), with the depth-2 tower (tower_2.h) as
// the exact back end.
//
// Circle C (centre c, radius rho, rational) and line l with chart
// x(t) = b0 + b1 sqrt(L) + t d (lines_2.h, d . b1 = 0):
//   |x(t) - c|^2 = rho^2  <=>  Ld t^2 + 2 beta t + gamma = 0,
//   Ld = |d|^2, beta = d.(b0 - c) (rational), gamma = g0 + g1 sqrt(L),
//   disc = beta^2 - Ld gamma in Q(sqrt(L)),
//   t = -beta / Ld +/- sqrt(disc) / Ld          (branch 0 is +).
// The point is x(t): both coordinates m + n sqrt(disc) with m in Q(sqrt(L)), n
// rational. disc < 0: no crossing; disc == 0: a tangency (double root, flips
// nothing); disc > 0: two transversal crossings.
//
// Piece line P and line Q (n_Q . x = k0 + k1 sqrt(L_Q)):
//   t (n_Q . d_P) = k0 + k1 sqrt(L_Q) - n_Q . b0 - (n_Q . b1) sqrt(L_P),
// a tower value with m in Q(sqrt(L_P)), n rational, radicand L_Q. n_Q . d_P == 0
// means parallel: coincident or a constant side, decided exactly.

#include "coverage_certificate_2/counters_2.h"
#include "coverage_certificate_2/kernel_2.h"
#include "coverage_certificate_2/lines_2.h"
#include "coverage_certificate_2/regions_2.h"
#include "coverage_certificate_2/tower_2.h"

#include <array>
#include <cstdint>
#include <memory>
#include <optional>

namespace compas_cgal::coverage_certificate {

/// Exact crossings of a circle and a line. For count == 1 (tangency) both
/// branches hold the same point.
struct CircleLineExact {
    std::uint8_t count = 0;  ///< 0, 1 (tangent, multiplicity 2) or 2 (transversal)
    std::array<Tower, 2> t, x, y;
};

[[nodiscard]] CircleLineExact exact_circle_line(const BoundaryCircle& circle, const LineSupport& line);

/// Certified enclosures of the crossings; `count` is set only when the
/// discriminant enclosure excludes zero (0 or 2 crossings), otherwise empty.
struct CircleLineFilter {
    std::optional<std::uint8_t> count;
    std::array<Interval, 2> t, x, y;
};

[[nodiscard]] CircleLineFilter filter_circle_line(const BoundaryCircle& circle, const LineSupport& line);

/// One (circle, line) crossing record: filter result now, exact towers on demand.
struct CircleLineRecord {
    std::uint32_t circle;
    std::uint32_t line;
    std::uint8_t count;
    std::array<Interval, 2> t, x, y;
    std::unique_ptr<CircleLineExact> exact;  ///< heap: a tower point is kilobytes, most records never need one
};

/// Build a record: the filter decides the count when it can (LAZY), the exact
/// crossings otherwise (or always in EXACT_ONLY).
[[nodiscard]] CircleLineRecord make_circle_line_record(std::uint32_t circle_id, const BoundaryCircle& circle,
                                                       std::uint32_t line_id, const LineSupport& line, FilterMode mode,
                                                       DecisionCounters& counters);

/// The exact crossings of a record (materialised on demand). Checks that every
/// filter enclosure contains its exact value (a contract on the filter).
///
/// Raises:
///     CoverageCertificateNativeError: an exact crossing lies outside its
///         certified enclosure, or the exact count disagrees with the filter.
const CircleLineExact& exact_of(CircleLineRecord& record, const BoundaryCircle& circle, const LineSupport& line,
                                DecisionCounters& counters);

/// The relation of line Q to piece line P.
enum class LineRelation : std::uint8_t { Crossing = 0, Parallel = 1, Coincident = 2 };

struct LineLineRecord {
    std::uint32_t other;
    LineRelation relation;
    int coincidence;          ///< +1/-1 normal orientation for a coincident pair
    CGAL::Sign parallel_side; ///< side of P relative to Q for a parallel pair
    Interval t;               ///< crossing parameter on P
    std::unique_ptr<Tower> exact_t;
};

/// Build a line-line record on piece line P (filtered parallelism, exact on doubt).
[[nodiscard]] LineLineRecord make_line_line_record(const LineSupport& piece, std::uint32_t other_id,
                                                   const LineSupport& other, bool same_support, FilterMode mode,
                                                   DecisionCounters& counters);

/// The exact crossing parameter of a Crossing record (materialised on demand).
const Tower& exact_t_of(LineLineRecord& record, const LineSupport& piece, const LineSupport& other,
                        DecisionCounters& counters);

}  // namespace compas_cgal::coverage_certificate
