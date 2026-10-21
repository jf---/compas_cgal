#pragma once

// Reporting adapter: exact values to binary64, one-way, never a decision
// (docs/exactness.md, "the deciding / reporting split").
//
// The Python oracle (benchmarks/held_declared_clearance.py) chooses which
// primitive to REPORT among equally decided ones by binary64 views of the
// distances. These views follow its formulas so that choice agrees:
//   * a rational converts to its nearest binary64, ties to even (Python's
//     float(Fraction)), from CGAL's enclosure (narrowed by exact bisection
//     where it is not tight) and one exact comparison with its midpoint;
//   * a one-root number converts as float(a) + float(b) * sqrt(float(c)),
//     evaluated without contraction into a fused multiply-add.

#include "held_declared_clearance_2/one_root_2.h"

namespace compas_cgal::held_declared_clearance {

/// The binary64 nearest `value`, ties to even.
///
/// Raises:
///     ReportingOverflowError: |value| rounds beyond the largest finite binary64.
[[nodiscard]] double binary64(const FT& value);

/// float(a) + float(b) * sqrt(float(c)), as the Python oracle evaluates it.
[[nodiscard]] double binary64(const Root& value);

}  // namespace compas_cgal::held_declared_clearance
