#include "engagement_2_bindings.h"

#include "engagement_2.h"
#include "stock_2.h"

#include <tuple>

#include <nanobind/stl/tuple.h>

namespace nb = nanobind;
using namespace nb::literals;

void register_engagement(nanobind::module_& m)
{
    nb::exception<RimSpanNormalizationError>(
        m, "RimSpanNormalizationError", PyExc_RuntimeError);
    nb::exception<NonFiniteEngagementInputError>(
        m, "NonFiniteEngagementInputError", PyExc_ValueError);
    nb::exception<NonPositiveEngagementToolRadiusError>(
        m, "NonPositiveEngagementToolRadiusError", PyExc_ValueError);
    nb::exception<InvalidEngagementCapRatioError>(
        m, "InvalidEngagementCapRatioError", PyExc_ValueError);
    nb::exception<InvalidEngagementGapRatioError>(
        m, "InvalidEngagementGapRatioError", PyExc_ValueError);
    nb::exception<NonFiniteSegmentCertificationInputError>(
        m, "NonFiniteSegmentCertificationInputError", PyExc_ValueError);
    nb::exception<NonPositiveSegmentCertificationToolRadiusError>(
        m,
        "NonPositiveSegmentCertificationToolRadiusError",
        PyExc_ValueError);
    nb::exception<InvalidSegmentCertificationCapError>(
        m, "InvalidSegmentCertificationCapError", PyExc_ValueError);
    nb::exception<NonFiniteMixedRadicalInputError>(
        m, "NonFiniteMixedRadicalInputError", PyExc_ValueError);
    nb::exception<InvalidMixedRadicalRootError>(
        m, "InvalidMixedRadicalRootError", PyExc_ValueError);

    m.def("cap_chord_ratio", &cap_chord_ratio, "cap_radians"_a);
    m.def("cap_chord_ratio_le", &cap_chord_ratio_le, "lhs"_a, "rhs"_a);

    // Nanobind boundary. cx, cy, tool_radius are measured/computed station data:
    // each double IS a rational and enters exact-land by exact injection
    // (Epeck::FT) -- no snapping, no tolerance at the seam. cap_chord_ratio is
    // the caller's exact rational surrogate 4*sin^2(cap/2) for the transcendental
    // cap (docs/exactness.md "Input semantics" and "boundary doctrine").
    m.def("engagement_at",
          [](const Stock2& stock, double cx, double cy, double tool_radius,
             double cap_chord_ratio, double gap_close_ratio) {
              EngagementSample s = engagement_at(stock, cx, cy, tool_radius,
                                                 cap_chord_ratio, gap_close_ratio);
              return std::make_tuple(s.total_tea, s.max_run_tea, s.cap_exceeded);
          },
          "stock"_a, "cx"_a, "cy"_a, "tool_radius"_a, "cap_chord_ratio"_a,
          "gap_close_ratio"_a = 0.0);

    // Nanobind boundary for the motion certificate: x0..y1, tool_radius and the
    // angular cap enter exact-land inside certify_segment_tea (validation + exact
    // chord-surrogate injection). Returns (max_tea, cap_certified, stations),
    // the CertifiedTea fields, matching the tuple style of engagement_at.
    m.def("certify_segment_tea",
          [](const Stock2& stock, double x0, double y0, double x1, double y1,
             double tool_radius, double cap_radians) {
              CertifiedTea c = certify_segment_tea(stock, x0, y0, x1, y1,
                                                   tool_radius, cap_radians);
              return std::make_tuple(c.max_tea, c.cap_certified, c.stations);
          },
          "stock"_a, "x0"_a, "y0"_a, "x1"_a, "y1"_a, "tool_radius"_a, "cap_radians"_a);

    // Test-only: exact sign of A + B*sqrt(alpha) + C*sqrt(beta) + D*sqrt(alpha*beta)
    // (returns -1/0/+1) so the cap predicate's core primitive is unit-tested
    // directly against high-precision references.
    m.def("_sign_mixed_radical",
          [](double a, double b, double c, double d, double alpha, double beta) {
              return static_cast<int>(
                  sign_mixed_radical_for_binding(a, b, c, d, alpha, beta));
          },
          "a"_a, "b"_a, "c"_a, "d"_a, "alpha"_a, "beta"_a);
}
