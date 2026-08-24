#pragma once

#include <stdexcept>

class Stock2;

// The engaged rim sub-arcs PARTITION the cutter circle, so their reported spans
// sum to at most one full turn and no assembled run can report more. engagement_at
// raises this when one does, which means a sub-arc span was normalized wrong.
//
// The invariant lives entirely on the REPORTING doubles -- no verdict consults
// them (see EngagementSample) -- but it is RAISED rather than CGAL_asserted
// because this project ships CMAKE_BUILD_TYPE Release, where CGAL_assertion
// compiles out and would leave the shipped build unguarded. The number it guards
// is the headline of every engagement report, and the failure mode it exists to
// catch is silent inflation, not a crash.
struct RimSpanNormalizationError : std::logic_error {
    using std::logic_error::logic_error;
};

class NonFiniteEngagementInputError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

class NonPositiveEngagementToolRadiusError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

class InvalidEngagementCapRatioError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

class InvalidEngagementGapRatioError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

class NonFiniteSegmentCertificationInputError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

class NonPositiveSegmentCertificationToolRadiusError
    : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

class InvalidSegmentCertificationCapError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

class NonFiniteMixedRadicalInputError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

class InvalidMixedRadicalRootError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

void validate_engagement_input_binary64(
    double cx,
    double cy,
    double tool_radius,
    double cap_chord_ratio,
    double gap_close_ratio);
void validate_segment_certification_input_binary64(
    double x0,
    double y0,
    double x1,
    double y1,
    double tool_radius,
    double cap_radians);
void validate_mixed_radical_input_binary64(
    double a,
    double b,
    double c,
    double d,
    double alpha,
    double beta);

// Convert the ergonomic angle cap exactly once at the native boundary. The
// returned binary64 value is subsequently injected into Epeck as a rational.
// Both the angle and the computed surrogate are validated; a positive angle
// whose squared sine underflows to zero is outside the representable contract.
double cap_chord_ratio(double cap_radians);

// Compare two valid binary64 cap surrogates after exact Epeck injection.
// Invalid values outside (0, 4] raise std::invalid_argument.
bool cap_chord_ratio_le(double lhs, double rhs);

// Tool-engagement-angle sample at one cutter station. Angles in radians and are
// REPORTING quantities (doubles produced for humans/statistics): total_tea sums
// the angular extent of every rim arc in contact with material, max_run_tea is
// the largest maximal contiguous engaged run. cap_exceeded is the only DECISION
// carried out here and it is decided EXACTLY on the exact arrangement, never
// from these reported doubles (see engagement_2.cpp).
//
// DECIDING / REPORTING SPLIT under gap-closure pessimism (see engagement_at):
// total_tea and max_run_tea always describe the TRUE runs (no gap closed);
// cap_exceeded is decided on the PESSIMISTIC runs (void gaps <= gap_close_ratio
// pre-absorbed). With the default gap_close_ratio == 0 the two coincide and
// every field is bit-for-bit the pre-pessimism result.
struct EngagementSample {
    double total_tea;
    double max_run_tea;
    bool cap_exceeded;
};

// Exact TEA query, read LOCALLY. Zone the cutter circle of radius tool_radius at
// (cx, cy) in the stock's OWN Gps arrangement (faces carry contained() =
// material) and harvest the engaged rim sub-arcs (the rim-in-material pieces),
// assemble maximal engaged runs by EXACT endpoint equality, and certify each run
// against the engagement cap. The query is O(cutter-crossings), not O(stock) --
// no whole-stock overlay (see engagement_2.cpp; the earlier overlay was validated
// to reproduce this exactly before it was removed).
//
// cap_chord_ratio is the dimensionless squared-chord surrogate for the angular
// cap: cap_chord_ratio = 4*sin^2(cap/2), which the CALLER computes as a double
// and which is injected exactly (Epeck::FT(double)). Contractually the cap is
// 0 < cap <= pi, so cap_chord_ratio lies in (0, 4]; values outside that range
// raise std::invalid_argument. The certificate is an exact statement about the
// exact rational threshold T = FT(cap_chord_ratio) * FT(tool_radius)^2; the
// sub-ulp gap between this rational surrogate and the transcendental angle cap
// is API semantics documented here, NOT an in-core correction constant.
//
// gap_close_ratio is the analogous squared-chord surrogate 4*sin^2(gamma/2) for
// a gap-closure angle gamma in [0, pi] (so gap_close_ratio lies in [0, 4]; out
// of range raises std::invalid_argument). GAP-CLOSURE PESSIMISM: before deciding
// the cap, every VOID gap between consecutive engaged runs whose angular span is
// <= gamma is absorbed, merging its two bounding runs (chained to a fixpoint) so
// the cap DECISION sees the PESSIMISTIC runs. This closes the certificate's merge
// hole: two runs a hair of void apart at a station are treated as already merged,
// so a merge that completes within a certifier step cannot slip past unseen. The
// closure test is the SAME exact orientation+chord predicate used for runs (a
// gap absorbed iff its span does not exceed gamma), applied to the gap's exact
// one-root endpoints -- no angle sums. The default 0.0 (gamma = 0) closes no gap:
// pessimistic runs == true runs, and every result is the pre-pessimism value
// bit-for-bit. Reported total_tea/max_run_tea always describe the TRUE runs.
//
// Raises RimSpanNormalizationError if the assembled runs report more engagement
// than a full turn, which the cutter circle cannot offer.
EngagementSample engagement_at(const Stock2& stock, double cx, double cy,
                               double tool_radius, double cap_chord_ratio,
                               double gap_close_ratio = 0.0);

// Result of certifying the engagement cap along one linear cutter motion.
// max_tea is the largest run TEA seen at any station visited -- a REPORTING
// double, best-effort over the sampled stations, never a decision input.
// cap_certified is the DECISION: true iff no cutter center on the motion can
// exceed the cap, established purely from EXACT station verdicts plus the
// analytic guard (see certify_segment_tea); false means a station violated the
// guarded cap and the margin could not be closed by refinement. stations counts
// the sub-intervals examined.
struct CertifiedTea {
    double max_tea;
    bool cap_certified;
    int stations;
};

// Certify TEA(P) <= cap_radians for EVERY cutter center P on the segment
// (x0,y0)->(x1,y1) with tool radius tool_radius, against the frozen stock.
//
// Method -- adaptive station sampling with a guarded exact test. Each station is
// measured by the EXACT engagement_at cap predicate, but against a GUARDED cap:
// at half-station-spacing d the station threshold is cap - guard(d), where guard
// is a conservative analytic TEA-growth bound carrying an explicit integer
// safety factor (derivation in engagement_2.cpp). Every center lies within its
// half-spacing of the nearer measured station, so two stations both under the
// guarded cap certify the whole span at TEA <= cap. If the guarded cap is
// exceeded -- or is non-positive, the spacing being too coarse to admit any
// guard -- the span is bisected; on reaching the spacing floor with the margin
// still open the motion is reported uncertified. The verdict is thus EXACT
// station predicates + the analytic guard ONLY; max_tea is reporting.
//
// BOUNDARY (docs/exactness.md, boundary doctrine): cap_radians is validated to
// (0, pi] and converted to exact chord surrogates here, at the one declared
// seam; every value crossing into the exact station test is injected exactly.
// Raises std::invalid_argument if cap_radians lies outside (0, pi].
CertifiedTea certify_segment_tea(const Stock2& stock, double x0, double y0,
                                 double x1, double y1, double tool_radius,
                                 double cap_radians);

// Binding-only diagnostic seam for the mixed-radical exact sign primitive.
// It carries only binary64 ingress and an integer result, keeping nanobind out
// of the deciding core translation unit.
int sign_mixed_radical_for_binding(
    double a,
    double b,
    double c,
    double d,
    double alpha,
    double beta);
