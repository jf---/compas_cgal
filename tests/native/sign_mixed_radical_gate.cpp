// Correctness gate for the mixed-radical sign predicate.
//
// `compas_cgal::exact::sign_mixed_radical` returns the exact sign of
//     a + b*sqrt(alpha) + c*sqrt(beta) + d*sqrt(alpha*beta)
// for rational a, b, c, d and rational alpha, beta >= 0. It is a deciding
// predicate on the exact cap certificate (engagement_2.cpp) and the certified
// audit station (audit_exact_station_2.cpp), which since 541b6ede both call
// this single definition. It therefore needs a direct native test of its own.
//
// This gate is the surviving half of the Task 1 equivalence gate (e54caf77).
// That gate had three arms: two textually distinct copies of the predicate and
// an independent oracle. 541b6ede merged the copies, which made the
// copy-vs-copy comparison a function compared to itself; the ORACLE arm is not
// tautological and is what is kept here, now measuring the merged definition.
// The corpus (83,850 inputs over four arms) and the branch accounting are
// unchanged from e54caf77.
//
// THE ORACLE. The expression is u + w*sqrt(beta) with u, w in Q(sqrt alpha),
// which IS the nested one-root number Sqrt_extension<CoordNT, FT>. Asking CGAL
// for that value's sign runs CGAL's own implementation, never this repository's,
// so agreement is evidence of correctness and not of a self-consistent mistake.
// ACDE_TAG is Tag_true (the arrangement's tag, and CoordNT's), which is what
// permits a root of zero -- the degenerate cases the predicate short-circuits.
#include "exact/one_root.h"

#include <CGAL/Sqrt_extension.h>
#include <CGAL/enum.h>
#include <CGAL/number_utils.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace exact = compas_cgal::exact;

namespace {

// The predicate's own vocabulary: exact::Rational is Epeck::FT and exact::OneRoot
// IS Gps_circle_segment_traits_2<Epeck>::Point_2::CoordNT, the identity
// exact_one_root_gate static_asserts.
using FT = exact::Rational;
using CoordNT = exact::OneRoot;

/// Raised when a gate check fails. Named so the gate obeys the same
/// named-exceptions-only rule as the code it exercises.
class GateCheckFailedError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// The project builds with NDEBUG, which erases <cassert>. A gate whose checks
// vanish reports a vacuous pass, so every decision here goes through a helper
// that is never compiled out.
void require(bool condition, const std::string& message)
{
    if (!condition) {
        throw GateCheckFailedError(message);
    }
}

// ---------------------------------------------------------------------------
// Independent oracle
// ---------------------------------------------------------------------------
using NestedOneRoot = CGAL::Sqrt_extension<CoordNT, FT, CGAL::Tag_true>;

// A one-root value carrying a ZERO root is representable under
// ACDE_TAG == Tag_true, but CGAL's own sign is NOT sound on it: with a0 == 0,
// Sqrt_extension::sign_ returns sign(a1) directly
// (CGAL/Sqrt_extension/Sqrt_extension_type.h:303), which is the sign of
// a1*sqrt(root) only when root > 0. At root == 0 the value IS a0 == 0 and the
// answer must be ZERO. Building this oracle naively reported NEGATIVE for
// a=-1, b=-1, c=1, d=-1, alpha=0, beta=1 -- an expression whose value is 0 --
// on 1315 of 83850 probes, every one of them a degenerate root.
//
// This is precisely why the predicate folds a degenerate root away with an
// exact FT compare before any extension is built (exact/one_root.cpp:65-77).
// The oracle folds it too, or it measures CGAL's unsound path instead of the
// predicate, and would then report the CORRECT predicate as red. Folding
// sqrt(0) == 0 is not the algorithm under test: the u/w decomposition and the
// magnitude comparison the predicate is checked on remain CGAL's
// implementation, never a transcription of the predicate's body.
CoordNT one_root(const FT& a0, const FT& a1, const FT& root)
{
    if (CGAL::is_zero(root)) return CoordNT(a0);
    return CoordNT(a0, a1, root);
}

CGAL::Sign oracle_sign(
    const FT& a,
    const FT& b,
    const FT& c,
    const FT& d,
    const FT& alpha,
    const FT& beta)
{
    const CoordNT u = one_root(a, b, alpha);
    const CoordNT w = one_root(c, d, alpha);
    if (CGAL::is_zero(beta)) {
        return CGAL::sign(u);   // u + w*sqrt(0) == u
    }
    return CGAL::sign(NestedOneRoot(u, w, beta));
}

// ---------------------------------------------------------------------------
// Branch accounting -- REPORTING ONLY
// ---------------------------------------------------------------------------
// This classifier labels which branch an input reaches. It never produces a
// sign and never decides whether the gate passes; it exists so the gate can
// prove it actually exercised every branch instead of passing vacuously.
enum class Branch : std::size_t {
    BOTH_ROOTS_RATIONAL = 0,
    BETA_RATIONAL,
    ALPHA_RATIONAL,
    EQUAL_ROOTS,
    W_ZERO,
    U_ZERO,
    LIKE_SIGNS,
    OPPOSITE_U_LARGER,
    OPPOSITE_W_LARGER,
    OPPOSITE_EQUAL,
    COUNT
};

const char* branch_text(Branch branch)
{
    switch (branch) {
        case Branch::BOTH_ROOTS_RATIONAL: return "both roots rational";
        case Branch::BETA_RATIONAL:       return "beta rational";
        case Branch::ALPHA_RATIONAL:      return "alpha rational";
        case Branch::EQUAL_ROOTS:         return "alpha == beta";
        case Branch::W_ZERO:              return "w == 0";
        case Branch::U_ZERO:              return "u == 0";
        case Branch::LIKE_SIGNS:          return "like signs";
        case Branch::OPPOSITE_U_LARGER:   return "opposite signs, |u| larger";
        case Branch::OPPOSITE_W_LARGER:   return "opposite signs, |w*sqrt(beta)| larger";
        case Branch::OPPOSITE_EQUAL:      return "opposite signs, equal magnitude";
        case Branch::COUNT:               break;
    }
    return "<unknown>";
}

Branch classify(
    const FT& a,
    const FT& b,
    const FT& c,
    const FT& d,
    const FT& alpha,
    const FT& beta)
{
    if (CGAL::is_zero(alpha) && CGAL::is_zero(beta)) return Branch::BOTH_ROOTS_RATIONAL;
    if (CGAL::is_zero(beta)) return Branch::BETA_RATIONAL;
    if (CGAL::is_zero(alpha)) return Branch::ALPHA_RATIONAL;
    if (alpha == beta) return Branch::EQUAL_ROOTS;

    const CoordNT u(a, b, alpha);
    const CoordNT w(c, d, alpha);
    const CGAL::Sign u_sign = CGAL::sign(u);
    const CGAL::Sign w_sign = CGAL::sign(w);
    if (w_sign == CGAL::ZERO) return Branch::W_ZERO;
    if (u_sign == CGAL::ZERO) return Branch::U_ZERO;
    if (u_sign == w_sign) return Branch::LIKE_SIGNS;
    switch (CGAL::compare(u * u, w * w * CoordNT(beta))) {
        case CGAL::LARGER:  return Branch::OPPOSITE_U_LARGER;
        case CGAL::SMALLER: return Branch::OPPOSITE_W_LARGER;
        default:            return Branch::OPPOSITE_EQUAL;
    }
}

// ---------------------------------------------------------------------------
// Probing
// ---------------------------------------------------------------------------
const char* sign_text(CGAL::Sign value)
{
    switch (value) {
        case CGAL::NEGATIVE: return "NEGATIVE";
        case CGAL::ZERO:     return "ZERO";
        case CGAL::POSITIVE: return "POSITIVE";
    }
    return "<unknown>";
}

// REPORTING ONLY: exact text of a rational, for a failure message. Never feeds
// a decision.
std::string exact_text(const FT& value)
{
    std::ostringstream out;
    out << CGAL::exact(value);
    return out.str();
}

struct ProbeLedger {
    std::size_t compared = 0;
    std::size_t mismatched = 0;
    std::array<std::size_t, static_cast<std::size_t>(Branch::COUNT)> branch_hits{};
    std::array<std::size_t, static_cast<std::size_t>(Branch::COUNT)> branch_misses{};
    std::string first_failure;
};

std::string describe(
    const FT& a,
    const FT& b,
    const FT& c,
    const FT& d,
    const FT& alpha,
    const FT& beta,
    CGAL::Sign predicate,
    CGAL::Sign oracle,
    Branch branch);

void probe(
    const FT& a,
    const FT& b,
    const FT& c,
    const FT& d,
    const FT& alpha,
    const FT& beta,
    ProbeLedger& ledger)
{
    const CGAL::Sign predicate = exact::sign_mixed_radical(a, b, c, d, alpha, beta);
    const CGAL::Sign oracle = oracle_sign(a, b, c, d, alpha, beta);
    const Branch branch = classify(a, b, c, d, alpha, beta);

    ++ledger.compared;
    ++ledger.branch_hits[static_cast<std::size_t>(branch)];

    if (predicate != oracle) {
        ++ledger.mismatched;
        ++ledger.branch_misses[static_cast<std::size_t>(branch)];
        if (ledger.first_failure.empty()) {
            ledger.first_failure =
                describe(a, b, c, d, alpha, beta, predicate, oracle, branch);
        }
    }
}

std::string describe(
    const FT& a,
    const FT& b,
    const FT& c,
    const FT& d,
    const FT& alpha,
    const FT& beta,
    CGAL::Sign predicate,
    CGAL::Sign oracle,
    Branch branch)
{
    std::ostringstream out;
    out << "sign disagreement on a=" << exact_text(a) << " b=" << exact_text(b)
        << " c=" << exact_text(c) << " d=" << exact_text(d)
        << " alpha=" << exact_text(alpha) << " beta=" << exact_text(beta)
        << " -- exact::sign_mixed_radical = " << sign_text(predicate)
        << ", CGAL nested-extension oracle = " << sign_text(oracle)
        << ", branch = " << branch_text(branch);
    return out.str();
}

FT rational(long numerator, long denominator)
{
    return FT(numerator) / FT(denominator);
}

// ---------------------------------------------------------------------------
// Arm 1 -- the prescribed exhaustive sign structure: 3^4 * 5^2 = 2025 inputs.
// alpha and beta range over 0..4, so alpha == 0, beta == 0, alpha == beta and
// perfect-square roots are all included.
// ---------------------------------------------------------------------------
void exhaustive_sign_structure(ProbeLedger& ledger)
{
    const std::array<FT, 3> coefficients{FT(-1), FT(0), FT(1)};
    for (const FT& a : coefficients) {
        for (const FT& b : coefficients) {
            for (const FT& c : coefficients) {
                for (const FT& d : coefficients) {
                    for (int alpha = 0; alpha <= 4; ++alpha) {
                        for (int beta = 0; beta <= 4; ++beta) {
                            probe(a, b, c, d, FT(alpha), FT(beta), ledger);
                        }
                    }
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Arm 2 -- the same grid widened so magnitudes cross: fractional and larger
// coefficients let |u| and |w*sqrt(beta)| overtake each other, and let u
// vanish over a perfect-square root (1 - (1/2)*sqrt(4) == 0).
// 7^4 * 5^2 = 60025 inputs.
// ---------------------------------------------------------------------------
void exhaustive_magnitude_grid(ProbeLedger& ledger)
{
    const std::array<FT, 7> coefficients{
        FT(-3), FT(-1), rational(-1, 2), FT(0), rational(1, 2), FT(1), FT(3)};
    for (const FT& a : coefficients) {
        for (const FT& b : coefficients) {
            for (const FT& c : coefficients) {
                for (const FT& d : coefficients) {
                    for (int alpha = 0; alpha <= 4; ++alpha) {
                        for (int beta = 0; beta <= 4; ++beta) {
                            probe(a, b, c, d, FT(alpha), FT(beta), ledger);
                        }
                    }
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Arm 3 -- inputs constructed to sit exactly ON the magnitude tie and on both
// sides of it, which is where an incorrect merge or an inexact comparison would
// show first. Choose w = c + d*sqrt(alpha) and a rational s > 0, set
// beta = s^2, and set u = -k*s*w. Then u^2 = k^2 * beta * w^2, so the exact
// comparison is EQUAL at k == 1 and strictly to one side otherwise.
// The k values are exact rational INPUTS chosen to straddle the tie, not
// tolerances: nothing in the gate or the predicate compares against them.
// ---------------------------------------------------------------------------
void magnitude_tie_neighbourhood(ProbeLedger& ledger)
{
    const std::array<int, 5> alphas{2, 3, 5, 6, 7};
    const std::array<FT, 6> tie_ratios{
        FT(1), FT(2), FT(3), rational(1, 2), rational(3, 2), rational(5, 3)};
    const std::array<std::pair<FT, FT>, 6> w_coefficients{
        std::pair<FT, FT>{FT(1), FT(0)},
        std::pair<FT, FT>{FT(0), FT(1)},
        std::pair<FT, FT>{FT(1), FT(1)},
        std::pair<FT, FT>{FT(1), FT(-1)},
        std::pair<FT, FT>{FT(2), FT(3)},
        std::pair<FT, FT>{FT(-1), FT(2)}};
    const std::array<FT, 5> straddle{
        rational(1, 2), rational(999, 1000), FT(1), rational(1001, 1000), FT(2)};

    for (int alpha_value : alphas) {
        const FT alpha(alpha_value);
        for (const FT& s : tie_ratios) {
            const FT beta = s * s;
            if (beta == alpha) continue;   // that is arm 1's EQUAL_ROOTS branch
            for (const std::pair<FT, FT>& w : w_coefficients) {
                for (const FT& k : straddle) {
                    const FT scale = k * s;
                    // Opposite signs: the branch the tie comparison decides.
                    probe(-scale * w.first, -scale * w.second, w.first, w.second,
                          alpha, beta, ledger);
                    // Like signs at the same magnitudes: the branch that returns
                    // without ever forming the magnitude comparison.
                    probe(scale * w.first, scale * w.second, w.first, w.second,
                          alpha, beta, ledger);
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Arm 4 -- randomised general-case rationals, fixed seed so a failure is
// reproducible. Roots are drawn non-negative (the predicate's precondition) and
// are occasionally forced to zero or to each other so the degenerate branches
// keep being visited with non-trivial coefficients.
// ---------------------------------------------------------------------------
constexpr std::uint32_t RANDOM_ARM_SEED = 0x5A17D00Du;
constexpr int RANDOM_ARM_PROBES = 20000;

void randomised_rationals(ProbeLedger& ledger)
{
    std::mt19937 generator(RANDOM_ARM_SEED);
    std::uniform_int_distribution<int> coefficient_numerator(-12, 12);
    std::uniform_int_distribution<int> root_numerator(0, 12);
    std::uniform_int_distribution<int> denominator(1, 7);
    std::uniform_int_distribution<int> degenerate(0, 7);

    for (int probe_index = 0; probe_index < RANDOM_ARM_PROBES; ++probe_index) {
        const FT a = rational(coefficient_numerator(generator), denominator(generator));
        const FT b = rational(coefficient_numerator(generator), denominator(generator));
        const FT c = rational(coefficient_numerator(generator), denominator(generator));
        const FT d = rational(coefficient_numerator(generator), denominator(generator));
        const FT alpha = rational(root_numerator(generator), denominator(generator));
        FT beta = rational(root_numerator(generator), denominator(generator));
        const int shape = degenerate(generator);
        if (shape == 0) beta = alpha;
        if (shape == 1) beta = FT(0);
        probe(a, b, c, d, alpha, beta, ledger);
    }
}

void report_and_require(const ProbeLedger& ledger)
{
    std::printf("  compared %zu inputs\n", ledger.compared);
    for (std::size_t index = 0; index < static_cast<std::size_t>(Branch::COUNT); ++index) {
        const Branch branch = static_cast<Branch>(index);
        std::printf("  branch %-42s %8zu\n", branch_text(branch), ledger.branch_hits[index]);
    }
    std::printf("  oracle disagreements %zu\n", ledger.mismatched);
    for (std::size_t index = 0; index < static_cast<std::size_t>(Branch::COUNT); ++index) {
        if (ledger.branch_misses[index] == 0) continue;
        std::printf(
            "  oracle disagreed in branch %-30s %8zu\n",
            branch_text(static_cast<Branch>(index)),
            ledger.branch_misses[index]);
    }

    // A gate that never reached a branch would pass without testing it, so
    // branch coverage is required before the agreement result is allowed to
    // mean anything.
    for (std::size_t index = 0; index < static_cast<std::size_t>(Branch::COUNT); ++index) {
        const Branch branch = static_cast<Branch>(index);
        require(
            ledger.branch_hits[index] > 0,
            std::string("no probe reached the branch: ") + branch_text(branch));
    }
    // The gate's contract: a disagreement with CGAL's own nested-extension sign
    // is a wrong answer from a deciding predicate on the exact certificate path.
    require(ledger.first_failure.empty(), ledger.first_failure);
    require(ledger.mismatched == 0, "oracle disagreement count is non-zero");
}

}  // namespace

int main()
{
    ProbeLedger ledger;
    try {
        exhaustive_sign_structure(ledger);
        exhaustive_magnitude_grid(ledger);
        magnitude_tie_neighbourhood(ledger);
        randomised_rationals(ledger);
        report_and_require(ledger);
    } catch (const std::exception& error) {
        std::printf("sign_mixed_radical_gate FAILED: %s\n", error.what());
        return 1;
    }
    std::printf("sign_mixed_radical_gate OK\n");
    return 0;
}
