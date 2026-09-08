// Equivalence gate for the two copies of the mixed-radical sign predicate.
//
// The predicate returns the exact sign of
//     a + b*sqrt(alpha) + c*sqrt(beta) + d*sqrt(alpha*beta)
// for rational a, b, c, d and rational alpha, beta >= 0. Two copies of it
// exist and they differ textually, so neither may be merged into the other
// until they are shown to return the same sign on the same input.
//
//   Copy A -- src/engagement_2.cpp:43, `sign_mixed_radical`, serving the
//             engagement cap certificate directly at :123 and :141.
//   Copy B -- src/audit_exact_station_2.cpp:22, `sign_mixed_radical_impl`,
//             published as `audit_sign_mixed_radical_exact` and reached by the
//             Python binding `_sign_mixed_radical`.
//
// Copy B MERGES two of copy A's branches:
//     A:  if (su == ZERO) return sw;      if (su == sw) return su;
//     B:  if (u_sign == ZERO || u_sign == w_sign) return w_sign;
// The merge is sound only if `u_sign == w_sign` really holds wherever the
// second disjunct fires, so that returning `w_sign` and returning `u_sign` name
// the same value. This gate settles that by measurement rather than by reading.
//
// HOW COPY A IS REACHED. Copy A sits in an anonymous namespace
// (src/engagement_2.cpp:22) and therefore has internal linkage: no declaration
// can name it. Its translation unit is compiled into this gate so that the
// SHIPPED function is the one called -- transcribing its body here would only
// prove a transcription equivalent to copy B, which is not the question. There
// is no ODR conflict: src/engagement_2.cpp belongs to the `_stock_2` nanobind
// module target (CMakeLists.txt:441) and is NOT a member of
// `continuous_tea_exact_core`, the only library this gate links, so this
// executable never links a second copy of that translation unit.
#include "engagement_2.cpp"  // NOLINT(bugprone-suspicious-include)

#include "audit_certification_2.h"

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
#include <vector>

// `FT` and `CoordNT` are the aliases of the included translation unit's
// anonymous namespace; unqualified lookup from this named namespace reaches
// them, and `sign_mixed_radical` with them.
namespace gate {

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
// The expression is u + w*sqrt(beta) with u, w in Q(sqrt alpha), which IS the
// nested one-root number Sqrt_extension<CoordNT, FT>. Asking CGAL for that
// value's sign runs CGAL's own implementation, not either copy, so a three-way
// agreement is evidence of correctness and not merely of a shared mistake.
// ACDE_TAG is Tag_true (the arrangement's tag, and CoordNT's), which is what
// permits a root of zero -- the degenerate cases both copies short-circuit.
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
// This is precisely why BOTH copies fold a degenerate root away with an exact
// FT compare before any extension is built (engagement_2.cpp:50-52). The
// oracle folds it too, or it measures CGAL's unsound path instead of the
// predicate. Folding sqrt(0) == 0 is not the algorithm under test: the u/w
// decomposition and the magnitude comparison the copies are checked on remain
// CGAL's implementation, never a transcription of copy A.
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
// prove it actually exercised the merged branch instead of passing vacuously.
enum class Branch : std::size_t {
    BOTH_ROOTS_RATIONAL = 0,
    BETA_RATIONAL,
    ALPHA_RATIONAL,
    EQUAL_ROOTS,
    W_ZERO,
    U_ZERO,             // copy A returns sw here; copy B folds it into w_sign
    LIKE_SIGNS,         // copy A returns su here; copy B returns w_sign
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
        case Branch::U_ZERO:              return "u == 0 (merged in copy B)";
        case Branch::LIKE_SIGNS:          return "like signs (merged in copy B)";
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
    std::size_t copy_mismatched = 0;
    std::size_t oracle_mismatched = 0;
    std::array<std::size_t, static_cast<std::size_t>(Branch::COUNT)> branch_hits{};
    std::array<std::size_t, static_cast<std::size_t>(Branch::COUNT)> oracle_misses{};
    std::string first_copy_failure;
    std::string first_oracle_failure;
};

std::string describe(
    const FT& a,
    const FT& b,
    const FT& c,
    const FT& d,
    const FT& alpha,
    const FT& beta,
    CGAL::Sign copy_a,
    CGAL::Sign copy_b,
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
    const CGAL::Sign copy_a = sign_mixed_radical(a, b, c, d, alpha, beta);
    const CGAL::Sign copy_b =
        audit_sign_mixed_radical_exact(a, b, c, d, alpha, beta);
    const CGAL::Sign oracle = oracle_sign(a, b, c, d, alpha, beta);
    const Branch branch = classify(a, b, c, d, alpha, beta);

    ++ledger.compared;
    ++ledger.branch_hits[static_cast<std::size_t>(branch)];

    if (copy_a != copy_b) {
        ++ledger.copy_mismatched;
        if (ledger.first_copy_failure.empty()) {
            ledger.first_copy_failure =
                describe(a, b, c, d, alpha, beta, copy_a, copy_b, oracle, branch);
        }
    }
    if (copy_a != oracle || copy_b != oracle) {
        ++ledger.oracle_mismatched;
        ++ledger.oracle_misses[static_cast<std::size_t>(branch)];
        if (ledger.first_oracle_failure.empty()) {
            ledger.first_oracle_failure =
                describe(a, b, c, d, alpha, beta, copy_a, copy_b, oracle, branch);
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
    CGAL::Sign copy_a,
    CGAL::Sign copy_b,
    CGAL::Sign oracle,
    Branch branch)
{
    std::ostringstream out;
    out << "sign disagreement on a=" << exact_text(a) << " b=" << exact_text(b)
        << " c=" << exact_text(c) << " d=" << exact_text(d)
        << " alpha=" << exact_text(alpha) << " beta=" << exact_text(beta)
        << " -- copy A (engagement_2.cpp) = " << sign_text(copy_a)
        << ", copy B (audit_exact_station_2.cpp) = " << sign_text(copy_b)
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
                    // Like signs at the same magnitudes: the merged branch.
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
    std::printf(
        "  copy disagreements %zu, oracle disagreements %zu\n",
        ledger.copy_mismatched,
        ledger.oracle_mismatched);
    for (std::size_t index = 0; index < static_cast<std::size_t>(Branch::COUNT); ++index) {
        if (ledger.oracle_misses[index] == 0) continue;
        std::printf(
            "  oracle disagreed in branch %-30s %8zu\n",
            branch_text(static_cast<Branch>(index)),
            ledger.oracle_misses[index]);
    }

    // A gate that never reached the merged branch would pass without testing
    // the one difference it exists to test, so branch coverage is required
    // before the agreement result is allowed to mean anything.
    for (std::size_t index = 0; index < static_cast<std::size_t>(Branch::COUNT); ++index) {
        const Branch branch = static_cast<Branch>(index);
        require(
            ledger.branch_hits[index] > 0,
            std::string("no probe reached the branch: ") + branch_text(branch));
    }
    // The gate's contract, in order of what a red result would mean: the two
    // copies disagreeing is a live wrong-answer defect in one lane; the oracle
    // disagreeing means both copies share a mistake CGAL does not.
    require(ledger.first_copy_failure.empty(), ledger.first_copy_failure);
    require(ledger.copy_mismatched == 0, "copy disagreement count is non-zero");
    require(ledger.first_oracle_failure.empty(), ledger.first_oracle_failure);
    require(ledger.oracle_mismatched == 0, "oracle disagreement count is non-zero");
}

}  // namespace gate

int main()
{
    gate::ProbeLedger ledger;
    try {
        gate::exhaustive_sign_structure(ledger);
        gate::exhaustive_magnitude_grid(ledger);
        gate::magnitude_tie_neighbourhood(ledger);
        gate::randomised_rationals(ledger);
        gate::report_and_require(ledger);
    } catch (const std::exception& error) {
        std::printf("sign_mixed_radical_equivalence_gate FAILED: %s\n", error.what());
        return 1;
    }
    std::printf("sign_mixed_radical_equivalence_gate OK\n");
    return 0;
}
