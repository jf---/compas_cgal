#pragma once

#include "stock_2.h"
#include "stock_local_2.h"

#include <stdexcept>

// ----------------------------------------------------------------------------
// The replay's ONE line-sweep model: what a tool of radius r removes moving
// along the segment A -> B, built in the circle-segment traits.
//
// The true swept region is the capsule A B (+) disk(r) = disk(A, r) U
// rect(A, B, r) U disk(B, r). The end disks are always exact (rational centre,
// rational radius). The rectangle's sides stand off the axis by r * n / L with
// L = |B - A|, rational EXACTLY when dx^2 + dy^2 is a square in Q -- decided
// exactly on the injected endpoints by is_pythagorean_direction:
//
//   DISK             A == B: the plain disk(A, r), exactly the swept region.
//   EXACT_CAPSULE    Pythagorean direction: half-width vector r * n / L, the
//                    union IS the capsule (equality).
//   CANONICAL_QUAD   otherwise: half-width vector h = s * perp(u) with
//                    u = (B - A) / max(|dx|, |dy|) the canonical direction of
//                    the class, s chosen in doubles and CHECKED exactly so |h|
//                    lies in a band about r. Only the side band is inexact.
//
// The end disks stay at the full radius r because the tool's endpoint positions
// are exact: a plunge, or the next motion's start, at A or B then meets nothing
// the model left behind (Ruling 0.2-endcaps). Collinear cuts of one direction
// class share h exactly, so a collinear re-cut changes the stock only inside the
// earlier cut's certified OVER region (its end disks poke into the band).
//
// Each quad is a certified BOUND on the true capsule:
//   UNDER  |h| in [(1 - f) r, r]  -> region subset of the capsule; what the
//          replay removes. Everything within (1 - f) r of the segment IS removed.
//   OVER   |h| in [r, (1 + f) r]  -> region contains the capsule; used only to
//          certify that a motion removes NOTHING.
// For DISK and EXACT_CAPSULE the two coincide with the swept region itself.
// ----------------------------------------------------------------------------

// Relative half-width band of a canonical quad's rectangle, as a fraction of the tool
// radius: the same 1e-4 * r budget CHAIN_SLACK_FRACTION grants the disk chain and
// the quad (stock_2.cpp), dominated by the generator's 2e-3 * r clearance. The
// construction aims at the band's middle (1 -/+ f/2) r, leaving f/2 = 5e-5 of
// relative headroom against ~2e-16 of double rounding in the scale.
inline constexpr double LINE_SWEEP_BAND_FRACTION = 1e-4;

// The tool radius of a line sweep is zero or negative: an argument fault at the
// double boundary (ValueError). A non-finite coordinate or radius reuses
// NonFiniteCapsuleInputError, the failure mode the quad already names.
class InvalidLineSweepRadiusError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// The rectangle's exact half-width left its certified band (or its double scale
// was not a positive finite number): a broken INTERNAL invariant (RuntimeError).
class LineSweepBandCertificateError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// The rectangle's exact corners are not counterclockwise: a broken internal
// invariant of the construction (RuntimeError), never an argument fault.
class DegenerateLineSweepRectangleError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// The constructed sweep region is not ONE simply connected polygon. A disk-
// rect-disk union always is, so this is a broken internal invariant (RuntimeError), raised
// before any boundary curve reaches the stock's arrangement.
class LineSweepRegionTopologyError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// Exact: is |along| rational, i.e. is dx^2 + dy^2 the square of a rational?
// Decided on the reduced fraction p/q of dx^2 + dy^2 through CGAL's number-type
// machinery (Fraction_traits::Decompose + Algebraic_structure_traits::Is_square
// on numerator and denominator); never through to_double.
bool is_pythagorean_direction(const EVector& along);

enum class LineSweepModel2 {
    DISK,
    EXACT_CAPSULE,
    CANONICAL_QUAD,
};

class LineSweepRegion2 {
public:
    // Validates at the double boundary (finite inputs, positive radius), decides
    // the model exactly, and certifies both rectangle bands exactly.
    static LineSweepRegion2 build(double x0, double y0, double x1, double y1, double radius);

    LineSweepModel2 model() const noexcept { return model_; }
    bool is_exact() const noexcept { return model_ != LineSweepModel2::CANONICAL_QUAD; }

    // Build the region IN PLACE: a Gps must never be returned by value (its
    // copy aliases the source's traits object; docs/exactness.md).
    void under_region_into(Gps& target) const;
    void over_region_into(Gps& target) const;

    // The UNDER region described for local depletion: its true boundary curves,
    // traversed with the surviving material on their left, and one strictly
    // interior seed (the exact midpoint of A and B; the region is connected).
    LocalRemovalRegion local_under_region() const;

private:
    LineSweepRegion2(EPoint start, EPoint end, Epeck::FT tool_radius,
                     LineSweepModel2 model, EVector under_half_width,
                     EVector over_half_width);

    EPoint start_;
    EPoint end_;
    Epeck::FT tool_radius_;
    LineSweepModel2 model_;
    EVector under_half_width_;
    EVector over_half_width_;
};
