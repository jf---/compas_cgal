#pragma once

#include "exact/rational.h"
#include "partition_certificate.h"

#include <string>

/// A rational in the station lane's canonical attestation form.
///
/// DERIVED VIEW, never a carrier. Since stage 2 the station source computes on
/// `compas_cgal::exact::Rational` and this type only renders those values for
/// the replay digest. Its `exact-rational-v1` framing is a frozen compatibility
/// contract, anchored to a byte literal in `exact_station_attestation_gate`.
class ExactRational2 {
public:
    /// Project an exact rational into the station lane's attestation form.
    ///
    /// The station lane's single exact-to-attestation-bytes projection site and
    /// its only caller of `compas_cgal::exact::to_canonical`.
    ///
    /// Args:
    ///     value: any exact rational.
    ///
    /// Returns:
    ///     Its reduced, positive-denominator canonical form.
    ///
    /// Raises:
    ///     compas_cgal::exact::UnreducedCanonicalRationalError: if the
    ///         decomposed denominator is not positive.
    static ExactRational2 project(
        const compas_cgal::exact::Rational& value);

    /// Decode decimal text into attestation form.
    ///
    /// LEGACY DECODER. Nothing in the pipeline calls it since stage 2; it is
    /// retained so `exact_station_attestation_gate` can keep proving `project`
    /// agrees with it byte for byte. Removing it is stage 6 and requires
    /// explicit user permission.
    ///
    /// Raises:
    ///     InvalidStationSourceError: if `text` is not one exact rational, or
    ///         has a zero denominator.
    static ExactRational2 build(
        const std::string& text);

    const std::string& numerator() const noexcept;
    const std::string& denominator() const noexcept;

    /// Decimal text, "n" when the denominator is 1 and "n/d" otherwise.
    std::string text() const;

    /// The frozen `exact-rational-v1` record. Changing this encoding
    /// invalidates every stored full-circle replay digest.
    std::string canonical_bytes() const;

private:
    ExactRational2(
        std::string numerator,
        std::string denominator);

    std::string numerator_;
    std::string denominator_;
};

/// The station source's derived attestation view.
///
/// Built once, in the source's constructor, so the frozen bytes are produced at
/// exactly one site and no consumer ever re-projects them. It is a separate
/// struct from the source because the two have different consumers and
/// different lifetimes: the source's fields feed predicates on every cell, the
/// attestation feeds the digest once.
struct StationAttestation2 {
    ExactRational2 center_x;
    ExactRational2 center_y;
    ExactRational2 tool_radius;
    ExactRational2 cap_chord_ratio;

    /// The frozen `station-event-source-v1` record over the four values above.
    std::string canonical_bytes;
};

class StationEventSource2 {
public:
    /// Build a station from exact values.
    ///
    /// Args:
    ///     center_x: the station centre's x coordinate.
    ///     center_y: the station centre's y coordinate.
    ///     tool_radius: the cutter radius; must be strictly positive.
    ///     cap_chord_ratio: the engagement cap surrogate; must lie in (0, 4].
    ///
    /// Raises:
    ///     InvalidStationSourceError: if either range constraint is violated.
    static StationEventSource2 build(
        const compas_cgal::exact::Rational& center_x,
        const compas_cgal::exact::Rational& center_y,
        const compas_cgal::exact::Rational& tool_radius,
        const compas_cgal::exact::Rational& cap_chord_ratio);

    /// CANONICAL STATE. Every consumer computes on these. No parsing anywhere.
    const compas_cgal::exact::Rational& center_x() const noexcept;
    const compas_cgal::exact::Rational& center_y() const noexcept;
    const compas_cgal::exact::Rational& tool_radius() const noexcept;
    const compas_cgal::exact::Rational& cap_chord_ratio() const noexcept;

    /// DERIVED VIEW of the canonical state, projected once at construction.
    const StationAttestation2& attestation() const noexcept;

    /// The frozen `station-event-source-v1` record.
    ///
    /// Surfaced on the source, not only on the attestation, because
    /// `circle_strata.cpp` embeds it in every full-circle cell decision record
    /// and does so once per parameter cell. It is the attestation's own bytes.
    const std::string& canonical_bytes() const noexcept;

private:
    StationEventSource2(
        compas_cgal::exact::Rational center_x,
        compas_cgal::exact::Rational center_y,
        compas_cgal::exact::Rational tool_radius,
        compas_cgal::exact::Rational cap_chord_ratio);

    compas_cgal::exact::Rational center_x_;
    compas_cgal::exact::Rational center_y_;
    compas_cgal::exact::Rational tool_radius_;
    compas_cgal::exact::Rational cap_chord_ratio_;
    StationAttestation2 attestation_;
};

class InvalidStationSourceError
    : public EventSubstrateError {
public:
    using EventSubstrateError::EventSubstrateError;
};
