"""Canonical row-level adjudication for generator measurement claims."""

from __future__ import annotations

from typing import TYPE_CHECKING
from typing import Literal

from tools.measurement_claim_units import Degrees
from tools.measurement_claim_units import Millimetres
from tools.measurement_claim_units import WorldMillimetres

if TYPE_CHECKING:
    from tools.measurement_claim_result import MC001EvidencePayload
    from tools.measurement_claim_result import MC002EvidencePayload
    from tools.measurement_claim_result import MC003EvidencePayload

StationDisposition = Literal["re-earned", "corrected"]

# Exact frozen-row baselines: centre/radii in mm, peaks in degrees, and audited Boolean material outcomes.
FROZEN_STATION_OCCURRENCE_COUNT = 1
FROZEN_STATION_CENTRE = (WorldMillimetres(18.482), WorldMillimetres(10.482))
FROZEN_MAXIMAL_RADIUS = Millimetres(0.5156)
FROZEN_COARSE_STEP = Millimetres(0.05)
FROZEN_RUNG_6_RADIUS = Millimetres(0.2156)
FROZEN_RUNG_6_PEAK = Degrees(61.3)
FROZEN_RUNG_6_CUTS_MATERIAL = True
FROZEN_RUNG_7_RADIUS = Millimetres(0.1656)
FROZEN_RUNG_7_PEAK = Degrees(5.9)
FROZEN_RUNG_7_CUTS_MATERIAL = False
FROZEN_REFINED_BAND_MIN_RADIUS = Millimetres(0.1719)
FROZEN_REFINED_BAND_MAX_RADIUS = Millimetres(0.2123)
FROZEN_FORCED_PEAK = Degrees(107.0)
FROZEN_RESCUED_PEAK = Degrees(59.0)


def adjudicate_station_claims(
    evidence_1: MC001EvidencePayload,
    evidence_2: MC002EvidencePayload,
    evidence_3: MC003EvidencePayload,
) -> tuple[tuple[StationDisposition, str], tuple[StationDisposition, str], tuple[StationDisposition, str]]:
    """Adjudicate the three frozen station rows without case-level aggregation."""
    centre = None if evidence_1["station_centre"] is None else tuple(evidence_1["station_centre"])
    mc001_matches = (evidence_1["occurrence_count"], centre, evidence_1["maximal_radius"]) == (
        FROZEN_STATION_OCCURRENCE_COUNT,
        FROZEN_STATION_CENTRE,
        FROZEN_MAXIMAL_RADIUS,
    )
    mc002_matches = (
        evidence_2["coarse_step"],
        evidence_2["rung_6_radius"],
        evidence_2["rung_6_peak"],
        evidence_2["rung_6_cuts_material"],
        evidence_2["rung_7_radius"],
        evidence_2["rung_7_peak"],
        evidence_2["rung_7_cuts_material"],
    ) == (
        FROZEN_COARSE_STEP,
        FROZEN_RUNG_6_RADIUS,
        FROZEN_RUNG_6_PEAK,
        FROZEN_RUNG_6_CUTS_MATERIAL,
        FROZEN_RUNG_7_RADIUS,
        FROZEN_RUNG_7_PEAK,
        FROZEN_RUNG_7_CUTS_MATERIAL,
    )
    mc003_matches = (
        evidence_3["rung_7_radius"],
        evidence_3["rung_7_peak"],
        evidence_3["rung_7_cuts_material"],
        evidence_3["refined_band_min_radius"],
        evidence_3["refined_band_max_radius"],
        evidence_3["forced_peak"],
        evidence_3["rescued_peak"],
    ) == (
        FROZEN_RUNG_7_RADIUS,
        FROZEN_RUNG_7_PEAK,
        FROZEN_RUNG_7_CUTS_MATERIAL,
        FROZEN_REFINED_BAND_MIN_RADIUS,
        FROZEN_REFINED_BAND_MAX_RADIUS,
        FROZEN_FORCED_PEAK,
        FROZEN_RESCUED_PEAK,
    )
    disposition_1: StationDisposition = "re-earned" if mc001_matches else "corrected"
    disposition_2: StationDisposition = "re-earned" if mc002_matches else "corrected"
    disposition_3: StationDisposition = "re-earned" if mc003_matches else "corrected"
    reasons = (
        "frozen station occurrence matches" if mc001_matches else "authenticated station occurrence differs",
        "frozen coarse-rung assertions match" if mc002_matches else "authenticated coarse-rung assertions differ",
        "frozen rung-7 and refined-band assertions match" if mc003_matches else "authenticated rung-7 or refined-band assertions differ",
    )
    return (disposition_1, reasons[0]), (disposition_2, reasons[1]), (disposition_3, reasons[2])
