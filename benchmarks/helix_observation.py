"""The survey's observation of a HELIX entry (Task 1.2): ends, bore, and measured exit.

A helix is not cap-sampled (E8): a bore into virgin stock is full immersion by
construction, scored by its ramp. `benchmarks.survey.survey_path` therefore
records it as a `HelixObservation`, never a `MotionQuality`, and the motion
after it is sampled as usual.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING
from typing import Tuple

from benchmarks.depletion import HelixEvent
from benchmarks.survey_station import EngagementSample
from benchmarks.survey_station import perpendicular
from benchmarks.survey_station import seam_station
from benchmarks.survey_station import swept_area
from compas_cgal import _coverage_2
from compas_cgal.helix import Helix
from compas_cgal.units import ExactMillimetre

if TYPE_CHECKING:
    from benchmarks.models import MachineModel


@dataclass(frozen=True)
class HelixObservation:
    """One helical entry: its bore, its ends, and the station it hands over at.

    A helix is not cap-sampled (E8): a bore into virgin stock is full
    immersion by construction, and it is scored by its ramp, not by
    engagement. So it is not a `MotionQuality`; what the junction and
    continuity reductions need of it is recorded here.

    Attributes:
        helix: The exact helix, for the ramp verdict (`ramp_within_limit`).
        index: Position of the operation in its toolpath.
        into_material: Whether the bore met the stock's interior. EXACT
            (`benchmarks.depletion.HelixEvent`).
        turns: Full turns of the descent.
        pitch: Descent per turn. EXACT.
        top_z: Start height. REPORTING.
        bottom_z: End height, the cut plane. REPORTING.
        length: Path length (`compas_cgal.helix.Helix.length`). REPORTING.
        curvature: Plan-view curvature ``1 / rho_h`` of the footprint, the
            cut-plane model's curvature of a loop on it (`MotionQuality.curvature`).
            REPORTING; read by the curvature-break check at the helix exit
            and by the feed-limited time (`benchmarks.quality.machine_outcome`).
        swept_area: The bore's area in closed form (`swept_area` of the
            footprint circle): the disk ``pi (rho_h + r)^2`` when ``rho_h <= r``,
            the annulus ``4 pi rho_h r`` when ``rho_h > r`` -- the region the
            replay subtracts (`benchmarks.depletion`).
        start: Cutter-centre XY at the start: the exact entry ``P`` as doubles.
        end: Cutter-centre XY at the end: ``P`` again (integer turns).
        start_tangent: Unit XY travel direction at the start, in the helix's
            turn sense.
        end_tangent: Unit XY travel direction at the end.
        exit_sample: The junction station (gate rule f): the footprint circle's
            seam, queried against the stock the helix started from with the
            FULL turn as prefix, i.e. against the bored stock (amendment M6).
            MEASURED, never assumed zero. It enters no per-motion statistic.
    """

    helix: Helix
    index: int
    into_material: bool
    turns: int
    pitch: ExactMillimetre
    top_z: float
    bottom_z: float
    length: float
    curvature: float
    swept_area: float
    start: Tuple[float, float]
    end: Tuple[float, float]
    start_tangent: Tuple[float, float]
    end_tangent: Tuple[float, float]
    exit_sample: EngagementSample

    def ramp_within_limit(self, machine: MachineModel) -> bool:
        """Whether the helix ramps within *machine*'s ramp cap: EXACT.

        Decided by `compas_cgal.helix.Helix.ramp_within` against the machine's
        exact surrogate `MachineModel.helix_pitch_per_radius`
        (``pitch <= exact(K) * rho_h`` over the rationals), never by a ramp
        angle in doubles.

        Args:
            machine: The machine whose ramp cap applies.

        Returns:
            True iff the helix is within the cap.
        """
        return self.helix.ramp_within(machine.helix_pitch_per_radius)


def measure_helix(event: HelixEvent, tool_radius: float, cap_ratio: float, centre_domain: "_coverage_2.CutterCentreDomain2") -> HelixObservation:
    """Observe one helix: its ends, its bore, and its measured exit station.

    Not cap-sampled (E8). The exit is the seam query of the survey's loop seam run
    on the footprint circle (`compas_cgal.helix.Helix.footprint_circle`) at the
    exact exit ``P``, against the stock the helix started from with the full
    turn as prefix: the bored stock the next motion meets (amendment M6).

    Args:
        event: The helix, its exact into-material verdict, and the pre-bore stock.
        tool_radius: Tool radius.
        cap_ratio: Exact squared-chord surrogate of the cap.
        centre_domain: Where a cutter of this radius may legally be centred.

    Returns:
        The observation.
    """
    helix = event.helix
    start, end = helix_endpoints(helix)
    start_tangent, end_tangent = helix_tangents(helix)
    length = float(helix.length())
    footprint = helix.footprint_circle()
    return HelixObservation(
        helix=helix,
        index=event.index,
        into_material=event.into_material,
        turns=helix.turns,
        pitch=helix.pitch,
        top_z=float(helix.top_z.value),
        bottom_z=float(helix.bottom_z.value),
        length=length,
        curvature=1.0 / float(footprint.radius),
        swept_area=swept_area(footprint, tool_radius),
        start=start,
        end=end,
        start_tangent=start_tangent,
        end_tangent=end_tangent,
        exit_sample=seam_station(event.stock_before, footprint, end, helix.clockwise, length, tool_radius, cap_ratio, centre_domain),
    )


def helix_endpoints(helix: Helix) -> Tuple[Tuple[float, float], Tuple[float, float]]:
    """Where a helix starts and ends in XY: its exact entry ``P``, both times.

    Args:
        helix: The helix.

    Returns:
        ``(start, end)``; equal, because the helix turns a whole number of times.
    """
    return (float(helix.entry_xy.x), float(helix.entry_xy.y)), (float(helix.exit_xy.x), float(helix.exit_xy.y))


def helix_tangents(helix: Helix) -> Tuple[Tuple[float, float], Tuple[float, float]]:
    """Unit XY travel directions where a helix starts and ends, in its turn sense.

    The plan-view tangent of the footprint circle at ``P``; the descent is
    along Z, which the cut-plane model drops.

    Args:
        helix: The helix.

    Returns:
        ``(start_tangent, end_tangent)``.
    """
    cx, cy = float(helix.footprint.centre.x), float(helix.footprint.centre.y)
    start, end = helix_endpoints(helix)
    turn = -1.0 if helix.clockwise else 1.0
    return perpendicular(start[0] - cx, start[1] - cy, turn), perpendicular(end[0] - cx, end[1] - cy, turn)
