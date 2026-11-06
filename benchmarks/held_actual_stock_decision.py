"""The actual-stock engagement decision of the refinement: one target on its circle view (design 2026-11-04, section 2.2).

The view C_t is the part polygon P minus the closed swept annuli of every
circle emitted before t, connectors left out. C_t contains the actual stock
S_t, and the 45-station decision is monotone under stock inclusion
(``src/swept_prefix_engagement_2.h``), so a pass here is a pass on S_t. The
engine is ``AnnulusView2``, its filtered path: ``circle_exceeded_filtered``
equals the exact verdict on the whole view.

One instance belongs to one emission: it appends what the emission appends, so
the view is the emitted list, except during a splice's trial of an adopted tail,
which it cuts back. Circle 0 has no view (ruling R3): it is the
``EntryCircle``, whose engagement the entry move owns, and no decision is ever
asked of it.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from benchmarks.errors import BenchmarkError
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_engagement_trace import ReadDisk
from benchmarks.held_material_domain import PartPolygon
from compas_cgal import _stock_2
from compas_cgal.motion import EngagementCap
from compas_cgal.pipeline.admissibility import loop_parameters
from compas_cgal.units import ToolRadius

# The certificate's decision: no gap closure (design section 2.2).
NO_GAP_CLOSURE = 0.0
# The emission index of the entry circle.
ENTRY_POSITION = 0


class ActualStockDecisionError(BenchmarkError):
    """A caller breached the stock decision's contract; never a refusal of a circle."""


class EntryCircleDecisionError(ActualStockDecisionError):
    """A decision was asked of emission index 0: its stock is P, and the entry move owns its engagement."""


class ActualStockViewMismatchError(ActualStockDecisionError):
    """A decision was asked at a position other than the next one: the view is out of step with the emission."""


@dataclass(frozen=True)
class EntryCircle:
    """Emission index 0, exempt by name: never decided, and in every later view."""

    circle: Figure5CounterclockwiseCircle


def read_disk(circle: Figure5CounterclockwiseCircle) -> ReadDisk:
    """The doubles a decision on ``circle`` injects, and its footprint."""
    return ReadDisk(float(circle.center.x), float(circle.center.y), float(circle.radius.value), float(circle.contact_point.x), float(circle.contact_point.y))


class ActualStockDecision:
    """The stock decision of one emission; the view holds exactly the emitted circles."""

    def __init__(self, part: PartPolygon, view: _stock_2.AnnulusView2, chord_ratio: float) -> None:
        self.part = part
        self._view = view
        self._chord_ratio = chord_ratio
        self._parameters = loop_parameters()

    @classmethod
    def open(cls, part: PartPolygon, tool_radius: ToolRadius, cap: EngagementCap, entry: EntryCircle, decided: Sequence[Figure5CounterclockwiseCircle]) -> ActualStockDecision:
        """A decision on P whose view holds ``entry`` and then ``decided``.

        The cap's ``chord_ratio`` is the native ``cap_chord_ratio`` the certificate uses.

        Raises:
            InvalidAnnulusViewPolygonError: P is not simple (from the engine).
        """
        decision = cls(part, _stock_2.AnnulusView2.open(part.stock_rows(), float(tool_radius.value)), cap.chord_ratio)
        decision.append(entry.circle)
        for circle in decided:
            decision.append(circle)
        return decision

    @property
    def size(self) -> int:
        """Circles in the view: the emission index of the next target."""
        return int(self._view.size)

    @property
    def counters(self) -> _stock_2.AnnulusFilterCounters2:
        """The engine's count of this view's decisions, a snapshot."""
        return self._view.filter_counters

    def append(self, circle: Figure5CounterclockwiseCircle) -> None:
        """The emission appended ``circle``."""
        self._view.append(float(circle.center.x), float(circle.center.y), float(circle.radius.value))

    def truncate(self, size: int) -> None:
        """Drop the circles from emission index ``size`` on, after a trial of an adopted tail."""
        self._view.truncate(size)

    def exceeded(self, target: ReadDisk, position: int) -> bool:
        """Exact: some station of ``target`` exceeds the cap on C_position.

        Raises:
            EntryCircleDecisionError: ``position`` is the entry circle's.
            ActualStockViewMismatchError: ``position`` is not the view's size.
        """
        if position == ENTRY_POSITION:
            raise EntryCircleDecisionError("the entry circle is exempt by name; its engagement is owned by the entry move")
        if position != self._view.size:
            raise ActualStockViewMismatchError(f"decision at position {position} on a view of {self._view.size} circles")
        return bool(
            self._view.circle_exceeded_filtered(0, position, target.cx, target.cy, target.rho, target.ex, target.ey, False, self._parameters, self._chord_ratio, NO_GAP_CLOSURE)
        )
