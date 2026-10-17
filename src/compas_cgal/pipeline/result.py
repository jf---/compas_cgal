"""The pipeline's result and its forced-decision records (Ruling 1.1-R15).

Every decision the pipeline takes AGAINST its own rules is a record, never a
warning and never silent:

* `DeclinedRun` -- guide stations that carry no loop: below the degeneracy
  floor, or the gap of a chain split at an unbridgeable radius step;
* `ForcedAdvance` -- an advance taken although no candidate passed every
  condition (the minimum one keeping the radius-step bound);
* `ForcedEntry` -- an entry from the clearance plane: the first cut (a helix,
  `EntryKind.HELIX`), a chain no admissible link reaches, or the re-entry after
  a split (plunges).

The ledger is the closed union of these three, built and checked in one place
(`ForcedDecisionLedger.build`). `OmittedLoop` (a loop not emitted because it
would remove nothing) and `LinkRecord` (a cut-plane link the walk emitted) are
reported beside the ledger, not in it: they are the rules working, not forced.
`RampOutRecord` (Task 1.2, brief N2) locates the first cut's helix, its pencil
members and its root loop in the operation stream, checked against it.

Every station index is in guide emission order.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from enum import Enum
from typing import Self

from compas.geometry import Circle

from compas_cgal.helix import Helix
from compas_cgal.pipeline.errors import LedgerInvariantError
from compas_cgal.pipeline.errors import RampOutRecordMismatchError
from compas_cgal.stock import SweepRemoval
from compas_cgal.toolpath import OperationType
from compas_cgal.toolpath import ToolpathOperation
from compas_cgal.toolpath import ToolpathResult
from compas_cgal.units import Point2
from compas_cgal.units import WorldXY


class DeclineCause(Enum):
    """Why a run of guide stations carries no loop."""

    FLOOR = "floor"
    RADIUS_STEP = "radius_step"


@dataclass(frozen=True)
class DeclinedRun:
    """A maximal run of stations on one chain that the walk leaves uncut.

    For `DeclineCause.FLOOR` the run is the stations below the degeneracy
    floor. For `DeclineCause.RADIUS_STEP` it is the closed range from the
    station the walk stopped at to the one it re-entered at: both ends are
    machined, the material between their annuli is not.

    Attributes:
        chain: The chain's index in the guide's emission order.
        first_station: First station of the run, in guide emission order.
        last_station: Last station of the run (inclusive), in guide emission order.
        cause: Why the run is declined.
    """

    chain: int
    first_station: int
    last_station: int
    cause: DeclineCause = DeclineCause.FLOOR


@dataclass(frozen=True)
class ForcedAdvance:
    """An advance taken past refusing conditions: the minimum one keeping the radius-step bound.

    Attributes:
        chain: The chain's index in the guide's emission order.
        reversed_walk: Whether the chain was walked against the guide's order.
        origin_station: Station advanced from, in guide emission order.
        station: Station advanced to, in guide emission order.
    """

    chain: int
    reversed_walk: bool
    origin_station: int
    station: int


class ForcedEntryReason(Enum):
    """Why the tool entered the cut plane from the clearance plane."""

    FIRST_CUT = "first_cut"
    NO_ADMISSIBLE_LINK = "no_admissible_link"
    RADIUS_STEP_SPLIT = "radius_step_split"


class EntryKind(Enum):
    """How the tool descended from the clearance plane to the cut plane."""

    PLUNGE = "plunge"
    HELIX = "helix"


@dataclass(frozen=True)
class ForcedEntry:
    """A plunge from the clearance plane onto a loop's exact entry.

    Attributes:
        path_index: Path index of the operations the entry begins.
        entry: The plunge foot, the loop's exact entry.
        reason: Why the entry was forced.
        into_material: Whether the plunge disk removes material, read from
            `Stock.circle_sweep_removal(x, y, 0, r)`, the gate's own plunge
            predicate.
    """

    path_index: int
    entry: Point2[WorldXY]
    reason: ForcedEntryReason
    into_material: bool
    entry_kind: EntryKind = EntryKind.PLUNGE


ForcedDecision = ForcedEntry | DeclinedRun | ForcedAdvance


@dataclass(frozen=True)
class ForcedDecisionLedger:
    """Every forced decision of one generated path, in the order taken.

    Attributes:
        records: The records.
    """

    records: tuple[ForcedDecision, ...]

    @classmethod
    def build(cls, records: Sequence[ForcedDecision]) -> Self:
        """Check the ledger invariant and freeze the records.

        Args:
            records: The records, in the order taken.

        Returns:
            The ledger.

        Raises:
            LedgerInvariantError: The records do not hold exactly one
                `ForcedEntryReason.FIRST_CUT` entry.
        """
        first_cuts = sum(1 for record in records if isinstance(record, ForcedEntry) and record.reason is ForcedEntryReason.FIRST_CUT)
        if first_cuts != 1:
            raise LedgerInvariantError(f"a generated path has exactly one FIRST_CUT entry; the ledger holds {first_cuts}.")
        return cls(records=tuple(records))

    @property
    def entries(self) -> tuple[ForcedEntry, ...]:
        """The forced entries."""
        return tuple(record for record in self.records if isinstance(record, ForcedEntry))

    @property
    def declined(self) -> tuple[DeclinedRun, ...]:
        """The declined runs."""
        return tuple(record for record in self.records if isinstance(record, DeclinedRun))

    @property
    def advances(self) -> tuple[ForcedAdvance, ...]:
        """The forced advances."""
        return tuple(record for record in self.records if isinstance(record, ForcedAdvance))


@dataclass(frozen=True)
class OmittedLoop:
    """A loop the walk did not emit because its exact sweep would remove nothing.

    Attributes:
        path_index: Path index of the chain walk it was skipped in.
        chain: The chain's index in the guide's emission order.
        station: The station, in guide emission order.
        verdict: The exact removal bracket that decided the skip.
    """

    path_index: int
    chain: int
    station: int
    verdict: SweepRemoval


class LinkKind(Enum):
    """The shape of a cut-plane link between two loops."""

    SAME_CIRCLE_ARC = "same_circle_arc"
    BELT = "belt"
    STRAIGHT = "straight"


@dataclass(frozen=True)
class LinkRecord:
    """A cut-plane link the walk emitted between two loops (in-chain bridges included).

    Attributes:
        path_index: Path index of the link's operations.
        kind: The link's shape.
        first_operation: Index of the link's first operation in the stream.
        operation_count: How many operations the link emitted.
    """

    path_index: int
    kind: LinkKind
    first_operation: int
    operation_count: int


@dataclass(frozen=True)
class RampOutRecord:
    """Where the first cut's helical entry sits in the operation stream (brief N2).

    Attributes:
        path_index: Path index of the entry's operations.
        helix_operation: Index of the ``HELIX`` operation.
        members: Indices of the pencil members (CUT circles), contiguous after
            the helix; empty when the first scan admitted the root.
        root_operation: Index of the root loop, right after the members.
    """

    path_index: int
    helix_operation: int
    members: range
    root_operation: int

    @classmethod
    def build(cls, operations: Sequence[ToolpathOperation], path_index: int, helix_operation: int, members: range, root_operation: int) -> Self:
        """Check the record against *operations* and freeze it.

        The helix is a ``HELIX`` operation; the members are a contiguous run of
        CUT circles right after it, and the root loop a circle right after them;
        every member and the root are entered exactly at the helix's entry
        ``P``, on *path_index*.

        Args:
            operations: The emitted stream.
            path_index: Path index of the entry.
            helix_operation: Index of the helix.
            members: Indices of the pencil members.
            root_operation: Index of the root loop.

        Returns:
            The record.

        Raises:
            RampOutRecordMismatchError: Any of the above does not hold.
        """
        if not 0 <= helix_operation < len(operations) or not 0 <= root_operation < len(operations):
            raise RampOutRecordMismatchError(f"record indices {helix_operation}, {root_operation} are outside the stream of {len(operations)} operations.")
        helix_op = operations[helix_operation]
        if helix_op.operation is not OperationType.HELIX or not isinstance(helix_op.geometry, Helix) or helix_op.path_index != path_index:
            raise RampOutRecordMismatchError(f"operation {helix_operation} is not the entry's HELIX on path {path_index}.")
        if members.step != 1 or members.start != helix_operation + 1 or members.stop != root_operation:
            raise RampOutRecordMismatchError(f"members {members!r} are not the contiguous run from the helix ({helix_operation}) to the root ({root_operation}).")
        for index in [*members, root_operation]:
            op = operations[index]
            is_cut = op.operation is OperationType.CUT or index == root_operation
            if not (is_cut and isinstance(op.geometry, Circle) and op.entry_xy == helix_op.entry_xy and op.path_index == path_index):
                raise RampOutRecordMismatchError(f"operation {index} is not a loop entered at the helix's entry {helix_op.entry_xy!r} on path {path_index}.")
        return cls(path_index=path_index, helix_operation=helix_operation, members=members, root_operation=root_operation)


@dataclass
class PipelineResult(ToolpathResult):
    """The generated path, its forced-decision ledger, and what the rules left out.

    Attributes:
        ledger: The forced decisions.
        omitted: Loops not emitted because they would remove nothing.
        links: Cut-plane links between loops.
        ramp_out: Where the first cut's helix, pencil and root loop sit in the stream.
    """

    ledger: ForcedDecisionLedger
    omitted: tuple[OmittedLoop, ...]
    links: tuple[LinkRecord, ...]
    ramp_out: RampOutRecord
