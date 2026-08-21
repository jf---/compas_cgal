"""The one record every benchmark run emits.

Generation and certification are timed SEPARATELY because they are different
businesses: generation is already competitive with the published state of the
art, certification is the thing the state of the art does not do. A single
combined number hides exactly the fact that matters.

`to_dict` / `from_dict` are the corpus's serialisation contract. `from_dict`
refuses any payload whose columns are not exactly the schema's, so a renamed or
dropped field fails at the seam instead of silently producing a record whose
missing measurement reads as a zero.

The cap is recorded as TWO columns because "over the cap" and "not proved under
the cap" are different measurements pointing in opposite directions, and a single
column named after either one misreports the other. `uncertified` over-counts by
construction, `truly_exceeding` under-counts by construction; see
`benchmarks.exceedance` for why neither substitutes for the other.
"""

from __future__ import annotations

from dataclasses import MISSING
from dataclasses import asdict
from dataclasses import dataclass
from dataclasses import fields
from typing import Any
from typing import Mapping

from benchmarks.errors import MalformedRecordError


@dataclass(frozen=True)
class MeasurementRecord:
    """One instance's measured outcome.

    Attributes:
        name: Instance name.
        family: Sweep family.
        params: Sweep coordinates, for plotting.
        tool_diameter: Cutter diameter used.
        tea_cap_deg: Engagement cap used.
        generate_seconds: Wall time of toolpath generation.
        certify_seconds: Wall time of the engagement audit.
        operations: Total toolpath operations.
        cut_operations: Operations that engaged material.
        stations: Total certifier stations, the refinement-depth cost proxy.
        max_tea_deg: Worst observed engagement angle in degrees.
        uncertified: Operations whose cap could NOT BE PROVED. Sound and
            conservative: an operation counts here whenever the certificate does
            not close, including when it was never measured because the growth
            guard could not close at that station density. It over-counts genuine
            violations and must never be reported as a violation count.
        truly_exceeding: Cut motions with at least one sampled cutter position
            where the exact predicate reports the cap exceeded. A sampled LOWER
            BOUND on true exceedance, never a certificate
            (`benchmarks.exceedance`).
        unresolved: Operations the certifier could not decide either way.
        arrangement_vertices_final: Arrangement vertices after the toolpath has
            been replayed onto a fresh stock.
        max_coordinate_digits: Longest exact rational after that replay; 0 when
            digit collection was disabled.
        error: Exception text when the instance failed, else None.
    """

    name: str
    family: str
    params: Mapping[str, float]
    tool_diameter: float
    tea_cap_deg: float
    generate_seconds: float
    certify_seconds: float
    operations: int
    cut_operations: int
    stations: int
    max_tea_deg: float
    uncertified: int
    truly_exceeding: int
    unresolved: int
    arrangement_vertices_final: int
    max_coordinate_digits: int
    error: str | None = None

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serialisable mapping.

        Returns:
            The record as a plain dict.
        """
        return asdict(self)

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "MeasurementRecord":
        """Rebuild a record from its serialised form, rejecting a mismatched schema.

        Args:
            data: A mapping produced by `to_dict`. Columns with a declared
                default may be omitted; every other column must be present, and
                no unknown column is accepted.

        Returns:
            The record.

        Raises:
            MalformedRecordError: The payload has unknown or missing columns.
        """
        declared = {f.name for f in fields(cls)}
        required = {f.name for f in fields(cls) if f.default is MISSING and f.default_factory is MISSING}
        present = set(data)
        unknown = sorted(present - declared)
        missing = sorted(required - present)
        if unknown or missing:
            raise MalformedRecordError(f"MeasurementRecord payload has unknown columns {unknown} and is missing {missing}.")
        return cls(**dict(data))

    @classmethod
    def failed(cls, name: str, family: str, params: Mapping[str, float], tool_diameter: float, tea_cap_deg: float, error: str) -> "MeasurementRecord":
        """Build the record for an instance that could not be measured.

        Args:
            name: Instance name.
            family: Sweep family.
            params: Sweep coordinates.
            tool_diameter: Cutter diameter.
            tea_cap_deg: Engagement cap.
            error: Exception text.

        Returns:
            A record with zeroed measurements and the error recorded.
        """
        return cls(
            name=name,
            family=family,
            params=dict(params),
            tool_diameter=tool_diameter,
            tea_cap_deg=tea_cap_deg,
            generate_seconds=0.0,
            certify_seconds=0.0,
            operations=0,
            cut_operations=0,
            stations=0,
            max_tea_deg=0.0,
            uncertified=0,
            truly_exceeding=0,
            unresolved=0,
            arrangement_vertices_final=0,
            max_coordinate_digits=0,
            error=error,
        )
