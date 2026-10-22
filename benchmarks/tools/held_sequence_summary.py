"""The numbers a sequence figure prints, written beside it as JSON so docs quote the generator, not a hand copy.

Every value is a reported double or count from the run that drew the figure;
no decision is made here. ``held_sequence_gallery_values`` reads these files
to fill the marked numbers on the gallery pages.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from benchmarks.errors import BenchmarkError
from benchmarks.tools.held_sequence_plot import CoverageVerdict

SUMMARY_VERSION = 2
SUMMARY_SUFFIX = ".summary.json"


class InvalidFigureSummaryError(BenchmarkError):
    """A figure summary contradicts itself or does not carry the version this module reads."""


@dataclass(frozen=True)
class CoverageSummary:
    """The certificate verdict as counts: exact witnesses and the display sites they group into."""

    covered: bool
    witnesses: int
    sites: int

    @classmethod
    def build(cls, *, covered: bool, witnesses: int, sites: int) -> CoverageSummary:
        """Validate the counts against the verdict.

        Raises:
            InvalidFigureSummaryError: a covered verdict with a witness or site,
                an uncovered one without, or more sites than witnesses.
        """
        if covered != (witnesses == 0) or covered != (sites == 0) or sites > witnesses or witnesses < 0:
            raise InvalidFigureSummaryError(f"covered={covered} with {witnesses} witnesses at {sites} sites.")
        return cls(covered, witnesses, sites)

    @classmethod
    def from_verdict(cls, verdict: CoverageVerdict, *, join_mm: float) -> CoverageSummary:
        """Count ``verdict``'s witnesses and its sites grouped at ``join_mm``, as the plot labels them."""
        return cls.build(covered=verdict.covered, witnesses=len(verdict.witnesses), sites=len(verdict.sites(join_mm)))

    def to_dict(self) -> dict[str, Any]:
        return {"covered": self.covered, "witnesses": self.witnesses, "sites": self.sites}

    @classmethod
    def from_dict(cls, raw: dict[str, Any]) -> CoverageSummary:
        return cls.build(covered=bool(raw["covered"]), witnesses=int(raw["witnesses"]), sites=int(raw["sites"]))


@dataclass(frozen=True)
class FigureSummary:
    """The metrics line and verdict of one emitted-path figure.

    ``over_cap`` and ``unresolved`` count successor pairs only; ``over_cap``
    and ``max_known_deg`` are binary64 reports of the predecessor model
    (``angle > cap`` on doubles), not the exact engagement gate, so at the cap
    they may disagree with it; ``seam_deg`` is
    the reported seam pair (None: open path; NaN: unresolved). ``coverage`` is
    None when no certificate was run. ``inserted_sources`` are the source
    indices a certificate-gated coverage repair re-inserted (None: no repair
    ran; empty: the repair needed none).
    """

    figure: str
    circles: int
    connectors: int
    over_cap: int
    unresolved: int
    max_known_deg: float
    seam_deg: float | None
    coverage: CoverageSummary | None
    inserted_sources: tuple[int, ...] | None = None

    @classmethod
    def build(
        cls,
        *,
        figure: str,
        circles: int,
        connectors: int,
        over_cap: int,
        unresolved: int,
        max_known_deg: float,
        seam_deg: float | None,
        coverage: CoverageSummary | None,
        inserted_sources: tuple[int, ...] | None = None,
    ) -> FigureSummary:
        """Validate the counts.

        Raises:
            InvalidFigureSummaryError: no circle, a connector count that is
                neither open (n - 1) nor closed (n), a negative count, or
                inserted sources that repeat or are negative.
        """
        if circles < 1 or connectors not in (circles - 1, circles) or over_cap < 0 or unresolved < 0:
            raise InvalidFigureSummaryError(f"{figure}: {circles} circles, {connectors} connectors, {over_cap} over cap, {unresolved} unresolved.")
        if inserted_sources is not None and (len(set(inserted_sources)) != len(inserted_sources) or any(source < 0 for source in inserted_sources)):
            raise InvalidFigureSummaryError(f"{figure}: inserted sources {inserted_sources} repeat or are negative.")
        return cls(figure, circles, connectors, over_cap, unresolved, max_known_deg, seam_deg, coverage, inserted_sources)

    def to_json(self) -> str:
        seam = None if self.seam_deg is None or math.isnan(self.seam_deg) else self.seam_deg
        return (
            json.dumps(
                {
                    "version": SUMMARY_VERSION,
                    "figure": self.figure,
                    "circles": self.circles,
                    "connectors": self.connectors,
                    "over_cap": self.over_cap,
                    "unresolved": self.unresolved,
                    "max_known_deg": self.max_known_deg,
                    "seam_deg": seam,
                    "seam_unresolved": self.seam_deg is not None and math.isnan(self.seam_deg),
                    "coverage": None if self.coverage is None else self.coverage.to_dict(),
                    "inserted_sources": None if self.inserted_sources is None else list(self.inserted_sources),
                },
                indent=2,
                allow_nan=False,
            )
            + "\n"
        )

    @classmethod
    def from_json(cls, text: str) -> FigureSummary:
        """Parse a summary written by :meth:`to_json`.

        Raises:
            InvalidFigureSummaryError: another version, or inconsistent counts.
        """
        raw = json.loads(text)
        if raw.get("version") != SUMMARY_VERSION:
            raise InvalidFigureSummaryError(f"summary version {raw.get('version')!r} is not {SUMMARY_VERSION}.")
        seam = math.nan if raw["seam_unresolved"] else (None if raw["seam_deg"] is None else float(raw["seam_deg"]))
        coverage = None if raw["coverage"] is None else CoverageSummary.from_dict(raw["coverage"])
        return cls.build(
            figure=str(raw["figure"]),
            circles=int(raw["circles"]),
            connectors=int(raw["connectors"]),
            over_cap=int(raw["over_cap"]),
            unresolved=int(raw["unresolved"]),
            max_known_deg=float(raw["max_known_deg"]),
            seam_deg=seam,
            coverage=coverage,
            inserted_sources=None if raw["inserted_sources"] is None else tuple(int(source) for source in raw["inserted_sources"]),
        )


def summary_path(figure: Path) -> Path:
    """Where the summary of the figure at ``figure`` lives: beside it, ``<stem>.summary.json``."""
    return figure.with_suffix(SUMMARY_SUFFIX)
