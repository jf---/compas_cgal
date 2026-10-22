"""Fill the marked numbers on the Held gallery pages from the figures' generated JSON; never hand-typed.

A page marks a generated value as ``<!-- gv:KEY -->value<!-- /gv -->``; this
tool replaces ``value`` with the one the figure run wrote, read from the
``*.summary.json`` sidecars (``held_sequence_summary``) and the
coverage-preserving reports (``held_coverage_preserving_plot``) in the images
directory. ``--check`` changes nothing and fails if any page is stale.

Keys are ``<figure>.<field>``; the figures and fields are the tables below.
A marked key with no generated value fails loud: a page never keeps a stale
number because its figure was not rendered.
"""

from __future__ import annotations

import argparse
import json
import re
from collections.abc import Callable
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from benchmarks.errors import BenchmarkError
from benchmarks.tools.held_sequence_summary import CoverageSummary
from benchmarks.tools.held_sequence_summary import FigureSummary
from benchmarks.tools.held_sequence_summary import summary_path

MARKER = re.compile(r"<!-- gv:(?P<key>[a-z0-9_.]+) -->(?P<value>.*?)<!-- /gv -->")
DEFAULT_IMAGES = Path("docs/assets/images")
GALLERY_PAGES = (Path("docs/held_toolpath_gallery.md"), Path("docs/held_motion_coverage.md"), Path("docs/segment_site_mat.md"))
# Emitted-path figures, by the key prefix the pages use.
PATH_FIGURES = {
    "fig5_refined": "held_figure5_engagement_refined.png",
    "fig5_current": "held_figure5_toolpath_current.png",
    "upper_lane": "held_figure8_upper_standard_refined.png",
    "upper_boundary": "held_figure8_upper_standard_refined_boundary.png",
    "skis": "held_figure8_crossed_skis_standard_refined_boundary.png",
    "monstera": "held_figure8_monstera_standard_refined_boundary_qualified_contacts_corner_approaches.png",
    "monstera_repaired": "held_figure8_monstera_standard_refined_boundary_qualified_contacts_corner_approaches_certified_repair.png",
}
# Coverage-preserving comparison reports, by key prefix.
COVERAGE_PRESERVING_REPORTS = {
    "fig5_covpres": "held_figure5_coverage_preserving.json",
    "upper_covpres": "held_figure8_upper_coverage_preserving.json",
}


class GalleryValuesError(BenchmarkError):
    """Base error of the gallery value filler."""


class UnknownGalleryValueError(GalleryValuesError):
    """A page marks a key that no rendered figure provides."""


class StaleGalleryValuesError(GalleryValuesError):
    """``--check``: a page's marked values differ from the generated ones."""


def _count(value: int) -> str:
    return f"{value:,}"


def _listing(values: tuple[int, ...]) -> str:
    """``744``, ``744 and 1273``, ``1, 2 and 3``; ``none`` for an empty tuple."""
    words = [str(value) for value in values]
    if not words:
        return "none"
    return words[0] if len(words) == 1 else f"{', '.join(words[:-1])} and {words[-1]}"


def _verdict(coverage: CoverageSummary | None) -> str:
    if coverage is None:
        return "not certified"
    return "CERTIFIED COVERED" if coverage.covered else "NOT COVERED"


def _path_values(summary: FigureSummary) -> dict[str, str]:
    values = {
        "circles": _count(summary.circles),
        "connectors": _count(summary.connectors),
        "over_cap": _count(summary.over_cap),
        "unresolved": _count(summary.unresolved),
        "max_deg": f"{summary.max_known_deg:.2f}",
        "verdict": _verdict(summary.coverage),
    }
    if summary.seam_deg is not None:
        values["seam_deg"] = f"{summary.seam_deg:.2f}"
    if summary.inserted_sources is not None:
        values["inserted_sources"] = _listing(summary.inserted_sources)
        values["inserted_count"] = _count(len(summary.inserted_sources))
    if summary.coverage is not None:
        values["witnesses"] = _count(summary.coverage.witnesses)
        values["sites"] = _count(summary.coverage.sites)
    return values


def _report_coverage(raw: dict[str, Any], field: str) -> CoverageSummary | None:
    return None if raw.get(field) is None else CoverageSummary.from_dict(raw[field])


def _coverage_preserving_values(raw: dict[str, Any]) -> dict[str, str]:
    values = {
        "before_circles": _count(int(raw["before_circles"])),
        "circles": _count(int(raw["circles"])),
        "connectors": _count(int(raw["connectors"])),
        "before_length_mm": f"{float(raw['before_length_mm']):,.1f}",
        "length_mm": f"{float(raw['length_mm']):,.1f}",
        "length_change": f"{float(raw['length_change_percent']):+.1f}%",
        "max_bound": f"{float(raw['max_reported_bound_degrees']):.3f}",
        "before_verdict": _verdict(_report_coverage(raw, "before_coverage")),
        "verdict": _verdict(_report_coverage(raw, "after_coverage")),
    }
    if raw.get("dense_circles") is not None:
        values["dense_circles"] = _count(int(raw["dense_circles"]))
    return values


def _collect(images: Path, files: dict[str, str], read: Callable[[Path], dict[str, str]]) -> dict[str, str]:
    values: dict[str, str] = {}
    for prefix, name in files.items():
        source = images / name
        if source.exists():
            values.update({f"{prefix}.{field}": value for field, value in read(source).items()})
    return values


def gallery_values(images: Path) -> dict[str, str]:
    """Every generated value available in ``images``, keyed ``<figure>.<field>``; figures not rendered are absent."""
    summaries = {prefix: summary_path(Path(name)).name for prefix, name in PATH_FIGURES.items()}
    values = _collect(images, summaries, lambda source: _path_values(FigureSummary.from_json(source.read_text())))
    values.update(_collect(images, COVERAGE_PRESERVING_REPORTS, lambda source: _coverage_preserving_values(json.loads(source.read_text()))))
    return values


def rewrite_marked_values(text: str, values: dict[str, str]) -> str:
    """``text`` with every marked value replaced by its generated one.

    Raises:
        UnknownGalleryValueError: a marked key has no generated value.
    """
    missing = sorted({match["key"] for match in MARKER.finditer(text)} - values.keys())
    if missing:
        raise UnknownGalleryValueError(f"no generated value for {', '.join(missing)}; render the figure first.")
    return MARKER.sub(lambda match: f"<!-- gv:{match['key']} -->{values[match['key']]}<!-- /gv -->", text)


def update_pages(pages: Sequence[Path], images: Path, *, check: bool) -> tuple[Path, ...]:
    """Rewrite (or, with ``check``, verify) the marked values of ``pages``; returns the pages that changed.

    Raises:
        UnknownGalleryValueError: a page marks a key with no generated value.
        StaleGalleryValuesError: ``check`` and some page is stale.
    """
    values = gallery_values(images)
    stale = []
    for page in pages:
        text = page.read_text()
        updated = rewrite_marked_values(text, values)
        if updated != text:
            stale.append(page)
            if not check:
                page.write_text(updated)
    if check and stale:
        raise StaleGalleryValuesError(f"stale generated values in {', '.join(map(str, stale))}; run held_sequence_gallery_values.")
    return tuple(stale)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--images", type=Path, default=DEFAULT_IMAGES, help="Directory holding the figures' summary and report JSON.")
    parser.add_argument("--check", action="store_true", help="Change nothing; fail if a page is stale.")
    parser.add_argument("pages", nargs="*", type=Path, default=list(GALLERY_PAGES))
    args = parser.parse_args()
    for page in update_pages(args.pages, args.images, check=args.check):
        print(f"updated {page}")


if __name__ == "__main__":
    main()
