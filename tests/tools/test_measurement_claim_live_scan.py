from __future__ import annotations

import pathlib
import re


PROJECT_ROOT = pathlib.Path(__file__).resolve().parents[2]
EXTRACTION_PATTERN = re.compile(r"# .*(MEASURED|[Mm]easured (on|at|against)|measures [0-9]|[0-9]+(\.[0-9]+)?[x×] (faster|slower))")
EXPECTED_RESIDUAL_HITS = (
    (
        "src/compas_cgal/engagement_radial_toolpath.py",
        "# MEASURED, on the 20x12 pocket at a 60 deg cap with the gate below in place:",
    ),
)


def test_frozen_extractor_finds_only_the_expected_live_residual_claim() -> None:
    hits: list[tuple[str, str]] = []
    for root_name in ("src/compas_cgal", "benchmarks"):
        root = PROJECT_ROOT / root_name
        for path in sorted(root.rglob("*.py")):
            relative = path.relative_to(PROJECT_ROOT).as_posix()
            for line in path.read_text(encoding="utf-8").splitlines():
                if "superseded" not in line and EXTRACTION_PATTERN.search(line) is not None:
                    hits.append((relative, line))
    assert tuple(hits) == EXPECTED_RESIDUAL_HITS
