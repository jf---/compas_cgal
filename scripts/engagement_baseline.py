#!/usr/bin/env python
"""Baseline engagement audit of the existing trochoid generator (SP2 reference).

For each benchmark pocket this generates a toolpath with the UNREGULATED
``trochoidal_mat_toolpath_circular``, replays it through the sound engagement
certifier ``audit_toolpath_engagement``, and records the tool-engagement angle
(TEA) each pocket actually sees against a 120 deg cap. The existing generator
has no engagement control, so the numbers are deliberately unflattering: this is
the reference the SP2 regulated generator must beat.

Outputs (into ``--out``, default ``docs/benchmarks``):
  * ``engagement_baseline.md``   BLUF-first report (worst pocket, worst TEA, totals)
  * ``engagement_baseline.json`` machine-readable per-pocket numbers

The audit depletes an exact ``Stock`` per replayed cut, so per-pocket cost is
O(n^2) in the cut count (the scalloped stock grows more complex with every loop --
see docs/superpowers/state/sp1-gate-c-analysis.md). The full 7-pocket run is ~1 hour
and is NOT regenerated automatically: this committed harness IS the deliverable, and
producing the reference ``.md``/``.json`` is deferred to SP2 (run
``python scripts/engagement_baseline.py`` once the depletion cost is addressed).
``--quick`` audits a tiny synthetic smoke pocket (~8 s) for the harness's own test;
``--pockets NAME[,NAME]`` selects a benchmark subset -- the same generate -> audit ->
report code path, fewer pockets.

Benchmark polygon coordinates are copied literally from ``tests/test_toolpath.py``;
scripts do not import test modules.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from dataclasses import dataclass
from pathlib import Path

from compas.geometry import Polygon
from compas_cgal.engagement import audit_toolpath_engagement
from compas_cgal.toolpath import trochoidal_mat_toolpath_circular

# --- Audit constants (fixed across all pockets) ----------------------------- #

TEA_CAP_DEG = 120.0
TEA_CAP_RAD = math.radians(TEA_CAP_DEG)
CLEARANCE_Z = 2.0

# --- Benchmark pockets (literal coordinates from tests/test_toolpath.py) ----- #

SQUARE = Polygon([[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0]])
IRREGULAR = Polygon(
    [
        (-1.91, 3.59, 0.0),
        (-5.53, -5.22, 0.0),
        (-0.39, -1.98, 0.0),
        (2.98, -5.51, 0.0),
        (4.83, -2.02, 0.0),
        (9.70, -3.63, 0.0),
        (12.23, 1.25, 0.0),
        (3.42, 0.66, 0.0),
        (2.92, 4.03, 0.0),
    ]
)
L_SHAPE = Polygon([(0, 0, 0), (6, 0, 0), (6, 2, 0), (2, 2, 0), (2, 8, 0), (0, 8, 0)])
STAR = Polygon(
    [
        (5.0, 0.0, 0.0),
        (1.5, 1.2, 0.0),
        (1.55, 4.76, 0.0),
        (-0.57, 1.9, 0.0),
        (-4.05, 2.94, 0.0),
        (-1.85, 0.0, 0.0),
        (-4.05, -2.94, 0.0),
        (-0.57, -1.9, 0.0),
        (1.55, -4.76, 0.0),
        (1.5, -1.2, 0.0),
    ]
)
KITE = Polygon([(5, 0, 0), (0, 5, 0), (-5, 0, 0), (0, -2.5, 0)])


def _dumbbell(pinch: float) -> Polygon:
    """20x10 pocket pinched to *pinch* at x=10 by two facing reflex notches."""
    h = pinch / 2.0
    return Polygon(
        [
            (0.0, 0.0, 0.0),
            (9.0, 0.0, 0.0),
            (10.0, 5.0 - h, 0.0),
            (11.0, 0.0, 0.0),
            (20.0, 0.0, 0.0),
            (20.0, 10.0, 0.0),
            (11.0, 10.0, 0.0),
            (10.0, 5.0 + h, 0.0),
            (9.0, 10.0, 0.0),
            (0.0, 10.0, 0.0),
        ]
    )


ISLAND_OUTER = Polygon([(0, 0, 0), (14, 0, 0), (14, 10, 0), (0, 10, 0)])
ISLAND_HOLE = Polygon([(5.5, 3.5, 0), (8.5, 3.5, 0), (8.5, 6.5, 0), (5.5, 6.5, 0)])


@dataclass(frozen=True)
class Pocket:
    """One benchmark pocket and the tool it is machined with.

    Tool diameter is 1.0 by default and 0.5 for the small-featured shapes
    (l_shape, star, kite), matching the test suite's per-pocket parameters.
    """

    name: str
    polygon: Polygon
    tool_diameter: float
    holes: list[Polygon] | None = None


POCKETS = [
    Pocket("square", SQUARE, 1.0),
    Pocket("irregular", IRREGULAR, 1.0),
    Pocket("l_shape", L_SHAPE, 0.5),
    Pocket("star", STAR, 0.5),
    Pocket("kite", KITE, 0.5),
    Pocket("dumbbell", _dumbbell(2.4), 1.0),
    Pocket("island", ISLAND_OUTER, 1.0, [ISLAND_HOLE]),
]

# Tiny synthetic smoke pocket for --quick: a 2.5 mm square (tool 1.0) yields a ~47-op
# trochoid that audits end-to-end in ~8 s (O(n^2) in cut count -> small n is cheap), so
# the subprocess test can drive the whole generate -> audit -> report path without any
# real benchmark's cost. NOT a benchmark; never part of the reference set.
SMOKE_POCKET = Pocket("smoke", Polygon([[0, 0, 0], [2.5, 0, 0], [2.5, 2.5, 0], [0, 2.5, 0]]), 1.0)


def audit_pocket(pocket: Pocket) -> dict:
    """Generate + audit one pocket; return its per-pocket record (audit wall-clock timed).

    Generation is negligible; the timed span is the audit alone (the expensive,
    stock-depleting replay), per the per-pocket wall-clock requirement.
    """
    result = trochoidal_mat_toolpath_circular(
        pocket.polygon,
        tool_diameter=pocket.tool_diameter,
        clearance_z=CLEARANCE_Z,
        holes=pocket.holes,
    )
    start = time.perf_counter()
    report = audit_toolpath_engagement(
        pocket.polygon,
        result,
        tool_diameter=pocket.tool_diameter,
        tea_cap=TEA_CAP_RAD,
        holes=pocket.holes,
    )
    wall_clock_s = time.perf_counter() - start

    return {
        "name": pocket.name,
        "tool_diameter": pocket.tool_diameter,
        "stepover": 0.4 * pocket.tool_diameter,
        "pitch": 0.75 * pocket.tool_diameter,
        "clearance_z": CLEARANCE_Z,
        "n_operations": len(report.operations),
        "engaged_ops": report.engaged_ops,
        "max_tea_deg": math.degrees(report.max_tea),
        "max_tea_rad": report.max_tea,
        "cap_violations": report.cap_violations,
        "stations": sum(op.stations for op in report.operations),
        "wall_clock_s": wall_clock_s,
    }


def build_report(records: list[dict]) -> dict:
    """Assemble the report-level aggregate (BLUF numbers) around the per-pocket records."""
    worst = max(records, key=lambda r: r["max_tea_deg"])
    return {
        "generator": "trochoidal_mat_toolpath_circular",
        "certifier": "audit_toolpath_engagement",
        "tea_cap_deg": TEA_CAP_DEG,
        "clearance_z": CLEARANCE_Z,
        "n_pockets": len(records),
        "worst_pocket": worst["name"],
        "worst_tea_deg": worst["max_tea_deg"],
        "total_cap_violations": sum(r["cap_violations"] for r in records),
        "total_stations": sum(r["stations"] for r in records),
        "total_wall_clock_s": sum(r["wall_clock_s"] for r in records),
        "pockets": records,
    }


def render_markdown(data: dict) -> str:
    """Render the BLUF-first markdown report from the aggregate data."""
    cap = data["tea_cap_deg"]
    lines = [
        "# Engagement baseline -- existing trochoid generator (SP2 reference)",
        "",
        (
            f"**BLUF.** The existing `trochoidal_mat_toolpath_circular` applies NO engagement "
            f"regulation. Audited against a {cap:.0f} deg tool-engagement-angle (TEA) cap across "
            f"{data['n_pockets']} pocket(s), the worst pocket is **{data['worst_pocket']}** with a "
            f"worst TEA of **{data['worst_tea_deg']:.1f} deg** (cap {cap:.0f} deg), "
            f"**{data['total_cap_violations']}** cap violations in total, and a total audit "
            f"wall-clock of **{data['total_wall_clock_s']:.1f} s**. "
            f"This baseline is the reference SP2 must beat."
        ),
        "",
        "## Per-pocket engagement",
        "",
        "| Pocket | Tool diameter | Worst TEA (deg) | Cap violations | Stations | Wall-clock (s) |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for r in data["pockets"]:
        lines.append(f"| {r['name']} | {r['tool_diameter']:.1f} | {r['max_tea_deg']:.1f} | " f"{r['cap_violations']} | {r['stations']} | {r['wall_clock_s']:.1f} |")
    lines.append(
        f"| **total** | | **{data['worst_tea_deg']:.1f}** (worst) | "
        f"**{data['total_cap_violations']}** | **{data['total_stations']}** | "
        f"**{data['total_wall_clock_s']:.1f}** |"
    )
    lines += [
        "",
        "## Method",
        "",
        (
            f"Each pocket is machined with `trochoidal_mat_toolpath_circular` (tool diameter per "
            f"pocket, `clearance_z={data['clearance_z']}`, generator defaults otherwise) and "
            f"replayed through `audit_toolpath_engagement` against a {cap:.0f} deg cap. The audit "
            f"depletes an exact `Stock` per cut and measures the TEA each motion sees *before* it "
            f"removes its own material; a cap violation is a motion the sound certifier cannot "
            f"prove stays under the cap. Worst TEA is a reporting quantity (it may exceed a half "
            f"turn for a fully immersed cut); the violation count is the certified decision. "
            f"Wall-clock is the per-pocket audit time (`time.perf_counter`); generation is "
            f"negligible. Regenerate with `python scripts/engagement_baseline.py`."
        ),
        "",
    ]
    return "\n".join(lines)


def select_pockets(args: argparse.Namespace) -> list[Pocket]:
    """Resolve the CLI selection to an ordered list of pockets (fail loud on unknown).

    ``--quick`` returns the tiny synthetic smoke pocket (harness self-test); otherwise the
    requested benchmark subset, or all benchmarks.
    """
    if args.quick:
        return [SMOKE_POCKET]
    by_name = {p.name: p for p in POCKETS}
    if args.pockets:
        names = [n.strip() for n in args.pockets.split(",") if n.strip()]
        unknown = [n for n in names if n not in by_name]
        if unknown:
            raise SystemExit(f"unknown pocket(s) {unknown}; valid names: {[p.name for p in POCKETS]}")
        return [by_name[n] for n in names]
    return list(POCKETS)


def main() -> None:
    parser = argparse.ArgumentParser(description="Baseline engagement audit of the existing generator (SP2 reference).")
    parser.add_argument("--out", default="docs/benchmarks", help="output directory for the .md and .json (default: docs/benchmarks)")
    parser.add_argument("--pockets", default=None, help="comma-separated subset of pocket names to audit")
    parser.add_argument("--quick", action="store_true", help="audit only the single fastest small pocket (fast subprocess test)")
    args = parser.parse_args()

    pockets = select_pockets(args)

    records = []
    for pocket in pockets:
        record = audit_pocket(pocket)
        records.append(record)
        print(
            f"{record['name']:<10} worst_tea={record['max_tea_deg']:6.1f}deg "
            f"viol={record['cap_violations']:<6} stations={record['stations']:<7} "
            f"wall={record['wall_clock_s']:7.1f}s",
            flush=True,
        )

    data = build_report(records)

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "engagement_baseline.json").write_text(json.dumps(data, indent=2) + "\n")
    (out_dir / "engagement_baseline.md").write_text(render_markdown(data))
    print(
        f"wrote {out_dir / 'engagement_baseline.md'} and {out_dir / 'engagement_baseline.json'} "
        f"({data['n_pockets']} pocket(s), worst={data['worst_pocket']} @ {data['worst_tea_deg']:.1f} deg)",
        flush=True,
    )


if __name__ == "__main__":
    main()
