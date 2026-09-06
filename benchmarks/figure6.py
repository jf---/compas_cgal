"""Reproduce Held & Pfeiffer's Figure 6: path length against engagement cap.

Figure 6 is the only directly reproducible published result in the paper. It plots
path length against maximum engagement angle for one pocket, and it is where the
paper's headline claim lives: an engagement-controlled path is substantially
SHORTER than a constant-spacing one at the same cap. Reproducing the comparison
turns "we are inside Held's range" from an assertion against a single
unreproducible timing figure into a like-for-like curve a reader can regenerate.

WHAT IS AND IS NOT REPRODUCED. The paper's three curves are its own standard,
contour-aware, and MATHSM generators on its own pocket; none of those are
available here, and no curve below is claimed to be one of them. What is
reproduced is the PROTOCOL and the axes on the reconstructed Figure 5 polygon
projection, with two generators of ours:

* **engagement-controlled** -- `engagement_controlled_toolpath`, whose advance is
  regulated by the exact engagement predicate against the depleting stock. Its cap
  is an input, so its length moves with the cap. This is the analogue of the
  paper's own contribution.
* **constant spacing (MATHSM protocol)** -- `trochoidal_mat_toolpath_circular`
  driven by `stepover`, with the cap met by brute-force search over spacing
  (`benchmarks.mathsm`). The generator has no cap input at all; the cap enters
  only through which trial is selected. This is the analogue of the paper's
  MATHSM baseline.

THE COMPARISON IS NOT SYMMETRIC, and the report says so on every row. The baseline
is SELECTED to meet the cap, so it meets it by construction. The controlled
generator is ASKED to meet the cap and may fail -- its bridge cuts between
machining circles are not regulated. Its measured maximum is therefore printed
beside its length, and a length ratio at a cap it misses is comparing a compliant
path against one that is not.

Engagement is read after the entry cut for both curves, for the reason given in
`benchmarks.pathmetrics`: the first circle after a plunge is a full slot for any
generator at any setting, and filtering on it would empty the comparison.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from typing import Dict
from typing import List
from typing import Optional
from typing import Sequence
from typing import Tuple

from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.mathsm import SPACING_SWEEP_TOOL_DIAMETERS
from benchmarks.mathsm import MathsmPoint
from benchmarks.mathsm import shortest_within_cap
from benchmarks.mathsm import sweep_spacing
from benchmarks.pathmetrics import CLEARANCE_Z_TOOL_DIAMETERS
from benchmarks.pathmetrics import REFERENCE_CAP_DEG
from benchmarks.pathmetrics import PathMetrics
from benchmarks.pathmetrics import demonstrated_exceedances_after_entry
from benchmarks.pathmetrics import measure_path
from benchmarks.spec import PocketSpec
from benchmarks.toolpath_coverage import require_toolpath_coverage
from compas_cgal.engagement_toolpath import engagement_controlled_toolpath
from compas_cgal.toolpath import ToolpathResult

MARKDOWN_NAME = "figure6.md"
JSON_NAME = "figure6.json"

# Caps swept across the figure, in degrees. The lower end is where neither
# generator can comply on any pocket entered from solid stock, and it is kept
# precisely so the report can say that; the upper end stays inside the exact
# kernel's contract of a half turn (`benchmarks.spec.MAX_CAP_DEG`).
FIGURE6_CAPS: Tuple[float, ...] = (20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0)

# Emitted verbatim under the conclusion. A reader meeting these two curves cold
# must not be able to mistake either for one of the paper's own.
CURVE_LEGEND = (
    "**controlled** is `engagement_controlled_toolpath`, asked for the cap on its "
    "row; its advance is regulated by the exact engagement predicate. "
    "**constant spacing** is `trochoidal_mat_toolpath_circular` at a fixed "
    "stepover, reproducing the MATHSM protocol: the generator has no cap input, so "
    "the cap enters only by selecting the shortest trial spacing that measured at "
    "or below it. Neither is one of the paper's curves. Both engagement figures "
    "are measured AFTER each chain's entry cut, which is a full slot for any "
    "generator entering solid stock and pins the raw maximum near a full turn "
    "regardless of the setting being swept. The exceedance column is an extra "
    "exact-predicate check carried for the controlled curve only; the baseline's "
    "compliance is already what selected it."
)

TABLE_HEADER = (
    "| cap (deg) | controlled length | controlled max TEA after entry (deg) | controlled meets cap | "
    "controlled exceedances after entry | spacing (tool diam.) | constant-spacing length | "
    "constant-spacing max TEA after entry (deg) | length ratio |"
)
TABLE_RULE = "| ---: | ---: | ---: | :--- | ---: | ---: | ---: | ---: | ---: |"

NO_BASELINE_CELL = "no compliant spacing"
EMPTY_CELL = "—"
YES_CELL = "yes"
NO_CELL = "no"


@dataclass(frozen=True)
class Figure6Point:
    """One cap value on the reproduction.

    Attributes:
        cap_deg: The engagement cap both curves were asked for, in degrees.
        controlled: Metrics of the engagement-controlled path generated at that
            cap.
        controlled_exceedances_after_entry: Non-entry cut motions of that path
            where the exact predicate fired against this cap. A sampled lower
            bound, never a certificate.
        mathsm: The shortest constant-spacing trial that met the cap, or None when
            no trial spacing did.
    """

    cap_deg: float
    controlled: PathMetrics
    controlled_exceedances_after_entry: int
    mathsm: Optional[MathsmPoint]

    @property
    def length_ratio(self) -> Optional[float]:
        """Controlled length as a multiple of the baseline's, or None without one."""
        if self.mathsm is None or self.mathsm.metrics.length <= 0.0:
            return None
        return self.controlled.length / self.mathsm.metrics.length

    @property
    def controlled_meets_cap(self) -> bool:
        """Whether the controlled path's measured maximum stayed at or under the cap.

        Measured after the entry cut, and measured -- not certified. It says no
        audited station exceeded the cap away from the entry.
        """
        return self.controlled.max_tea_after_entry_deg <= self.cap_deg


def reference_pocket(tea_cap_deg: float = REFERENCE_CAP_DEG) -> PocketSpec:
    """The pocket the reproduction is reported on.

    One definition, so the artifact, the docs page, and every test describe the
    same geometry.

    Args:
        tea_cap_deg: Cap recorded on the returned spec. Irrelevant to the figure,
            which replaces it per point, and defaulted to the measurement
            reference so an unswept use of this pocket measures rather than tests.

    Returns:
        The spec.
    """
    reconstructed = load_held_reference_case("figure5").pocket_spec()
    return PocketSpec.build(
        name=reconstructed.name,
        family=reconstructed.family,
        polygon=reconstructed.polygon,
        tool_diameter=reconstructed.tool_diameter,
        tea_cap_deg=tea_cap_deg,
        holes=reconstructed.holes,
        params=reconstructed.params,
    )


def controlled_path(spec: PocketSpec, cap_deg: float) -> ToolpathResult:
    """Generate the engagement-controlled path for one cap.

    The single place the figure fixes this generator's parameters, so a point, a
    test, and the artifact all describe the same path. Clearance matches
    `benchmarks.mathsm.constant_spacing_path`, so the two curves' lengths differ by
    their cutting, not by their retract height.

    Args:
        spec: The pocket and tool.
        cap_deg: The engagement cap in degrees.

    Returns:
        The generated toolpath.

    Warns:
        UnavoidableEngagementWarning: Raised through from the generator when it
            had to emit circles at positions the predicate refuses. Deliberately
            not suppressed: it is the generator's own honest report of where the
            cap could not be honoured.
    """
    return engagement_controlled_toolpath(
        spec.polygon,
        tool_diameter=spec.tool_diameter,
        tea_cap_deg=cap_deg,
        holes=list(spec.holes),
        clearance_z=CLEARANCE_Z_TOOL_DIAMETERS * spec.tool_diameter,
    )


@dataclass(frozen=True)
class Figure6Run:
    """Everything one reproduction measured.

    Attributes:
        spec: The pocket and tool both curves ran on.
        points: One point per cap, in plotting order.
        trials: Every constant-spacing trial the baseline was selected from, in
            sweep order. Carried because the shape of this sweep is the evidence
            that the baseline search has to be brute force.
    """

    spec: PocketSpec
    points: Tuple[Figure6Point, ...]
    trials: Tuple[MathsmPoint, ...]


def run_figure6(spec: PocketSpec, caps: Sequence[float] = FIGURE6_CAPS, spacings_tool_diameters: Sequence[float] = SPACING_SWEEP_TOOL_DIAMETERS) -> Figure6Run:
    """Sweep the baseline once, then measure both curves at every cap.

    The constant-spacing trials are generated and measured ONCE, before the cap
    loop, because their engagement is a property of their spacing alone
    (`benchmarks.pathmetrics`). Re-measuring them per cap would cost the whole
    sweep again and would let the baseline's x-values move with the cap axis they
    are plotted against.

    Args:
        spec: The pocket and tool. Its own cap is replaced per point.
        caps: Engagement caps in degrees, in plotting order.
        spacings_tool_diameters: Trial spacings for the baseline, as multiples of
            the tool diameter.

    Returns:
        The completed run.

    Raises:
        InvalidCapError: A cap lies outside the exact kernel's contract.
    """
    trials = sweep_spacing(spec, spacings_tool_diameters)
    return Figure6Run(spec=spec, points=tuple(figure6_points(spec, caps, trials)), trials=tuple(trials))


def figure6_points(spec: PocketSpec, caps: Sequence[float], trials: Sequence[MathsmPoint]) -> List[Figure6Point]:
    """Measure the controlled curve at every cap against an already-swept baseline.

    Exact design-pocket coverage is mandatory before accepting each path's
    metrics. Unsupported motions and nonempty residuals raise; engagement
    failures remain reported independently on coverage-qualified paths.

    Args:
        spec: The pocket and tool. Its own cap is replaced per point.
        caps: Engagement caps in degrees, in plotting order.
        trials: Constant-spacing trials from `benchmarks.mathsm.sweep_spacing`.

    Returns:
        One point per cap, in input order.

    Raises:
        InvalidCapError: A cap lies outside the exact kernel's contract.
    """
    points: List[Figure6Point] = []
    for cap_deg in caps:
        at_cap = PocketSpec.build(
            name=spec.name,
            family=spec.family,
            polygon=spec.polygon,
            tool_diameter=spec.tool_diameter,
            tea_cap_deg=cap_deg,
            holes=spec.holes,
            params=spec.params,
        )
        result = controlled_path(at_cap, cap_deg)
        require_toolpath_coverage(at_cap, result)
        points.append(
            Figure6Point(
                cap_deg=cap_deg,
                controlled=measure_path(at_cap, result),
                controlled_exceedances_after_entry=demonstrated_exceedances_after_entry(at_cap, result),
                mathsm=shortest_within_cap(trials, cap_deg),
            )
        )
    return points


def render_figure6_markdown(run: Figure6Run) -> str:
    """Render the reproduction, conclusion first.

    The ordering is built into this function rather than left to whoever edits the
    output, so every re-run leads with the answer.

    Args:
        run: The completed reproduction.

    Returns:
        The markdown document.
    """
    lines: List[str] = [f"# Figure 6 reproduction — path length against engagement cap ({run.spec.name})", ""]
    lines += _conclusion(run.points, run.spec)
    lines += [CURVE_LEGEND, ""]
    lines += ["## Per-cap comparison", "", TABLE_HEADER, TABLE_RULE]
    for point in run.points:
        lines.append(_row(point))
    lines.append("")
    if run.trials:
        lines += _trial_table(run.trials)
    return "\n".join(lines)


def _conclusion(points: Sequence[Figure6Point], spec: PocketSpec) -> List[str]:
    """The answer, stated before any evidence.

    Args:
        points: Measured points.
        spec: The pocket both curves ran on.

    Returns:
        The conclusion paragraph as markdown lines, blank-line terminated.
    """
    if not points:
        return ["No cap was measured, so the reproduction produced no comparison.", ""]

    ratios = [p.length_ratio for p in points if p.length_ratio is not None]
    missing = [p for p in points if p.mathsm is None]
    met = [p for p in points if p.controlled_meets_cap]
    tool = f"tool {spec.tool_diameter:g}"

    if not ratios:
        return [
            f"No trial spacing met any of the {len(points)} caps swept on {spec.name} ({tool}), so this run carries no length "
            f"comparison; the constant-spacing baseline's best measured engagement is above every cap requested, and the "
            f"sweep needs widening before the figure says anything.",
            "",
        ]

    mean_ratio = sum(ratios) / len(ratios)
    verdict = "shorter" if mean_ratio < 1.0 else "longer"
    lead = (
        f"Across the {len(ratios)} of {len(points)} caps where a constant-spacing path met the cap at all, the "
        f"engagement-controlled path is on average **{mean_ratio:.2f}x** the baseline's length — **{verdict}** — on "
        f"{spec.name} ({tool})."
    )

    # The only rows nobody can argue with: both paths were measured under the cap
    # they are filed under. Everything else compares a baseline selected for
    # compliance against a path merely asked for it.
    like_for_like = [p.length_ratio for p in points if p.length_ratio is not None and p.controlled_meets_cap]
    if like_for_like:
        lead += (
            f" At the **{len(like_for_like)}** cap(s) where BOTH paths met the cap it is **{sum(like_for_like) / len(like_for_like):.2f}x**, "
            f"which is the fully like-for-like comparison in this table."
        )
    else:
        lead += " No cap had both paths meet it, so no row here is a fully like-for-like comparison."

    if missing:
        floor = min(p.controlled.max_tea_after_entry_deg for p in points)
        lead += (
            f" The remaining **{len(missing)}** cap(s) have no compliant spacing at all: no trial spacing reached them, and "
            f"neither did the controlled generator, whose lowest measured maximum over the sweep is {floor:.1f} deg."
        )

    never_fired = [p for p in points if p.controlled_exceedances_after_entry == 0]
    lead += (
        f" The controlled path stayed under its own cap away from the entry at **{len(met)} of {len(points)}** caps, and the "
        f"exact predicate never fired away from an entry at **{len(never_fired)} of {len(points)}**."
    )
    if len(met) < len(points):
        lead += (
            " Where it did not, the ratio on that row compares a baseline selected for compliance against a controlled path "
            "that is not compliant, and the measured maximum in the row says by how much."
        )
    return [lead, ""]


def _row(point: Figure6Point) -> str:
    """One table row, with the baseline cells collapsed when no trial complied.

    Args:
        point: The measured point.

    Returns:
        The markdown row.
    """
    cells = [
        f"{point.cap_deg:.0f}",
        f"{point.controlled.length:.1f}",
        f"{point.controlled.max_tea_after_entry_deg:.1f}",
        YES_CELL if point.controlled_meets_cap else NO_CELL,
        str(point.controlled_exceedances_after_entry),
    ]
    if point.mathsm is None:
        cells += [EMPTY_CELL, NO_BASELINE_CELL, EMPTY_CELL, EMPTY_CELL]
    else:
        ratio = point.length_ratio
        cells += [
            f"{point.mathsm.spacing_tool_diameters:.3f}",
            f"{point.mathsm.metrics.length:.1f}",
            f"{point.mathsm.metrics.max_tea_after_entry_deg:.1f}",
            EMPTY_CELL if ratio is None else f"{ratio:.2f}",
        ]
    return "| " + " | ".join(cells) + " |"


def _trial_table(trials: Sequence[MathsmPoint]) -> List[str]:
    """The complete constant-spacing sweep used by exhaustive selection.

    Args:
        trials: Every trial, in sweep order.

    Returns:
        Markdown lines, blank-line terminated.
    """
    return (
        [
            "## Constant-spacing trials",
            "",
            "The fixed sweep is reported in full below. Baseline selection evaluates every compliant trial and chooses the shortest without assuming spacing orders engagement.",
            "",
            "| spacing (tool diam.) | length | cut motions | max TEA after entry (deg) | raw max TEA (deg) |",
            "| ---: | ---: | ---: | ---: | ---: |",
        ]
        + [
            f"| {t.spacing_tool_diameters:.3f} | {t.metrics.length:.1f} | {t.metrics.cut_motions} | {t.metrics.max_tea_after_entry_deg:.1f} | {t.metrics.max_tea_deg:.1f} |"
            for t in trials
        ]
        + [""]
    )


def figure6_payload(run: Figure6Run) -> Dict[str, Any]:
    """The reproduction as a JSON-serialisable mapping.

    Args:
        run: The completed reproduction.

    Returns:
        The payload.
    """
    spec = run.spec
    return {
        "pocket": {
            "name": spec.name,
            "family": spec.family,
            "tool_diameter": spec.tool_diameter,
            "params": dict(spec.params),
            "holes": [[[point.x, point.y, point.z] for point in ring.points] for ring in spec.holes],
        },
        "engagement_measured_at_cap_deg": REFERENCE_CAP_DEG,
        "points": [
            {
                "cap_deg": p.cap_deg,
                "controlled_length": p.controlled.length,
                "controlled_cut_motions": p.controlled.cut_motions,
                "controlled_entry_cuts": p.controlled.entry_cuts,
                "controlled_max_tea_deg": p.controlled.max_tea_deg,
                "controlled_max_tea_after_entry_deg": p.controlled.max_tea_after_entry_deg,
                "controlled_exceedances_after_entry": p.controlled_exceedances_after_entry,
                "controlled_meets_cap": p.controlled_meets_cap,
                "mathsm_spacing_tool_diameters": None if p.mathsm is None else p.mathsm.spacing_tool_diameters,
                "mathsm_length": None if p.mathsm is None else p.mathsm.metrics.length,
                "mathsm_max_tea_after_entry_deg": None if p.mathsm is None else p.mathsm.metrics.max_tea_after_entry_deg,
                "length_ratio": p.length_ratio,
            }
            for p in run.points
        ],
        "spacing_trials": [
            {
                "spacing_tool_diameters": t.spacing_tool_diameters,
                "length": t.metrics.length,
                "cut_motions": t.metrics.cut_motions,
                "entry_cuts": t.metrics.entry_cuts,
                "max_tea_deg": t.metrics.max_tea_deg,
                "max_tea_after_entry_deg": t.metrics.max_tea_after_entry_deg,
            }
            for t in run.trials
        ],
    }


def write_figure6(run: Figure6Run, out_dir: Path) -> Tuple[Path, Path]:
    """Write the markdown and JSON artifacts.

    Args:
        run: The completed reproduction.
        out_dir: Destination directory; created when missing.

    Returns:
        ``(markdown_path, json_path)``.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    md_path = out_dir / MARKDOWN_NAME
    json_path = out_dir / JSON_NAME
    md_path.write_text(render_figure6_markdown(run), encoding="utf-8")
    json_path.write_text(json.dumps(figure6_payload(run), indent=2), encoding="utf-8")
    return md_path, json_path
