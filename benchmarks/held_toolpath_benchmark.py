"""Time the Held toolpath generators and gate their results against a recorded baseline.

Five concrete workloads call the existing pipelines stage by stage: guide,
inward offset, initial placement, engagement refinement, connectors. Each run
records the toolpath's counts and lengths with the time, so a regression in
geometry is caught as surely as a regression in speed. The exact coverage
replay is validation, not generation, and is timed on its own clock only when
asked for.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import math
import pathlib
import time
from collections.abc import Callable
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

from benchmarks.errors import BenchmarkError
from benchmarks.held_boundary_order_spacing import refine_boundary_ordered_path
from benchmarks.held_exact_motion_coverage import require_full_motion_coverage
from benchmarks.held_figure5_boundary_path import Figure5BoundaryTransition
from benchmarks.held_figure5_boundary_path import Figure5CounterclockwiseCircle
from benchmarks.held_figure5_boundary_path import _side_lengths
from benchmarks.held_figure5_boundary_path import _transition
from benchmarks.held_figure5_engagement_refinement import refine_figure5_engagement
from benchmarks.held_figure5_path import _paper_candidate
from benchmarks.held_figure5_path import build_figure5_approximate_path
from benchmarks.held_figure5_path import build_held_reference_path
from benchmarks.held_figure5_raw_guide import build_figure5_raw_guide
from benchmarks.held_figure5_raw_guide import build_held_reference_raw_guide
from benchmarks.held_reference_cases import HeldReferenceCase
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.held_reference_figures import figure7_inward_offset
from benchmarks.held_standard_placement import StandardPlacementFragmentationError
from benchmarks.held_standard_placement import maximum_predecessor_engagement
from benchmarks.performance_record import REPOSITORY
from benchmarks.performance_record import BuildIdentity
from benchmarks.performance_record import MachineIdentity
from benchmarks.performance_record import build_identity
from benchmarks.performance_record import machine_identity
from benchmarks.performance_record import time_factor
from benchmarks.tools.held_toolpath_progress import _containment_report
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import WorldXY

SCHEMA_VERSION = 1
BASELINE_PATH = REPOSITORY / "benchmarks/results/held_toolpath_baseline.json"
# Reporting doubles of exact lengths and angles: round-off, never geometry.
REPORTED_LENGTH_TOLERANCE_MM = 1e-6
REPORTED_ANGLE_TOLERANCE_DEG = 1e-9


class ToolpathBenchmarkError(BenchmarkError):
    """Base error of the toolpath benchmark."""


class BaselineSchemaError(ToolpathBenchmarkError):
    """A recorded benchmark does not carry the schema this module reads."""


class MissingCeilingError(ToolpathBenchmarkError):
    """A gated workload has no absolute ceiling."""


class UnknownWorkloadError(ToolpathBenchmarkError):
    """The requested workload is not one of the concrete generators."""


@dataclass(frozen=True)
class ToolpathMetrics:
    circles: int
    connectors: int
    path_length_mm: float
    max_engagement_deg: float
    containment_rejections: int
    unresolved_engagements: int = 0


@dataclass(frozen=True)
class WorkloadRun:
    workload: str
    metrics: ToolpathMetrics
    stage_seconds: dict[str, float]
    generation_seconds: float
    validation_seconds: float | None = None


@dataclass(frozen=True)
class ToolpathBenchmark:
    schema_version: int
    identity: MachineIdentity
    build: BuildIdentity
    runs: dict[str, WorkloadRun]

    def to_json(self) -> str:
        return json.dumps(dataclasses.asdict(self), indent=2, sort_keys=True) + "\n"

    @classmethod
    def from_json(cls, text: str) -> ToolpathBenchmark:
        raw: dict[str, Any] = json.loads(text)
        if raw.get("schema_version") != SCHEMA_VERSION:
            raise BaselineSchemaError(f"schema_version {raw.get('schema_version')!r} is not {SCHEMA_VERSION}.")
        runs = {
            name: WorkloadRun(
                workload=run["workload"],
                metrics=ToolpathMetrics(**run["metrics"]),
                stage_seconds=dict(run["stage_seconds"]),
                generation_seconds=run["generation_seconds"],
                validation_seconds=run["validation_seconds"],
            )
            for name, run in raw["runs"].items()
        }
        return cls(schema_version=raw["schema_version"], identity=MachineIdentity(**raw["identity"]), build=BuildIdentity(**raw["build"]), runs=runs)


@dataclass(frozen=True)
class Finding:
    workload: str
    kind: str
    current: float
    bound: float
    message: str


class _Stages:
    """Stage clocks of one generation; every stage is disjoint and the total is their sum."""

    def __init__(self) -> None:
        self.seconds: dict[str, float] = {}
        self._started = 0.0

    def start(self) -> None:
        self._started = time.perf_counter()

    def stop(self, stage: str) -> None:
        self.seconds[stage] = time.perf_counter() - self._started


def _metrics(case: HeldReferenceCase, circles: tuple[Figure5CounterclockwiseCircle, ...], transitions: tuple[Figure5BoundaryTransition, ...], boundary: tuple[Point2[WorldXY], ...]) -> ToolpathMetrics:
    length = sum(2 * math.pi * float(circle.radius.value) for circle in circles)
    for transition in transitions:
        samples = transition.samples
        length += sum(math.dist((float(a.x), float(a.y)), (float(b.x), float(b.y))) for a, b in zip(samples, samples[1:]))
    max_engagement = 0.0
    unresolved = 0
    for first, second in zip(circles, circles[1:]):
        try:
            angle = math.degrees(float(maximum_predecessor_engagement(_paper_candidate(first), _paper_candidate(second), case.tool_radius)))
        except StandardPlacementFragmentationError:
            unresolved += 1
            continue
        max_engagement = max(max_engagement, angle)
    rejections, _ = _containment_report(circles, boundary)
    return ToolpathMetrics(
        circles=len(circles),
        connectors=len(transitions),
        path_length_mm=length,
        max_engagement_deg=max_engagement,
        containment_rejections=rejections,
        unresolved_engagements=unresolved,
    )


def _connectors(boundary: tuple[Point2[WorldXY], ...], circles: tuple[Figure5CounterclockwiseCircle, ...]) -> tuple[Figure5BoundaryTransition, ...]:
    lengths = _side_lengths(boundary)
    return tuple(_transition(boundary, lengths, a, b) for a, b in zip(circles, circles[1:]))


def _finish(name: str, case: HeldReferenceCase, circles, transitions, boundary, stages: _Stages, validate: bool) -> WorkloadRun:
    metrics = _metrics(case, circles, transitions, boundary)
    validation: float | None = None
    if validate:
        started = time.perf_counter()
        require_full_motion_coverage(case, circles, transitions, pathlib.Path("build/held-toolpath-benchmark") / f"{name}.exact-coverage.json")
        validation = time.perf_counter() - started
    return WorkloadRun(workload=name, metrics=metrics, stage_seconds=dict(stages.seconds), generation_seconds=sum(stages.seconds.values()), validation_seconds=validation)


def figure5_refined(validate: bool = False) -> WorkloadRun:
    """Figure 5 standard-model draft: publisher-start polygon proposal, corner-aware engagement refinement."""
    case = load_held_reference_case("figure5")
    stages = _Stages()
    stages.start(); guide = build_figure5_raw_guide(case); stages.stop("guide")
    stages.start(); components = figure7_inward_offset(case).components; stages.stop("offset")
    stages.start(); result = build_figure5_approximate_path(guide, components); stages.stop("initial")
    cap = EngagementCap.build(math.radians(float(case.tea_cap)))
    stages.start(); refined = refine_figure5_engagement(result.path.circles, components[0], case.tool_radius, cap); stages.stop("refine")
    stages.start(); transitions = _connectors(components[0], refined.circles); stages.stop("connectors")
    return _finish("figure5_refined", case, refined.circles, transitions, components[0], stages, validate)


def _figure8(name: str, case_name: str, *, spacing: str, qualify_contacts: bool, corner_approaches: bool, validate: bool) -> WorkloadRun:
    case = load_held_reference_case(case_name)
    stages = _Stages()
    stages.start(); guide = build_held_reference_raw_guide(case); stages.stop("guide")
    stages.start(); components = figure7_inward_offset(case).components; stages.stop("offset")
    stages.start(); result = build_held_reference_path(guide, components, qualify_contacts=qualify_contacts); stages.stop("initial")
    cap = EngagementCap.build(math.radians(float(case.tea_cap)))
    stages.start()
    if spacing == "boundary":
        _, refined = refine_boundary_ordered_path(result, components[0], case.tool_radius, cap, corner_approaches=corner_approaches)
    else:
        refined = refine_figure5_engagement(result.path.circles, components[0], case.tool_radius, cap, corner_approaches=corner_approaches)
    stages.stop("refine")
    stages.start(); transitions = _connectors(components[0], refined.circles); stages.stop("connectors")
    return _finish(name, case, refined.circles, transitions, components[0], stages, validate)


def figure8_upper_lane(validate: bool = False) -> WorkloadRun:
    return _figure8("figure8_upper_lane", "figure8_upper", spacing="lane", qualify_contacts=False, corner_approaches=False, validate=validate)


def figure8_upper_boundary(validate: bool = False) -> WorkloadRun:
    return _figure8("figure8_upper_boundary", "figure8_upper", spacing="boundary", qualify_contacts=False, corner_approaches=False, validate=validate)


def figure8_crossed_skis_boundary(validate: bool = False) -> WorkloadRun:
    return _figure8("figure8_crossed_skis_boundary", "figure8_crossed_skis", spacing="boundary", qualify_contacts=False, corner_approaches=False, validate=validate)


def figure8_monstera_boundary_qualified_corners(validate: bool = False) -> WorkloadRun:
    return _figure8("figure8_monstera_boundary_qualified_corners", "figure8_monstera", spacing="boundary", qualify_contacts=True, corner_approaches=True, validate=validate)


WORKLOADS: dict[str, Callable[..., WorkloadRun]] = {
    "figure5_refined": figure5_refined,
    "figure8_upper_lane": figure8_upper_lane,
    "figure8_upper_boundary": figure8_upper_boundary,
    "figure8_crossed_skis_boundary": figure8_crossed_skis_boundary,
    "figure8_monstera_boundary_qualified_corners": figure8_monstera_boundary_qualified_corners,
}

# Absolute generation-time ceilings per workload, seconds: three times the
# generation measured on the Apple M1 Max on 2026-09-08 (20.6, 88.2, 87.9,
# 34.9 and 436 s), where all but a fraction of a second is the initial
# placement. Workloads under a minute gate in tests/benchmarks; the others in
# benchmarks/acceptance, which pytest does not collect.
ABSOLUTE_CEILINGS_S: dict[str, float] = {
    "figure5_refined": 62.0,
    "figure8_upper_lane": 265.0,
    "figure8_upper_boundary": 264.0,
    "figure8_crossed_skis_boundary": 105.0,
    "figure8_monstera_boundary_qualified_corners": 1310.0,
}
FAST_WORKLOADS: tuple[str, ...] = ("figure5_refined", "figure8_crossed_skis_boundary")
SLOW_WORKLOADS: tuple[str, ...] = tuple(name for name in WORKLOADS if name not in FAST_WORKLOADS)


def measure(names: tuple[str, ...] = tuple(WORKLOADS), validate: bool = False) -> ToolpathBenchmark:
    unknown = [name for name in names if name not in WORKLOADS]
    if unknown:
        raise UnknownWorkloadError(f"unknown workloads: {unknown}; known: {list(WORKLOADS)}")
    return ToolpathBenchmark(
        schema_version=SCHEMA_VERSION,
        identity=machine_identity(),
        build=build_identity(),
        runs={name: WORKLOADS[name](validate=validate) for name in names},
    )


def compare(current: ToolpathBenchmark, baseline: ToolpathBenchmark, ceilings: Mapping[str, float] = ABSOLUTE_CEILINGS_S) -> list[Finding]:
    """Every way the current run differs in geometry from, or is slower than, the baseline."""
    findings: list[Finding] = []
    factor, regime = time_factor(current.identity, baseline.identity)
    for name, recorded in baseline.runs.items():
        present = current.runs.get(name)
        if present is None:
            findings.append(Finding(name, "missing_workload", 0.0, 0.0, f"{name}: not measured in the current run."))
            continue
        if name not in ceilings:
            raise MissingCeilingError(f"{name} has no absolute ceiling; add it to ABSOLUTE_CEILINGS_S.")
        was, now = recorded.metrics, present.metrics
        exact = (("circles", now.circles, was.circles), ("connectors", now.connectors, was.connectors))
        for kind, current_value, recorded_value in exact:
            if current_value != recorded_value:
                findings.append(Finding(name, kind, current_value, recorded_value, f"{name}: {kind} {current_value} differs from baseline {recorded_value}."))
        if abs(now.path_length_mm - was.path_length_mm) > REPORTED_LENGTH_TOLERANCE_MM:
            findings.append(Finding(name, "path_length", now.path_length_mm, was.path_length_mm, f"{name}: path length {now.path_length_mm:.6f} mm differs from baseline {was.path_length_mm:.6f} mm."))
        if abs(now.max_engagement_deg - was.max_engagement_deg) > REPORTED_ANGLE_TOLERANCE_DEG:
            findings.append(Finding(name, "max_engagement", now.max_engagement_deg, was.max_engagement_deg, f"{name}: max engagement {now.max_engagement_deg:.6f} deg differs from baseline {was.max_engagement_deg:.6f} deg."))
        if now.containment_rejections != was.containment_rejections:
            findings.append(Finding(name, "containment_rejections", now.containment_rejections, was.containment_rejections, f"{name}: containment rejections {now.containment_rejections} differ from baseline {was.containment_rejections}."))
        if now.unresolved_engagements != was.unresolved_engagements:
            findings.append(Finding(name, "unresolved_engagements", now.unresolved_engagements, was.unresolved_engagements, f"{name}: unresolved engagements {now.unresolved_engagements} differ from baseline {was.unresolved_engagements}."))
        bound = recorded.generation_seconds * factor
        if present.generation_seconds > bound:
            findings.append(Finding(name, "generation_time", present.generation_seconds, bound, f"{name}: generation {present.generation_seconds:.2f} s > {factor}x baseline {recorded.generation_seconds:.2f} s ({regime})."))
        if present.generation_seconds > ceilings[name]:
            findings.append(Finding(name, "absolute_ceiling", present.generation_seconds, ceilings[name], f"{name}: generation {present.generation_seconds:.2f} s exceeds the {ceilings[name]} s ceiling."))
    return findings


def write_json(benchmark: ToolpathBenchmark, path: pathlib.Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(benchmark.to_json(), encoding="utf-8")


def read_json(path: pathlib.Path) -> ToolpathBenchmark:
    return ToolpathBenchmark.from_json(path.read_text(encoding="utf-8"))


def summary_lines(benchmark: ToolpathBenchmark) -> list[str]:
    lines = [f"{benchmark.identity.cpu_brand} · {benchmark.build.commit[:10]}{' (dirty)' if benchmark.build.dirty else ''}"]
    for name, run in benchmark.runs.items():
        m = run.metrics
        stages = " ".join(f"{stage} {seconds:.2f}" for stage, seconds in run.stage_seconds.items())
        validation = "" if run.validation_seconds is None else f"  validation {run.validation_seconds:.1f} s"
        lines.append(f"{name:<44} {run.generation_seconds:8.2f} s [{stages}]  circles {m.circles:>5}  connectors {m.connectors:>5}  length {m.path_length_mm:12.3f} mm  max eng {m.max_engagement_deg:8.4f} deg  rejections {m.containment_rejections}{validation}")
    return lines


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--workloads", nargs="+", default=list(WORKLOADS), choices=list(WORKLOADS))
    parser.add_argument("--validate", action="store_true", help="also run the exact coverage replay on its own clock")
    args = parser.parse_args()
    benchmark = measure(tuple(args.workloads), validate=args.validate)
    write_json(benchmark, args.output)
    print("\n".join(summary_lines(benchmark)))
    print(f"RECORD {args.output}")


if __name__ == "__main__":
    main()
