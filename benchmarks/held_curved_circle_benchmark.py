"""Reproducible timing of the native curved-circle query on the prepared Held pockets.

One job: measure, record with identity, and compare against a recorded
baseline. Plotting lives in `held_curved_circle_performance_plot`; the
construction diagrams live in `held_native_curved_circle_plot`.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import math
import pathlib
import platform
import statistics
import subprocess
import time
from dataclasses import dataclass
from datetime import datetime
from datetime import timezone
from typing import Any

from benchmarks.errors import BenchmarkError
from benchmarks.held_native_curve_import import import_held_boundary
from benchmarks.held_reference_cases import CANONICAL_CASE_NAMES
from benchmarks.held_reference_cases import load_held_reference_case
from compas_cgal import __version__ as PACKAGE_VERSION
from compas_cgal import _coverage_2 as native

SCHEMA_VERSION = 1
DEFAULT_REPEATS = 3
SAMPLE_PARAMETER = 0.5
REPOSITORY = pathlib.Path(__file__).resolve().parents[1]
BASELINE_PATH = REPOSITORY / "benchmarks/results/held_curved_circle_baseline.json"

# Regression policy. Run-to-run spread of a case median on one quiet machine
# was measured under 1.5x on 2026-09-08; 3x and 5x leave room for a loaded
# box while catching the order-of-magnitude regressions this lane produced
# (100x to 1000x). Across machines only order of magnitude is comparable, so
# the factor is 10x and the absolute ceilings carry the rest: 0.25 s per query
# is the corpus witness budget, and a 20 ms case median is fifty times the
# recorded medians on the M1 Max, above any CI runner slowdown seen so far.
SAME_MACHINE_MEDIAN_FACTOR = 3.0
SAME_MACHINE_MAX_FACTOR = 5.0
CROSS_MACHINE_FACTOR = 10.0
ABSOLUTE_PIECE_CEILING_S = 0.25
ABSOLUTE_MEDIAN_CEILING_S = 0.02

# Named geometry stops are recorded outcomes of a case, never process failures.
GEOMETRY_STOPS = (
    native.InvalidNativeBoundaryCurveError,
    native.InvalidNativeBoundaryChainError,
    native.ReachableDomainConstructionError,
    native.InvalidBoundaryCircleContactError,
    native.BoundaryContactConstructionError,
    native.CoverageTransitionError,
    native.InvalidCoverageGeometryError,
    native.InvalidNativeBoundaryMedialInputError,
    native.NoPositiveNativeBoundaryCircleError,
    native.NativeBoundaryMedialConstructionError,
)


class CurvedCircleBenchmarkError(BenchmarkError):
    """Base error of the curved-circle benchmark."""


class BaselineSchemaError(CurvedCircleBenchmarkError):
    """A recorded benchmark does not carry the schema this module reads."""


class EmptyCaseBenchmarkError(CurvedCircleBenchmarkError):
    """A case with no completed query has no timing statistics."""


class BuildIdentityError(CurvedCircleBenchmarkError):
    """The repository identity could not be read from git."""


@dataclass(frozen=True)
class MachineIdentity:
    machine: str
    system: str
    cpu_brand: str
    python: str

    def matches(self, other: MachineIdentity) -> bool:
        """Same hardware class: the same-machine factors apply."""
        return (self.machine, self.system, self.cpu_brand) == (other.machine, other.system, other.cpu_brand)


@dataclass(frozen=True)
class BuildIdentity:
    commit: str
    dirty: bool
    package_version: str
    recorded_at: str


@dataclass(frozen=True)
class PieceTiming:
    piece_index: int
    kind: str
    seconds: float  # median of samples
    samples: tuple[float, ...]


@dataclass(frozen=True)
class CaseBenchmark:
    case: str
    native_pieces: int
    completed_queries: int
    stop_piece_index: int | None
    stop_type: str | None
    import_seconds: float
    pieces: tuple[PieceTiming, ...]

    def _seconds(self) -> list[float]:
        if not self.pieces:
            raise EmptyCaseBenchmarkError(f"{self.case}: no completed query to summarise.")
        return sorted(piece.seconds for piece in self.pieces)

    @property
    def median_seconds(self) -> float:
        return statistics.median(self._seconds())

    @property
    def p95_seconds(self) -> float:
        ordered = self._seconds()
        return ordered[math.ceil(0.95 * len(ordered)) - 1]

    @property
    def max_seconds(self) -> float:
        return self._seconds()[-1]


@dataclass(frozen=True)
class CurvedCircleBenchmark:
    schema_version: int
    identity: MachineIdentity
    build: BuildIdentity
    repeats: int
    parameter: float
    cases: dict[str, CaseBenchmark]

    def to_json(self) -> str:
        return json.dumps(dataclasses.asdict(self), indent=2, sort_keys=True) + "\n"

    @classmethod
    def from_json(cls, text: str) -> CurvedCircleBenchmark:
        raw: dict[str, Any] = json.loads(text)
        if raw.get("schema_version") != SCHEMA_VERSION:
            raise BaselineSchemaError(f"schema_version {raw.get('schema_version')!r} is not {SCHEMA_VERSION}.")
        cases = {
            name: CaseBenchmark(
                case=case["case"],
                native_pieces=case["native_pieces"],
                completed_queries=case["completed_queries"],
                stop_piece_index=case["stop_piece_index"],
                stop_type=case["stop_type"],
                import_seconds=case["import_seconds"],
                pieces=tuple(
                    PieceTiming(piece_index=p["piece_index"], kind=p["kind"], seconds=p["seconds"], samples=tuple(p["samples"]))
                    for p in case["pieces"]
                ),
            )
            for name, case in raw["cases"].items()
        }
        return cls(
            schema_version=raw["schema_version"],
            identity=MachineIdentity(**raw["identity"]),
            build=BuildIdentity(**raw["build"]),
            repeats=raw["repeats"],
            parameter=raw["parameter"],
            cases=cases,
        )


@dataclass(frozen=True)
class RegressionFinding:
    case: str
    kind: str
    current: float
    bound: float
    message: str


def cpu_brand() -> str:
    """Marketing name of the CPU, or 'unknown' where the platform does not publish one."""
    system = platform.system()
    if system == "Darwin":
        return subprocess.run(["sysctl", "-n", "machdep.cpu.brand_string"], capture_output=True, text=True, check=True).stdout.strip()
    if system == "Linux":
        for line in pathlib.Path("/proc/cpuinfo").read_text(encoding="utf-8").splitlines():
            if line.startswith("model name"):
                return line.split(":", 1)[1].strip()
        return "unknown"
    return platform.processor() or "unknown"


def machine_identity() -> MachineIdentity:
    return MachineIdentity(machine=platform.machine(), system=platform.system(), cpu_brand=cpu_brand(), python=platform.python_version())


def build_identity(repository: pathlib.Path = REPOSITORY) -> BuildIdentity:
    try:
        commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=repository, capture_output=True, text=True, check=True).stdout.strip()
        status = subprocess.run(["git", "status", "--porcelain"], cwd=repository, capture_output=True, text=True, check=True).stdout
    except (subprocess.CalledProcessError, FileNotFoundError) as error:
        raise BuildIdentityError(f"git identity unavailable in {repository}: {error}") from error
    return BuildIdentity(
        commit=commit,
        dirty=bool(status.strip()),
        package_version=PACKAGE_VERSION,
        recorded_at=datetime.now(timezone.utc).isoformat(timespec="seconds"),
    )


def measure_case(case_name: str, repeats: int = DEFAULT_REPEATS) -> CaseBenchmark:
    """Query every native piece `repeats` times; a geometry stop ends the case and is recorded."""
    if repeats < 1:
        raise CurvedCircleBenchmarkError("repeats must be positive.")
    case = load_held_reference_case(case_name)
    started = time.perf_counter()
    owner = import_held_boundary(case)
    import_seconds = time.perf_counter() - started
    tool = float(case.tool_radius.value)
    primitives = list(owner.cycle.primitives)
    pieces: list[PieceTiming] = []
    stop_index: int | None = None
    stop_type: str | None = None
    for index, primitive in enumerate(primitives):
        samples: list[float] = []
        try:
            for _ in range(repeats):
                query_started = time.perf_counter()
                owner.circle_on_piece(index, SAMPLE_PARAMETER, tool)
                samples.append(time.perf_counter() - query_started)
        except GEOMETRY_STOPS as stop:
            stop_index, stop_type = index, type(stop).__name__
            break
        pieces.append(PieceTiming(piece_index=index, kind=primitive.kind, seconds=statistics.median(samples), samples=tuple(samples)))
    return CaseBenchmark(
        case=case_name,
        native_pieces=len(primitives),
        completed_queries=len(pieces),
        stop_piece_index=stop_index,
        stop_type=stop_type,
        import_seconds=import_seconds,
        pieces=tuple(pieces),
    )


def measure_all(repeats: int = DEFAULT_REPEATS, case_names: tuple[str, ...] = tuple(CANONICAL_CASE_NAMES)) -> CurvedCircleBenchmark:
    return CurvedCircleBenchmark(
        schema_version=SCHEMA_VERSION,
        identity=machine_identity(),
        build=build_identity(),
        repeats=repeats,
        parameter=SAMPLE_PARAMETER,
        cases={name: measure_case(name, repeats) for name in case_names},
    )


def compare(current: CurvedCircleBenchmark, baseline: CurvedCircleBenchmark) -> list[RegressionFinding]:
    """Every way the current run is worse than, or behaves differently from, the baseline."""
    findings: list[RegressionFinding] = []
    same_machine = current.identity.matches(baseline.identity)
    median_factor = SAME_MACHINE_MEDIAN_FACTOR if same_machine else CROSS_MACHINE_FACTOR
    max_factor = SAME_MACHINE_MAX_FACTOR if same_machine else CROSS_MACHINE_FACTOR
    regime = "same machine" if same_machine else f"cross-machine ({current.identity.cpu_brand} vs {baseline.identity.cpu_brand})"
    for name, recorded in baseline.cases.items():
        present = current.cases.get(name)
        if present is None:
            findings.append(RegressionFinding(name, "missing_case", 0.0, 0.0, f"{name}: not measured in the current run."))
            continue
        if present.completed_queries != recorded.completed_queries:
            findings.append(RegressionFinding(name, "completed", present.completed_queries, recorded.completed_queries,
                                              f"{name}: completed {present.completed_queries} queries, baseline {recorded.completed_queries}; behaviour changed, re-record deliberately if intended."))
        if (present.stop_piece_index, present.stop_type) != (recorded.stop_piece_index, recorded.stop_type):
            findings.append(RegressionFinding(name, "stop", present.stop_piece_index or -1, recorded.stop_piece_index or -1,
                                              f"{name}: stop {present.stop_type} at {present.stop_piece_index}, baseline {recorded.stop_type} at {recorded.stop_piece_index}."))
        if not present.pieces or not recorded.pieces:
            continue
        if present.max_seconds > ABSOLUTE_PIECE_CEILING_S:
            findings.append(RegressionFinding(name, "absolute_piece_ceiling", present.max_seconds, ABSOLUTE_PIECE_CEILING_S,
                                              f"{name}: slowest query {present.max_seconds:.4f} s exceeds the {ABSOLUTE_PIECE_CEILING_S} s ceiling."))
        if present.median_seconds > ABSOLUTE_MEDIAN_CEILING_S:
            findings.append(RegressionFinding(name, "absolute_median_ceiling", present.median_seconds, ABSOLUTE_MEDIAN_CEILING_S,
                                              f"{name}: median query {present.median_seconds:.4f} s exceeds the {ABSOLUTE_MEDIAN_CEILING_S} s ceiling."))
        median_bound = recorded.median_seconds * median_factor
        if present.median_seconds > median_bound:
            findings.append(RegressionFinding(name, "median", present.median_seconds, median_bound,
                                              f"{name}: median {present.median_seconds:.5f} s > {median_factor}x baseline {recorded.median_seconds:.5f} s ({regime})."))
        max_bound = recorded.max_seconds * max_factor
        if present.max_seconds > max_bound:
            findings.append(RegressionFinding(name, "max", present.max_seconds, max_bound,
                                              f"{name}: slowest {present.max_seconds:.5f} s > {max_factor}x baseline {recorded.max_seconds:.5f} s ({regime})."))
    return findings


def write_json(benchmark: CurvedCircleBenchmark, path: pathlib.Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(benchmark.to_json(), encoding="utf-8")


def read_json(path: pathlib.Path) -> CurvedCircleBenchmark:
    return CurvedCircleBenchmark.from_json(path.read_text(encoding="utf-8"))


def summary_lines(benchmark: CurvedCircleBenchmark) -> list[str]:
    lines = [f"{benchmark.identity.cpu_brand} · {benchmark.build.commit[:10]}{' (dirty)' if benchmark.build.dirty else ''} · repeats {benchmark.repeats}"]
    for name, case in benchmark.cases.items():
        stop = "none" if case.stop_type is None else f"{case.stop_type} at piece {case.stop_piece_index}"
        if case.pieces:
            lines.append(f"{name:<22} {case.completed_queries:>3}/{case.native_pieces:<3} median {case.median_seconds * 1000:7.3f} ms  p95 {case.p95_seconds * 1000:7.3f} ms  max {case.max_seconds * 1000:7.3f} ms  stop {stop}")
        else:
            lines.append(f"{name:<22} {case.completed_queries:>3}/{case.native_pieces:<3} no completed query  stop {stop}")
    return lines


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=pathlib.Path, required=True, help="JSON record to write")
    parser.add_argument("--repeats", type=int, default=DEFAULT_REPEATS)
    args = parser.parse_args()
    benchmark = measure_all(repeats=args.repeats)
    write_json(benchmark, args.output)
    print("\n".join(summary_lines(benchmark)))
    print(f"RECORD {args.output}")


if __name__ == "__main__":
    main()
