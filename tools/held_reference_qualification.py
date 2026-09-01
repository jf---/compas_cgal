"""Qualify the existing generator on the four Held reference pockets."""

from __future__ import annotations

import math
import time
from collections.abc import Callable
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Literal
from typing import Protocol
from typing import TypeAlias

from benchmarks.errors import EmptyToolpathError
from benchmarks.errors import MalformedRecordError
from benchmarks.held_reference_cases import CANONICAL_CASE_NAMES
from benchmarks.held_reference_cases import load_all_held_reference_cases
from benchmarks.runner import generate_toolpath
from benchmarks.spec import PocketSpec
from compas_cgal.adaptive.errors import UncertifiedZeroGuideEdgeError
from compas_cgal.toolpath import DegeneratePrimitiveError

EXPECTED_CASES = CANONICAL_CASE_NAMES
DEFAULT_REPORT_PATH = Path("docs/benchmarks/held_reference_qualification.md")

QualificationStatus = Literal["passed", "failed"]


class GeneratedPath(Protocol):
    """Minimal existing-generator result consumed by qualification."""

    @property
    def operations(self) -> Sequence[object]: ...


Generator: TypeAlias = Callable[[PocketSpec], GeneratedPath]
Clock: TypeAlias = Callable[[], float]

# These are the existing, closed product failures qualification may preserve.
# Everything else is a broken run and deliberately propagates.
RECORDED_PRODUCT_FAILURES = (
    UncertifiedZeroGuideEdgeError,
    EmptyToolpathError,
    DegeneratePrimitiveError,
)


@dataclass(frozen=True)
class QualificationOutcome:
    """One serial generator outcome with its corpus context."""

    case_name: str
    primitive_count: int
    projection_count: int
    status: QualificationStatus
    operation_count: int
    generation_seconds: float
    failure_type: str | None
    failure_message: str | None

    def __post_init__(self) -> None:
        if self.case_name not in EXPECTED_CASES:
            raise MalformedRecordError(f"Unknown qualification case: {self.case_name!r}.")
        counts = (self.primitive_count, self.projection_count, self.operation_count)
        if any(count < 0 for count in counts):
            raise MalformedRecordError("Qualification counts must be non-negative.")
        if not math.isfinite(self.generation_seconds) or self.generation_seconds < 0.0:
            raise MalformedRecordError("Qualification generation time must be finite and non-negative.")
        if self.status == "passed":
            if self.operation_count == 0 or self.failure_type is not None or self.failure_message is not None:
                raise MalformedRecordError("A passed qualification needs operations and no failure.")
        elif self.operation_count != 0 or not self.failure_type or not self.failure_message:
            raise MalformedRecordError("A failed qualification needs one named failure and no operations.")

    @classmethod
    def passed(
        cls,
        *,
        case_name: str,
        primitive_count: int,
        projection_count: int,
        operation_count: int,
        generation_seconds: float,
    ) -> QualificationOutcome:
        """Build a successful qualification outcome."""
        return cls(
            case_name=case_name,
            primitive_count=primitive_count,
            projection_count=projection_count,
            status="passed",
            operation_count=operation_count,
            generation_seconds=generation_seconds,
            failure_type=None,
            failure_message=None,
        )

    @classmethod
    def failed(
        cls,
        *,
        case_name: str,
        primitive_count: int,
        projection_count: int,
        generation_seconds: float,
        error: Exception,
    ) -> QualificationOutcome:
        """Build an outcome for one approved named product failure."""
        return cls(
            case_name=case_name,
            primitive_count=primitive_count,
            projection_count=projection_count,
            status="failed",
            operation_count=0,
            generation_seconds=generation_seconds,
            failure_type=type(error).__name__,
            failure_message=str(error),
        )


def qualify_cases(
    generator: Generator,
    *,
    clock: Clock = time.perf_counter,
) -> tuple[QualificationOutcome, ...]:
    """Run the injected generator serially on exactly four canonical cases."""
    cases = load_all_held_reference_cases()
    names = tuple(case.name for case in cases)
    if names != EXPECTED_CASES:
        raise MalformedRecordError(f"Qualification cases {names!r} do not match {EXPECTED_CASES!r}.")

    outcomes: list[QualificationOutcome] = []
    for case in cases:
        started = clock()
        try:
            generated = generator(case.pocket_spec())
            operation_count = len(generated.operations)
            if operation_count == 0:
                raise EmptyToolpathError(f"{case.name}: existing generator returned no operations.")
        except RECORDED_PRODUCT_FAILURES as error:
            outcomes.append(
                QualificationOutcome.failed(
                    case_name=case.name,
                    primitive_count=len(case.boundary.primitives),
                    projection_count=case.projection_vertex_count,
                    generation_seconds=clock() - started,
                    error=error,
                )
            )
            continue

        outcomes.append(
            QualificationOutcome.passed(
                case_name=case.name,
                primitive_count=len(case.boundary.primitives),
                projection_count=case.projection_vertex_count,
                operation_count=operation_count,
                generation_seconds=clock() - started,
            )
        )
    return tuple(outcomes)


def render_qualification_markdown(outcomes: Sequence[QualificationOutcome]) -> str:
    """Render qualification outcomes as concise product evidence."""
    names = tuple(outcome.case_name for outcome in outcomes)
    if names != EXPECTED_CASES:
        raise MalformedRecordError(f"Report cases {names!r} do not match {EXPECTED_CASES!r}.")

    lines = [
        "# Held reference qualification",
        "",
        "This product-evidence run exercises the existing generator serially. It does not compare generator performance or establish Held parity.",
        "",
        "| case | primitives | projection vertices | status | operations | generation seconds | named failure |",
        "| --- | ---: | ---: | --- | ---: | ---: | --- |",
    ]
    for outcome in outcomes:
        failure = "—"
        if outcome.failure_type is not None:
            failure = f"`{outcome.failure_type}`: {outcome.failure_message}"
        lines.append(
            f"| {outcome.case_name} | {outcome.primitive_count} | {outcome.projection_count} "
            f"| {outcome.status} | {outcome.operation_count} | {outcome.generation_seconds:.6f} | {failure} |"
        )
    lines.extend(
        [
            "",
            "A pass means only that this reconstructed input produced a non-empty path in this run.",
            "",
        ]
    )
    return "\n".join(lines)


def write_qualification_report(
    outcomes: Sequence[QualificationOutcome],
    path: Path = DEFAULT_REPORT_PATH,
) -> Path:
    """Write the rendered qualification report to *path*."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(render_qualification_markdown(outcomes), encoding="utf-8")
    return path


def main() -> None:
    """Run the real serial qualification and write its documentation evidence."""
    path = write_qualification_report(qualify_cases(generate_toolpath))
    print(path)


if __name__ == "__main__":
    main()
