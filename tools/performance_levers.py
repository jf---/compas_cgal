"""The lever mapping: which profiled code each performance lever replaces, and its share per workload.

`benchmarks/results/performance_levers.json` (committed) names, per lever, the
functions whose time the lever replaces. A cell is the profiled share of a
workload's time spent in that code: an upper bound on the lever's saving
there. Every name is checked twice: it must be defined in the source file the
mapping names, and it must be in every profile the mapping applies it to. A
rename therefore breaks the figure instead of silently zeroing a cell.
"""

from __future__ import annotations

import ast
import json
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any
from typing import Mapping

from tools.performance_stats import FunctionKey
from tools.performance_stats import ProfileStats
from tools.performance_stats import function_key
from tools.performance_stats import is_within
from tools.performance_workloads import ROOT
from tools.performance_workloads import WORKLOADS_BY_KEY

LEVERS_PATH = ROOT / "benchmarks" / "results" / "performance_levers.json"
SCHEMA = "performance-levers/1"
PERCENT = 100.0


class LeverStatus(Enum):
    """How far a lever is: wired into generation, built beside it, proposed, or in design."""

    WIRED = "wired"
    BUILT = "built"
    PROPOSED = "proposed"
    DESIGN = "design"


class TimeKind(Enum):
    """Which pstats time of a function counts."""

    CUMULATIVE = "cumulative"
    SELF = "self"


class MalformedLeverMappingError(Exception):
    """The lever mapping does not follow its schema."""


class UnknownLeverFunctionError(Exception):
    """A mapped function is not defined in the source file the mapping names."""


class LeverDoubleCountError(Exception):
    """Within one lever, one function's counted time lies inside another's."""


@dataclass(frozen=True)
class LeverFunction:
    """One function whose time a lever replaces.

    Attributes:
        module: Repository-relative path of the defining file.
        function: The function name.
        time: Cumulative (Python subtree) or self time (native wrapper).
        workloads: Keys of the workloads in which this function's time counts.
    """

    module: str
    function: str
    time: TimeKind
    workloads: tuple[str, ...]


@dataclass(frozen=True)
class Lever:
    """One performance lever.

    Attributes:
        key: JSON key.
        label: Row label.
        status: How far it is.
        functions: The code it replaces.
        estimated: Workload key -> why the profile cannot attribute the cell exactly.
    """

    key: str
    label: str
    status: LeverStatus
    functions: tuple[LeverFunction, ...]
    estimated: Mapping[str, str]

    def applies_to(self, workload: str) -> bool:
        """Whether any of the lever's code runs in *workload*."""
        return any(workload in function.workloads for function in self.functions)

    @classmethod
    def build(cls, data: Mapping[str, Any]) -> Lever:
        """Validate one lever entry of the mapping.

        Args:
            data: The JSON object.

        Returns:
            The lever.

        Raises:
            MalformedLeverMappingError: An unknown status, time kind or workload, or an
                estimate for a workload the lever does not apply to.
            UnknownLeverFunctionError: A function is not defined in its module.
        """
        try:
            key, label = str(data["key"]), str(data["label"])
            status = LeverStatus(data["status"])
            functions = tuple(
                LeverFunction(str(entry["module"]), str(entry["function"]), TimeKind(entry["time"]), tuple(str(w) for w in entry["workloads"])) for entry in data["functions"]
            )
            estimated = {str(workload): str(reason) for workload, reason in data["estimated"].items()}
        except (KeyError, ValueError, TypeError, AttributeError) as error:
            raise MalformedLeverMappingError(f"lever {data.get('key')!r}: {error}") from error
        for function in functions:
            unknown = [workload for workload in function.workloads if workload not in WORKLOADS_BY_KEY]
            if unknown:
                raise MalformedLeverMappingError(f"{key}: {function.function} names unknown workloads {unknown}")
            _require_defined(function)
        lever = cls(key=key, label=label, status=status, functions=functions, estimated=estimated)
        stray = [workload for workload in estimated if not lever.applies_to(workload)]
        if stray:
            raise MalformedLeverMappingError(f"{key}: estimate for workloads it does not apply to: {stray}")
        return lever


def _require_defined(function: LeverFunction) -> None:
    path = ROOT / function.module
    if not path.is_file():
        raise UnknownLeverFunctionError(f"{function.module} does not exist")
    tree = ast.parse(path.read_text(), filename=str(path))
    names = {node.name for node in ast.walk(tree) if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))}
    if function.function not in names:
        raise UnknownLeverFunctionError(f"{function.module} defines no function {function.function!r}; was it renamed?")


def load_levers(path: Path = LEVERS_PATH) -> tuple[Lever, ...]:
    """Read and validate the lever mapping.

    Args:
        path: The mapping file.

    Returns:
        The levers, in row order.

    Raises:
        MalformedLeverMappingError: Wrong schema or a malformed lever.
    """
    document = json.loads(path.read_text())
    if document.get("schema") != SCHEMA:
        raise MalformedLeverMappingError(f"{path}: schema {document.get('schema')!r}, expected {SCHEMA!r}")
    levers = tuple(Lever.build(entry) for entry in document["levers"])
    keys = [lever.key for lever in levers]
    if len(set(keys)) != len(keys):
        raise MalformedLeverMappingError(f"{path}: duplicate lever keys {keys}")
    return levers


def lever_share(lever: Lever, workload: str, stats: ProfileStats) -> float | None:
    """The share of *workload*'s profiled time in the code *lever* replaces.

    Args:
        lever: The lever.
        workload: The workload key.
        stats: That workload's profile.

    Returns:
        Percent of the profiled total, or ``None`` when none of the lever's code runs in the workload.

    Raises:
        MissingProfileFunctionError: A function the mapping applies here is not in the profile.
        LeverDoubleCountError: One counted function lies inside another's counted time.
    """
    counted: list[tuple[FunctionKey, TimeKind]] = [(function_key(stats, f.function, f.module), f.time) for f in lever.functions if workload in f.workloads]
    if not counted:
        return None
    for i, (outer, outer_time) in enumerate(counted):
        for j, (inner, _) in enumerate(counted):
            if i == j:
                continue
            if inner == outer or (outer_time is TimeKind.CUMULATIVE and is_within(stats, inner, outer)):
                raise LeverDoubleCountError(f"{lever.key} in {workload}: {inner[2]} lies inside the counted subtree of {outer[2]}")
    seconds = sum(stats.stats[key][3] if time is TimeKind.CUMULATIVE else stats.stats[key][2] for key, time in counted)
    return PERCENT * seconds / stats.total_tt
