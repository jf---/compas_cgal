"""Fail-loud lookups into cProfile statistics.

A bucket or table row names the profiled functions it reads. A name that is
not in the profile (code renamed or removed) or that matches more than one
function raises, so a refactor breaks the figure instead of silently zeroing
a bucket.
"""

from __future__ import annotations

import pstats
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping
from typing import Protocol

# pstats key: (filename, first line, function name); builtins use ("~", 0, "<built-in method ...>").
FunctionKey = tuple[str, int, str]
# pstats value: (primitive calls, total calls, self seconds, cumulative seconds, callers).
FunctionRow = tuple[int, int, float, float, Mapping[FunctionKey, object]]


class MissingProfileFunctionError(Exception):
    """A function a bucket or row reads is not in the profile."""


class AmbiguousProfileFunctionError(Exception):
    """A function name matches more than one profiled function."""


class ProfileStats(Protocol):
    """What the bucket and table code reads from a `pstats.Stats`."""

    @property
    def stats(self) -> Mapping[FunctionKey, FunctionRow]:
        """Per-function rows."""
        ...

    @property
    def total_tt(self) -> float:
        """Total profiled seconds."""
        ...


@dataclass(frozen=True)
class LoadedStats:
    """A `ProfileStats` read from a ``.prof`` file.

    Attributes:
        stats: Per-function rows.
        total_tt: Total profiled seconds.
    """

    stats: Mapping[FunctionKey, FunctionRow]
    total_tt: float


def load_stats(path: Path) -> LoadedStats:
    """Read a cProfile file.

    Args:
        path: The ``.prof`` file.

    Returns:
        Its statistics.
    """
    raw = pstats.Stats(str(path))
    return LoadedStats(stats=dict(raw.stats), total_tt=float(raw.total_tt))  # type: ignore[attr-defined]


def module_name(filename: str) -> str:
    """The basename of a pstats filename (``"~"`` for builtins)."""
    return filename.rsplit("/", 1)[-1]


def _module_matches(filename: str, module: str) -> bool:
    if "/" in module:
        return filename == module or filename.endswith("/" + module)
    return module_name(filename) == module


def function_key(stats: ProfileStats, name: str, module: str | None = None) -> FunctionKey:
    """The pstats key of the one profiled function called *name* (in *module*, if given).

    Args:
        stats: The profile.
        name: The function name as pstats records it.
        module: The defining file, as a basename (``"stock.py"``) or a
            repository-relative path (``"src/compas_cgal/stock.py"``).

    Returns:
        Its pstats key.

    Raises:
        MissingProfileFunctionError: No function matches.
        AmbiguousProfileFunctionError: More than one function matches.
    """
    hits = [key for key in stats.stats if key[2] == name and (module is None or _module_matches(key[0], module))]
    where = f"{module}:{name}" if module is not None else name
    if not hits:
        raise MissingProfileFunctionError(f"{where} is not in the profile; was it renamed?")
    if len(hits) > 1:
        raise AmbiguousProfileFunctionError(f"{where} matches {len(hits)} profiled functions; name the module.")
    return hits[0]


def function_row(stats: ProfileStats, name: str, module: str | None = None) -> FunctionRow:
    """The pstats row of the one profiled function called *name* (see `function_key`)."""
    return stats.stats[function_key(stats, name, module)]


def is_within(stats: ProfileStats, inner: FunctionKey, outer: FunctionKey) -> bool:
    """Whether *inner* is *outer* or is called, directly or transitively, from *outer*.

    Args:
        stats: The profile.
        inner: The candidate descendant.
        outer: The candidate ancestor.

    Returns:
        True when *outer* is reachable from *inner* along caller edges.
    """
    seen: set[FunctionKey] = set()
    pending = [inner]
    while pending:
        key = pending.pop()
        if key == outer:
            return True
        if key in seen:
            continue
        seen.add(key)
        callers = stats.stats[key][4] if key in stats.stats else {}
        pending.extend(callers)
    return False


def cumulative_s(stats: ProfileStats, name: str, module: str | None = None) -> float:
    """Cumulative seconds of the one function called *name* (see `function_row`)."""
    return function_row(stats, name, module)[3]


def call_count(stats: ProfileStats, name: str, module: str | None = None) -> int:
    """Total calls of the one function called *name* (see `function_row`)."""
    return function_row(stats, name, module)[1]
