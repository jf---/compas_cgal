"""The performance baseline and the regression check against it.

Current entries: per record version (record path, record SHA-256) only the
entries read by the newest adapter version count (``current_entries``); older
ones stay in the file, superseded, never rewritten.

Baseline (orchestrator ruling on RV40 I3, 2026-10-23): per (workload, stage,
machine), the BEST quiet entry, the one with the smallest gated minimum (ties:
the later instant, then history order). A slower result never becomes the
baseline by being recorded, so a regression cannot creep in 9 % steps. To
accept a slower baseline, append a ``Rebaseline`` naming the accepted quiet
entry: then only the key's quiet entries measured at or after that entry
(compared as instants, ``PerfEntry.when``) are candidates, and the best of
them is the baseline. The latest rebaseline per key (by its ``recorded_at``)
governs. Indicative entries never set a baseline.

Self-exclusion (RV40 I2): ``check`` leaves every history entry read from the
same record bytes as a checked entry (same ``record_sha256``) out of the
baseline, so checking a record after ingesting it compares it with the
history before it, never with itself.

Check: each entry of a fresh record is compared with its key's baseline on
the gated clock (process CPU when both have it, wall when neither does) by its
minimum, since interference only adds time. A QUIET entry whose minimum
exceeds ``REGRESSION_RATIO_MAX`` times the baseline's fails the check
(``PerformanceRegressionError``); an indicative entry over the threshold is
reported and never fails, since its load alone can explain it. The decision is
``current > REGRESSION_RATIO_MAX * baseline`` in binary64; the threshold is a
policy ratio, not a geometric tolerance.
"""

from __future__ import annotations

from collections.abc import Iterable
from collections.abc import Mapping
from collections.abc import Sequence
from collections.abc import Set as AbstractSet
from dataclasses import dataclass
from datetime import datetime
from datetime import timezone
from typing import Optional

from tools.perf_history_entry import ClockStatistic
from tools.perf_history_entry import PerfEntry
from tools.perf_history_entry import parse_instant
from tools.perf_history_rebaseline import Rebaseline
from tools.perf_history_store import History

# Largest admitted ratio of a quiet run's gated minimum to its baseline's (dimensionless); above it
# the check fails. Derived from the backfilled history (2026-10-23, 52 repeated-run clocks): the
# widest process-CPU max/min spread of one record's runs is 1.027 (coverage_certificate_omega_star,
# 3 runs, 1-minute load ~29); the minimum of N runs varies less than max/min, so 1.10 is ~3.7x the
# widest measured CPU spread. Loaded wall spreads (to 1.70, one outlier 22.6 at load 175) never gate:
# only quiet entries do. No quiet repeated run exists yet; re-derive once one does.
REGRESSION_RATIO_MAX = 1.10

# Clock statistics whose minimum is comparable: an observed minimum over one or more runs.
_COMPARABLE: Mapping[ClockStatistic, str] = {
    ClockStatistic.RUNS: "observed minimum",
    ClockStatistic.SINGLE: "observed minimum",
    ClockStatistic.REPORTED_MEDIAN: "reported median",
    ClockStatistic.SUM_OF_SINGLES: "sum of single runs",
    ClockStatistic.SUM_OF_MEDIANS: "sum of medians",
}

Key = tuple[str, str, str]

# Ordering stand-in for an entry with no recorded time: the earliest possible instant.
_UNTIMED = datetime.min.replace(tzinfo=timezone.utc)


class GatedMetricMismatchError(Exception):
    """An entry and its baseline gate on different clocks (CPU vs wall)."""


class ClockStatisticMismatchError(Exception):
    """An entry's gated clock and its baseline's are different statistics (e.g. a minimum vs a median)."""


class PerformanceRegressionError(Exception):
    """A quiet entry's gated minimum exceeds ``REGRESSION_RATIO_MAX`` times its baseline's."""

    def __init__(self, regressions: Sequence[Comparison]) -> None:
        lines = [comparison.describe() for comparison in regressions]
        super().__init__(f"{len(regressions)} quiet regression(s) over {REGRESSION_RATIO_MAX}x baseline:\n" + "\n".join(lines))
        self.regressions = tuple(regressions)


def current_entries(entries: Sequence[PerfEntry]) -> tuple[PerfEntry, ...]:
    """The entries of the newest adapter version per record version, in history order."""
    newest: dict[tuple[str, str], int] = {}
    for entry in entries:
        version = (entry.source.record_path, entry.source.record_sha256)
        newest[version] = max(newest.get(version, entry.source.adapter_version), entry.source.adapter_version)
    return tuple(entry for entry in entries if entry.source.adapter_version == newest[(entry.source.record_path, entry.source.record_sha256)])


def _instant(entry: PerfEntry) -> datetime:
    when = entry.when
    return _UNTIMED if when is None else when


def _latest_rebaselines(rebaselines: Sequence[Rebaseline]) -> dict[Key, Rebaseline]:
    latest: dict[Key, tuple[tuple[datetime, int], Rebaseline]] = {}
    for index, rebaseline in enumerate(rebaselines):
        order = (parse_instant(rebaseline.recorded_at), index)
        current = latest.get(rebaseline.key)
        if current is None or order >= current[0]:
            latest[rebaseline.key] = (order, rebaseline)
    return {key: rebaseline for key, (_order, rebaseline) in latest.items()}


def baselines(history: History, exclude_records: AbstractSet[str] = frozenset()) -> dict[Key, PerfEntry]:
    """The baseline per (workload, stage, machine): the best quiet current entry, from the latest rebaseline on.

    Args:
        history: The history.
        exclude_records: ``record_sha256`` values whose entries are left out (the records under check).

    Raises:
        RebaselineTargetError: a rebaseline names an entry the history does not hold.
    """
    entries = current_entries(history.entries)
    windows = {key: _instant(rebaseline.target(history.entries)) for key, rebaseline in _latest_rebaselines(history.rebaselines).items()}
    best: dict[Key, PerfEntry] = {}
    for entry in entries:
        if not entry.quiet or entry.source.record_sha256 in exclude_records:
            continue
        start = windows.get(entry.key)
        if start is not None and (entry.when is None or entry.when < start):
            continue
        current = best.get(entry.key)
        # The smallest minimum wins; among equal minima, the later instant, then the later line.
        if (
            current is None
            or entry.gated_clock.minimum < current.gated_clock.minimum
            or (entry.gated_clock.minimum == current.gated_clock.minimum and _instant(entry) >= _instant(current))
        ):
            best[entry.key] = entry
    return best


@dataclass(frozen=True)
class Comparison:
    """One entry against its key's baseline.

    Attributes:
        entry: The compared entry.
        baseline: Its key's baseline, or None when the key has no quiet entry.
        over_threshold: Whether its gated minimum exceeds ``REGRESSION_RATIO_MAX`` times the baseline's.
    """

    entry: PerfEntry
    baseline: Optional[PerfEntry]
    over_threshold: bool

    @property
    def regressed(self) -> bool:
        """Over the threshold and quiet: a failing verdict."""
        return self.over_threshold and self.entry.quiet

    @property
    def ratio(self) -> Optional[float]:
        """Gated minimum over the baseline's (reporting), None without a baseline or with a zero baseline."""
        if self.baseline is None or self.baseline.gated_clock.minimum == 0.0:
            return None
        return self.entry.gated_clock.minimum / self.baseline.gated_clock.minimum

    def describe(self) -> str:
        workload, stage, machine = self.entry.key
        state = "quiet" if self.entry.quiet else f"indicative ({self.entry.indicative_reason})"
        head = f"{workload} / {stage} [{machine}] {self.entry.gated_metric} {self.entry.gated_clock.minimum:.4f} s, {state}"
        if self.baseline is None:
            return f"{head}: no quiet baseline"
        ratio = "n/a" if self.ratio is None else f"{self.ratio:.3f}x"
        verdict = "REGRESSION" if self.regressed else "over threshold (indicative, no verdict)" if self.over_threshold else "ok"
        commit = self.baseline.provenance.commit or self.baseline.source.committed_in
        return f"{head}: baseline {self.baseline.gated_clock.minimum:.4f} s ({commit}, {self.baseline.when_text}), {ratio}: {verdict}"


def compare(entry: PerfEntry, baseline: Optional[PerfEntry]) -> Comparison:
    """Compare ``entry`` with ``baseline`` on the gated clock's minimum.

    Raises:
        GatedMetricMismatchError: one gates on CPU and the other on wall.
        ClockStatisticMismatchError: their gated clocks are different statistics.
    """
    if baseline is None:
        return Comparison(entry, None, False)
    if entry.gated_metric != baseline.gated_metric:
        raise GatedMetricMismatchError(f"{entry.key}: the entry gates on {entry.gated_metric}, its baseline on {baseline.gated_metric}.")
    current, reference = entry.gated_clock, baseline.gated_clock
    if _COMPARABLE[current.statistic] != _COMPARABLE[reference.statistic]:
        raise ClockStatisticMismatchError(
            f"{entry.key}: the entry's {entry.gated_metric} is a {_COMPARABLE[current.statistic]}, its baseline's a {_COMPARABLE[reference.statistic]}."
        )
    return Comparison(entry, baseline, current.minimum > REGRESSION_RATIO_MAX * reference.minimum)


def check(history: History, entries: Iterable[PerfEntry]) -> tuple[Comparison, ...]:
    """Compare every entry with its baseline in ``history``, leaving the entries' own records out of the baseline.

    Raises:
        PerformanceRegressionError: a quiet entry is over the threshold (all such entries listed).
        GatedMetricMismatchError, ClockStatisticMismatchError: an entry is not comparable with its baseline.
    """
    checked = tuple(entries)
    reference = baselines(history, frozenset(entry.source.record_sha256 for entry in checked))
    comparisons = tuple(compare(entry, reference.get(entry.key)) for entry in checked)
    regressions = [comparison for comparison in comparisons if comparison.regressed]
    if regressions:
        raise PerformanceRegressionError(regressions)
    return comparisons
