"""Performance history CLI: ``ingest``, ``backfill``, ``check`` and ``render`` over ``benchmarks/results/perf_history.jsonl``.

- ``ingest <record.json>...``: append the entries of each record as it is on
  disk (``tools.perf_history_ingest.file_entries``); duplicates are skipped.
- ``backfill``: append the entries of every committed version of every timing
  record on this branch, oldest first (``backfill_entries``).
- ``check <record.json>...``: compare each record's entries with the baseline
  (``tools.perf_history_baseline.check``); exits non-zero on a quiet regression.
  It does not append.
- ``rebaseline --entry <entry_sha256> --reason <text>``: accept that quiet
  entry as its key's baseline from now on, though slower (an appended
  ``rebaseline`` line; ``tools.perf_history_rebaseline``).
- ``render``: rewrite ``benchmarks/results/perf_history.md``.

Every record path must be a known timing family; anything else raises.
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from datetime import datetime
from datetime import timezone
from pathlib import Path
from typing import Optional

from tools.perf_history_baseline import PerformanceRegressionError
from tools.perf_history_baseline import check
from tools.perf_history_entry import PerfEntry
from tools.perf_history_ingest import backfill_entries
from tools.perf_history_ingest import file_entries
from tools.perf_history_ingest import record_versions
from tools.perf_history_rebaseline import Rebaseline
from tools.perf_history_render import write_markdown
from tools.perf_history_store import HISTORY_PATH
from tools.perf_history_store import AppendResult
from tools.perf_history_store import append_entries
from tools.perf_history_store import append_rebaseline
from tools.perf_history_store import load_history

# Exit status of ``check`` on a quiet regression.
REGRESSION_EXIT_STATUS = 1


def _report(result: AppendResult) -> None:
    quiet = sum(1 for entry in result.added if entry.quiet)
    print(f"appended {len(result.added)} entries ({quiet} quiet, {len(result.added) - quiet} indicative); {len(result.duplicates)} already present.")
    for entry in result.added:
        if not entry.quiet:
            print(f"  indicative: {entry.measurement.workload} / {entry.measurement.stage}: {entry.indicative_reason}")


def ingest(paths: Sequence[Path], history: Path = HISTORY_PATH) -> AppendResult:
    """Append the entries of every record in ``paths``; all are adapted before any is written."""
    entries: list[PerfEntry] = []
    for path in paths:
        entries.extend(file_entries(path))
    return append_entries(entries, history)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("--history", type=Path, default=HISTORY_PATH)
    commands = parser.add_subparsers(dest="command", required=True)
    ingest_parser = commands.add_parser("ingest", help="append the entries of records on disk")
    ingest_parser.add_argument("records", nargs="+", type=Path)
    commands.add_parser("backfill", help="append every committed version of every timing record")
    check_parser = commands.add_parser("check", help="compare records with the baseline; non-zero exit on a quiet regression")
    check_parser.add_argument("records", nargs="+", type=Path)
    rebaseline_parser = commands.add_parser("rebaseline", help="accept a quiet entry as its key's baseline, though slower")
    rebaseline_parser.add_argument("--entry", required=True, help="the accepted entry's entry_sha256")
    rebaseline_parser.add_argument("--reason", required=True, help="why the slower baseline is accepted")
    commands.add_parser("render", help="rewrite benchmarks/results/perf_history.md")
    arguments = parser.parse_args(argv)
    history_path: Path = arguments.history
    if arguments.command == "ingest":
        _report(ingest(arguments.records, history_path))
    elif arguments.command == "backfill":
        result = backfill_entries(record_versions())
        print(f"read {len(result.timing_versions)} timing-record versions; {len(result.non_timing_versions)} non-timing record versions not ingested:")
        for version in result.non_timing_versions:
            print(f"  {version.rev[:8]} {version.path}")
        _report(append_entries(result.entries, history_path))
    elif arguments.command == "check":
        entries = [entry for path in arguments.records for entry in file_entries(path)]
        try:
            comparisons = check(load_history(history_path), entries)
        except PerformanceRegressionError as error:
            print(error, file=sys.stderr)
            return REGRESSION_EXIT_STATUS
        for comparison in comparisons:
            print(comparison.describe())
    elif arguments.command == "rebaseline":
        accepted = Rebaseline.build(load_history(history_path), arguments.entry, arguments.reason, datetime.now(timezone.utc).isoformat(timespec="seconds"))
        append_rebaseline(accepted, history_path)
        print(f"rebaselined {accepted.key} on {accepted.entry_sha256}: {accepted.reason}")
    else:
        write_markdown(history_path, history_path.with_suffix(".md"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
