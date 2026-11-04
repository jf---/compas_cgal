"""Circle-view probe: per sampled circle, the max-run engagement on the actual stock S_t and on the circle view C_t.

Port of the 2026-11-04 design probe (``build/actual-stock-design/probe_circles_only.py``)
onto the oracle's typed views. Both sides are the last-16 certificate's exact
FULL-window stock on the part polygon P (``Last16Oracle.view(i, None)``):
``SweepView.ACTUAL`` (every prior circle and connector not proven exactly to
miss the window) and ``SweepView.CIRCLES`` (the same without connectors). C_t
contains S_t, so the circle view's max run is >= the actual one at every
station; the gap is how much looser the circle view is (design risk 4: up to
21 degrees on the samples). The verdict counts are the D4 statistic on the
sample; the D4 gate itself is ``benchmarks.held_actual_stock_oracle``.

Degrees are REPORTING (``to_double`` of the engine's max run); the exceeded
flags are the engine's exact decisions.

Usage (``build/pyrun.sh -m tools.actual_stock_circle_view_probe ...``):

    WORKLOAD [--stride K] [--out DIR]     # -> DIR/<workload>.circle_view.stride<K>.{json,md}
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from collections.abc import Mapping
from collections.abc import Sequence
from dataclasses import asdict
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from benchmarks.held_actual_stock_last16 import Last16Oracle
from benchmarks.held_actual_stock_last16 import SweepView
from benchmarks.held_actual_stock_oracle import differential_workload
from benchmarks.performance_record import REPOSITORY
from tools.actual_stock_profile import ReplayWorkload

SCHEMA = "actual-stock-circle-view-probe/1"
DEFAULT_OUT = REPOSITORY / "build" / "actual-stock-u0"
# The design probe's strides: 51 / 55 / 74 samples on Figure 5 / skis / Monstera.
DEFAULT_STRIDE = 40
# Quantiles of the gap table.
P50, P90 = 0.5, 0.9


class CircleViewProbeError(Exception):
    """Base class of the circle-view probe failures."""


class MalformedProbeRecordError(CircleViewProbeError):
    """A probe record breaks its schema."""


@dataclass(frozen=True)
class ProbeRow:
    """One sampled circle on its two full-window stocks."""

    circle: int
    sweeps: int  # sweeps of the actual full window (circles + connectors)
    circles: int  # of which circles: the circle view's sweeps
    actual_deg: float
    circle_view_deg: float
    actual_exceeded: bool
    circle_view_exceeded: bool


def _quantile(values: Sequence[float], q: float) -> float:
    ordered = sorted(values)
    return ordered[min(len(ordered) - 1, int(q * (len(ordered) - 1)))] if ordered else 0.0


@dataclass(frozen=True)
class ProbeRecord:
    """The probe over one workload."""

    workload: str
    stride: int
    rows: tuple[ProbeRow, ...]
    cpu_s: float
    load: tuple[float, float, float]

    @property
    def gaps_deg(self) -> tuple[float, ...]:
        return tuple(r.circle_view_deg - r.actual_deg for r in self.rows)

    @property
    def gap_max_deg(self) -> float:
        return max(self.gaps_deg, default=0.0)

    @property
    def circle_view_only_refusals(self) -> tuple[int, ...]:
        """Sampled circles exceeded on C_t and not on S_t (the D4 statistic on the sample)."""
        return tuple(r.circle for r in self.rows if r.circle_view_exceeded and not r.actual_exceeded)

    def to_json(self) -> dict[str, object]:
        return {"schema": SCHEMA, "workload": self.workload, "stride": self.stride, "rows": [asdict(r) for r in self.rows], "cpu_s": self.cpu_s, "load": list(self.load)}

    @classmethod
    def from_json(cls, raw: Mapping[str, object]) -> ProbeRecord:
        """Raises:
        MalformedProbeRecordError: wrong schema or missing / ill-typed fields.
        """
        if raw.get("schema") != SCHEMA:
            raise MalformedProbeRecordError(f"schema {raw.get('schema')!r}, expected {SCHEMA!r}")
        try:
            rows = raw["rows"]
            load = raw["load"]
            if not isinstance(rows, list) or not isinstance(load, list) or len(load) != 3:
                raise MalformedProbeRecordError("rows / load are not lists of the recorded shape")
            return cls(
                workload=str(raw["workload"]),
                stride=int(str(raw["stride"])),
                rows=tuple(ProbeRow(**row) for row in rows),
                cpu_s=float(str(raw["cpu_s"])),
                load=(float(load[0]), float(load[1]), float(load[2])),
            )
        except (KeyError, TypeError) as error:
            raise MalformedProbeRecordError(f"probe record: {error}") from error


def _max_run_deg(oracle: Last16Oracle, index: int) -> tuple[int, float, bool]:
    view, stock = oracle.view(index, None)
    samples = oracle.workload.prefixes(stock, index)
    return view.sweeps, math.degrees(max(run for _start, run, _e in samples)), any(e for *_x, e in samples)


def probe(workload: ReplayWorkload, stride: int) -> ProbeRecord:
    """Every ``stride``-th circle of ``workload`` on its actual and circle-view full windows."""
    start = time.process_time()
    actual, circle_view = Last16Oracle(workload), Last16Oracle(workload, sweeps=SweepView.CIRCLES)
    rows = []
    for index in range(0, len(workload.circles), stride):
        sweeps, actual_deg, actual_exceeded = _max_run_deg(actual, index)
        circles, circle_deg, circle_exceeded = _max_run_deg(circle_view, index)
        rows.append(ProbeRow(index, sweeps, circles, actual_deg, circle_deg, actual_exceeded, circle_exceeded))
    return ProbeRecord(workload.name, stride, tuple(rows), time.process_time() - start, os.getloadavg())


def render_markdown(record: ProbeRecord) -> str:
    gaps = record.gaps_deg
    later = [r.circle_view_deg for r in record.rows if r.circle > 0]
    return (
        f"`{record.workload}` circle-view probe, every {record.stride}th circle: n {len(record.rows)}, "
        f"exceeded actual {sum(r.actual_exceeded for r in record.rows)}, circle view {sum(r.circle_view_exceeded for r in record.rows)}; "
        f"circle-view-only refusals: {list(record.circle_view_only_refusals)}\n\n"
        f"gap (circle view - actual) deg p50 {_quantile(gaps, P50):.2f} p90 {_quantile(gaps, P90):.2f} max {record.gap_max_deg:.2f}; "
        f"max circle-view run after the entry {max(later, default=0.0):.2f} deg; CPU {record.cpu_s:.1f} s, load {record.load[0]:.0f}\n"
    )


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0] if __doc__ else None)
    parser.add_argument("workload")
    parser.add_argument("--stride", type=int, default=DEFAULT_STRIDE)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args(argv)
    record = probe(differential_workload(args.workload), args.stride)
    args.out.mkdir(parents=True, exist_ok=True)
    stem = f"{record.workload}.circle_view.stride{record.stride}"
    (args.out / f"{stem}.json").write_text(json.dumps(record.to_json()) + "\n", encoding="utf-8")
    text = render_markdown(record)
    (args.out / f"{stem}.md").write_text(text, encoding="utf-8")
    sys.stdout.write(text)


if __name__ == "__main__":
    main()
