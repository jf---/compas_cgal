from __future__ import annotations

import dataclasses

from benchmarks.families.analytic import rectangle
from benchmarks.instrument import probe_digits
from benchmarks.instrument import probe_size
from benchmarks.runner import run_corpus
from benchmarks.runner import run_spec
from benchmarks.spec import PocketSpec
from compas_cgal.stock import Stock

# The engagement audit's contract is a cap in (0, pi]; 270 degrees is outside it,
# so this instance is guaranteed to raise inside the timed region. Built by
# `dataclasses.replace` rather than `PocketSpec.build` precisely to bypass the
# factory's validation -- an unmeasurable instance is what the sweep must survive.
UNCERTIFIABLE_CAP_DEG = 270.0


def _pocket() -> PocketSpec:
    """The smallest instance that still produces a many-operation toolpath."""
    return rectangle(width=8.0, height=6.0, tool_diameter=2.0, tea_cap_deg=120.0)


def test_run_spec_populates_both_phases() -> None:
    record = run_spec(_pocket(), collect_digits=True)
    assert record.error is None
    assert record.generate_seconds > 0.0
    assert record.certify_seconds > 0.0
    assert record.operations > 0
    assert record.cut_operations > 0
    assert record.stations > 0
    assert record.max_coordinate_digits > 0


def test_run_spec_reports_diagnostics_from_the_depleted_stock() -> None:
    # Both kernel diagnostics are read after the toolpath has been replayed onto a
    # fresh stock. Reading them on an undepleted stock would pin them to the virgin
    # values forever, making the two fields incapable of showing the growth they
    # exist to expose.
    spec = _pocket()
    virgin = Stock(spec.polygon, list(spec.holes))
    record = run_spec(spec, collect_digits=True)
    assert record.arrangement_vertices_final > probe_size(virgin).vertices
    assert record.max_coordinate_digits > probe_digits(virgin).max_digits


def test_run_spec_without_digits_leaves_the_field_zero() -> None:
    record = run_spec(_pocket(), collect_digits=False)
    assert record.max_coordinate_digits == 0
    # The arrangement size is a pure counter read, so it stays populated even when
    # the exact-evaluation probe is switched off.
    assert record.arrangement_vertices_final > 0


def test_run_corpus_records_a_failure_without_aborting_the_sweep() -> None:
    good = _pocket()
    bad = dataclasses.replace(good, name="uncertifiable_cap", tea_cap_deg=UNCERTIFIABLE_CAP_DEG)
    records = run_corpus([bad, good], collect_digits=False)
    assert [r.name for r in records] == ["uncertifiable_cap", good.name]
    assert records[0].error is not None
    assert records[0].error.startswith("InvalidEngagementCapError")
    assert records[0].certify_seconds == 0.0
    assert records[1].error is None
    assert records[1].operations > 0


def test_run_spec_records_both_cap_columns() -> None:
    """The runner fills both cap columns, and they are not each other.

    `uncertified` is the audit's could-not-prove count; `truly_exceeding` is the
    sampled demonstration that the exact predicate fired. The soundness relation
    between them holds in one direction only, which is exactly why recording a
    single column under either name misreports the other.
    """
    record = run_spec(_pocket(), collect_digits=False)
    assert record.truly_exceeding > 0
    assert record.truly_exceeding <= record.uncertified
