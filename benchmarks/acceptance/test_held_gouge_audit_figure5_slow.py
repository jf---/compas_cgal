"""Slow-tier acceptance: the CURRENT published Figure 5 path gouges the analytic pocket.

Not collected by pytest (``testpaths = ["tests"]``); run deliberately with
``python -m pytest benchmarks/acceptance/test_held_gouge_audit_figure5_slow.py``.
It generates the path live (about 25 s) and audits it exactly (under a second).

This test records a KNOWN DEFECT, not a goal: Omega is offset from the
polygonal projection, so at concave arcs the chords lie in material and every
circle tangent to Omega gouges by up to the chord sagitta. When Omega becomes
the exact disk erosion, this test must fail and be replaced by
``require_gouge_free`` passing.

The counts are exact decisions. The 2026-10-17 Decimal audit reported 214 / 598
circles and 138 / 533 connectors (> 1e-9 mm / all). The exact audit agrees on
214 and 138 and on every circle's depth; it decides 595, not 598, gouging
circles because three of the Decimal "gouges" (circles 1843, 1844, 1912, margin
1e-79 mm against line 27) are exact ties with ``dist == r``: touching, allowed.
The 381 circle and 395 connector gouges at or below 1e-9 mm are against lines
and are exact too: Omega's straight edges sit a rounding amount (>= 2e-17 mm)
closer than ``r`` to the wall.

Re-pinned 2026-10-21 to the current generator, whose path closes with the seam
connector (last circle back to the first): 2010 circles and 2010 connectors;
592 / 211 gouging circles and 531 / 136 gouging connectors (exact / > 1e-9 mm),
78 and 38 touching. The seam connector itself is clear. The counts above are
the 2026-10-17 path (2013 circles, open, 2012 connectors).

Re-pinned 2026-10-22 (user go) after the exact NOTHING pass joined generation:
it removes circles 980, 987, 995, 1017, 1036 and 1045, each a gouging circle
deeper than the 1e-9 mm split, so 2004 circles and 2004 connectors; 586 / 205
gouging circles (592 - 6, 211 - 6), 78 touching; the connectors' gouges are
unchanged (531 / 136, 38 touching).
"""

from __future__ import annotations

import pytest

from benchmarks.held_gouge_audit import HeldPathGougesBoundaryError
from benchmarks.held_gouge_audit import require_gouge_free
from benchmarks.held_gouge_workloads import audit_generated
from benchmarks.held_gouge_workloads import generate


def test_published_figure5_path_gouges_concave_arcs_by_the_chord_sagitta() -> None:
    report = audit_generated(generate("figure5_refined"))
    assert (report.circles, report.connectors) == (2004, 2004)
    assert (len(report.gouging_circles), report.above_split("circle"), len(report.touching_circles)) == (586, 205, 78)
    assert (len(report.gouging_connectors), report.above_split("connector"), len(report.touching_connectors)) == (531, 136, 38)
    geometric = [item for item in report.gouges if item.depth_mm > 1e-9]
    assert {item.primitive_kind for item in geometric} == {"concave-arc"}
    assert {item.primitive_kind for item in report.gouges if item not in geometric} == {"line"}
    with pytest.raises(HeldPathGougesBoundaryError, match="586 circles and 531 connectors gouge"):
        require_gouge_free(report)
