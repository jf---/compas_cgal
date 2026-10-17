"""A HELIX entry arms the entry-cut window like a plunge (Task 1.2, GAP 3).

`entry_cut_indices` names the first material contact after each entry, so the
post-entry engagement statistics can set it aside. A helix is the entry; if it
did not arm, the first pencil member would count as post-entry engagement.
"""

from __future__ import annotations

import dataclasses

import pytest

from benchmarks.errors import UnmeasurableOperationLengthError
from benchmarks.pathmetrics import entry_cut_indices
from benchmarks.pathmetrics import path_length
from tests.benchmarks.helix_fixtures import emitted
from tests.benchmarks.helix_fixtures import first_member
from tests.benchmarks.helix_fixtures import make_helix
from tests.benchmarks.helix_fixtures import result


def test_a_helix_arms_the_entry_cut_window() -> None:
    helix = make_helix()
    assert entry_cut_indices(result([emitted(helix), first_member(helix), first_member(helix)])) == frozenset({1})


def test_path_length_counts_the_helix_by_its_own_length() -> None:
    helix = make_helix()
    assert path_length(result([emitted(helix)])) == float(helix.length())


def test_path_length_refuses_geometry_outside_the_union_by_name() -> None:
    stranger = dataclasses.replace(emitted(make_helix()), geometry=object())
    with pytest.raises(UnmeasurableOperationLengthError):
        path_length(result([stranger]))
