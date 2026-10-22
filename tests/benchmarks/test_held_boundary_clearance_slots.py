"""``BoundaryClearanceSite`` carries no per-instance ``__dict__``: its derived values are slot-backed or recomputed.

One site is built per distance-admissible hypothesis (120 k on Monstera); a
``cached_property`` gave each one a ``__dict__``. Required: no ``__dict__``,
the exact values cached in slots on first read, the reporting doubles equal to
the correctly rounded exact values, frozen fields, the generated ``repr``
unchanged, and value equality, hashing, copying and pickling unaffected by
which derived values have been read.
"""

from __future__ import annotations

import copy
import dataclasses
import pickle
from fractions import Fraction

import pytest

from compas_cgal.held_boundary_clearance import BoundaryClearanceSite
from compas_cgal.held_boundary_clearance import _site

# A native record: segment 3, parameter 0x1/0x3, footpoint (0x5, -0x7) and
# squared distance 0xb on the grid exponent -2 over L = 3.
_RECORD = (3, "1", "3", "5", "-7", "b", -2)


def _fresh() -> BoundaryClearanceSite:
    return _site(_RECORD)


def test_a_site_has_no_instance_dict() -> None:
    assert not hasattr(_fresh(), "__dict__")


def test_the_exact_values_are_cached_once_read() -> None:
    site = _fresh()
    assert site.footpoint == (Fraction(5, 12), Fraction(-7, 12))
    assert site.squared_distance == Fraction(11, 48)
    assert site.footpoint is site.footpoint
    assert site.squared_distance is site.squared_distance


def test_the_reporting_doubles_are_the_rounded_exact_values() -> None:
    site = _fresh()
    assert site.footpoint_doubles == (float(Fraction(5, 12)), float(Fraction(-7, 12)))
    assert site.squared_distance_double == float(Fraction(11, 48))


def test_the_fields_stay_frozen() -> None:
    with pytest.raises(dataclasses.FrozenInstanceError):
        _fresh().segment = 4  # type: ignore[misc]


def test_the_repr_prints_the_fields_only() -> None:
    site = _fresh()
    site.footpoint  # noqa: B018 - fill the slot cache
    assert repr(site) == "BoundaryClearanceSite(segment=3, parameter=Fraction(1, 3), foot_x=5, foot_y=-7, distance=11, length_squared=3, exponent=-2)"


def test_equality_and_hashing_ignore_the_cache_state() -> None:
    read, unread = _fresh(), _fresh()
    read.footpoint  # noqa: B018 - fill the slot cache
    read.squared_distance  # noqa: B018
    assert read == unread and hash(read) == hash(unread)
    scaled = BoundaryClearanceSite(3, Fraction(1, 3), 10, -14, 22, 6, -2)
    assert scaled == unread and hash(scaled) == hash(unread)


@pytest.mark.parametrize("read_first", [False, True])
def test_copies_and_pickles_round_trip(read_first: bool) -> None:
    site = _fresh()
    if read_first:
        site.footpoint  # noqa: B018 - fill the slot cache
    for twin in (copy.copy(site), copy.deepcopy(site), pickle.loads(pickle.dumps(site))):
        assert twin == site and repr(twin) == repr(site)
        assert twin.footpoint == site.footpoint and twin.squared_distance == site.squared_distance
