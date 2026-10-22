"""Vector-unit factories validate once and keep the pre-change values and errors (lever 2, step 3)."""

from __future__ import annotations

import math
from collections.abc import Callable
from collections.abc import Sequence
from fractions import Fraction
from typing import Any

import numpy as np
import pytest
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from compas_cgal import units
from compas_cgal.errors import InvalidUnitValueError
from compas_cgal.units import Direction3
from compas_cgal.units import Point2
from compas_cgal.units import Point3
from compas_cgal.units import Vector2
from compas_cgal.units import WorldXY


def _oracle_finite(value: object, name: str) -> float:
    """Verbatim pre-step-3 ``units._finite``."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise InvalidUnitValueError(f"{name} must be a finite real number.")
    try:
        numeric = float(value)
    except OverflowError:
        raise InvalidUnitValueError(f"{name} exceeds the binary64 range.") from None
    if not math.isfinite(numeric):
        raise InvalidUnitValueError(f"{name} must be finite.")
    return numeric


def _oracle_components(args: tuple[object, ...], name: str) -> tuple[float, float]:
    """Verbatim pre-step-3 ``units._components``."""
    if len(args) == 1 and isinstance(args[0], Sequence) and not isinstance(args[0], (str, bytes)):
        values = args[0]
        if len(values) != 2:
            raise InvalidUnitValueError(f"{name} requires exactly two components.")
        return _oracle_finite(values[0], f"{name}.x"), _oracle_finite(values[1], f"{name}.y")
    if len(args) == 2:
        return _oracle_finite(args[0], f"{name}.x"), _oracle_finite(args[1], f"{name}.y")
    raise InvalidUnitValueError(f"{name} requires two scalars or one two-component sequence.")


def _oracle_components3(args: tuple[object, ...], name: str) -> tuple[float, float, float]:
    """Verbatim pre-step-3 ``units._components3``."""
    if len(args) == 1 and isinstance(args[0], Sequence) and not isinstance(args[0], (str, bytes)):
        values = args[0]
        if len(values) != 3:
            raise InvalidUnitValueError(f"{name} requires exactly three components.")
        return (
            _oracle_finite(values[0], f"{name}.x"),
            _oracle_finite(values[1], f"{name}.y"),
            _oracle_finite(values[2], f"{name}.z"),
        )
    if len(args) == 3:
        return (
            _oracle_finite(args[0], f"{name}.x"),
            _oracle_finite(args[1], f"{name}.y"),
            _oracle_finite(args[2], f"{name}.z"),
        )
    raise InvalidUnitValueError(f"{name} requires three scalars or one three-component sequence.")


def _oracle_build2(cls: Any, *args: object) -> object:
    """Verbatim pre-step-3 ``Point2.build`` / ``Vector2.build`` (with the oracle ``__post_init__`` checks)."""
    x, y = _oracle_components(args, cls.__name__)
    return cls(x, y)


def _oracle_build3(cls: Any, *args: object) -> object:
    """Verbatim pre-step-3 ``Point3.build`` / ``Direction3.build``."""
    x, y, z = _oracle_components3(args, cls.__name__)
    return cls(x, y, z)


class SemanticPoint(Point2[WorldXY]):
    """A subclass: its error messages carry its own name through ``build``."""


class SemanticVector(Vector2[WorldXY]):
    """A subclass of the other two-component unit."""


def _outcome(call: Callable[[], object]) -> tuple[object, ...]:
    try:
        value = call()
    except InvalidUnitValueError as error:
        return ("error", type(error), str(error))
    components = tuple(getattr(value, name) for name in ("x", "y", "z") if hasattr(value, name))
    return ("value", type(value), repr(value), tuple(type(component) for component in components))


_SCALARS = st.one_of(
    st.floats(allow_nan=True, allow_infinity=True),
    st.floats(allow_nan=True, allow_infinity=True).map(np.float64),
    st.integers(min_value=-(2**70), max_value=2**70),
    st.integers(min_value=10**309, max_value=10**310),
    st.booleans(),
    st.fractions(max_denominator=8),
    st.just(None),
    st.just("1.0"),
    st.just(b"1"),
)

_TWO_CLASSES = (Point2, Vector2, SemanticPoint, SemanticVector, Point2[WorldXY])
_THREE_CLASSES = (Point3, Direction3, Point3[WorldXY])


@settings(max_examples=600, deadline=None)
@given(
    cls=st.sampled_from(_TWO_CLASSES),
    args=st.one_of(
        st.tuples(_SCALARS, _SCALARS),
        st.tuples(st.lists(_SCALARS, min_size=0, max_size=3)),
        st.tuples(st.lists(_SCALARS, min_size=2, max_size=2).map(tuple)),
        st.tuples(_SCALARS),
        st.tuples(_SCALARS, _SCALARS, _SCALARS),
        st.tuples(st.just("ab")),
    ),
)
def test_two_component_build_matches_oracle(cls: Any, args: tuple[object, ...]) -> None:
    assert _outcome(lambda: cls.build(*args)) == _outcome(lambda: _oracle_build2(cls, *args))


@settings(max_examples=400, deadline=None)
@given(
    cls=st.sampled_from(_THREE_CLASSES),
    args=st.one_of(
        st.tuples(_SCALARS, _SCALARS, _SCALARS),
        st.tuples(st.lists(_SCALARS, min_size=0, max_size=4)),
        st.tuples(_SCALARS, _SCALARS),
    ),
)
def test_three_component_build_matches_oracle(cls: Any, args: tuple[object, ...]) -> None:
    assert _outcome(lambda: cls.build(*args)) == _outcome(lambda: _oracle_build3(cls, *args))


@settings(max_examples=400, deadline=None)
@given(value=_SCALARS)
def test_finite_matches_oracle(value: object) -> None:
    assert _outcome(lambda: units._finite(value, "v")) == _outcome(lambda: _oracle_finite(value, "v"))


@pytest.mark.parametrize("cls", [Point2, Vector2, Point2[WorldXY]])
def test_two_float_build_validates_once(cls: Any, monkeypatch: pytest.MonkeyPatch) -> None:
    """Finite doubles go straight to the constructor; its ``__post_init__`` is the one validation."""

    def _refuse(*_args: object) -> None:
        raise AssertionError("build re-validated through _components")

    monkeypatch.setattr(units, "_components", _refuse)

    assert cls.build(1.5, -2.25) == cls(1.5, -2.25)


@pytest.mark.parametrize("cls", [Point3, Direction3])
def test_three_float_build_validates_once(cls: Any, monkeypatch: pytest.MonkeyPatch) -> None:
    def _refuse(*_args: object) -> None:
        raise AssertionError("build re-validated through _components3")

    monkeypatch.setattr(units, "_components3", _refuse)

    assert cls.build(1.5, -2.25, 0.5) == cls(1.5, -2.25, 0.5)


def test_fraction_component_is_still_refused() -> None:
    with pytest.raises(InvalidUnitValueError, match=r"^Point2\.x must be a finite real number\.$"):
        Point2.build(Fraction(1, 3), 0.0)  # type: ignore[call-overload]


def _base_name(cls: Any) -> str:
    """The class whose ``__post_init__`` names the components (subclasses inherit it)."""
    origin = getattr(cls, "__origin__", cls)
    return str(next(base.__name__ for base in origin.__mro__ if base.__name__ in {"Point2", "Vector2", "Point3", "Direction3"}))


def _oracle_construct(cls: Any, *values: object) -> object:
    """Pre-step-3 ``__post_init__``: every component through ``_finite``, in order, with the base class's name."""
    name = _base_name(cls)
    return cls(*(_oracle_finite(value, f"{name}.{axis}") for value, axis in zip(values, "xyz")))


@settings(max_examples=600, deadline=None)
@given(cls=st.sampled_from(_TWO_CLASSES), x=_SCALARS, y=_SCALARS)
def test_two_component_constructor_matches_oracle(cls: Any, x: object, y: object) -> None:
    assert _outcome(lambda: cls(x, y)) == _outcome(lambda: _oracle_construct(cls, x, y))


@settings(max_examples=400, deadline=None)
@given(cls=st.sampled_from(_THREE_CLASSES), x=_SCALARS, y=_SCALARS, z=_SCALARS)
def test_three_component_constructor_matches_oracle(cls: Any, x: object, y: object, z: object) -> None:
    assert _outcome(lambda: cls(x, y, z)) == _outcome(lambda: _oracle_construct(cls, x, y, z))
