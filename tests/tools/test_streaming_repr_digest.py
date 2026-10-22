"""The streamed ``repr`` digest: byte-identical to ``sha256(repr(value).encode())``, without the string."""

from __future__ import annotations

import dataclasses
import enum
import hashlib
import tracemalloc
from collections import namedtuple
from dataclasses import dataclass
from dataclasses import field
from fractions import Fraction
from typing import Any
from typing import ClassVar

import pytest
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from benchmarks.tools import held_declared_clearance_differential as differential
from tools import held_path_digest
from tools.streaming_repr_digest import MEMORY_ADDRESS
from tools.streaming_repr_digest import AddressGuard
from tools.streaming_repr_digest import ReprCarriesMemoryAddressError
from tools.streaming_repr_digest import repr_sha256

# Number of small dataclasses in the peak-memory case: their ``repr`` is about 9 MB.
_LARGE_COUNT = 200_000
# Ceiling on the traced peak while digesting the large case, bytes: the streamed
# digest holds one chunk of text (well under a megabyte) at a time, the ``repr``
# digest holds the whole string, its parts and its encoding (tens of megabytes).
_STREAMED_PEAK_CEILING = 3 * 1024 * 1024


def _reference(value: object) -> str:
    return hashlib.sha256(repr(value).encode()).hexdigest()


@dataclass(frozen=True)
class _Point:
    x: float
    y: float


@dataclass
class _Mutable:
    items: list[Any]
    label: str = "a"
    hidden: int = field(default=7, repr=False)
    shared: ClassVar[int] = 3


@dataclass(frozen=True)
class _CustomRepr:
    value: int

    def __repr__(self) -> str:
        return f"<custom {self.value}>"


@dataclass(frozen=True)
class _Child(_Point):
    z: float = 0.0


@dataclass(frozen=True, repr=False)
class _InheritedRepr(_Point):
    z: float = 0.0


class _PlainSubclass(_Point):
    pass


@dataclass(frozen=True)
class _Slotted:
    __slots__ = ("a", "b")
    a: Fraction
    b: tuple[Any, ...]


class _Colour(enum.Enum):
    RED = "red"


_Pair = namedtuple("_Pair", "left right")


class _AddressRepr:
    """Default ``object.__repr__``: carries the memory address."""


def _local_dataclass() -> type:
    @dataclass(frozen=True)
    class Local:
        value: object

    return Local


_Local = _local_dataclass()

_LEAVES = st.one_of(
    st.none(),
    st.booleans(),
    st.integers(),
    st.floats(allow_nan=True, allow_infinity=True),
    st.fractions(),
    st.text(max_size=8),
    st.sampled_from([_Colour.RED, -0.0, 1e-300, 2.0**-1074, b"\x00\xff", "é∂\n'\""]),
)


def _compose(children: st.SearchStrategy[Any]) -> st.SearchStrategy[Any]:
    return st.one_of(
        st.lists(children, max_size=4),
        st.lists(children, max_size=4).map(tuple),
        st.tuples(children),
        st.builds(_Point, children, children),
        st.builds(_Mutable, st.lists(children, max_size=3), st.text(max_size=4), st.integers()),
        st.builds(_CustomRepr, st.integers()),
        st.builds(_Child, children, children, children),
        st.builds(_InheritedRepr, children, children, children),
        st.builds(_PlainSubclass, children, children),
        st.builds(_Slotted, st.fractions(), st.lists(children, max_size=3).map(tuple)),
        st.builds(_Pair, children, children),
        st.builds(_Local, children),
        st.dictionaries(st.integers(), children, max_size=3),
        st.frozensets(st.integers(), max_size=3),
    )


_VALUES = st.recursive(_LEAVES, _compose, max_leaves=40)


@settings(max_examples=400, deadline=None)
@given(_VALUES)
def test_the_streamed_digest_is_the_sha256_of_the_repr(value: object) -> None:
    assert repr_sha256(value) == _reference(value)


@settings(max_examples=100, deadline=None)
@given(_VALUES, st.integers(min_value=1, max_value=64))
def test_the_chunk_size_never_changes_the_digest(value: object, chunk_parts: int) -> None:
    assert repr_sha256(value, chunk_parts=chunk_parts) == _reference(value)


@pytest.mark.parametrize(
    "value",
    [
        (),
        (1,),
        [],
        [()],
        ((),),
        _Point(1.0, (2,)),
        _Child(1, 2),
        _InheritedRepr(1, 2),
        _PlainSubclass(1, 2),
        _Mutable([1, [2]]),
        _Local(_Local(None)),
    ],
)
def test_structural_corner_cases(value: object) -> None:
    assert repr_sha256(value) == _reference(value)


def test_self_reference_prints_as_repr_does() -> None:
    loop: list[Any] = [1]
    loop.append(loop)
    holder = _Mutable([])
    holder.items.append(holder)
    for value in (loop, (loop,), holder, [holder, loop]):
        assert repr_sha256(value) == _reference(value)


@pytest.mark.parametrize("chunk_parts", [1, 2, 3, 5, 8, 13, 1 << 20])
def test_a_memory_address_is_refused_at_every_chunk_boundary(chunk_parts: int) -> None:
    with pytest.raises(ReprCarriesMemoryAddressError):
        repr_sha256((_Point(0.0, 0.0), [_AddressRepr()]), chunk_parts=chunk_parts)


@pytest.mark.parametrize("chunk_parts", [1, 2, 3, 7, 1 << 20])
def test_an_address_split_across_leaves_is_refused(chunk_parts: int) -> None:
    @dataclass(frozen=True)
    class Text:
        head: object
        tail: object

    class Raw:
        def __init__(self, text: str) -> None:
            self.text = text

        def __repr__(self) -> str:
            return self.text

    value = [Raw("<thing at 0x1f"), Raw("a>")]
    assert " at 0x1f, a>" in repr(value)  # the separator breaks the match in the full text too
    assert repr_sha256(value, chunk_parts=chunk_parts) == _reference(value)
    joined = (Raw("x at 0x"), Raw("ab>"))
    with pytest.raises(ReprCarriesMemoryAddressError):
        repr_sha256(Raw(repr(joined).replace(", ", "")), chunk_parts=chunk_parts)
    with pytest.raises(ReprCarriesMemoryAddressError):
        repr_sha256(Text(Raw("x"), Raw("<y at 0xAB>")), chunk_parts=chunk_parts)


def _traced_peak(digest: Any, value: object) -> int:
    tracemalloc.start()
    try:
        tracemalloc.reset_peak()
        before = tracemalloc.get_traced_memory()[0]
        digest(value)
        return tracemalloc.get_traced_memory()[1] - before
    finally:
        tracemalloc.stop()


def _large() -> list[_Point]:
    return [_Point(index * 0.1, (index, Fraction(index, 7))) for index in range(_LARGE_COUNT)]


@pytest.mark.parametrize(
    "digest",
    [repr_sha256, held_path_digest.canonical_digest, differential._digest],
    ids=["streaming", "held_path_digest", "differential"],
)
def test_digesting_never_materialises_the_repr(digest: Any) -> None:
    value = _large()
    assert digest(value) == _reference(value)
    assert _traced_peak(digest, value) < _STREAMED_PEAK_CEILING


def test_the_digest_tools_keep_their_own_error_model() -> None:
    with pytest.raises(held_path_digest.NonCanonicalReprError):
        held_path_digest.canonical_digest([_AddressRepr()])
    with pytest.raises(differential.NonCanonicalReprError):
        differential._digest([_AddressRepr()])


def test_dataclass_fields_order_and_repr_flags_are_honoured() -> None:
    value = _Mutable([1.5], "b", 9)
    assert [f.name for f in dataclasses.fields(value) if f.repr] == ["items", "label"]
    assert repr_sha256(value) == _reference(value)


_ADDRESS_ALPHABET = st.text(alphabet=" at0xAbF9>z", max_size=40)


@settings(max_examples=500, deadline=None)
@given(_ADDRESS_ALPHABET, st.lists(st.integers(min_value=0, max_value=40), max_size=6))
def test_the_address_guard_finds_a_match_exactly_when_the_whole_text_has_one(text: str, cuts: list[int]) -> None:
    bounds = sorted({0, len(text), *(min(cut, len(text)) for cut in cuts)})
    guard = AddressGuard()
    refused = False
    try:
        for start, end in zip(bounds, bounds[1:]):
            guard.scan(text[start:end])
    except ReprCarriesMemoryAddressError:
        refused = True
    assert refused == (MEMORY_ADDRESS.search(text) is not None)
