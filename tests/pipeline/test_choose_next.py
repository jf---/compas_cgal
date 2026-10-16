"""The next-chain choice checks links lazily in key order and chooses exactly what the eager rule chooses.

The key is ``(not linked, link length, pending position)``, first offered
orientation on ties. Length and position are link geometry; only ``linked``
needs the in-process admissibility, so the candidates are checked in
``(length, position)`` order and the first admissible one is the minimum.
"""

from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

import pytest
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

from compas_cgal.pipeline import generate
from compas_cgal.pipeline.generate import NextChain
from compas_cgal.pipeline.generate import WalkEnd
from compas_cgal.pipeline.linking import NoLink


@dataclass(frozen=True)
class FakeLink:
    """A scripted link: its ordering length and whether its pieces are admissible."""

    length: float
    admissible: bool


@dataclass(frozen=True)
class Offer:
    """One offered orientation: no link at all, or a scripted one."""

    link: FakeLink | None


def _fake_chain(offer: Offer) -> SimpleNamespace:
    chain = SimpleNamespace(first=0, entry_station=SimpleNamespace(clockwise=True), offer=offer)
    chain.loops = [chain]  # `link_between` receives the first loop; the fake reads the offer from it
    return chain


def _eager_choice(pending: list[tuple[SimpleNamespace, ...]]) -> tuple[int, SimpleNamespace, FakeLink | None]:
    """The pre-round-2 rule, verbatim in its key and tie-break: every link checked, strict ``<`` over the key."""
    best: tuple[tuple[bool, float, int], tuple[int, SimpleNamespace, FakeLink | None]] | None = None
    for position, offered in enumerate(pending):
        for chain in offered:
            link = chain.offer.link
            linked = link is not None and link.admissible
            key = (not linked, link.length if link is not None else 0.0, position)
            if best is None or key < best[0]:
                best = (key, (position, chain, link if linked else None))
    assert best is not None
    return best[1]


def _install(monkeypatch: pytest.MonkeyPatch, checked: list[FakeLink]) -> None:
    def scripted_link_between(_current_loop: object, first: object, _current_cw: bool, _cw: bool) -> object:
        offer: Offer = first.offer  # type: ignore[attr-defined]
        return NoLink() if offer.link is None else offer.link

    def link_admissible(_path: object, link: FakeLink, _clockwise: bool) -> bool:
        checked.append(link)
        return link.admissible

    monkeypatch.setattr(generate, "link_between", scripted_link_between)
    monkeypatch.setattr(generate, "_link_admissible", link_admissible)
    monkeypatch.setattr(generate, "_link_length", lambda link: link.length)


def _as_pending(offers: list[list[Offer]]) -> list[tuple[SimpleNamespace, ...]]:
    return [tuple(_fake_chain(offer) for offer in offered) for offered in offers]


LENGTHS = st.sampled_from([0.0, 1.0, 2.5, 2.5, 7.0])  # few values: ties are the point
OFFERS = st.builds(Offer, st.one_of(st.none(), st.builds(FakeLink, LENGTHS, st.booleans())))


@settings(max_examples=400, deadline=None)
@given(st.lists(st.lists(OFFERS, min_size=1, max_size=2), min_size=1, max_size=6))
def test_the_lazy_choice_is_the_eager_choice(offers: list[list[Offer]]) -> None:
    with pytest.MonkeyPatch.context() as monkeypatch:
        checked: list[FakeLink] = []
        _install(monkeypatch, checked)
        pending = _as_pending(offers)
        current = WalkEnd(loop="current-loop", clockwise=True, resume=None)  # type: ignore[arg-type]
        choice = generate.choose_next(generate.PathState.__new__(generate.PathState), current, pending)  # type: ignore[arg-type]
    position, chain, link = _eager_choice(pending)
    assert isinstance(choice, NextChain)
    assert (choice.pending, choice.chain, choice.link) == (position, chain, link)
    assert choice.chain is chain  # the same orientation object, ties included


def test_link_checks_stop_at_the_first_admissible_candidate_in_key_order(monkeypatch: pytest.MonkeyPatch) -> None:
    """Counting fake: the short refused link, then the admissible one; the longer links are never checked."""
    checked: list[FakeLink] = []
    _install(monkeypatch, checked)
    short_refused = FakeLink(1.0, admissible=False)
    winner = FakeLink(2.0, admissible=True)
    pending = _as_pending([[Offer(FakeLink(9.0, True))], [Offer(None), Offer(winner)], [Offer(short_refused)], [Offer(FakeLink(5.0, True))]])
    current = WalkEnd(loop="current-loop", clockwise=True, resume=None)  # type: ignore[arg-type]
    choice = generate.choose_next(generate.PathState.__new__(generate.PathState), current, pending)  # type: ignore[arg-type]
    assert choice.link is winner and choice.pending == 1
    assert checked == [short_refused, winner]


def test_no_admissible_link_checks_every_link_once_and_enters_the_shortest(monkeypatch: pytest.MonkeyPatch) -> None:
    checked: list[FakeLink] = []
    _install(monkeypatch, checked)
    links = [FakeLink(3.0, False), FakeLink(1.0, False)]
    pending = _as_pending([[Offer(links[0])], [Offer(links[1])]])
    current = WalkEnd(loop="current-loop", clockwise=True, resume=None)  # type: ignore[arg-type]
    choice = generate.choose_next(generate.PathState.__new__(generate.PathState), current, pending)  # type: ignore[arg-type]
    assert choice.link is None and choice.pending == 1
    assert sorted(checked, key=lambda link: link.length) == sorted(links, key=lambda link: link.length)
