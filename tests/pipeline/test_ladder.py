"""The generic furthest-first ladder (Ruling 1.1-R13, amendments M4/D2): order, skip, no fallback."""

from __future__ import annotations

from compas.geometry import Polygon

from compas_cgal.pipeline.ladder import NoAdmissibleCandidate
from compas_cgal.pipeline.ladder import furthest_first
from compas_cgal.stock import Stock

RECT_12X8 = Polygon([[-6.0, -4.0, 0.0], [6.0, -4.0, 0.0], [6.0, 4.0, 0.0], [-6.0, 4.0, 0.0]])


def test_the_ladder_evaluates_conditions_in_order_and_stops_at_the_first_refusal() -> None:
    calls: list[str] = []

    def radius(_stock: Stock, _origin: int, candidate: int) -> bool:
        calls.append(f"r{candidate}")
        return candidate != 6

    def cap(_stock: Stock, _origin: int, candidate: int) -> bool:
        calls.append(f"c{candidate}")
        return True

    chosen = furthest_first(Stock(RECT_12X8), 0, (6, 5, 4, 3, 2, 1), (radius, cap))
    # 6 is refused by the first condition, so the second is never asked about it;
    # the scan moves on to 5 (it does not stop), which passes both.
    assert chosen == 5
    assert calls == ["r6", "r5", "c5"]


def test_the_ladder_returns_the_furthest_pass_under_a_non_monotone_pattern() -> None:
    # Admissible {2, 5} in a window reaching 6: a bisection on [1, 6] returns 2.
    seen: list[int] = []

    def admissible(_stock: Stock, _origin: int, candidate: int) -> bool:
        seen.append(candidate)
        return candidate in {2, 5}

    assert furthest_first(Stock(RECT_12X8), 0, (6, 5, 4, 3, 2, 1), (admissible,)) == 5
    assert seen[0] == 6


def test_the_ladder_never_takes_a_forced_fallback() -> None:
    def refuse(_stock: Stock, _origin: int, _candidate: int) -> bool:
        return False

    chosen = furthest_first(Stock(RECT_12X8), 0, (3, 2, 1), (refuse,))
    assert chosen == NoAdmissibleCandidate(scanned=3)
