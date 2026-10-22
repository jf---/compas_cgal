"""``CentreDomain.unreachable`` is attributed on first access, once, to the value the eager attribution gives.

Generation never reads the attribution, so :meth:`CentreDomain.build` must not
pay for it; every reader still sees the same exact value.
"""

from __future__ import annotations

from fractions import Fraction

import pytest

import benchmarks.held_exact_centre_domain as centre_domain
from benchmarks.held_exact_centre_domain import CentreDomain
from tests.benchmarks.test_held_exact_centre_domain import concave_arc_pocket
from tests.benchmarks.test_held_exact_centre_domain import l_shape_pocket

# Tool radius (mm) of the fixtures.
TOOL_RADIUS = Fraction(1)


@pytest.mark.parametrize("pocket", [concave_arc_pocket, l_shape_pocket])
def test_unreachable_is_attributed_on_first_access_only(monkeypatch: pytest.MonkeyPatch, pocket: object) -> None:
    boundary = pocket()  # type: ignore[operator]
    attribute = centre_domain._attribute_unreachable
    calls: list[int] = []

    def spy(*args: object) -> object:
        calls.append(1)
        return attribute(*args)  # type: ignore[arg-type]

    monkeypatch.setattr(centre_domain, "_attribute_unreachable", spy)
    domain = CentreDomain.build(boundary, TOOL_RADIUS)
    assert not calls
    first = domain.unreachable
    assert len(calls) == 1
    assert domain.unreachable is first
    assert len(calls) == 1
    eager = attribute(boundary, domain.supports, domain.components, domain.native)
    assert first == eager
