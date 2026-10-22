"""Exact Held centre domain Omega* (docs/held_exact_centre_domain.md); CORE-free.

Mirrors src/held_centre_domain_2_bindings.cpp. Rationals cross as
``(numerator, denominator)`` decimal text, one-root numbers as
``(rational, coefficient, radicand)`` triples of those. No float crosses
except the binary64 query coordinates of ``project_binary64``, injected exactly.
"""

from collections.abc import Sequence
from typing import TypeAlias

_RationalWire: TypeAlias = tuple[str, str]
_OneRootWire: TypeAlias = tuple[_RationalWire, _RationalWire, _RationalWire]
_PointWire: TypeAlias = tuple[_OneRootWire, _OneRootWire]
# (kind, sense, rationals): a line is (0, 0, [sx, sy, ex, ey, nx, ny]); an arc is
# (1, sense, [sx, sy, ex, ey, cx, cy, radius]).
_PrimitiveWire: TypeAlias = tuple[int, int, Sequence[_RationalWire]]
# (circular, support, sense, source, target); support is (a, b, c) of a*x + b*y + c = 0
# for a line, (centre_x, centre_y, radius_squared) for a circle.
_PieceWire: TypeAlias = tuple[bool, tuple[_RationalWire, _RationalWire, _RationalWire], int, _PointWire, _PointWire]
# (ambiguous, kind, piece, squared_distance, x, y); ambiguous is 0 or 1.
_NearestWire: TypeAlias = tuple[int, int, int, _OneRootWire, _OneRootWire, _OneRootWire]
# (status, kind, piece); status 0 projected, 1 ambiguous, 2 beyond gap.
_DecisionWire: TypeAlias = tuple[int, int, int]

class InvalidCentreDomainInputError(ValueError): ...
class BrokenCentreDomainChainError(RuntimeError): ...
class AmbiguousOuterCcbError(RuntimeError): ...
class BoundedFaceWithoutOuterCcbError(RuntimeError): ...
class UnboundedFaceOuterCcbError(RuntimeError): ...

def last_exact_fallbacks() -> int:
    """Queries of the last ``project_binary64`` batch the interval stage could not decide."""

class CentreDomain2:
    @staticmethod
    def build(primitives: Sequence[_PrimitiveWire], tool_radius: _RationalWire) -> CentreDomain2: ...
    def components(self) -> list[list[list[_PieceWire]]]:
        """Per component, the outer cycle (CCW) then its holes (CW), each a list of pieces."""
    @property
    def member_count(self) -> int: ...
    @property
    def construction_seconds(self) -> float: ...
    @property
    def refiled_outer_ccbs(self) -> int:
        """Outer CCBs the Gps difference misfiled and the build re-filed as holes (reporting only)."""
    def locate(self, x: _RationalWire, y: _RationalWire) -> int:
        """``1`` interior, ``0`` boundary, ``-1`` exterior."""
    def nearest(self, x: _RationalWire, y: _RationalWire) -> _NearestWire: ...
    def project_binary64(self, xs: Sequence[float], ys: Sequence[float], gap: _RationalWire) -> list[_DecisionWire]: ...
