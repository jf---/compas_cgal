"""Native placement batch at VERTEX contacts against the Python ``place_local``, by repr.

The generation A/B (57,743 placements on fig5, skis and Monstera) holds no
vertex placement: on an eroded pocket the only kept vertex contacts are probes
exactly on the common normal of a G1 junction. These fixtures decide such
contacts with ``CentreDomain.decide_contacts`` (not hand-built
``DecidedContact``s) and place them in one native batch:

- the L-shape's junction arc: a rational G1 vertex at an arc's source and one at
  a line's source (the contacts the Ω* path keeps);
- a shallow notch whose junction arc cuts the offset top line: one-root convex
  vertices at an arc's source and at a line's source (refused by the path, but
  the placement contract is the same and exercises one-root source points).
"""

from __future__ import annotations

from fractions import Fraction

from benchmarks.held_centre_domain_supports import OffsetCircleSupport
from benchmarks.held_centre_domain_supports import OffsetLineSupport
from benchmarks.held_declared_boundary import DeclaredBoundary
from benchmarks.held_declared_boundary import DeclaredLine
from benchmarks.held_declared_boundary import PythagoreanDirection
from benchmarks.held_exact_centre_domain import CentreDomain
from benchmarks.held_exact_centre_domain import DecidedContact
from benchmarks.held_omega_circles import NoContainedCircleError
from benchmarks.held_omega_circles import place_local
from benchmarks.held_omega_circles import place_local_batch
from benchmarks.held_omega_sites import OmegaLoop
from benchmarks.held_omega_sites import OmegaSite
from tests.benchmarks.test_held_exact_centre_domain import _Pocket
from tests.benchmarks.test_held_exact_centre_domain import _point
from tests.benchmarks.test_held_exact_centre_domain import l_shape_pocket
from tests.benchmarks.test_held_omega_native_placement import GAP
from tests.benchmarks.test_held_omega_native_placement import GUIDE_RADII
from tests.benchmarks.test_held_omega_native_placement import _provenance

F = Fraction
R = F(1)
# Notch side length (mm): binary64, and short enough (< 1/4 with the 15/17,
# 8/17 directions) that the side offsets leave Ω* and the junction arc of
# radius r about the tip cuts the offset top line y = 9 directly.
NOTCH_SIDE = F(17, 128)
# Half-angle parameters of the notch sides: directions (-15/17, -8/17) down, (-15/17, 8/17) up.
NOTCH_DOWN_T, NOTCH_UP_T = F(1, 4), F(-1, 4)
# Top segments either side of the notch (mm): 9.75 + 30/128 + 641/64 = 20.
NOTCH_RIGHT_TOP, NOTCH_LEFT_TOP = F(39, 4), F(641, 64)
# Probe grid around a vertex (mm): binary64 offsets on a 2^-5 lattice within ±0.5.
GRID_STEP = 2.0**-5
GRID_HALF_STEPS = 16
# Offsets (mm) along the L-shape's G1 normals, either side of the vertex.
NORMAL_OFFSETS = (-0.75, -0.5, -0.25, 0.25, 0.5, 0.75)


def shallow_notch_pocket() -> DeclaredBoundary:
    """20 x 10 rectangle whose top edge carries a 1/16 mm deep notch with its tip at (10.1328125, 9.9375)."""
    pocket = _Pocket(_point(0, 0)).line(0, 20).line(1, 10).line(2, NOTCH_RIGHT_TOP)
    for t in (NOTCH_DOWN_T, NOTCH_UP_T):
        side = DeclaredLine.build(pocket.cursor, PythagoreanDirection.build(t, 2), NOTCH_SIDE, len(pocket.primitives))
        pocket.primitives.append(side)
        pocket.cursor = side.end
    return pocket.line(2, NOTCH_LEFT_TOP).line(3, 10).close("shallow_notch")


def _loop(boundary: DeclaredBoundary) -> OmegaLoop:
    return OmegaLoop.build(CentreDomain.build(boundary, R))


def _vertex_contacts(loop: OmegaLoop, probes: list[tuple[float, float]]) -> list[DecidedContact]:
    decisions = loop.domain.decide_contacts(probes, GAP).contacts
    return [decision for decision in decisions if decision is not None and decision.is_vertex]


def _reference(loop: OmegaLoop, contact: DecidedContact, radius: float) -> str:
    try:
        return repr(place_local(OmegaSite.decided(loop, contact), radius, _provenance(0)))
    except NoContainedCircleError:
        return "NoContainedCircleError"


def _batched_one(loop: OmegaLoop, contact: DecidedContact, radius: float) -> str:
    try:
        (circle,) = place_local_batch(loop, [contact], [radius], [_provenance(0)])
    except NoContainedCircleError:
        return "NoContainedCircleError"
    return repr(circle)


def _assert_batch_equals_reference(loop: OmegaLoop, contacts: list[DecidedContact]) -> int:
    """Every contact at every guide radius, singly and (where all place) in one batch; returns the placed count."""
    placed = 0
    for radius in GUIDE_RADII:
        expected = [_reference(loop, contact, radius) for contact in contacts]
        assert [_batched_one(loop, contact, radius) for contact in contacts] == expected, radius
        if "NoContainedCircleError" not in expected:
            batch = place_local_batch(loop, contacts, [radius] * len(contacts), [_provenance(0)] * len(contacts))
            assert [repr(circle) for circle in batch] == expected, radius
        placed += sum(text != "NoContainedCircleError" for text in expected)
    return placed


def _pieces(loop: OmegaLoop, contacts: list[DecidedContact]) -> set[tuple[int, bool]]:
    return {(contact.piece.index, isinstance(contact.piece.support, OffsetCircleSupport)) for contact in contacts}


def test_g1_junction_vertices_place_as_place_local() -> None:
    loop = _loop(l_shape_pocket())
    # Junction arc about (10, 10): source (10, 9) on the normal x = 10, target (9, 10) on y = 10.
    probes = [(10.0, 9.0 + offset) for offset in NORMAL_OFFSETS] + [(9.0 + offset, 10.0) for offset in NORMAL_OFFSETS]
    contacts = _vertex_contacts(loop, probes)
    assert len(contacts) == len(probes)
    kinds = {is_circle for _, is_circle in _pieces(loop, contacts)}
    assert kinds == {True, False}
    assert all(not loop.convex_sources[contact.piece.index] for contact in contacts)
    assert all(contact.piece.source.x.is_rational and contact.piece.source.y.is_rational for contact in contacts)
    assert _assert_batch_equals_reference(loop, contacts) > 0


def test_one_root_vertices_place_as_place_local() -> None:
    loop = _loop(shallow_notch_pocket())
    vertices = [piece for piece in loop.pieces if not (piece.source.x.is_rational and piece.source.y.is_rational)]
    assert {isinstance(piece.support, OffsetCircleSupport) for piece in vertices} == {True, False}
    assert all(isinstance(piece.support, (OffsetCircleSupport, OffsetLineSupport)) for piece in vertices)
    lattice = [step * GRID_STEP for step in range(-GRID_HALF_STEPS, GRID_HALF_STEPS + 1)]
    contacts: list[DecidedContact] = []
    for piece in vertices:
        x, y = float(piece.source.x), float(piece.source.y)
        found = [contact for contact in _vertex_contacts(loop, [(x + dx, y + dy) for dx in lattice for dy in lattice]) if contact.piece.index == piece.index]
        assert found, piece.index
        contacts.extend(found)
    assert _assert_batch_equals_reference(loop, contacts) > 0
