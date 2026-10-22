"""Ω* output is a function of the input alone: no hash-seed or cache-state dependence.

Two defects this pins:

- owner sets of ∂Ω* pieces, vertex disks and unreachable stretches were plain
  ``frozenset``s of :class:`BoundaryFeature`, whose ``repr`` and iteration order
  follow the per-process string hash seed, so cross-process digests of
  identical Ω* paths differed (now :class:`OwnerSet`, canonical order);
- :class:`OmegaSite` compared by identity, so two identical sites were unequal.
"""

from __future__ import annotations

import os
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

import pytest

from benchmarks.held_centre_domain_supports import BoundaryFeature
from benchmarks.held_centre_domain_supports import BoundaryFeatureKind
from benchmarks.held_centre_domain_supports import OwnerSet
from benchmarks.held_centre_domain_supports import canonical_owners
from benchmarks.held_centre_domain_supports import owner_order_key
from benchmarks.held_exact_centre_domain import CentreDomain
from benchmarks.held_omega_sites import OmegaLoop
from benchmarks.held_omega_sites import OmegaSite
from tests.benchmarks.test_held_exact_centre_domain import concave_arc_pocket
from tests.benchmarks.test_held_exact_centre_domain import l_shape_pocket
from tests.benchmarks.test_held_omega_path import GAP
from tests.benchmarks.test_held_omega_path import _arc_contacts

REPO = Path(__file__).resolve().parents[2]
# Distinct string hash seeds: with k two-owner sets each seed permutes every set
# independently with probability 1/2, so four seeds expose a seed-dependent repr.
HASH_SEEDS = ("1", "2", "3", "4")
# Child process: SHA-256 of the repr of each fixture's Ω* (pieces, vertex disks,
# unreachable stretches) and of the contained Ω* circles on the concave fixture.
DIGEST_SCRIPT = """
import hashlib
from fractions import Fraction
from benchmarks.held_declared_clearance import DeclaredClearance
from benchmarks.held_exact_centre_domain import CentreDomain
from tests.benchmarks.test_held_exact_centre_domain import concave_arc_pocket, l_shape_pocket
from tests.benchmarks.test_held_omega_path import _arc_contacts, _omega_circles
parts = []
for boundary in (l_shape_pocket(), concave_arc_pocket()):
    domain = CentreDomain.build(boundary, Fraction(1))
    parts.append(repr((domain.components, domain.vertex_disks, domain.unreachable)))
boundary = concave_arc_pocket()
domain = CentreDomain.build(boundary, Fraction(1))
parts.append(repr(_omega_circles(domain, DeclaredClearance.build(boundary, Fraction(1)), _arc_contacts())))
print(" ".join(hashlib.sha256(part.encode()).hexdigest() for part in parts))
"""


def _digests(seed: str) -> str:
    env = {**os.environ, "PYTHONHASHSEED": seed}
    result = subprocess.run([sys.executable, "-c", DIGEST_SCRIPT], cwd=REPO, env=env, capture_output=True, text=True, check=True)
    return result.stdout.strip()


def test_omega_repr_digests_do_not_depend_on_the_hash_seed() -> None:
    digests = {seed: _digests(seed) for seed in HASH_SEEDS}
    assert len(set(digests.values())) == 1, digests


def test_owner_sets_iterate_and_print_in_canonical_order() -> None:
    for boundary in (l_shape_pocket(), concave_arc_pocket()):
        domain = CentreDomain.build(boundary, Fraction(1))
        pieces = [piece for component in domain.components for cycle in (component.outer, *component.holes) for piece in cycle]
        owner_sets = [*(piece.owners for piece in pieces), *(disk.owners for disk in domain.vertex_disks), *(stretch.owners for stretch in domain.unreachable)]
        assert owner_sets
        for owners in owner_sets:
            assert isinstance(owners, OwnerSet)
            ordered = sorted(frozenset(owners), key=owner_order_key)
            assert list(owners) == ordered
            assert repr(owners) == f"OwnerSet({tuple(ordered)!r})"


def test_owner_set_keeps_set_semantics() -> None:
    line, arc, junction = BoundaryFeature(BoundaryFeatureKind.LINE, 5), BoundaryFeature(BoundaryFeatureKind.ARC, 1), BoundaryFeature(BoundaryFeatureKind.JUNCTION, 0)
    owners = canonical_owners([junction, arc, line, arc])
    assert list(owners) == [line, arc, junction]
    assert owners == frozenset({line, arc, junction})
    assert hash(owners) == hash(frozenset({line, arc, junction}))
    assert list(canonical_owners(owners | {BoundaryFeature(BoundaryFeatureKind.LINE, 2)})) == [BoundaryFeature(BoundaryFeatureKind.LINE, 2), line, arc, junction]


def _two_identical_sites() -> tuple[OmegaSite, OmegaSite]:
    domain = CentreDomain.build(concave_arc_pocket(), Fraction(1))
    loop = OmegaLoop.build(domain)
    decision = domain.decide_contacts(_arc_contacts()[:1], GAP).contacts[0]
    assert decision is not None
    return OmegaSite.decided(loop, decision), OmegaSite.decided(loop, decision)


def test_identical_sites_are_equal_whatever_their_point_cache() -> None:
    first, second = _two_identical_sites()
    assert first == second
    first.point()
    assert first == second
    assert second == first
    second.point()
    assert first == second


def test_sites_are_unhashable_like_their_exact_parts() -> None:
    first, _ = _two_identical_sites()
    with pytest.raises(TypeError):
        hash(first)


def test_sites_at_different_positions_are_unequal() -> None:
    domain = CentreDomain.build(concave_arc_pocket(), Fraction(1))
    loop = OmegaLoop.build(domain)
    first, second = (OmegaSite.decided(loop, decision) for decision in domain.decide_contacts(_arc_contacts()[:2], GAP).contacts if decision is not None)
    assert first != second
