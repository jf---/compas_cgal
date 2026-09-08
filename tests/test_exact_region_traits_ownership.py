"""Every ExactRegion2 must keep alive the traits object its arrangement reads.

CGAL's `Gps_on_surface_base_2` copy constructor allocates a fresh `Traits_2` for
the copy but builds the copy's arrangement as `Aos_2(*(ps.m_arr))`, and
`Arrangement_on_surface_2::assign` propagates a BORROWED traits pointer verbatim
(`m_geom_traits = arr.m_own_traits ? new Traits_adaptor_2 : arr.m_geom_traits`).
Every `Gps` arrangement borrows, so a copied set reads the traits of the ROOT set
it descends from, and the copy's own traits is never used.

`ExactRegion2` used to adopt its set through `make_shared<const ReachSet>(...)`.
`ReachSet` has no move constructor, so that adopted by copy, and the by-value
parameter -- the root the copy's arrangement points at -- was destroyed on the
return. Every region held a freed traits pointer from the moment it was built.

That is undefined behaviour, but not a crash today: `oriented_side` reaches point
location through `Arr_walk_along_line_point_location`, which only stores the
pointer and calls stateless functor accessors through it. The sibling defect on
`Stock2` DID crash, because `Stock2::contains` builds an
`Arr_trapezoid_ric_point_location`, whose constructor copy-constructs a
`Td_traits` out of the arrangement's traits. The difference between silence and a
SIGSEGV is one locator choice, so this suite asserts the OWNERSHIP INVARIANT
rather than a symptom: a region must reach its traits through a `ReachSet` it
keeps alive.

Answers are asserted alongside the invariant, so a change that satisfies the
audit by rebuilding geometry fails here too.
"""

import numpy as np
import pytest

from compas_cgal import _containment_2
from compas_cgal import _coverage_2

SQUARE = np.array(
    [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [10.0, 10.0, 0.0], [0.0, 10.0, 0.0]],
    dtype=np.float64,
)
ISLAND = np.array(
    [[4.0, 4.0, 0.0], [6.0, 4.0, 0.0], [6.0, 6.0, 0.0], [4.0, 6.0, 0.0]],
    dtype=np.float64,
)
TOOL_RADIUS = 0.5
# Inside the ring, inside the island void, outside, and hard against two edges.
PROBES = ((2.0, 2.0), (5.0, 5.0), (-1.0, 5.0), (11.0, 5.0), (0.0, 0.0))


def owns_traits(region: _coverage_2.ExactRegion2) -> bool:
    return region.arrangement_traits_are_owned_for_audit()


def probe(region: _coverage_2.ExactRegion2) -> tuple[bool, ...]:
    return tuple(region.contains(x, y) for x, y in PROBES)


@pytest.fixture
def domain() -> _coverage_2.ReachableDomain2:
    return _coverage_2.ReachableDomain2(SQUARE, [ISLAND], TOOL_RADIUS)


def test_design_region_from_polygon_owns_its_traits() -> None:
    design = _coverage_2.ExactRegion2.from_polygon(SQUARE, [ISLAND])
    assert owns_traits(design)
    assert probe(design) == (True, False, False, False, True)


def test_clone_owns_its_traits() -> None:
    design = _coverage_2.ExactRegion2.from_polygon(SQUARE, [ISLAND])
    clone = design.clone()
    del design
    assert owns_traits(clone)
    assert probe(clone) == (True, False, False, False, True)


def test_reachable_domain_regions_own_their_traits(domain: _coverage_2.ReachableDomain2) -> None:
    regions = {
        "design": domain.design_region(),
        "center": domain.center_domain(),
        "material": domain.reachable_material(),
        "residual": domain.unreachable_residual(),
    }
    del domain
    assert {name: owns_traits(region) for name, region in regions.items()} == {
        "design": True,
        "center": True,
        "material": True,
        "residual": True,
    }


def test_precleared_coverage_regions_own_their_traits(domain: _coverage_2.ReachableDomain2) -> None:
    coverage = _coverage_2.Coverage2(domain.reachable_material(), 1.0, 1.0, 0.4)
    assert owns_traits(coverage.residual())
    assert owns_traits(coverage.accumulated_sweeps())


def test_swept_coverage_regions_own_their_traits(domain: _coverage_2.ReachableDomain2) -> None:
    coverage = _coverage_2.Coverage2(domain.reachable_material(), 1.0, 1.0, 0.4)
    coverage.add_segment_sweep(1.0, 1.0, 8.0, 1.0, 0.4)
    coverage.add_full_circle_sweep(2.0, 8.0, 3.0, 8.0, 0.4)
    coverage.add_disk_sweep(8.0, 8.0, 0.4)
    assert owns_traits(coverage.residual())
    assert owns_traits(coverage.accumulated_sweeps())
    assert coverage.residual_component_count() == 2


def test_uncut_coverage_regions_own_their_traits() -> None:
    """`from_uncut` copies the target's set without operating on it, so nothing
    rebuilds the arrangement on traits of its own -- the case the sibling audit
    named."""
    design = _coverage_2.ExactRegion2.from_polygon(SQUARE, [ISLAND])
    coverage = _coverage_2.Coverage2.from_uncut(design)
    del design
    assert owns_traits(coverage.residual())
    assert owns_traits(coverage.accumulated_sweeps())
    assert probe(coverage.residual()) == (True, False, False, False, True)


def test_uncut_coverage_regions_own_their_traits_after_a_sweep() -> None:
    """The first sweep after `from_uncut` joins onto an EMPTY accumulated set,
    which is the branch where CGAL assigns the other arrangement wholesale."""
    design = _coverage_2.ExactRegion2.from_polygon(SQUARE, [ISLAND])
    coverage = _coverage_2.Coverage2.from_uncut(design)
    del design
    coverage.add_segment_sweep(1.0, 1.0, 8.0, 1.0, 0.4)
    assert owns_traits(coverage.accumulated_sweeps())
    assert owns_traits(coverage.residual())


def test_remaining_material_owns_its_traits(domain: _coverage_2.ReachableDomain2) -> None:
    material = domain.reachable_material()
    remaining = _coverage_2.remaining_material(
        material,
        np.array([[2.0, 8.0, 1.0]], dtype=np.float64),
        np.array([[1.0, 1.0, 8.0, 1.0]], dtype=np.float64),
        np.array([[8.0, 8.0]], dtype=np.float64),
        0.4,
    )
    del domain
    assert owns_traits(remaining)
    assert remaining.is_subset_of(material)


def test_remaining_material_without_motions_owns_its_traits(domain: _coverage_2.ReachableDomain2) -> None:
    """With no motions nothing is subtracted, so no boolean operation rebuilds
    the arrangement -- the un-operated copy the sibling audit flagged."""
    material = domain.reachable_material()
    remaining = _coverage_2.remaining_material(
        material,
        np.zeros((0, 3), dtype=np.float64),
        np.zeros((0, 4), dtype=np.float64),
        np.zeros((0, 2), dtype=np.float64),
        0.4,
    )
    del domain
    assert owns_traits(remaining)
    assert remaining.exactly_equals(material)


def test_native_boundary_design_region_owns_its_traits() -> None:
    boundary = _coverage_2.NativeBoundary2(
        [
            _coverage_2.NativeBoundaryCurve2.line((0.0, 0.0), (10.0, 0.0)),
            _coverage_2.NativeBoundaryCurve2.line((10.0, 0.0), (10.0, 10.0)),
            _coverage_2.NativeBoundaryCurve2.line((10.0, 10.0), (0.0, 10.0)),
            _coverage_2.NativeBoundaryCurve2.line((0.0, 10.0), (0.0, 0.0)),
        ]
    )
    design = boundary.design_region()
    del boundary
    assert owns_traits(design)
    assert design.contains(5.0, 5.0)


def test_orphaned_regions_still_answer_containment(domain: _coverage_2.ReachableDomain2) -> None:
    """Every producer at once, with each parent dropped before the query and the
    allocator churned in between, so a region that reads freed traits has to read
    reused bytes."""
    material = domain.reachable_material()
    coverage = _coverage_2.Coverage2(material, 1.0, 1.0, 0.4)
    coverage.add_segment_sweep(1.0, 1.0, 8.0, 1.0, 0.4)
    regions = [
        _coverage_2.ExactRegion2.from_polygon(SQUARE, [ISLAND]),
        domain.design_region(),
        domain.center_domain(),
        material,
        domain.unreachable_residual(),
        coverage.residual(),
        coverage.accumulated_sweeps(),
    ]
    del coverage
    del domain
    ballast = [bytearray(32) for _ in range(8192)]
    assert all(owns_traits(region) for region in regions)
    assert [probe(region) for region in regions] == [
        (True, False, False, False, True),
        (True, False, False, False, True),
        (True, False, False, False, False),
        (True, False, False, False, False),
        (False, False, False, False, True),
        (True, False, False, False, False),
        (False, False, False, False, False),
    ]
    del ballast


def test_containment_evaluation_survives_orphaned_regions(domain: _coverage_2.ReachableDomain2) -> None:
    """`containment_2.cpp` reads a region's arrangement through `oriented_side`."""
    design = domain.design_region()
    center = domain.center_domain()
    del domain
    authority = bytes(range(32))
    record = _containment_2.evaluate_exact_segment_containment(
        design, center, authority, 1.0, 1.0, 8.0, 1.0, 0.4
    )
    assert owns_traits(design)
    assert owns_traits(center)
    assert record.contained is True
    assert record.guide_anchor_in_center_domain is True
