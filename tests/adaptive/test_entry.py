import hashlib
from dataclasses import replace

import pytest
from compas.geometry import Polygon

from compas_cgal.adaptive.canonical import CanonicalRingV1
from compas_cgal.adaptive.entry import BoreProcessIdentity
from compas_cgal.adaptive.entry import PreclearedEntry
from compas_cgal.adaptive.entry import QualifiedBore
from compas_cgal.adaptive.errors import GougeContainmentError
from compas_cgal.adaptive.errors import InvalidContainmentCertificateError
from compas_cgal.adaptive.errors import InvalidDepletionWitnessError
from compas_cgal.adaptive.errors import InvalidPreclearedEntryError
from compas_cgal.adaptive.errors import InvalidStockAreaError
from compas_cgal.adaptive.motion import ExactCircleMotion
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.operation import ApproachOperation
from compas_cgal.adaptive.operation import PlungeOperation
from compas_cgal.adaptive.reachable_domain import ReachableDomain
from compas_cgal.adaptive.stock_area import Stock2Area
from compas_cgal.adaptive.units import ClearanceZ
from compas_cgal.adaptive.units import CutPlane
from compas_cgal.adaptive.units import CutZ
from compas_cgal.adaptive.units import EntryRadius
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.adaptive.units import Vector2
from compas_cgal.adaptive.units import WorldXY
from compas_cgal.stock import Stock

OUTER = (
    (0.0, 0.0),
    (12.0, 0.0),
    (12.0, 12.0),
    (0.0, 12.0),
)
ISLAND = (
    (5.0, 5.0),
    (7.0, 5.0),
    (7.0, 7.0),
    (5.0, 7.0),
)


def _ring(
    points: tuple[tuple[float, float], ...],
    *,
    outer: bool,
) -> CanonicalRingV1:
    framed = tuple(Point2[WorldXY].build(x, y) for x, y in points)
    return CanonicalRingV1.build_outer(framed) if outer else CanonicalRingV1.build_hole(framed)


def _domain(
    *,
    holes: tuple[tuple[tuple[float, float], ...], ...] = (),
    tool_radius: float = 0.5,
) -> ReachableDomain:
    return ReachableDomain.build(
        design_boundary=_ring(OUTER, outer=True),
        holes=tuple(_ring(hole, outer=False) for hole in holes),
        tool_radius=ToolRadius.build(tool_radius),
    )


def _cut_plane(*, cut_z: float = -2.0, clearance_z: float = 5.0) -> CutPlane:
    return CutPlane.build(
        CutZ.build(cut_z),
        ClearanceZ.build(clearance_z),
    )


def _qualification(
    cut_plane: CutPlane,
    *,
    evidence: bytes = b"qualified-bore-metrology-v1",
) -> QualifiedBore:
    return QualifiedBore.build(
        cut_plane=cut_plane,
        process_identity=BoreProcessIdentity(b"predrill-cycle-revision-7"),
        evidence_bytes=evidence,
    )


def _entry(
    *,
    domain: ReachableDomain | None = None,
    center: Point2[WorldXY] | None = None,
    radius: float = 2.0,
    tool_radius: float = 0.5,
    cut_plane: CutPlane | None = None,
    evidence: bytes = b"qualified-bore-metrology-v1",
) -> PreclearedEntry:
    owned_domain = domain if domain is not None else _domain(tool_radius=tool_radius)
    owned_cut_plane = cut_plane if cut_plane is not None else _cut_plane()
    return PreclearedEntry.build(
        reachable_domain=owned_domain,
        center=center if center is not None else Point2[WorldXY].build(3.0, 3.0),
        radius=EntryRadius.build(radius),
        tool_radius=ToolRadius.build(tool_radius),
        cut_plane=owned_cut_plane,
        qualified_bore=_qualification(
            owned_cut_plane,
            evidence=evidence,
        ),
    )


def _stock(*, holes: tuple[tuple[tuple[float, float], ...], ...] = ()) -> Stock:
    return Stock(
        Polygon([[x, y, 0.0] for x, y in OUTER]),
        [Polygon([[x, y, 0.0] for x, y in hole]) for hole in holes],
    )


def test_precleared_entry_binds_exact_disk_bore_and_canonical_vertical_path() -> None:
    entry = _entry()

    assert entry.containment.center_in_center_domain
    assert entry.containment.disk_contained
    assert entry.qualified_bore.cut_plane is entry.cut_plane
    assert hashlib.sha256(entry.qualified_bore.evidence_bytes).digest() == entry.qualified_bore.evidence_digest
    assert type(entry.approach) is ApproachOperation
    assert entry.approach.position == entry.center
    assert entry.approach.clearance_z == entry.cut_plane.clearance_z
    assert type(entry.plunge) is PlungeOperation
    assert entry.plunge.position == entry.center
    assert entry.plunge.clearance_z == entry.cut_plane.clearance_z
    assert entry.plunge.cut_z == entry.cut_plane.cut_z
    assert len(entry.digest) == 32

    changed_process = _entry(evidence=b"qualified-bore-metrology-v2")
    changed_radius = _entry(radius=2.25)
    assert changed_process.digest != entry.digest
    assert changed_radius.digest != entry.digest


def test_entry_center_disk_island_and_bore_interval_are_independent_contracts() -> None:
    domain = _domain(tool_radius=1.0)
    cut_plane = _cut_plane()

    with pytest.raises(InvalidPreclearedEntryError, match="launch"):
        _entry(
            domain=domain,
            radius=1.0,
            tool_radius=1.0,
            cut_plane=cut_plane,
        )

    with pytest.raises(GougeContainmentError, match="entry disk"):
        _entry(
            domain=domain,
            center=Point2[WorldXY].build(1.0, 1.0),
            radius=1.25,
            tool_radius=1.0,
            cut_plane=cut_plane,
        )

    island_domain = _domain(holes=(ISLAND,), tool_radius=0.5)
    with pytest.raises(GougeContainmentError, match="entry disk"):
        _entry(
            domain=island_domain,
            center=Point2[WorldXY].build(6.0, 6.0),
            radius=1.0,
        )

    other_cut_plane = _cut_plane(cut_z=-3.0)
    with pytest.raises(InvalidPreclearedEntryError, match="qualified interval"):
        PreclearedEntry.build(
            reachable_domain=_domain(),
            center=Point2[WorldXY].build(3.0, 3.0),
            radius=EntryRadius.build(2.0),
            tool_radius=ToolRadius.build(0.5),
            cut_plane=cut_plane,
            qualified_bore=_qualification(other_cut_plane),
        )


def test_entry_depletion_is_first_exactly_once_and_path_moves_no_stock() -> None:
    entry = _entry()
    area = Stock2Area.build(_stock())
    initial = area.raw

    assert entry.approach
    assert entry.plunge
    assert area.raw.exactly_equals(initial)

    witness = area.deplete(entry)

    assert witness.motion is entry
    assert witness.policy is None
    assert witness.center_parameters == ()
    assert witness.parent_lineage == ()
    assert area.lineage == (witness,)
    assert not area.raw.contains(entry.center.x, entry.center.y)
    assert not area.raw.exactly_equals(initial)

    with pytest.raises(InvalidStockAreaError, match="exactly once"):
        area.deplete(entry)

    arbitrary_lateral = ExactSegmentMotion.build(
        entry.center,
        Point2[WorldXY].build(3.25, 3.0),
    )
    with pytest.raises(InvalidStockAreaError, match="three depletion inputs"):
        Stock2Area.build(_stock()).deplete(arbitrary_lateral)

    with pytest.raises(InvalidStockAreaError, match="supported exact depletion"):
        Stock2Area.build(_stock()).deplete(entry.approach)
    with pytest.raises(InvalidStockAreaError, match="supported exact depletion"):
        Stock2Area.build(_stock()).deplete(entry.plunge)


def test_entry_depletion_witness_rejects_every_motion_witness_field() -> None:
    entry = _entry()
    witness = Stock2Area.build(_stock()).deplete(entry)

    with pytest.raises(InvalidDepletionWitnessError, match="policy"):
        replace(witness, policy=object())  # type: ignore[arg-type]
    with pytest.raises(InvalidDepletionWitnessError, match="no sampled"):
        replace(witness, center_parameters=(object(),))  # type: ignore[arg-type]
    with pytest.raises(InvalidDepletionWitnessError, match="tool radius"):
        replace(witness, tool_radius=ToolRadius.build(0.25))
    with pytest.raises(InvalidDepletionWitnessError, match="strategy"):
        replace(witness, native_strategy_version=b"wrong-entry-strategy-v0")
    with pytest.raises(InvalidDepletionWitnessError, match="first"):
        replace(witness, parent_lineage=(b"\x01" * 32,))  # type: ignore[arg-type]


def test_first_full_circle_must_fit_declared_entry_void() -> None:
    entry = _entry()
    area = Stock2Area.build(_stock())
    area.deplete(entry)
    first_circle = ExactCircleMotion.build(
        entry.center,
        Vector2[WorldXY].build(1.0, 0.0),
        False,
    )

    certificate = entry.certify_first_circle(first_circle)

    assert certificate.motion is first_circle
    assert certificate.entry_digest == entry.digest
    assert not area.raw.contains(
        entry.center.x + first_circle.phase_vector.x,
        entry.center.y + first_circle.phase_vector.y,
    )

    outside = ExactCircleMotion.build(
        entry.center,
        Vector2[WorldXY].build(1.75, 0.0),
        False,
    )
    with pytest.raises(GougeContainmentError, match="entry disk"):
        entry.certify_first_circle(outside)


def test_entry_certificates_cannot_be_relabelled_to_another_authority() -> None:
    entry = _entry()
    other_domain = _domain(holes=(ISLAND,))

    with pytest.raises(
        InvalidContainmentCertificateError,
        match="native record",
    ):
        replace(
            entry.containment,
            domain_digest=other_domain.certificate.digest,
        )

    motion = ExactCircleMotion.build(
        entry.center,
        Vector2[WorldXY].build(1.0, 0.0),
        False,
    )
    launch = entry.certify_first_circle(motion)

    with pytest.raises(
        InvalidContainmentCertificateError,
        match="native record",
    ):
        replace(launch, entry_digest=b"\x00" * 32)
