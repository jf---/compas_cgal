"""`PostSetup.build` and its parts: every invariant raises its own named error."""

from fractions import Fraction

import pytest

from compas_cgal.post.errors import ArcToleranceBelowQuantumError
from compas_cgal.post.errors import ControllerResolutionError
from compas_cgal.post.errors import DegenerateEnvelopeError
from compas_cgal.post.errors import DepthExceedsStockError
from compas_cgal.post.errors import InexactQuantityError
from compas_cgal.post.errors import InvalidPlungeAllowanceError
from compas_cgal.post.errors import NonCentreCuttingPlungeError
from compas_cgal.post.errors import NonPositiveQuantityError
from compas_cgal.post.machine import CutterSpec
from compas_cgal.post.machine import CuttingParameters
from compas_cgal.post.machine import Dialect
from compas_cgal.post.machine import PostSetup
from compas_cgal.post.machine import StockBlank
from compas_cgal.post.machine import TravelEnvelope
from compas_cgal.post.testcut_defaults import build_testcut_defaults
from tests.post.support import TEST_ENVELOPE
from tests.post.support import rebuild


def test_defaults_carry_every_question(setup: PostSetup) -> None:
    assert setup.arc_tolerance == Fraction(2, 1000)
    assert setup.cutter.diameter == 6 and setup.cutter.flutes == 3 and setup.cutter.centre_cutting
    assert setup.scale == 3
    assert setup.axial_depth == Fraction(3, 2)
    assert (setup.cutting.spindle_rpm, setup.cutting.cut_feed, setup.cutting.plunge_feed) == (10000, 600, 100)
    assert setup.clearance_above_top == 12
    assert setup.envelope is TEST_ENVELOPE
    assert setup.xy_shift == (Fraction(18), Fraction(12))
    assert setup.cut_z == Fraction(-3, 2) and setup.clearance_z == 12
    assert setup.dialect is Dialect.LINUXCNC
    assert setup.controller_resolution == Fraction(1, 10000)
    assert setup.spindle_dwell == 2
    assert setup.stock.min_corner == (Fraction(-12), Fraction(-13))


def test_envelope_has_no_default() -> None:
    with pytest.raises(TypeError, match="envelope"):
        build_testcut_defaults(pocket_min_corner=(Fraction(0), Fraction(0)), accepted_plunges=1)  # type: ignore[call-arg]


def test_degenerate_envelope_raises() -> None:
    with pytest.raises(DegenerateEnvelopeError):
        TravelEnvelope.build(x=(Fraction(0), Fraction(0)), y=(Fraction(0), Fraction(1)), z=(Fraction(0), Fraction(1)))


def test_float_quantity_raises() -> None:
    with pytest.raises(InexactQuantityError):
        CutterSpec.build(diameter=6.0, flutes=3, centre_cutting=True)  # type: ignore[arg-type]


@pytest.mark.parametrize("field", ["spindle_rpm", "cut_feed", "plunge_feed"])
def test_non_positive_cutting_parameter_raises(field: str) -> None:
    values = {"spindle_rpm": 10000, "cut_feed": 600, "plunge_feed": 100, field: 0}
    with pytest.raises(NonPositiveQuantityError):
        CuttingParameters.build(**values)


def test_non_positive_scale_raises(setup: PostSetup) -> None:
    with pytest.raises(NonPositiveQuantityError):
        rebuild(setup, scale=Fraction(0))


def test_scope_quantum_cannot_hold_the_controller_tolerance(setup: PostSetup) -> None:
    # 0.001 mm quantum: worst-case emitted radius gap 2*sqrt(2)*q = 0.00283 mm > 0.002 mm.
    with pytest.raises(ArcToleranceBelowQuantumError):
        rebuild(setup, decimals=3)


def test_plunge_needs_centre_cutting_tool(setup: PostSetup) -> None:
    with pytest.raises(NonCentreCuttingPlungeError):
        rebuild(setup, cutter=CutterSpec.build(diameter=Fraction(6), flutes=3, centre_cutting=False))


def test_depth_beyond_stock_raises(setup: PostSetup) -> None:
    with pytest.raises(DepthExceedsStockError):
        rebuild(setup, stock=StockBlank.build(material="wax", min_corner=(Fraction(0), Fraction(0)), size_x=Fraction(60), size_y=Fraction(50), size_z=Fraction(1)))


def test_negative_plunge_allowance_raises(setup: PostSetup) -> None:
    with pytest.raises(InvalidPlungeAllowanceError):
        rebuild(setup, accepted_plunges=-1)


def test_controller_that_would_round_emitted_digits_raises(setup: PostSetup) -> None:
    # 5 digits (1e-5 mm) on a 1e-4 mm controller: the controller would round the fifth digit.
    with pytest.raises(ControllerResolutionError):
        rebuild(setup, decimals=5)
    rebuild(setup, decimals=5, controller_resolution=Fraction(1, 10**5))


def test_non_positive_dwell_raises(setup: PostSetup) -> None:
    with pytest.raises(NonPositiveQuantityError):
        rebuild(setup, spindle_dwell=Fraction(0))
