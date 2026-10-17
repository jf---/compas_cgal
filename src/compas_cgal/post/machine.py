"""The machine and job setup of one post run: every scope question as an explicit, validated field.

`PostSetup.build` owns the error model. The test-cut answers live in exactly one
place, `compas_cgal.post.testcut_defaults`, and every value there is a
placeholder pending the user's answer to the numbered scope question
(`testcut-scope.md`, "Questions for the user"). The travel envelope (Q10) has
no default.
"""

from dataclasses import dataclass
from enum import Enum
from fractions import Fraction
from typing import Self

from compas_cgal.helix_units import HelixPitchPerRadius
from compas_cgal.post.errors import ArcToleranceBelowQuantumError
from compas_cgal.post.errors import ControllerResolutionError
from compas_cgal.post.errors import DegenerateEnvelopeError
from compas_cgal.post.errors import DepthExceedsStockError
from compas_cgal.post.errors import InexactQuantityError
from compas_cgal.post.errors import InvalidPlungeAllowanceError
from compas_cgal.post.errors import NonCentreCuttingPlungeError
from compas_cgal.post.errors import NonPositiveQuantityError
from compas_cgal.post.errors import OffGridValueError
from compas_cgal.post.errors import RampStartAboveClearanceError
from compas_cgal.post.grid import quantum
from compas_cgal.post.quantities import DecimalPlaces
from compas_cgal.post.quantities import MillimetrePerMinute
from compas_cgal.post.quantities import ProgramMillimetre
from compas_cgal.post.quantities import RevolutionsPerMinute
from compas_cgal.post.quantities import ScaleFactor
from compas_cgal.post.quantities import Seconds

# Worst-case gap between an emitted arc's start and end radii, in output quanta,
# squared. Each emitted point is off its exact position by at most q/sqrt(2)
# (half a step per axis). The centre is written as I/J from the emitted start,
# so the centre the controller sees is the quantised centre; r_start moves by
# |d_start| + |d_centre| and r_end by |d_end| + |d_centre|, so
# |r_start - r_end| <= 4 * q / sqrt(2) = 2*sqrt(2)*q, i.e. 8 q^2 squared.
ARC_RADIUS_GAP_QUANTA_SQUARED = 8

# The stock top in the work frame, mm: Z0 by the Q9 convention (the cut plane is -axial_depth below it).
STOCK_TOP_Z = ProgramMillimetre(Fraction(0))


def _exact(value: object, name: str) -> Fraction:
    if not isinstance(value, Fraction):
        raise InexactQuantityError(f"{name} must be a Fraction, got {type(value).__name__}.")
    return value


def _positive(value: object, name: str) -> Fraction:
    exact = _exact(value, name)
    if exact <= 0:
        raise NonPositiveQuantityError(f"{name} must be positive, got {exact}.")
    return exact


def _positive_count(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise InexactQuantityError(f"{name} must be an int, got {type(value).__name__}.")
    if value <= 0:
        raise NonPositiveQuantityError(f"{name} must be positive, got {value}.")
    return value


class Dialect(Enum):
    """Q1: controller dialect. LinuxCNC is the only one the writer maps.

    Not Fanuc: Fanuc has no G91.1 (its I/J are always incremental), and its
    input resolution may round the fourth mm decimal, which would void the
    arc-radius certificate. Any other member makes the writer raise.
    """

    LINUXCNC = "linuxcnc"


class ArcCentreMode(Enum):
    """Q2: arc plane and centre addressing. Only G17 with incremental I/J (G91.1)."""

    XY_INCREMENTAL_IJ = "G17 G91.1"


class WorkOffset(Enum):
    """Q9: the work coordinate system the program runs in."""

    G54 = "G54"


@dataclass(frozen=True)
class CutterSpec:
    """Q4: the cutter.

    Attributes:
        diameter: Cutting diameter, mm.
        flutes: Flute count.
        centre_cutting: Whether the end face cuts (required for any plunge into material).
    """

    diameter: ProgramMillimetre
    flutes: int
    centre_cutting: bool

    @classmethod
    def build(cls, diameter: Fraction, flutes: int, centre_cutting: bool) -> Self:
        """Validate and freeze.

        Raises:
            InexactQuantityError: *diameter* is not a Fraction or *flutes* not an int.
            NonPositiveQuantityError: Either is not positive.
        """
        return cls(ProgramMillimetre(_positive(diameter, "cutter diameter")), _positive_count(flutes, "flute count"), bool(centre_cutting))


@dataclass(frozen=True)
class StockBlank:
    """Q6: the stock blank.

    Attributes:
        material: Material name (reporting only).
        min_corner: Where the blank's min XY corner sits in the work frame, mm.
        size_x: Blank length along X, mm.
        size_y: Blank length along Y, mm.
        size_z: Blank thickness, mm.
    """

    material: str
    min_corner: tuple[ProgramMillimetre, ProgramMillimetre]
    size_x: ProgramMillimetre
    size_y: ProgramMillimetre
    size_z: ProgramMillimetre

    @classmethod
    def build(cls, material: str, min_corner: tuple[Fraction, Fraction], size_x: Fraction, size_y: Fraction, size_z: Fraction) -> Self:
        """Validate and freeze.

        Raises:
            InexactQuantityError: A corner coordinate or size is not a Fraction.
            NonPositiveQuantityError: A size is not positive.
        """
        return cls(
            material,
            (ProgramMillimetre(_exact(min_corner[0], "stock min X")), ProgramMillimetre(_exact(min_corner[1], "stock min Y"))),
            ProgramMillimetre(_positive(size_x, "stock X")),
            ProgramMillimetre(_positive(size_y, "stock Y")),
            ProgramMillimetre(_positive(size_z, "stock Z")),
        )

    @property
    def x_range(self) -> tuple[ProgramMillimetre, ProgramMillimetre]:
        """(min, max) X of the blank in the work frame."""
        return self.min_corner[0], ProgramMillimetre(self.min_corner[0] + self.size_x)

    @property
    def y_range(self) -> tuple[ProgramMillimetre, ProgramMillimetre]:
        """(min, max) Y of the blank in the work frame."""
        return self.min_corner[1], ProgramMillimetre(self.min_corner[1] + self.size_y)


@dataclass(frozen=True)
class CuttingParameters:
    """Q8 (speeds and feeds).

    Attributes:
        spindle_rpm: Spindle speed.
        cut_feed: Feed of every XY motion on the cut plane.
        plunge_feed: Feed of every plunge.
    """

    spindle_rpm: RevolutionsPerMinute
    cut_feed: MillimetrePerMinute
    plunge_feed: MillimetrePerMinute

    @classmethod
    def build(cls, spindle_rpm: int, cut_feed: int, plunge_feed: int) -> Self:
        """Validate and freeze.

        Raises:
            InexactQuantityError: A value is not an int.
            NonPositiveQuantityError: A value is not positive.
        """
        return cls(
            RevolutionsPerMinute(_positive_count(spindle_rpm, "spindle rpm")),
            MillimetrePerMinute(_positive_count(cut_feed, "cut feed")),
            MillimetrePerMinute(_positive_count(plunge_feed, "plunge feed")),
        )


@dataclass(frozen=True)
class TravelEnvelope:
    """Q10: the machine's travel box in the work frame, mm, closed on every side."""

    x: tuple[ProgramMillimetre, ProgramMillimetre]
    y: tuple[ProgramMillimetre, ProgramMillimetre]
    z: tuple[ProgramMillimetre, ProgramMillimetre]

    @classmethod
    def build(cls, x: tuple[Fraction, Fraction], y: tuple[Fraction, Fraction], z: tuple[Fraction, Fraction]) -> Self:
        """Validate and freeze.

        Args:
            x: (min, max) along X.
            y: (min, max) along Y.
            z: (min, max) along Z.

        Raises:
            InexactQuantityError: A bound is not a Fraction.
            DegenerateEnvelopeError: An axis has min >= max.
        """
        axes: list[tuple[ProgramMillimetre, ProgramMillimetre]] = []
        for name, (low, high) in (("x", x), ("y", y), ("z", z)):
            low_exact, high_exact = _exact(low, f"envelope {name} min"), _exact(high, f"envelope {name} max")
            if low_exact >= high_exact:
                raise DegenerateEnvelopeError(f"envelope {name}: min {low_exact} must be below max {high_exact}.")
            axes.append((ProgramMillimetre(low_exact), ProgramMillimetre(high_exact)))
        return cls(axes[0], axes[1], axes[2])


@dataclass(frozen=True)
class PostSetup:
    """Everything one post run needs beyond the operations; one field per scope question.

    Construct with `build` (validates) or `compas_cgal.post.testcut_defaults.build_testcut_defaults`.

    Attributes:
        dialect: Q1.
        controller_resolution: Q1, the smallest coordinate step the controller reads without rounding, mm.
        arc_mode: Q2.
        arc_tolerance: Q3, the controller's arc radius tolerance, mm.
        cutter: Q4.
        scale: Q5, source units to mm.
        stock: Q6.
        axial_depth: Q7, depth of the single Z level below the stock top (Z0), mm.
        cutting: Q8, speeds and feeds.
        spindle_dwell: Q8, the wait after spindle-on before any motion, s.
        clearance_above_top: Q8, the rapid plane above the stock top, mm.
        work_offset: Q9, work coordinate system.
        xy_shift: Q9, where the source origin lands in the work frame, mm (applied after scaling).
        envelope: Q10, required.
        decimals: Output quantum digits (scope, not a user question).
        accepted_plunges: Plunges the user accepts (the ledger has no acceptance field; audit section 4).
        ramp_feed: Q8, feed of every block of a HELIX entry (ramp turns and the flat closing turn).
        ramp_cap: The machine ramp-angle cap as its exact surrogate ``K = 2*pi*tan(alpha)``
            (`compas_cgal.helix_units`): an emitted turn of descent ``D`` on radius ``r`` is
            within the cap iff ``D <= exact(K) * r``. Certified per full turn (exact: each
            turn sweeps exactly 2 pi); a single half-turn block can be steeper by at most
            ``q / (pi r)`` in slope (its descent is within ``q`` of half the turn's), about
            1e-5 at r = 3 mm against tan(2 deg) ~ 0.035.
        ramp_start_above_top: Q8, where a HELIX entry's ramp starts above the stock top, mm: the tool
            feeds (G1, plunge feed, in air) from clearance down to it, then ramps from it to the cut plane.
    """

    dialect: Dialect
    controller_resolution: ProgramMillimetre
    arc_mode: ArcCentreMode
    arc_tolerance: ProgramMillimetre
    cutter: CutterSpec
    scale: ScaleFactor
    stock: StockBlank
    axial_depth: ProgramMillimetre
    cutting: CuttingParameters
    spindle_dwell: Seconds
    clearance_above_top: ProgramMillimetre
    work_offset: WorkOffset
    xy_shift: tuple[ProgramMillimetre, ProgramMillimetre]
    envelope: TravelEnvelope
    decimals: DecimalPlaces
    accepted_plunges: int
    ramp_feed: MillimetrePerMinute
    ramp_cap: HelixPitchPerRadius
    ramp_start_above_top: ProgramMillimetre

    @property
    def ramp_start_z(self) -> ProgramMillimetre:
        """Where a HELIX entry's ramp starts in the work frame, strictly between the stock top and clearance."""
        return ProgramMillimetre(STOCK_TOP_Z + self.ramp_start_above_top)

    @property
    def cut_z(self) -> ProgramMillimetre:
        """The cut plane in the work frame (Z0 = stock top)."""
        return ProgramMillimetre(STOCK_TOP_Z - self.axial_depth)

    @property
    def clearance_z(self) -> ProgramMillimetre:
        """The rapid plane in the work frame."""
        return self.clearance_above_top

    @classmethod
    def build(
        cls,
        *,
        dialect: Dialect,
        controller_resolution: Fraction,
        arc_mode: ArcCentreMode,
        arc_tolerance: Fraction,
        cutter: CutterSpec,
        scale: Fraction,
        stock: StockBlank,
        axial_depth: Fraction,
        cutting: CuttingParameters,
        spindle_dwell: Fraction,
        clearance_above_top: Fraction,
        work_offset: WorkOffset,
        xy_shift: tuple[Fraction, Fraction],
        envelope: TravelEnvelope,
        decimals: int,
        accepted_plunges: int,
        ramp_feed: int,
        ramp_cap: HelixPitchPerRadius,
        ramp_start_above_top: Fraction,
    ) -> Self:
        """Validate every field and the invariants between them.

        Raises:
            InexactQuantityError: A length is not a Fraction, a count or feed not an int, or the ramp cap not a `HelixPitchPerRadius`.
            NonPositiveQuantityError: A resolution, tolerance, scale, depth, dwell, clearance, ramp feed, ramp start or digit count is not positive.
            RampStartAboveClearanceError: ``ramp_start_above_top >= clearance_above_top``.
            OffGridValueError: ``spindle_dwell`` or ``ramp_start_above_top`` is not on the output grid.
            ControllerResolutionError: The output quantum is not a whole multiple of ``controller_resolution``.
            ArcToleranceBelowQuantumError: ``arc_tolerance**2 <= 8 * quantum**2``.
            DepthExceedsStockError: ``axial_depth >= stock.size_z``.
            InvalidPlungeAllowanceError: ``accepted_plunges < 0``.
            NonCentreCuttingPlungeError: Plunges accepted with a non-centre-cutting cutter.
        """
        tolerance = _positive(arc_tolerance, "arc tolerance")
        digits = DecimalPlaces(_positive_count(decimals, "decimal places"))
        if tolerance * tolerance <= ARC_RADIUS_GAP_QUANTA_SQUARED * quantum(digits) ** 2:
            raise ArcToleranceBelowQuantumError(f"quantum {quantum(digits)} mm can open an arc's radii by 2*sqrt(2)*q, beyond the {tolerance} mm tolerance.")
        resolution = _positive(controller_resolution, "controller resolution")
        if (quantum(digits) / resolution).denominator != 1:
            raise ControllerResolutionError(f"quantum {quantum(digits)} mm is not a multiple of the controller's {resolution} mm: it would round emitted digits.")
        dwell = _positive(spindle_dwell, "spindle dwell")
        if (dwell / quantum(digits)).denominator != 1:
            raise OffGridValueError(f"dwell {dwell} s is not on the 10**-{digits} grid.")
        depth = _positive(axial_depth, "axial depth")
        if depth >= stock.size_z:
            raise DepthExceedsStockError(f"depth {depth} mm reaches through the {stock.size_z} mm blank.")
        if isinstance(accepted_plunges, bool) or not isinstance(accepted_plunges, int):
            raise InexactQuantityError("accepted plunges must be an int.")
        if accepted_plunges < 0:
            raise InvalidPlungeAllowanceError(f"accepted plunges must be >= 0, got {accepted_plunges}.")
        if not isinstance(ramp_cap, HelixPitchPerRadius):
            raise InexactQuantityError(f"ramp cap must be a HelixPitchPerRadius, got {type(ramp_cap).__name__}.")
        clearance = _positive(clearance_above_top, "clearance above top")
        ramp_start = _positive(ramp_start_above_top, "ramp start above top")
        if ramp_start >= clearance:
            raise RampStartAboveClearanceError(f"ramp start {ramp_start} mm above the top is not below the {clearance} mm clearance.")
        if (ramp_start / quantum(digits)).denominator != 1:
            raise OffGridValueError(f"ramp start {ramp_start} mm is not on the 10**-{digits} grid.")
        if accepted_plunges > 0 and not cutter.centre_cutting:
            raise NonCentreCuttingPlungeError(f"{accepted_plunges} plunge(s) accepted with a non-centre-cutting cutter.")
        return cls(
            dialect=dialect,
            controller_resolution=ProgramMillimetre(resolution),
            arc_mode=arc_mode,
            arc_tolerance=ProgramMillimetre(tolerance),
            cutter=cutter,
            scale=ScaleFactor(_positive(scale, "scale")),
            stock=stock,
            axial_depth=ProgramMillimetre(depth),
            cutting=cutting,
            spindle_dwell=Seconds(dwell),
            clearance_above_top=ProgramMillimetre(clearance),
            work_offset=work_offset,
            xy_shift=(ProgramMillimetre(_exact(xy_shift[0], "x shift")), ProgramMillimetre(_exact(xy_shift[1], "y shift"))),
            envelope=envelope,
            decimals=digits,
            accepted_plunges=accepted_plunges,
            ramp_feed=MillimetrePerMinute(_positive_count(ramp_feed, "ramp feed")),
            ramp_cap=ramp_cap,
            ramp_start_above_top=ProgramMillimetre(ramp_start),
        )
