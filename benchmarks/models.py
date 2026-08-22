"""Calibration coefficients: the line between what we computed and what we assumed.

Everything in `benchmarks.quality` that reads only the toolpath and the exact
stock is OURS -- reproducible from this repository, defensible on its own. The
moment a number depends on a workpiece material or a machine's dynamics it stops
being a property of the path and becomes a property of a CALIBRATION, and it is
only as good as that calibration. These two frozen models are where that
dependency is made visible: a metric needing one takes it as an argument, and
raises a NAMED error when it is absent rather than substituting a guessed
coefficient.

The defaults below are documented, cited to a concept, and DEFAULTS -- not
measurements of any real machine or alloy in this project. A figure or a report
built on them must say so. `MaterialModel.build` and `MachineModel.build`
validate every coefficient, so an out-of-range value fails at construction rather
than propagating into a plausible-looking number.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from benchmarks.errors import InvalidMachineModelError
from benchmarks.errors import InvalidMaterialModelError

# Feed per tooth for the default material model, in millimetres. A mid-range
# value for a carbide endmill in aluminium; it is the SCALE that converts the
# geometry-derived chip-thickness RATIO into a length, and nothing else in this
# module depends on it.
DEFAULT_FEED_PER_TOOTH_MM = 0.05

# Cutting-edge radius of a sharp carbide endmill, in millimetres. Ground carbide
# edges are commonly quoted at 3-10 um; 5 um is the middle of that band.
DEFAULT_EDGE_RADIUS_MM = 0.005

# Minimum undeformed chip thickness as a fraction of the cutting-edge radius.
# Below this the edge ploughs and rubs instead of shearing a chip, which wears
# the tool faster than a heavier cut does. The literature places the ratio
# between roughly 0.05 and 0.3 depending on material and rake; 0.2 is mid-band
# and gives h_min ~ 1 um at the default edge radius.
DEFAULT_MIN_CHIP_FRACTION_OF_EDGE_RADIUS = 0.2

# Taylor tool-life exponent n in V * T^n = C. Carbide tooling is usually quoted
# at 0.2-0.4; 0.25 is the conventional mid-value for carbide in steel.
DEFAULT_TAYLOR_EXPONENT = 0.25

# Taylor constant C, in metres per minute, for V in m/min and T in minutes. A
# round mid-range figure for carbide; it fixes only the SCALE of a life estimate
# and never the shape of its dependence on speed.
DEFAULT_TAYLOR_CONSTANT_M_PER_MIN = 250.0

# Programmed feed rate along the path, in millimetres per minute. The speed the
# machine would hold if no curvature or acceleration limit bound it.
DEFAULT_FEED_RATE_MM_PER_MIN = 2000.0

# Rapid-traverse rate, in millimetres per minute, for the non-cutting moves.
DEFAULT_RAPID_RATE_MM_PER_MIN = 15000.0

# Path-normal acceleration limit, in millimetres per second squared. This is what
# bounds feed through a curve: v <= sqrt(a_max / kappa). 2000 mm/s^2 is a
# mid-range figure for a rigid small-format machining centre.
DEFAULT_MAX_ACCELERATION_MM_PER_S2 = 2000.0

# Jerk limit, in millimetres per second cubed, used to charge a cost to each
# tangent discontinuity where the machine cannot pass through at speed.
DEFAULT_MAX_JERK_MM_PER_S3 = 50000.0

# Spindle speed in revolutions per minute. With the tool diameter this fixes the
# surface cutting speed V that a Taylor life estimate is a function of; 10000 rpm
# is a mid-range figure for a small-format machining centre on a small cutter.
DEFAULT_SPINDLE_RPM = 10000.0


@dataclass(frozen=True)
class MaterialModel:
    """Workpiece and cutting-edge coefficients.

    Attributes:
        feed_per_tooth_mm: Programmed feed per tooth, the scale that turns the
            geometry-derived chip-thickness ratio into a length.
        edge_radius_mm: Cutting-edge radius of the tool.
        min_chip_fraction_of_edge_radius: Ratio below which the edge rubs
            instead of cutting.
        taylor_exponent: The exponent n in ``V * T^n = C``.
        taylor_constant_m_per_min: The constant C, for V in m/min, T in minutes.
    """

    feed_per_tooth_mm: float
    edge_radius_mm: float
    min_chip_fraction_of_edge_radius: float
    taylor_exponent: float
    taylor_constant_m_per_min: float

    @property
    def min_chip_thickness_mm(self) -> float:
        """The rubbing floor ``h_min``, as a length."""
        return self.min_chip_fraction_of_edge_radius * self.edge_radius_mm

    @classmethod
    def build(
        cls,
        feed_per_tooth_mm: float = DEFAULT_FEED_PER_TOOTH_MM,
        edge_radius_mm: float = DEFAULT_EDGE_RADIUS_MM,
        min_chip_fraction_of_edge_radius: float = DEFAULT_MIN_CHIP_FRACTION_OF_EDGE_RADIUS,
        taylor_exponent: float = DEFAULT_TAYLOR_EXPONENT,
        taylor_constant_m_per_min: float = DEFAULT_TAYLOR_CONSTANT_M_PER_MIN,
    ) -> "MaterialModel":
        """Validate every coefficient and return a frozen model.

        Args:
            feed_per_tooth_mm: Feed per tooth; must be finite and positive.
            edge_radius_mm: Edge radius; must be finite and positive.
            min_chip_fraction_of_edge_radius: Rubbing ratio; must lie in (0, 1].
            taylor_exponent: Taylor n; must lie in (0, 1).
            taylor_constant_m_per_min: Taylor C; must be finite and positive.

        Returns:
            The validated model.

        Raises:
            InvalidMaterialModelError: A coefficient is NaN or out of range.
        """
        _require_positive(feed_per_tooth_mm, "feed_per_tooth_mm", InvalidMaterialModelError)
        _require_positive(edge_radius_mm, "edge_radius_mm", InvalidMaterialModelError)
        _require_positive(taylor_constant_m_per_min, "taylor_constant_m_per_min", InvalidMaterialModelError)
        if not math.isfinite(min_chip_fraction_of_edge_radius) or not (0.0 < min_chip_fraction_of_edge_radius <= 1.0):
            raise InvalidMaterialModelError(f"min_chip_fraction_of_edge_radius must lie in (0, 1], got {min_chip_fraction_of_edge_radius!r}.")
        if not math.isfinite(taylor_exponent) or not (0.0 < taylor_exponent < 1.0):
            raise InvalidMaterialModelError(f"taylor_exponent must lie in (0, 1), got {taylor_exponent!r}.")
        return cls(
            feed_per_tooth_mm=feed_per_tooth_mm,
            edge_radius_mm=edge_radius_mm,
            min_chip_fraction_of_edge_radius=min_chip_fraction_of_edge_radius,
            taylor_exponent=taylor_exponent,
            taylor_constant_m_per_min=taylor_constant_m_per_min,
        )


@dataclass(frozen=True)
class MachineModel:
    """Machine feed and dynamics limits.

    Attributes:
        feed_rate_mm_per_min: Programmed feed along the path.
        rapid_rate_mm_per_min: Traverse rate for non-cutting moves.
        max_acceleration_mm_per_s2: Path-normal acceleration limit, which bounds
            feed through a curve by ``v <= sqrt(a_max / kappa)``.
        max_jerk_mm_per_s3: Jerk limit, used to cost each tangent break.
        spindle_rpm: Spindle speed, which with the tool diameter fixes the
            surface cutting speed a Taylor life estimate depends on.
    """

    feed_rate_mm_per_min: float
    rapid_rate_mm_per_min: float
    max_acceleration_mm_per_s2: float
    max_jerk_mm_per_s3: float
    spindle_rpm: float

    @property
    def feed_mm_per_s(self) -> float:
        """The programmed feed in millimetres per second."""
        return self.feed_rate_mm_per_min / 60.0

    @property
    def rapid_mm_per_s(self) -> float:
        """The rapid rate in millimetres per second."""
        return self.rapid_rate_mm_per_min / 60.0

    def curvature_limited_feed_mm_per_s(self, curvature: float) -> float:
        """The largest feed the acceleration limit allows at *curvature*.

        The standard normal-acceleration bound ``v <= sqrt(a_max / kappa)``,
        clamped to the programmed feed. A straight move (``kappa == 0``) is
        unbounded by it and returns the programmed feed.

        Args:
            curvature: Path curvature in reciprocal millimetres; non-negative.

        Returns:
            The feed in millimetres per second.
        """
        if curvature <= 0.0:
            return self.feed_mm_per_s
        return min(self.feed_mm_per_s, math.sqrt(self.max_acceleration_mm_per_s2 / curvature))

    @classmethod
    def build(
        cls,
        feed_rate_mm_per_min: float = DEFAULT_FEED_RATE_MM_PER_MIN,
        rapid_rate_mm_per_min: float = DEFAULT_RAPID_RATE_MM_PER_MIN,
        max_acceleration_mm_per_s2: float = DEFAULT_MAX_ACCELERATION_MM_PER_S2,
        max_jerk_mm_per_s3: float = DEFAULT_MAX_JERK_MM_PER_S3,
        spindle_rpm: float = DEFAULT_SPINDLE_RPM,
    ) -> "MachineModel":
        """Validate every limit and return a frozen model.

        Args:
            feed_rate_mm_per_min: Programmed feed; must be finite and positive.
            rapid_rate_mm_per_min: Rapid rate; must be finite and positive.
            max_acceleration_mm_per_s2: Acceleration limit; finite and positive.
            max_jerk_mm_per_s3: Jerk limit; finite and positive.
            spindle_rpm: Spindle speed; must be finite and positive.

        Returns:
            The validated model.

        Raises:
            InvalidMachineModelError: A limit is NaN or non-positive.
        """
        _require_positive(feed_rate_mm_per_min, "feed_rate_mm_per_min", InvalidMachineModelError)
        _require_positive(rapid_rate_mm_per_min, "rapid_rate_mm_per_min", InvalidMachineModelError)
        _require_positive(max_acceleration_mm_per_s2, "max_acceleration_mm_per_s2", InvalidMachineModelError)
        _require_positive(max_jerk_mm_per_s3, "max_jerk_mm_per_s3", InvalidMachineModelError)
        _require_positive(spindle_rpm, "spindle_rpm", InvalidMachineModelError)
        return cls(
            feed_rate_mm_per_min=feed_rate_mm_per_min,
            rapid_rate_mm_per_min=rapid_rate_mm_per_min,
            max_acceleration_mm_per_s2=max_acceleration_mm_per_s2,
            max_jerk_mm_per_s3=max_jerk_mm_per_s3,
            spindle_rpm=spindle_rpm,
        )


def _require_positive(value: float, name: str, error: type) -> None:
    """Raise *error* unless *value* is finite and strictly positive.

    Args:
        value: The coefficient to check.
        name: Its attribute name, for the message.
        error: The named exception class to raise.

    Raises:
        BenchmarkError: The subclass named by *error*, when the value is NaN,
            infinite, zero, or negative.
    """
    if not math.isfinite(value) or value <= 0.0:
        raise error(f"{name} must be finite and positive, got {value!r}.")
