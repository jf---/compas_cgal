"""Units of the helical entry: the ramp angle, its exact rational pitch surrogate, and the declared entry ramp.

The ramp-angle cap is transcendental user intent (``tan``). Per the exactness
doctrine it is not enforced transcendentally: `HelixPitchPerRadius` is the
surrogate ``K = 2*pi*tan(alpha)`` computed ONCE as a double and injected
exactly (`HelixPitchPerRadius.exact`). Every ramp decision is an exact statement
about ``K``; the few-ulp gap between ``K`` and ``2*pi*tan(alpha)`` is this API's
semantic, not an in-core correction.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction
from typing import Self

# Default machine ramp-angle cap, degrees (plan D1). PLACEHOLDER until a
# datasheet value for the gate tool is supplied; nothing derives it.
DEFAULT_MAX_RAMP_ANGLE_DEG = 2.0

# Default declared entry ramp, tool radii (USER 2026-10-16 Q4): one tool
# diameter. PLACEHOLDER declared intent, like the cap; nothing derives it.
DEFAULT_ENTRY_RAMP_TOOL_RADII = 2.0

# Open upper bound of a ramp angle, degrees: tan diverges at a vertical descent,
# which is a plunge, not a ramp.
VERTICAL_DESCENT_DEG = 90.0


class InvalidRampAngleError(ValueError):
    """A ramp angle is not a finite real in the open interval (0, 90) degrees."""


class InvalidHelixPitchError(ValueError):
    """A helix pitch-per-radius surrogate is not a finite positive real."""


class InvalidEntryRampLengthError(ValueError):
    """A declared entry ramp length is not a finite positive real number of tool radii."""


def _real(value: object) -> float | None:
    """The value as a float when it is a non-bool real number, else None."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    return float(value)


@dataclass(frozen=True)
class RampAngleDegrees:
    """Machine ramp-angle cap, degrees, in the open interval (0, 90).

    Attributes:
        value: The angle, degrees. REPORTING and surrogate construction only;
            no decision reads it (decisions read `HelixPitchPerRadius.exact`).
    """

    value: float

    def __post_init__(self) -> None:
        numeric = _real(self.value)
        if numeric is None or not math.isfinite(numeric) or not 0.0 < numeric < VERTICAL_DESCENT_DEG:
            raise InvalidRampAngleError(f"ramp angle must be finite and in (0, {VERTICAL_DESCENT_DEG}) deg; got {self.value!r}.")
        object.__setattr__(self, "value", numeric)

    @classmethod
    def build(cls, value: float) -> Self:
        """Validate a ramp angle.

        Args:
            value: Angle, degrees.

        Returns:
            The typed angle.

        Raises:
            InvalidRampAngleError: *value* is not finite or not in (0, 90).
        """
        return cls(value)


@dataclass(frozen=True)
class HelixPitchPerRadius:
    """Exact surrogate ``K`` of the ramp cap: descent per turn per footprint radius.

    A helix of footprint radius ``rho`` and pitch ``p`` (mm of descent per turn)
    ramps at ``atan(p / (2*pi*rho))``; it is within a cap ``alpha`` iff
    ``p <= 2*pi*tan(alpha) * rho``. ``K`` is that factor as one double, and the
    API contract is ``p <= exact(K) * rho`` decided over the rationals. The gap
    between ``exact(K)`` and ``2*pi*tan(alpha)`` is the rounding of the double
    expression ``2.0 * pi * tan(radians(alpha))``: ``radians``, ``tan`` and the
    products each round, so a few ulps of ``K``.

    Attributes:
        value: ``K``, dimensionless (mm of descent per turn per mm of radius).
    """

    value: float

    def __post_init__(self) -> None:
        numeric = _real(self.value)
        if numeric is None or not math.isfinite(numeric) or numeric <= 0.0:
            raise InvalidHelixPitchError(f"helix pitch per radius must be finite and positive; got {self.value!r}.")
        object.__setattr__(self, "value", numeric)

    @classmethod
    def build(cls, value: float) -> Self:
        """Validate a surrogate given directly as a double.

        Args:
            value: ``K``, dimensionless.

        Returns:
            The typed surrogate.

        Raises:
            InvalidHelixPitchError: *value* is not finite or not positive.
        """
        return cls(value)

    @classmethod
    def from_ramp_angle(cls, angle: RampAngleDegrees) -> Self:
        """The surrogate of a ramp-angle cap, ``2*pi*tan(alpha)`` as one double.

        Args:
            angle: The ramp-angle cap.

        Returns:
            The typed surrogate.

        Raises:
            InvalidHelixPitchError: The double rounds to a non-positive or
                non-finite value (unreachable for angles in (0, 90)).
        """
        return cls(2.0 * math.pi * math.tan(math.radians(angle.value)))

    @property
    def exact(self) -> Fraction:
        """``K`` injected exactly: the rational equal to the double `value`."""
        return Fraction(self.value)


# The default machine ramp surrogate: `DEFAULT_MAX_RAMP_ANGLE_DEG` through
# `HelixPitchPerRadius.from_ramp_angle` (a PLACEHOLDER, as that angle is).
DEFAULT_HELIX_PITCH_PER_RADIUS = HelixPitchPerRadius.from_ramp_angle(RampAngleDegrees(DEFAULT_MAX_RAMP_ANGLE_DEG))


@dataclass(frozen=True)
class EntryRampLength:
    """Declared entry ramp: the path length, in tool radii, over which engagement may rise from the seam to the cap.

    DECLARED USER INTENT, like the cap (USER 2026-10-16 Q4); the default
    `DEFAULT_ENTRY_RAMP_TOOL_RADII` (one tool diameter) is a placeholder value.
    The rise limit it implies at an entry distance ``s_e`` is
    ``theta_e = cap * s_e / L``; the pipeline evaluates it at the gate's
    station 1 of the smallest legal loop (``s_e = 2 pi r / 45``), so
    ``theta_e = cap * pi / 45`` at the default
    (`compas_cgal.pipeline.admissibility.entry_rise_ratio`).

    Attributes:
        tool_radii: ``L / r``, a finite positive real.
    """

    tool_radii: float

    def __post_init__(self) -> None:
        numeric = _real(self.tool_radii)
        if numeric is None or not math.isfinite(numeric) or numeric <= 0.0:
            raise InvalidEntryRampLengthError(f"entry ramp length must be a finite positive number of tool radii; got {self.tool_radii!r}.")
        object.__setattr__(self, "tool_radii", numeric)

    @classmethod
    def build(cls, tool_radii: float) -> Self:
        """Validate a declared entry ramp length.

        Args:
            tool_radii: ``L / r``.

        Returns:
            The typed length.

        Raises:
            InvalidEntryRampLengthError: *tool_radii* is not finite or not positive.
        """
        return cls(tool_radii)

    @classmethod
    def default(cls) -> Self:
        """The placeholder declared ramp, `DEFAULT_ENTRY_RAMP_TOOL_RADII` (one tool diameter)."""
        return cls(DEFAULT_ENTRY_RAMP_TOOL_RADII)
