"""One named exception per failure mode of the test-cut post."""


class PostError(Exception):
    """Base of every post failure: the program is refused, never repaired."""


# --- setup (`PostSetup.build` and its parts) ---------------------------------


class InexactQuantityError(PostError):
    """A setup quantity is not an exact `Fraction` (or `int` count): doubles enter only at the mapping boundary."""


class NonPositiveQuantityError(PostError):
    """A length, scale, feed, speed, count or tolerance that must be positive is not."""


class DegenerateEnvelopeError(PostError):
    """A travel-envelope axis has min >= max."""


class ArcToleranceBelowQuantumError(PostError):
    """The output quantum can move an arc's radii apart by more than the controller arc tolerance."""


class NonCentreCuttingPlungeError(PostError):
    """Plunges are accepted but the cutter is not centre-cutting."""


class DepthExceedsStockError(PostError):
    """The axial depth reaches through the stock blank."""


class InvalidPlungeAllowanceError(PostError):
    """The accepted plunge count is negative."""


class ControllerResolutionError(PostError):
    """The controller's input resolution would round the emitted digits (the output quantum is not a multiple of it)."""


# --- mapping (operations -> program) ----------------------------------------


class UnsupportedOperationError(PostError):
    """An `OperationType`, or an operation/geometry pairing, the post has no mapping for."""


class UnmappedHeightError(PostError):
    """An operation height is neither the source cut plane nor the source clearance plane."""


class NonVerticalMoveError(PostError):
    """A PLUNGE or RETRACT moves in XY."""


class MissingStartTangentError(PostError):
    """A circular operation carries no start tangent to cross-check its `clockwise` flag against."""


class WindingDisagreementError(PostError):
    """A circular operation's `clockwise` flag disagrees with the turn sense of its own start tangent."""


class MissingCircleEntryError(PostError):
    """A full circle or a helix carries no exact `entry_xy`."""


class HelixOperationMismatchError(PostError):
    """A HELIX operation's `clockwise` or `entry_xy` is not its own `Helix`'s."""


class DiscontinuousPathError(PostError):
    """An operation does not start where the emitted program left the tool."""


class MissingEntryPlungeError(PostError):
    """The operation stream does not start with a PLUNGE or a HELIX entry."""


class MissingFinalRetractError(PostError):
    """The operation stream does not end with a RETRACT."""


class ToolDiameterMismatchError(PostError):
    """The setup's cutter is not the source tool scaled by the setup's scale."""


class UnsupportedDialectError(PostError):
    """A controller dialect the writer has no mapping for."""


class UnsupportedArcModeError(PostError):
    """An arc plane / centre mode the writer has no mapping for."""


class OffGridValueError(PostError):
    """A value to be written is not on the output decimal grid."""


# --- parsing (text -> program) ----------------------------------------------


class UnsupportedWordError(PostError):
    """A word outside the emitted dialect."""


class MalformedProgramError(PostError):
    """Dialect words in a structure the writer never emits."""


class RoundTripMismatchError(PostError):
    """The written text does not parse back to the program that was written."""


# --- checks (program -> verdict) --------------------------------------------


class ArcRadiusMismatchError(PostError):
    """An arc's start and end radii differ by more than the controller arc tolerance."""


class DegenerateArcError(PostError):
    """An arc of zero radius, or one that ends where it starts (a controller-dependent full turn)."""


class OutOfEnvelopeError(PostError):
    """A motion point or an arc's swept extent leaves the travel envelope."""


class UnsafeRapidError(PostError):
    """A rapid below clearance outside the named exemptions (declared approach, vertical retract)."""


class OffPlaneFeedError(PostError):
    """A feed motion that is neither an XY cut on the cut plane, a descending plunge onto it, a Z-only descent in air, nor a helical block descending between the planes."""


class PocketOutsideStockError(PostError):
    """A cut-plane motion's tool disc, or an arc's swept tool extent, leaves the stock blank in XY."""


class NonCentreCuttingHelixError(PostError):
    """A HELIX entry is posted for a cutter that is not centre-cutting (ramping cuts with the end face; fail closed pending the user)."""


class RampStartAboveClearanceError(PostError):
    """The helix ramp start is not below the clearance plane."""


class InexactProgramHelixError(PostError):
    """A program-unit helix radius is not exactly a double, so `compas_cgal.helix.ramp_turns` cannot decide its turn count exactly."""


class UnacceptedPlungeError(PostError):
    """More plunges than the setup accepts."""


class HelixStructureError(PostError):
    """A run of helical blocks is not whole same-sense turns about one centre ending on the cut plane."""


class HelixDescentMismatchError(PostError):
    """A helical block's Z descent departs from the run's uniform pitch by more than one output quantum."""


class HelixRampAngleError(PostError):
    """An emitted helical turn descends more than the setup's ramp-cap surrogate allows on its radius."""
