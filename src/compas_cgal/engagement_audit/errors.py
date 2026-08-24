"""Named failures at authoritative engagement-audit boundaries."""


class EngagementAuditError(ValueError):
    """Base class for engagement-audit contract failures."""


class InvalidBuildIdentityError(EngagementAuditError):
    """Build provenance is incomplete or noncanonical."""


class DuplicateBuildComponentError(InvalidBuildIdentityError):
    """Two build components claim the same stable domain."""


class InvalidMeasuredOperationAuditError(EngagementAuditError):
    """A measured operation omits or corrupts native evidence."""


class InvalidMotionVerdictError(InvalidMeasuredOperationAuditError):
    """A native motion verdict is outside the closed three-way domain."""


class InvalidNonEngagingOperationAuditError(EngagementAuditError):
    """A non-engaging operation carries an invalid geometric reason."""


class InvalidAuthenticatedLateralOperationError(EngagementAuditError):
    """An authenticated lateral operation omits valid source motion identity."""


class InvalidAuthenticatedPlungeOperationError(EngagementAuditError):
    """An authenticated plunge omits its exact opaque native motion."""


class InvalidAuthenticatedNonEngagingOperationError(EngagementAuditError):
    """An authenticated non-engaging operation carries unsupported native geometry."""


class InvalidEngagementAuditInputError(EngagementAuditError):
    """An engagement-audit input is incomplete or noncanonical."""


class InvalidAuditDepletionPolicyError(InvalidEngagementAuditInputError):
    """The audit input does not own one exact depletion policy."""


class InvalidNativeAuditPolicyError(InvalidEngagementAuditInputError):
    """The native audit policy rejected an authoritative input."""


class InconsistentEngagementCapSurrogateError(InvalidNativeAuditPolicyError):
    """The supplied cap surrogate differs from the native cap observation."""


class InvalidNativeAuditRequestIdentityError(InvalidEngagementAuditInputError):
    """Native stock, policy, or motion identity construction failed."""


class InvalidAuditDecisionLimitsError(InvalidEngagementAuditInputError):
    """The audit input does not own one exact sealed decision-limits value."""


class InvalidAuditSpatialFloorError(InvalidAuditDecisionLimitsError):
    """The audit decision spatial floor is not a positive finite length."""


class InvalidAuditDecisionDepthError(InvalidAuditDecisionLimitsError):
    """The audit decision depth is outside its sealed non-negative range."""


class InvalidAuditDecisionNodeLimitError(InvalidAuditDecisionLimitsError):
    """The audit decision node limit is outside its sealed positive range."""


class EmptyToolpathAuditError(InvalidEngagementAuditInputError):
    """An audit was requested for an empty operation stream."""


class InvalidAuditOperationError(InvalidEngagementAuditInputError):
    """A legacy toolpath operation cannot enter the authenticated stream."""


class NonFiniteAuditGeometryError(InvalidAuditOperationError):
    """Toolpath geometry contains a nonfinite binary64 value."""


class UnsupportedAuditGeometryError(InvalidEngagementAuditInputError):
    """Toolpath geometry has no authoritative audit semantics."""


class MultipleCutPlaneError(UnsupportedAuditGeometryError):
    """A lateral operation lies on a plane other than the declared cut plane."""


class ContradictoryOperationRoleError(InvalidEngagementAuditInputError):
    """An operation label contradicts its geometry-derived role."""


class ContradictoryOperationOrientationError(InvalidEngagementAuditInputError):
    """An operation orientation contradicts its native arc traversal."""


class InvalidAuditPlaneError(InvalidEngagementAuditInputError):
    """The declared clearance plane is not strictly above the cut plane."""
