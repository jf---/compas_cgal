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
