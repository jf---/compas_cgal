class InvalidUnitValueError(ValueError):
    """A typed unit or framed coordinate violates its numeric invariant."""


class InvalidEngagementCapError(ValueError):
    """An engagement cap cannot be represented by the native exact boundary."""


class DegenerateSegmentMotionError(ValueError):
    """An exact segment motion has no progress."""


class DegenerateCircleMotionError(ValueError):
    """An exact circle motion has no nonzero phase vector."""


class InvalidCandidatePolicyError(ValueError):
    """A candidate lattice is not finite and deterministic."""


class InvalidNeckPolicyError(ValueError):
    """A neck classification or effective-cap mapping is incomplete or invalid."""


class InvalidDepletionPolicyError(ValueError):
    """A depletion construction has an invalid bound."""


class InvalidTraversalPolicyError(ValueError):
    """A traversal order or forward window is invalid."""


class InvalidCutDirectionPolicyError(ValueError):
    """A cut-direction intent or material-side decision is invalid."""


class CanonicalEncodingError(ValueError):
    """A value cannot be represented by the frozen canonical byte grammar."""


class ArtifactIdentityError(ValueError):
    """An identity input or decision record violates its structural contract."""


class InvalidComponentIdentityError(ArtifactIdentityError):
    """A component identity omits or corrupts a required build input."""


class InvalidBoundaryVertexIdentityError(ArtifactIdentityError):
    """A structural boundary-vertex identity is invalid or ambiguous."""


class InvalidEffectiveCapDecisionError(ArtifactIdentityError):
    """An effective-cap decision is not a closed, exact policy result."""


class InvalidTraversalDecisionError(ArtifactIdentityError):
    """A traversal decision does not encode one legal cursor transition."""


class InvalidOperationIdentityError(ArtifactIdentityError):
    """A canonical operation omits or contradicts an owned decision."""
