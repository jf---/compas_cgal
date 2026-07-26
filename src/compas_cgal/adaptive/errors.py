from compas_cgal import _containment_2
from compas_cgal import _coverage_2
from compas_cgal import _stock_2

ExactDepletionConstructionError = _stock_2.ExactDepletionConstructionError
ExactDepletionCenterLimitError = _stock_2.ExactDepletionCenterLimitError
ReachableDomainConstructionError = _coverage_2.ReachableDomainConstructionError
InvalidReachableDomainInputError = _coverage_2.InvalidReachableDomainInputError
ReachableArrangementTopologyError = _coverage_2.ReachableArrangementTopologyError
PocketNotMachinableError = _coverage_2.PocketNotMachinableError
ReachableMaterialContainmentError = _coverage_2.ReachableMaterialContainmentError
InvalidCoverageGeometryError = _coverage_2.InvalidCoverageGeometryError
CoverageTransitionError = _coverage_2.CoverageTransitionError
ContainmentConstructionError = _containment_2.ContainmentConstructionError


class InvalidUnitValueError(ValueError):
    """A typed unit or framed coordinate violates its numeric invariant."""


class InvalidEngagementCapError(ValueError):
    """An engagement cap cannot be represented by the native exact boundary."""


class InvalidMotionCertificateError(ValueError):
    """A motion-certifier request violates its typed local contract."""


class EngagementCapExceededError(RuntimeError):
    """An exact motion oracle proves that the effective cap is exceeded."""


class UnresolvedMotionEventError(RuntimeError):
    """An exact motion oracle cannot reconstruct a complete event partition."""


class DegenerateSegmentMotionError(ValueError):
    """An exact segment motion has no progress."""


class DegenerateCircleMotionError(ValueError):
    """An exact circle motion has no nonzero phase vector."""


class InvalidCandidatePolicyError(ValueError):
    """A candidate lattice is not finite and deterministic."""


class InvalidMiddleCurveSpanError(ValueError):
    """A candidate span does not advance between two owned MAT cursors."""


class InvalidMathsmProposalError(ValueError):
    """A one-sided MATHSM proposal contradicts its MAT/site inputs."""


class InvalidCandidateLatticeError(ValueError):
    """A finite candidate cell omits or contradicts bound decision state."""


class InvalidContainmentCertificateError(ValueError):
    """A native exact-containment record contradicts its typed motion."""


class GougeContainmentError(RuntimeError):
    """An exact cutter sweep is not contained in the design pocket."""


class InvalidPreclearedEntryError(ValueError):
    """A declared precleared bore cannot launch the certified programme."""


class InvalidNeckPolicyError(ValueError):
    """A neck classification or effective-cap mapping is incomplete or invalid."""


class InvalidMedialAxisProjectionError(ValueError):
    """A native MAT projection contradicts its exact topology owner."""


class InvalidNeckEvidenceError(ValueError):
    """Native neck evidence or exact classification is structurally inconsistent."""


class InvalidNeckPassageTransitionError(ValueError):
    """A neck passage decision does not advance its owned orientation once."""


class TerminalNeckPassageError(RuntimeError):
    """A terminal oriented neck passage cannot propose another cut."""


class InvalidDepletionPolicyError(ValueError):
    """A depletion construction has an invalid bound."""


class InvalidDepletionTraceError(ValueError):
    """A native exact-depletion trace violates its structural contract."""


class InvalidDepletionWitnessError(ValueError):
    """A depletion witness omits or corrupts an immutable identity input."""


class InvalidStockAreaError(ValueError):
    """A stock-area owner received an aliased or invalid stock lifetime."""


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


class InvalidInputIdentityError(ArtifactIdentityError):
    """A Phase-1 input identity omits or cross-wires a validated input."""


class InvalidBoundaryVertexIdentityError(ArtifactIdentityError):
    """A structural boundary-vertex identity is invalid or ambiguous."""


class InvalidEffectiveCapDecisionError(ArtifactIdentityError):
    """An effective-cap decision is not a closed, exact policy result."""


class InvalidTraversalDecisionError(ArtifactIdentityError):
    """A traversal decision does not encode one legal cursor transition."""


class InvalidOperationIdentityError(ArtifactIdentityError):
    """A canonical operation omits or contradicts an owned decision."""


class InvalidReachableDomainCertificateError(ValueError):
    """Reachable-domain reconstruction evidence violates its exact contract."""


class InvalidCoverageCertificateError(ValueError):
    """An exact coverage certificate violates its owned structural contract."""


class InvalidCoverageSweepError(ValueError):
    """A native coverage transition diverges from its exact Python motion."""


class IncompletePocketCoverageError(RuntimeError):
    """Exact reachable material remains after the owned sweep lineage."""
