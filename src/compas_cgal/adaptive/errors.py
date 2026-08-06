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


class InvalidMatSamplingPolicyError(ValueError):
    """MAT proposal sampling omits or corrupts a cursor-affecting bound."""


class InvalidMiddleCurveSpanError(ValueError):
    """A candidate span does not advance between two owned MAT cursors."""


class InvalidMiddleCurveCursorError(ValueError):
    """A derived cursor cannot prove one exact candidate-span lineage."""


class InvalidMathsmProposalError(ValueError):
    """A one-sided MATHSM proposal contradicts its MAT/site inputs."""


class InvalidCandidateLatticeError(ValueError):
    """A finite candidate cell omits or contradicts bound decision state."""


class InvalidZeroGuideCandidateError(ValueError):
    """A zero-guide link candidate contradicts its exact proof or span state."""


class UncertifiedZeroGuideEdgeError(ValueError):
    """A link-only candidate was requested for an edge without exact proof."""


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


class InvalidZeroGuideCertificateError(ValueError):
    """A zero-guide proof record contradicts its MAT or native owner."""


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


class InvalidTraversalGraphError(ValueError):
    """An exact MAT traversal graph or directed route is inconsistent."""


class InvalidCausalNeckTransitError(ValueError):
    """A causal transition does not bind two owned exact neck sides."""


class AmbiguousNeckSideError(ValueError):
    """An exact route edge does not resolve one unique certified neck side."""


class OverlappingNeckTransitError(ValueError):
    """One operation would cross multiple independently active necks."""


class InvalidMatTraversalStateError(ValueError):
    """Global MAT traversal state violates its exact graph authority."""


class StaleTraversalCursorError(ValueError):
    """A candidate does not begin at the active exact global cursor."""


class TerminalTraversalCursorError(RuntimeError):
    """A candidate attempts to advance an already terminal MAT edge."""


class NonterminalMatTraversalError(RuntimeError):
    """Global MAT traversal is queried or switched before terminal state."""


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


class ReplayInputMismatchError(ValueError):
    """Submitted replay inputs differ from the authenticated input root."""


class ReplayGrammarError(ValueError):
    """An operation sequence violates the closed Phase-1 grammar."""


class ReplayContinuityError(ValueError):
    """A lateral operation breaks exact phase or cut-plane continuity."""


class ReplayCutDirectionError(ValueError):
    """A circle orientation contradicts its material side and cut policy."""


class ReplayCandidateError(ValueError):
    """A recorded circle is not one unique finite-lattice candidate."""


class ReplayZeroGuideCandidateError(ValueError):
    """An advancing segment lacks one unique fresh zero-guide candidate."""


class ReplayPairingError(ValueError):
    """A link does not belong exactly to its following circle candidate."""


class ReplayEffectiveCapError(ValueError):
    """A recorded cap decision contradicts fresh neck passage state."""


class ReplayTraversalError(ValueError):
    """A recorded cursor transition contradicts the fresh MAT graph."""


class InvalidReplayTraceError(ValueError):
    """Fresh replay evidence is incomplete, cross-wired, or out of order."""


class InvalidReplayCertificateError(ValueError):
    """A replay certificate omits or contradicts recomputed proof state."""


class InvalidBoundaryVertexIdentityError(ArtifactIdentityError):
    """A structural boundary-vertex identity is invalid or ambiguous."""


class InvalidEffectiveCapDecisionError(ArtifactIdentityError):
    """An effective-cap decision is not a closed, exact policy result."""


class InvalidTraversalDecisionError(ArtifactIdentityError):
    """A traversal decision does not encode one legal cursor transition."""


class InvalidOperationIdentityError(ArtifactIdentityError):
    """A canonical operation omits or contradicts an owned decision."""


class InvalidAdvanceSegmentOperationError(ArtifactIdentityError):
    """An advancing segment violates its one-motion traversal contract."""


class InvalidReachableDomainCertificateError(ValueError):
    """Reachable-domain reconstruction evidence violates its exact contract."""


class InvalidCoverageCertificateError(ValueError):
    """An exact coverage certificate violates its owned structural contract."""


class InvalidCoverageSweepError(ValueError):
    """A native coverage transition diverges from its exact Python motion."""


class IncompletePocketCoverageError(RuntimeError):
    """Exact reachable material remains after the owned sweep lineage."""


class InvalidGenerationStateError(ValueError):
    """An authoritative generation snapshot violates cross-owner chronology."""


class CandidateStateMismatchError(ValueError):
    """A candidate transition does not begin at the authoritative state."""


class InvalidCandidateTransactionError(ValueError):
    """Accepted candidate evidence is malformed, cross-wired, or reordered."""


class StaleCandidateTransactionError(RuntimeError):
    """A transaction parent digest no longer names authoritative state."""


class InvalidZeroGuideTransactionError(ValueError):
    """A zero-guide transaction contradicts its candidate or proof lineage."""


class StaleZeroGuideTransactionError(RuntimeError):
    """A zero-guide transaction no longer names its physical parent or cursor."""


class InvalidInitialCandidateEvaluatorError(ValueError):
    """An entry-launch evaluator has foreign or incomplete root authority."""


class InitialCandidateStateMismatchError(ValueError):
    """An entry candidate does not begin at the seeded global MAT state."""


class InitialCandidatePhaseError(ValueError):
    """An entry candidate phase differs from the qualified bore center."""


class InvalidInitialCandidateTransactionError(ValueError):
    """Entry-launch evidence is malformed, cross-wired, or out of order."""


class StaleInitialCandidateTransactionError(RuntimeError):
    """An entry transaction no longer names the authoritative traversal parent."""


class InvalidTraversalCommitError(ValueError):
    """A continuation commit cross-wires its physical or global state axis."""


class InvalidCandidateFamilyError(ValueError):
    """A finite candidate family is foreign, duplicated, or out of order."""


class EngagementCapInfeasibleError(RuntimeError):
    """Every non-neck candidate is proved to exceed its engagement cap."""


class NeckTooTightError(RuntimeError):
    """Every candidate for one causal neck crossing is proved cap-infeasible."""


class NoFeasibleCandidateError(RuntimeError):
    """A finite candidate family is exhausted by mixed exact rejection."""


class CandidateSelectionError(ValueError):
    """Accepted transactions cannot form one deterministic winner set."""
