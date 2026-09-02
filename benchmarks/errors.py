from __future__ import annotations


class BenchmarkError(Exception):
    """Base class for every benchmark-corpus failure mode."""


class NonPositiveToolError(BenchmarkError):
    """A tool diameter was zero, negative, or NaN."""


class InvalidCapError(BenchmarkError):
    """An engagement cap fell outside the exact kernel's contract of (0, 180] degrees."""


class DegeneratePocketError(BenchmarkError):
    """A generated pocket has zero area, fewer than three vertices, or no room for the tool."""


class PocketNotSimpleError(BenchmarkError):
    """A pocket boundary self-intersects, so it is not a valid general polygon."""


class UnreplayableOperationError(BenchmarkError):
    """A toolpath operation lies outside the cut-plane depletion model (e.g. a 3D ramp)."""


class UnsampleableMotionError(BenchmarkError):
    """A cut motion's geometry carries no cutter-centre path to sample positions from."""


class MalformedRecordError(BenchmarkError):
    """A serialised measurement does not carry exactly the columns the schema declares."""


class InvalidSideCountError(BenchmarkError):
    """A family was asked for a polygon with fewer than three sides."""


class InvalidArcRatioError(BenchmarkError):
    """An arc fraction fell outside the closed unit interval it is defined on."""


class ImpassableNeckError(BenchmarkError):
    """A neck was narrower than the tool, so no toolpath can reach through it."""


class UnpinchedChannelError(BenchmarkError):
    """A neck was as wide as the pocket it pinches, so the instance carries no neck."""


class MissingSweepParameterError(BenchmarkError):
    """A record handed to a sweep analysis does not carry that sweep's parameter."""


class InvalidDecimalsError(BenchmarkError):
    """A coordinate precision was negative, so it names no rounding."""


class InvalidIslandCountError(BenchmarkError):
    """An island grid was asked for a non-positive number of rows or columns."""


class CrowdedIslandGridError(BenchmarkError):
    """An island grid leaves channels too narrow for the tool to machine."""


class UnmeasurableOperationLengthError(BenchmarkError):
    """A toolpath operation's geometry is not a primitive whose length is defined."""


class UnknownCorpusError(BenchmarkError):
    """A corpus was requested by a name no family answers to."""


class MissingExternalDirectoryError(BenchmarkError):
    """The external corpus was requested without the operator-supplied profile directory."""


class UnplottableGeometryError(BenchmarkError):
    """A toolpath operation carries a primitive with no drawable path."""


class UnplottableBoundaryError(BenchmarkError):
    """A pocket boundary has too few points to close into a ring."""


class UnknownOperationClassError(BenchmarkError):
    """A toolpath operation names an operation class the drawing has no mark for."""


class UnknownColourModeError(BenchmarkError):
    """A drawing was asked for a colour mode no encoder answers to."""


class EmptyToolpathError(BenchmarkError):
    """A toolpath carries no operations, so there is nothing to draw."""


class EmptyComparisonError(BenchmarkError):
    """A figure comparing toolpaths was given no panels to compare."""


class MissingEngagementDataError(BenchmarkError):
    """Engagement colouring was asked for without the measurements it draws."""


class EngagementLengthMismatchError(BenchmarkError):
    """A drawing was given a number of engagement measurements other than one per operation."""


class MissingToolDiameterError(BenchmarkError):
    """A swept-tool envelope was asked for without a tool diameter to sweep."""


class InvalidBandCountError(BenchmarkError):
    """A colour ramp was asked to split a range into fewer than one band."""


class AmbiguousPanelError(BenchmarkError):
    """A drawing with several panels was asked for its single panel."""


class InvalidGridResolutionError(BenchmarkError):
    """A coverage grid was requested with fewer than one sample along an axis."""


class CoarseCoverageGridError(BenchmarkError):
    """A coverage grid's cell is too coarse to resolve residue at the instance's tool size."""


class InvalidMotionSampleCountError(BenchmarkError):
    """A cut motion was to be probed at fewer than one cutter position."""


class EmptyReachableRegionError(BenchmarkError):
    """No coverage sample landed inside the tool-reachable region, so no coverage fraction is defined."""


class ZeroLengthToolpathError(BenchmarkError):
    """A toolpath's operations sum to zero length, so no length fraction is defined."""


class MissingMaterialModelError(BenchmarkError):
    """A metric that depends on material coefficients was asked for without a `MaterialModel`."""


class MissingMachineModelError(BenchmarkError):
    """A metric that depends on machine limits was asked for without a `MachineModel`."""


class InvalidMaterialModelError(BenchmarkError):
    """A material coefficient was NaN, non-positive, or outside the range its definition allows."""


class InvalidMachineModelError(BenchmarkError):
    """A machine limit was NaN, non-positive, or outside the range its definition allows."""


class NoCornerDefectError(BenchmarkError):
    """A figure that annotates the corner defect found no degenerate loop to annotate."""


class InvalidPublishedPrimitiveError(BenchmarkError):
    """A published vector primitive is non-finite or degenerate."""


class DisconnectedPublishedBoundaryError(BenchmarkError):
    """Published boundary primitives do not form one closed cycle."""


class UnresolvedPublishedCurveError(BenchmarkError):
    """A published cubic cannot be reconstructed inside its fidelity bound."""


class InvalidReferenceProjectionError(BenchmarkError):
    """A polygon projection violates its declared chord-deviation contract."""


class InvalidReferenceReconstructionError(BenchmarkError):
    """A reference reconstruction is empty, disconnected, or carries an invalid proof bound."""


class UnknownHeldReferenceCaseError(BenchmarkError):
    """A requested Held-Pfeiffer case has no committed corpus document."""


class MalformedHeldReferenceCaseError(BenchmarkError):
    """A Held-Pfeiffer document violates its closed schema or geometric evidence."""


class UnsupportedHeldReferenceVersionError(BenchmarkError):
    """A Held-Pfeiffer document uses an unsupported schema version."""


class UnsupportedPdfBoundaryOperatorError(BenchmarkError):
    """A selected publisher path uses an unsupported drawing operator."""


class MissingPublishedToolCircleError(BenchmarkError):
    """A figure crop contains no unambiguous depicted tool circle."""


class AmbiguousPublishedBoundaryError(BenchmarkError):
    """A figure crop contains more than one valid boundary selection."""


class MissingPublishedBoundaryMarkerError(BenchmarkError):
    """A selected boundary endpoint has no unique published marker."""


class InvalidReferenceOverlayError(BenchmarkError):
    """A Held reference overlay lacks valid, case-matched semantic geometry."""


class MissingFigureAxesError(BenchmarkError):
    """A publisher plot lacks the axes required for geometric registration."""


class EmptyFigureColourSamplesError(BenchmarkError):
    """A publisher panel contains no coloured tool-centre samples."""


class InvalidFigureInwardOffsetError(BenchmarkError):
    """A certified projection has no valid one-radius comparator component."""


class AmbiguousFigurePanelRegistrationError(BenchmarkError):
    """Publisher axes admit no unique panel-to-world registration."""


class InvalidHeldOperationSnapshotError(BenchmarkError):
    """A snapshotted operation contains malformed geometric or motion data."""


class InvalidHeldPathEvidenceError(BenchmarkError):
    """Typed values or operation coverage in Held path evidence are invalid."""


class InvalidHeldPathReportContextError(BenchmarkError):
    """Held report invocation metadata is empty, non-UTC, or malformed."""


class ContradictoryEngagementEvidenceError(BenchmarkError):
    """Guarded replay and sampled exact-predicate evidence contradict."""


class ContradictoryPathQualityEvidenceError(BenchmarkError):
    """Path-quality attribution does not reduce to its aggregate value."""


class MutatedHeldToolpathError(BenchmarkError):
    """A path-replay consumer changed the generated operation stream."""


class HeldPathNotEligibleForPostQualificationError(BenchmarkError):
    """A complete characterization has at least one open Phase 1 criterion."""


class UnexpectedHeldPathCaseError(BenchmarkError):
    """Characterization received a Held case other than Figure 5."""
