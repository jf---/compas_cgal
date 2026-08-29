"""Named failure modes for measurement-claim production and validation."""

from tools.measurement_artifact import MeasurementArtifactError


class MeasurementClaimError(MeasurementArtifactError): ...


class UnknownMeasurementClaimCaseError(MeasurementClaimError): ...


class InvalidMeasurementClaimConfigError(MeasurementClaimError): ...


class ProbeInstrumentationContractError(MeasurementClaimError): ...


class InvalidMeasurementClaimPayloadError(MeasurementClaimError): ...


class MeasurementClaimChildError(MeasurementClaimError): ...
