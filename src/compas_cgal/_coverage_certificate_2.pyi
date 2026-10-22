from collections.abc import Sequence
from typing import TypedDict

RationalText = tuple[str, str]
PointText = tuple[RationalText, RationalText]
EdgeRecord = tuple[int, int, int, int, bool, list[int], list[int], RationalText, PointText, PointText, bool]
WitnessRecord = tuple[int, int, int, int, RationalText, PointText]
ConnectorRow = tuple[int, int, RationalText, int, RationalText]

class InvalidCentreDomainError(ValueError): ...
class InvalidCertificateInputError(ValueError): ...
class CoverageCertificateNativeError(RuntimeError): ...
class NonFiniteBinary64Error(ValueError): ...
class StartWitnessSearchError(CoverageCertificateNativeError): ...
class SignRuleContractError(CoverageCertificateNativeError): ...
class SeparationBoundExceededError(CoverageCertificateNativeError): ...
class CoincidentTargetSupportError(CoverageCertificateNativeError): ...
class ConnectorBoundaryCoverageError(CoverageCertificateNativeError): ...
class BoundaryChainError(CoverageCertificateNativeError): ...

def contract_probe_start_witness(cx: float, cy: float, radius: float, coincident: bool) -> tuple[int, int]: ...
def contract_probe_sign_check(claimed: int, lo: float, hi: float) -> None: ...
def contract_probe_sign_rational_two(u: RationalText, v: RationalText, p: RationalText, w: RationalText, q: RationalText) -> int: ...
def contract_probe_rational_between(
    lo_m: RationalText, lo_n: RationalText, lo_d: RationalText, hi_m: RationalText, hi_n: RationalText, hi_d: RationalText
) -> RationalText: ...
def contract_probe_vertex_disk(x0: RationalText, x1: RationalText, y0: RationalText, y1: RationalText, delta: RationalText) -> tuple[bool, bool]: ...

class CircleCoreOutput(TypedDict):
    exposed: list[EdgeRecord]
    counters: dict[str, int]
    region_count: int
    piece_count: int
    elapsed_s: float
    phase_s: list[float]

def validate_centre_domain(outer: Sequence[tuple[float, float]], holes: Sequence[Sequence[tuple[float, float]]]) -> None: ...
def run_circle_core(
    centres: Sequence[tuple[float, float]],
    path_radii: Sequence[float],
    tool_radius: float,
    outer: Sequence[tuple[float, float]],
    holes: Sequence[Sequence[tuple[float, float]]],
    vertex_disks: Sequence[Sequence[RationalText]],
    exact_only: bool,
    audit: bool,
) -> CircleCoreOutput: ...

class CertificateOutput(TypedDict):
    exposed: list[WitnessRecord]
    undecided: int
    counters: dict[str, int]
    region_count: int
    piece_count: int
    edge_count: int
    elapsed_s: float
    phase_s: list[float]

def run_certificate(
    centres: Sequence[tuple[float, float]],
    path_radii: Sequence[float],
    tool_radius: float,
    outer: Sequence[tuple[float, float]],
    holes: Sequence[Sequence[tuple[float, float]]],
    connectors: Sequence[ConnectorRow],
    vertex_disks: Sequence[Sequence[RationalText]],
    exact_only: bool,
    audit: bool,
) -> CertificateOutput: ...

PrimitiveRow = tuple[bool, int, list[RationalText]]
PieceRow = tuple[bool, int, list[RationalText], list[RationalText], list[RationalText]]
ExactPieceRow = tuple[int, bool, list[RationalText], list[RationalText]]

def run_exact_target_certificate(
    centres: Sequence[tuple[float, float]],
    path_radii: Sequence[float],
    tool_radius: float,
    outer: Sequence[tuple[float, float]],
    holes: Sequence[Sequence[tuple[float, float]]],
    connectors: Sequence[ConnectorRow],
    boundary: Sequence[PrimitiveRow],
    pieces: Sequence[PieceRow],
    vertex_disks: Sequence[Sequence[RationalText]],
    exact_only: bool,
    audit: bool,
    exact_connectors: Sequence[Sequence[ExactPieceRow]] = ...,
) -> CertificateOutput: ...
def contract_probe_probe_step(w: PointText, d: PointText, circles: Sequence[tuple[RationalText, RationalText, RationalText]], exact_only: bool) -> RationalText:
    """The side probe's step h at ``w`` along ``d`` against rational support circles (centre x, y, squared radius)."""

def contract_probe_exact_target_contains(
    tool_radius: float,
    boundary: Sequence[PrimitiveRow],
    pieces: Sequence[PieceRow],
    vertex_disks: Sequence[Sequence[RationalText]],
    points: Sequence[PointText],
    exact_only: bool,
) -> list[bool]: ...
