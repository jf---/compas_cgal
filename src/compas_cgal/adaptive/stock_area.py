import hashlib
from dataclasses import dataclass
from typing import overload

from compas_cgal import _stock_2
from compas_cgal.adaptive.canonical import ExactCenterParameterV1
from compas_cgal.adaptive.canonical import canonical_depletion_witness_bytes
from compas_cgal.adaptive.errors import InvalidDepletionTraceError
from compas_cgal.adaptive.errors import InvalidDepletionWitnessError
from compas_cgal.adaptive.errors import InvalidStockAreaError
from compas_cgal.adaptive.identity import IdentityDigest
from compas_cgal.adaptive.motion import ExactCircleMotion
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.policy import DepletionPolicy
from compas_cgal.adaptive.units import ToolRadius
from compas_cgal.stock import Stock


@dataclass(frozen=True)
class DepletionWitness:
    motion: ExactSegmentMotion | ExactCircleMotion
    policy: DepletionPolicy
    tool_radius: ToolRadius
    center_parameters: tuple[ExactCenterParameterV1, ...]
    native_strategy_version: bytes
    parent_lineage: tuple[IdentityDigest, ...]

    def __post_init__(self) -> None:
        if type(self.motion) not in (ExactSegmentMotion, ExactCircleMotion):
            raise InvalidDepletionWitnessError("depletion witness requires one exact motion type.")
        if type(self.policy) is not DepletionPolicy:
            raise InvalidDepletionWitnessError("depletion witness requires an exact depletion policy.")
        if type(self.tool_radius) is not ToolRadius:
            raise InvalidDepletionWitnessError("depletion witness requires an exact typed tool radius.")
        if type(self.center_parameters) is not tuple or not self.center_parameters:
            raise InvalidDepletionWitnessError("depletion witness requires nonempty structural center parameters.")
        if any(type(parameter) is not ExactCenterParameterV1 for parameter in self.center_parameters):
            raise InvalidDepletionWitnessError("depletion witness center parameters must use exact closed values.")
        if type(self.native_strategy_version) is not bytes or not self.native_strategy_version:
            raise InvalidDepletionWitnessError("depletion witness requires a nonempty native strategy version.")
        if type(self.parent_lineage) is not tuple or any(type(digest) is not bytes or len(digest) != 32 for digest in self.parent_lineage):
            raise InvalidDepletionWitnessError("depletion witness parent lineage must contain exact SHA-256 digests.")
        try:
            self.canonical_bytes
        except ValueError as error:
            raise InvalidDepletionWitnessError(str(error)) from None

    @property
    def canonical_bytes(self) -> bytes:
        return canonical_depletion_witness_bytes(
            motion=self.motion,
            policy=self.policy,
            tool_radius=self.tool_radius,
            center_parameters=self.center_parameters,
            native_strategy_version=self.native_strategy_version,
            parent_lineage=tuple(bytes(digest) for digest in self.parent_lineage),
        )

    @property
    def digest(self) -> IdentityDigest:
        return IdentityDigest(hashlib.sha256(self.canonical_bytes).digest())


class Stock2Area:
    def __init__(
        self,
        raw: "_stock_2.Stock2",
        lineage: tuple[DepletionWitness, ...],
    ) -> None:
        if not isinstance(raw, _stock_2.Stock2):
            raise InvalidStockAreaError("stock area requires one owned native Stock2.")
        if type(lineage) is not tuple or any(type(witness) is not DepletionWitness for witness in lineage):
            raise InvalidStockAreaError("stock area lineage must contain exact depletion witnesses.")
        self._raw = raw
        self._lineage = lineage

    @classmethod
    def build(cls, stock: Stock) -> "Stock2Area":
        if type(stock) is not Stock:
            raise InvalidStockAreaError("stock area factory requires exact Stock.")
        return cls(stock.raw.clone(), ())

    @property
    def raw(self) -> "_stock_2.Stock2":
        return self._raw

    @property
    def lineage(self) -> tuple[DepletionWitness, ...]:
        return self._lineage

    def fork(self) -> "Stock2Area":
        return Stock2Area(self._raw.clone(), self._lineage)

    @overload
    def deplete(
        self,
        motion: ExactSegmentMotion,
        tool_radius: ToolRadius,
        policy: DepletionPolicy,
    ) -> DepletionWitness: ...

    @overload
    def deplete(
        self,
        motion: ExactCircleMotion,
        tool_radius: ToolRadius,
        policy: DepletionPolicy,
    ) -> DepletionWitness: ...

    def deplete(
        self,
        motion: ExactSegmentMotion | ExactCircleMotion,
        tool_radius: ToolRadius,
        policy: DepletionPolicy,
    ) -> DepletionWitness:
        if type(tool_radius) is not ToolRadius:
            raise InvalidStockAreaError("stock depletion requires an exact typed tool radius.")
        if type(policy) is not DepletionPolicy:
            raise InvalidStockAreaError("stock depletion requires an exact depletion policy.")
        trial = self._raw.clone()
        if type(motion) is ExactSegmentMotion:
            trace = trial.subtract_exact_segment(
                motion.start.x,
                motion.start.y,
                motion.end.x,
                motion.end.y,
                tool_radius.value,
                policy.chord_bound.value,
                policy.center_count_limit,
            )
        elif type(motion) is ExactCircleMotion:
            trace = trial.subtract_exact_full_circle(
                motion.center.x,
                motion.center.y,
                motion.phase_vector.x,
                motion.phase_vector.y,
                motion.clockwise,
                tool_radius.value,
                policy.chord_bound.value,
                policy.center_count_limit,
            )
        else:
            raise InvalidStockAreaError("stock depletion requires one exact motion type.")

        center_parameters = _validated_center_parameters(
            trace,
            motion,
            tool_radius,
            policy,
        )
        witness = DepletionWitness(
            motion,
            policy,
            tool_radius,
            center_parameters,
            trace.strategy_version,
            tuple(witness.digest for witness in self._lineage),
        )
        self._raw = trial
        self._lineage = self._lineage + (witness,)
        return witness


def _validated_center_parameters(
    trace: "_stock_2.DepletionTrace",
    motion: ExactSegmentMotion | ExactCircleMotion,
    tool_radius: ToolRadius,
    policy: DepletionPolicy,
) -> tuple[ExactCenterParameterV1, ...]:
    if not (
        trace.exact_incidence
        and trace.exact_parameters_in_range
        and trace.exact_anchors_present
        and trace.exact_removal_radius_valid
        and trace.exact_chord_bound_holds
        and (not trace.cyclic or trace.exact_seam_chord_bound_holds)
    ):
        raise InvalidDepletionTraceError("native exact-depletion trace failed an exact structural predicate.")
    try:
        parameters = tuple(
            ExactCenterParameterV1.build(
                chart=chart,
                numerator=numerator,
                denominator=denominator,
            )
            for chart, numerator, denominator in trace.center_parameters
        )
    except (TypeError, ValueError):
        raise InvalidDepletionTraceError("native exact-depletion trace contains malformed structural parameters.") from None
    if len(parameters) != trace.center_count or not parameters:
        raise InvalidDepletionTraceError("native exact-depletion trace center count is inconsistent.")
    if trace.center_count > policy.center_count_limit:
        raise InvalidDepletionTraceError("native exact-depletion trace exceeds the owned center-count policy.")
    if trace.max_chord != policy.chord_bound.value:
        raise InvalidDepletionTraceError("native exact-depletion trace changed the owned chord bound.")
    if trace.removal_radius != tool_radius.value:
        raise InvalidDepletionTraceError("native exact-depletion trace changed the owned tool radius.")
    if type(trace.strategy_version) is not bytes or not trace.strategy_version:
        raise InvalidDepletionTraceError("native exact-depletion trace has no closed strategy version.")
    if type(motion) is ExactSegmentMotion and (trace.cyclic or any(parameter.chart != -1 for parameter in parameters)):
        raise InvalidDepletionTraceError("native exact segment trace has circular structure.")
    if type(motion) is ExactCircleMotion and (not trace.cyclic or any(parameter.chart == -1 for parameter in parameters)):
        raise InvalidDepletionTraceError("native exact circle trace has segment structure.")
    return parameters
