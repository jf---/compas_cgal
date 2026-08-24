"""Unit-bearing sealed limits for native audit decisions."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Self

from compas_cgal import _stock_2
from compas_cgal.adaptive.units import Millimetre
from compas_cgal.engagement_audit.errors import InvalidAuditDecisionDepthError
from compas_cgal.engagement_audit.errors import InvalidAuditDecisionLimitsError
from compas_cgal.engagement_audit.errors import InvalidAuditDecisionNodeLimitError
from compas_cgal.engagement_audit.errors import InvalidAuditSpatialFloorError


@dataclass(frozen=True)
class AuditSpatialFloor:
    """Positive spatial-resolution length in millimetres."""

    value: Millimetre

    def __post_init__(self) -> None:
        if isinstance(self.value, bool) or not isinstance(self.value, (int, float)):
            raise InvalidAuditSpatialFloorError("audit spatial floor must be a finite length.")
        numeric = float(self.value)
        if not math.isfinite(numeric) or numeric <= 0.0:
            raise InvalidAuditSpatialFloorError("audit spatial floor must be finite and positive.")
        object.__setattr__(self, "value", Millimetre(numeric))

    @classmethod
    def build(cls, value: float) -> Self:
        return cls(Millimetre(value))


@dataclass(frozen=True, init=False)
class AuditDecisionLimits:
    """Validated Python value paired with one opaque native limits value."""

    spatial_floor: AuditSpatialFloor
    max_depth: int
    max_nodes: int
    _native: _stock_2.AuditDecisionLimits2

    def __init__(self, *args: object, **kwargs: object) -> None:
        raise InvalidAuditDecisionLimitsError("AuditDecisionLimits must be created by AuditDecisionLimits.build().")

    @classmethod
    def build(
        cls,
        *,
        spatial_floor_mm: float,
        max_depth: int,
        max_nodes: int,
    ) -> Self:
        spatial_floor = AuditSpatialFloor.build(spatial_floor_mm)
        if type(max_depth) is not int or max_depth < 0:
            raise InvalidAuditDecisionDepthError("audit max depth must be an exact non-negative integer.")
        if type(max_nodes) is not int or max_nodes <= 0:
            raise InvalidAuditDecisionNodeLimitError("audit max nodes must be an exact positive integer.")
        try:
            native = _stock_2.build_audit_decision_limits(
                float(spatial_floor.value),
                max_depth,
                max_nodes,
            )
        except (
            _stock_2.AuditDecisionLimitsNonFiniteInputError,
            _stock_2.AuditSquaredSpatialFloorError,
        ) as error:
            raise InvalidAuditSpatialFloorError(str(error)) from error
        except _stock_2.AuditDecisionDepthLimitError as error:
            raise InvalidAuditDecisionDepthError(str(error)) from error
        except _stock_2.AuditDecisionNodeLimitError as error:
            raise InvalidAuditDecisionNodeLimitError(str(error)) from error
        instance = object.__new__(cls)
        object.__setattr__(instance, "spatial_floor", spatial_floor)
        object.__setattr__(instance, "max_depth", max_depth)
        object.__setattr__(instance, "max_nodes", max_nodes)
        object.__setattr__(instance, "_native", native)
        return instance

    @property
    def native(self) -> _stock_2.AuditDecisionLimits2:
        return self._native

    @property
    def canonical_bytes(self) -> bytes:
        return self._native.canonical_bytes
