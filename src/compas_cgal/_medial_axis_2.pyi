from collections.abc import Sequence
from typing import TypeAlias

import numpy as np
import numpy.typing as npt

Float64Array: TypeAlias = npt.NDArray[np.float64]
Int64Array: TypeAlias = npt.NDArray[np.int64]
MedialAxisResult: TypeAlias = tuple[
    Float64Array,
    Int64Array,
    Int64Array,
    Int64Array,
    Int64Array,
    Int64Array,
    Int64Array,
    Int64Array,
    Int64Array,
    Float64Array,
    Float64Array,
    Float64Array,
    Int64Array,
    Int64Array,
    Float64Array,
    tuple[bytes, ...],
    Int64Array,
    Int64Array,
    bytes,
    bytes,
]

class MedialAxisConstructionError(RuntimeError): ...
class InvalidReachableDomainInputError(MedialAxisConstructionError): ...
class ReachableArrangementTopologyError(MedialAxisConstructionError): ...
class PocketNotMachinableError(MedialAxisConstructionError): ...
class ReachableMaterialContainmentError(MedialAxisConstructionError): ...
class InvalidMatSamplingPolicyError(MedialAxisConstructionError): ...
class ConicSamplingLimitError(MedialAxisConstructionError): ...
class UnsupportedCanonicalMatLShapeGraphError(MedialAxisConstructionError): ...

def segment_site_medial_axis(
    vertices: Float64Array,
    holes: Sequence[Float64Array],
    tool_radius: float,
    station_spacing: float,
    max_sagitta: float,
    max_refinement_depth: int,
) -> MedialAxisResult: ...
