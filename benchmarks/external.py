"""Load third-party 2D profiles as reference problems.

How often the certifier answers "could not prove it" on geometry this project did
NOT author is the number that decides whether the certifier is a product. Every
other family in this corpus is authored here, and an author who also chooses the
test cases cannot measure that.

NOTHING IS VENDORED. Profile datasets are large and separately licensed, so the
operator prepares a directory and points this loader at it. Each file is
`{"name": str, "points": [[x, y], ...]}` -- a closed ring, with or without its
first vertex repeated at the end. Converting a licensed dataset (Autodesk's
Fusion 360 Gallery sketch profiles, say) into that shape is the operator's step
and is deliberately outside this repository.

REJECTIONS ARE REPORTED, NEVER SWALLOWED. A ring the corpus cannot machine -- too
small for the tool, self-crossing, below the vertex floor -- is real dataset noise
rather than a certifier finding, so it does not become an instance. But dropping
it quietly would bias the very resolution rate this family exists to measure,
toward the easy geometry, invisibly. Every rejection is returned alongside the
accepted specs with the reason attached.

A MALFORMED FILE IS DIFFERENT and raises: it is an error in the operator's
conversion step, not a property of any geometry, and continuing past it would
measure a corpus nobody chose. So is an empty directory, which would otherwise
report a flawless run over zero instances.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from typing import List
from typing import Sequence
from typing import Tuple

from compas.geometry import Polygon
from compas.tolerance import TOL

from benchmarks.errors import BenchmarkError
from benchmarks.errors import DegeneratePocketError
from benchmarks.errors import PocketNotSimpleError
from benchmarks.spec import PocketSpec

# Family name every external instance is filed under, so a report row is
# attributable to third-party geometry at a glance.
EXTERNAL_FAMILY = "external"

# Below this many vertices a profile is a triangle or a sliver: legal geometry,
# but not a pocket anyone machines, and its presence would dilute the resolution
# rate with instances nobody cares about.
MIN_PROFILE_VERTICES = 4

# Coordinates per vertex. Profiles are planar rings; a third coordinate would mean
# the operator exported something other than a 2D sketch profile.
COORDINATES_PER_VERTEX = 2


class ExternalCorpusError(BenchmarkError):
    """An external corpus directory is missing, unreadable, or holds a malformed profile."""


class EmptyExternalCorpusError(ExternalCorpusError):
    """An external corpus directory holds no profile files at all."""


@dataclass(frozen=True)
class RejectedProfile:
    """A profile the corpus cannot machine, kept so the rejection is visible.

    Attributes:
        name: The profile's name, as the file declared it.
        reason: Why it was rejected, carrying the exception type so a reader can
            tell a sliver from a self-crossing ring without opening the file.
    """

    name: str
    reason: str


@dataclass(frozen=True)
class ExternalCorpus:
    """Everything one operator-supplied directory yielded.

    Attributes:
        directory: The directory that was read.
        specs: Instances the corpus will measure, sorted by name.
        rejected: Profiles that could not become instances, sorted by name. Never
            empty-by-omission: a caller that ignores this is choosing to.
    """

    directory: Path
    specs: Tuple[PocketSpec, ...]
    rejected: Tuple[RejectedProfile, ...]


def profile_from_points(name: str, points: Sequence[Sequence[float]], tool_diameter: float, tea_cap_deg: float) -> PocketSpec:
    """Build a spec from a raw closed ring.

    A first vertex repeated at the end is dropped: it encodes the same ring, and
    left in place it becomes a zero-length edge that the simplicity test reads as
    a crossing, throwing out a perfectly machinable profile.

    Args:
        name: Instance name.
        points: Ring vertices as ``[x, y]`` pairs, in order.
        tool_diameter: Cutter diameter.
        tea_cap_deg: Engagement cap in degrees.

    Returns:
        The spec.

    Raises:
        DegeneratePocketError: The ring has too few vertices or too little area to
            admit the tool.
        PocketNotSimpleError: The ring self-intersects.
        NonPositiveToolError: The diameter is not finite and positive.
        InvalidCapError: The cap is outside the exact kernel's contract.
    """
    ring = [(float(p[0]), float(p[1])) for p in points]
    if len(ring) > 1 and TOL.is_allclose(ring[0], ring[-1]):
        ring = ring[:-1]
    polygon = Polygon([[x, y, 0.0] for x, y in ring])
    return PocketSpec.build(
        name=name,
        family=EXTERNAL_FAMILY,
        polygon=polygon,
        tool_diameter=tool_diameter,
        tea_cap_deg=tea_cap_deg,
        params={"vertices": float(len(ring))},
    )


def load_profiles(directory: Path, tool_diameter: float, tea_cap_deg: float, min_vertices: int = MIN_PROFILE_VERTICES) -> ExternalCorpus:
    """Load every profile in *directory*, reporting the ones the corpus cannot use.

    Args:
        directory: Directory of profile JSON files.
        tool_diameter: Cutter diameter.
        tea_cap_deg: Engagement cap in degrees.
        min_vertices: Profiles with fewer vertices are rejected.

    Returns:
        The corpus, with accepted and rejected profiles both sorted by name.

    Raises:
        ExternalCorpusError: The directory is missing or a file is malformed.
        EmptyExternalCorpusError: The directory holds no profile files.
    """
    directory = Path(directory)
    if not directory.is_dir():
        raise ExternalCorpusError(f"External corpus directory does not exist: {directory}")
    paths = sorted(directory.glob("*.json"))
    if not paths:
        raise EmptyExternalCorpusError(f"External corpus directory holds no *.json profiles: {directory}")

    specs: List[PocketSpec] = []
    rejected: List[RejectedProfile] = []
    for path in paths:
        name, points = _read_profile(path)
        if len(points) < min_vertices:
            rejected.append(RejectedProfile(name=name, reason=f"{len(points)} vertices, below the floor of {min_vertices}"))
            continue
        try:
            specs.append(profile_from_points(name, points, tool_diameter, tea_cap_deg))
        except (DegeneratePocketError, PocketNotSimpleError) as exc:
            # Dataset noise, not a certifier finding -- but recorded, because a
            # silently thinner corpus measures a resolution rate nobody chose.
            rejected.append(RejectedProfile(name=name, reason=f"{type(exc).__name__}: {exc}"))
    specs.sort(key=lambda s: s.name)
    rejected.sort(key=lambda r: r.name)
    return ExternalCorpus(directory=directory, specs=tuple(specs), rejected=tuple(rejected))


def _read_profile(path: Path) -> Tuple[str, List[List[float]]]:
    """Parse one profile file, failing loudly on anything the schema does not allow.

    Args:
        path: The file to read.

    Returns:
        ``(name, points)``.

    Raises:
        ExternalCorpusError: The file is unreadable, is not JSON, lacks a required
            key, or carries vertices that are not numeric ``[x, y]`` pairs.
    """
    try:
        payload: Any = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ExternalCorpusError(f"Malformed profile {path}: {type(exc).__name__}: {exc}") from exc
    if not isinstance(payload, dict) or "name" not in payload or "points" not in payload:
        raise ExternalCorpusError(f"Malformed profile {path}: expected an object with 'name' and 'points'.")
    raw_points = payload["points"]
    if not isinstance(raw_points, list):
        raise ExternalCorpusError(f"Malformed profile {path}: 'points' is {type(raw_points).__name__!r}, expected a list of [x, y] pairs.")
    points: List[List[float]] = []
    for index, vertex in enumerate(raw_points):
        if not isinstance(vertex, (list, tuple)) or len(vertex) != COORDINATES_PER_VERTEX:
            raise ExternalCorpusError(f"Malformed profile {path}: vertex {index} is {vertex!r}, expected exactly {COORDINATES_PER_VERTEX} coordinates.")
        try:
            points.append([float(vertex[0]), float(vertex[1])])
        except (TypeError, ValueError) as exc:
            raise ExternalCorpusError(f"Malformed profile {path}: vertex {index} is {vertex!r}, whose coordinates are not numbers.") from exc
    return str(payload["name"]), points
