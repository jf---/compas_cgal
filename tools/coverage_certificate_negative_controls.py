"""Frozen negative-control inputs of the coverage certificate, and their witness digest.

Two uncovered runs gated the certificate's performance levers (a lever was kept
only if it left the witness set unchanged):

- ``figure5_omega_drop_600_899``: the Figure 5 Ω* path (exact target T(Ω*),
  exact ∂Ω* connectors) with circles 600–899 dropped;
- ``polygonal_monstera``: the Monstera path on its polygonal target (polygon
  connectors).

The generated paths drift with generator work, so their inputs are frozen here
as tracked gzip JSON (``tests/fixtures/coverage_certificate/negative_controls/``):
circles and polygon vertices as hex binary64, connector sites and one-root
connector ends as exact decimal fractions. Figure 5's exact target is rebuilt
from the reference case (it is deterministic). ``freeze`` writes the fixtures
from captured ``run_certificate`` inputs (a pickle of ``{"args", "kwargs"}``);
``load`` rebuilds the call; ``witness_digest`` is the digest the lever A/B used:
sha256 over the sorted (kind, region, x, y) of every witness, first 16 hex.

Usage: ``python -m tools.coverage_certificate_negative_controls freeze NAME CAPTURE.pkl``.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import pickle
from collections.abc import Sequence
from fractions import Fraction
from pathlib import Path
from typing import Any

from benchmarks.held_coverage_target import exact_coverage_target
from benchmarks.held_declared_boundary_emission import declared_boundary_from_reference
from benchmarks.held_exact_centre_domain import CentreDomain as ExactCentreDomain
from benchmarks.held_reference_cases import load_held_reference_case
from compas_cgal.coverage_certificate.exact_connectors import ExactConnectorPiece
from compas_cgal.coverage_certificate.exact_connectors import ExactConnectorSweep
from compas_cgal.coverage_certificate.exact_target import ExactCentreDomainTarget
from compas_cgal.coverage_certificate.exact_target import OneRootPoint
from compas_cgal.coverage_certificate.native_certificate import CertificateReport
from compas_cgal.coverage_certificate.sweeps import BoundaryConnectorSweep
from compas_cgal.coverage_certificate.sweeps import FullCircleSweep
from compas_cgal.coverage_certificate.target import CentreBoundarySite
from compas_cgal.coverage_certificate.target import CentreDomain
from compas_cgal.coverage_certificate.target import OneRootMillimetre
from compas_cgal.units import GuideRadius
from compas_cgal.units import Point2
from compas_cgal.units import ToolRadius
from compas_cgal.units import WorldXY

REPOSITORY = Path(__file__).resolve().parents[1]
FIXTURES = REPOSITORY / "tests" / "fixtures" / "coverage_certificate" / "negative_controls"
SCHEMA = "coverage-certificate-negative-control/1"
# Reference case whose declared boundary defines the exact target, per control.
EXACT_TARGET_CASES = {"figure5_omega_drop_600_899": "figure5"}
# Circles dropped from the captured Figure 5 path (index range, half-open).
FIGURE5_DROP = (600, 900)
# Hex digits of the witness digest.
DIGEST_HEX = 16


class NegativeControlFixtureError(Exception):
    """A frozen fixture does not rebuild the captured inputs exactly, or names an unknown control."""


def _fraction(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def _parse(text: str) -> Fraction:
    numerator, denominator = text.split("/")
    return Fraction(int(numerator), int(denominator))


def _point_json(point: OneRootPoint) -> list[str]:
    return [_fraction(point.x.rational), _fraction(point.x.coefficient), _fraction(point.y.rational), _fraction(point.y.coefficient), _fraction(point.radicand)]


def _point(row: Sequence[str]) -> OneRootPoint:
    xr, xc, yr, yc, radicand = (_parse(value) for value in row)
    return OneRootPoint.build(OneRootMillimetre(xr, xc), OneRootMillimetre(yr, yc), radicand)


def _site_json(site: CentreBoundarySite) -> list[object]:
    return [site.loop, site.segment, _fraction(site.parameter)]


def _site(row: Sequence[Any]) -> CentreBoundarySite:
    return CentreBoundarySite.build(int(row[0]), int(row[1]), _parse(row[2]))


def _encode(name: str, args: Sequence[Any], kwargs: dict[str, Any]) -> dict[str, Any]:
    domain, tool, circles = args[0], args[1], args[2]
    connectors = args[3] if len(args) > 3 else kwargs.get("connectors", ())
    if domain.holes:
        raise NegativeControlFixtureError("a frozen negative control carries no domain holes")
    exact = kwargs.get("exact_connectors", ())
    return {
        "schema": SCHEMA,
        "name": name,
        "tool_radius": float(tool.value).hex(),
        "domain": [[float(p.x).hex(), float(p.y).hex()] for p in domain.outer],
        "circles": [[float(c.centre.x).hex(), float(c.centre.y).hex(), float(c.path_radius.value).hex()] for c in circles],
        "connectors": [[_site_json(c.start), _site_json(c.end)] for c in connectors],
        "exact_connectors": [[[piece.piece, _point_json(piece.start), _point_json(piece.end), piece.full] for piece in sweep.pieces] for sweep in exact],
    }


def _exact_target(name: str) -> ExactCentreDomainTarget | None:
    if name not in EXACT_TARGET_CASES:
        return None
    case = load_held_reference_case(EXACT_TARGET_CASES[name])
    declared = declared_boundary_from_reference(case.name, case.boundary)
    return exact_coverage_target(declared, ExactCentreDomain.build(declared, Fraction(float(case.tool_radius.value))))


def load(name: str) -> tuple[tuple[Any, ...], dict[str, Any]]:
    """The frozen ``run_certificate`` call ``(args, kwargs)`` of control ``name``.

    Raises:
        NegativeControlFixtureError: An unknown control or a fixture of another schema.
    """
    path = FIXTURES / f"{name}.json.gz"
    if not path.exists():
        raise NegativeControlFixtureError(f"no frozen negative control {name!r} in {FIXTURES}")
    data = json.loads(gzip.decompress(path.read_bytes()))
    if data.get("schema") != SCHEMA or data.get("name") != name:
        raise NegativeControlFixtureError(f"{path} is not a {SCHEMA} fixture of {name!r}")
    domain = CentreDomain.build([Point2[WorldXY].build(float.fromhex(x), float.fromhex(y)) for x, y in data["domain"]])
    tool = ToolRadius.build(float.fromhex(data["tool_radius"]))
    circles = [FullCircleSweep.build(Point2[WorldXY].build(float.fromhex(x), float.fromhex(y)), GuideRadius.build(float.fromhex(r))) for x, y, r in data["circles"]]
    connectors = [BoundaryConnectorSweep.build(domain, _site(start), _site(end)) for start, end in data["connectors"]]
    kwargs: dict[str, Any] = {}
    target = _exact_target(name)
    if target is not None:
        kwargs["exact_target"] = target
        kwargs["exact_connectors"] = [
            ExactConnectorSweep.build([ExactConnectorPiece.build(target, int(index), _point(start), _point(end), full=bool(full)) for index, start, end, full in sweep])
            for sweep in data["exact_connectors"]
        ]
    args: tuple[Any, ...] = (domain, tool, circles, connectors) if connectors else (domain, tool, circles)
    return args, kwargs


def witness_digest(report: CertificateReport) -> str:
    """sha256 over the sorted (kind, region, x, y) of every witness, first ``DIGEST_HEX`` hex digits."""
    witnesses = sorted((e.piece_kind.name, e.region_index, str(e.uncovered[0]), str(e.uncovered[1])) for e in report.exposed)
    return hashlib.sha256(repr(witnesses).encode()).hexdigest()[:DIGEST_HEX]


def freeze(name: str, capture: Path) -> Path:
    """Write the fixture of ``name`` from a captured call and check that :func:`load` rebuilds it exactly.

    Raises:
        NegativeControlFixtureError: The rebuilt call differs from the capture.
    """
    captured = pickle.loads(capture.read_bytes())
    args, kwargs = list(captured["args"]), dict(captured["kwargs"])
    if name == "figure5_omega_drop_600_899":
        low, high = FIGURE5_DROP
        args[2] = [c for i, c in enumerate(args[2]) if not low <= i < high]
    FIXTURES.mkdir(parents=True, exist_ok=True)
    path = FIXTURES / f"{name}.json.gz"
    path.write_bytes(gzip.compress(json.dumps(_encode(name, args, kwargs), separators=(",", ":")).encode(), mtime=0))
    rebuilt_args, rebuilt_kwargs = load(name)
    if repr(tuple(rebuilt_args)) != repr(tuple(a for a in args if not (isinstance(a, list) and not a))) or repr(rebuilt_kwargs) != repr(kwargs):
        raise NegativeControlFixtureError(f"{path} does not rebuild the captured call exactly")
    return path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)
    frozen = sub.add_parser("freeze")
    frozen.add_argument("name")
    frozen.add_argument("capture", type=Path)
    args = parser.parse_args()
    print(freeze(args.name, args.capture))


if __name__ == "__main__":
    main()
