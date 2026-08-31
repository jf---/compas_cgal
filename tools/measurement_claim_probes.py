"""Run and authenticate the fixed generator measurement-claim batch."""

from __future__ import annotations

import argparse
import datetime
import json
import pathlib
import sys
from typing import Dict
from typing import List
from typing import Optional
from typing import Sequence
from typing import Tuple
from typing import cast

from tools.measurement_artifact import ArtifactKind
from tools.measurement_artifact import IdentityVersion
from tools.measurement_artifact import MeasurementArtifactCollisionError
from tools.measurement_artifact import build_envelope
from tools.measurement_artifact import capture_clean_source
from tools.measurement_artifact import input_identity_sha256
from tools.measurement_artifact import publication_stage
from tools.measurement_artifact import publish_stage
from tools.measurement_artifact import write_envelope
from tools.measurement_claim_advance import run_advance_case
from tools.measurement_claim_artifact_validation import validate_benchmark_claim_artifact
from tools.measurement_claim_artifact_validation import validate_claim_artifact
from tools.measurement_claim_benchmark_identity import ARTIFACT_KIND as BENCHMARK_ARTIFACT_KIND
from tools.measurement_claim_benchmark_identity import INPUT_VERSION as BENCHMARK_INPUT_VERSION
from tools.measurement_claim_benchmark_identity import MARKDOWN_NAME
from tools.measurement_claim_benchmark_identity import PAYLOAD_NAME as BENCHMARK_PAYLOAD_NAME
from tools.measurement_claim_benchmark_identity import RAW_JSON_NAME
from tools.measurement_claim_benchmark_identity import RESULT_VERSION as BENCHMARK_RESULT_VERSION
from tools.measurement_claim_benchmark_semantic_input import benchmark_payload_semantic_input
from tools.measurement_claim_benchmark_semantic_input import benchmark_semantic_input
from tools.measurement_claim_errors import InvalidMeasurementClaimConfigError
from tools.measurement_claim_errors import InvalidMeasurementClaimPayloadError
from tools.measurement_claim_figure6 import run_figure6_case
from tools.measurement_claim_identity import ARTIFACT_KIND
from tools.measurement_claim_identity import INPUT_VERSION
from tools.measurement_claim_identity import PAYLOAD_NAME
from tools.measurement_claim_identity import RESULT_VERSION
from tools.measurement_claim_ledger import validate_ledger_evidence
from tools.measurement_claim_payload import compose_generator_payload
from tools.measurement_claim_radial import run_radial_case
from tools.measurement_claim_schema import GeneratorCase
from tools.measurement_claim_schema import GeneratorCasePayload
from tools.measurement_claim_semantic_input import generator_semantic_input
from tools.measurement_claim_task6_ledger import render_ledger_evidence

GENERATOR_CASE_ORDER: Tuple[GeneratorCase, ...] = (
    "radial-station",
    "radial-subdivisions",
    "radial-floor",
    "radial-margin",
    "advance-placement",
    "advance-probe-count",
)
BENCHMARK_CASE = "figure6-spacing"


def _selected_cases(selected: Sequence[str]) -> Tuple[GeneratorCase, ...]:
    if not selected:
        raise InvalidMeasurementClaimConfigError("at least one generator measurement-claim case is required")
    if len(set(selected)) != len(selected):
        raise InvalidMeasurementClaimConfigError("generator measurement-claim cases must be unique")
    unknown = tuple(case for case in selected if case not in GENERATOR_CASE_ORDER)
    if unknown:
        raise InvalidMeasurementClaimConfigError(f"unknown generator measurement-claim cases: {unknown!r}")
    ordered = tuple(case for case in GENERATOR_CASE_ORDER if case in selected)
    if tuple(selected) != ordered:
        raise InvalidMeasurementClaimConfigError("generator measurement-claim cases must follow the fixed order")
    return ordered


def _run_cases(selected: Sequence[str]) -> List[GeneratorCasePayload]:
    cases = _selected_cases(selected)
    results: List[GeneratorCasePayload] = []
    for case in cases:
        if case == "radial-station":
            result = run_radial_case("radial-station")
        elif case == "radial-subdivisions":
            result = run_radial_case("radial-subdivisions")
        elif case == "radial-floor":
            result = run_radial_case("radial-floor")
        elif case == "radial-margin":
            result = run_radial_case("radial-margin")
        elif case == "advance-placement":
            result = run_advance_case("advance-placement")
        else:
            result = run_advance_case("advance-probe-count")
        results.append(result)
    return results


def _now_utc() -> datetime.datetime:
    return datetime.datetime.now(datetime.timezone.utc)


def _payload_bytes(payload: object) -> bytes:
    return (json.dumps(payload, indent=2, allow_nan=False) + "\n").encode("utf-8")


def _input_sha256(envelope: Dict[str, object]) -> str:
    identity = cast(Dict[str, object], envelope["input_identity"])
    digest = identity["sha256"]
    if type(digest) is not str or len(digest) != 64:
        raise InvalidMeasurementClaimPayloadError("common envelope returned an invalid input identity")
    return digest


def _produce_generator(repository: pathlib.Path, results_root: pathlib.Path) -> pathlib.Path:
    source = capture_clean_source(repository)
    started = _now_utc()
    cases = _run_cases(GENERATOR_CASE_ORDER)
    payload = compose_generator_payload(source.commit, cases)
    payload_bytes = _payload_bytes(payload)
    finished = _now_utc()
    command = [sys.executable, "-m", "tools.measurement_claim_probes", "run-generator", "--all"]
    envelope = build_envelope(
        artifact_kind=ArtifactKind(ARTIFACT_KIND),
        source=source,
        started=started,
        finished=finished,
        argv=command,
        input_version=IdentityVersion(INPUT_VERSION),
        input_payload=generator_semantic_input(payload),
        result_version=IdentityVersion(RESULT_VERSION),
        payloads={PAYLOAD_NAME: payload_bytes},
    )
    logical_name = f"{started.date().isoformat()}-{str(source.commit)[:12]}-generator-{_input_sha256(envelope)[:12]}"
    final = results_root / logical_name
    if final.exists() or final.is_symlink():
        raise MeasurementArtifactCollisionError(f"measurement result already exists: {final}")
    with publication_stage(results_root, logical_name) as stage:
        (stage / PAYLOAD_NAME).write_bytes(payload_bytes)
        write_envelope(stage, envelope)
        staged = validate_claim_artifact(stage)
        published = publish_stage(source=source, stage=stage, final=final)
        validated = validate_claim_artifact(published)
    if staged[3] != validated[3]:
        raise InvalidMeasurementClaimPayloadError("stage and final validation produced different canonical artifact paths")
    final_payload, final_envelope, final_started, final_directory = validated
    render_ledger_evidence(
        final_payload,
        final_envelope,
        started=final_started,
        artifact_directory=final_directory,
    )
    return published


def _produce_benchmark(repository: pathlib.Path, results_root: pathlib.Path) -> pathlib.Path:
    source = capture_clean_source(repository)
    started = _now_utc()
    semantic_input = benchmark_semantic_input(source.commit)
    input_digest = input_identity_sha256(version=IdentityVersion(BENCHMARK_INPUT_VERSION), payload=semantic_input)
    logical_name = f"{started.date().isoformat()}-{str(source.commit)[:12]}-benchmark-{str(input_digest)[:12]}"
    final = results_root / logical_name
    if final.exists() or final.is_symlink():
        raise MeasurementArtifactCollisionError(f"measurement result already exists: {final}")
    with publication_stage(results_root, logical_name) as stage:
        payload, executed_argv = run_figure6_case(stage, source_commit=source.commit)
        if benchmark_payload_semantic_input(payload) != semantic_input:
            raise InvalidMeasurementClaimPayloadError("Figure-6 result semantic input differs from the precomputed input")
        payload_bytes = _payload_bytes(payload)
        (stage / BENCHMARK_PAYLOAD_NAME).write_bytes(payload_bytes)
        raw_payloads = {
            MARKDOWN_NAME: (stage / MARKDOWN_NAME).read_bytes(),
            RAW_JSON_NAME: (stage / RAW_JSON_NAME).read_bytes(),
            BENCHMARK_PAYLOAD_NAME: payload_bytes,
        }
        envelope = build_envelope(
            artifact_kind=ArtifactKind(BENCHMARK_ARTIFACT_KIND),
            source=source,
            started=started,
            finished=_now_utc(),
            argv=executed_argv,
            input_version=IdentityVersion(BENCHMARK_INPUT_VERSION),
            input_payload=semantic_input,
            result_version=IdentityVersion(BENCHMARK_RESULT_VERSION),
            payloads=raw_payloads,
        )
        if _input_sha256(envelope) != input_digest:
            raise InvalidMeasurementClaimPayloadError("Figure-6 envelope input digest differs from its precomputed digest")
        write_envelope(stage, envelope)
        staged = validate_benchmark_claim_artifact(stage)
        published = publish_stage(source=source, stage=stage, final=final)
        validated = validate_benchmark_claim_artifact(published)
    if staged[3] != validated[3]:
        raise InvalidMeasurementClaimPayloadError("benchmark stage and final validation produced different canonical artifact paths")
    return published


def _run_benchmark_case(case: str, repository: pathlib.Path, results_root: pathlib.Path) -> pathlib.Path:
    if case != BENCHMARK_CASE:
        raise InvalidMeasurementClaimConfigError(f"unknown benchmark measurement-claim case: {case!r}")
    return _produce_benchmark(repository, results_root)


def main(argv: Optional[Sequence[str]] = None) -> int:
    """List or run the fixed generator measurement-claim cases."""
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("list", help="list generator claim cases")
    run_parser = subparsers.add_parser("run-generator", help="run the fixed generator claim batch")
    selection = run_parser.add_mutually_exclusive_group(required=True)
    selection.add_argument("--all", action="store_true", help="run the canonical authenticated batch")
    selection.add_argument("--case", action="append", default=[], help="run one or more cases in canonical order as diagnostics")
    benchmark_parser = subparsers.add_parser("run-benchmark", help="run the fixed Figure-6 benchmark claim batch")
    benchmark_parser.add_argument("--case", choices=(BENCHMARK_CASE,), required=True)
    validate_parser = subparsers.add_parser("validate", help="validate one generator claim artifact")
    validate_parser.add_argument("artifact", type=pathlib.Path)
    ledger_parser = subparsers.add_parser("validate-ledger", help="validate the ledger against generator and benchmark artifacts")
    ledger_parser.add_argument("--ledger", type=pathlib.Path, default=pathlib.Path("docs/measurement_claims.md"))
    ledger_parser.add_argument("artifacts", type=pathlib.Path, nargs=2)
    arguments = parser.parse_args(argv)
    if arguments.command == "list":
        for case in GENERATOR_CASE_ORDER:
            print(case)
        return 0
    if arguments.command == "run-generator":
        if arguments.all:
            repository = pathlib.Path.cwd().resolve()
            result = _produce_generator(repository, repository / "benchmarks" / "measurement_claim_results")
            print(result.relative_to(repository).as_posix())
            return 0
        selected = _selected_cases(arguments.case)
        print(json.dumps(_run_cases(selected), indent=2, allow_nan=False))
        return 0
    if arguments.command == "validate":
        _, _, _, canonical = validate_claim_artifact(arguments.artifact)
        print(canonical)
        return 0
    if arguments.command == "run-benchmark":
        repository = pathlib.Path.cwd().resolve()
        result = _run_benchmark_case(arguments.case, repository, repository / "benchmarks" / "measurement_claim_results")
        print(result.relative_to(repository).as_posix())
        return 0
    if arguments.command == "validate-ledger":
        artifacts = cast(Tuple[pathlib.Path, pathlib.Path], tuple(arguments.artifacts))
        validate_ledger_evidence(arguments.ledger, artifacts)
        print(arguments.ledger)
        return 0
    raise InvalidMeasurementClaimConfigError(f"unsupported measurement-claim command: {arguments.command!r}")


if __name__ == "__main__":
    raise SystemExit(main())
