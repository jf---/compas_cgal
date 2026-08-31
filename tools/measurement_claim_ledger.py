"""Joint generator and benchmark measurement-claim ledger consumer."""

from __future__ import annotations

import pathlib
from typing import Dict
from typing import Tuple

from tools.measurement_claim_artifact_validation import validate_benchmark_claim_artifact
from tools.measurement_claim_artifact_validation import validate_claim_artifact
from tools.measurement_claim_benchmark_history import validate_benchmark_history
from tools.measurement_claim_benchmark_schema import BenchmarkClaimPayload
from tools.measurement_claim_benchmark_validation import render_benchmark_ledger_evidence
from tools.measurement_claim_schema import GeneratorClaimPayload
from tools.measurement_claim_task6_ledger import COMPLETE_STATUS as COMPLETE_STATUS
from tools.measurement_claim_task6_ledger import FROZEN_SOURCE_COMMIT as FROZEN_SOURCE_COMMIT
from tools.measurement_claim_task6_ledger import IMMUTABLE_COLUMNS_SHA256 as IMMUTABLE_COLUMNS_SHA256
from tools.measurement_claim_task6_ledger import InvalidMeasurementClaimLedgerError as InvalidMeasurementClaimLedgerError
from tools.measurement_claim_task6_ledger import LedgerRow as LedgerRow
from tools.measurement_claim_task6_ledger import render_ledger_evidence as render_generator_ledger_evidence
from tools.measurement_claim_task6_ledger import render_ledger_evidence as render_ledger_evidence
from tools.measurement_claim_task6_ledger import validate_ledger_structure as validate_ledger_structure
from tools.measurement_claim_task6_ledger import validate_task6_source_correction as validate_task6_source_correction
from tools.measurement_claim_task6_ledger import validate_task6_source_lineage as validate_task6_source_lineage

ClaimArtifactDirectories = Tuple[pathlib.Path, pathlib.Path]


def _require_claim_union(generator_payload: GeneratorClaimPayload, benchmark_payload: BenchmarkClaimPayload) -> None:
    claims = [*generator_payload["claims"], *benchmark_payload["claims"]]
    identifiers = [claim["claim_id"] for claim in claims]
    expected = [f"MC-{index:03d}" for index in range(1, 15)]
    if identifiers != expected:
        raise InvalidMeasurementClaimLedgerError("artifact claim union must be the exact disjoint MC-001..MC-014 sequence")


def _compare_rows(
    rows: Tuple[LedgerRow, ...],
    generator_payload: GeneratorClaimPayload,
    benchmark_payload: BenchmarkClaimPayload,
    generator_evidence: Dict[str, str],
    benchmark_evidence: Dict[str, str],
) -> None:
    if len(rows) != 14:
        raise InvalidMeasurementClaimLedgerError(f"joint ledger acceptance requires exactly 14 rows, got {len(rows)}")
    _require_claim_union(generator_payload, benchmark_payload)
    claims = [*generator_payload["claims"], *benchmark_payload["claims"]]
    evidence = {**generator_evidence, **benchmark_evidence}
    for index, (row, claim) in enumerate(zip(rows, claims), start=1):
        claim_id = f"MC-{index:03d}"
        if row["ordinal"] != f"{index:03d}" or row["claim_id"] != claim_id:
            raise InvalidMeasurementClaimLedgerError(f"joint ledger row {index:03d} identity is not canonical")
        if row["disposition"] != claim["disposition"]:
            raise InvalidMeasurementClaimLedgerError(f"{claim_id} disposition differs from authenticated evidence")
        if row["evidence"] != evidence[claim_id]:
            raise InvalidMeasurementClaimLedgerError(f"{claim_id} evidence is not byte-equal to authenticated rendering")
    if any(row["disposition"] == "pending" for row in rows):
        raise InvalidMeasurementClaimLedgerError("joint ledger must contain zero pending rows")


def validate_ledger_evidence(ledger: pathlib.Path, artifact_directories: ClaimArtifactDirectories) -> None:
    """Authenticate the exact ordered generator/benchmark pair and ledger."""
    if type(artifact_directories) is not tuple or len(artifact_directories) != 2 or any(not isinstance(path, pathlib.Path) for path in artifact_directories):
        raise InvalidMeasurementClaimLedgerError("joint ledger acceptance requires an exact two-member tuple of paths")
    generator_directory, benchmark_directory = artifact_directories
    rows = validate_ledger_structure(ledger)
    generator_payload, generator_envelope, generator_started, generator_artifact = validate_claim_artifact(generator_directory)
    mc007 = next(claim for claim in generator_payload["claims"] if claim["claim_id"] == "MC-007")
    generator_repository = generator_directory.resolve().parent.parent.parent
    validate_task6_source_correction(
        generator_repository,
        str(generator_payload["source_correction_commit"]),
        mc007_disposition=mc007["disposition"],
    )
    validate_task6_source_lineage(
        generator_repository,
        str(generator_payload["source_correction_commit"]),
        str(generator_payload["source_commit"]),
        mc007_disposition=mc007["disposition"],
    )
    benchmark_payload, benchmark_envelope, benchmark_started, benchmark_artifact = validate_benchmark_claim_artifact(benchmark_directory)
    benchmark_repository = benchmark_directory.resolve().parent.parent.parent
    if generator_repository != benchmark_repository:
        raise InvalidMeasurementClaimLedgerError("generator and benchmark artifacts must belong to the same repository root")
    validate_benchmark_history(benchmark_repository, benchmark_payload)
    _require_claim_union(generator_payload, benchmark_payload)
    generator_evidence = dict(
        render_generator_ledger_evidence(
            generator_payload,
            generator_envelope,
            started=generator_started,
            artifact_directory=generator_artifact,
        )
    )
    benchmark_evidence = dict(
        render_benchmark_ledger_evidence(
            benchmark_payload,
            benchmark_envelope,
            started=benchmark_started,
            artifact_directory=benchmark_artifact,
        )
    )
    _compare_rows(
        rows,
        generator_payload,
        benchmark_payload,
        generator_evidence,
        benchmark_evidence,
    )
