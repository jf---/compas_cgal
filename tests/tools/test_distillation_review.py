from pathlib import Path

import pytest

from tools.distillation_review import DistillationArtifactError
from tools.distillation_review import validate_review


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def _manifest(
    tmp_path: Path,
    *,
    duplicate_variant: bool = False,
    pass_1: str = "complete",
    pass_2: str = "complete",
    pass_3: str = "complete",
) -> Path:
    header = "VARIANT\tREF\tWORKTREE\tPATH\tFAMILY\tSCOPE\tSAME_AS\tPASS_1\tPASS_2\tPASS_3\tDRIFT\tNOTE\n"
    row = f"frontier:src/stock_2.cpp\tfrontier\t/worktree/frontier\tsrc/stock_2.cpp\tstock\tincluded\t\t{pass_1}\t{pass_2}\t{pass_3}\tclean\tCAP-0001\n"
    return _write(
        tmp_path / "manifest.tsv",
        header + row + (row if duplicate_variant else ""),
    )


def _complete_manifest(tmp_path: Path) -> Path:
    return _manifest(tmp_path)


def _capabilities(tmp_path: Path) -> Path:
    return _write(tmp_path / "capabilities.md", "## CAP-0001 — Exact stock\n")


def _findings(
    tmp_path: Path,
    *,
    capability: str = "",
    status: str = "proposed",
) -> Path:
    header = "ID\tPASS\tFAMILY\tVARIANT\tLOCATION\tQUOTE\tCAPABILITY\tISSUE\tCONSUMER\tLOSS_IF_CHANGED\tPROPOSED_CONDENSATION\tORACLE\tCOST\tFALSIFIER\tEVIDENCE\tSTATUS\n"
    if not capability:
        return _write(tmp_path / "findings.tsv", header)
    row = f"F-0001\t1\tstock\tfrontier:src/stock_2.cpp\tstock_2.cpp:1\tquote\t{capability}\tissue\tconsumer\tloss\tproposal\toracle\tcheap\tfalsifier\tevidence\t{status}\n"
    return _write(tmp_path / "findings.tsv", header + row)


def _surgery(
    tmp_path: Path,
    *,
    capability: str = "CAP-0001",
    disposition: str = "KEEP",
    retained_owner: str = "Exact stock",
    valuable_nucleus: str = "exact predicate",
    consumer_contracts: str = "Stock consumer",
    loss_argument: str = "loss argument",
    oracle: str = "focused test",
    falsifier: str = "A smaller implementation proves consumer equivalence",
    residual_risk: str = "none",
) -> Path:
    row = f"| {capability} | {disposition} | {retained_owner} | {valuable_nucleus} | {consumer_contracts} | {loss_argument} | {oracle} | {falsifier} | {residual_risk} |\n"
    return _write(
        tmp_path / "surgery.md",
        "| Capability | Disposition | Retained owner | Valuable nucleus | "
        "Consumer contracts | Loss argument | Oracle | Falsifier | Residual risk |\n"
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- |\n" + row,
    )


def test_manifest_requires_unique_variant_and_three_pass_states(tmp_path: Path) -> None:
    manifest = _manifest(tmp_path, duplicate_variant=True)
    with pytest.raises(DistillationArtifactError, match="duplicate variant"):
        validate_review(manifest, _capabilities(tmp_path), _findings(tmp_path), None)


@pytest.mark.parametrize("state", ["not-started", "complete", "stale"])
def test_manifest_accepts_allowed_pass_states(tmp_path: Path, state: str) -> None:
    validate_review(
        _manifest(tmp_path, pass_1=state),
        _capabilities(tmp_path),
        _findings(tmp_path),
        None,
    )


def test_manifest_rejects_unknown_pass_state(tmp_path: Path) -> None:
    with pytest.raises(DistillationArtifactError, match="PASS_2.*in-progress"):
        validate_review(
            _manifest(tmp_path, pass_2="in-progress"),
            _capabilities(tmp_path),
            _findings(tmp_path),
            None,
        )


@pytest.mark.parametrize("status", ["proposed", "verified", "rejected", "queued-c2", "jelle-c3"])
def test_findings_accept_spec_statuses(tmp_path: Path, status: str) -> None:
    validate_review(
        _complete_manifest(tmp_path),
        _capabilities(tmp_path),
        _findings(tmp_path, capability="CAP-0001", status=status),
        None,
    )


def test_findings_reject_unknown_status(tmp_path: Path) -> None:
    with pytest.raises(DistillationArtifactError, match="STATUS.*accepted"):
        validate_review(
            _complete_manifest(tmp_path),
            _capabilities(tmp_path),
            _findings(tmp_path, capability="CAP-0001", status="accepted"),
            None,
        )


def test_findings_require_existing_capability(tmp_path: Path) -> None:
    with pytest.raises(DistillationArtifactError, match="unknown capability CAP-9999"):
        validate_review(
            _complete_manifest(tmp_path),
            _capabilities(tmp_path),
            _findings(tmp_path, capability="CAP-9999"),
            None,
        )


def test_complete_surgery_requires_three_reads_per_scoped_variant(tmp_path: Path) -> None:
    manifest = _manifest(tmp_path, pass_3="not-started")
    surgery = _surgery(tmp_path)
    with pytest.raises(DistillationArtifactError, match="three complete readings"):
        validate_review(manifest, _capabilities(tmp_path), _findings(tmp_path), surgery)


def test_surgery_requires_existing_capability(tmp_path: Path) -> None:
    with pytest.raises(DistillationArtifactError, match="unknown capability CAP-9999"):
        validate_review(
            _complete_manifest(tmp_path),
            _capabilities(tmp_path),
            _findings(tmp_path),
            _surgery(tmp_path, capability="CAP-9999"),
        )


def test_condense_requires_nucleus_and_consumer_contracts(tmp_path: Path) -> None:
    surgery = _surgery(tmp_path, disposition="CONDENSE", valuable_nucleus="")
    with pytest.raises(DistillationArtifactError, match="CONDENSE.*valuable nucleus"):
        validate_review(_complete_manifest(tmp_path), _capabilities(tmp_path), _findings(tmp_path), surgery)


@pytest.mark.parametrize("disposition", ["ABSORB", "REMOVE"])
def test_reductions_require_loss_oracle_falsifier_and_retained_owner(tmp_path: Path, disposition: str) -> None:
    surgery = _surgery(tmp_path, disposition=disposition, falsifier="")
    with pytest.raises(DistillationArtifactError, match=rf"{disposition}.*falsifier"):
        validate_review(_complete_manifest(tmp_path), _capabilities(tmp_path), _findings(tmp_path), surgery)


def test_unknown_accepts_named_missing_evidence_and_falsifier(tmp_path: Path) -> None:
    validate_review(
        _complete_manifest(tmp_path),
        _capabilities(tmp_path),
        _findings(tmp_path),
        _surgery(
            tmp_path,
            disposition="UNKNOWN",
            falsifier="Run controller backplot",
            residual_risk="MISSING_EVIDENCE: controller backplot",
        ),
    )


def test_unknown_rejects_unnamed_missing_evidence(tmp_path: Path) -> None:
    with pytest.raises(DistillationArtifactError, match="UNKNOWN.*MISSING_EVIDENCE"):
        validate_review(
            _complete_manifest(tmp_path),
            _capabilities(tmp_path),
            _findings(tmp_path),
            _surgery(tmp_path, disposition="UNKNOWN", residual_risk="unknown"),
        )
