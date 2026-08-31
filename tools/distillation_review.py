from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


class DistillationArtifactError(Exception):
    """Raised when a distillation review artifact violates its structure."""


MANIFEST_COLUMNS = (
    "VARIANT",
    "REF",
    "WORKTREE",
    "PATH",
    "FAMILY",
    "SCOPE",
    "SAME_AS",
    "PASS_1",
    "PASS_2",
    "PASS_3",
    "DRIFT",
    "NOTE",
)
FINDING_COLUMNS = (
    "ID",
    "PASS",
    "FAMILY",
    "VARIANT",
    "LOCATION",
    "QUOTE",
    "CAPABILITY",
    "ISSUE",
    "CONSUMER",
    "LOSS_IF_CHANGED",
    "PROPOSED_CONDENSATION",
    "ORACLE",
    "COST",
    "FALSIFIER",
    "EVIDENCE",
    "STATUS",
)
SURGERY_COLUMNS = (
    "CAPABILITY",
    "DISPOSITION",
    "RETAINED_OWNER",
    "VALUABLE_NUCLEUS",
    "CONSUMER_CONTRACTS",
    "LOSS_ARGUMENT",
    "ORACLE",
    "FALSIFIER",
    "RESIDUAL_RISK",
)
PASS_STATES = frozenset({"not-started", "complete", "stale"})
SCOPES = frozenset({"included", "excluded"})
DRIFT_STATES = frozenset({"clean", "stale"})
FINDING_STATUSES = frozenset({"proposed", "verified", "rejected", "queued-c2", "jelle-c3"})
DISPOSITIONS = frozenset({"KEEP", "CONDENSE", "ABSORB", "QUARANTINE", "REMOVE", "UNKNOWN"})
CAPABILITY_HEADING = re.compile(r"^##\s+(CAP-\d{4})\s+—\s+\S", re.MULTILINE)

STATE_DIR = Path("docs/superpowers/state")
DEFAULT_MANIFEST = STATE_DIR / "2026-08-31-distillation-manifest.tsv"
DEFAULT_CAPABILITIES = STATE_DIR / "2026-08-31-distillation-capabilities.md"
DEFAULT_FINDINGS = STATE_DIR / "2026-08-31-distillation-findings.tsv"
DEFAULT_SURGERY = STATE_DIR / "2026-08-31-distillation-surgery.md"


def _read_tsv(path: Path, required: tuple[str, ...]) -> list[dict[str, str]]:
    try:
        stream = path.open(encoding="utf-8", newline="")
    except OSError as error:
        raise DistillationArtifactError(f"cannot read {path}: {error}") from error

    with stream:
        reader = csv.DictReader(stream, delimiter="\t")
        fields = tuple(reader.fieldnames or ())
        missing = [column for column in required if column not in fields]
        if missing:
            raise DistillationArtifactError(f"{path} missing required columns: {', '.join(missing)}")
        rows: list[dict[str, str]] = []
        for line, raw_row in enumerate(reader, start=2):
            if None in raw_row:
                raise DistillationArtifactError(f"{path}:{line} has extra TSV fields")
            row: dict[str, str] = {}
            for column in required:
                value = raw_row[column]
                if value is None:
                    raise DistillationArtifactError(f"{path}:{line} has no value for {column}")
                row[column] = value.strip()
            rows.append(row)
        return rows


def _require_allowed(
    path: Path,
    rows: list[dict[str, str]],
    column: str,
    allowed: frozenset[str],
) -> None:
    for line, row in enumerate(rows, start=2):
        value = row[column]
        if value not in allowed:
            raise DistillationArtifactError(f"{path}:{line} {column} has disallowed value {value!r}")


def _read_capabilities(path: Path) -> set[str]:
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as error:
        raise DistillationArtifactError(f"cannot read {path}: {error}") from error
    capabilities = set(CAPABILITY_HEADING.findall(text))
    if not capabilities:
        raise DistillationArtifactError(f"{path} has no capability headings")
    return capabilities


def _markdown_cells(line: str) -> list[str]:
    return [cell.strip() for cell in line.strip().strip("|").split("|")]


def _column_name(cell: str) -> str:
    return re.sub(r"[^A-Z0-9]+", "_", cell.upper()).strip("_")


def _read_surgery(path: Path) -> list[dict[str, str]]:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as error:
        raise DistillationArtifactError(f"cannot read {path}: {error}") from error

    for index, line in enumerate(lines):
        columns = tuple(_column_name(cell) for cell in _markdown_cells(line))
        if columns != SURGERY_COLUMNS:
            continue
        if index + 1 >= len(lines):
            break
        separator = _markdown_cells(lines[index + 1])
        if len(separator) != len(columns) or any(re.fullmatch(r":?-{3,}:?", cell) is None for cell in separator):
            raise DistillationArtifactError(f"{path} has malformed surgery table separator")
        rows: list[dict[str, str]] = []
        for line_number, row_line in enumerate(lines[index + 2 :], start=index + 3):
            if not row_line.lstrip().startswith("|"):
                break
            cells = _markdown_cells(row_line)
            if len(cells) != len(columns):
                raise DistillationArtifactError(f"{path}:{line_number} has {len(cells)} surgery fields; expected {len(columns)}")
            rows.append(dict(zip(columns, cells, strict=True)))
        return rows
    raise DistillationArtifactError(f"{path} missing required surgery table columns")


def _require_capability(path: Path, line: int, capability: str, capabilities: set[str]) -> None:
    if capability not in capabilities:
        raise DistillationArtifactError(f"{path}:{line} references unknown capability {capability or '<empty>'}")


def _validate_manifest(path: Path) -> list[dict[str, str]]:
    rows = _read_tsv(path, MANIFEST_COLUMNS)
    _require_allowed(path, rows, "SCOPE", SCOPES)
    _require_allowed(path, rows, "DRIFT", DRIFT_STATES)
    for pass_column in ("PASS_1", "PASS_2", "PASS_3"):
        _require_allowed(path, rows, pass_column, PASS_STATES)

    seen: set[str] = set()
    for line, row in enumerate(rows, start=2):
        variant = row["VARIANT"]
        if variant in seen:
            raise DistillationArtifactError(f"{path}:{line} duplicate variant {variant}")
        seen.add(variant)
    return rows


def _validate_findings(path: Path, capabilities: set[str]) -> list[dict[str, str]]:
    rows = _read_tsv(path, FINDING_COLUMNS)
    _require_allowed(path, rows, "STATUS", FINDING_STATUSES)
    for line, row in enumerate(rows, start=2):
        _require_capability(path, line, row["CAPABILITY"], capabilities)
    return rows


def _require_fields(
    path: Path,
    line: int,
    row: dict[str, str],
    disposition: str,
    fields: tuple[str, ...],
) -> None:
    for field in fields:
        if not row[field]:
            label = field.lower().replace("_", " ")
            raise DistillationArtifactError(f"{path}:{line} {disposition} requires {label}")


def _validate_surgery(path: Path, rows: list[dict[str, str]], capabilities: set[str]) -> None:
    for line, row in enumerate(rows, start=3):
        capability = row["CAPABILITY"]
        _require_capability(path, line, capability, capabilities)
        disposition = row["DISPOSITION"]
        if disposition not in DISPOSITIONS:
            raise DistillationArtifactError(f"{path}:{line} DISPOSITION has disallowed value {disposition!r}")
        if disposition == "CONDENSE":
            _require_fields(
                path,
                line,
                row,
                disposition,
                ("VALUABLE_NUCLEUS", "CONSUMER_CONTRACTS"),
            )
        elif disposition in {"ABSORB", "REMOVE"}:
            _require_fields(
                path,
                line,
                row,
                disposition,
                ("RETAINED_OWNER", "LOSS_ARGUMENT", "ORACLE", "FALSIFIER"),
            )
        elif disposition == "UNKNOWN":
            _require_fields(path, line, row, disposition, ("FALSIFIER",))
            if not row["RESIDUAL_RISK"].startswith("MISSING_EVIDENCE:"):
                raise DistillationArtifactError(f"{path}:{line} UNKNOWN requires named MISSING_EVIDENCE")
            if not row["RESIDUAL_RISK"].removeprefix("MISSING_EVIDENCE:").strip():
                raise DistillationArtifactError(f"{path}:{line} UNKNOWN requires named MISSING_EVIDENCE")


def _validated_counts(manifest: Path, capabilities: Path, findings: Path, surgery: Path | None) -> tuple[int, int, int, int]:
    manifest_rows = _validate_manifest(manifest)
    capability_ids = _read_capabilities(capabilities)
    finding_rows = _validate_findings(findings, capability_ids)
    if surgery is None:
        return len(manifest_rows), len(capability_ids), len(finding_rows), 0

    incomplete = [row["VARIANT"] for row in manifest_rows if row["SCOPE"] == "included" and any(row[column] != "complete" for column in ("PASS_1", "PASS_2", "PASS_3"))]
    if incomplete:
        raise DistillationArtifactError("surgery validation requires three complete readings for every scoped variant: " + ", ".join(incomplete))
    surgery_rows = _read_surgery(surgery)
    _validate_surgery(surgery, surgery_rows, capability_ids)
    return len(manifest_rows), len(capability_ids), len(finding_rows), len(surgery_rows)


def validate_review(manifest: Path, capabilities: Path, findings: Path, surgery: Path | None) -> None:
    """Validate structural contracts for the distillation review artifacts."""
    _validated_counts(manifest, capabilities, findings, surgery)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Validate distillation review artifacts")
    parser.add_argument("command", choices=("validate",))
    args = parser.parse_args(argv)
    if args.command != "validate":
        return 2

    surgery = DEFAULT_SURGERY if DEFAULT_SURGERY.exists() else None
    counts = _validated_counts(DEFAULT_MANIFEST, DEFAULT_CAPABILITIES, DEFAULT_FINDINGS, surgery)
    print(f"validated variants={counts[0]} capabilities={counts[1]} findings={counts[2]} surgery={counts[3]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
