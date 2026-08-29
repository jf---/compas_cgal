"""Validate the frozen measurement-claim ledger identity and state."""

from __future__ import annotations

import ast
import datetime
import hashlib
import io
import json
import pathlib
import re
import shlex
import subprocess
import tokenize
from typing import Dict
from typing import Literal
from typing import Mapping
from typing import Optional
from typing import Sequence
from typing import TypedDict
from typing import cast

from tools.measurement_artifact import ValidatedEnvelope
from tools.measurement_claim_result import GeneratorCasePayload
from tools.measurement_claim_result import GeneratorClaimPayload
from tools.measurement_claim_result import GeneratorClaimRecord
from tools.measurement_claim_result import ValidatedArtifactDirectory
from tools.measurement_claim_result import ValidatedArtifactStartedUtc
from tools.measurement_claim_result import validate_claim_artifact
from tools.measurement_claim_result import validate_generator_payload

Disposition = Literal["pending", "re-earned", "corrected", "historical", "deleted", "not-a-claim"]


class LedgerRow(TypedDict):
    """One parsed measurement-claim ledger row."""

    ordinal: str
    claim_id: str
    extracted_location: str
    stable_anchor: str
    anchor_match: str
    disposition: Disposition
    evidence: str


class _FrozenHit(TypedDict):
    path: str
    line: int
    text: str
    stable_anchor: str
    anchor_match: str


class _SourceRegion(TypedDict):
    raw: bytes
    byte_start: int
    byte_end: int
    line_start: int
    line_end: int


DISPOSITIONS = ("pending", "re-earned", "corrected", "historical", "deleted", "not-a-claim")
IN_AUDIT_STATUS = "> **status: in audit — {done}/14 Task-5 extractor rows dispositioned**"
COMPLETE_STATUS = "> **status: complete — 14/14 Task-5 extractor rows dispositioned**"
FROZEN_SOURCE_COMMIT = "eec665c1df1cd8d1e98dd9dd1001b5984e17a703"
IMMUTABLE_COLUMNS_SHA256 = "e0761be0cce4f4bb5dfcc4335ec3b54ff5751c7e1454870c2273f2147250883e"

_RADIAL_SOURCE = "src/compas_cgal/engagement_radial_toolpath.py"
_ADVANCE_SOURCE = "src/compas_cgal/engagement_toolpath.py"
_SOURCE_PATHS = (_RADIAL_SOURCE, _ADVANCE_SOURCE)
_SOURCE_FILE_SHA256 = {
    _RADIAL_SOURCE: "4a12b0d8a404a7355271eafb10baa864e33d64411a60ca3195ccf85bd1a19af0",
    _ADVANCE_SOURCE: "6bc5095fd7853949c4c6f32bcfbe7b2d85ffd038b82e6ce8cf6d259e224117cf",
}
_SOURCE_REGION_SHA256 = {
    "radial-module-which-circle": "1ab67dae4cbc8dea240b7f7a572be932eed99022427e7d773bb2442918ec7ccc",
    "radius-ladder-subdivisions": "e3898ee926c10c0dcac29a92d4811aa2f98ef1dc0fa9591b54f5b1e99a769a34",
    "radius-ladder-refinement-margin": "ba810f64945c4da67a1e994b2bab241dbba0bbb3d24241c9430f36178b1d0461",
    "gentlest-rung-peak": "19357c059c7b7f67247485b1613c8b043077ac5be332e16e4b9c8e0bde89b4b5",
    "least-bad-rung-double": "4aeee76e5df523410d1f0adc00d26fd6207d3c2f6b52937f9eceae80d6dfc311",
    "regulation-cap-angle": "5ecd9cc249e9138e4dc9a7339ae7fd6cf33bc76432fab2c6a756c5347f45ca4a",
    "measured-peak-reporting": "39fa70015f93d09c420d337aa80dbca94dbfe2dad4b86b72060c9698c5b96680",
    "loop-probe-count": "a913fa0d7ce1c67304c72c1751e445c702798c6a6a887464f37d93394b6ad2ad",
    "radius-ladder-floor-steps": "a3b0954dde586e1abd2d3399c81c05d70e2eb748c1db247f538b945217ea8128",
}
_SOURCE_REGION_PATH = {
    "radial-module-which-circle": _RADIAL_SOURCE,
    "radius-ladder-subdivisions": _RADIAL_SOURCE,
    "radius-ladder-refinement-margin": _RADIAL_SOURCE,
    "gentlest-rung-peak": _RADIAL_SOURCE,
    "least-bad-rung-double": _RADIAL_SOURCE,
    "regulation-cap-angle": _ADVANCE_SOURCE,
    "measured-peak-reporting": _ADVANCE_SOURCE,
    "loop-probe-count": _ADVANCE_SOURCE,
    "radius-ladder-floor-steps": _RADIAL_SOURCE,
}
_MANDATORY_SOURCE_REGIONS = (
    "radial-module-which-circle",
    "radius-ladder-subdivisions",
    "radius-ladder-refinement-margin",
    "gentlest-rung-peak",
    "least-bad-rung-double",
    "regulation-cap-angle",
    "measured-peak-reporting",
    "loop-probe-count",
)
_FLOOR_SOURCE_REGION = "radius-ladder-floor-steps"
_LEADING_COMMENT_ASSIGNMENTS = {
    "radius-ladder-subdivisions": "RADIUS_LADDER_SUBDIVISIONS",
    "radius-ladder-refinement-margin": "RADIUS_LADDER_REFINEMENT_MARGIN",
    "loop-probe-count": "LOOP_PROBE_COUNT",
    "radius-ladder-floor-steps": "RADIUS_LADDER_FLOOR_STEPS",
}
_DOCSTRING_REGION_BOUNDARIES = {
    "radial-module-which-circle": (
        b"emitted, counted, and warned about rather than hidden.\n\n",
        b"\nIT IS NOT FREE,",
    ),
    "gentlest-rung-peak": (
        b'    """The mildest of a station\'s already-refused candidate radii, with its measurement.\n\n',
        b"\n    Attributes:\n",
    ),
    "least-bad-rung-double": (
        b"    `tests/test_engagement_radial_toolpath.py` pins that state.\n\n",
        b"\n    CANDIDATES are",
    ),
    "regulation-cap-angle": (
        b"            `_cap_surrogate` -- the only form of the cap that reaches a predicate.\n",
        b"        tool_radius: Tool radius in model units.\n",
    ),
    "measured-peak-reporting": (
        b'    """Largest engaged-run angle REPORTED over this machining circle\'s evaluated positions.\n\n',
        b"\n    The whole ring is evaluated",
    ),
}
_DOCSTRING_REGION_OWNERS: Mapping[str, Optional[str]] = {
    "radial-module-which-circle": None,
    "gentlest-rung-peak": "_GentlestRung",
    "least-bad-rung-double": "_least_bad_rung",
    "regulation-cap-angle": "_Regulation",
    "measured-peak-reporting": "_measured_peak_engagement",
}
_COMMENT_SOURCE_REGIONS = frozenset(_LEADING_COMMENT_ASSIGNMENTS)
_MC007_DISPOSITIONS = ("re-earned", "corrected", "historical")
_RAW_SOURCE_DIFF = re.compile(r"^:100644 100644 [0-9a-f]+ [0-9a-f]+ M\t(.+)$")
_DIFF_FILE = re.compile(r"^diff --git a/(.+) b/(.+)$")
_DIFF_HUNK = re.compile(r"^@@ -(\d+)(?:,(\d+))? \+(\d+)(?:,(\d+))? @@")
_OBJECT_ID = re.compile(r"(?:[0-9a-f]{40}|[0-9a-f]{64})\Z")

_ROW_COUNT = 14
_FILE_COUNT = 5
_LEDGER_HEADER = "| ordinal | claim ID | extracted location | stable anchor | anchor match | disposition | evidence / change |"
_LEDGER_SEPARATOR = "| --- | --- | --- | --- | --- | --- | --- |"
_SECTION_HEADINGS = (
    "## Scope and limits",
    "## Frozen extraction",
    "## Disposition semantics",
    "## Ledger",
    "## Frozen matched lines",
)
_SCOPE_REQUIRED = (
    "This ledger freezes the Task-5 Python-comment regex population. It is not a",
    "claim of repository-wide invariant-I2 closure. Tasks 6 and 7 may change only",
    "the `disposition` and `evidence / change` cells; all identity columns and",
    "frozen matched lines remain immutable.",
)
_FROZEN_METADATA_REQUIRED = (
    f"Extraction source commit: `{FROZEN_SOURCE_COMMIT}`",
    "",
    "Population: 14 line hits across 5 files",
)
_EXTRACTOR_REQUIRED = (
    "```bash",
    'grep -rnE "# .*(MEASURED|[Mm]easured (on|at|against)|measures [0-9]|[0-9]+(\\.[0-9]+)?[x×] (faster|slower))" \\',
    "  src/compas_cgal benchmarks --include='*.py' | grep -v superseded",
    "```",
)
_ORDER_REQUIRED = (
    "Canonical order is `src/compas_cgal` before `benchmarks`, then `LC_ALL=C`",
    "relative-path order, numeric source line, and raw hit as the final tie-breaker.",
)
_DRIFT_REQUIRED = (
    "Execution scan (2026-08-28 UTC): exact match — 14/14 frozen hits present;",
    "0 changed, 0 missing, 0 new. Normalized stream SHA-256:",
    "`b78ac690bceda37fdebcdcdeb9731a8b213e7d946a897be45181cb95b1100f66`.",
    "This scan is drift evidence only and does not redefine the frozen population.",
)
_DISPOSITION_SEMANTICS_REQUIRED = (
    "- `pending` — final adjudication has not happened.",
    ("- `re-earned` — an authenticated rerun confirmed every material assertion and records the exact command, configuration, artifact, result digest, and full input commit."),
    ("- `corrected` — the assertion was wrong or incomplete; authenticated evidence and the source correction commit are recorded."),
    ("- `historical` — the original configuration cannot be reconstructed; the source is explicitly labelled and names the missing identity or configuration."),
    "- `deleted` — the assertion was removed while its frozen identity remains here.",
    ("- `not-a-claim` — semantic review found a regex false positive and records an explicit rationale."),
    "",
    "A case-level result does not automatically disposition every mapped row. Partial",
    "reproduction cannot become `re-earned`; every material assertion is adjudicated",
    "row by row. `reproduced` is not a ledger disposition.",
)
_EXTRACTION_PATTERN = re.compile(r"# .*(MEASURED|[Mm]easured (on|at|against)|measures [0-9]|[0-9]+(\.[0-9]+)?[x×] (faster|slower))")


class InvalidMeasurementClaimLedgerError(RuntimeError):
    """The measurement-claim ledger violates its frozen structural contract."""


def _git(repository: pathlib.Path, *arguments: str) -> bytes:
    command = ["git", "--no-replace-objects", "-C", str(repository), *arguments]
    try:
        completed = subprocess.run(command, check=False, capture_output=True)
    except OSError as exc:
        raise InvalidMeasurementClaimLedgerError(f"Git command could not start: {shlex.join(command)}") from exc
    if completed.returncode != 0:
        detail = completed.stderr.decode("utf-8", errors="replace").strip()
        raise InvalidMeasurementClaimLedgerError(f"Git command failed ({completed.returncode}): {shlex.join(command)}: {detail}")
    return completed.stdout


def _repository_root() -> pathlib.Path:
    start = pathlib.Path(__file__).resolve().parent
    raw = _git(start, "rev-parse", "--show-toplevel")
    try:
        root = pathlib.Path(raw.decode("utf-8").strip()).resolve()
    except UnicodeDecodeError as exc:
        raise InvalidMeasurementClaimLedgerError(f"Git repository path is not UTF-8: {start}") from exc
    if not root.is_dir():
        raise InvalidMeasurementClaimLedgerError(f"Git repository path is not a directory: {root}")
    return root


def _decode_source_blob(path: str, raw: bytes) -> str:
    try:
        return raw.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise InvalidMeasurementClaimLedgerError(f"Task-6 source is not UTF-8: {path}") from exc


def _source_region(raw: bytes, byte_start: int, byte_end: int, *, name: str) -> _SourceRegion:
    region = raw[byte_start:byte_end]
    if not region or not region.endswith(b"\n"):
        raise InvalidMeasurementClaimLedgerError(f"Task-6 source region {name} must be non-empty and newline-terminated")
    line_start = raw.count(b"\n", 0, byte_start) + 1
    return {
        "raw": region,
        "byte_start": byte_start,
        "byte_end": byte_end,
        "line_start": line_start,
        "line_end": line_start + region.count(b"\n") - 1,
    }


def _leading_comment_region(path: str, raw: bytes, *, name: str, assignment: str) -> _SourceRegion:
    source = _decode_source_blob(path, raw)
    try:
        tree = ast.parse(source, filename=path)
    except SyntaxError as exc:
        raise InvalidMeasurementClaimLedgerError(f"Task-6 source cannot be parsed while locating {name}: {path}") from exc
    matches = [node.lineno for node in tree.body if _assignment_name(node) == assignment]
    if len(matches) != 1:
        raise InvalidMeasurementClaimLedgerError(f"Task-6 source assignment {assignment} for {name} must occur exactly once in {path}")
    lines = raw.splitlines(keepends=True)
    assignment_index = matches[0] - 1
    start_index = assignment_index - 1
    while start_index >= 0 and (not lines[start_index].strip() or lines[start_index].lstrip().startswith(b"#")):
        start_index -= 1
    start_index += 1
    while start_index < assignment_index and not lines[start_index].strip():
        start_index += 1
    if start_index == assignment_index or any(not line.lstrip().startswith(b"#") for line in lines[start_index:assignment_index] if line.strip()):
        raise InvalidMeasurementClaimLedgerError(f"Task-6 source region {name} is not one complete leading comment in {path}")
    byte_start = sum(len(line) for line in lines[:start_index])
    byte_end = sum(len(line) for line in lines[:assignment_index])
    return _source_region(raw, byte_start, byte_end, name=name)


def _docstring_region(path: str, raw: bytes, *, name: str, before: bytes, after: bytes) -> _SourceRegion:
    if raw.count(before) != 1:
        raise InvalidMeasurementClaimLedgerError(f"Task-6 immutable docstring start boundary for {name} must occur exactly once in {path}")
    byte_start = raw.index(before) + len(before)
    byte_end = raw.find(after, byte_start)
    if byte_end < 0:
        raise InvalidMeasurementClaimLedgerError(f"Task-6 immutable docstring end boundary for {name} is missing in {path}")
    region = _source_region(raw, byte_start, byte_end, name=name)
    source = _decode_source_blob(path, raw)
    try:
        tree = ast.parse(source, filename=path)
    except SyntaxError as exc:
        raise InvalidMeasurementClaimLedgerError(f"Task-6 source cannot be parsed while binding {name}: {path}") from exc
    owner_name = _DOCSTRING_REGION_OWNERS[name]
    if owner_name is None:
        owners: list[ast.Module | ast.ClassDef | ast.FunctionDef | ast.AsyncFunctionDef] = [tree]
    else:
        owners = [node for node in ast.walk(tree) if isinstance(node, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == owner_name]
    if len(owners) != 1 or not owners[0].body:
        raise InvalidMeasurementClaimLedgerError(f"Task-6 docstring owner for {name} must occur exactly once in {path}")
    expression = owners[0].body[0]
    if not isinstance(expression, ast.Expr) or not isinstance(expression.value, ast.Constant) or not isinstance(expression.value.value, str):
        raise InvalidMeasurementClaimLedgerError(f"Task-6 region {name} is not inside its expected owner docstring")
    if expression.end_lineno is None or expression.end_col_offset is None:
        raise InvalidMeasurementClaimLedgerError(f"Task-6 docstring owner location is incomplete for {name}")
    lines = raw.splitlines(keepends=True)
    owner_start = sum(len(line) for line in lines[: expression.lineno - 1]) + expression.col_offset
    owner_end = sum(len(line) for line in lines[: expression.end_lineno - 1]) + expression.end_col_offset
    if not (owner_start <= region["byte_start"] < region["byte_end"] <= owner_end):
        raise InvalidMeasurementClaimLedgerError(f"Task-6 region {name} was relocated outside its expected owner docstring")
    return region


def _task6_source_region_spans(path: str, raw: bytes) -> Dict[str, _SourceRegion]:
    if path not in _SOURCE_PATHS:
        raise InvalidMeasurementClaimLedgerError(f"Task-6 source path is not allowlisted: {path}")
    regions: Dict[str, _SourceRegion] = {}
    for name, region_path in _SOURCE_REGION_PATH.items():
        if region_path != path:
            continue
        if name in _LEADING_COMMENT_ASSIGNMENTS:
            regions[name] = _leading_comment_region(
                path,
                raw,
                name=name,
                assignment=_LEADING_COMMENT_ASSIGNMENTS[name],
            )
        else:
            before, after = _DOCSTRING_REGION_BOUNDARIES[name]
            regions[name] = _docstring_region(path, raw, name=name, before=before, after=after)
    return regions


def _task6_source_regions(path: str, raw: bytes) -> Dict[str, bytes]:
    """Resolve exact Task-6-owned source bytes for baseline attestation tests."""
    return {name: region["raw"] for name, region in _task6_source_region_spans(path, raw).items()}


def _source_sentinel(name: str) -> bytes:
    marker = name.upper().replace("-", "_").encode("ascii")
    if name in _COMMENT_SOURCE_REGIONS:
        return b"# __TASK6_OWNED_" + marker + b"__\n"
    return b"__TASK6_OWNED_" + marker + b"__\n"


def _normalize_source(raw: bytes, regions: Mapping[str, _SourceRegion], applicable: frozenset[str]) -> bytes:
    normalized = raw
    owned = [(name, region) for name, region in regions.items() if name in applicable]
    for name, region in sorted(owned, key=lambda item: item[1]["byte_start"], reverse=True):
        normalized = normalized[: region["byte_start"]] + _source_sentinel(name) + normalized[region["byte_end"] :]
    return normalized


def _source_token_pairs(path: str, raw: bytes) -> tuple[tuple[int, str], ...]:
    try:
        return tuple((token.type, token.string) for token in tokenize.tokenize(io.BytesIO(raw).readline))
    except (IndentationError, SyntaxError, tokenize.TokenError) as exc:
        raise InvalidMeasurementClaimLedgerError(f"Task-6 normalized token stream is invalid for {path}") from exc


def _strip_leading_docstrings(node: ast.AST) -> None:
    for child in ast.iter_child_nodes(node):
        _strip_leading_docstrings(child)
    if not isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)) or not node.body:
        return
    first = node.body[0]
    if isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant) and isinstance(first.value.value, str):
        del node.body[0]


def _docstring_stripped_ast_dump(path: str, raw: bytes) -> str:
    source = _decode_source_blob(path, raw)
    try:
        tree = ast.parse(source, filename=path, type_comments=True)
    except SyntaxError as exc:
        raise InvalidMeasurementClaimLedgerError(f"Task-6 source AST is invalid for {path}") from exc
    _strip_leading_docstrings(tree)
    return ast.dump(tree, annotate_fields=True, include_attributes=False)


def _commit_parents(repository: pathlib.Path, commit: str) -> tuple[str, ...]:
    if _OBJECT_ID.fullmatch(commit) is None:
        raise InvalidMeasurementClaimLedgerError("Task-6 source commit must be one full Git object ID")
    raw = _git(repository, "cat-file", "commit", commit)
    hasher = hashlib.sha1() if len(commit) == 40 else hashlib.sha256()
    hasher.update(b"commit " + str(len(raw)).encode("ascii") + b"\0" + raw)
    if hasher.hexdigest() != commit:
        raise InvalidMeasurementClaimLedgerError("Task-6 source raw commit identity differs from its object ID")
    header, separator, _ = raw.partition(b"\n\n")
    if not separator:
        raise InvalidMeasurementClaimLedgerError("Task-6 source commit object lacks its header separator")
    parent_headers = [line[len(b"parent ") :] for line in header.splitlines() if line.startswith(b"parent ")]
    try:
        parents = tuple(parent.decode("ascii") for parent in parent_headers)
    except UnicodeDecodeError as exc:
        raise InvalidMeasurementClaimLedgerError("Task-6 source ancestry is not ASCII") from exc
    if any(_OBJECT_ID.fullmatch(parent) is None or len(parent) != len(commit) for parent in parents):
        raise InvalidMeasurementClaimLedgerError("Task-6 source parent is not one same-format full Git object ID")
    return parents


def _commit_parent(repository: pathlib.Path, correction_commit: str) -> str:
    parents = _commit_parents(repository, correction_commit)
    if len(parents) != 1:
        raise InvalidMeasurementClaimLedgerError("Task-6 source correction commit must have exactly one parent")
    return parents[0]


def _source_blob(repository: pathlib.Path, commit: str, path: str) -> bytes:
    return _git(repository, "show", f"{commit}:{path}")


def _validate_source_diff_entries(repository: pathlib.Path, parent: str, correction_commit: str) -> None:
    raw = _git(repository, "diff-tree", "--no-commit-id", "--raw", "-r", "--no-renames", parent, correction_commit)
    try:
        lines = raw.decode("ascii").splitlines()
    except UnicodeDecodeError as exc:
        raise InvalidMeasurementClaimLedgerError("Task-6 source correction raw diff is not ASCII") from exc
    paths: list[str] = []
    for line in lines:
        match = _RAW_SOURCE_DIFF.fullmatch(line)
        if match is None:
            raise InvalidMeasurementClaimLedgerError("Task-6 source correction must contain exactly two ordinary 100644 modified paths")
        paths.append(match.group(1))
    if len(paths) != len(_SOURCE_PATHS) or set(paths) != set(_SOURCE_PATHS):
        raise InvalidMeasurementClaimLedgerError(f"Task-6 source correction must contain exactly two ordinary modified paths: {_SOURCE_PATHS!r}")


def _containing_regions(
    path: str,
    line_start: int,
    line_count: int,
    regions: Mapping[str, _SourceRegion],
    applicable: frozenset[str],
) -> frozenset[str]:
    containing: set[str] = set()
    for name, region in regions.items():
        if name not in applicable or _SOURCE_REGION_PATH[name] != path:
            continue
        if line_count == 0:
            inside = region["line_start"] - 1 <= line_start <= region["line_end"]
        else:
            line_end = line_start + line_count - 1
            inside = region["line_start"] <= line_start and line_end <= region["line_end"]
        if inside:
            containing.add(name)
    return frozenset(containing)


def _validate_source_diff_hunks(
    repository: pathlib.Path,
    parent: str,
    correction_commit: str,
    baseline_regions: Mapping[str, Mapping[str, _SourceRegion]],
    candidate_regions: Mapping[str, Mapping[str, _SourceRegion]],
    applicable: frozenset[str],
) -> None:
    raw = _git(
        repository,
        "diff",
        "--unified=0",
        "--no-ext-diff",
        "--no-textconv",
        parent,
        correction_commit,
        "--",
        *_SOURCE_PATHS,
    )
    try:
        lines = raw.decode("utf-8").splitlines()
    except UnicodeDecodeError as exc:
        raise InvalidMeasurementClaimLedgerError("Task-6 zero-context source diff is not UTF-8") from exc
    current_path: Optional[str] = None
    touched: set[str] = set()
    for line in lines:
        file_match = _DIFF_FILE.fullmatch(line)
        if file_match is not None:
            old_path, new_path = file_match.groups()
            if old_path != new_path or old_path not in _SOURCE_PATHS:
                raise InvalidMeasurementClaimLedgerError("Task-6 source diff path is outside the exact two-path allowlist")
            current_path = old_path
            continue
        hunk_match = _DIFF_HUNK.match(line)
        if hunk_match is None:
            continue
        if current_path is None:
            raise InvalidMeasurementClaimLedgerError("Task-6 source diff hunk has no owning allowlisted path")
        old_start, old_count_text, new_start, new_count_text = hunk_match.groups()
        old_count = 1 if old_count_text is None else int(old_count_text)
        new_count = 1 if new_count_text is None else int(new_count_text)
        old_owned = _containing_regions(current_path, int(old_start), old_count, baseline_regions[current_path], applicable)
        new_owned = _containing_regions(current_path, int(new_start), new_count, candidate_regions[current_path], applicable)
        owners = old_owned & new_owned
        if len(owners) != 1:
            raise InvalidMeasurementClaimLedgerError(f"Task-6 zero-context diff hunk is outside one applicable owned region in {current_path}: {line}")
        touched.update(owners)
    missing = set(_MANDATORY_SOURCE_REGIONS) - touched
    if missing:
        raise InvalidMeasurementClaimLedgerError(f"Task-6 source diff lacks mandatory owned-region coverage: {sorted(missing)!r}")


def validate_task6_source_correction(
    repository: pathlib.Path,
    correction_commit: str,
    *,
    mc007_disposition: str,
) -> None:
    """Prove one committed Task-6 correction changes only its exact prose allowlist."""
    if mc007_disposition not in _MC007_DISPOSITIONS:
        raise InvalidMeasurementClaimLedgerError(f"Task-6 MC-007 disposition is invalid for source correction: {mc007_disposition!r}")
    parent = _commit_parent(repository, correction_commit)
    _validate_source_diff_entries(repository, parent, correction_commit)
    baseline_blobs = {path: _source_blob(repository, parent, path) for path in _SOURCE_PATHS}
    candidate_blobs = {path: _source_blob(repository, correction_commit, path) for path in _SOURCE_PATHS}
    for path, raw in baseline_blobs.items():
        observed = hashlib.sha256(raw).hexdigest()
        expected = _SOURCE_FILE_SHA256[path]
        if observed != expected:
            raise InvalidMeasurementClaimLedgerError(f"Task-6 baseline raw source SHA-256 differs for {path}: expected {expected}, got {observed}")

    baseline_regions = {path: _task6_source_region_spans(path, raw) for path, raw in baseline_blobs.items()}
    candidate_regions = {path: _task6_source_region_spans(path, raw) for path, raw in candidate_blobs.items()}
    for name, expected in _SOURCE_REGION_SHA256.items():
        path = _SOURCE_REGION_PATH[name]
        observed = hashlib.sha256(baseline_regions[path][name]["raw"]).hexdigest()
        if observed != expected:
            raise InvalidMeasurementClaimLedgerError(f"Task-6 baseline source-region SHA-256 differs for {name}: expected {expected}, got {observed}")

    for name in _MANDATORY_SOURCE_REGIONS:
        path = _SOURCE_REGION_PATH[name]
        if baseline_regions[path][name]["raw"] == candidate_regions[path][name]["raw"]:
            raise InvalidMeasurementClaimLedgerError(f"Task-6 mandatory source region was not corrected: {name}")
    floor_path = _SOURCE_REGION_PATH[_FLOOR_SOURCE_REGION]
    floor_changed = baseline_regions[floor_path][_FLOOR_SOURCE_REGION]["raw"] != candidate_regions[floor_path][_FLOOR_SOURCE_REGION]["raw"]
    if floor_changed and mc007_disposition != "corrected":
        raise InvalidMeasurementClaimLedgerError("Task-6 radius-ladder-floor-steps may change only when MC-007 is corrected")

    applicable = frozenset((*_MANDATORY_SOURCE_REGIONS, _FLOOR_SOURCE_REGION)) if mc007_disposition == "corrected" else frozenset(_MANDATORY_SOURCE_REGIONS)
    for path in _SOURCE_PATHS:
        normalized_baseline = _normalize_source(baseline_blobs[path], baseline_regions[path], applicable)
        normalized_candidate = _normalize_source(candidate_blobs[path], candidate_regions[path], applicable)
        if _source_token_pairs(path, normalized_baseline) != _source_token_pairs(path, normalized_candidate):
            raise InvalidMeasurementClaimLedgerError(f"Task-6 protected token stream differs outside owned source regions: {path}")
        if _docstring_stripped_ast_dump(path, baseline_blobs[path]) != _docstring_stripped_ast_dump(path, candidate_blobs[path]):
            raise InvalidMeasurementClaimLedgerError(f"Task-6 docstring-stripped AST differs in executable source: {path}")

    _validate_source_diff_hunks(
        repository,
        parent,
        correction_commit,
        baseline_regions,
        candidate_regions,
        applicable,
    )


def validate_task6_source_lineage(
    repository: pathlib.Path,
    correction_commit: str,
    execution_commit: str,
    *,
    mc007_disposition: str,
) -> None:
    """Prove strict raw ancestry and preserved Task-6 correction-region bytes."""
    if mc007_disposition not in _MC007_DISPOSITIONS:
        raise InvalidMeasurementClaimLedgerError(f"Task-6 MC-007 disposition is invalid for source lineage: {mc007_disposition!r}")
    if correction_commit == execution_commit:
        raise InvalidMeasurementClaimLedgerError("Task-6 correction must be a strict ancestor of the distinct execution commit")
    if len(correction_commit) != len(execution_commit):
        raise InvalidMeasurementClaimLedgerError("Task-6 correction and execution commits must use the same object-ID format")
    _commit_parents(repository, correction_commit)
    pending = list(_commit_parents(repository, execution_commit))
    visited: set[str] = set()
    while pending:
        candidate = pending.pop()
        if candidate == correction_commit:
            break
        if candidate not in visited:
            visited.add(candidate)
            pending.extend(_commit_parents(repository, candidate))
    else:
        raise InvalidMeasurementClaimLedgerError("Task-6 source correction is not a raw-parent ancestor of the execution commit")

    correction_parent = _commit_parent(repository, correction_commit)
    floor_path = _SOURCE_REGION_PATH[_FLOOR_SOURCE_REGION]
    parent_floor = _task6_source_region_spans(floor_path, _source_blob(repository, correction_parent, floor_path))[_FLOOR_SOURCE_REGION]["raw"]
    correction_floor = _task6_source_region_spans(floor_path, _source_blob(repository, correction_commit, floor_path))[_FLOOR_SOURCE_REGION]["raw"]
    applicable = frozenset((*_MANDATORY_SOURCE_REGIONS, _FLOOR_SOURCE_REGION)) if parent_floor != correction_floor else frozenset(_MANDATORY_SOURCE_REGIONS)
    for path in _SOURCE_PATHS:
        corrected = _task6_source_region_spans(path, _source_blob(repository, correction_commit, path))
        executed = _task6_source_region_spans(path, _source_blob(repository, execution_commit, path))
        for name in applicable:
            if _SOURCE_REGION_PATH[name] == path and corrected[name]["raw"] != executed[name]["raw"]:
                raise InvalidMeasurementClaimLedgerError(f"Task-6 protected correction region differs at execution: {name}")


def _decode_frozen_blob(repository: pathlib.Path, path: str) -> str:
    raw = _git(repository, "show", f"{FROZEN_SOURCE_COMMIT}:{path}")
    try:
        return raw.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise InvalidMeasurementClaimLedgerError(f"frozen source blob is not UTF-8: {path}") from exc


def _frozen_python_paths(repository: pathlib.Path) -> tuple[str, ...]:
    raw = _git(
        repository,
        "ls-tree",
        "-r",
        "--name-only",
        FROZEN_SOURCE_COMMIT,
        "--",
        "src/compas_cgal",
        "benchmarks",
    )
    try:
        paths = raw.decode("utf-8").splitlines()
    except UnicodeDecodeError as exc:
        raise InvalidMeasurementClaimLedgerError("frozen Git path list is not UTF-8") from exc
    selected = [path for path in paths if path.endswith(".py")]
    return tuple(sorted(selected, key=lambda path: (0 if path.startswith("src/compas_cgal/") else 1, path)))


def _assignment_name(node: ast.stmt) -> Optional[str]:
    target: Optional[ast.expr] = None
    if isinstance(node, ast.Assign) and len(node.targets) == 1:
        target = node.targets[0]
    elif isinstance(node, ast.AnnAssign):
        target = node.target
    if isinstance(target, ast.Name):
        return target.id
    return None


def _module_name(path: str) -> str:
    without_suffix = path[:-3]
    if without_suffix.startswith("src/"):
        without_suffix = without_suffix[len("src/") :]
    return without_suffix.replace("/", ".")


def _anchor_for_hit(path: str, source: str, line: int) -> str:
    try:
        tree = ast.parse(source, filename=path)
    except SyntaxError as exc:
        raise InvalidMeasurementClaimLedgerError(f"frozen source cannot be parsed for stable anchors: {path}") from exc
    assignments = [(node.lineno, name) for node in tree.body if (name := _assignment_name(node)) is not None and node.lineno > line]
    if not assignments:
        raise InvalidMeasurementClaimLedgerError(f"frozen hit has no following module assignment anchor: {path}:{line}")
    assignment_line, name = min(assignments)
    intervening = source.splitlines()[line : assignment_line - 1]
    if any(text.strip() and not text.lstrip().startswith("#") for text in intervening):
        raise InvalidMeasurementClaimLedgerError(f"frozen hit is not a leading comment for {name}: {path}:{line}")
    return f"{_module_name(path)}.{name}"


def _frozen_hits(repository: pathlib.Path) -> tuple[_FrozenHit, ...]:
    hits: list[_FrozenHit] = []
    for path in _frozen_python_paths(repository):
        source = _decode_frozen_blob(repository, path)
        for line_number, text in enumerate(source.splitlines(), start=1):
            if "superseded" in text or _EXTRACTION_PATTERN.search(text) is None:
                continue
            hits.append(
                {
                    "path": path,
                    "line": line_number,
                    "text": text,
                    "stable_anchor": _anchor_for_hit(path, source, line_number),
                    "anchor_match": "",
                }
            )
    hits.sort(
        key=lambda hit: (
            0 if hit["path"].startswith("src/compas_cgal/") else 1,
            hit["path"],
            hit["line"],
            hit["text"],
        )
    )
    if len(hits) != _ROW_COUNT or len({hit["path"] for hit in hits}) != _FILE_COUNT:
        raise InvalidMeasurementClaimLedgerError(
            f"frozen extraction must produce {_ROW_COUNT} hits across {_FILE_COUNT} files, got {len(hits)} hits across {len({hit['path'] for hit in hits})} files"
        )
    totals = {anchor: sum(hit["stable_anchor"] == anchor for hit in hits) for anchor in {hit["stable_anchor"] for hit in hits}}
    seen: dict[str, int] = {}
    for hit in hits:
        anchor = hit["stable_anchor"]
        seen[anchor] = seen.get(anchor, 0) + 1
        hit["anchor_match"] = f"{seen[anchor]}/{totals[anchor]}"
    return tuple(hits)


def _read_ledger(ledger: pathlib.Path) -> str:
    try:
        return ledger.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError) as exc:
        raise InvalidMeasurementClaimLedgerError(f"measurement claim ledger is not readable UTF-8: {ledger}") from exc


def _unique_line_index(lines: Sequence[str], value: str, *, field: str) -> int:
    indices = [index for index, line in enumerate(lines) if line == value]
    if len(indices) != 1:
        raise InvalidMeasurementClaimLedgerError(f"ledger {field} must occur exactly once: {value}")
    return indices[0]


def _section(lines: Sequence[str], heading: str, next_heading: Optional[str]) -> tuple[str, ...]:
    start = _unique_line_index(lines, heading, field="section heading") + 1
    end = len(lines) if next_heading is None else _unique_line_index(lines, next_heading, field="section heading")
    return tuple(lines[start:end])


def _block_index(section: Sequence[str], block: Sequence[str], *, field: str) -> int:
    width = len(block)
    indices = [index for index in range(len(section) - width + 1) if tuple(section[index : index + width]) == tuple(block)]
    if len(indices) != 1:
        raise InvalidMeasurementClaimLedgerError(f"ledger {field} must occur exactly once in its owning section")
    return indices[0]


def _validate_page_frame(text: str) -> tuple[str, ...]:
    lines = tuple(text.splitlines())
    opening = (
        "# Measurement-claim ledger",
        "",
        lines[2] if len(lines) > 2 else "",
        "",
        "- Opened (UTC): `2026-08-28`",
        "- Programme: coherence Wave 1 backlog A",
        "",
        "## Scope and limits",
    )
    if len(lines) < len(opening) or tuple(lines[: len(opening)]) != opening or not opening[2].startswith("> **status:"):
        raise InvalidMeasurementClaimLedgerError("ledger opening status and metadata layout is malformed")

    actual_headings = tuple(line for line in lines if line.startswith("## "))
    if actual_headings != _SECTION_HEADINGS:
        raise InvalidMeasurementClaimLedgerError(f"ledger section headings must be exact and ordered: expected {_SECTION_HEADINGS!r}, got {actual_headings!r}")

    scope = _section(lines, "## Scope and limits", "## Frozen extraction")
    _block_index(scope, _SCOPE_REQUIRED, field="scope and I2 limitation")

    extraction = _section(lines, "## Frozen extraction", "## Disposition semantics")
    extraction_blocks = (
        ("frozen extraction metadata", _FROZEN_METADATA_REQUIRED),
        ("canonical extractor", _EXTRACTOR_REQUIRED),
        ("canonical order", _ORDER_REQUIRED),
        ("separate drift report", _DRIFT_REQUIRED),
    )
    positions = [_block_index(extraction, block, field=field) for field, block in extraction_blocks]
    if positions != sorted(positions):
        raise InvalidMeasurementClaimLedgerError("ledger extractor, canonical order, and drift report are out of order")

    semantics = _section(lines, "## Disposition semantics", "## Ledger")
    _block_index(semantics, _DISPOSITION_SEMANTICS_REQUIRED, field="disposition semantics")
    return lines


def _unwrap_code(value: str, *, row: int, field: str) -> str:
    if len(value) < 3 or not value.startswith("`") or not value.endswith("`") or "`" in value[1:-1]:
        raise InvalidMeasurementClaimLedgerError(f"ledger row {row:03d} {field} must be one exact inline-code span")
    return value[1:-1]


def _parse_table(lines: Sequence[str]) -> tuple[LedgerRow, ...]:
    header_index = _unique_line_index(lines, _LEDGER_HEADER, field="ledger table heading")
    if header_index + 1 >= len(lines) or lines[header_index + 1] != _LEDGER_SEPARATOR:
        raise InvalidMeasurementClaimLedgerError("ledger table separator is missing or malformed")
    body: list[str] = []
    for line in lines[header_index + 2 :]:
        if not line.startswith("|"):
            break
        body.append(line)
    if len(body) != _ROW_COUNT:
        raise InvalidMeasurementClaimLedgerError(f"ledger table must contain {_ROW_COUNT} body rows, got {len(body)}")
    rows: list[LedgerRow] = []
    for expected, line in enumerate(body, start=1):
        if not line.endswith("|"):
            raise InvalidMeasurementClaimLedgerError(f"ledger row {expected:03d} is not pipe-delimited")
        cells = [cell.strip() for cell in line[1:-1].split("|")]
        if len(cells) != 7:
            raise InvalidMeasurementClaimLedgerError(f"ledger row {expected:03d} column count must be 7, got {len(cells)}")
        ordinal, claim_id, location_cell, anchor_cell, anchor_match, disposition_text, evidence = cells
        expected_ordinal = f"{expected:03d}"
        expected_id = f"MC-{expected:03d}"
        if ordinal != expected_ordinal or claim_id != expected_id:
            raise InvalidMeasurementClaimLedgerError(f"ledger row {expected:03d} must use ordinal {expected_ordinal} and claim ID {expected_id}")
        location = _unwrap_code(location_cell, row=expected, field="extracted location")
        anchor_suffix = " / leading comment"
        if not anchor_cell.endswith(anchor_suffix):
            raise InvalidMeasurementClaimLedgerError(f"ledger row {expected:03d} stable anchor suffix is malformed")
        anchor = _unwrap_code(anchor_cell[: -len(anchor_suffix)], row=expected, field="stable anchor")
        if disposition_text not in DISPOSITIONS:
            raise InvalidMeasurementClaimLedgerError(f"ledger row {expected:03d} disposition is empty or unknown: {disposition_text!r}")
        disposition = cast(Disposition, disposition_text)
        if disposition == "pending" and evidence != "—":
            raise InvalidMeasurementClaimLedgerError(f"ledger row {expected:03d} pending disposition must have em-dash evidence")
        if disposition != "pending" and (not evidence or evidence == "—"):
            raise InvalidMeasurementClaimLedgerError(f"ledger row {expected:03d} terminal disposition requires evidence")
        rows.append(
            {
                "ordinal": ordinal,
                "claim_id": claim_id,
                "extracted_location": location,
                "stable_anchor": anchor,
                "anchor_match": anchor_match,
                "disposition": disposition,
                "evidence": evidence,
            }
        )
    return tuple(rows)


def _immutable_digest(rows: Sequence[LedgerRow]) -> str:
    immutable = [
        {
            "ordinal": row["ordinal"],
            "claim_id": row["claim_id"],
            "extracted_location": row["extracted_location"],
            "stable_anchor": row["stable_anchor"],
            "anchor_match": row["anchor_match"],
        }
        for row in rows
    ]
    raw = json.dumps(immutable, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _validate_rows(rows: Sequence[LedgerRow], hits: Sequence[_FrozenHit]) -> None:
    if len({row["claim_id"] for row in rows}) != _ROW_COUNT:
        raise InvalidMeasurementClaimLedgerError("ledger contains a duplicate claim ID")
    if len({row["extracted_location"] for row in rows}) != _ROW_COUNT:
        raise InvalidMeasurementClaimLedgerError("ledger contains a duplicate extracted location")
    for index, (row, hit) in enumerate(zip(rows, hits), start=1):
        expected = (
            f"{index:03d}",
            f"MC-{index:03d}",
            f"{hit['path']}:{hit['line']}",
            hit["stable_anchor"],
            hit["anchor_match"],
        )
        actual = (
            row["ordinal"],
            row["claim_id"],
            row["extracted_location"],
            row["stable_anchor"],
            row["anchor_match"],
        )
        if actual != expected:
            raise InvalidMeasurementClaimLedgerError(f"ledger row {index:03d} immutable identity differs from frozen source: expected {expected!r}, got {actual!r}")
    digest = _immutable_digest(rows)
    if digest != IMMUTABLE_COLUMNS_SHA256:
        raise InvalidMeasurementClaimLedgerError(f"ledger immutable-column digest differs: expected {IMMUTABLE_COLUMNS_SHA256}, got {digest}")


def _validate_frozen_blocks(lines: Sequence[str], hits: Sequence[_FrozenHit]) -> None:
    expected_headings = tuple(f"### MC-{index:03d}" for index in range(1, _ROW_COUNT + 1))
    actual_headings = tuple(line for line in lines if line.startswith("### MC-"))
    if actual_headings != expected_headings:
        raise InvalidMeasurementClaimLedgerError(
            f"ledger claim headings must be exact and confined to the frozen section: expected {expected_headings!r}, got {actual_headings!r}"
        )
    expected_section: list[str] = [""]
    for index, hit in enumerate(hits, start=1):
        expected_section.extend((f"### MC-{index:03d}", "", "```text", hit["text"], "```"))
        if index != _ROW_COUNT:
            expected_section.append("")
    actual_section = _section(lines, "## Frozen matched lines", None)
    if actual_section != tuple(expected_section):
        raise InvalidMeasurementClaimLedgerError("ledger frozen blocks must be exact and wholly contained by the Frozen matched lines section")


def _validate_status(lines: Sequence[str], rows: Sequence[LedgerRow]) -> None:
    status_lines = [line for line in lines if line.startswith("> **status:")]
    done = sum(row["disposition"] != "pending" for row in rows)
    expected = COMPLETE_STATUS if done == _ROW_COUNT else IN_AUDIT_STATUS.format(done=done)
    if status_lines != [expected]:
        raise InvalidMeasurementClaimLedgerError(f"ledger status must equal observed disposition count byte-for-byte: expected {expected!r}, got {status_lines!r}")


def _artifact_text(
    payload: GeneratorClaimPayload,
    envelope: ValidatedEnvelope,
    *,
    started: ValidatedArtifactStartedUtc,
    artifact_directory: ValidatedArtifactDirectory,
) -> str:
    if type(started) is not datetime.datetime or started.tzinfo is None or started.utcoffset() != datetime.timedelta(0):
        raise InvalidMeasurementClaimLedgerError("Task-6 artifact started value must be one timezone-aware UTC datetime")
    if started > envelope.finished:
        raise InvalidMeasurementClaimLedgerError("Task-6 artifact started value must not follow the authenticated finish")
    if payload["source_commit"] != envelope.commit:
        raise InvalidMeasurementClaimLedgerError("Task-6 payload source commit differs from the authenticated artifact commit")
    expected = pathlib.PurePosixPath(
        "benchmarks",
        "measurement_claim_results",
        f"{started.date().isoformat()}-{str(envelope.commit)[:12]}-generator-{str(envelope.input_sha256)[:12]}",
    )
    if type(artifact_directory) is not pathlib.PurePosixPath or artifact_directory != expected:
        raise InvalidMeasurementClaimLedgerError(f"Task-6 artifact directory must equal the authenticated repository-relative path: {expected}")
    text = artifact_directory.as_posix()
    if any(token in text for token in ("|", "\r", "\n", "\u2028", "\u2029")):
        raise InvalidMeasurementClaimLedgerError("Task-6 artifact directory is not Markdown-table-safe")
    return text


def _canonical_json(value: object, *, field: str) -> str:
    try:
        return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False)
    except (TypeError, ValueError) as exc:
        raise InvalidMeasurementClaimLedgerError(f"Task-6 {field} is not canonical finite JSON") from exc


def _render_validated_ledger_evidence(
    payload: GeneratorClaimPayload,
    envelope: ValidatedEnvelope,
    *,
    started: ValidatedArtifactStartedUtc,
    artifact_directory: ValidatedArtifactDirectory,
) -> Dict[str, str]:
    artifact = _artifact_text(payload, envelope, started=started, artifact_directory=artifact_directory)
    cases: Dict[str, GeneratorCasePayload] = {case["case"]: case for case in payload["cases"]}
    rendered: Dict[str, str] = {}
    for claim in payload["claims"]:
        claim_id = claim["claim_id"]
        case = cases[claim["case"]]
        cell = (
            f"artifact={artifact}; claim={claim_id}; input={envelope.input_sha256}; result={envelope.result_sha256}; "
            f"source={payload['source_commit']}; source_correction={payload['source_correction_commit']}; "
            f"case={claim['case']}; disposition={claim['disposition']}; reason={claim['reason']}; "
            f"config={_canonical_json(case['config'], field=f'{claim_id} config')}; "
            f"evidence={_canonical_json(claim['evidence'], field=f'{claim_id} evidence')}; "
            f"selection={_canonical_json(claim['selection_decision_provenance'], field=f'{claim_id} selection')}; "
            "continuous_certificate=null"
        )
        if any(token in cell for token in ("|", "\r", "\n", "\u2028", "\u2029")):
            raise InvalidMeasurementClaimLedgerError(f"Task-6 {claim_id} evidence must be Markdown-safe and one physical line")
        rendered[claim_id] = cell
    return rendered


def render_ledger_evidence(
    payload: object,
    envelope: ValidatedEnvelope,
    *,
    started: ValidatedArtifactStartedUtc,
    artifact_directory: ValidatedArtifactDirectory,
) -> Mapping[str, str]:
    """Render the ten authenticated generator-claim ledger evidence cells."""
    validated = validate_generator_payload(payload)
    return _render_validated_ledger_evidence(
        validated,
        envelope,
        started=started,
        artifact_directory=artifact_directory,
    )


def _validate_task6_ledger_rows(
    rows: Sequence[LedgerRow],
    payload: GeneratorClaimPayload,
    envelope: ValidatedEnvelope,
    *,
    started: ValidatedArtifactStartedUtc,
    artifact_directory: ValidatedArtifactDirectory,
) -> None:
    if len(rows) != _ROW_COUNT:
        raise InvalidMeasurementClaimLedgerError(f"Task-6 ledger acceptance requires exactly {_ROW_COUNT} rows")
    evidence = _render_validated_ledger_evidence(
        payload,
        envelope,
        started=started,
        artifact_directory=artifact_directory,
    )
    claims: Dict[str, GeneratorClaimRecord] = {claim["claim_id"]: claim for claim in payload["claims"]}
    for index, row in enumerate(rows, start=1):
        claim_id = f"MC-{index:03d}"
        if row["ordinal"] != f"{index:03d}" or row["claim_id"] != claim_id:
            raise InvalidMeasurementClaimLedgerError(f"Task-6 ledger row {index:03d} identity is not canonical")
        if index <= 10:
            claim = claims[claim_id]
            if row["disposition"] != claim["disposition"]:
                raise InvalidMeasurementClaimLedgerError(f"Task-6 {claim_id} disposition differs from the authenticated payload")
            if row["evidence"] != evidence[claim_id]:
                raise InvalidMeasurementClaimLedgerError(f"Task-6 {claim_id} evidence is not byte-equal to the authenticated rendering")
        elif row["disposition"] != "pending" or row["evidence"] != "—":
            raise InvalidMeasurementClaimLedgerError(f"Task-6 {claim_id} must remain pending with em-dash evidence")


def validate_ledger_evidence(ledger: pathlib.Path, artifact_directories: Sequence[pathlib.Path]) -> None:
    """Authenticate one Task-6 artifact and require its exact ledger projection."""
    if len(artifact_directories) != 1:
        raise InvalidMeasurementClaimLedgerError("Task-6 ledger acceptance requires exactly one authenticated artifact")
    rows = validate_ledger_structure(ledger)
    payload, envelope, started, artifact_directory = validate_claim_artifact(artifact_directories[0])
    mc007 = next(claim for claim in payload["claims"] if claim["claim_id"] == "MC-007")
    repository = artifact_directories[0].absolute().parent.parent.parent
    validate_task6_source_correction(
        repository,
        str(payload["source_correction_commit"]),
        mc007_disposition=mc007["disposition"],
    )
    validate_task6_source_lineage(
        repository,
        str(payload["source_correction_commit"]),
        str(payload["source_commit"]),
        mc007_disposition=mc007["disposition"],
    )
    _validate_task6_ledger_rows(
        rows,
        payload,
        envelope,
        started=started,
        artifact_directory=artifact_directory,
    )


def validate_ledger_structure(ledger: pathlib.Path) -> tuple[LedgerRow, ...]:
    """Validate the sole Markdown claim ledger against its frozen Git source.

    Args:
        ledger: Markdown ledger to parse.

    Returns:
        The ordered, validated ledger rows.

    Raises:
        InvalidMeasurementClaimLedgerError: The page, Git source, or frozen
            identity violates the contract.
    """
    text = _read_ledger(ledger)
    lines = _validate_page_frame(text)
    rows = _parse_table(lines)
    hits = _frozen_hits(_repository_root())
    _validate_rows(rows, hits)
    _validate_frozen_blocks(lines, hits)
    _validate_status(lines, rows)
    return rows
