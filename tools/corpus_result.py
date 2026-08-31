"""Validate and select committed benchmark-corpus result bundles."""

from __future__ import annotations

import datetime
import pathlib
import re

from benchmarks.cli import AGGREGATE_CORPUS_NAMES
from benchmarks.cli import DEFAULT_CAP_DEG
from benchmarks.cli import DEFAULT_TOOL_DIAMETER
from benchmarks.errors import MalformedRecordError
from benchmarks.measurement import MeasurementRecord
from benchmarks.report import JSON_NAME
from benchmarks.report import MARKDOWN_NAME
from tools import measurement_artifact
from tools import measurement_claim_json
from tools.measurement_artifact import ArtifactKind
from tools.measurement_artifact import InvalidMeasurementEnvelopeError
from tools.measurement_artifact import MeasurementArtifactError

ARTIFACT_KIND = ArtifactKind("benchmark-corpus-result/v1")
INPUT_VERSION = "benchmark-corpus-input/v1"
RESULT_VERSION = "benchmark-corpus-result/v1"
_RESULT_NAME = re.compile(r"(?P<date>\d{4}-\d{2}-\d{2})-(?P<commit>[0-9a-f]{12})\Z")


class NoMeasuredResultError(MeasurementArtifactError):
    """A corpus result root has no visible immutable bundle."""


class InvalidMeasuredResultError(MeasurementArtifactError):
    """A visible corpus result violates its payload contract."""


def _repository_for_results_root(results_root: pathlib.Path) -> pathlib.Path:
    root = results_root.resolve()
    if root.name != "results" or root.parent.name != "benchmarks":
        raise InvalidMeasuredResultError(f"corpus results root must be <repository>/benchmarks/results: {root}")
    candidate = root.parents[1]
    try:
        repository = measurement_artifact._resolve_repository(candidate)
    except measurement_artifact.MeasurementGitError as exc:
        raise InvalidMeasuredResultError(f"corpus result is not inside a Git worktree: {root}") from exc
    if repository != candidate:
        raise InvalidMeasuredResultError(f"corpus results root is not owned by its declared repository: {root}")
    return repository


def _read_stamp(result: pathlib.Path) -> dict[str, object]:
    try:
        value = measurement_artifact._decode_strict_json_bytes((result / measurement_artifact.STAMP_NAME).read_bytes(), field="stamp")
    except InvalidMeasurementEnvelopeError as exc:
        raise InvalidMeasuredResultError(f"invalid corpus stamp: {result}") from exc
    if type(value) is not dict:
        raise InvalidMeasuredResultError(f"corpus stamp must be an object: {result}")
    return value


def _expected_input_payload() -> dict[str, object]:
    return {
        "name": "all",
        "tool_diameter": DEFAULT_TOOL_DIAMETER,
        "cap_deg": DEFAULT_CAP_DEG,
        "collect_digits": True,
        "aggregate_corpora": list(AGGREGATE_CORPUS_NAMES),
    }


def _validate_input_payload(value: object, *, result: pathlib.Path) -> None:
    expected_keys = {"name", "tool_diameter", "cap_deg", "collect_digits", "aggregate_corpora"}
    if type(value) is not dict or set(value) != expected_keys:
        raise InvalidMeasuredResultError(f"corpus input payload shape is invalid: {result}")
    if (
        type(value["name"]) is not str
        or type(value["tool_diameter"]) is not float
        or type(value["cap_deg"]) is not float
        or type(value["collect_digits"]) is not bool
        or type(value["aggregate_corpora"]) is not list
        or any(type(name) is not str for name in value["aggregate_corpora"])
        or value != _expected_input_payload()
    ):
        raise InvalidMeasuredResultError(f"corpus input payload does not match committed all-corpus defaults: {result}")


def _validate_argv(value: object, *, logical_name: str, result: pathlib.Path) -> None:
    if type(value) is not list or len(value) != 8 or any(type(argument) is not str or not argument for argument in value):
        raise InvalidMeasuredResultError(f"corpus argv shape is invalid: {result}")
    if value[1:7] != ["-m", "benchmarks.cli", "corpus", "--name", "all", "--out"]:
        raise InvalidMeasuredResultError(f"corpus argv is not the fixed all-corpus command: {result}")
    output = pathlib.PurePosixPath(value[7])
    expected_parent = pathlib.PurePosixPath("benchmarks/results")
    if ".." in output.parts or output.parent != expected_parent or not output.name.startswith(f".{logical_name}.stage-"):
        raise InvalidMeasuredResultError(f"corpus argv does not name its hidden repository-relative stage: {result}")


def _validate_corpus_payloads(result: pathlib.Path) -> None:
    markdown_path = result / MARKDOWN_NAME
    json_path = result / JSON_NAME
    for path in (markdown_path, json_path):
        if not path.is_file() or path.is_symlink():
            raise InvalidMeasuredResultError(f"corpus payload must be a regular non-symlink file: {path}")
        if path.stat().st_size == 0:
            raise InvalidMeasuredResultError(f"corpus payload is empty: {path}")
    try:
        markdown = markdown_path.read_text(encoding="utf-8")
    except UnicodeDecodeError as exc:
        raise InvalidMeasuredResultError(f"corpus Markdown is not UTF-8: {markdown_path}") from exc
    if not markdown.startswith("# Benchmark corpus result"):
        raise InvalidMeasuredResultError(f"corpus Markdown has the wrong heading: {markdown_path}")
    records = measurement_claim_json.decode_strict(json_path.read_bytes(), str(json_path), InvalidMeasuredResultError)
    if type(records) is not list or not records:
        raise InvalidMeasuredResultError(f"corpus JSON must be a non-empty exact list: {json_path}")
    for row in records:
        if type(row) is not dict:
            raise InvalidMeasuredResultError(f"corpus JSON rows must be exact objects: {json_path}")
        try:
            MeasurementRecord.from_dict(row)
        except MalformedRecordError as exc:
            raise InvalidMeasuredResultError(f"corpus row violates the exact-column schema: {json_path}") from exc


def _validate_result(result: pathlib.Path, *, repository: pathlib.Path, logical_name: str) -> datetime.datetime:
    match = _RESULT_NAME.fullmatch(logical_name)
    if match is None:
        raise InvalidMeasuredResultError(f"invalid corpus result directory name: {logical_name}")
    try:
        validated = measurement_artifact.validate_envelope(
            result,
            logical_name=logical_name,
            artifact_kind=ARTIFACT_KIND,
            repository=repository,
        )
    except InvalidMeasurementEnvelopeError as exc:
        raise InvalidMeasuredResultError(f"invalid corpus measurement envelope: {result}") from exc
    stamp = _read_stamp(result)
    if match.group("commit") != str(validated.commit)[:12]:
        raise InvalidMeasuredResultError(f"corpus directory commit prefix disagrees with stamp: {result}")
    try:
        started = datetime.datetime.fromisoformat(str(stamp["started"]))
    except (KeyError, ValueError) as exc:
        raise InvalidMeasuredResultError(f"corpus start timestamp is invalid: {result}") from exc
    if match.group("date") != started.date().isoformat():
        raise InvalidMeasuredResultError(f"corpus directory date disagrees with start timestamp: {result}")
    input_identity = stamp.get("input_identity")
    result_identity = stamp.get("result_identity")
    if type(input_identity) is not dict or input_identity.get("version") != INPUT_VERSION:
        raise InvalidMeasuredResultError(f"corpus input identity version is invalid: {result}")
    if type(result_identity) is not dict or result_identity.get("version") != RESULT_VERSION:
        raise InvalidMeasuredResultError(f"corpus result identity version is invalid: {result}")
    _validate_input_payload(input_identity.get("payload"), result=result)
    _validate_argv(stamp.get("argv"), logical_name=logical_name, result=result)
    _validate_corpus_payloads(result)
    return validated.finished


def validate_result(result: pathlib.Path) -> datetime.datetime:
    """Validate one authenticated corpus bundle and return its finish time."""
    if result.is_symlink():
        raise InvalidMeasuredResultError(f"corpus result directory must not be a symlink: {result}")
    resolved = result.resolve()
    return _validate_result(resolved, repository=_repository_for_results_root(resolved.parent), logical_name=resolved.name)


def latest_result(results: pathlib.Path) -> pathlib.Path:
    """Return the latest valid visible result by finish time then name."""
    root = results.resolve()
    if not root.is_dir() or root.is_symlink():
        raise NoMeasuredResultError(f"corpus result root holds no run directory: {root}")
    visible = [child for child in root.iterdir() if not child.name.startswith(".")]
    if not visible:
        raise NoMeasuredResultError(f"corpus result root holds no run directory: {root}")
    repository = _repository_for_results_root(root)
    candidates: list[tuple[datetime.datetime, str, pathlib.Path]] = []
    for child in visible:
        if not child.is_dir() or child.is_symlink():
            raise InvalidMeasuredResultError(f"visible corpus result child is not a real directory: {child}")
        finished = _validate_result(child.resolve(), repository=repository, logical_name=child.name)
        candidates.append((finished, child.name, child.resolve()))
    return max(candidates)[2]
