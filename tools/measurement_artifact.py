"""Content-addressed measurement provenance and atomic bundle publication.

This module authenticates reporting artifacts. It does not certify geometric
decisions or interpret any benchmark-family payload.
"""

from __future__ import annotations

import contextlib
import datetime
import errno
import hashlib
import json
import pathlib
import platform
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from types import MappingProxyType
from typing import Iterator
from typing import Mapping
from typing import NewType
from typing import Optional
from typing import Sequence

GitObjectId = NewType("GitObjectId", str)
Sha256Hex = NewType("Sha256Hex", str)
ArtifactKind = NewType("ArtifactKind", str)
IdentityVersion = NewType("IdentityVersion", str)

ENVELOPE_VERSION = "measurement-artifact-envelope/v1"
BUILD_IDENTITY_VERSION = "measurement-build-identity/v1"
STAMP_NAME = "stamp.json"
ENVELOPE_KEYS = (
    "envelope_version",
    "artifact_kind",
    "commit",
    "dirty",
    "python",
    "platform",
    "started",
    "finished",
    "argv",
    "build_identity",
    "input_identity",
    "result_identity",
)

_OBJECT_ID = re.compile(r"(?:[0-9a-f]{40}|[0-9a-f]{64})\Z")
_SHA256 = re.compile(r"[0-9a-f]{64}\Z")


class MeasurementArtifactError(RuntimeError):
    """Base failure for authenticated measurement artifacts."""


class MeasurementGitError(MeasurementArtifactError):
    """Git could not provide required source provenance."""


class InvalidGitObjectIdError(MeasurementArtifactError):
    """A Git object identity is abbreviated or non-canonical."""


class DirtyMeasurementTreeError(MeasurementArtifactError):
    """Measurement input contains worktree or index changes."""


class MeasurementInputChangedError(MeasurementArtifactError):
    """The source commit changed while a measurement was running."""


class MeasurementArtifactCollisionError(MeasurementArtifactError):
    """The intended immutable result directory already exists."""


class InvalidMeasurementEnvelopeError(MeasurementArtifactError):
    """A measurement envelope or authenticated payload is malformed."""


def _validate_object_id(value: object) -> GitObjectId:
    if type(value) is not str or _OBJECT_ID.fullmatch(value) is None:
        raise InvalidGitObjectIdError(f"Git object ID must be 40 or 64 lowercase hexadecimal characters: {value!r}")
    return GitObjectId(value)


def _validate_sha256(value: object, *, field: str) -> Sha256Hex:
    if type(value) is not str or _SHA256.fullmatch(value) is None:
        raise InvalidMeasurementEnvelopeError(f"{field} must be 64 lowercase hexadecimal characters")
    return Sha256Hex(value)


def _validate_repository_path(repository: pathlib.Path) -> pathlib.Path:
    if not repository.is_absolute() or repository != repository.resolve():
        raise MeasurementArtifactError(f"repository must be an absolute resolved path: {repository}")
    return repository


def _validate_utc(value: datetime.datetime, *, field: str) -> datetime.datetime:
    if value.tzinfo is None or value.utcoffset() != datetime.timedelta(0):
        raise InvalidMeasurementEnvelopeError(f"{field} must be timezone-aware UTC")
    return value


@dataclass(frozen=True)
class SourceSnapshot:
    """Committed build input captured before a measurement."""

    repository: pathlib.Path
    commit: GitObjectId
    pixi_lock_sha256: Sha256Hex

    def __post_init__(self) -> None:
        object.__setattr__(self, "repository", _validate_repository_path(self.repository))
        object.__setattr__(self, "commit", _validate_object_id(self.commit))
        object.__setattr__(self, "pixi_lock_sha256", _validate_sha256(self.pixi_lock_sha256, field="pixi_lock_sha256"))

    @classmethod
    def build(cls, *, repository: pathlib.Path, commit: str, pixi_lock_sha256: str) -> SourceSnapshot:
        """Validate and build a source snapshot."""
        return cls(repository=repository, commit=GitObjectId(commit), pixi_lock_sha256=Sha256Hex(pixi_lock_sha256))


@dataclass(frozen=True)
class ValidatedEnvelope:
    """Identity values reconstructed from an authenticated bundle."""

    finished: datetime.datetime
    commit: GitObjectId
    input_sha256: Sha256Hex
    result_sha256: Sha256Hex
    payload_sha256: Mapping[str, Sha256Hex]

    def __post_init__(self) -> None:
        object.__setattr__(self, "finished", _validate_utc(self.finished, field="finished"))
        object.__setattr__(self, "commit", _validate_object_id(self.commit))
        object.__setattr__(self, "input_sha256", _validate_sha256(self.input_sha256, field="input_sha256"))
        object.__setattr__(self, "result_sha256", _validate_sha256(self.result_sha256, field="result_sha256"))
        copied = {name: _validate_sha256(digest, field=f"payload_sha256[{name!r}]") for name, digest in self.payload_sha256.items()}
        object.__setattr__(self, "payload_sha256", MappingProxyType(copied))

    @classmethod
    def build(
        cls,
        *,
        finished: datetime.datetime,
        commit: str,
        input_sha256: str,
        result_sha256: str,
        payload_sha256: Mapping[str, str],
    ) -> ValidatedEnvelope:
        """Validate and build a consumed envelope summary."""
        return cls(
            finished=finished,
            commit=GitObjectId(commit),
            input_sha256=Sha256Hex(input_sha256),
            result_sha256=Sha256Hex(result_sha256),
            payload_sha256={name: Sha256Hex(digest) for name, digest in payload_sha256.items()},
        )


def _git(repository: pathlib.Path, *arguments: str) -> bytes:
    command = ["git", "-C", str(repository), *arguments]
    try:
        completed = subprocess.run(command, check=False, capture_output=True)
    except OSError as exc:
        raise MeasurementGitError(f"Git invocation failed in {repository}: {' '.join(arguments)}") from exc
    if completed.returncode != 0:
        detail = completed.stderr.decode("utf-8", errors="replace").strip()
        raise MeasurementGitError(f"Git failed in {repository}: {' '.join(arguments)}: {detail}")
    return completed.stdout


def _resolve_repository(start: pathlib.Path) -> pathlib.Path:
    probe = start.resolve()
    while not probe.exists() and probe.parent != probe:
        probe = probe.parent
    output = _git(probe, "rev-parse", "--show-toplevel")
    try:
        repository = pathlib.Path(output.decode("utf-8").strip()).resolve()
    except UnicodeDecodeError as exc:
        raise MeasurementGitError(f"Git toplevel is not UTF-8 in {probe}") from exc
    return _validate_repository_path(repository)


def _status(repository: pathlib.Path, *, owned_stage: Optional[pathlib.Path] = None) -> bytes:
    arguments = ["status", "--porcelain=v1", "-z", "--untracked-files=all", "--", "."]
    if owned_stage is not None:
        stage = _require_descendant(owned_stage, repository, field="owned_stage")
        relative = stage.relative_to(repository).as_posix()
        arguments.append(f":(exclude,top){relative}")
    return _git(repository, *arguments)


def _require_descendant(path: pathlib.Path, parent: pathlib.Path, *, field: str) -> pathlib.Path:
    resolved = path.resolve()
    root = parent.resolve()
    try:
        relative = resolved.relative_to(root)
    except ValueError as exc:
        raise InvalidMeasurementEnvelopeError(f"{field} escapes {root}: {resolved}") from exc
    if relative == pathlib.Path("."):
        raise InvalidMeasurementEnvelopeError(f"{field} must be below {root}")
    return resolved


def capture_clean_source(repository: pathlib.Path) -> SourceSnapshot:
    """Capture full HEAD and committed lock identity from a clean repository."""
    requested = _validate_repository_path(repository)
    resolved = _resolve_repository(requested)
    if resolved != requested:
        raise MeasurementGitError(f"measurement repository must be the exact Git toplevel: {requested}")
    if _status(resolved):
        raise DirtyMeasurementTreeError(f"measurement repository is dirty: {resolved}")
    commit_text = _git(resolved, "rev-parse", "HEAD^{commit}").decode("ascii").strip()
    commit = _validate_object_id(commit_text)
    lock_bytes = _git(resolved, "show", f"{commit}:pixi.lock")
    return SourceSnapshot.build(repository=resolved, commit=commit, pixi_lock_sha256=hashlib.sha256(lock_bytes).hexdigest())


def require_source_unchanged(source: SourceSnapshot, *, owned_stage: pathlib.Path) -> None:
    """Require unchanged HEAD and cleanliness excluding one exact owned stage."""
    stage = _require_descendant(owned_stage, source.repository, field="owned_stage")
    current = _validate_object_id(_git(source.repository, "rev-parse", "HEAD^{commit}").decode("ascii").strip())
    if current != source.commit:
        raise MeasurementInputChangedError(f"measurement HEAD changed from {source.commit} to {current}")
    if _status(source.repository, owned_stage=stage):
        raise DirtyMeasurementTreeError(f"measurement repository changed outside owned stage: {source.repository}")


@contextlib.contextmanager
def publication_stage(results_root: pathlib.Path, logical_name: str) -> Iterator[pathlib.Path]:
    """Yield one hidden sibling stage and clean it on ordinary failure."""
    if not logical_name or pathlib.PurePath(logical_name).name != logical_name or logical_name.startswith("."):
        raise InvalidMeasurementEnvelopeError(f"invalid logical artifact name: {logical_name!r}")
    if results_root.is_symlink():
        raise InvalidMeasurementEnvelopeError(f"results_root must not be a symlink: {results_root}")
    root = results_root.resolve()
    try:
        repository = _resolve_repository(root.parent)
    except MeasurementGitError as exc:
        raise InvalidMeasurementEnvelopeError(f"results_root is not inside a Git worktree: {root}") from exc
    _require_descendant(root, repository, field="results_root")
    root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f".{logical_name}.stage-", dir=root) as temporary:
        yield pathlib.Path(temporary).resolve()


def _canonical_json_bytes(value: object) -> bytes:
    try:
        encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
    except (TypeError, ValueError) as exc:
        raise InvalidMeasurementEnvelopeError("identity payload is not strict JSON") from exc
    return encoded


def _reject_json_constant(value: str) -> object:
    raise ValueError(f"non-standard JSON constant {value}")


def _decode_strict_json_bytes(data: bytes, *, field: str) -> object:
    try:
        return json.loads(data.decode("utf-8"), parse_constant=_reject_json_constant)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise InvalidMeasurementEnvelopeError(f"{field} is not strict UTF-8 JSON") from exc


def _sha256(data: bytes) -> Sha256Hex:
    return Sha256Hex(hashlib.sha256(data).hexdigest())


def _identity(version: str, fields: Mapping[str, object]) -> dict[str, object]:
    if type(version) is not str or not version:
        raise InvalidMeasurementEnvelopeError("identity version must be a non-empty string")
    body: dict[str, object] = {"version": version, **fields}
    return {**body, "sha256": _sha256(_canonical_json_bytes(body))}


def _isoformat_utc(value: datetime.datetime, *, field: str) -> str:
    return _validate_utc(value, field=field).isoformat(timespec="microseconds")


def _payload_digest_map(payloads: Mapping[str, bytes]) -> dict[str, dict[str, Sha256Hex]]:
    if not payloads:
        raise InvalidMeasurementEnvelopeError("result identity requires at least one payload")
    result: dict[str, dict[str, Sha256Hex]] = {}
    for name in sorted(payloads):
        payload = payloads[name]
        if type(name) is not str or name in ("", ".", "..") or pathlib.PurePosixPath(name).name != name or name == STAMP_NAME:
            raise InvalidMeasurementEnvelopeError(f"payload name must be one plain relative file name: {name!r}")
        if type(payload) is not bytes:
            raise InvalidMeasurementEnvelopeError(f"payload {name!r} must be exact bytes")
        result[name] = {"sha256": _sha256(payload)}
    return result


def build_envelope(
    *,
    artifact_kind: ArtifactKind,
    source: SourceSnapshot,
    started: datetime.datetime,
    finished: datetime.datetime,
    argv: Sequence[str],
    input_version: IdentityVersion,
    input_payload: Mapping[str, object],
    result_version: IdentityVersion,
    payloads: Mapping[str, bytes],
) -> dict[str, object]:
    """Build a canonical identity envelope over source, input, and payload bytes."""
    started_text = _isoformat_utc(started, field="started")
    finished_text = _isoformat_utc(finished, field="finished")
    if finished < started:
        raise InvalidMeasurementEnvelopeError("finished precedes started")
    if type(artifact_kind) is not str or not artifact_kind:
        raise InvalidMeasurementEnvelopeError("artifact_kind must be a non-empty string")
    if not argv or any(type(argument) is not str or not argument for argument in argv):
        raise InvalidMeasurementEnvelopeError("argv must contain non-empty strings")
    canonical_input = _decode_strict_json_bytes(_canonical_json_bytes(dict(input_payload)), field="input payload")
    if type(canonical_input) is not dict:
        raise InvalidMeasurementEnvelopeError("input payload must be a JSON object")
    payload_digests = _payload_digest_map(payloads)
    return {
        "envelope_version": ENVELOPE_VERSION,
        "artifact_kind": artifact_kind,
        "commit": source.commit,
        "dirty": False,
        "python": sys.version,
        "platform": platform.platform(),
        "started": started_text,
        "finished": finished_text,
        "argv": list(argv),
        "build_identity": _identity(
            BUILD_IDENTITY_VERSION,
            {"commit": source.commit, "pixi_lock_sha256": source.pixi_lock_sha256},
        ),
        "input_identity": _identity(input_version, {"payload": canonical_input}),
        "result_identity": _identity(result_version, {"payloads": payload_digests}),
    }


def write_envelope(stage: pathlib.Path, envelope: Mapping[str, object]) -> pathlib.Path:
    """Write one UTF-8 envelope into an existing stage."""
    if not stage.is_dir() or stage.is_symlink():
        raise InvalidMeasurementEnvelopeError(f"stage is not a real directory: {stage}")
    destination = stage / STAMP_NAME
    destination.write_text(json.dumps(dict(envelope), indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return destination


def _exact_dict(value: object, keys: Sequence[str], *, field: str, ordered: bool = False) -> dict[str, object]:
    if type(value) is not dict:
        raise InvalidMeasurementEnvelopeError(f"{field} must be an exact JSON object")
    actual = tuple(value)
    if (ordered and actual != tuple(keys)) or (not ordered and set(actual) != set(keys)):
        raise InvalidMeasurementEnvelopeError(f"{field} keys are {actual}, expected {tuple(keys)}")
    return value


def _parse_utc(value: object, *, field: str) -> datetime.datetime:
    if type(value) is not str:
        raise InvalidMeasurementEnvelopeError(f"{field} must be a UTC timestamp string")
    try:
        parsed = datetime.datetime.fromisoformat(value)
    except ValueError as exc:
        raise InvalidMeasurementEnvelopeError(f"{field} is not ISO-8601") from exc
    _validate_utc(parsed, field=field)
    if parsed.isoformat(timespec="microseconds") != value:
        raise InvalidMeasurementEnvelopeError(f"{field} must include canonical UTC microseconds")
    return parsed


def _read_regular_file(path: pathlib.Path, *, field: str) -> bytes:
    if not path.is_file() or path.is_symlink():
        raise InvalidMeasurementEnvelopeError(f"{field} must be a regular non-symlink file: {path}")
    return path.read_bytes()


def _validate_envelope_impl(
    result: pathlib.Path,
    *,
    logical_name: str,
    artifact_kind: ArtifactKind,
    repository: pathlib.Path,
) -> ValidatedEnvelope:
    """Reconstruct every common identity and exact payload-byte digest."""
    resolved_repository = _validate_repository_path(repository)
    if _resolve_repository(resolved_repository) != resolved_repository:
        raise InvalidMeasurementEnvelopeError(f"repository is not the exact Git toplevel: {resolved_repository}")
    if result.is_symlink():
        raise InvalidMeasurementEnvelopeError(f"result must not be a symlink: {result}")
    resolved_result = _require_descendant(result, resolved_repository, field="result")
    if not resolved_result.is_dir() or resolved_result.is_symlink():
        raise InvalidMeasurementEnvelopeError(f"result is not a real directory: {resolved_result}")
    if type(logical_name) is not str or logical_name in ("", ".", "..") or pathlib.PurePath(logical_name).name != logical_name:
        raise InvalidMeasurementEnvelopeError("logical_name must be one plain non-empty name")
    stamp = _decode_strict_json_bytes(_read_regular_file(resolved_result / STAMP_NAME, field="stamp"), field="stamp")
    envelope = _exact_dict(stamp, ENVELOPE_KEYS, field="stamp", ordered=True)
    if envelope["envelope_version"] != ENVELOPE_VERSION:
        raise InvalidMeasurementEnvelopeError("unknown envelope_version")
    if envelope["artifact_kind"] != artifact_kind:
        raise InvalidMeasurementEnvelopeError("artifact_kind does not match consumer")
    commit = _validate_object_id(envelope["commit"])
    if envelope["dirty"] is not False:
        raise InvalidMeasurementEnvelopeError("dirty must be exactly false")
    for field in ("python", "platform"):
        if type(envelope[field]) is not str or not envelope[field]:
            raise InvalidMeasurementEnvelopeError(f"{field} must be a non-empty string")
    started = _parse_utc(envelope["started"], field="started")
    finished = _parse_utc(envelope["finished"], field="finished")
    if finished < started:
        raise InvalidMeasurementEnvelopeError("finished precedes started")
    argv = envelope["argv"]
    if type(argv) is not list or not argv or any(type(argument) is not str or not argument for argument in argv):
        raise InvalidMeasurementEnvelopeError("argv must be a non-empty exact string list")

    build = _exact_dict(envelope["build_identity"], ("version", "commit", "pixi_lock_sha256", "sha256"), field="build_identity", ordered=True)
    if build["version"] != BUILD_IDENTITY_VERSION or build["commit"] != commit:
        raise InvalidMeasurementEnvelopeError("build identity version or commit mismatch")
    stamped_lock = _validate_sha256(build["pixi_lock_sha256"], field="build_identity.pixi_lock_sha256")
    try:
        committed_lock = _git(resolved_repository, "show", f"{commit}:pixi.lock")
    except MeasurementGitError as exc:
        raise InvalidMeasurementEnvelopeError(f"stamped commit cannot provide pixi.lock: {commit}") from exc
    if _sha256(committed_lock) != stamped_lock:
        raise InvalidMeasurementEnvelopeError("committed pixi.lock digest mismatch")
    expected_build = _identity(BUILD_IDENTITY_VERSION, {"commit": commit, "pixi_lock_sha256": stamped_lock})
    if build != expected_build:
        raise InvalidMeasurementEnvelopeError("build identity digest mismatch")

    input_identity = _exact_dict(envelope["input_identity"], ("version", "payload", "sha256"), field="input_identity", ordered=True)
    input_version = input_identity["version"]
    if type(input_version) is not str or not input_version or type(input_identity["payload"]) is not dict:
        raise InvalidMeasurementEnvelopeError("input identity version/payload is malformed")
    canonical_input = _decode_strict_json_bytes(_canonical_json_bytes(input_identity["payload"]), field="input payload")
    expected_input = _identity(input_version, {"payload": canonical_input})
    if input_identity != expected_input:
        raise InvalidMeasurementEnvelopeError("input identity digest mismatch")

    result_identity = _exact_dict(envelope["result_identity"], ("version", "payloads", "sha256"), field="result_identity", ordered=True)
    result_version = result_identity["version"]
    payload_entries = result_identity["payloads"]
    if type(result_version) is not str or not result_version or type(payload_entries) is not dict or not payload_entries:
        raise InvalidMeasurementEnvelopeError("result identity version/payloads is malformed")
    expected_names = {STAMP_NAME}
    payload_sha256: dict[str, str] = {}
    for name, entry_value in payload_entries.items():
        if type(name) is not str or pathlib.PurePosixPath(name).name != name or name == STAMP_NAME:
            raise InvalidMeasurementEnvelopeError(f"invalid result payload name: {name!r}")
        entry = _exact_dict(entry_value, ("sha256",), field=f"result_identity.payloads[{name!r}]", ordered=True)
        digest = _validate_sha256(entry["sha256"], field=f"payload {name!r} digest")
        actual = _sha256(_read_regular_file(resolved_result / name, field=f"payload {name!r}"))
        if actual != digest:
            raise InvalidMeasurementEnvelopeError(f"payload digest mismatch: {name}")
        expected_names.add(name)
        payload_sha256[name] = digest
    actual_names = {child.name for child in resolved_result.iterdir()}
    if actual_names != expected_names:
        raise InvalidMeasurementEnvelopeError(f"payload set mismatch: {sorted(actual_names)} != {sorted(expected_names)}")
    expected_result = _identity(result_version, {"payloads": {name: {"sha256": payload_sha256[name]} for name in sorted(payload_sha256)}})
    if result_identity != expected_result:
        raise InvalidMeasurementEnvelopeError("result identity digest mismatch")
    return ValidatedEnvelope.build(
        finished=finished,
        commit=commit,
        input_sha256=str(input_identity["sha256"]),
        result_sha256=str(result_identity["sha256"]),
        payload_sha256=payload_sha256,
    )


def validate_envelope(
    result: pathlib.Path,
    *,
    logical_name: str,
    artifact_kind: ArtifactKind,
    repository: pathlib.Path,
) -> ValidatedEnvelope:
    """Reconstruct identities, normalizing failures at the artifact boundary."""
    try:
        return _validate_envelope_impl(
            result,
            logical_name=logical_name,
            artifact_kind=artifact_kind,
            repository=repository,
        )
    except InvalidGitObjectIdError as exc:
        raise InvalidMeasurementEnvelopeError(f"{result.resolve()}: commit: {exc}") from exc
    except InvalidMeasurementEnvelopeError as exc:
        cause: BaseException = exc
        if isinstance(exc.__cause__, MeasurementGitError):
            cause = exc.__cause__
        raise InvalidMeasurementEnvelopeError(f"{result.resolve()}: {exc}") from cause


def publish_stage(*, source: SourceSnapshot, stage: pathlib.Path, final: pathlib.Path) -> pathlib.Path:
    """Publish one complete stage under the operational single-writer contract."""
    if stage.is_symlink():
        raise InvalidMeasurementEnvelopeError(f"stage must not be a symlink: {stage}")
    if final.exists() or final.is_symlink():
        raise MeasurementArtifactCollisionError(f"measurement result already exists: {final}")
    resolved_stage = _require_descendant(stage, source.repository, field="stage")
    resolved_final = _require_descendant(final, source.repository, field="final")
    if resolved_stage.parent != resolved_final.parent:
        raise InvalidMeasurementEnvelopeError("stage and final must be siblings")
    if not resolved_stage.name.startswith(f".{resolved_final.name}.stage-"):
        raise InvalidMeasurementEnvelopeError("stage logical name does not match final name")
    if not resolved_stage.is_dir() or resolved_stage.is_symlink():
        raise InvalidMeasurementEnvelopeError(f"stage is not a real directory: {resolved_stage}")
    require_source_unchanged(source, owned_stage=resolved_stage)
    if resolved_final.exists() or resolved_final.is_symlink():
        raise MeasurementArtifactCollisionError(f"measurement result appeared before publication: {resolved_final}")
    try:
        published = resolved_stage.rename(resolved_final)
    except FileExistsError as exc:
        raise MeasurementArtifactCollisionError(f"measurement result appeared during publication: {resolved_final}") from exc
    except OSError as exc:
        if exc.errno in (errno.EEXIST, errno.ENOTEMPTY):
            raise MeasurementArtifactCollisionError(f"measurement result appeared during publication: {resolved_final}") from exc
        raise
    return published.resolve()
