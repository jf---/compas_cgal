"""Git authentication for benchmark-claim history commits."""

from __future__ import annotations

import pathlib
import re
import shlex
import subprocess
from typing import Dict
from typing import Set

from tools.measurement_claim_benchmark_schema import BenchmarkClaimPayload
from tools.measurement_claim_errors import InvalidMeasurementClaimPayloadError

_OBJECT_ID = re.compile(r"(?:[0-9a-f]{40}|[0-9a-f]{64})\Z")


def _git(repository: pathlib.Path, *arguments: str) -> bytes:
    command = ["git", "--no-replace-objects", "-C", str(repository), *arguments]
    try:
        completed = subprocess.run(command, check=False, capture_output=True)
    except OSError as exc:
        raise InvalidMeasurementClaimPayloadError(f"benchmark history Git command could not start: {shlex.join(command)}") from exc
    if completed.returncode != 0:
        detail = completed.stderr.decode("utf-8", errors="replace").strip()
        raise InvalidMeasurementClaimPayloadError(f"benchmark history Git command failed: {shlex.join(command)}: {detail}")
    return completed.stdout


def _touched_paths(repository: pathlib.Path, commit: str) -> Set[str]:
    if _OBJECT_ID.fullmatch(commit) is None:
        raise InvalidMeasurementClaimPayloadError("benchmark history commit must be one full lowercase Git object ID")
    _git(repository, "cat-file", "commit", commit)
    raw = _git(repository, "diff-tree", "--root", "--no-commit-id", "--name-only", "-r", commit)
    try:
        return set(raw.decode("utf-8").splitlines())
    except UnicodeDecodeError as exc:
        raise InvalidMeasurementClaimPayloadError(f"benchmark history commit path list is not UTF-8: {commit}") from exc


def validate_benchmark_history(repository: pathlib.Path, payload: BenchmarkClaimPayload) -> None:
    """Require every recorded history commit to exist and touch its source file.

    Args:
        repository: Absolute repository root containing the artifact.
        payload: Validated benchmark claim payload.

    Raises:
        InvalidMeasurementClaimPayloadError: The repository, commit, claim source,
            or commit-to-source relationship is invalid.
    """
    if not repository.is_absolute() or repository != repository.resolve() or not repository.is_dir():
        raise InvalidMeasurementClaimPayloadError(f"benchmark history repository must be one absolute resolved directory: {repository}")
    paths_by_commit: Dict[str, Set[str]] = {}
    for claim in payload["claims"]:
        commit = str(claim["history_commit"])
        source, separator, line = claim["source"].rpartition(":")
        if not separator or not line.isdigit() or pathlib.PurePosixPath(source).as_posix() != source:
            raise InvalidMeasurementClaimPayloadError(f"{claim['claim_id']}: claim source must be one repository-relative path and line")
        if commit not in paths_by_commit:
            paths_by_commit[commit] = _touched_paths(repository, commit)
        touched = paths_by_commit[commit]
        if source not in touched:
            raise InvalidMeasurementClaimPayloadError(f"{claim['claim_id']}: history commit {commit} does not touch claim source {source}")
