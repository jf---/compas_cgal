from __future__ import annotations

import hashlib
import json
import pathlib
import subprocess
from typing import Any


REPOSITORY = pathlib.Path(__file__).resolve().parents[2]
ARTIFACT_DIRECTORY = REPOSITORY / "benchmarks/measurement_claim_history/2026-08-29-53135e04390e"
MANIFEST = ARTIFACT_DIRECTORY / "manifest.json"
PARENT_COMMIT = "b531d215e6ca06741d040b070ba43164f61abd58"
CORRECTION_COMMIT = "53135e04390e84bf69aa74dc4d0c1ce6ca308eb4"
RADIAL_SOURCE = "src/compas_cgal/engagement_radial_toolpath.py"
ADVANCE_SOURCE = "src/compas_cgal/engagement_toolpath.py"
SOURCE_PATHS = [RADIAL_SOURCE, ADVANCE_SOURCE]
MANIFEST_KEYS = {
    "schema_version",
    "status",
    "parent_commit",
    "correction_commit",
    "source_paths",
    "patch_file",
    "patch_sha256",
    "claim_sources",
}


class DuplicateManifestKeyError(RuntimeError):
    """A JSON object repeats a key and therefore has ambiguous provenance."""


def _reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise DuplicateManifestKeyError(f"duplicate manifest key: {key}")
        result[key] = value
    return result


def _load_manifest() -> dict[str, Any]:
    loaded = json.loads(MANIFEST.read_text(encoding="utf-8"), object_pairs_hook=_reject_duplicate_keys)
    assert isinstance(loaded, dict)
    return loaded


def _git(*arguments: str) -> bytes:
    return subprocess.run(
        ["git", "--no-replace-objects", *arguments],
        cwd=REPOSITORY,
        check=True,
        capture_output=True,
    ).stdout


def _git_object_hexdigest(object_id: str, canonical_object: bytes) -> str:
    if len(object_id) == 40:
        return hashlib.sha1(canonical_object).hexdigest()
    if len(object_id) == 64:
        return hashlib.sha256(canonical_object).hexdigest()
    raise AssertionError(f"unsupported Git object ID length: {len(object_id)}")


def _authenticated_commit(object_id: str) -> bytes:
    commit = _git("cat-file", "commit", object_id)
    canonical_object = f"commit {len(commit)}\0".encode() + commit
    assert _git_object_hexdigest(object_id, canonical_object) == object_id
    return commit


def test_historical_assertions_match_authenticated_correction_diff() -> None:
    manifest = _load_manifest()
    assert set(manifest) == MANIFEST_KEYS
    assert manifest["schema_version"] == "measurement-claim-history/v1"
    assert manifest["status"] == "superseded-source-assertions"
    assert manifest["parent_commit"] == PARENT_COMMIT
    assert manifest["correction_commit"] == CORRECTION_COMMIT
    assert manifest["patch_file"] == "source-correction.patch"

    correction = _authenticated_commit(CORRECTION_COMMIT)
    parent_headers = [line.removeprefix(b"parent ").decode("ascii") for line in correction.partition(b"\n\n")[0].splitlines() if line.startswith(b"parent ")]
    assert parent_headers == [PARENT_COMMIT]
    _authenticated_commit(PARENT_COMMIT)

    source_paths = manifest["source_paths"]
    assert source_paths == SOURCE_PATHS
    assert len(source_paths) == len(set(source_paths))

    patch = (ARTIFACT_DIRECTORY / manifest["patch_file"]).read_bytes()
    expected_patch = _git(
        "diff",
        "--binary",
        PARENT_COMMIT,
        CORRECTION_COMMIT,
        "--",
        *SOURCE_PATHS,
    )
    assert patch == expected_patch
    assert hashlib.sha256(patch).hexdigest() == manifest["patch_sha256"]

    claim_sources = manifest["claim_sources"]
    expected_claim_ids = [f"MC-{ordinal:03d}" for ordinal in range(1, 11)]
    assert [record["claim_id"] for record in claim_sources] == expected_claim_ids
    assert all(set(record) == {"claim_id", "source_path"} for record in claim_sources)
    assert len({record["claim_id"] for record in claim_sources}) == len(claim_sources)
    assert [record["source_path"] for record in claim_sources] == [RADIAL_SOURCE] * 8 + [ADVANCE_SOURCE] * 2
