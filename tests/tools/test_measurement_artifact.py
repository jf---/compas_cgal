from __future__ import annotations

import datetime
import hashlib
import importlib
import importlib.util
import json
import pathlib
import subprocess
from typing import Any
from typing import Optional

import pytest


UTC = datetime.timezone.utc
COMMIT_40 = "a" * 40
COMMIT_64 = "b" * 64
DIGEST = "c" * 64
REPLACEMENT_LOCK = b"replacement dependency graph\n"


def _module() -> Any:
    return importlib.import_module("tools.measurement_artifact")


def _git(repository: pathlib.Path, *arguments: str) -> bytes:
    return subprocess.run(
        ["git", "-C", str(repository), *arguments],
        check=True,
        capture_output=True,
    ).stdout


def _repository(tmp_path: pathlib.Path) -> pathlib.Path:
    repository = tmp_path / "repository"
    repository.mkdir(parents=True)
    _git(repository, "init", "-q")
    (repository / "pixi.lock").write_bytes(b"locked dependency graph\n")
    _git(repository, "add", "pixi.lock")
    _git(
        repository,
        "-c",
        "user.name=Jelle Feringa",
        "-c",
        "user.email=jelleferinga@gmail.com",
        "commit",
        "-qm",
        "lock",
    )
    return repository.resolve()


def _replace_head_commit(repository: pathlib.Path) -> str:
    original = _git(repository, "rev-parse", "HEAD^{commit}").decode("ascii").strip()
    (repository / "pixi.lock").write_bytes(REPLACEMENT_LOCK)
    _git(repository, "add", "pixi.lock")
    _git(
        repository,
        "-c",
        "user.name=Jelle Feringa",
        "-c",
        "user.email=jelleferinga@gmail.com",
        "commit",
        "-qm",
        "replacement lock",
    )
    replacement = _git(repository, "rev-parse", "HEAD^{commit}").decode("ascii").strip()
    _git(repository, "replace", original, replacement)
    _git(repository, "update-ref", "HEAD", original)
    return original


def _payloads() -> dict[str, bytes]:
    return {"benchmark_report.json": b"[]\n", "benchmark_report.md": b"# Benchmark corpus result\n"}


def _envelope(
    module: Any,
    source: Any,
    payloads: dict[str, bytes],
    *,
    finished: Optional[datetime.datetime] = None,
) -> dict[str, object]:
    started = datetime.datetime(2026, 8, 28, 12, 0, tzinfo=UTC)
    return module.build_envelope(
        artifact_kind=module.ArtifactKind("benchmark-corpus-result/v1"),
        source=source,
        started=started,
        finished=finished or started + datetime.timedelta(seconds=1),
        argv=("python", "-m", "benchmarks.cli"),
        input_version=module.IdentityVersion("benchmark-corpus-input/v1"),
        input_payload={"b": 2, "a": 1},
        result_version=module.IdentityVersion("benchmark-corpus-result/v1"),
        payloads=payloads,
    )


def _write_result(
    module: Any,
    repository: pathlib.Path,
    logical_name: str,
    *,
    payloads: Optional[dict[str, bytes]] = None,
) -> tuple[pathlib.Path, Any]:
    source = module.capture_clean_source(repository)
    result = repository / "benchmarks" / "results" / logical_name
    result.mkdir(parents=True)
    actual_payloads = payloads or _payloads()
    for name, payload in actual_payloads.items():
        (result / name).write_bytes(payload)
    module.write_envelope(result, _envelope(module, source, actual_payloads))
    return result, source


def test_measurement_artifact_module_exists() -> None:
    assert importlib.util.find_spec("tools.measurement_artifact") is not None


@pytest.mark.parametrize("commit", [COMMIT_40, COMMIT_64])
def test_source_snapshot_accepts_full_lowercase_git_object_ids(tmp_path: pathlib.Path, commit: str) -> None:
    module = _module()
    snapshot = module.SourceSnapshot.build(repository=tmp_path.resolve(), commit=commit, pixi_lock_sha256=DIGEST)
    assert snapshot.commit == commit


@pytest.mark.parametrize("commit", ["a" * 12, "A" * 40, "g" * 40, "a" * 41])
def test_source_snapshot_rejects_noncanonical_git_object_ids(tmp_path: pathlib.Path, commit: str) -> None:
    module = _module()
    with pytest.raises(module.InvalidGitObjectIdError):
        module.SourceSnapshot.build(repository=tmp_path.resolve(), commit=commit, pixi_lock_sha256=DIGEST)


def test_source_snapshot_rejects_relative_repository() -> None:
    module = _module()
    with pytest.raises(module.MeasurementArtifactError):
        module.SourceSnapshot.build(repository=pathlib.Path("relative"), commit=COMMIT_40, pixi_lock_sha256=DIGEST)


def test_validated_envelope_copies_payload_digest_mapping() -> None:
    module = _module()
    digests = {"result.json": DIGEST}
    validated = module.ValidatedEnvelope.build(
        finished=datetime.datetime(2026, 8, 28, tzinfo=UTC),
        commit=COMMIT_40,
        input_sha256=DIGEST,
        result_sha256=DIGEST,
        payload_sha256=digests,
    )
    digests["result.json"] = "d" * 64
    assert validated.payload_sha256["result.json"] == DIGEST
    with pytest.raises(TypeError):
        validated.payload_sha256["new"] = DIGEST


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("commit", "a" * 12),
        ("pixi_lock_sha256", "x" * 64),
    ],
)
def test_source_snapshot_raw_constructor_is_bypass_safe(tmp_path: pathlib.Path, field: str, value: str) -> None:
    module = _module()
    values = {"repository": tmp_path.resolve(), "commit": COMMIT_40, "pixi_lock_sha256": DIGEST}
    values[field] = value
    with pytest.raises(module.MeasurementArtifactError):
        module.SourceSnapshot(**values)


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("finished", datetime.datetime(2026, 8, 28)),
        ("commit", "a" * 12),
        ("input_sha256", "x" * 64),
        ("result_sha256", "x" * 64),
        ("payload_sha256", {"result": "x" * 64}),
    ],
)
def test_validated_envelope_raw_constructor_is_bypass_safe(field: str, value: object) -> None:
    module = _module()
    values: dict[str, object] = {
        "finished": datetime.datetime(2026, 8, 28, tzinfo=UTC),
        "commit": COMMIT_40,
        "input_sha256": DIGEST,
        "result_sha256": DIGEST,
        "payload_sha256": {"result": DIGEST},
    }
    values[field] = value
    with pytest.raises(module.MeasurementArtifactError):
        module.ValidatedEnvelope(**values)


def test_capture_clean_source_uses_full_head_and_committed_lock(tmp_path: pathlib.Path) -> None:
    module = _module()
    repository = _repository(tmp_path)
    source = module.capture_clean_source(repository)
    assert source.repository == repository
    assert source.commit == _git(repository, "rev-parse", "HEAD^{commit}").decode().strip()
    assert source.pixi_lock_sha256 == hashlib.sha256(b"locked dependency graph\n").hexdigest()


def test_capture_clean_source_does_not_bind_original_commit_to_replacement_tree(tmp_path: pathlib.Path) -> None:
    module = _module()
    repository = _repository(tmp_path)
    _replace_head_commit(repository)
    assert _git(repository, "status", "--porcelain=v1") == b""

    with pytest.raises(module.DirtyMeasurementTreeError):
        module.capture_clean_source(repository)


def test_capture_clean_source_rejects_untracked_input(tmp_path: pathlib.Path) -> None:
    module = _module()
    repository = _repository(tmp_path)
    (repository / "untracked.txt").write_text("input", encoding="utf-8")
    with pytest.raises(module.DirtyMeasurementTreeError):
        module.capture_clean_source(repository)


@pytest.mark.parametrize("damage", ["tracked", "staged"])
def test_capture_clean_source_rejects_tracked_and_staged_changes(tmp_path: pathlib.Path, damage: str) -> None:
    module = _module()
    repository = _repository(tmp_path)
    (repository / "pixi.lock").write_text("changed\n", encoding="utf-8")
    if damage == "staged":
        _git(repository, "add", "pixi.lock")
    with pytest.raises(module.DirtyMeasurementTreeError):
        module.capture_clean_source(repository)


def test_capture_clean_source_rejects_non_git_missing_lock_and_nested_path(tmp_path: pathlib.Path) -> None:
    module = _module()
    plain = tmp_path / "plain"
    plain.mkdir()
    with pytest.raises(module.MeasurementGitError):
        module.capture_clean_source(plain.resolve())

    repository = _repository(tmp_path / "missing")
    _git(repository, "rm", "pixi.lock")
    _git(repository, "-c", "user.name=Jelle Feringa", "-c", "user.email=jelleferinga@gmail.com", "commit", "-qm", "remove lock")
    with pytest.raises(module.MeasurementGitError):
        module.capture_clean_source(repository)

    nested = repository / "nested"
    nested.mkdir()
    with pytest.raises(module.MeasurementGitError):
        module.capture_clean_source(nested.resolve())


def test_build_envelope_separates_build_input_and_result_identities(tmp_path: pathlib.Path) -> None:
    module = _module()
    source = module.SourceSnapshot.build(repository=tmp_path.resolve(), commit=COMMIT_40, pixi_lock_sha256=DIGEST)
    first = _envelope(module, source, {"result": b"one"})
    reordered = module.build_envelope(
        artifact_kind=module.ArtifactKind("benchmark-corpus-result/v1"),
        source=source,
        started=datetime.datetime(2026, 8, 28, 12, 0, tzinfo=UTC),
        finished=datetime.datetime(2026, 8, 28, 12, 0, 1, tzinfo=UTC),
        argv=("python", "-m", "benchmarks.cli"),
        input_version=module.IdentityVersion("benchmark-corpus-input/v1"),
        input_payload={"a": 1, "b": 2},
        result_version=module.IdentityVersion("benchmark-corpus-result/v1"),
        payloads={"result": b"one"},
    )
    changed_result = _envelope(module, source, {"result": b"two"})
    assert first["build_identity"] == reordered["build_identity"]
    assert first["input_identity"] == reordered["input_identity"]
    assert first["result_identity"] == reordered["result_identity"]
    assert first["input_identity"] == changed_result["input_identity"]
    assert first["result_identity"] != changed_result["result_identity"]


def test_build_envelope_uses_specified_canonical_sha256_bytes(tmp_path: pathlib.Path) -> None:
    module = _module()
    source = module.SourceSnapshot.build(repository=tmp_path.resolve(), commit=COMMIT_40, pixi_lock_sha256=DIGEST)
    envelope = _envelope(module, source, {"z": b"last", "a": b"first"})
    build_body = {"version": "measurement-build-identity/v1", "commit": COMMIT_40, "pixi_lock_sha256": DIGEST}
    input_body = {"version": "benchmark-corpus-input/v1", "payload": {"a": 1, "b": 2}}
    result_body = {
        "version": "benchmark-corpus-result/v1",
        "payloads": {
            "a": {"sha256": hashlib.sha256(b"first").hexdigest()},
            "z": {"sha256": hashlib.sha256(b"last").hexdigest()},
        },
    }

    def digest(value: object) -> str:
        raw = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
        return hashlib.sha256(raw).hexdigest()

    assert envelope["build_identity"]["sha256"] == digest(build_body)
    assert envelope["input_identity"]["sha256"] == digest(input_body)
    assert envelope["result_identity"]["sha256"] == digest(result_body)


def test_precomputed_input_digest_equals_envelope_identity(tmp_path: pathlib.Path) -> None:
    module = _module()
    source = module.SourceSnapshot.build(repository=tmp_path.resolve(), commit=COMMIT_40, pixi_lock_sha256=DIGEST)
    payload = {"source": COMMIT_40, "semantic_command": ["-m", "benchmarks.cli", "figure6"]}
    expected = module.input_identity_sha256(version=module.IdentityVersion("benchmark-input/v1"), payload=payload)
    envelope = module.build_envelope(
        artifact_kind=module.ArtifactKind("benchmark/v1"),
        source=source,
        started=datetime.datetime(2026, 8, 28, 12, 0, tzinfo=UTC),
        finished=datetime.datetime(2026, 8, 28, 12, 0, 1, tzinfo=UTC),
        argv=("python", "-m", "benchmarks.cli", "figure6"),
        input_version=module.IdentityVersion("benchmark-input/v1"),
        input_payload=payload,
        result_version=module.IdentityVersion("benchmark-result/v1"),
        payloads={"result.json": b"{}\n"},
    )
    assert envelope["input_identity"]["sha256"] == expected


def test_build_envelope_rejects_nan_input(tmp_path: pathlib.Path) -> None:
    module = _module()
    source = module.SourceSnapshot.build(repository=tmp_path.resolve(), commit=COMMIT_40, pixi_lock_sha256=DIGEST)
    with pytest.raises(module.InvalidMeasurementEnvelopeError):
        module.build_envelope(
            artifact_kind=module.ArtifactKind("benchmark-corpus-result/v1"),
            source=source,
            started=datetime.datetime(2026, 8, 28, tzinfo=UTC),
            finished=datetime.datetime(2026, 8, 28, 0, 0, 1, tzinfo=UTC),
            argv=("python",),
            input_version=module.IdentityVersion("input/v1"),
            input_payload={"invalid": float("nan")},
            result_version=module.IdentityVersion("result/v1"),
            payloads={"result": b"payload"},
        )


@pytest.mark.parametrize("version", ["", 1, None])
def test_build_envelope_rejects_non_string_or_empty_identity_versions(tmp_path: pathlib.Path, version: object) -> None:
    module = _module()
    source = module.SourceSnapshot.build(repository=tmp_path.resolve(), commit=COMMIT_40, pixi_lock_sha256=DIGEST)
    with pytest.raises(module.InvalidMeasurementEnvelopeError):
        module.build_envelope(
            artifact_kind=module.ArtifactKind("benchmark-corpus-result/v1"),
            source=source,
            started=datetime.datetime(2026, 8, 28, tzinfo=UTC),
            finished=datetime.datetime(2026, 8, 28, 0, 0, 1, tzinfo=UTC),
            argv=("python",),
            input_version=version,
            input_payload={},
            result_version=module.IdentityVersion("result/v1"),
            payloads={"result": b"payload"},
        )


@pytest.mark.parametrize("name", ["", ".", "..", "nested/result", "stamp.json"])
def test_build_envelope_rejects_unsafe_payload_names(tmp_path: pathlib.Path, name: str) -> None:
    module = _module()
    source = module.SourceSnapshot.build(repository=tmp_path.resolve(), commit=COMMIT_40, pixi_lock_sha256=DIGEST)
    with pytest.raises(module.InvalidMeasurementEnvelopeError):
        _envelope(module, source, {name: b"payload"})


@pytest.mark.parametrize("mutation", ["raw", "missing", "extra"])
def test_validate_envelope_rejects_payload_set_or_byte_mutation(tmp_path: pathlib.Path, mutation: str) -> None:
    module = _module()
    repository = _repository(tmp_path)
    result, _ = _write_result(module, repository, "2026-08-28-" + _git(repository, "rev-parse", "--short=12", "HEAD").decode().strip())
    if mutation == "raw":
        (result / "benchmark_report.json").write_bytes(b"[1]\n")
    elif mutation == "missing":
        (result / "benchmark_report.json").unlink()
    else:
        (result / "extra.txt").write_text("extra", encoding="utf-8")
    with pytest.raises(module.InvalidMeasurementEnvelopeError):
        module.validate_envelope(
            result,
            logical_name=result.name,
            artifact_kind=module.ArtifactKind("benchmark-corpus-result/v1"),
            repository=repository,
        )


def test_validate_envelope_rejects_nonstandard_json_constant(tmp_path: pathlib.Path) -> None:
    module = _module()
    repository = _repository(tmp_path)
    result, _ = _write_result(module, repository, "2026-08-28-" + _git(repository, "rev-parse", "--short=12", "HEAD").decode().strip())
    stamp = (result / "stamp.json").read_text(encoding="utf-8").replace('"dirty": false', '"dirty": NaN')
    (result / "stamp.json").write_text(stamp, encoding="utf-8")
    with pytest.raises(module.InvalidMeasurementEnvelopeError):
        module.validate_envelope(
            result,
            logical_name=result.name,
            artifact_kind=module.ArtifactKind("benchmark-corpus-result/v1"),
            repository=repository,
        )


def _mutate_stamp(result: pathlib.Path, mutate: Any) -> None:
    stamp = json.loads((result / "stamp.json").read_text(encoding="utf-8"))
    mutate(stamp)
    (result / "stamp.json").write_text(json.dumps(stamp, indent=2, allow_nan=False) + "\n", encoding="utf-8")


@pytest.mark.parametrize(
    ("field", "mutate"),
    [
        ("keys", lambda stamp: stamp.update(extra=True)),
        ("envelope_version", lambda stamp: stamp.__setitem__("envelope_version", "unknown/v1")),
        ("artifact_kind", lambda stamp: stamp.__setitem__("artifact_kind", "wrong/v1")),
        ("dirty", lambda stamp: stamp.__setitem__("dirty", True)),
        ("commit", lambda stamp: stamp.__setitem__("commit", "a" * 12)),
        ("started", lambda stamp: stamp.__setitem__("started", "2026-08-28T12:00:00")),
        ("finished", lambda stamp: stamp.__setitem__("finished", "2026-08-28T11:00:00.000000+00:00")),
        ("argv", lambda stamp: stamp.__setitem__("argv", "python")),
        ("build.sha256", lambda stamp: stamp["build_identity"].__setitem__("sha256", "d" * 64)),
        ("input.version", lambda stamp: stamp["input_identity"].__setitem__("version", "")),
        ("input.sha256", lambda stamp: stamp["input_identity"].__setitem__("sha256", "d" * 64)),
        ("result.version", lambda stamp: stamp["result_identity"].__setitem__("version", "")),
        ("result.sha256", lambda stamp: stamp["result_identity"].__setitem__("sha256", "d" * 64)),
        ("payload.sha256", lambda stamp: stamp["result_identity"]["payloads"]["benchmark_report.json"].__setitem__("sha256", "x" * 64)),
    ],
)
def test_validate_envelope_rejects_adversarial_stamp_grammar(
    tmp_path: pathlib.Path,
    field: str,
    mutate: Any,
) -> None:
    module = _module()
    repository = _repository(tmp_path)
    name = "2026-08-28-" + _git(repository, "rev-parse", "--short=12", "HEAD").decode().strip()
    result, _ = _write_result(module, repository, name)
    _mutate_stamp(result, mutate)
    with pytest.raises(module.InvalidMeasurementEnvelopeError) as caught:
        module.validate_envelope(result, logical_name=name, artifact_kind=module.ArtifactKind("benchmark-corpus-result/v1"), repository=repository)
    assert str(result) in str(caught.value)
    assert field.split(".")[0] in str(caught.value)


def test_validate_envelope_normalizes_unavailable_stamped_commit_with_cause(tmp_path: pathlib.Path) -> None:
    module = _module()
    repository = _repository(tmp_path)
    name = "2026-08-28-" + _git(repository, "rev-parse", "--short=12", "HEAD").decode().strip()
    result, _ = _write_result(module, repository, name)
    unavailable = "f" * 40

    def mutate(stamp: dict[str, object]) -> None:
        stamp["commit"] = unavailable
        build = stamp["build_identity"]
        build["commit"] = unavailable
        body = {"version": build["version"], "commit": unavailable, "pixi_lock_sha256": build["pixi_lock_sha256"]}
        build["sha256"] = hashlib.sha256(json.dumps(body, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()

    _mutate_stamp(result, mutate)
    with pytest.raises(module.InvalidMeasurementEnvelopeError) as caught:
        module.validate_envelope(result, logical_name=name, artifact_kind=module.ArtifactKind("benchmark-corpus-result/v1"), repository=repository)
    assert isinstance(caught.value.__cause__, module.MeasurementGitError)


def test_validate_envelope_rejects_committed_lock_digest_mismatch(tmp_path: pathlib.Path) -> None:
    module = _module()
    repository = _repository(tmp_path)
    name = "2026-08-28-" + _git(repository, "rev-parse", "--short=12", "HEAD").decode().strip()
    result, _ = _write_result(module, repository, name)

    def mutate(stamp: dict[str, object]) -> None:
        build = stamp["build_identity"]
        build["pixi_lock_sha256"] = "d" * 64
        body = {"version": build["version"], "commit": build["commit"], "pixi_lock_sha256": build["pixi_lock_sha256"]}
        build["sha256"] = hashlib.sha256(json.dumps(body, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()

    _mutate_stamp(result, mutate)
    with pytest.raises(module.InvalidMeasurementEnvelopeError):
        module.validate_envelope(result, logical_name=name, artifact_kind=module.ArtifactKind("benchmark-corpus-result/v1"), repository=repository)


def test_validate_envelope_rejects_lock_content_from_replacement_commit(tmp_path: pathlib.Path) -> None:
    module = _module()
    repository = _repository(tmp_path)
    original = _replace_head_commit(repository)
    source = module.SourceSnapshot.build(
        repository=repository,
        commit=original,
        pixi_lock_sha256=hashlib.sha256(REPLACEMENT_LOCK).hexdigest(),
    )
    name = f"2026-08-28-{original[:12]}"
    result = repository / "benchmarks" / "results" / name
    result.mkdir(parents=True)
    payloads = _payloads()
    for payload_name, payload in payloads.items():
        (result / payload_name).write_bytes(payload)
    module.write_envelope(result, _envelope(module, source, payloads))

    with pytest.raises(module.InvalidMeasurementEnvelopeError):
        module.validate_envelope(
            result,
            logical_name=name,
            artifact_kind=module.ArtifactKind("benchmark-corpus-result/v1"),
            repository=repository,
        )


def test_require_source_unchanged_excludes_only_owned_stage(tmp_path: pathlib.Path) -> None:
    module = _module()
    repository = _repository(tmp_path)
    source = module.capture_clean_source(repository)
    stage = repository / "benchmarks" / "results" / ".owned.stage-1"
    stage.mkdir(parents=True)
    (stage / "payload").write_text("owned", encoding="utf-8")
    module.require_source_unchanged(source, owned_stage=stage)
    (repository / "sibling.txt").write_text("not owned", encoding="utf-8")
    with pytest.raises(module.DirtyMeasurementTreeError):
        module.require_source_unchanged(source, owned_stage=stage)


def test_require_source_unchanged_rejects_competing_stage_under_result_root(tmp_path: pathlib.Path) -> None:
    module = _module()
    repository = _repository(tmp_path)
    source = module.capture_clean_source(repository)
    owned = repository / "benchmarks" / "results" / ".owned.stage-1"
    competing = repository / "benchmarks" / "results" / ".competing.stage-1"
    owned.mkdir(parents=True)
    competing.mkdir()
    (competing / "payload").write_text("competing", encoding="utf-8")
    with pytest.raises(module.DirtyMeasurementTreeError):
        module.require_source_unchanged(source, owned_stage=owned)


def test_require_source_unchanged_rejects_changed_head(tmp_path: pathlib.Path) -> None:
    module = _module()
    repository = _repository(tmp_path)
    source = module.capture_clean_source(repository)
    (repository / "tracked.txt").write_text("next", encoding="utf-8")
    _git(repository, "add", "tracked.txt")
    _git(repository, "-c", "user.name=Jelle Feringa", "-c", "user.email=jelleferinga@gmail.com", "commit", "-qm", "next")
    stage = repository / "benchmarks" / "results" / ".owned.stage-1"
    stage.mkdir(parents=True)
    with pytest.raises(module.MeasurementInputChangedError):
        module.require_source_unchanged(source, owned_stage=stage)


def test_require_source_unchanged_rejects_stage_outside_repository(tmp_path: pathlib.Path) -> None:
    module = _module()
    repository = _repository(tmp_path)
    source = module.capture_clean_source(repository)
    outside = tmp_path / "outside"
    outside.mkdir()
    with pytest.raises(module.InvalidMeasurementEnvelopeError):
        module.require_source_unchanged(source, owned_stage=outside)


def test_publication_stage_is_hidden_and_cleans_after_failure(tmp_path: pathlib.Path) -> None:
    module = _module()
    repository = _repository(tmp_path)
    results = repository / "benchmarks" / "results"
    with pytest.raises(RuntimeError):
        with module.publication_stage(results, "2026-08-28-abcdef012345") as stage:
            assert stage.parent == results.resolve()
            assert stage.name.startswith(".2026-08-28-abcdef012345.stage-")
            raise RuntimeError("producer failed")
    assert not results.exists() or list(results.iterdir()) == []


@pytest.mark.parametrize("logical_name", ["../escape", "nested/name", "..", ".hidden"])
def test_publication_stage_rejects_unsafe_logical_names(tmp_path: pathlib.Path, logical_name: str) -> None:
    module = _module()
    repository = _repository(tmp_path)
    with pytest.raises(module.InvalidMeasurementEnvelopeError):
        with module.publication_stage(repository / "benchmarks" / "results", logical_name):
            pass


def test_publication_stage_rejects_root_outside_or_symlinked_outside_repository(tmp_path: pathlib.Path) -> None:
    module = _module()
    repository = _repository(tmp_path)
    outside = tmp_path / "outside"
    outside.mkdir()
    with pytest.raises(module.InvalidMeasurementEnvelopeError):
        with module.publication_stage(outside, "2026-08-28-abcdef012345"):
            pass
    linked = repository / "benchmarks" / "results"
    linked.parent.mkdir()
    linked.symlink_to(outside, target_is_directory=True)
    with pytest.raises(module.InvalidMeasurementEnvelopeError):
        with module.publication_stage(linked, "2026-08-28-abcdef012345"):
            pass


def test_publish_stage_renames_once_and_preserves_bundle(tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch) -> None:
    module = _module()
    repository = _repository(tmp_path)
    source = module.capture_clean_source(repository)
    results = repository / "benchmarks" / "results"
    rename_calls: list[tuple[pathlib.Path, pathlib.Path]] = []
    original_rename = pathlib.Path.rename

    def recording_rename(stage: pathlib.Path, final: pathlib.Path) -> pathlib.Path:
        rename_calls.append((stage, final))
        return original_rename(stage, final)

    monkeypatch.setattr(pathlib.Path, "rename", recording_rename)
    with module.publication_stage(results, "2026-08-28-abcdef012345") as stage:
        (stage / "payload").write_text("complete", encoding="utf-8")
        final = results / "2026-08-28-abcdef012345"
        assert module.publish_stage(source=source, stage=stage, final=final) == final.resolve()
    assert len(rename_calls) == 1
    assert (final / "payload").read_text(encoding="utf-8") == "complete"


def test_publish_stage_rejects_and_preserves_existing_nonempty_destination(tmp_path: pathlib.Path) -> None:
    module = _module()
    repository = _repository(tmp_path)
    source = module.capture_clean_source(repository)
    results = repository / "benchmarks" / "results"
    with module.publication_stage(results, "2026-08-28-abcdef012345") as stage:
        (stage / "payload").write_text("new", encoding="utf-8")
        final = results / "2026-08-28-abcdef012345"
        final.mkdir()
        (final / "payload").write_text("existing", encoding="utf-8")
        with pytest.raises(module.MeasurementArtifactCollisionError):
            module.publish_stage(source=source, stage=stage, final=final)
        assert (final / "payload").read_text(encoding="utf-8") == "existing"


def test_publish_stage_treats_existing_final_symlink_as_collision(tmp_path: pathlib.Path) -> None:
    module = _module()
    repository = _repository(tmp_path)
    source = module.capture_clean_source(repository)
    results = repository / "benchmarks" / "results"
    outside = tmp_path / "outside"
    outside.mkdir()
    with module.publication_stage(results, "2026-08-28-abcdef012345") as stage:
        final = results / "2026-08-28-abcdef012345"
        final.symlink_to(outside, target_is_directory=True)
        with pytest.raises(module.MeasurementArtifactCollisionError):
            module.publish_stage(source=source, stage=stage, final=final)
        assert final.is_symlink()


def test_publish_stage_rejects_mismatched_stage_and_final_roots(tmp_path: pathlib.Path) -> None:
    module = _module()
    repository = _repository(tmp_path)
    source = module.capture_clean_source(repository)
    with module.publication_stage(repository / "benchmarks" / "results", "2026-08-28-abcdef012345") as stage:
        outside_final = repository / "benchmarks" / "measurement_claim_results" / "2026-08-28-abcdef012345"
        with pytest.raises(module.InvalidMeasurementEnvelopeError):
            module.publish_stage(source=source, stage=stage, final=outside_final)


def test_publish_stage_rejects_final_name_that_does_not_match_stage(tmp_path: pathlib.Path) -> None:
    module = _module()
    repository = _repository(tmp_path)
    source = module.capture_clean_source(repository)
    results = repository / "benchmarks" / "results"
    with module.publication_stage(results, "2026-08-28-abcdef012345") as stage:
        with pytest.raises(module.InvalidMeasurementEnvelopeError):
            module.publish_stage(source=source, stage=stage, final=results / "2026-08-28-different0000")
