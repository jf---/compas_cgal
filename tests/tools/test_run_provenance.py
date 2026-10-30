"""Run provenance: what tree, tool, environment, machine and time a timing record came from, validated on load."""

from __future__ import annotations

import dataclasses
import hashlib
import subprocess
from datetime import datetime
from datetime import timezone
from pathlib import Path

import pytest

from tools.pipeline_pencil_cost_results import PIPELINE_PENCIL_COST_PATH
from tools.pipeline_pencil_cost_results import PIPELINE_PENCIL_COST_TABLE_PATH
from tools.pipeline_time_gate import PIPELINE_TIME_GATE_TABLE_PATH
from tools.pipeline_timing_1_1_base import PIPELINE_TIMING_1_1_SESSION_PATH
from tools.pipeline_timing_results import PIPELINE_TIMING_PATH
from tools.pipeline_timing_results import PIPELINE_TIMING_TABLE_PATH
from tools.run_provenance import GATE_OUTPUT_RECORDS
from tools.run_provenance import MalformedProvenanceError
from tools.run_provenance import RunProvenance
from tools.run_provenance import TreeChangedDuringRunError
from tools.run_provenance import file_set_sha256
from tools.run_provenance import interpreter_identity_here
from tools.run_provenance import machine_identity
from tools.run_provenance import require_unchanged_tree
from tools.run_provenance import tree_identity
from tools.run_provenance import utc_now

ROOT = Path(__file__).resolve().parents[2]
REV = "7b001d473d6e94b79fc8ff3cb852ecde1ed67733"
SHA = "a" * 64


def provenance(**overrides: object) -> RunProvenance:
    default = RunProvenance(
        rev=REV,
        dirty=False,
        patch_sha256=None,
        timing_tool_sha256=SHA,
        pixi_lock_sha256=SHA,
        environment="default",
        python_version="3.12.13",
        stock_2_sha256=SHA,
        native_sha256=SHA,
        native_sources_sha256=SHA,
        host="host.local",
        cpu_model="Apple M1 Max",
        started_utc="2026-10-22T18:00:00+00:00",
        ended_utc="2026-10-22T18:10:00+00:00",
    )
    return dataclasses.replace(default, **overrides)


def test_provenance_round_trips_through_its_json() -> None:
    record = provenance()
    assert RunProvenance.build(record.to_json()) == record


def test_a_dirty_record_carries_its_patch_digest() -> None:
    record = provenance(dirty=True, patch_sha256=SHA)
    assert RunProvenance.build(record.to_json()) == record


@pytest.mark.parametrize(
    "overrides",
    [
        {"dirty": True, "patch_sha256": None},
        {"dirty": False, "patch_sha256": SHA},
        {"rev": "7b001d47"},
        {"pixi_lock_sha256": "not-a-hash"},
        {"host": ""},
        {"started_utc": "2026-10-22 18:00"},
        {"started_utc": "2026-10-22T18:00:00"},
        {"ended_utc": "2026-10-22T17:00:00+00:00"},
    ],
)
def test_a_malformed_provenance_is_refused(overrides: dict[str, object]) -> None:
    with pytest.raises(MalformedProvenanceError):
        RunProvenance.build(provenance(**overrides).to_json())


def test_a_provenance_with_a_missing_field_is_refused() -> None:
    raw = provenance().to_json()
    del raw["host"]
    with pytest.raises(MalformedProvenanceError):
        RunProvenance.build(raw)


def test_the_utc_times_parse_to_aware_datetimes() -> None:
    record = provenance()
    assert record.started == datetime(2026, 10, 22, 18, 0, tzinfo=timezone.utc)
    assert (record.ended - record.started).total_seconds() == 600


def test_utc_now_is_an_aware_utc_iso_string() -> None:
    assert datetime.fromisoformat(utc_now()).utcoffset() is not None


def test_the_tool_digest_covers_paths_and_bytes(tmp_path: Path) -> None:
    first, second = tmp_path / "a.py", tmp_path / "b.py"
    first.write_text("one")
    second.write_text("two")
    digest = file_set_sha256(tmp_path, (Path("a.py"), Path("b.py")))
    assert digest == file_set_sha256(tmp_path, (Path("b.py"), Path("a.py")))
    second.write_text("three")
    assert file_set_sha256(tmp_path, (Path("a.py"), Path("b.py"))) != digest


def test_tree_identity_names_the_full_rev_and_digests_a_dirty_tree(tmp_path: Path) -> None:
    def git(*args: str) -> None:
        subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t", *args], cwd=tmp_path, check=True, capture_output=True)

    git("init", "-q")
    (tmp_path / "f.txt").write_text("one")
    git("add", "f.txt")
    git("commit", "-q", "-m", "one")
    rev = subprocess.run(["git", "rev-parse", "HEAD"], cwd=tmp_path, check=True, capture_output=True, text=True).stdout.strip()
    assert tree_identity(tmp_path) == (rev, False, None)
    (tmp_path / "f.txt").write_text("two")
    _rev, dirty, patch = tree_identity(tmp_path)
    assert dirty is True and patch is not None and len(patch) == len(SHA)
    (tmp_path / "new.txt").write_text("untracked")
    _rev, _dirty, with_untracked = tree_identity(tmp_path)
    assert with_untracked != patch


def test_the_excluded_records_are_exactly_the_gate_tools_default_outputs() -> None:
    assert set(GATE_OUTPUT_RECORDS) == {
        PIPELINE_TIMING_PATH,
        PIPELINE_TIMING_TABLE_PATH,
        PIPELINE_TIMING_1_1_SESSION_PATH,
        PIPELINE_PENCIL_COST_PATH,
        PIPELINE_PENCIL_COST_TABLE_PATH,
        PIPELINE_TIME_GATE_TABLE_PATH,
    }


def test_writing_the_gate_outputs_leaves_the_tree_identity_unchanged(tmp_path: Path) -> None:
    def git(*args: str) -> None:
        subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t", *args], cwd=tmp_path, check=True, capture_output=True)

    git("init", "-q")
    tracked = (PIPELINE_TIMING_PATH, PIPELINE_TIMING_TABLE_PATH, Path("tools") / "measured.py")
    for path in tracked:
        (tmp_path / path).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / path).write_text("committed")
    git("add", *(path.as_posix() for path in tracked))
    git("commit", "-q", "-m", "one")
    (tmp_path / "tools" / "measured.py").write_text("edited before the run")
    at_timing_end = tree_identity(tmp_path)
    assert at_timing_end[1] is True
    for record in GATE_OUTPUT_RECORDS:
        (tmp_path / record).write_text(f"written after the end read: {record}")
    assert tree_identity(tmp_path) == at_timing_end
    (tmp_path / "tools" / "measured.py").write_text("edited after the run")
    assert tree_identity(tmp_path) != at_timing_end
    (tmp_path / "tools" / "measured.py").write_text("edited before the run")
    (tmp_path / "tools" / "new.py").write_text("an untracked file outside the records")
    assert tree_identity(tmp_path) != at_timing_end


def test_a_tree_that_differs_only_in_the_gate_outputs_is_clean(tmp_path: Path) -> None:
    def git(*args: str) -> None:
        subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t", *args], cwd=tmp_path, check=True, capture_output=True)

    git("init", "-q")
    (tmp_path / PIPELINE_TIMING_PATH).parent.mkdir(parents=True)
    (tmp_path / PIPELINE_TIMING_PATH).write_text("committed")
    git("add", PIPELINE_TIMING_PATH.as_posix())
    git("commit", "-q", "-m", "one")
    rev = subprocess.run(["git", "rev-parse", "HEAD"], cwd=tmp_path, check=True, capture_output=True, text=True).stdout.strip()
    for record in GATE_OUTPUT_RECORDS:
        (tmp_path / record).write_text("rewritten")
    assert tree_identity(tmp_path) == (rev, False, None)


class _Repository:
    """A git repository at *root* with an export of its HEAD under its ignored ``build/`` (the landing export layout)."""

    def __init__(self, root: Path) -> None:
        self.root = root
        self.export = root / "build" / "land" / "tree"
        self.git("init", "-q")
        (root / ".gitignore").write_text("build/\n__pycache__/\n")
        (root / "tools").mkdir()
        (root / "tools" / "measured.py").write_text("committed")
        self.git("add", ".gitignore", "tools/measured.py")
        self.git("commit", "-q", "-m", "one")
        self.export.mkdir(parents=True)
        archive = subprocess.run(["git", "archive", "HEAD"], cwd=root, check=True, capture_output=True).stdout
        subprocess.run(["tar", "-x", "-C", str(self.export)], input=archive, check=True)

    def git(self, *args: str) -> None:
        subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t", *args], cwd=self.root, check=True, capture_output=True)

    def head(self) -> str:
        return subprocess.run(["git", "rev-parse", "HEAD"], cwd=self.root, check=True, capture_output=True, text=True).stdout.strip()


def test_an_export_is_identified_by_its_content_never_as_its_clean_parent(tmp_path: Path) -> None:
    repository = _Repository(tmp_path)
    rev, dirty, digest = tree_identity(repository.export)
    assert rev == repository.head() and dirty is True and digest is not None and len(digest) == len(SHA)
    assert tree_identity(repository.root) == (repository.head(), False, None)


def test_a_commit_in_the_containing_repository_leaves_an_unchanged_export_unchanged(tmp_path: Path) -> None:
    repository = _Repository(tmp_path)
    at_start = tree_identity(repository.export)
    (tmp_path / "other.txt").write_text("landed elsewhere")
    repository.git("add", "other.txt")
    repository.git("commit", "-q", "-m", "two")
    at_end = tree_identity(repository.export)
    assert at_end[0] == repository.head() != at_start[0] and at_end[2] == at_start[2]
    require_unchanged_tree(at_start, at_end, "the export")


def test_a_content_change_in_an_export_is_refused(tmp_path: Path) -> None:
    repository = _Repository(tmp_path)
    at_start = tree_identity(repository.export)
    (repository.export / "tools" / "measured.py").write_text("edited during the run")
    with pytest.raises(TreeChangedDuringRunError):
        require_unchanged_tree(at_start, tree_identity(repository.export), "the export")
    (repository.export / "tools" / "measured.py").write_text("committed")
    (repository.export / "tools" / "new.py").write_text("added during the run")
    with pytest.raises(TreeChangedDuringRunError):
        require_unchanged_tree(at_start, tree_identity(repository.export), "the export")


def test_an_export_ignores_its_gate_outputs_and_ignored_files(tmp_path: Path) -> None:
    repository = _Repository(tmp_path)
    at_start = tree_identity(repository.export)
    for record in GATE_OUTPUT_RECORDS:
        (repository.export / record).parent.mkdir(parents=True, exist_ok=True)
        (repository.export / record).write_text("written by the run")
    (repository.export / "tools" / "__pycache__").mkdir()
    (repository.export / "tools" / "__pycache__" / "measured.cpython-312.pyc").write_bytes(b"bytecode")
    assert tree_identity(repository.export) == at_start


def test_a_worktree_commit_still_changes_the_worktree_identity(tmp_path: Path) -> None:
    repository = _Repository(tmp_path)
    at_start = tree_identity(repository.root)
    (tmp_path / "other.txt").write_text("committed during the run")
    repository.git("add", "other.txt")
    repository.git("commit", "-q", "-m", "two")
    with pytest.raises(TreeChangedDuringRunError):
        require_unchanged_tree(at_start, tree_identity(repository.root), "the worktree")


def test_a_worktree_and_an_export_of_the_same_content_never_share_an_identity(tmp_path: Path) -> None:
    """Same files, two kinds of digest: a dirty worktree's diff digest and an export's content digest never compare equal."""
    repository = _Repository(tmp_path)
    (tmp_path / "tools" / "measured.py").write_text("edited")
    (repository.export / "tools" / "measured.py").write_text("edited")
    worktree, export = tree_identity(repository.root), tree_identity(repository.export)
    assert worktree[:2] == export[:2] and worktree[2] != export[2]
    with pytest.raises(TreeChangedDuringRunError):
        require_unchanged_tree(worktree, export, "mixed")


def test_this_machine_and_interpreter_are_identified() -> None:
    host, cpu = machine_identity()
    assert host and cpu
    environment, python_version, stock_2 = interpreter_identity_here()
    assert environment and python_version
    assert len(stock_2) == len(hashlib.sha256().hexdigest())
