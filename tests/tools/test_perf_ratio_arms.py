"""Arm setup and the per-arm proofs of the ratio harness, on fake arm trees only."""

from __future__ import annotations

import dataclasses
import os
import subprocess
import sys
from collections.abc import Mapping
from collections.abc import Sequence
from pathlib import Path

import pytest

from tests.tools.perf_ratio_fake_arms import FAKE_NATIVE
from tests.tools.perf_ratio_fake_arms import fake_arm
from tests.tools.perf_ratio_fake_arms import fake_arm_tree
from tests.tools.perf_ratio_fake_arms import fake_bracket
from tests.tools.perf_ratio_fake_arms import fake_interpreter
from tests.tools.perf_ratio_fake_arms import fake_origin
from tests.tools.perf_ratio_fake_arms import fake_tracker
from tests.tools.perf_ratio_fake_arms import git
from tools import perf_ratio_arms
from tools.held_native_identity import NativeBuildStaleError
from tools.perf_ratio_arms import PYPROJECT
from tools.perf_ratio_arms import Arm
from tools.perf_ratio_arms import ArmProbeFailedError
from tools.perf_ratio_arms import ArmProbeNoOutputError
from tools.perf_ratio_arms import HeadCheckoutFailedError
from tools.perf_ratio_arms import HeadRevUnknownError
from tools.perf_ratio_arms import SameArmError
from tools.perf_ratio_arms import arm_python
from tools.perf_ratio_arms import bracket_arm
from tools.perf_ratio_arms import checkout_head
from tools.perf_ratio_arms import interpreter_prefix_of
from tools.perf_ratio_arms import live_bracket
from tools.perf_ratio_arms import open_arm
from tools.perf_ratio_arms import rebuild_arm
from tools.perf_ratio_arms import require_distinct_arms
from tools.perf_ratio_arms import require_settled
from tools.perf_ratio_record import ArmName
from tools.perf_ratio_record import CpuSeconds
from tools.perf_ratio_workloads import WORKLOAD_SETS
from tools.perf_ratio_workloads import Family
from tools.perf_ratio_workloads import child_command
from tools.perf_ratio_workloads import parse_child
from tools.pipeline_timing_1_1_base import ForeignLibraryError
from tools.pipeline_timing_1_1_base import LibraryOrigin
from tools.pipeline_timing_1_1_base import MeasuredTreeError
from tools.pipeline_timing_1_1_base import Runner
from tools.pipeline_timing_1_1_base import base_child_env
from tools.run_native_identity import NativeFinisher
from tools.run_native_identity import NativeIdentity
from tools.run_native_identity import NativeTracker
from tools.run_provenance import TreeChangedDuringRunError
from tools.run_provenance import TreeIdentity
from tools.run_provenance import sha256_file
from tools.run_provenance import tree_identity

# CPU seconds and output digest the fake arms' stub children report.
FAKE_CPU_S = 1.5
FAKE_DIGEST = "fake-digest"


def _tree(path: Path) -> tuple[Path, str]:
    path.mkdir()
    return path, fake_arm_tree(path, cpu_seconds=FAKE_CPU_S, digest=FAKE_DIGEST)


def _two_arms(tmp_path: Path) -> tuple[Arm, Arm]:
    ref_root, ref_rev = _tree(tmp_path / "ref")
    head_root, head_rev = _tree(tmp_path / "head")
    return fake_arm(ArmName.REF, ref_root, ref_rev), fake_arm(ArmName.HEAD, head_root, head_rev)


def _head(root: Path) -> str:
    return git(root, "rev-parse", "HEAD")


def test_arm_python_names_the_arm_manifest(tmp_path: Path) -> None:
    python = arm_python(tmp_path)
    index = python.index("--manifest-path")
    assert python[index + 1] == str(tmp_path / "pyproject.toml")
    assert python[-1] == "python"


def test_arm_python_never_installs_into_the_arm_environment(tmp_path: Path) -> None:
    assert "--no-install" in arm_python(tmp_path)


def test_arms_sharing_a_build_dir_are_refused(tmp_path: Path) -> None:
    ref, head = _two_arms(tmp_path)
    with pytest.raises(SameArmError, match="editable build"):
        require_distinct_arms(ref, dataclasses.replace(head, editable_build_dir=ref.editable_build_dir))


def test_arms_sharing_an_interpreter_prefix_are_refused(tmp_path: Path) -> None:
    ref, head = _two_arms(tmp_path)
    with pytest.raises(SameArmError, match="interpreter prefix"):
        require_distinct_arms(ref, dataclasses.replace(head, interpreter_prefix=ref.interpreter_prefix))


def test_one_tree_as_both_arms_is_refused(tmp_path: Path) -> None:
    ref, head = _two_arms(tmp_path)
    with pytest.raises(SameArmError, match="root"):
        require_distinct_arms(ref, dataclasses.replace(head, root=ref.root))


def test_a_symlinked_alias_of_one_tree_is_refused(tmp_path: Path) -> None:
    ref, head = _two_arms(tmp_path)
    alias = tmp_path / "alias"
    alias.symlink_to(ref.root)
    with pytest.raises(SameArmError, match="root"):
        require_distinct_arms(ref, dataclasses.replace(head, root=alias))


def test_distinct_fake_arms_pass(tmp_path: Path) -> None:
    ref, head = _two_arms(tmp_path)
    assert require_distinct_arms(ref, head) is None


def test_a_foreign_library_is_refused(tmp_path: Path) -> None:
    # The test interpreter's compas_cgal resolves to this worktree, outside the fake tree.
    root, rev = _tree(tmp_path / "head")
    with pytest.raises(ForeignLibraryError):
        open_arm(ArmName.HEAD, root, (sys.executable,), rev)


def test_open_arm_probes_the_arm_with_its_own_python_tree_and_pinned_env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    root, rev = _tree(tmp_path / "head")
    python = ("arm-python", "--flag")
    build_dir = root / "build" / "editable"
    prefix = root / ".env"
    calls: list[tuple[str, tuple[str, ...], Path, dict[str, str]]] = []

    def build_dir_of(py: Sequence[str], cwd: Path, env: Mapping[str, str]) -> Path:
        calls.append(("build", tuple(py), cwd, dict(env)))
        return build_dir

    def origin_of(py: Sequence[str], cwd: Path, env: Mapping[str, str]) -> LibraryOrigin:
        calls.append(("origin", tuple(py), cwd, dict(env)))
        return fake_origin(root)

    def prefix_of(py: Sequence[str], cwd: Path, env: Mapping[str, str]) -> Path:
        calls.append(("prefix", tuple(py), cwd, dict(env)))
        return prefix

    monkeypatch.setenv("PYTHONPATH", "/orchestrator/tree")
    monkeypatch.setenv("SKBUILD_EDITABLE_SKIP", "/orchestrator/build")
    monkeypatch.setattr(perf_ratio_arms, "editable_build_dir_of", build_dir_of)
    monkeypatch.setattr(perf_ratio_arms, "library_origin_of", origin_of)
    monkeypatch.setattr(perf_ratio_arms, "interpreter_prefix_of", prefix_of)

    arm = open_arm(ArmName.HEAD, root, python, rev)

    assert [name for name, *_ in calls] == ["build", "origin", "prefix"]
    assert all(py == python and cwd == root for _, py, cwd, _ in calls)
    unpinned, pinned, prefix_env = (env for *_, env in calls)
    assert unpinned["PYTHONPATH"] == str(root) and "SKBUILD_EDITABLE_SKIP" not in unpinned
    assert pinned == prefix_env == base_child_env(root, os.environ, build_dir)
    assert pinned["SKBUILD_EDITABLE_SKIP"] == str(build_dir) and pinned["PYTHONPATH"] == str(root)
    assert arm == Arm(name=ArmName.HEAD, root=root, python=python, env=pinned, editable_build_dir=build_dir, interpreter_prefix=prefix, rev=rev)


def test_interpreter_prefix_of_reads_the_interpreters_own_prefix(tmp_path: Path) -> None:
    assert interpreter_prefix_of((sys.executable,), tmp_path, os.environ) == Path(sys.prefix).resolve()


def test_a_tree_changed_mid_run_is_refused(tmp_path: Path) -> None:
    root, rev = _tree(tmp_path / "head")
    arm = fake_arm(ArmName.HEAD, root, rev)
    bracket = fake_bracket(arm, trees=iter([TreeIdentity(rev, False, None), TreeIdentity(rev, True, "ab")]))
    with pytest.raises(TreeChangedDuringRunError):
        bracket()()


@pytest.mark.parametrize("moved", ["another-rev", "dirty"])
def test_a_bracket_started_off_the_arms_clean_rev_is_refused(tmp_path: Path, moved: str) -> None:
    # open_arm proved the tree clean at arm.rev; a tree that moved or got dirty since is refused at the start.
    root, rev = _tree(tmp_path / "head")
    arm = fake_arm(ArmName.HEAD, root, rev)
    start = TreeIdentity("f" * len(rev), False, None) if moved == "another-rev" else TreeIdentity(rev, True, "ab")
    with pytest.raises(MeasuredTreeError, match="head arm"):
        fake_bracket(arm, trees=iter([start, start]))()


def test_a_library_that_moves_mid_run_is_refused(tmp_path: Path) -> None:
    root, rev = _tree(tmp_path / "head")
    arm = fake_arm(ArmName.HEAD, root, rev)
    outside = tmp_path / "elsewhere"
    origins = iter([LibraryOrigin(outside / "compas_cgal" / "__init__.py", outside / "_stock_2.so")])
    trees = iter([TreeIdentity(rev, False, None)] * 2)
    bracket = bracket_arm(arm, read_tree=lambda: next(trees), track_native=fake_tracker(), origin_of=lambda: next(origins), interpreter_of=fake_interpreter)
    finish = bracket()
    with pytest.raises(ForeignLibraryError):
        finish()


def test_the_bracket_records_both_native_identities_and_the_lock(tmp_path: Path) -> None:
    root, rev = _tree(tmp_path / "head")
    arm = fake_arm(ArmName.HEAD, root, rev)
    identity = fake_bracket(arm, trees=iter([TreeIdentity(rev, False, None)] * 2))()()
    assert identity.native_sha256 == FAKE_NATIVE.native_sha256
    assert identity.native_sources_sha256 == FAKE_NATIVE.native_sources_sha256
    assert identity.pixi_lock_sha256 == sha256_file(root / "pixi.lock")
    assert (identity.rev, identity.dirty, identity.patch_sha256) == (rev, False, None)
    assert identity.editable_build_dir == str(arm.editable_build_dir)
    assert identity.interpreter_prefix == str(arm.interpreter_prefix)
    assert identity.compas_cgal_file == str(fake_origin(root).compas_cgal_file)


def test_an_unsettled_build_is_refused_before_timing() -> None:
    def stale() -> NativeIdentity:
        raise NativeBuildStaleError("a source is newer than the build")

    def tracker() -> NativeFinisher:
        return stale

    with pytest.raises(NativeBuildStaleError):
        require_settled(tracker)


def test_a_dirty_head_tree_is_refused_before_checkout(tmp_path: Path) -> None:
    root, first = _tree(tmp_path / "head")
    (root / "extra.txt").write_text("second\n")
    git(root, "add", "extra.txt")
    git(root, "commit", "-q", "-m", "second")
    before = _head(root)
    (root / "untracked.txt").write_text("dirt\n")
    with pytest.raises(MeasuredTreeError):
        checkout_head(root, first)
    assert _head(root) == before


def test_an_unknown_head_rev_is_refused_and_leaves_the_tree_alone(tmp_path: Path) -> None:
    root, rev = _tree(tmp_path / "head")
    with pytest.raises(HeadRevUnknownError):
        checkout_head(root, "0" * 40)
    assert tree_identity(root) == TreeIdentity(rev, False, None)


def test_a_failed_head_checkout_is_refused_by_name_quoting_git(tmp_path: Path) -> None:
    root, first = _tree(tmp_path / "head")
    (root / "extra.txt").write_text("second\n")
    git(root, "add", "extra.txt")
    git(root, "commit", "-q", "-m", "second")
    # A held index lock makes ``git checkout`` fail while the tree still reads clean.
    (root / ".git" / "index.lock").write_text("")
    with pytest.raises(HeadCheckoutFailedError, match="index.lock"):
        checkout_head(root, first)


def test_a_failed_interpreter_probe_is_refused_by_name_quoting_its_stderr(tmp_path: Path) -> None:
    failing = ("sh", "-c", "echo no interpreter here >&2; exit 3", "sh")
    with pytest.raises(ArmProbeFailedError, match="no interpreter here"):
        interpreter_prefix_of(failing, tmp_path, dict(os.environ))


def test_a_silent_interpreter_probe_is_refused_by_name(tmp_path: Path) -> None:
    silent = ("sh", "-c", "echo only a warning >&2", "sh")
    with pytest.raises(ArmProbeNoOutputError, match="only a warning"):
        interpreter_prefix_of(silent, tmp_path, dict(os.environ))


def test_checkout_moves_a_clean_head_to_the_rev(tmp_path: Path) -> None:
    root, first = _tree(tmp_path / "head")
    (root / "extra.txt").write_text("second\n")
    git(root, "add", "extra.txt")
    git(root, "commit", "-q", "-m", "second")
    assert _head(root) != first
    assert checkout_head(root, first[:12]) == first
    assert tree_identity(root) == TreeIdentity(first, False, None)


def test_the_arm_rebuild_pins_the_orchestrators_build(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    root, rev = _tree(tmp_path / "head")
    arm = fake_arm(ArmName.HEAD, root, rev)
    orchestrator_build = tmp_path / "orchestrator" / "build"
    runs: list[tuple[Path, dict[str, str], tuple[str, ...]]] = []

    def recording_run_in_tree(cwd: Path, env: Mapping[str, str]) -> Runner:
        def run(command: Sequence[str]) -> str:
            runs.append((cwd, dict(env), tuple(command)))
            return ""

        return run

    monkeypatch.setattr(perf_ratio_arms, "run_in_tree", recording_run_in_tree)
    command = rebuild_arm(arm, orchestrator_build)
    expected = ("pixi", "run", "--manifest-path", str(root / PYPROJECT), "_editable-rebuild")
    assert command == expected
    assert [(cwd, cmd) for cwd, _, cmd in runs] == [(root, expected)]
    assert runs[0][1]["SKBUILD_EDITABLE_SKIP"] == str(orchestrator_build)
    assert runs[0][1]["PYTHONPATH"] == str(root)


def test_live_bracket_probes_the_arms_own_interpreter(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    root, rev = _tree(tmp_path / "head")
    arm = fake_arm(ArmName.HEAD, root, rev)
    seen: dict[str, tuple[object, ...]] = {}

    def track_native_in(python: Sequence[str], cwd: Path, env: Mapping[str, str]) -> NativeTracker:
        seen["native"] = (tuple(python), cwd, env)
        return fake_tracker()

    def library_origin_of(python: Sequence[str], cwd: Path, env: Mapping[str, str]) -> LibraryOrigin:
        seen["origin"] = (tuple(python), cwd, env)
        return fake_origin(root)

    def interpreter_identity_of(python: Sequence[str], cwd: Path, env: dict[str, str]) -> tuple[str, str, str]:
        seen["interpreter"] = (tuple(python), cwd, env)
        return fake_interpreter()

    def read_tree(path: Path) -> TreeIdentity:
        seen["tree"] = (path,)
        return TreeIdentity(rev, False, None)

    monkeypatch.setattr(perf_ratio_arms, "track_native_in", track_native_in)
    monkeypatch.setattr(perf_ratio_arms, "library_origin_of", library_origin_of)
    monkeypatch.setattr(perf_ratio_arms, "interpreter_identity_of", interpreter_identity_of)
    monkeypatch.setattr(perf_ratio_arms, "tree_identity", read_tree)

    identity = live_bracket(arm)()()

    assert seen["tree"] == (root,)
    for name in ("native", "origin", "interpreter"):
        assert seen[name] == (arm.python, root, dict(arm.env)), name
    assert identity.rev == rev


def test_fake_arm_children_report_the_given_cpu_and_digest(tmp_path: Path) -> None:
    arms = _two_arms(tmp_path)
    one_per_family = {workload.family: workload for workload in WORKLOAD_SETS["default"]}
    assert set(one_per_family) == set(Family)
    for workload in one_per_family.values():
        runs = []
        for arm in arms:
            command = child_command(workload, arm.python, CpuSeconds(0.0))
            completed = subprocess.run(list(command), cwd=arm.root, env=dict(arm.env), capture_output=True, text=True, check=True)
            runs.append(parse_child(workload, completed.stdout))
        assert [run.cpu_seconds for run in runs] == [FAKE_CPU_S, FAKE_CPU_S], workload.name
        # Same digest argument, same output: no tree path leaks into the digest.
        assert runs[0].digest == runs[1].digest, workload.name
    for arm in arms:
        assert tree_identity(arm.root) == TreeIdentity(arm.rev, False, None)
