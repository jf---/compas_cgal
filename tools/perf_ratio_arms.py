"""The two arms of a ratio session: setup, distinctness, HEAD checkout, rebuild and the per-arm bracket.

An arm is a worktree run by its own interpreter (`arm_python`: ``pixi run`` with
the arm's own manifest, never the orchestrator's). `open_arm` proves the tree is
the expected clean commit and that the arm's interpreter loads its library from
inside the tree; `require_distinct_arms` refuses two arms that share a tree,
an editable build or an interpreter prefix (which would time one environment
twice and give a silent 1.0). `bracket_arm` brackets a session per arm: the
tree, library and native build must be the same at the end as at the start,
and the end yields the arm's `ArmIdentity`. Both arms are identified the same
way; HEAD is never special-cased in-process.

The generic tree, library and rebuild proofs are reused from
``tools.pipeline_timing_1_1_base``; their error texts name the 1.1 tree, but
their logic is tree-agnostic.
"""

from __future__ import annotations

import os
import subprocess
from collections.abc import Mapping
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

from tools.perf_ratio_record import ArmIdentity
from tools.perf_ratio_record import ArmName
from tools.pipeline_timing_1_1_base import LibraryOrigin
from tools.pipeline_timing_1_1_base import MeasuredTreeError
from tools.pipeline_timing_1_1_base import base_child_env
from tools.pipeline_timing_1_1_base import editable_build_dir_of
from tools.pipeline_timing_1_1_base import library_origin_of
from tools.pipeline_timing_1_1_base import rebuild_base_build
from tools.pipeline_timing_1_1_base import require_library_in_tree
from tools.pipeline_timing_1_1_base import require_measured_tree
from tools.pipeline_timing_1_1_base import run_in_tree
from tools.run_native_identity import NativeTracker
from tools.run_native_identity import track_native_in
from tools.run_provenance import TreeReader
from tools.run_provenance import interpreter_identity_of
from tools.run_provenance import require_unchanged_tree
from tools.run_provenance import sha256_file
from tools.run_provenance import tree_identity

# The pinned reference commit every ratio is taken against.
PINNED_REF_COMMIT = "15fbed70ebb016bdd8ef49453e04606d03e0b1fe"
ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REF_WORKTREE = ROOT.parent / "compas_cgal_prs-ratio-ref"
DEFAULT_HEAD_WORKTREE = ROOT.parent / "compas_cgal_prs-ratio-head"
PIXI_LOCK = Path("pixi.lock")
PYPROJECT = Path("pyproject.toml")
# Prints the interpreter's environment prefix; two arms sharing it run one environment.
INTERPRETER_PREFIX_SNIPPET = "import sys; print(sys.prefix)"

# Starts an arm's bracket; the returned callable ends it and yields the arm's identity.
ArmBracket = Callable[[], Callable[[], ArmIdentity]]


class SameArmError(Exception):
    """The two arms share a tree, an editable build or an interpreter prefix: the session would time one arm twice."""


class HeadRevUnknownError(Exception):
    """The requested HEAD revision names no commit in the HEAD worktree's repository."""


class HeadCheckoutFailedError(Exception):
    """``git checkout`` of a known commit in the clean HEAD worktree exited non-zero."""


class ArmProbeFailedError(Exception):
    """The arm interpreter's prefix probe exited non-zero."""


class ArmProbeNoOutputError(Exception):
    """The arm interpreter's prefix probe succeeded but printed no prefix."""


@dataclass(frozen=True)
class Arm:
    """One arm of a ratio session, as `open_arm` proved it.

    Attributes:
        name: Which arm.
        root: The arm's worktree.
        python: The arm interpreter's command prefix.
        env: The pinned environment every arm child runs with (`base_child_env`).
        editable_build_dir: The arm's editable native build directory.
        interpreter_prefix: ``sys.prefix`` of the arm interpreter, resolved.
        rev: Full commit id the arm's tree was proved clean at.
    """

    name: ArmName
    root: Path
    python: tuple[str, ...]
    env: Mapping[str, str]
    editable_build_dir: Path
    interpreter_prefix: Path
    rev: str


def arm_python(root: Path) -> tuple[str, ...]:
    """The interpreter command of the arm at *root*: ``pixi run`` naming the arm's own manifest.

    Naming the manifest keeps the orchestrator's pixi environment from resolving
    for the arm when the variables of its activation leak into the child.
    """
    return ("pixi", "run", "--no-install", "--manifest-path", str(root / PYPROJECT), "python")


def interpreter_prefix_of(python: Sequence[str], root: Path, env: Mapping[str, str]) -> Path:
    """``sys.prefix`` (resolved) of the interpreter *python* run in *root* with *env*.

    Raises:
        ArmProbeFailedError: The probe exited non-zero (its stderr is quoted).
        ArmProbeNoOutputError: The probe printed nothing on stdout (its stderr is quoted).
    """
    command = [*python, "-c", INTERPRETER_PREFIX_SNIPPET]
    completed = subprocess.run(command, cwd=root, env=dict(env), capture_output=True, text=True, check=False)
    if completed.returncode != 0:
        raise ArmProbeFailedError(f"the prefix probe {command} in {root} exited {completed.returncode}: {completed.stderr.strip()!r}")
    lines = completed.stdout.strip().splitlines()
    if not lines:
        raise ArmProbeNoOutputError(f"the prefix probe {command} in {root} printed no prefix; stderr: {completed.stderr.strip()!r}")
    return Path(lines[-1]).resolve()


def open_arm(name: ArmName, root: Path, python: Sequence[str], expected_rev: str) -> Arm:
    """Prove the arm at *root* and return it.

    Order: the tree is *expected_rev* and clean; the editable build directory is
    read unpinned; every later probe runs pinned to it; the library resolves
    inside the tree; then the interpreter prefix is read.

    Args:
        name: Which arm.
        root: The arm's worktree.
        python: The arm interpreter's command prefix (`arm_python` for a real arm).
        expected_rev: The full commit id the tree must be at.

    Returns:
        The proved arm.

    Raises:
        MeasuredTreeError: The tree is another commit, or modified.
        ForeignLibraryError: The arm's library or editable build lies outside *root*.
    """
    command = tuple(python)
    require_measured_tree(root, expected_rev)
    editable_build_dir = editable_build_dir_of(command, root, base_child_env(root, os.environ, None))
    env = base_child_env(root, os.environ, editable_build_dir)
    require_library_in_tree(library_origin_of(command, root, env), root, editable_build_dir)
    prefix = interpreter_prefix_of(command, root, env)
    return Arm(name=name, root=root, python=command, env=env, editable_build_dir=editable_build_dir, interpreter_prefix=prefix, rev=expected_rev)


def require_distinct_arms(ref: Arm, head: Arm) -> None:
    """Refuse two arms that share a resolved root, editable build directory or interpreter prefix.

    Raises:
        SameArmError: Any of the three is equal for both arms.
    """
    shared = [
        label
        for label, ref_path, head_path in (
            ("root", ref.root, head.root),
            ("editable build", ref.editable_build_dir, head.editable_build_dir),
            ("interpreter prefix", ref.interpreter_prefix, head.interpreter_prefix),
        )
        if ref_path.resolve() == head_path.resolve()
    ]
    if shared:
        raise SameArmError(f"the ref and head arms share their {', '.join(shared)}: the session would time one arm twice.")


def checkout_head(root: Path, rev: str) -> str:
    """Move the clean HEAD worktree at *root* to the commit *rev* names, detached, and return its full id.

    Raises:
        MeasuredTreeError: The tree is modified (refused before anything moves), or is not at the commit afterwards.
        HeadRevUnknownError: *rev* names no commit (the tree is left alone).
        HeadCheckoutFailedError: ``git checkout`` itself failed (its stderr is quoted).
    """
    require_measured_tree(root, tree_identity(root).rev)
    resolved = subprocess.run(["git", "-C", str(root), "rev-parse", "--verify", "--quiet", f"{rev}^{{commit}}"], capture_output=True, text=True, check=False)
    if resolved.returncode != 0:
        raise HeadRevUnknownError(f"{rev!r} names no commit in {root}.")
    full = resolved.stdout.strip()
    checkout = subprocess.run(["git", "-C", str(root), "checkout", "-q", "--detach", full], capture_output=True, text=True, check=False)
    if checkout.returncode != 0:
        raise HeadCheckoutFailedError(f"git checkout --detach {full} in {root} exited {checkout.returncode}: {checkout.stderr.strip()!r}")
    require_measured_tree(root, full)
    return full


def rebuild_arm(arm: Arm, orchestrator_build_dir: Path) -> tuple[str, ...]:
    """Bring *arm*'s editable build up to date with its own rebuild task, and return the command run.

    The rebuild's environment pins the orchestrator's own build directory
    (*orchestrator_build_dir*), so the orchestrator's build is never rebuilt and
    the arm's is not skipped.

    Raises:
        BaseRebuildFailedError: The arm's rebuild task failed.
    """
    return rebuild_base_build(arm.root, run_in_tree(arm.root, base_child_env(arm.root, os.environ, orchestrator_build_dir)))


def require_settled(track_native: NativeTracker) -> None:
    """Run one full native bracket before timing, so a stale or changing build is refused up front.

    Raises:
        NativeIdentityError: The build is stale, changed, or not compiled from the tree (a subclass propagates).
    """
    track_native()()


def bracket_arm(
    arm: Arm,
    *,
    read_tree: TreeReader,
    track_native: NativeTracker,
    origin_of: Callable[[], LibraryOrigin],
    interpreter_of: Callable[[], tuple[str, str, str]],
) -> ArmBracket:
    """The per-arm bracket around a session: start reads the tree and starts the native bracket; finish proves and identifies.

    Args:
        arm: The proved arm.
        read_tree: Reads the arm's tree identity.
        track_native: The arm's native bracket.
        origin_of: Reads where the arm's library resolves.
        interpreter_of: Reads ``(environment, python_version, stock_2_sha256)`` of the arm interpreter.

    Returns:
        The bracket. Starting it raises `MeasuredTreeError` when the tree is not
        *arm*'s proved clean rev; finishing it raises `ForeignLibraryError`,
        `TreeChangedDuringRunError` or a `NativeIdentityError` and otherwise
        returns the arm's `ArmIdentity`.
    """
    label = f"the {arm.name.value} arm"

    def start() -> Callable[[], ArmIdentity]:
        at_start = read_tree()
        if at_start.rev != arm.rev or at_start.dirty:
            raise MeasuredTreeError(f"{label} at {arm.root} reads {at_start} at its bracket start, not the clean {arm.rev} it was opened at.")
        finish_native = track_native()

        def finish() -> ArmIdentity:
            origin = origin_of()
            require_library_in_tree(origin, arm.root, arm.editable_build_dir)
            at_end = read_tree()
            require_unchanged_tree(at_start, at_end, label)
            native = finish_native()
            interpreter = interpreter_of()
            return ArmIdentity.of(
                tree=at_end,
                native=native,
                interpreter=interpreter,
                pixi_lock_sha256=sha256_file(arm.root / PIXI_LOCK),
                compas_cgal_file=str(origin.compas_cgal_file),
                stock_2_origin=str(origin.stock_2_origin),
                editable_build_dir=str(arm.editable_build_dir),
                interpreter_prefix=str(arm.interpreter_prefix),
            )

        return finish

    return start


def live_bracket(arm: Arm) -> ArmBracket:
    """`bracket_arm` over the real readers, each run with *arm*'s own interpreter, tree and pinned environment."""
    return bracket_arm(
        arm,
        read_tree=lambda: tree_identity(arm.root),
        track_native=track_native_in(arm.python, arm.root, arm.env),
        origin_of=lambda: library_origin_of(arm.python, arm.root, arm.env),
        interpreter_of=lambda: interpreter_identity_of(arm.python, arm.root, dict(arm.env)),
    )
