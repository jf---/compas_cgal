"""Fake ratio arms for unit tests: tmp git trees whose tool modules are stubs.

Each stub prints its real tool's result marker (imported here, written into
the stub) with a fixed CPU and digest, so the harness's child commands and
parsers run unchanged against a fake tree and get exact numbers back. No
stub imports ``compas_cgal`` or this worktree's ``tools``.
"""

from __future__ import annotations

import os
import platform
import subprocess
import sys
from collections.abc import Iterator
from pathlib import Path

from tools.coverage_certificate_budget_gate import BUDGET_WORKER_MARKER
from tools.held_generation_clock import WORKER_MARKER
from tools.perf_ratio_arms import PIXI_LOCK
from tools.perf_ratio_arms import Arm
from tools.perf_ratio_arms import ArmBracket
from tools.perf_ratio_arms import bracket_arm
from tools.perf_ratio_record import ArmName
from tools.perf_ratio_workloads import CERTIFICATE_DIGEST_FIELDS
from tools.perf_ratio_workloads import CERTIFICATE_WORKLOAD
from tools.pipeline_timing import RESULT_LINE_PREFIX
from tools.pipeline_timing_1_1_base import LibraryOrigin
from tools.pipeline_timing_1_1_base import base_child_env
from tools.run_native_identity import NativeFinisher
from tools.run_native_identity import NativeIdentity
from tools.run_native_identity import NativeTracker
from tools.run_provenance import TreeIdentity

# The native identity the in-memory tracker reports (two distinct 64-hex digests).
FAKE_NATIVE = NativeIdentity("a" * 64, "b" * 64)
# (environment, python_version, stock_2_sha256) the fake interpreter reader reports.
FAKE_ENVIRONMENT = "default"
FAKE_STOCK_2_SHA256 = "c" * 64
# Operation count the GT stub reports (any fixed integer).
FAKE_GT_OPERATIONS = 1
# Ignored by the fake trees: bytecode the stub children write, the arm's build and interpreter dirs.
FAKE_GITIGNORE = "__pycache__/\nbuild/\n.env/\n"
_GIT_IDENTITY = ("-c", "user.name=fake-arm", "-c", "user.email=fake-arm@example.invalid")

_GENERATION_CLOCK_STUB = """\
import json
import sys

MARKER = {marker!r}
CPU_SECONDS = {cpu!r}
DIGEST = {digest!r}


def main(argv):
    clock = {{"stage_cpu_seconds": {{"stub": CPU_SECONDS}}, "stage_wall_seconds": {{"stub": CPU_SECONDS}}, "exit_cpu_seconds": 0.0, "exit_wall_seconds": 0.0}}
    print(MARKER + json.dumps({{"digest": {{"digest": DIGEST, "argv": list(argv)}}, "clock": clock}}), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
"""

_CERTIFICATE_STUB = """\
import json
import sys

MARKER = {marker!r}
ROW = {row!r}


def main(argv):
    print(MARKER + json.dumps(ROW), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
"""

# Reads time.perf_counter() exactly twice through the module-level `time`, as the real tool does,
# so the harness's GT clock stand-in sees its expected read count; reports the fixed CPU.
_PIPELINE_TIMING_STUB = """\
import json
import sys
import time

RESULT_LINE_PREFIX = {marker!r}
CPU_SECONDS = {cpu!r}
DIGEST = {digest!r}
OPERATIONS = {operations!r}


def main(argv):
    start = time.perf_counter()
    end = time.perf_counter()
    del start, end
    print(RESULT_LINE_PREFIX + json.dumps({{"generation_s": CPU_SECONDS, "operations": OPERATIONS, "digest": DIGEST + ":" + argv[-1]}}), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
"""


def git(root: Path, *args: str) -> str:
    """Run git in *root* with a fixed test identity and return its stripped stdout."""
    return subprocess.run(["git", *_GIT_IDENTITY, *args], cwd=root, capture_output=True, text=True, check=True).stdout.strip()


def fake_arm_tree(root: Path, *, cpu_seconds: float, digest: str) -> str:
    """Make *root* a one-commit git tree whose tool modules are stubs, and return its full rev.

    Args:
        root: An existing empty directory.
        cpu_seconds: The CPU every stub child reports.
        digest: The output digest every stub child reports (combined with its workload argument).

    Returns:
        The full commit id of the tree's only commit.
    """
    git(root, "init", "-q")
    tools = root / "tools"
    tools.mkdir()
    (tools / "__init__.py").write_text("")
    (tools / "held_generation_clock.py").write_text(_GENERATION_CLOCK_STUB.format(marker=WORKER_MARKER, cpu=cpu_seconds, digest=digest))
    row: dict[str, object] = {name: digest for name in CERTIFICATE_DIGEST_FIELDS}
    row.update(workload=CERTIFICATE_WORKLOAD, certificate_cpu_seconds=cpu_seconds, certificate_wall_seconds=cpu_seconds)
    (tools / "coverage_certificate_budget_gate.py").write_text(_CERTIFICATE_STUB.format(marker=BUDGET_WORKER_MARKER, row=row))
    (tools / "pipeline_timing.py").write_text(_PIPELINE_TIMING_STUB.format(marker=RESULT_LINE_PREFIX, cpu=cpu_seconds, digest=digest, operations=FAKE_GT_OPERATIONS))
    (root / PIXI_LOCK).write_text(f"# fake arm lock\nversion: 6\ndigest: {digest}\n")
    (root / ".gitignore").write_text(FAKE_GITIGNORE)
    git(root, "add", "-A")
    git(root, "commit", "-q", "-m", "fake arm")
    return git(root, "rev-parse", "HEAD")


def fake_arm(name: ArmName, root: Path, rev: str) -> Arm:
    """An `Arm` over the fake tree *root* at *rev*, run by this test interpreter."""
    build_dir = root / "build" / "editable"
    return Arm(
        name=name,
        root=root,
        python=(sys.executable,),
        env=base_child_env(root, os.environ, build_dir),
        editable_build_dir=build_dir,
        interpreter_prefix=root / ".env",
        rev=rev,
    )


def fake_origin(root: Path) -> LibraryOrigin:
    """A library origin inside *root* (resolved, as `library_origin_of` reports it)."""
    tree = root.resolve()
    return LibraryOrigin(compas_cgal_file=tree / "src" / "compas_cgal" / "__init__.py", stock_2_origin=tree / "build" / "editable" / "_stock_2.so")


def fake_tracker(identity: NativeIdentity = FAKE_NATIVE) -> NativeTracker:
    """An in-memory native tracker whose bracket always reports *identity*."""

    def start() -> NativeFinisher:
        return lambda: identity

    return start


def fake_interpreter() -> tuple[str, str, str]:
    """(environment, python_version, stock_2_sha256) of a fake arm interpreter."""
    return FAKE_ENVIRONMENT, platform.python_version(), FAKE_STOCK_2_SHA256


def fake_bracket(arm: Arm, *, trees: Iterator[TreeIdentity]) -> ArmBracket:
    """`bracket_arm` over *arm* with tree reads from *trees*, an in-memory tracker and an origin inside the tree."""
    return bracket_arm(arm, read_tree=lambda: next(trees), track_native=fake_tracker(), origin_of=lambda: fake_origin(arm.root), interpreter_of=fake_interpreter)
