"""The gate ratchet: no pinned gate cell may worsen, and touching the gate must improve one.

Semantics (Ruling 0.3-C2), for a parent document B and a child document H:

* WORSENS when any gated cell of a case both carry is larger in H than in B, or
  when a case of B is missing from H. Checked on every commit, whatever it
  touches, so a kernel or generator regression fails even though neither is a
  gate path.
* IMPROVES when some gated cell of a shared case is smaller in H, or when a NEW
  case of H (a new generator column) beats, on some cell, the minimum over the
  existing generators for its (cap, pocket) and exceeds that minimum on no cell
  (option i). A new case that does not is a stagnation, not a regression.
* A commit that changes a path under `GATE_PATHS` must improve; one that does
  not need not (Ruling 0.3-C3).
* A different rules version re-baselines only the criteria whose definition
  changed (Ruling 0.3-bump-worsens). H must then carry a ``rebaseline`` record
  naming B's rules version, the approval (the ledger line) and the changed
  gated criteria; without it the ratchet refuses. Never silent. Every unchanged
  criterion is still checked for worsening, and every case of B must still be
  in H, exactly as for a same-version pair. A gated criterion of H that B does
  not gate must be listed. Every criterion B gates that H does not (removed,
  or demoted to ``reported``) must be listed in the record's ``retired``, and
  ``retired`` may name only such criteria. A recorded bump counts as the
  commit's improvement (it re-baselines or retires at least one criterion, the
  record being validated non-empty), so a bump that touches the gate need not
  also improve a cell. Under one rules version the gated criteria may not
  change at all.
* A parent with no document is a baseline: the child passes.

`check` walks every commit of a range and compares the document at each commit
with the document at its parent (`git show`), which costs no measurement;
`fresh` re-measures the head (`tools.gate_measure`) and requires every gated
cell and every reported integer to match exactly (Ruling 0.3-CI-fresh: reported
floats are never compared, because their last digit is not guaranteed to agree
across platforms).
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path
from typing import Sequence

from tools.gate_vectors import CASE_KEY_SEPARATOR
from tools.gate_vectors import GATE_VECTORS_PATH
from tools.gate_vectors import GateVectors
from tools.gate_vectors import MalformedGateVectorsError
from tools.gate_vectors import load_gate_vectors
from tools.gate_vectors import parse_gate_vectors

# Prefixes of the repository paths whose change "touches the gate" (Ruling
# 0.3-C3): the Python gate modules and the pipeline generator. C++ kernels are
# out of scope; their regressions are still caught, because worsening is
# checked on every commit. `tests/tools/test_gate_ratchet.py` holds this list to
# the import closure of the gate.
GATE_PATHS: tuple[str, ...] = (
    "src/compas_cgal/pipeline/",
    "benchmarks/survey.py",
    "benchmarks/depletion.py",
    "benchmarks/quality.py",
    "benchmarks/quality_observations.py",
    "benchmarks/coverage.py",
    "benchmarks/gate.py",
    "benchmarks/spec.py",
    "src/compas_cgal/replay_classification.py",
    "src/compas_cgal/engagement.py",
    "src/compas_cgal/stock.py",
)


# The range `check` starts after: the environment's GATE_BASE (CI sets it to the
# pushed-from or pull-request base commit), else the parent of HEAD.
GATE_BASE_ENV = "GATE_BASE"
DEFAULT_GATE_BASE = "HEAD^"

# When the range base is not an ancestor of the head -- the pushed-from commit
# of a force-push after a rebase -- the range starts at the merge base of the
# head with this ref instead (CI sets it to origin/main). Unset, an unreachable
# base is a named failure.
GATE_FALLBACK_REF_ENV = "GATE_FALLBACK_REF"


class GateVectorRegressionError(Exception):
    """A gated cell worsened, or a pinned case was dropped."""


class GateVectorStagnationError(Exception):
    """The commit touches the gate and improves no gated cell."""


class GateVectorsMissingError(Exception):
    """A commit that had gate vectors removed them."""


class GateRulesChangedError(Exception):
    """The rules version changed without a recorded re-baseline approval."""


class UnacknowledgedRetiredCriterionError(GateRulesChangedError):
    """A rules bump stopped gating a criterion the parent gated without listing it as retired."""


class SpuriousRetiredCriterionError(GateRulesChangedError):
    """A rules bump lists as retired a criterion the parent did not gate."""


class GateVectorsStaleError(Exception):
    """The committed vectors do not match a fresh measurement."""


class GitQueryError(Exception):
    """A git query the ratchet depends on failed."""


class UnreachableGateBaseError(Exception):
    """The range base is not an ancestor of the head (a force-pushed range)."""


@dataclass(frozen=True)
class RatchetVerdict:
    """What one parent-to-child comparison found.

    Attributes:
        improved: The cells that got better, as ``case/criterion`` strings, and
            the new cases that count as improvements.
        rebaselined: The criteria a recorded rules bump re-baselined, in record
            order; empty for a same-version pair.
        retired: The criteria a recorded rules bump stopped gating, in record
            order; empty for a same-version pair.
    """

    improved: tuple[str, ...]
    rebaselined: tuple[str, ...]
    retired: tuple[str, ...] = ()


def _cap_and_pocket(case: str) -> tuple[str, str]:
    cap_id, _generator, pocket = case.split(CASE_KEY_SEPARATOR)
    return cap_id, pocket


def _new_case_improves(base: GateVectors, head: GateVectors, case: str, compared: Sequence[str]) -> bool:
    """Option (i): a new column improves iff it beats the per-(cap, pocket) minimum somewhere and loses nowhere."""
    siblings = [key for key in base.cases if _cap_and_pocket(key) == _cap_and_pocket(case)]
    if not siblings:
        return False
    beats = False
    for name in compared:
        minimum = min(base.value(key, name) for key in siblings)
        value = head.value(case, name)
        if value > minimum:
            return False
        beats = beats or value < minimum
    return beats


def check_pair(base: GateVectors | None, head: GateVectors, *, touches_gate: bool) -> RatchetVerdict:
    """Judge one commit's gate vectors against its parent's.

    Args:
        base: The parent's vectors, or None when the parent has none.
        head: The commit's vectors.
        touches_gate: Whether the commit changes a path under `GATE_PATHS`.

    Returns:
        The improvements found.

    Raises:
        GateVectorRegressionError: A compared gated cell worsened or a case was
            dropped.
        GateVectorStagnationError: The gate was touched and nothing improved.
        GateRulesChangedError: The rules version changed without a recorded
            re-baseline from the parent's version, a bump left a newly gated
            criterion unlisted, or the gated criteria changed without a bump.
        UnacknowledgedRetiredCriterionError: A bump stopped gating a criterion
            the parent gated without listing it as retired.
        SpuriousRetiredCriterionError: A bump lists as retired a criterion the
            parent did not gate.
    """
    if base is None:
        return RatchetVerdict(improved=(), rebaselined=())
    rebaselined = _rebaselined_criteria(base, head)
    retired = _retired_criteria(base, head)
    base_gated = {criterion.name for criterion in base.gated}
    compared = [criterion.name for criterion in head.gated if criterion.name not in rebaselined]
    unlisted = [name for name in compared if name not in base_gated]
    if unlisted:
        raise GateRulesChangedError(f"rules {base.rules_version!r} -> {head.rules_version!r} gates {unlisted} newly; list them in the rebaseline record's criteria.")

    dropped = sorted(set(base.cases) - set(head.cases))
    worse: list[str] = [f"{case} (dropped)" for case in dropped]
    improved: list[str] = []
    for case in base.cases:
        if case not in head.cases:
            continue
        for name in compared:
            before, after = base.value(case, name), head.value(case, name)
            if after > before:
                worse.append(f"{case}/{name}: {before} -> {after}")
            elif after < before:
                improved.append(f"{case}/{name}: {before} -> {after}")
    if worse:
        raise GateVectorRegressionError("gate cells worsened:\n  " + "\n  ".join(worse))
    improved.extend(f"{case} (new case beats the existing minimum)" for case in head.cases if case not in base.cases and _new_case_improves(base, head, case, compared))
    if touches_gate and not improved and not rebaselined and not retired:
        raise GateVectorStagnationError("the commit touches the gate and improves no gated cell.")
    return RatchetVerdict(improved=tuple(improved), rebaselined=rebaselined, retired=retired)


def _retired_criteria(base: GateVectors, head: GateVectors) -> tuple[str, ...]:
    """The criteria a bump stops gating, exactly as its record acknowledges them.

    Called after `_rebaselined_criteria`, so under one rules version the gated
    criteria are already equal and nothing is retired.

    Args:
        base: The parent's vectors.
        head: The commit's vectors.

    Returns:
        The retired criterion names, in record order.

    Raises:
        UnacknowledgedRetiredCriterionError: The parent gates a criterion the
            commit does not, and the record does not list it as retired.
        SpuriousRetiredCriterionError: The record retires a criterion the
            parent did not gate.
    """
    if head.rules_version == base.rules_version:
        return ()
    if head.rebaseline is None:
        raise GateRulesChangedError(f"rules version {base.rules_version!r} -> {head.rules_version!r} needs a recorded re-baseline.")
    head_gated = {criterion.name for criterion in head.gated}
    base_gated = [criterion.name for criterion in base.gated]
    unacknowledged = [name for name in base_gated if name not in head_gated and name not in head.rebaseline.retired]
    if unacknowledged:
        raise UnacknowledgedRetiredCriterionError(
            f"rules {base.rules_version!r} -> {head.rules_version!r} stops gating {unacknowledged}; list them in the rebaseline record's 'retired'."
        )
    spurious = [name for name in head.rebaseline.retired if name not in base_gated]
    if spurious:
        raise SpuriousRetiredCriterionError(f"the rebaseline record retires {spurious}, which rules {base.rules_version!r} did not gate.")
    return head.rebaseline.retired


def _rebaselined_criteria(base: GateVectors, head: GateVectors) -> tuple[str, ...]:
    """The criteria the pair re-baselines: none under one rules version, the record's list across a bump.

    Args:
        base: The parent's vectors.
        head: The commit's vectors.

    Returns:
        The re-baselined criterion names, in record order.

    Raises:
        GateRulesChangedError: A bump without a record from the parent's
            version, or a gated-criteria change without a bump.
    """
    if head.rules_version == base.rules_version:
        if [criterion.name for criterion in head.gated] != [criterion.name for criterion in base.gated]:
            raise GateRulesChangedError("the gated criteria changed without a rules-version bump.")
        return ()
    if head.rebaseline is None or head.rebaseline.from_version != base.rules_version:
        raise GateRulesChangedError(
            f"rules version {base.rules_version!r} -> {head.rules_version!r} needs a recorded re-baseline; record it as rebaseline "
            f"{{'from': {base.rules_version!r}, 'approval': '<user-approved ledger line>', 'criteria': [<changed gated criteria>], "
            "'retired': [<criteria no longer gated>]}} in the vectors."
        )
    return head.rebaseline.criteria


def touches_gate(paths: Sequence[str]) -> bool:
    """Whether any changed path lies under `GATE_PATHS`.

    Args:
        paths: Repository-relative paths, as ``git diff --name-only`` prints them.

    Returns:
        True when the change touches the gate.
    """
    return any(path.startswith(prefix) for path in paths for prefix in GATE_PATHS)


def _git(repo: Path, *args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(["git", *args], cwd=repo, capture_output=True, text=True, check=False)


def _git_out(repo: Path, *args: str) -> str:
    done = _git(repo, *args)
    if done.returncode != 0:
        raise GitQueryError(f"git {' '.join(args)} failed ({done.returncode}): {done.stderr.strip()}")
    return done.stdout


def _parent(repo: Path, commit: str) -> str | None:
    """The first parent of *commit*, or None for a root commit."""
    parents = _git_out(repo, "rev-list", "--parents", "-n", "1", commit).split()[1:]
    return parents[0] if parents else None


def _vectors_at(repo: Path, commit: str | None) -> GateVectors | None:
    """The vectors committed at *commit*; None for no commit or no file there."""
    if commit is None:
        return None
    listing = _git_out(repo, "ls-tree", "--name-only", commit, "--", GATE_VECTORS_PATH.as_posix()).split()
    if not listing:
        return None
    return parse_gate_vectors(_git_out(repo, "show", f"{commit}:{GATE_VECTORS_PATH.as_posix()}"))


def resolve_base(repo: Path, base: str, head: str, fallback_ref: str | None) -> str:
    """The commit a ratchet range starts after.

    Args:
        repo: The repository root.
        base: The requested base.
        head: The range head.
        fallback_ref: A ref whose merge base with *head* replaces a base that
            is not an ancestor of *head*, or None to refuse.

    Returns:
        *base*, or the merge base with *fallback_ref*.

    Raises:
        UnreachableGateBaseError: *base* is unknown or not an ancestor of *head*
            and no fallback is given.
        GitQueryError: A git query failed.
    """
    known = _git(repo, "cat-file", "-e", f"{base}^{{commit}}").returncode == 0
    if known and _git(repo, "merge-base", "--is-ancestor", base, head).returncode == 0:
        return base
    if fallback_ref is None:
        raise UnreachableGateBaseError(f"base {base!r} is not an ancestor of {head!r} (force-pushed range?); set {GATE_FALLBACK_REF_ENV} to ratchet from the merge base.")
    merge_base = _git_out(repo, "merge-base", fallback_ref, head).strip()
    print(f"gate ratchet: base {base[:10]} is not an ancestor of {head}; ratcheting from merge-base({fallback_ref}, {head}) = {merge_base[:10]}")
    return merge_base


def check_range(repo: Path, base: str, head: str) -> list[tuple[str, RatchetVerdict]]:
    """Ratchet every commit of ``base..head`` against its parent.

    Args:
        repo: The repository root.
        base: The commit the range starts after.
        head: The last commit of the range.

    Returns:
        Each commit with its verdict, oldest first.

    Raises:
        GateVectorsMissingError: A commit removed the vectors its parent had.
        GateVectorRegressionError: See `check_pair`.
        GateVectorStagnationError: See `check_pair`.
        GateRulesChangedError: See `check_pair`.
        GitQueryError: A git query failed.
    """
    verdicts: list[tuple[str, RatchetVerdict]] = []
    for commit in _git_out(repo, "rev-list", "--reverse", f"{base}..{head}").split():
        parent = _parent(repo, commit)
        parent_vectors = _vectors_at(repo, parent)
        vectors = _vectors_at(repo, commit)
        if vectors is None:
            if parent_vectors is not None:
                raise GateVectorsMissingError(f"commit {commit} removed {GATE_VECTORS_PATH}.")
            continue
        # A root commit's change set is its whole tree.
        changed = (_git_out(repo, "ls-tree", "-r", "--name-only", commit) if parent is None else _git_out(repo, "diff", "--name-only", parent, commit)).split()
        verdicts.append((commit, check_pair(parent_vectors, vectors, touches_gate=touches_gate(changed))))
    return verdicts


def stale_cells(committed: GateVectors, measured: GateVectors) -> list[str]:
    """Gated cells and reported integers whose committed value differs from a fresh measurement.

    Args:
        committed: The vectors in the tree.
        measured: A fresh measurement.

    Returns:
        ``case/criterion: committed != measured`` lines, empty when fresh.
    """
    stale: list[str] = []
    if set(committed.cases) != set(measured.cases):
        stale.append(f"case sets differ: {sorted(set(committed.cases) ^ set(measured.cases))}")
    for case in sorted(set(committed.cases) & set(measured.cases)):
        for criterion in committed.freshness_checked:
            before = committed.cases[case][criterion.name]
            after = measured.cases[case].get(criterion.name)
            if after is None or Decimal(before) != Decimal(after):
                stale.append(f"{case}/{criterion.name}: committed {before} != measured {after}")
    return stale


def main(argv: Sequence[str] | None = None) -> int:
    """Command-line entry: ``check``, ``check-files``, ``fresh``.

    Args:
        argv: Arguments, defaulting to ``sys.argv[1:]``.

    Returns:
        The process exit status.
    """
    parser = argparse.ArgumentParser(prog="python -m tools.gate_ratchet", description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)
    check = sub.add_parser("check", help="ratchet every commit of BASE..HEAD")
    check.add_argument("--base", default=os.environ.get(GATE_BASE_ENV, DEFAULT_GATE_BASE), help=f"default: ${GATE_BASE_ENV}, else {DEFAULT_GATE_BASE}")
    check.add_argument("--head", default="HEAD")
    check.add_argument("--fresh", action="store_true", help="also re-measure the working tree's vectors")
    check.add_argument("--fallback-ref", default=os.environ.get(GATE_FALLBACK_REF_ENV), help=f"default: ${GATE_FALLBACK_REF_ENV}; merge-base fallback for an unreachable base")
    files = sub.add_parser("check-files", help="ratchet two vector files")
    files.add_argument("--base", type=Path, required=True)
    files.add_argument("--head", type=Path, required=True)
    files.add_argument("--touches-gate", action="store_true")
    fresh = sub.add_parser("fresh", help="re-measure and compare with the committed vectors")
    fresh.add_argument("--vectors", type=Path, default=GATE_VECTORS_PATH)
    args = parser.parse_args(argv)
    repo = Path.cwd()
    try:
        if args.command == "check-files":
            verdict = check_pair(load_gate_vectors(args.base), load_gate_vectors(args.head), touches_gate=args.touches_gate)
            print(f"gate ratchet: pass ({len(verdict.improved)} improved)")
            return 0
        if args.command == "check":
            for commit, verdict in check_range(repo, resolve_base(repo, args.base, args.head, args.fallback_ref), args.head):
                rebaselined = f"re-baselined {', '.join(verdict.rebaselined)}; " if verdict.rebaselined else ""
                retired = f"retired {', '.join(verdict.retired)}; " if verdict.retired else ""
                kind = f"{rebaselined}{retired}{len(verdict.improved)} improved"
                print(f"{commit[:10]}: pass ({kind})")
            if not args.fresh:
                return 0
            return _fresh(GATE_VECTORS_PATH)
        return _fresh(args.vectors)
    except (
        GateVectorRegressionError,
        GateVectorStagnationError,
        GateRulesChangedError,
        GateVectorsMissingError,
        GateVectorsStaleError,
        MalformedGateVectorsError,
        GitQueryError,
        UnreachableGateBaseError,
    ) as error:
        print(f"gate ratchet: FAIL {type(error).__name__}: {error}", file=sys.stderr)
        return 1


def _fresh(path: Path) -> int:
    from tools.gate_measure import measure_gate_vectors

    stale = stale_cells(load_gate_vectors(path), measure_gate_vectors())
    if stale:
        raise GateVectorsStaleError(f"{path} is stale:\n  " + "\n  ".join(stale))
    print(f"gate vectors fresh: {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
