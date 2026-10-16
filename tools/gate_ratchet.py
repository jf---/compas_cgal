"""The gate ratchet: no pinned gate cell may worsen, and touching the gate must improve one.

Semantics (Ruling 0.3-C2), for a parent document B and a child document H:

* WORSENS when any gated cell of a case both carry is larger in H than in B, or
  when a case of B is missing from H. Checked on every commit, whatever it
  touches, so a kernel or generator regression fails even though neither is a
  gate path.
* IMPROVES when some gated cell of a shared case is smaller in H, or when a NEW
  case of H (a new generator column) beats, on some cell, the minimum over the
  existing generators for its (cap, pocket) and exceeds that minimum on no cell
  (option i), or when H appends an admission entry (Ruling 0.3-C2-amend).
  Across a rules bump, a re-baselined criterion's minimum is read from H's
  pre-existing cases, the only cells measured under the same definition.
* Every NEW case of H must pass option (i) or be admitted by an entry H
  appends to ``admissions``; otherwise the commit fails, whatever it touches
  and whatever else improves (Ruling 0.3-C2-amend). An appended entry's cases
  must be new. On EVERY commit, each admitted case of H (under any entry)
  whose generator names a reference must be no worse than the reference
  generator's case at the same (cap, pocket) on every gated cell of H,
  re-baselined criteria included (user ruling 2026-10-15); a null reference
  makes the admitted cells the baseline on the ruling's authority. B's entries must be H's first entries,
  unchanged: admissions are append-only.
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

from tools.gate_admissions import Admission
from tools.gate_admissions import MalformedAdmissionError
from tools.gate_case_keys import CASE_KEY_SEPARATOR
from tools.gate_case_keys import cap_and_pocket
from tools.gate_case_keys import case_generator
from tools.gate_vectors import GATE_VECTORS_PATH
from tools.gate_vectors import SCHEMA
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
    # Ruling 1.1-closure.
    "src/compas_cgal/incident_loop.py",
    "src/compas_cgal/toolpath.py",
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


class UnadmittedCaseError(Exception):
    """A new case neither beats the existing minimum (option i) nor is covered by an admission this commit adds."""


class LegacyHeadVectorsError(Exception):
    """The commit's vectors are not the current schema; only a parent may be read as an older one."""


class AdmissionRecordAlteredError(Exception):
    """An admission entry of the parent is missing or changed in the child (admissions are append-only)."""


class AdmissionAboveReferenceError(Exception):
    """An admitted case is worse than its reference generator's case on some gated cell."""


class MissingReferenceCaseError(Exception):
    """An admitted case's reference generator has no case at the same (cap, pocket)."""


class AdmissionOfExistingCaseError(Exception):
    """An admission this commit adds names a case the parent already pins."""


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


def _new_case_improves(base: GateVectors, head: GateVectors, case: str, compared: Sequence[str]) -> bool:
    """Option (i): a new column improves iff it beats the per-(cap, pocket) minimum somewhere and loses nowhere.

    The minimum of a compared criterion is over the parent's cells; that of a
    criterion a bump re-baselines is over the child's cells of the same
    pre-existing cases, measured under the new definition.
    """
    siblings = [key for key in base.cases if cap_and_pocket(key) == cap_and_pocket(case) and key in head.cases]
    if not siblings:
        return False
    beats = False
    for criterion in head.gated:
        name = criterion.name
        source = base if name in compared else head
        minimum = min(source.value(key, name) for key in siblings)
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
        AdmissionRecordAlteredError: An admission of the parent is missing or
            changed in the child.
        AdmissionOfExistingCaseError: An appended admission names a case the
            parent pins.
        MissingReferenceCaseError: An admitted case's reference generator has
            no case at its (cap, pocket).
        AdmissionAboveReferenceError: An admitted case is worse than its
            reference on a gated cell.
        UnadmittedCaseError: A new case neither passes option (i) nor is
            admitted by an appended entry.
        LegacyHeadVectorsError: The commit's document is not the current
            schema; only a parent may be read as ``gate-vectors/1``.
    """
    if head.schema != SCHEMA:
        raise LegacyHeadVectorsError(f"the commit's gate vectors are {head.schema!r}; write them as {SCHEMA!r} (only a parent may be an older schema).")
    _check_references(head)
    if base is None:
        return RatchetVerdict(improved=(), rebaselined=())
    appended = _appended_admissions(base, head)
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
    admitted = {case: admission.ruling for admission in appended for case in admission.cases}
    unadmitted: list[str] = []
    for case in head.cases:
        if case in base.cases or case in admitted:
            continue
        if _new_case_improves(base, head, case, compared):
            improved.append(f"{case} (new case beats the existing minimum)")
        else:
            unadmitted.append(case)
    if unadmitted:
        raise UnadmittedCaseError(
            f"new cases {unadmitted} neither beat the existing minimum for their (cap, pocket) somewhere while losing nowhere (option i), "
            "nor are admitted by an entry this commit appends to 'admissions'."
        )
    improved.extend(f"{case} (admitted under {ruling})" for case, ruling in admitted.items())
    if touches_gate and not improved and not rebaselined and not retired:
        raise GateVectorStagnationError("the commit touches the gate and improves no gated cell.")
    return RatchetVerdict(improved=tuple(improved), rebaselined=rebaselined, retired=retired)


def _appended_admissions(base: GateVectors, head: GateVectors) -> tuple[Admission, ...]:
    """The admissions the child appends, each naming only new cases.

    Args:
        base: The parent's vectors.
        head: The commit's vectors.

    Returns:
        The appended entries, in append order.

    Raises:
        AdmissionRecordAlteredError: The parent's entries are not the child's
            first entries, unchanged.
        AdmissionOfExistingCaseError: An appended entry names a case the
            parent pins.
    """
    kept = head.admissions[: len(base.admissions)]
    if kept != base.admissions:
        altered = [admission.ruling for position, admission in enumerate(base.admissions) if position >= len(kept) or kept[position] != admission]
        raise AdmissionRecordAlteredError(f"admissions are append-only: the parent's entries {altered} are missing or changed.")
    appended = head.admissions[len(base.admissions) :]
    for admission in appended:
        pinned = [case for case in admission.cases if case in base.cases]
        if pinned:
            raise AdmissionOfExistingCaseError(f"admission {admission.ruling} names {pinned}, which the parent already pins; pre-existing cells are ratcheted, never admitted.")
    return appended


def _check_references(head: GateVectors) -> None:
    """The reference rule on every commit: each admitted case of the child, under every entry, against its reference.

    Between rules bumps the ratchet implies it; across a bump it is what keeps
    a re-baselined criterion from lifting an admitted column above its
    reference (user ruling 2026-10-15). Both sides are the child's cells.

    Raises:
        MissingReferenceCaseError: See `_check_reference`.
        AdmissionAboveReferenceError: See `_check_reference`.
    """
    for admission in head.admissions:
        for case in admission.cases:
            _check_reference(head, admission, case)


def _check_reference(head: GateVectors, admission: Admission, case: str) -> None:
    """The reference rule: an admitted case is no worse than its reference on any gated cell of the child.

    Raises:
        MissingReferenceCaseError: The reference generator has no case at
            *case*'s (cap, pocket).
        AdmissionAboveReferenceError: *case* exceeds the reference on a gated
            cell.
    """
    reference = admission.reference[case_generator(case)]
    if reference is None:
        return
    cap_id, pocket = cap_and_pocket(case)
    reference_case = CASE_KEY_SEPARATOR.join((cap_id, reference, pocket))
    if reference_case not in head.cases:
        raise MissingReferenceCaseError(f"admission {admission.ruling}: {case} references {reference!r}, which has no case {reference_case}.")
    above = [
        f"{case}/{criterion.name}: {head.value(case, criterion.name)} > {reference_case} {head.value(reference_case, criterion.name)}"
        for criterion in head.gated
        if head.value(case, criterion.name) > head.value(reference_case, criterion.name)
    ]
    if above:
        raise AdmissionAboveReferenceError(f"admission {admission.ruling}: admitted cells exceed their reference:\n  " + "\n  ".join(above))


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
        UnadmittedCaseError: See `check_pair`.
        AdmissionRecordAlteredError: See `check_pair`.
        AdmissionOfExistingCaseError: See `check_pair`.
        MissingReferenceCaseError: See `check_pair`.
        AdmissionAboveReferenceError: See `check_pair`.
        LegacyHeadVectorsError: See `check_pair`.
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
        UnadmittedCaseError,
        LegacyHeadVectorsError,
        AdmissionRecordAlteredError,
        AdmissionOfExistingCaseError,
        MissingReferenceCaseError,
        AdmissionAboveReferenceError,
        MalformedAdmissionError,
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

    committed = load_gate_vectors(path)
    stale = stale_cells(committed, measure_gate_vectors(admissions=committed.admissions))
    if stale:
        raise GateVectorsStaleError(f"{path} is stale:\n  " + "\n  ".join(stale))
    print(f"gate vectors fresh: {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
