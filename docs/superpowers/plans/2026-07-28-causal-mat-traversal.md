# Causal MAT Traversal and Covered Generation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> superpowers:executing-plans to implement this plan task-by-task. Steps use
> checkbox (`- [ ]`) syntax for tracking.

**Goal:** Derive neck scope from exact global MAT side history, traverse every
supported MAT edge, generate a terminal covered operation stream, and close it
with fresh replay.

**Architecture:** An additive native topology projection preserves exact neck
loci and every separating-cut partition. `MatTraversalState` owns global
directed graph cursors and causal sides; `GenerationState` continues to own
the physical stock/coverage boundary. A distinct entry-circle bootstrap and
the existing atomic link/circle engine are cross-bound by immutable traversal
commits.

**Tech Stack:** C++20, CGAL exact kernels, nanobind, Python 3.12, COMPAS framed
geometry, frozen typed dataclasses, CCAN, SHA-256, pytest-xdist,
pytest-testmon, Ruff, strict mypy, strict MkDocs.

## Global constraints

- Exact topology comes from the retained native evidence owner. No Python
  evidence parser, coordinate adjacency, sampled side test, or epsilon.
- Preserve the frozen 20-field numeric projection; add exact owner
  properties and cross-validate them against the existing cut union.
- `GenerationState` and `MatTraversalState` remain separate lifetimes.
- One authoritative candidate-evaluation engine owns containment,
  certification, depletion, coverage, and witness construction.
- No `max_passes`, skip, xfail, fallback, partial certificate, or silent
  branch drop.
- Every test has a geometric/proof-boundary docstring.
- Every pytest command uses `-n auto`; affected tests use `--testmon`.
- Update developer MkDocs and the Held comparison at each completed coherent
  stage.

---

### Task 1: Exact neck topology projection

**Status (2026-07-29): complete and published at `d289f79`.**

**Files**

- Modify: `src/segment_site_mat_bundle.h`
- Modify: `src/segment_site_mat_bundle.cpp`
- Modify: `src/segment_site_neck_evidence.h`
- Modify: `src/segment_site_neck_evidence.cpp`
- Modify: `src/segment_site_neck_evidence_bytes.cpp`
- Modify: `src/medial_axis_2.cpp`
- Modify: `src/compas_cgal/_medial_axis_2.pyi`
- Modify: `src/compas_cgal/adaptive/neck.py`
- Create: `src/compas_cgal/adaptive/neck_topology.py`
- Test: `tests/native/task9_neck_evidence_gate.cpp`
- Test: `tests/native/task9_mat_numeric_table_gate.cpp`
- Test: `tests/adaptive/test_neck.py`

**Interfaces**

- Additive `SegmentSiteMedialAxis.neck_location_tags`.
- Additive `SegmentSiteMedialAxis.neck_location_edge_ids`.
- Additive `SegmentSiteMedialAxis.neck_location_node_ids`.
- Additive `SegmentSiteMedialAxis.neck_parameter_root_ids`.
- Additive `SegmentSiteMedialAxis.neck_cut_edge_partitions`.
- Typed four-variant Python neck locus.
- `NeckSide` and complete `ClassifiedNeck.sides`.

- [x] **Step 1: Write native RED retention tests**

Require exact evidence construction to retain all location fields and nested
partitions after `SegmentSiteMatBundle2::build()`. Native synthetic gates cover
all four location variants and malformed endpoint projection. The L gate
proves both plateau necks have three canonical partitions and that reversal
produces the same records.

- [x] **Step 2: Write Python RED projection tests**

Require the generic native projection to cover all four synthetic location
variants, then require the production L owner to cross the binding with:

- two plateau loci;
- three sides per neck;
- side unions equal the established projected cut unions;
- side IDs and complete `ClassifiedNeck` bytes survive input reversal; and
- deleted, reordered, duplicated, unknown-edge, and cross-wired partitions
  fail `InvalidNeckEvidenceError`.

- [x] **Step 3: Implement retained native topology**

Derive binding fields directly from each retained `MatNeckEvidenceV1`. Keep
the exact evidence records or their complete derived topology inside
`SegmentSiteMatBundle2`; never decode its canonical bytes.

- [x] **Step 4: Implement typed Python topology**

Construct the closed neck-locus union, canonical `NeckSide` values, and
cross-validated `ClassifiedNeck` fields. Keep the existing cut union as an
additive compatibility fact during validation; do not use it as traversal
authority.

- [x] **Step 5: GREEN and publish the coherent stage**

```bash
pixi run task9-mat-compile-gate
pixi run pytest tests/adaptive/test_neck.py -n auto --testmon -q
pixi run lint
pixi run types-adaptive
pixi run docs
```

Commit: `feat(mat): expose exact neck sides`

---

### Task 2: Directed global traversal ledger

**Status (2026-07-29): complete and published at `6ac2050`.**

**Files**

- Create: `src/compas_cgal/adaptive/traversal.py`
- Create: `src/compas_cgal/adaptive/traversal_graph.py`
- Modify: `src/compas_cgal/adaptive/candidates.py`
- Modify: `src/compas_cgal/adaptive/errors.py`
- Test: `tests/adaptive/test_traversal.py`
- Modify: `tests/adaptive/typecheck/consumer_contract.py`

**Interfaces**

- `DirectedEdgeCursor.build(...) -> DirectedEdgeCursor`
- `CausalNeckTransit.build(...) -> CausalNeckTransit`
- `MatTraversalState.seed(...) -> MatTraversalState`
- `MatTraversalState.advance(candidate) -> MatTraversalState`
- `MatTraversalState.activate_next(...) -> MatTraversalState`
- `MatTraversalState.require_terminal() -> None`
- `TraversalSampleIndex.build(...) -> TraversalSampleIndex`
- `TraversalGraph.from_axis(...) -> TraversalGraph`

- [x] **Step 1: Write RED state contracts**

Test one cursor per edge, canonical state identity, exact owner binding,
single-cursor advancement, stale cursor rejection, terminal-edge rejection,
alias-free immutable transitions, and deterministic component/branch order.

- [x] **Step 2: Write RED graph-route contracts**

Use exact synthetic chain, branch, cycle, and multi-component topologies.
Require transitions to use shared node IDs, retain both cycle incidences, and
visit every edge. Mutate a node ID while preserving reporting coordinates and
require loud failure.

- [x] **Step 3: Write RED causal-side contracts**

Require unique side initialization, a three-partition plateau transition,
canonical forward/reverse orientation, passage-state lookup, and exact
NoNeckScope outside a transit. Reject union-only inference, ambiguous sides,
overlapping active necks, and a relabelled partition.

- [x] **Step 4: Implement immutable traversal**

Keep route discovery, per-edge cursors, visited incidences, and neck sides in
one graph-lifetime value. Hash the MAT certificate digest plus complete
canonical state. Do not include mutable stock or coverage.

- [x] **Step 5: GREEN and document**

```bash
pixi run pytest tests/adaptive/test_traversal.py -n auto --testmon -q
pixi run lint
pixi run types-adaptive
pixi run docs
```

Commit: `feat(adaptive): add causal traversal`

---

### Task 3: Bidirectional finite-lattice spans

**Status (2026-07-29): complete locally; publication pending.**

**Files**

- Modify: `src/compas_cgal/adaptive/candidates.py`
- Modify: `src/compas_cgal/adaptive/replay.py`
- Test: `tests/adaptive/test_candidates.py`
- Test: `tests/adaptive/test_replay.py`
- Test: `tests/adaptive/test_traversal.py`

- [x] **Step 1: Write RED reverse-span tests**

Mirror a real L line edge and P-S parabola span. Require reverse progress,
geometry, native/derived cursor identity, deterministic enumeration, terminal
limit, and forward/reverse identity separation. Reject a direction reversal
inside one derived continuation.

- [x] **Step 2: Implement inferred span direction**

Allow ordered native cursor pairs in either ordinal direction. Retain the
direction in `DerivedCandidateCursor`, compute the next legal limit ordinal,
and leave positive progress and MATHSM geometry unchanged.

- [x] **Step 3: Extend fresh unique reconstruction**

Replay searches the directed window owned by traversal state. It must uniquely
reconstruct old forward fixtures and new reverse fixtures with no
orientation-based inference.

- [x] **Step 4: GREEN and document**

```bash
pixi run pytest tests/adaptive/test_candidates.py \
  tests/adaptive/test_replay.py \
  tests/adaptive/test_traversal.py -n auto --testmon -q
pixi run lint
pixi run types-adaptive
pixi run docs
```

Commit: `feat(adaptive): traverse mat both ways`

---

### Task 4: Entry-circle bootstrap

**Files**

- Create: `src/compas_cgal/adaptive/bootstrap.py`
- Modify: `src/compas_cgal/adaptive/errors.py`
- Test: `tests/adaptive/test_acceptance.py`

**Interfaces**

- `InitialCandidateEvaluator.build(...) -> InitialCandidateEvaluator`
- `evaluate(traversal, candidate) -> InitialCandidateTransaction`
- `commit(traversal, transaction) -> tuple[GenerationState,
  MatTraversalState]`

- [ ] **Step 1: Write RED launch tests**

Require exact phase equality with entry center, full sweep inside the qualified
entry, design containment, motion certification against post-entry stock,
certify-before-deplete chronology, coverage addition, one cursor advance,
causal passage advance only on commit, and deterministic evidence.

- [ ] **Step 2: Write RED atomicity and mutation tests**

Reject a zero-length synthetic link, circle outside the entry disk, foreign
entry/domain/tool/cut depth, stale traversal parent, deplete-before-certify,
cross-wired circle witness, and a launch that advances two cursors.

- [ ] **Step 3: Implement isolated launch and independent commit**

Build pristine stock, apply entry depletion once, build empty coverage, then
evaluate the candidate circle on forks. Commit repeats the complete launch and
requires byte-identical transaction evidence.

- [ ] **Step 4: GREEN and document**

```bash
pixi run pytest tests/adaptive/test_acceptance.py -n auto --testmon -q
pixi run lint
pixi run types-adaptive
pixi run docs
```

Commit: `feat(adaptive): certify entry launch`

---

### Task 5: Global continuation and finite search

**Files**

- Modify: `src/compas_cgal/adaptive/transaction.py`
- Create: `src/compas_cgal/adaptive/generator.py`
- Modify: `src/compas_cgal/adaptive/errors.py`
- Test: `tests/adaptive/test_transaction.py`
- Test: `tests/adaptive/test_traversal.py`
- Test: `tests/adaptive/test_acceptance.py`

**Interfaces**

- One internal candidate-trial engine accepting an authenticated active cursor.
- `TraversalCommit.build(...) -> TraversalCommit`
- `generate_exact_adaptive(...) -> GenerationResult`

- [ ] **Step 1: RED branch-switch transaction tests**

Require the same Task 12 proof order when the active global edge differs from
the last physical edge. Bind physical and traversal parents, advance exactly
one cursor, and reject stale/cross-wired parents. Existing same-edge Task 12
tests stay unchanged and GREEN.

- [ ] **Step 2: RED finite-search tests**

Require invariant ordered search, first-feasible winner, no lower-ranked
evaluation afterward, exact continuation after gouge/cap failures, immediate
propagation of unresolved events, and named exhaustion errors. Count one
candidate-family build per attempted span.

- [ ] **Step 3: Implement one deciding trial engine**

Refactor only enough to let both Task 12 and the traversal consumer supply an
authenticated cursor. Keep exactly one containment/certification/depletion/
coverage/witness path.

- [ ] **Step 4: Implement route orchestration**

Seed with the launch transaction, activate deterministic directed edges,
derive scope/cap from causal side state, enumerate each forward window once,
commit the first feasible candidate, and continue until every component is
terminal. No recursion or iteration bound substitutes for the well-founded
cursor measure.

- [ ] **Step 5: GREEN and document**

```bash
pixi run pytest tests/adaptive/test_transaction.py \
  tests/adaptive/test_traversal.py \
  tests/adaptive/test_acceptance.py -n auto --testmon -q
pixi run lint
pixi run types-adaptive
pixi run docs
```

Commit: `feat(adaptive): generate covered traversal`

---

### Task 6: Tractable fixture, terminal coverage, and fresh replay

**Files**

- Create: `tests/adaptive/fixtures/tractable_pocket.json`
- Modify: `src/compas_cgal/adaptive/generator.py`
- Modify: `src/compas_cgal/adaptive/replay.py`
- Modify: `src/compas_cgal/adaptive/replay_trace.py`
- Test: `tests/adaptive/test_acceptance.py`
- Test: `tests/adaptive/test_replay.py`

- [ ] **Step 1: Commit exact fixture inputs**

Record polygon, holes, cut plane, tool, entry evidence identity, MAT sampling,
candidate, neck, depletion, traversal, and cut-direction policies. No expected
certificate bytes are hand-authored.

- [ ] **Step 2: RED complete acceptance**

Require approach/plunge, at least one segment and two circles, complete witness
bijection, terminal traversal, exact empty reachable residual, separately
bound unreachable residual, and fresh causal replay with zero cap violations.

- [ ] **Step 3: RED mutations**

Kill empty, entry-only, first-circle-only, dropped-final-branch, relabelled
side, reversed transit, stale traversal parent, nonterminal cursor, nonempty
residual, and certify-after-deplete histories.

- [ ] **Step 4: Implement terminal seal**

After traversal terminality, call `CoverageLedger.require_complete()`, rebuild
the complete traversal and physical state through fresh replay, re-derive each
causal neck scope, and build content-addressed `GenerationResult`.

- [ ] **Step 5: GREEN**

```bash
pixi run pytest tests/adaptive/test_acceptance.py \
  tests/adaptive/test_replay.py -n auto --testmon -q
```

Commit: `feat(adaptive): seal covered generation`

---

### Task 7: Verification, performance baseline, docs, and publication

**Files**

- Modify: `docs/segment_site_mat.md`
- Modify: `docs/continuous_engagement.md`
- Modify:
  `docs/superpowers/plans/2026-07-24-exact-certified-adaptive-trochoidal-phase1.md`
- Modify: `docs/superpowers/plans/2026-07-28-causal-mat-traversal.md`

- [ ] **Step 1: Structural performance counters**

Gate one MAT/inventory build, one candidate-family build per attempted span,
no post-acceptance lower-ranked trials, one winner replay per commit, and one
terminal fresh replay.

- [ ] **Step 2: Bounded runtime measurement**

Measure warm Release medians for generation and replay on the committed
fixture. Report build/oracle/candidate counts and wall time. Do not call this
Held parity; matched Fig. 5 measurement remains Task 16.

- [ ] **Step 3: Full verification**

```bash
pixi run format-adaptive
pixi run lint
pixi run types-adaptive
pixi run affected
pixi run pytest tests/adaptive -n auto -q
pixi run docs
git diff --check
```

- [ ] **Step 4: ETH review and regression repair**

Audit exact authority, ownership, error models, file responsibility, strict
typing, proof chronology, mutation coverage, and documentation accuracy.
Repair every finding before completion.

- [ ] **Step 5: Update contemporaneous MkDocs**

Document implemented behavior, topology, three-partition insight, causal scope,
entry bootstrap, traversal/coverage/replay evidence, measured runtime,
limitations, and the Held–Pfeiffer comparison. Mark only verified checklist
items complete.

- [ ] **Step 6: Commit, cancel superseded CI, push, and verify**

Commit each remaining coherent stage as Jelle Feringa. Before every push, list
the branch's runs with explicit repository, cancel every `in_progress` run,
push without force, and require remote branch SHA to equal local `HEAD`.
