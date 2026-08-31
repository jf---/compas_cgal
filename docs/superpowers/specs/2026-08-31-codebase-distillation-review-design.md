# Codebase Distillation Review Design

> **status: approved** — approved in conversation 2026-08-31.

## Decision

Review every scoped first-party code variant three complete times before
proposing any condensation. Preserve demonstrated mathematical and engineering
capability; reduce duplicate authority, accidental orchestration, and machinery
that cannot justify itself through a consumer, invariant, counterexample, or
measured outcome.

The review produces a surgery proposal. It does not alter production code,
tests, build configuration, branches, or worktrees. Any later condensation is a
separate add-and-validate programme requiring explicit approval before removal.

## Goal

Produce an evidence-backed codebase map in which:

- every exceptional capability has a named preserved home;
- every deciding domain truth has one proposed authority;
- every competing implementation has an explicit relationship to that
  authority;
- every proposed condensation states the valuable nucleus and the contracts
  that must survive;
- every proposed removal has a loss argument, falsifier, and consumer oracle;
  and
- the retained capabilities compose into one honest route toward a matched
  Held–Pfeiffer result and downstream machine evidence.

The target is fewer independent concepts and competing authorities, not a line
count or file-count reduction.

## Non-goals

This review does not:

- implement, repair, delete, move, or rename product code;
- choose the future generator algorithm before all three readings;
- treat Git history as a substitute for extracting and naming valuable ideas;
- call unconsumed code worthless merely because it is not currently composed;
- weaken exact predicates, proof boundaries, or named unresolved states for a
  simpler narrative;
- build a new identity, fingerprinting, receipt, or review-pinning system;
- claim superiority over Held from internal sophistication; or
- require machine evidence for a mathematical capability whose claim ends at a
  lower consumer boundary.

## Review universe

### Included

The scope is the union of differing first-party textual code variants across
all ratified active local branches, remote-only active branches without a local
equivalent, and dirty or untracked code in live worktrees:

- C++ sources and headers;
- Python sources and binding stubs;
- tests, fixtures, and executable examples;
- benchmarks and measurement code;
- scripts and developer tools;
- CMake, Pixi, CI, task, and code-generating configuration;
- schemas consumed by executable code; and
- behavior-bearing code snippets in documentation.

The initial lineage candidates are:

- `main`;
- `jf/toolpath-redesign`;
- `codex/held-reference-corpus`;
- `codex/auditor-convergence-sdd`;
- `codex/sdd-coherence`;
- `codex/exact-certified-adaptive-phase1-t9-zero-guide`;
- `perf/exact-rational-representation`; and
- any additional branch found by the scope ratification whose unique
  first-party code is not represented above.

Branch age and naming never establish relevance. Exact textual duplicates may
share one review row only after direct Git comparison. Textually different
variants remain separate until Pass 3, even if they appear semantically
equivalent.

### Excluded

- compiled objects, wheels, caches, and build directories;
- vendored third-party source with no local modification;
- generated images, PDFs, and benchmark result dumps;
- dependency lockfiles;
- prose without executable snippets; and
- every historical revision of every file.

Named historical lineages receive a separate deleted-capability salvage sweep
when Pass 1 finds evidence that a capability may exist only in history. The
review never implies that all historical revisions were read.

## Review stability

Reviewers read the same worktree state. Edits pause while a pass is active.
Ordinary branch names, paths, status output, and Git comparisons document the
state; they do not become an application-level review identity.

If a scoped variant changes during a pass:

1. mark that variant's affected pass record stale;
2. identify dependent capability and finding entries;
3. reread the changed variant in the affected pass;
4. revisit only dependent later-pass judgments; and
5. record the refreshed evidence in the progress ledger.

No source worktree or branch is modified, rebased, merged, deleted, or cleaned
during the review. User-owned untracked content is preserved.

## Preservation unit: capability

Files are not the preservation unit. One file may contain several capabilities;
one capability may span several files and branches.

Each capability entry records:

- ordinary capability ID and name;
- locations and variants;
- mathematical claim;
- product claim, if any;
- invariants and exactness boundary;
- named failure states;
- consumers and public boundary;
- counterexamples and oracles;
- measured value and workload;
- uniqueness and competing implementations;
- replacement cost;
- Held relevance;
- Buchli relevance;
- proposed nucleus;
- surrounding machinery;
- missing evidence; and
- provisional disposition after Pass 3 only.

## The three readings

### Pass 1 — capability archaeology

Read every scoped variant bottom-up in forward dependency order. Start from
pathological tests and hard fixtures, trace into the deciding implementation,
then trace outward to consumers.

For every meaningful symbol or cohesive algorithm, answer:

- What observable capability exists?
- What theorem, invariant, or engineering idea does it embody?
- Which rare topology or numerical case does it preserve?
- Which tests can distinguish its presence from its absence?
- Which measurements establish its cost or value?
- Who consumes it today?
- What would be expensive or impossible to reconstruct?
- What is explicitly incomplete?

No `KEEP`, `CONDENSE`, `ABSORB`, `QUARANTINE`, `REMOVE`, or `UNKNOWN`
classification is allowed in Pass 1.

Output: the complete capability ledger and the first manifest coverage column.

### Pass 2 — authority and product truth

Read every scoped variant again, top-down in reverse dependency order. Trace
outputs back toward kernels and identify who decides each domain truth.

For every capability family, answer:

- Which implementation is authoritative, and for which exact claim?
- Do multiple layers decide the same fact?
- Is complexity mathematical or created by ownership and lifecycle machinery?
- Does the proof or validation apparatus have a consumer?
- Is a valuable kernel hidden behind accidental orchestration?
- Is an assertion tested only through a copy of its own logic?
- What smaller boundary could carry the same capability honestly?

#### Held lens

Trace the complete competitive path:

```text
pocket -> MAT and clearance -> candidate motion -> stock mutation
       -> engagement decision -> coverage -> motion operations -> G-code
       -> independent parse/backplot -> comparable machining outcome
```

Classify every break as absent, present but uncomposed, hidden behind excessive
machinery, asserted but unverified, or mathematically verified without a
product result. Stronger local mathematics never implies a stronger completed
toolpath.

#### Buchli lens

Ask whether local capability survives as a complete observable system:

- Does output survive downstream controller semantics?
- Are entry, retraction, continuity, feed, acceleration, and latency represented
  at the claim boundary?
- Does validation observe the same artifact consumed downstream?
- Is timing measured end to end?
- Can failures be reproduced and diagnosed?
- Does the evidence stop at the correct layer?

Evidence escalates proportionally:

```text
contract -> independent geometry -> parsed G-code
         -> controller simulation -> machine execution
```

Pass 2 records authorities, duplicates, essential complexity, accidental
complexity, and candidate condensation seams. It does not assign final
dispositions.

### Pass 3 — condensation falsification

Read every scoped variant a third time, grouped by semantic family across
branches. Prosecute every proposed simplification on behalf of the exceptional
code.

For each candidate reduction, state:

> Capability X can be preserved by nucleus Y, consumer contracts Z, and
> counterexamples Q; machinery M is not semantically required.

Then try to falsify it:

- find hidden consumers and ABI commitments;
- search branch-only improvements;
- compare named failure states and unresolved distinctions;
- run exactness, differential, metamorphic, counterexample, and bounded
  performance oracles;
- test whether two owners have different lifetimes or evolution rates;
- test whether apparent equivalence fails on a legitimate input; and
- state what valuable option would be lost.

Only Pass 3 may assign a disposition.

## Dispositions

`UNKNOWN` is the safe default.

| Disposition | Required evidence |
| --- | --- |
| `KEEP` | Unique validated capability or expensive mathematical result with a coherent owner |
| `CONDENSE` | Valuable nucleus is separable and named contracts preserve its semantics |
| `ABSORB` | A natural existing owner has matching semantics, lifetime, consumers, and failure distinctions |
| `QUARANTINE` | Valuable unresolved research should remain recoverable but must not own a production claim |
| `REMOVE` | No unique capability, consumer, counterexample, failure distinction, or measured advantage survives comparison |
| `UNKNOWN` | A named missing oracle or experiment prevents a safe judgment |

Size, age, internal status, lack of current consumers, or self-contained tests
are insufficient evidence for `REMOVE`.

## Semantic families

Every variant belongs to one primary review family while cross-family
dependencies remain explicit:

1. general COMPAS/CGAL bindings and public compatibility;
2. motion vocabulary and legacy toolpath generation;
3. stock, depletion, containment, and coverage;
4. point, guarded, and continuous engagement geometry;
5. medial axis, segment-site geometry, reachability, and neck topology;
6. adaptive candidates, traversal, transactions, replay, and generation;
7. benchmarks, quality metrics, comparison, and visualization;
8. build, bindings, typing, CI, scripts, tools, and executable documentation.

The family map may split a family when two best reviewers detect disjoint
failure surfaces. A split is recorded as a scope correction, not silently
introduced work.

## Adversarial panel

Named personas are internal public-record-conditioned instruments, never
attributions. Findings are restated in the project's voice and accepted only
through evidence.

Standing global lenses:

- Martin Held — direct CAM comparator: MAT/Voronoi fidelity, engagement
  regulation, stock awareness, coverage, path quality, and fair comparison;
- Jonas Buchli — end-to-end systems closure, observability, resource limits,
  controller contact, and proportional evidence;
- John Ousterhout — cognitive load, change amplification, shallow modules, and
  dependencies whose cost exceeds their hidden capability;
- Sylvain Pion — exact computation, predicate/construction discipline,
  degeneracy, robustness, and performance without weakened truth;
- Peter Smid — controller-explicit G-code, operator diagnosability, and
  downstream machine cost; and
- one null reviewer with no persona or domain briefing.

Before dispatch, `adversarial-panel-construction` ratifies the failure surface,
public anchors, opponent, cost-bearer, resource role, interface role, and null
control. Thin or confabulated named lenses become de-named commitment lenses.

Reviewers receive the artifact and finding contract but are not told which
failure classes they are expected to own. Capacity-limited dispatch occurs in
blind waves; later reviewers receive no earlier output. The panel generates
verification proposals. It never votes on dispositions.

After Pass 1, each semantic family gets its own panel against that family's
failure surface. Held and Buchli remain standing only where their detection
functions are relevant; family specialists may replace other seats.

## Oracles

### Baseline

Each executable lineage gets a truthful, sequential baseline using its own
policy-declared Pixi tasks. Known red product criteria remain red and are
recorded; they are never laundered into regression success. Native editable
builds never run concurrently.

### Liveness

Before a green oracle can exonerate code, a reversible seeded defect in an
isolated audit worktree must demonstrate that the exact gate can detect its
claimed failure class. Seeded source changes are never committed and are fully
reversed before the readings begin.

At minimum validate detection of:

- a weakened exact decision;
- a lost frame or unit distinction;
- malformed motion semantics;
- incorrect stock depletion;
- reachable residual stock;
- cap-violating motion; and
- a planner/consumer disagreement.

An absent G-code round-trip, controller simulation, or machine test is recorded
as a missing oracle, not replaced by an inferred product claim.

### Routing

- C0: compiler, type checker, linter, existing tests — always verify.
- C1: focused characterization, differential test, property test, or bounded
  benchmark — verify when plausible or high-impact.
- C2: sanitizer, fuzzer, stress, controller, large benchmark, or machine run —
  verify above the severity threshold; otherwise queue explicitly.
- C3: architecture, ownership, and explanatory complexity — Jelle adjudicates.

One reviewer finding is sufficient to run a cheap oracle. Agreement signals a
missing characterization test; it does not establish truth.

## Review artifacts

The review creates exactly five durable artifacts:

1. `.superpowers/sdd/2026-08-31-codebase-distillation-review/progress.md`
   — one active step, last evidence, drift, blockers, and next command;
2. `docs/superpowers/state/2026-08-31-distillation-manifest.tsv`
   — ordinary variant labels, family, scope decision, and three pass states;
3. `docs/superpowers/state/2026-08-31-distillation-capabilities.md`
   — the capability ledger;
4. `docs/superpowers/state/2026-08-31-distillation-findings.tsv`
   — adversarial findings, oracles, falsifiers, and adjudication; and
5. `docs/superpowers/state/2026-08-31-distillation-surgery.md`
   — final authority graph and disposition proposal.

A small validator checks artifact shape, allowed states, referential integrity,
three-pass coverage, and that no final reduction lacks its required fields. It
does not decide code value or create another evidence authority.

## Finding contract

Every finding records:

```text
ID | PASS | FAMILY | VARIANT | LOCATION | QUOTE | CAPABILITY | ISSUE
CONSUMER | LOSS_IF_CHANGED | PROPOSED_CONDENSATION | ORACLE | COST
FALSIFIER | EVIDENCE | STATUS
```

Allowed statuses are `proposed`, `verified`, `rejected`, `queued-c2`, and
`jelle-c3`.

## Acceptance

The review is complete only when:

- every scoped textual variant has three valid pass records;
- no changed worktree state remains unread;
- every exceptional capability has a named retained, absorbed, or quarantined
  home;
- every deciding domain truth has one proposed authority;
- every competing implementation has an explicit relationship to that
  authority;
- every `CONDENSE` names its nucleus and preserved consumer contracts;
- every `ABSORB` and `REMOVE` has a loss argument, oracle, and falsifier;
- all C0/C1 findings are verified or rejected;
- every remaining C2 is listed with severity and cost;
- all C3 judgments are presented to Jelle;
- a fresh adversarial synthesis round yields no new verified disposition
  change;
- strict documentation and artifact validation pass; and
- Jelle approves the surgery ledger.

Completion authorizes planning the surgery. It does not authorize product-code
changes or removal.
