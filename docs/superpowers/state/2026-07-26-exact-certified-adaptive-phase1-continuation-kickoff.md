# Exact-certified adaptive Phase 1 continuation — T9 exact MAT

I'm starting **Task 9 continuation — exact segment-site Voronoi medial axis**
of the **exact-certified adaptive Phase 1** project
(`/private/tmp/compas_cgal_prs-exact-certified-adaptive-phase1-t9`, branch
`codex/exact-certified-adaptive-phase1-t9`).
`T9-endpoint-binding-237ff85` is the last published, independently reviewed
hard checkpoint. Local commit `2cab98c` is a green but unreviewed point-graph
checkpoint; three uncommitted files form a failed cocircular experiment.

## 1. What I want you to do

Resume the governing written plan with
`superpowers:executing-plans` and
`superpowers:subagent-driven-development`. Use one same-model/same-effort
implementer per shared file set and separate read-only reviewers.

First resolve the three-file failed experiment without losing local commit
`2cab98c`. Independently review and publish that focused point-only production
slice if it survives. Then continue the remaining Task 9 verticals in strict
proof-dependency order. Keep Tasks 10–16 frozen until the complete Task 9
contract is green. Do not advance Task 8 while its Task 6 dependency is open.

Each implementation slice has a maximum wall-clock budget of 30 minutes.
The budget never permits breaking dependency order, weakening exactness, or
exposing a falsely generic API.

## 2. Read first, in dependency order

1. `/private/tmp/compas_cgal_prs-exact-certified-adaptive-phase1-t9/CLAUDE.md`
   — read completely; binding repository policy.
2. `/private/tmp/compas_cgal_prs-exact-certified-adaptive-phase1/docs/superpowers/state/2026-07-26-exact-certified-adaptive-phase1-resume.md`
   — current multi-worktree state, accepted evidence, failed experiment, and
   immediate recovery protocol.
3. `/private/tmp/compas_cgal_prs-exact-certified-adaptive-phase1/.superpowers/sdd/2026-07-24-exact-certified-adaptive-trochoidal-phase1/progress.md`
   — authoritative SDD ledger and dependency status.
4. `docs/superpowers/plans/2026-07-24-exact-certified-adaptive-trochoidal-phase1.md`
   — read Task 9 completely, plus the Task 10 boundary. The plan is governing.
5. `docs/exactness.md` and `docs/dev/exact-cgal-idioms.md` — exactness doctrine
   and mandatory idiomatic CGAL route.
6. `src/segment_site_mat.{h,cpp}`,
   `src/segment_site_parameterization.{h,cpp}`,
   `src/segment_site_clipping.{h,cpp}`,
   `src/segment_site_endpoint_binding.{h,cpp}`, and
   `src/segment_site_provenance.{h,cpp}` — current native substrate.
7. `tests/native/task9_mat_compile_gate.cpp`,
   `tests/native/task9_point_graph_gate.cpp`, and
   `tests/native/task9_segment_limiter_gate.cpp` — current exact contract
   gates.
8. Official locked CGAL 6.0.1 headers under
   `external/cgal/include/CGAL/`, especially
   `Segment_Delaunay_graph_2.h`,
   `Segment_Delaunay_graph_storage_traits_with_info_2.h`,
   `Segment_Delaunay_graph_adaptation_policies_2.h`, and
   `Voronoi_diagram_2.h`. Read the authoritative backend before designing an
   adaptor path.

Do not rely on prior chat or agent claims. Reproduce every current-state claim
from the files, branch tips, and gates above.

## 3. Boundary contract (= your input)

Published `237ff85` supplies exact live parabola endpoint binding:

```cpp
MatParameterEndpoint2 bind_point_limiter_parabola_endpoint(
    const MatExactPointSiteSource2& focus,
    const MatExactOpenSegmentSource2& segment,
    const MatExactPointSiteSource2& segment_source,
    const MatExactPointSiteSource2& segment_target,
    const MatExactPointSiteSource2& limiter,
    const SegmentSiteVoronoi2& voronoi,
    const SegmentSiteVoronoi2::Halfedge_handle& halfedge);

MatParameterEndpoint2 bind_segment_limiter_parabola_endpoint(
    const MatExactPointSiteSource2& focus,
    const MatExactOpenSegmentSource2& segment,
    const MatExactPointSiteSource2& segment_source,
    const MatExactPointSiteSource2& segment_target,
    const MatExactOpenSegmentSource2& limiter,
    const MatExactPointSiteSource2& limiter_source,
    const MatExactPointSiteSource2& limiter_target,
    const SegmentSiteVoronoi2& voronoi,
    const SegmentSiteVoronoi2::Halfedge_handle& halfedge);
```

Local commit `2cab98c` adds this point-only production boundary:

```cpp
struct NormalizedPointSource2 {
    std::string stable_site_id;
    CORE::BigRat x;
    CORE::BigRat y;
};

MatExactGraph2 exact_point_site_graph(
    const std::vector<NormalizedPointSource2>& points,
    const MatDomainPolygonWithHoles2& domain,
    const CORE::BigRat& radius_squared);
```

At committed `2cab98c`, the native build/run gate passes. The current
uncommitted tree adds `DegeneratePointSiteTopologyError` plus a cocircular
square expectation; it builds but the native executable exits 1. That
experiment is not part of the accepted boundary.

The current graph support matrix is:

```text
P-P LINE/RAY/SEGMENT:
  exact rational parameterization, D clipping, radius clipping, graph append

nonincident P-S PARABOLA:
  exact source parameterization and exact endpoint binders exist
  but binders, algebraic bounds, and real-radius clipping are not wired into
  production graph construction

S-S:
  no production feature-domain, endpoint, clearance, or graph path
```

## 4. Deliverables

Complete only plan Task 9 in this continuation:

1. Diagnose the current cocircular gate failure. Either restore clean
   `2cab98c` or commit a proved fix-forward; never tune a fixture to trigger a
   preconceived outcome.
2. Independently review the point-only production graph, fix all
   Critical/Important findings, commit focused, push to the user fork, and
   verify the remote SHA.
3. Wire nonincident P-S halfedge endpoint ownership into graph construction
   using the accepted binders.
4. Add algebraic clipping bounds and the true exact
   `distance_to_defining_site^2 >= radius_squared` parabola predicate.
5. Implement S-S feature domains and one unified segment-Delaunay/Voronoi
   adaptor traversal with stable site info or a proved-bijective exact
   post-map.
6. Build degeneracy-normalized node CSR from every incident halfedge and
   preserve every exact endpoint provenance source.
7. Create `src/segment_site_neck.{h,cpp}` with the four exact
   `NeckEvidenceV1` variants and deterministic separating-cut ownership.
8. Create `src/segment_site_mat_sampling.{h,cpp}` for proposal-only sampling;
   exact topology and certificates must not depend on samples.
9. Create `src/medial_axis_2.cpp`,
   `src/compas_cgal/_medial_axis_2.pyi`, and
   `tests/adaptive/test_medial_axis.py`; return the fixed 20-field tuple from
   the governing plan and bind the Task 3 `C_r` digest.

Preserve all existing spike paths until the new unified path is independently
validated. Removal or replacement of an old path requires user permission.

## 5. Exit gates

The governing plan's Task 9 GREEN gates are quoted verbatim:

> ```bash
> pixi run format-adaptive
> pixi run lint
> pixi run pytest tests/adaptive/test_medial_axis.py -n auto --testmon -q
> pixi run pytest tests/test_toolpath.py -n auto -q
> ```
>
> Commit: `feat(mat): add exact segment Voronoi graph`

Before those final gates, every native slice must also pass:

```bash
pixi run task9-mat-compile-gate
git diff --check
```

Run Ruff and strict mypy for every Python slice. Every pytest invocation uses
`-n auto`; affected tests use `--testmon`.

If an exact primitive, adaptor contract, topology identity, or performance
gate slips, write
`docs/superpowers/state/2026-07-26-task9-gate-analysis.md` first. State the
root cause, exact counterexample, and options before implementing a changed
course. Never substitute sampling or approximation.

## 6. Decoupling rule

Task 9 depends only on accepted Tasks 3 and 5, so continue it independently of
the Task 6 fixture decision. Task 10 consumes the complete Task 9 certificate
and is frozen.

The T9 branch currently diverges from the primary Phase 1 branch. Do not merge,
rebase, or cherry-pick while the three-file experiment is unresolved. After
T9 is green and reviewed, plan the linear integration explicitly; never merge
implicitly.

Native fixtures must be bounded, exact, and analytic. Proposal samples may
change under refinement, but topology, provenance, neck ownership, and
certificate bytes may not.

## 7. Load-bearing conventions

- Exact-kernel decisions use CGAL predicates and exact number types; no
  epsilon, deflation, ULP nudge, finite probes, or `to_double`.
- Task 9 is a CGAL adaptor, not a hand-written replacement MAT.
- Never use `Parabola_2::side_of_parabola`,
  `Parabola_segment_2::compute_k`, `generate_points`, or streamed parabolas
  for topology, clipping, identity, or evidence.
- A quartic root remains an `AlgebraicRootIdV1`; do not coerce it into a
  `MatKernel::Point_2`.
- Use stable site/root/curve identities only. Handles, addresses, iterator
  order, and doubles never enter identity.
- Fail loudly with a named exception for every unsupported or malformed
  producer state. No fallback behavior or parallel calling conventions.
- Use Pixi exclusively. Never invoke pip, conda, Homebrew, or a venv.
- Never edit a reference test merely to make it pass. Never skip or xfail.
- Use `apply_patch` for edits. Preserve old paths until the replacement passes
  and the user permits removal.
- Commit focused slices proactively on this agent branch with author and
  committer `Jelle Feringa <jelleferinga@gmail.com>`, concise message, no
  attribution.
- Never mutate `main`, `master`, `jf/toolpath-redesign`, or the protected
  checkout.
- Before a push, cancel in-progress runs on the exact branch with explicit
  `gh --repo jf---/compas_cgal`; verify the pushed SHA with
  `git ls-remote --heads origin`.

## 8. Environment booby traps

- Locked native stack: CGAL 6.0.1, repository Boost/Eigen sources, CORE with
  Boost backend, GMP-disabled flags.
- `cmake -S` for the Task 9 target currently costs roughly 80–92 seconds;
  the native executable runs in roughly 1.1 seconds. Do not misdiagnose
  configure latency as algorithmic runtime.
- Do not run concurrent editable/native builds.
- `SegmentSiteVoronoi2::Halfedge` access is sided:
  call `left()` only when `has_source()` and `right()` only when
  `has_target()`. Violating this can terminate a Release build.
- Direct `CGAL::is_square(CORE::BigRat)` is not supported in the locked
  version. Decompose with `Fraction_traits<CORE::BigRat>::Decompose`, then use
  `CGAL::is_square` on exact integer numerator and denominator.
- Square-free factor coefficients plus factor-local root ordinal are the
  authoritative algebraic identity. Never use whole-polynomial ordinal after
  factorization.
- The current native gate combines many boolean predicates and reports only
  process exit status. Isolate the new cocircular predicate before drawing a
  conclusion from exit 1.
- GitHub HTTPS API calls timed out twice in the prior session. SSH remote
  checks and pushes worked. Retry boundedly; do not omit mandatory CI-state
  checks when a remote branch already exists.

## 9. Suggested file structure

```text
compas_cgal_prs-exact-certified-adaptive-phase1-t9/
├── CLAUDE.md                                      # READ
├── docs/
│   ├── exactness.md                               # READ
│   ├── dev/exact-cgal-idioms.md                  # READ
│   └── superpowers/
│       ├── plans/2026-07-24-exact-certified-adaptive-trochoidal-phase1.md  # READ
│       └── state/2026-07-26-exact-certified-adaptive-phase1-continuation-kickoff.md  # THIS FILE
├── src/
│   ├── segment_site_mat.{h,cpp}                  # MODIFY
│   ├── segment_site_parameterization.{h,cpp}     # MODIFY
│   ├── segment_site_clipping.{h,cpp}             # MODIFY
│   ├── segment_site_endpoint_binding.{h,cpp}     # MODIFY
│   ├── segment_site_provenance.{h,cpp}           # MODIFY
│   ├── segment_site_mat_sampling.{h,cpp}         # CREATE
│   ├── segment_site_neck.{h,cpp}                 # CREATE
│   ├── medial_axis_2.cpp                         # CREATE
│   └── compas_cgal/_medial_axis_2.pyi           # CREATE
├── tests/
│   ├── native/task9_mat_compile_gate.cpp         # MODIFY
│   ├── native/task9_point_graph_gate.cpp         # MODIFY
│   ├── native/task9_segment_limiter_gate.cpp     # MODIFY
│   └── adaptive/test_medial_axis.py               # CREATE
└── CMakeLists.txt                                 # MODIFY
```

## 10. Risk + budget

The governing plan provides no separate Task 9 wall-clock budget. The user
sets a maximum of 30 minutes per independently green initial slice.

The primary technical risk is not polynomial root isolation; that substrate is
green. It is unified exact topology ownership across P-P, P-S, and S-S duals,
especially algebraic endpoint equality, degeneracy normalization, and stable
provenance.

Use this performance axiom before any profiler:

```text
For a fixed normalized site graph, work should scale with emitted/incident
Voronoi primitives and isolated roots, not with repeated full graph builds or
all-pairs rematching.
```

First record structural counts. Profile only if a bounded representative
workload violates that axiom after duplicate construction/rematching is
removed.

## 11. Reviews and publication

- Review every focused commit independently against the governing Task 9
  contract and `docs/dev/exact-cgal-idioms.md`.
- Critical/Important findings require a separate fix-forward commit and fresh
  rereview.
- A green native binary is necessary but not sufficient evidence.
- Published reviewed checkpoints must have local HEAD, remote branch SHA, and
  clean/dirty state reported separately.
- Current reviewed remote T9 tip:
  `237ff85cc2fadb8b90775c2485d6bbe0601559b0`.
- Current local committed T9 tip:
  `2cab98c5f7fb5148bf26db39ff7e4135fa9f63de`.
- Current local tip is not reviewed or published.

## 12. Start

Run this literal first command from the active T9 worktree:

```bash
git status --short --branch --untracked-files=all
```

Then read the three-file diff, run the committed-vs-WIP native predicate in
isolation, and decide the cocircular experiment only from exact observed
topology. Do not start another vertical until the tree is clean or a proved
fix-forward commit is green and reviewed.

End of prompt.
