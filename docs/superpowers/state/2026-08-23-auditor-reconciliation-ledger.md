# Auditor reconciliation ledger

**Recorded:** 2026-08-23
**Integration source:** `073a0f7f833da440fe323a0a046ad31550579f82`
**Certifier source:** `73d53729564851d69eb1c28d6c2f1e68350cd20f`
**Merge base:** `1860167929e50f38fdae8d67ef77e9c967364f1c`
**Integration divergence:** 254 commits
**Certifier divergence:** 36 commits

## Source topology

The convergence worktree is
`/Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-auditor-convergence-sdd`
on `codex/auditor-convergence-sdd`. At the topology freeze it was clean at
`233527a35bffccaf364ccc4a741a64241d1c8234`.

The integration source worktree is
`/Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-exact-certified-adaptive-phase1-t9-zero-guide`
on `codex/exact-certified-adaptive-phase1-t9-zero-guide`. Its pre-execution
status was clean at `073a0f7f833da440fe323a0a046ad31550579f82`.

The certifier source worktree is `/Users/jelle/Code/CADCAM/compas_cgal_prs` on
`jf/toolpath-redesign`. Its pre-execution status contained the unrelated
tracked modification `docs/examples/example_isolines.py` and unrelated
untracked path `0` at `73d53729564851d69eb1c28d6c2f1e68350cd20f`.
Those two paths are excluded from every comparison and are never staged,
modified, or cleaned by this programme.

Read-only ref verification from the convergence worktree produced:

```text
073a0f7f833da440fe323a0a046ad31550579f82 refs/heads/codex/exact-certified-adaptive-phase1-t9-zero-guide
73d53729564851d69eb1c28d6c2f1e68350cd20f refs/heads/jf/toolpath-redesign
merge-base 1860167929e50f38fdae8d67ef77e9c967364f1c
left-right 254 36
```

No Git mutation is performed against either source worktree or source ref.

## Commit dispositions

`git cherry -v 073a0f7 73d5372` proves exact patch equivalence for three rows;
all other dispositions use inspected symbols, tests, and current integration
history rather than commit subjects alone.

| Source commit | Subject | Symbols/files | Disposition | Integration evidence | Target plan |
| --- | --- | --- | --- | --- | --- |
| `33cbcb6d589a6d65f75c558f5d1a87f84f7fde51` | zone-query measured A/B docs | `docs/engagement_zone_query.md`, `mkdocs.yml` | unrelated | Historical performance note; no auditor contract or executable evidence | none |
| `23a1de0a2f78fe4e077ff5055b35c50807cccc84` | validated `PocketSpec` | `PocketSpec.build`, `test_spec.py` | equivalent | Patch-equivalent to `a1c0101`; current `tests/benchmarks/test_spec.py` owns invariants | P0 |
| `6f404354f4879778f87bcb3c6310e15adbf4b37e` | analytic pocket families | `disk`, `rectangle`, `stadium`, `arc_channel` | equivalent | Patch-equivalent to `87ea587`; current `test_analytic.py` covers closed forms | P0 |
| `4c058e339a39ca4e0bed22e88bcff20922bdee18` | congruence invariants | rational rotation/translation tests | equivalent | Patch-equivalent to `e4495b6`; current `test_congruence.py` retains controls | P0 |
| `ebf0dcbd080ab53f4ff18517f6c3db507bc6b0f5` | arrangement/rational-digit probes | `Stock.arrangement_stats`, `coordinate_digits` | superseded | `7c1da88` plus `ce04adc` expose both probes and bind them into benchmark measurement; `test_instrument.py` covers both | P4 |
| `83b36ab1f5ef0da77ac6087410fb0a650393d49e` | initial Pixi manifest | `baseline`, `affected`, environment tests | superseded | `d42c688` established locked Pixi; current manifest also owns docs, native gates, schema, testmon, Ruff, and mypy | P4 |
| `8b1d431e8993de7bd4db233e7cbaa60ac8eb4ac7` | Pixi strict-config/editable helper | `scripts/pytest_editable.py` | superseded | Current no-wrapper `_editable-rebuild` and `SKBUILD_EDITABLE_SKIP` tasks provide one explicit rebuild and named task boundaries | P4 |
| `e333c46f8c9dda8e380f4b2b3e20dd555d5b1853` | scoped lint gate | `lint`, owned Python cleanup | superseded | Current `lint = "ruff check src/compas_cgal tests"`; P4 expands enforcement without importing legacy edits | P4 |
| `2af221cd50b88ad443d3823a3dcd75ed1ea73f00` | strict-config comment correction | `pyproject.toml` comment | unrelated | Commentary for the retired `scripts/pytest_editable.py` path | none |
| `d90e18b8c69e7ef4b428000652b13467dabb302a` | bound legacy `mat_scale` | legacy `trochoidal_toolpath`, gouge tests | unrelated | P3 regulates the `engagement_*toolpath` family and its typed override ledger, not the legacy Epick generator | none |
| `2568453b2ea5ca5ee052fd3bf0a2073459761043` | repair editable rebuild | `_editable-rebuild` task | superseded | Current task imports every native module before pytest and was exercised by the P0 environment build | P4 |
| `b0094cb61d1558ecb4f19d683301b4d71074e249` | certify unlinked traverse | `link_paths=False` legacy generator | unrelated | New authoritative audit consumes an explicit ordered operation stream; it never infers omitted legacy traverses | none |
| `446daca5e3c2a7a4009ba3b3105d7168b75c3c28` | lift unlinked traverse | legacy `clearance_z` path | unrelated | Typed approach/link/plunge/retract operations already carry cut/clearance planes in `adaptive/operation.py` | none |
| `b2fa087e6b8c87aab56d1723d465b95eb4805d0d` | discontinuous tessellation | legacy `tessellate_toolpath` | unrelated | Auditor P1 classifies canonical operations directly; no polyline reconstruction is an authority | none |
| `34bddcc0253040948bddd7a87af81950c64d24f0` | output junction dedup | `OUTPUT_DEDUP_TOL` | unrelated | Rendering/tessellation tolerance is outside exact audit identity and verdict semantics | none |
| `2f983d8b84cbaa5a1f8d26b4487487cedc82015a` | exact engagement seam validation | `require_finite`, `require_positive_tool_radius` | dependent | Current legacy `engagement_2` lacks the complete native seam tests; P1 must retain these guards around the ported certifier | P1 Task 4 |
| `832765b053410754b1b129dd0fbb3dfffe2d3306` | exact Stock seam validation | `exact_boundary.h`, all `Stock2` entries | dependent | Current integration only pins selected annulus nonfinite inputs; port needs the complete exact-injection boundary | P1 Task 4 |
| `5c489fbd522a12d2e9d12c3c070880596d2e21ef` | bounded disk-chain count | `chain_intervals`, `MAX_CHAIN_INTERVALS` | dependent | No bounded-chain symbol or unbuildable-chain control exists at `073a0f7`; arc replay depends on loud finite allocation bounds | P1 Task 4 |
| `92cedda24e95ca0b7cb435e351ff93525f34442b` | legacy toolpath parameter contracts | `radial_clearance`, radius/plane checks | unrelated | P3 factories own generator-specific parameter invariants; legacy API patch is not a dependency | none |
| `728b6c8215a6242246dd3d2a6ce4a9954ef3a8e0` | one compiled growth guard | `tea_growth_bound`, `tea_guard` | dependent | Establishes the source certifier lineage that the swept-annulus repair replaces; never imported as certification authority | P1 Task 4 |
| `ebb627d3379e8fd0c92c72abe57b5586187dde11` | growth-bound falsifier | `test_growth_bound.py` harness | dependent | Integration lacks this independent legacy falsifier; its liveness controls are prerequisites to trusting the repair | P1 Task 4 |
| `eee1c35bde6d0c9cd28f4361ebbff6d2c6f0d890` | annular-rib false certificate | `test_false_certificate.py` | required | C2 negative control absent on integration; it proves a station-green segment can contain a cap violation | P1 Task 4 |
| `3c270cb5f3279c95a531a9d99ac543ca0802e783` | spiral-rib generalization | spiral probe/alignment controls | dependent | Prevents the annular-rib repair from overfitting one radius/alignment | P1 Task 4 |
| `c3bdef8419c437a4606f9b3031d16d8b426197c2` | separate witness liveness/verdict | falsifier assertions | dependent | Required evidence discipline: a dead witness may not make a negative control green | P1 Task 4 |
| `2e2422b58d2a6a03fc394520d190531af3b33b2c` | exact swept-annulus segment guard | `interior_run_within_cap` | required | Integration segment oracle is different; adaptive arc support still lacks this proved interior guard | P1 Task 4 |
| `03866fccd2b54b7034b9ab63c984cbe1332ff31c` | rotation-invariant interior bound | exact direction spans | dependent | Removes the axis-aligned-box dependence exposed by rotated falsifiers; needed with the swept guard | P1 Task 4 |
| `e446d035f5e8bcc2f5a93ec73288ac0520c28a9f` | witness statement cleanup | false-certificate docstring | dependent | Preserves witness provenance without asserting an incidental station count | P1 Task 4 |
| `3db92d6445b15e34a55d5de768fede6c9d5a6335` | sub-ulp rim reporting repair | `max_run_tea` span normalization | superseded | `62aea8f` independently forbids promoting a rim sub-arc to a full turn; `test_a_full_turn_report_means_a_buried_rim` is the stronger control | P1 Task 4 |
| `f227dab2fe04af8d8956f623c2839579c0cd7142` | TEA certificate explanation | `docs/engagement_certificate.md` | dependent | The negative-control derivation accompanies the port, but P1 documentation must describe the new authoritative package | P1 Task 5 |
| `387bb7efa36783d8149916cde1f7a2b50afc373f` | correct rim-repair status | certificate docs | superseded | Integration executable full-turn controls already encode the repaired status | P1 Task 5 |
| `5c1d492a9478de5feece1dfcaefbca5d30fcb130` | earlier ETH remediation plan | findings/spec documents | superseded | Approved `233527a` SDD binds the current 254/36 topology and C1-C7 closure | P0 |
| `148a49ad223c8aa4e620890ee70f5afc11defe74` | release shared-root guard | `test_release_build_still_enforces_the_shared_root_precondition` | required | Current `engagement_2.cpp` relies on a CGAL shared-root assumption without this release negative control | P1 Task 4 |
| `48b9e08fab687989d81c253d5b710e64c2537a52` | slack/full-turn preconditions | `<= pi` guard, overflow-safe derivation | dependent | Required by the final arc/segment report assembly and source-side ULP budget | P1 Task 4 |
| `c90a4fef4a4c1da3e98beae75033edfe91fe2619` | machining-quality task track | benchmark plan/policy | superseded | P3/P4 plans separate instrumentation from product gates and bind generator/audit identities | P3/P4 |
| `f369b3df9fd2d2645f02826781463796f93b7434` | reviews and Figure 6 artifacts | review documents, SVG/PDF | unrelated | Review findings informed the approved SDD; generated artifacts are not source-code dependencies | none |
| `73d53729564851d69eb1c28d6c2f1e68350cd20f` | adaptive arc certification | `_certify_arc_engagement`, native recursive arc proof | required | C2 remains unimplemented at `073a0f7`; P1 adds typed partial arcs and ports this dependency-closed proof | P1 Tasks 2/4 |

## Required dependency closure

The source-side proof-bearing transfer order is:

1. `2f983d8`, `832765b`, `5c489fb` — exact input seams and bounded stock
   construction.
2. `728b6c8`, `ebb627d`, `eee1c35`, `3c270cb`, `c3bdef8` — falsified
   predecessor, live negative controls, and verdict separation.
3. `2e2422b`, `03866fc`, `e446d03` — swept-annulus proof and
   rotation/witness closure.
4. `148a49a`, `48b9e08` — release shared-root and full-turn/slack
   preconditions.
5. `73d5372` — adaptive arc recursion, only after the prerequisites above.

`3db92d6` is not transferred because integration commit `62aea8f` has the
stronger rim-reporting contract. `f227dab` is documentation input, rewritten
for the new package after executable proof transfer.

### Exact seam validation

- Source: `require_finite`, `require_positive_tool_radius`, and
  `exact_boundary.h` from `2f983d8`/`832765b`.
- Integration: typed adaptive factories validate many Python inputs, but the
  native legacy engagement/stock entry points do not have the complete source
  boundary.
- Difference: nonfinite or nonphysical binary64 values can reach exact
  injection through uncovered entry points.
- Negative control: source tests for nonfinite segment endpoints, centers,
  boundary vertices, and mixed-radical inputs.
- Positive control: ordinary finite, positive stock and engagement calls.
- Destination: P1 Task 4.

### Bounded disk chain

- Source: `chain_intervals` and `MAX_CHAIN_INTERVALS` from `5c489fb`.
- Integration: no bounded-chain symbol or unbuildable-chain test exists.
- Difference: finite input alone does not bound allocation/work.
- Negative control: `test_subtract_capsule_refuses_an_unbuildable_chain` and
  `test_subtract_arc_sweep_refuses_an_unbuildable_chain`.
- Positive control: `test_ordinary_sweeps_are_far_under_the_chain_limit`.
- Destination: P1 Task 4.

### Compiled growth guard

- Source: `tea_growth_bound`/`tea_guard` from `728b6c8`, followed by its
  falsification and replacement.
- Integration: exact adaptive segment/circle oracles exist, but the legacy
  Python audit still owns a separate sampled arc path.
- Difference: the analytic growth lemma is not sound certification authority;
  it remains useful only as predecessor evidence.
- Negative control: `test_bound_holds_against_a_void_smaller_than_the_tool`
  and the later rib counterexamples.
- Positive control: straight-wall orientation sweep.
- Destination: P1 Task 4; never expose the old guard as an authority.

### False-certificate liveness

- Source: `eee1c35`, `3c270cb`, and `c3bdef8`.
- Integration: `tests/test_false_certificate.py` and
  `tests/test_growth_bound.py` are absent.
- Difference: the current tree cannot demonstrate that its negative witnesses
  are live independently of the verdict under test.
- Negative control: annular rib, machined rib, sector rib, spiral rib, and
  alignment sweeps.
- Positive control: motion clear of the rib and motion whose station lands on
  the violation.
- Destination: P1 Task 4 before implementation transfer.

### Swept-annulus segment certificate

- Source: `interior_run_within_cap` from `2e2422b`.
- Integration: `continuous_tea_2` owns a newer exact segment event oracle, but
  no typed partial-arc audit path consumes a dependency-closed equivalent.
- Difference: station verdicts alone do not prove the continuous interior.
- Negative control: a station-green motion whose interior crosses the annular
  rib.
- Positive control: wall-following motions at every orientation.
- Destination: adapt behind P1's narrow native audit seam; do not replace the
  newer segment oracle.

### Rotation-invariant bound

- Source: exact direction spans from `03866fc`.
- Integration: no equivalent source-side arc interior bound is wired.
- Difference: an axis-aligned bounding argument is orientation-dependent.
- Negative control: rationally rotated rib/alignment cases.
- Positive control: identical verdict across the complete orientation sweep.
- Destination: P1 Task 4 with the swept-annulus guard.

### Rim reporting

- Source: sub-ULP repair `3db92d6`.
- Integration: `62aea8f` independently forbids promoting a rim sub-arc to a
  full turn.
- Difference: none requiring transfer; integration is stricter.
- Negative control: a near-zero contact arc at large coordinate scale.
- Positive control: `test_a_full_turn_report_means_a_buried_rim`.
- Destination: reuse the integration implementation and regression in P1.

### Release shared-root guard and full-turn theorem

- Source: `148a49a` and `48b9e08`.
- Integration: `engagement_2.cpp` documents the CGAL shared-root assumption,
  but lacks the release-build negative control; full-turn reporting has a
  stronger integration regression.
- Difference: a debug-only assertion is not a release proof boundary, and the
  slack derivation needs its explicit `<= pi` domain.
- Negative control:
  `test_release_build_still_enforces_the_shared_root_precondition` plus an
  over-domain span.
- Positive control: legitimate shared-root consecutive arcs and a truly buried
  full rim.
- Destination: P1 Task 4.

### Adaptive arc certification

- Source: native recursive arc proof and `_certify_arc_engagement` wiring from
  `73d5372`.
- Integration: `ExactArcMotion` is absent; the legacy Python arc audit is not a
  native three-verdict authoritative boundary.
- Difference: C2 remains open for partial arcs.
- Negative control: annular/spiral rib arc whose sampled stations appear safe.
- Positive control: constant-stock arc with a known below-cap maximum.
- Destination: P1 Tasks 2 and 4 after the full dependency closure above.

## Baseline

All commands ran against production tree `233527a35bffccaf364ccc4a741a64241d1c8234`.
The only working-tree changes were this ledger and the SDD progress document.

| Gate | Exit | Result | Wall time | Classification |
| --- | ---: | --- | ---: | --- |
| requested false-certificate/growth/audit group | 5 | no tests ran; first two files absent | 101.66 s | missing executable evidence |
| `tests/test_engagement_audit.py` | 0 | 15 passed | 90.05 s | green legacy audit baseline |
| adaptive replay/segment/circle group | 0 | 53 passed | 101.22 s | green focused baseline |
| `types-adaptive` | 1 | 6 errors in 3 files | 4.64 s | known typing regression baseline |
| `lint` | 0 | all checks passed | 0.39 s | green |
| strict docs | 0 | built successfully | 5.82 s | green |
| all benchmark tests | 1 | 259 passed, 6 product failures | 145.24 s | instrument green; product red |
| isolated quality test | 1 | 36 passed, same 6 product failures | 87.03 s | deliberately red product gate |
| full `baseline` | 130 | interrupted after 20:41; summary at 19:25 was 1,315 passed, 10 failed, 38 warnings, one CPU-bound test outstanding | 1,241.37 s | known failures plus unresolved extreme tail |

The six strict-mypy errors are:

- unused ignores in `adaptive/replay.py` and `adaptive/bootstrap.py`;
- two `int | None` arguments passed where route indices require `int`;
- `ExactCircleMotion.end` accessed through a segment/circle union;
- `source_commit_digest` accessed through a traversal/retrace union.

The four adaptive full-suite failures are:

- `tests/adaptive/test_route_retrace_generator.py::test_route_retrace_derivation_rejects_unsupported_source_scope[terminal]`;
- `tests/adaptive/test_route_retrace_generator.py::test_continuation_rejects_missing_retrace_commit`;
- `tests/adaptive/test_generator.py::test_real_active_family_stops_at_unresolved_exact_event`;
- `tests/adaptive/test_generator.py::test_task13f_full_continuation`.

The six product failures are every `engagement_controlled` and
`radius_regulated` combination for `rect_12x8`, `rect_20x12`, and `L_shape` in
`test_the_generated_path_is_worth_running`. The reports consistently identify
uncut stock, degenerate loops, redundant operations, cap exceedances, large
engagement steps, and tangent breaks; `rect_20x12` additionally reports
slotting. They are product-gate failures, not instrument failures.

The full run did not produce a new failing node ID. It did prove an unresolved
performance outcome: one exact-geometry test remained CPU-bound after the
suite had reached its final tail. P4 must identify and budget or decompose that
test before CI enforcement; this result is not called green.

## P0 acceptance

- All 36 certifier-side commits have one allowed disposition; none is
  unadjudicated.
- Required/dependent transfers have an explicit source order, negative
  control, positive control, and destination task.
- P0 changed documentation only.
- `git diff --check` passed.
- `pixi run -e docs docs` passed under strict MkDocs.
- Focused audit and adaptive oracle/replay tests passed.
- Known mypy, adaptive-continuation, quality-gate, missing-falsifier, and
  full-suite-tail outcomes are preserved as negative evidence, not hidden.
- Source refs remained fixed at `073a0f7` and `73d5372`.

P0 is accepted as a truthful reconciliation baseline. It does not claim the
integration tree is release-green.
