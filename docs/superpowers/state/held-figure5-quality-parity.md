# Held Figure 5 quality parity

Date: 2026-09-02

## Decision

The additive attributed reducer is **equivalent** to the current quality judge
for every gated field on the synthetic machinery, invariant-generated paths,
and all twelve `QUALITY_GATE_CASES`. The repository's declared seventeen-red
acceptance condition is now deterministic and satisfied: the full manifest run
found the twelve quality reds, four exact adaptive reds, and the pinned benchmark
translation-invariance defect, with no unexpected failures.

Production still calls `measure_quality()` and its existing private reducers.
No reducer was changed or removed, and `measure_quality()` was not redirected.
Removing the old reducers remains an explicit decision for Jelle; Task 6 must
not start without that approval.

## Compared contract

The parity helper compares each existing gate input with its typed observation
and pins the observation threshold to the existing gate constant:

| Existing field | Typed observation | Required |
| --- | --- | ---: |
| `elementary.uncut_fraction` | `uncut_fraction.measured` | 0 |
| `elementary.gouging_motions` | `gouging_motions.measured` | 0 |
| `elementary.unsafe_rapids` | `unsafe_rapids.measured` | 0 |
| `elementary.continuity_breaks` | `continuity_breaks.measured` | 0 |
| `elementary.zero_length_motions` | `zero_length_motions.measured` | 0 |
| `elementary.degenerate_loops` | `degenerate_loops.measured` | 0 |
| `elementary.redundant_operations` | `redundant_operations.measured` | 0 |
| `cut.cap_exceedances` | `cap_exceedances.measured` | 0 |
| `cut.slotting_motions` | `slotting_motions.measured` | 0 |
| `cut.max_engagement_step_deg` | `max_engagement_step.measured` | case cap |
| `cut.max_loop_radius_step` | `max_loop_radius_step.measured` | 2 tool radii |
| `speed.tangent_breaks` | `tangent_breaks.measured` | 0 |

For each of the nine count observations, the test independently checks that the
measured count equals the complete attributed operation or pair tuple length.
Uncut stock remains spatial aggregate evidence and has no invented operation.
For every synthetic, invariant-generated, and gate invocation, a test-side
oracle independently reconstructs adjacent engagement candidates and loop
runs from the raw survey and operation snapshot. It resets loop succession at
rapid and path-chain boundaries, retains the first occurrence of a tied
maximum, and requires the attributed maximum to name that exact pair; a zero
maximum must carry `None`. It does not call the attributed reducer. The
synthetic distinct-maximum case additionally pins the engagement winner to
operation pair `2 -> 3` and the loop-radius winner to operation pair `1 -> 3`.
The completed full manifest run also executes the Task 4 first-tie contracts, which pin
engagement pair `0 -> 1`, loop pair `0 -> 2`, every over-threshold pair, and
the no-pair result for a zero engagement maximum.

## Commands and results

### Task 5A.4 final transition evidence (2026-09-04)

- Green transition oracle: `12 passed in 167.02s`; every unchanged case used the
  shared pre-verdict evaluator and matched its reviewed ordered violation vector.
- Protected red wrapper: exactly 12 failures in 162.54s, all at the final
  `assert not violations`; exact JUnit IDs and reviewed messages matched.
- Exact-source mutation coverage: `12 passed in 2.81s`; the complementary
  same-cardinality wrong-source rejection gate passed `47 passed in 2.34s`.
- Full repository reconciliation: `2735 passed, 17 failed`; `red set == manifest,
  both directions`.
- Native stability: three independent runs each completed with `21 passed, 4
  failed`; each JUnit report matched the same four exact adaptive IDs and no
  process error occurred.
- Production routing remains unchanged at `measure_quality()`; Task 6 has not
  started.

```text
pixi run pytest -- tests/benchmarks/test_quality.py tests/benchmarks/test_quality_invariants.py -k 'not test_the_generated_path_is_worth_running and not test_moving_the_pocket_across_the_table_changes_no_metric' -n auto --testmon -q
```

Fix-round-1 result: exit 0, `43 passed in 11.75s`.

```text
pixi run pytest -- tests/benchmarks/test_quality.py::test_the_generated_path_is_worth_running -n auto -q
```

Fix-round-1 result: required exit 1, `12 failed, 12 warnings in 170.58s`.
Every cell reached the unchanged final assertion after field, exact winning-pair
parity, and the local spies had proved exactly one `PathSurvey` and one
`CoverageEstimate`. There were no parity, spy, snapshot, collection,
configuration, interruption, or infrastructure failures.

The exact quality-red membership and unchanged criterion messages were:

| Test ID | Violated criteria (`measured <= required`) |
| --- | --- |
| `cap-120-default-engagement_controlled-rect_12x8` | uncut `0.008041 <= 0`; degenerate `4 <= 0`; redundant `10 <= 0`; cap `5 <= 0`; engagement step `329.49 <= 120`; tangent `24 <= 0` |
| `cap-120-default-radius_regulated-rect_12x8` | uncut `0.008041 <= 0`; degenerate `4 <= 0`; redundant `10 <= 0`; cap `5 <= 0`; engagement step `329.49 <= 120`; tangent `24 <= 0` |
| `cap-120-default-engagement_controlled-rect_20x12` | uncut `0.002843 <= 0`; degenerate `4 <= 0`; redundant `10 <= 0`; cap `9 <= 0`; slotting `4 <= 0`; engagement step `337.54 <= 120`; tangent `32 <= 0` |
| `cap-120-default-radius_regulated-rect_20x12` | uncut `0.002843 <= 0`; degenerate `4 <= 0`; redundant `10 <= 0`; cap `9 <= 0`; slotting `4 <= 0`; engagement step `337.54 <= 120`; tangent `32 <= 0` |
| `cap-120-default-engagement_controlled-L_shape` | uncut `0.006930 <= 0`; degenerate `5 <= 0`; redundant `10 <= 0`; cap `6 <= 0`; engagement step `325.42 <= 120`; tangent `22 <= 0` |
| `cap-120-default-radius_regulated-L_shape` | uncut `0.006930 <= 0`; degenerate `5 <= 0`; redundant `10 <= 0`; cap `6 <= 0`; engagement step `325.42 <= 120`; tangent `22 <= 0` |
| `cap-40-attribution-engagement_controlled-rect_12x8` | uncut `0.003186 <= 0`; degenerate `80 <= 0`; redundant `26 <= 0`; cap `145 <= 0`; slotting `40 <= 0`; engagement step `356.18 <= 40`; tangent `184 <= 0` |
| `cap-40-attribution-radius_regulated-rect_12x8` | uncut `0.005462 <= 0`; degenerate `124 <= 0`; redundant `311 <= 0`; cap `54 <= 0`; slotting `41 <= 0`; engagement step `349.89 <= 40`; tangent `798 <= 0` |
| `cap-40-attribution-engagement_controlled-rect_20x12` | uncut `0.000335 <= 0`; degenerate `112 <= 0`; redundant `62 <= 0`; cap `185 <= 0`; slotting `108 <= 0`; engagement step `357.19 <= 40`; tangent `400 <= 0` |
| `cap-40-attribution-radius_regulated-rect_20x12` | uncut `0.000669 <= 0`; degenerate `308 <= 0`; redundant `490 <= 0`; cap `125 <= 0`; slotting `116 <= 0`; engagement step `357.19 <= 40`; tangent `1384 <= 0` |
| `cap-40-attribution-engagement_controlled-L_shape` | uncut `0.006059 <= 0`; degenerate `40 <= 0`; redundant `15 <= 0`; cap `221 <= 0`; slotting `15 <= 0`; engagement step `355.04 <= 40`; tangent `92 <= 0` |
| `cap-40-attribution-radius_regulated-L_shape` | uncut `0.002916 <= 0`; degenerate `236 <= 0`; redundant `499 <= 0`; cap `229 <= 0`; slotting `62 <= 0`; engagement step `347.52 <= 40`; tangent `956 <= 0` |

```text
pixi run red-manifest
```

Fresh result: exit 1 after a complete suite run, `16 failed, 2613 passed, 44
warnings in 765.17s`. Actual reds were all twelve quality IDs above plus:

- `tests/adaptive/test_generator.py::test_task13f_full_continuation`
- `tests/adaptive/test_generator.py::test_real_active_family_stops_at_unresolved_exact_event`
- `tests/adaptive/test_route_retrace_generator.py::test_route_retrace_derivation_rejects_unsupported_source_scope[terminal]`
- `tests/adaptive/test_route_retrace_generator.py::test_continuation_rejects_missing_retrace_commit`

The manifest checker reported:

```text
expected-red-went-green: tests\.benchmarks\.test_quality_invariants::test_moving_the_pocket_across_the_table_changes_no_metric
```

That sixteen-red result was a probabilistic miss, not evidence that the defect
was absent or that the new reducer changed a result. The translation-invariance
cell is deliberately excluded from the focused parity GREEN command and was
already unexpectedly green in the pre-Task-1 baseline (`16 failed, 2481
passed`).

Task 5A.1 retained the generated property and added one explicit dyadic witness.
Before the witness, fixed Hypothesis seeds 1 and 2 each exited 0 with `1 passed`.
After the witness, both commands exited 1 at the protected count assertion:

```text
AssertionError: cut.cap_exceedances moved with the pocket
assert 1.0 == 2.0
```

The unchanged manifest-owned node ID therefore now has deterministic failure
evidence independent of generated examples. This establishes its declared-red
membership witness; it does not replace the fresh exact seventeen-red full-run
gate in Task 5A.4.

## Blocker and resumption history

The first live cell revealed that generated vertical `PLUNGE`/`RETRACT`
operations use zero tangent sentinels. Task 2 correction `81417c5` maps only
those sentinels to absent typed tangents and retains strict lateral validation.
The next run revealed that the radius-regulated generator returns the declared
`RadialToolpathResult` subtype. Task 2 correction `0ea3fd5` accepts
`ToolpathResult` subtypes while retaining exact operation validation. Both
corrections were independently reviewed before Task 5 resumed.

The intervening 12-cell run was invalid: six cells reached the expected quality
assertion and six stopped at subtype rejection. Its following full-suite attempt
segfaulted in adaptive native execution. Neither interrupted run is counted as
acceptance evidence. The final runs above started from scratch and completed.
