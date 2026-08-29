# Wave-1 Task-6 Historical Assertions Design

> **status: approved** — design approved 2026-08-29. This specification repairs
> the traceability gap identified at the `7506a47` interim review without
> extending the accepted measurement-result schema.

## Purpose

Commit an immutable artifact containing every source assertion removed by
`53135e0`, including assertions later shown false, incomplete, or
non-reconstructible. The artifact records provenance only. It does not certify
the assertions and cannot satisfy a measurement-claim ledger row by itself.

Task 6 remains open until a separate canonical v2 measurement artifact is
committed and accepted by the ledger validator.

## Decisions

### Preserve the exact correction diff

The historical artifact is the exact byte output of:

```bash
git --no-replace-objects diff --binary \
  b531d215e6ca06741d040b070ba43164f61abd58 \
  53135e04390e84bf69aa74dc4d0c1ce6ca308eb4 -- \
  src/compas_cgal/engagement_radial_toolpath.py \
  src/compas_cgal/engagement_toolpath.py
```

The deleted lines preserve every removed number in its original prose and table
context. Keeping the diff avoids hand-transcribing historical values into a
second semantic schema.

### Bind the diff with a minimal manifest

The sibling `manifest.json` contains exactly:

- schema version `measurement-claim-history/v1`;
- the raw parent and correction Git object IDs;
- the two ordered source paths;
- the patch filename and SHA-256;
- one mapping from each `MC-001` through `MC-010` to its source path;
- status `superseded-source-assertions`, which cannot be confused with an
  accepted measurement disposition.

The artifact directory is:

```text
benchmarks/measurement_claim_history/2026-08-29-53135e04390e/
```

### Verify, do not generalize

A focused contract test recomputes the correction parent from the raw commit
object, recomputes the exact diff, verifies the patch bytes and SHA-256, rejects
extra or missing manifest keys, and proves MC-001 through MC-010 map exactly
once. No production CLI, Pixi task, reusable producer, alternate ledger
renderer, or new measurement-result schema is added.

The test verifies the committed artifact as a repository consumer boundary. It
does not introduce a second acceptance path for current measurements.

### Resolve the `8` versus `12` conflict by provenance

The Wave-1 Task-9 plan records the historical `20x12` baseline as
`126.1/8/3236`; the source assertion removed by `53135e0` records
`126.1/12/3236`. Their configurations are not authenticated well enough to
prove equality. Both values are retained with their distinct origins and
neither is promoted as a reconstructed fact. Task 9 must disclose the conflict
instead of silently choosing one.

## Wave-1 scope ruling

The existing Task-6 framework is retained through Task 11 because it is the
only implemented path that rejects semantic misadjudication, raw-history
substitution, protected-source drift, and manually mismatched ledger cells.
Its scope freezes at the ten generator claims. Task 7 follows its approved
benchmark contract and must not generalize the generator schema.

After Task 11, a separate explicitly approved cleanup may remove producer-only
modules and tests while retaining the committed artifacts, canonical validator,
source-lineage checks, and ledger acceptance boundary. No removal occurs during
Wave 1.

## Fleet-file ruling

Commit `681f045` remains on the Wave-1 branch because reverting it restores a
twelfth unmanifested red. It changes only the diagnostic text introduced by
fleet commit `e678fe0`, restoring the pre-existing `representable chord ratio`
contract. The deviation and absorption rule belong in
`docs/superpowers/state/wave1-gate-analysis.md`: drop `681f045` during S1
absorption only if the fleet has landed a patch-equivalent fix; otherwise
surface it for explicit fleet acceptance.

## Completion conditions

Task 6 closes only when all of the following hold:

1. The historical artifact and its contract test are committed and green.
2. The rejected v1 bundle is preserved outside the accepted-result root and
   labelled rejected.
3. Exactly one canonical v2 result exists under
   `benchmarks/measurement_claim_results/`.
4. MC-001 through MC-010 are rendered from that v2 result; MC-011 through
   MC-014 remain pending.
5. `test_real_ledger_has_authenticated_task6_acceptance` passes.
6. The focused generator tests, strict mypy, Ruff, red manifest, plan-header
   lint, and source-lineage gates pass.
7. The Wave-1 plan header records Task 6 only in the commit that lands this
   accepted state.
