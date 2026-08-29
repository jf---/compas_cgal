# Wave-1 gate analysis

> **status: deviation rulings accepted** — recorded 2026-08-29 for Task 6.

## Fleet diagnostic prerequisite: keep `681f045`

Before `681f0453848925910d081c6cb202fe790deb40b6`, the Wave-1 affected
suite command was:

```bash
pixi run affected
```

It exited 1 with `12 failed, 1502 passed, 38 warnings in 252.86s`. The
twelfth red, outside the eleven-entry manifest, was
`tests.adaptive.test_motion::test_segment_certifier_reuses_native_cap_conversion`.
`pixi run red-manifest` reproduced the same 12 reds, and the direct diagnostic
command:

```bash
pixi run python -m tools.red_manifest build/junit-baseline.xml
```

reported exactly:

```text
unexpected-red: tests.adaptive.test_motion::test_segment_certifier_reuses_native_cap_conversion
```

Fleet commit `e678fe04a7819363d7d6cd315cec532dd021b177` centralized cap
conversion in `src/audit_policy_2.cpp` but changed the stable diagnostic from
`representable chord ratio` to `representable chord surrogate`. Commit
`681f045` changes only that reporting text to `representable chord ratio
surrogate`; exception type, predicate, values, and routing remain unchanged.
The exact focused consumer-boundary command is now green:

```bash
pixi run pytest -- tests/adaptive/test_motion.py::test_segment_certifier_reuses_native_cap_conversion --testmon -n auto -q
```

Fresh Task-2 result: `1 passed in 1.82s`. The related no-deselection run passed 72
tests, and the serialized `pixi run audit-native` gate passed.

Ruling: keep `681f045`. Reverting it would deliberately restore a red outside
`docs/red_manifest.json`, violating I3's definition that every unmanifested
red is a defect. This is a compatibility repair for a fleet-introduced
diagnostic regression, not suppression or reclassification of a red.

At S1 absorption, drop `681f045` only if the fleet has already landed a
patch-equivalent repair: the authoritative audit-policy diagnostic retains the
contiguous `representable chord ratio` contract, the exception semantics and
single call path are unchanged, and the focused command above passes without
this commit. Otherwise surface `681f045` to the fleet for explicit acceptance;
do not silently discard it during conflict resolution.

## Task-6 authenticated-claims framework: retain, freeze, retire later

The prose ledger alone could not prevent a reviewer or later edit from pairing
a plausible disposition with the wrong payload, substituting raw historical
provenance for current evidence, accepting protected-source drift, or manually
mismatching a ledger cell. The Task-6 framework makes those cases executable
consumer-boundary failures through authenticated input/result identities,
semantic adjudication, source-lineage checks, canonical validation, and
renderer-owned ledger evidence.

That guarantee has a measured cost. The following inventory command over the
15 producer/validator modules and five directly owned test modules reports
7,757 lines across 20 files:

```bash
wc -l tools/measurement_claim*.py \
  tests/tools/test_measurement_claim_{advance,ledger,probes,radial,result}.py
```

The accepted cost includes maintenance and review load plus fleet-rebase
surface in shared configuration such as `pyproject.toml`.

Ruling: retain the framework through Task 11, but freeze the generator scope
at exactly MC-001 through MC-010. Task 7 continues under its approved benchmark
contract and must not generalize the generator schema, producer family, or
ledger renderer.

After Task 11, a separate explicitly approved cleanup may retire producer-only
modules and their tests. It must retain the committed artifacts, canonical
validator, source-lineage checks, and ledger acceptance boundary. No framework
removal occurs during Wave 1.

## Historical `8` versus `12` conflict

The Task-9 plan records the historical `20x12` baseline as
`126.1/8/3236`; the source assertion removed by `53135e0` records
`126.1/12/3236`. Their configurations are insufficiently identified to prove
that they describe the same run. They are distinct historical assertions with
distinct provenance. Preserve both and disclose the conflict; neither value is
promoted as reconstructed truth.
