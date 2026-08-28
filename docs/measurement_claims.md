# Measurement-claim ledger

> **status: in audit — 0/14 Task-5 extractor rows dispositioned**

- Opened (UTC): `2026-08-28`
- Programme: coherence Wave 1 backlog A

## Scope and limits

This ledger freezes the Task-5 Python-comment regex population. It is not a
claim of repository-wide invariant-I2 closure. Tasks 6 and 7 may change only
the `disposition` and `evidence / change` cells; all identity columns and
frozen matched lines remain immutable.

## Frozen extraction

Extraction source commit: `eec665c1df1cd8d1e98dd9dd1001b5984e17a703`

Population: 14 line hits across 5 files

The historical extractor was:

```bash
grep -rnE "# .*(MEASURED|[Mm]easured (on|at|against)|measures [0-9]|[0-9]+(\.[0-9]+)?[x×] (faster|slower))" \
  src/compas_cgal benchmarks --include='*.py' | grep -v superseded
```

Canonical order is `src/compas_cgal` before `benchmarks`, then `LC_ALL=C`
relative-path order, numeric source line, and raw hit as the final tie-breaker.

Execution scan (2026-08-28 UTC): exact match — 14/14 frozen hits present;
0 changed, 0 missing, 0 new. Normalized stream SHA-256:
`b78ac690bceda37fdebcdcdeb9731a8b213e7d946a897be45181cb95b1100f66`.
This scan is drift evidence only and does not redefine the frozen population.

## Disposition semantics

- `pending` — final adjudication has not happened.
- `re-earned` — an authenticated rerun confirmed every material assertion and records the exact command, configuration, artifact, result digest, and full input commit.
- `corrected` — the assertion was wrong or incomplete; authenticated evidence and the source correction commit are recorded.
- `historical` — the original configuration cannot be reconstructed; the source is explicitly labelled and names the missing identity or configuration.
- `deleted` — the assertion was removed while its frozen identity remains here.
- `not-a-claim` — semantic review found a regex false positive and records an explicit rationale.

A case-level result does not automatically disposition every mapped row. Partial
reproduction cannot become `re-earned`; every material assertion is adjudicated
row by row. `reproduced` is not a ledger disposition.

## Ledger

| ordinal | claim ID | extracted location | stable anchor | anchor match | disposition | evidence / change |
| --- | --- | --- | --- | --- | --- | --- |
| 001 | MC-001 | `src/compas_cgal/engagement_radial_toolpath.py:193` | `compas_cgal.engagement_radial_toolpath.RADIUS_LADDER_SUBDIVISIONS` / leading comment | 1/6 | pending | — |
| 002 | MC-002 | `src/compas_cgal/engagement_radial_toolpath.py:194` | `compas_cgal.engagement_radial_toolpath.RADIUS_LADDER_SUBDIVISIONS` / leading comment | 2/6 | pending | — |
| 003 | MC-003 | `src/compas_cgal/engagement_radial_toolpath.py:195` | `compas_cgal.engagement_radial_toolpath.RADIUS_LADDER_SUBDIVISIONS` / leading comment | 3/6 | pending | — |
| 004 | MC-004 | `src/compas_cgal/engagement_radial_toolpath.py:201` | `compas_cgal.engagement_radial_toolpath.RADIUS_LADDER_SUBDIVISIONS` / leading comment | 4/6 | pending | — |
| 005 | MC-005 | `src/compas_cgal/engagement_radial_toolpath.py:207` | `compas_cgal.engagement_radial_toolpath.RADIUS_LADDER_SUBDIVISIONS` / leading comment | 5/6 | pending | — |
| 006 | MC-006 | `src/compas_cgal/engagement_radial_toolpath.py:225` | `compas_cgal.engagement_radial_toolpath.RADIUS_LADDER_SUBDIVISIONS` / leading comment | 6/6 | pending | — |
| 007 | MC-007 | `src/compas_cgal/engagement_radial_toolpath.py:257` | `compas_cgal.engagement_radial_toolpath.RADIUS_LADDER_FLOOR_STEPS` / leading comment | 1/1 | pending | — |
| 008 | MC-008 | `src/compas_cgal/engagement_radial_toolpath.py:284` | `compas_cgal.engagement_radial_toolpath.RADIUS_LADDER_REFINEMENT_MARGIN` / leading comment | 1/1 | pending | — |
| 009 | MC-009 | `src/compas_cgal/engagement_toolpath.py:120` | `compas_cgal.engagement_toolpath.LOOP_PROBE_COUNT` / leading comment | 1/2 | pending | — |
| 010 | MC-010 | `src/compas_cgal/engagement_toolpath.py:141` | `compas_cgal.engagement_toolpath.LOOP_PROBE_COUNT` / leading comment | 2/2 | pending | — |
| 011 | MC-011 | `benchmarks/gate.py:58` | `benchmarks.gate.GATE_LARGE_RECT_WIDTH` / leading comment | 1/1 | pending | — |
| 012 | MC-012 | `benchmarks/gate.py:66` | `benchmarks.gate.L_ARM_TOOL_DIAMETERS` / leading comment | 1/1 | pending | — |
| 013 | MC-013 | `benchmarks/mathsm.py:47` | `benchmarks.mathsm.SPACING_SWEEP_TOOL_DIAMETERS` / leading comment | 1/1 | pending | — |
| 014 | MC-014 | `benchmarks/quality.py:150` | `benchmarks.quality.IMMERSION_STEADY_BAND_FRACTION` / leading comment | 1/1 | pending | — |

## Frozen matched lines

### MC-001

```text
# Measured on the 20x12 pocket at a 60 deg cap, station (18.482, 10.482), maximal
```

### MC-002

```text
# radius 0.5156, coarse step 0.05: rung 6 (radius 0.2156) measures 61.3 deg and
```

### MC-003

```text
# still cuts, rung 7 (radius 0.1656) measures 5.9 deg and cuts NOTHING -- one step
```

### MC-004

```text
# WHY THIS COUNT: MEASURED, NOT DERIVED -- the same footing as `LOOP_PROBE_COUNT`,
```

### MC-005

```text
# Measured on the 20x12 pocket at a 60 deg cap, tool diameter 2.0, over the 244
```

### MC-006

```text
# AND DO NOT KEEP THIS WHILE DROPPING `_least_bad_rung`. Measured on 6x4 at a
```

### MC-007

```text
# MEASURED, on the 20x12 pocket at a 60 deg cap with the gate below in place:
```

### MC-008

```text
# MEASURED, on the 6x4 pocket at a 40 deg cap -- the hard case, a pocket three tool
```

### MC-009

```text
# is FALSE. Measured on a 20x12 pocket, 2 mm tool, by walking every machining
```

### MC-010

```text
# WHY THIS COUNT: MEASURED CONVERGENCE, NOT A DERIVATION. There is a geometric
```

### MC-011

```text
# The pocket the corner defect was measured on. Ten by six tool diameters.
```

### MC-012

```text
# W > 4r, hence an arm STRICTLY WIDER THAN TWO TOOL DIAMETERS. Measured at
```

### MC-013

```text
# spacing stops helping -- on the reference pocket 0.025 measures 131.14 degrees
```

### MC-014

```text
# measured at the same cap, which is also what makes two generators comparable.
```
