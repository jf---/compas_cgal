# Spec — ETH-grade audit findings, `jf/toolpath-redesign`

**Date:** 2026-08-19
**Audited head:** `33cbcb6` (+ uncommitted `docs/examples/example_isolines.py`)
**Scope:** `src/{toolpath,stock_2,engagement_2}.{cpp,h}`, `src/compas_cgal/{toolpath,stock,engagement,isolines}.py`,
`tests/`, `mkdocs.yml`, `pyproject.toml`, `.github/workflows/`. Branch
`codex/exact-certified-adaptive-phase1-t9` is **out of scope**.
**Baseline:** 119/119 tests pass; `ruff check` clean; `ruff format --check` fails on one file.

This document is the **spec** for `docs/superpowers/plans/2026-08-19-eth-audit-remediation.md`.
Every finding below carries reproduced evidence. Numbers were measured on this branch with the
built extensions (`PYTHONPATH=src`), not inferred from reading.

---

## Verdict

**CONDITIONAL.** The exact-kernel core is publication-grade. Three defects would each fail peer
review, and two of them attack the word the project rests on — *certified*.

The exact algebra is **correct, including the part that is easy to get wrong**:
`GpsTraits::Point_2::CoordNT` is `Sqrt_extension<NT, NT, Tag_true, …>` (verified in
`external/cgal/include/CGAL/Arr_geometry_traits/Circle_segment_2.h:46`), so `ACDE_TAG::value == true`
and `Sqrt_extension::compare`'s default `in_same_extension = !ACDE_TAG::value` is **false** — the
cross-root comparison takes the exact branch
(`external/cgal/include/CGAL/Sqrt_extension/Sqrt_extension_type.h:531`). `sign_mixed_radical`
decomposes into `sign(u)`, `sign(w)`, `compare(u², β·w²)` and never forms cross-root arithmetic.
Run assembly merges on exact point equality, not angular tolerance. The deciding/reporting split is
held consistently. This is real work and the remediation must not damage it.

---

## Critical

### C1 — `mat_scale > 1` silently emits gouging toolpaths under an unconditional guarantee

`src/compas_cgal/toolpath.py:282` and `:379` expose `mat_scale: float = 1.0` as public, documented
API ("Scale factor applied to the clearance-derived available radius") with **no stated range and
no validation**. `validate_toolpath_params` (`src/toolpath.cpp:602`) validates seven parameters and
not this one.

`RadiusModel::radius_from_clearance` (`src/toolpath.cpp:205`) computes
`available = mat_scale · max(0, clearance − R − rc)`. A trochoid circle of radius `r` reaches
`r + R` from its centre, so gouge-freedom requires `mat_scale ≤ 1`. Only bridges (`certify_bridges`)
and leads (`certified_lead`) are certified — **the trochoid circles, the dominant primitive, are
never certified**; they are safe purely by this construction.

Measured on `_dumbbell(2.4)` with tool ⌀1.0:

| `mat_scale` | max cut-circle penetration past the wall | max tool-**centre** excursion outside the pocket |
| --- | --- | --- |
| 1.0 | 0.0000 | −0.0010 (safe) |
| 1.5 | **+2.2485** | **+0.3999** |
| 2.0 | +4.4980 | +0.2596 |
| 3.0 | +8.9970 | +0.1876 |

`tests/test_toolpath.py::test_no_gouge_circular_primitives_dumbbell` asserts exactly the right
invariant — only at the default.

**Required contract:** `mat_scale` is validated to `(0, 1]` at the seam, with a named error, and the
existing gouge tests are parametrised over it.

### C2 — `tea_growth_bound` is not an upper bound; nothing in the suite would notice

`src/engagement_2.cpp:462`. Clause (a) claims a rim/boundary crossing drifts by at most
`2·asin(d/2r)` under centre travel `d`. Near tangency the true drift is `√(2d/r)`:

| `d` (r = 1) | true drift of one crossing | clause (a) bound | ratio |
| --- | --- | --- | --- |
| 1e-2 | 0.141539 | 0.010000 | 14.2× |
| 1e-4 | 0.014142 | 0.000100 | 141.4× |
| 1e-6 | 0.001414 | 0.000001 | 1414.2× |

The composite `GROWTH(d) = 4·asin(d/2r) + 2·acos(1 − d/r)` survives on a **straight** wall, but its
slack collapses to 1.007× at `d = 1e-4`. It does **not** survive on **concave circular** boundary
arcs. For a void disk of radius ρ > r, the exact emerged run half-angle at internal tangency
`s₀ = ρ − r` after travel `δ` is

```
psi**2 = (2*delta*rho + delta**2) / (r * (rho - r + delta))          # exact, not asymptotic
run_growth = 2*psi  ->  2*sqrt(2*delta*rho / (r*(rho - r)))          # small-delta limit
```

so `true / guard = ½·√(ρ/(ρ − r))`, unbounded as ρ → r⁺. Measured at `d = 1e-4`, `r = 0.5`,
verified two independent ways (closed form vs. 4 000 000-sample discrete rim occupancy, agreeing to
4 significant figures):

| void radius ρ (as `r_guide + R`) | true growth | guard `2·GROWTH(d)` | ratio |
| --- | --- | --- | --- |
| ρ − r = 1.00 R | 0.039999 | 0.056969 | 0.70× (holds) |
| ρ − r = 0.30 R | 0.058870 | 0.056969 | **1.03×** |
| ρ − r = 0.10 R | 0.093770 | 0.056969 | **1.65×** |
| ρ − r = 0.01 R | 0.283086 | 0.056969 | **4.97×** |

For a **convex** feature (island of radius ρ) the same derivation gives
`2·√(2δρ / (r(ρ + r)))`, which is *below* the current clause (b) and reduces to it at ρ → ∞. So the
existing bound is correct for straight walls and convex material, and wrong for concave voids.

**Reachability — the part that sets the priority.** Every `Stock2::subtract_*` removes a union of
disks of the *subtraction* radius, so within `audit_toolpath_engagement` the void arcs always have
ρ = r exactly (the query radius equals the radius that cut the stock). A 12 000-sample random probe
against a stock depleted by 81 real toolpath operations found **no violation**:

| `d` | worst observed growth | `GROWTH(d)` | verdict |
| --- | --- | --- | --- |
| 1e-2 | 0.354246 | 0.440670 | holds (1.24×) |
| 1e-3 | 0.108296 | 0.130512 | holds (1.21×) |
| 1e-4 | 0.019760 | 0.040401 | holds (2.04×) |

The violating regime is reachable through the **public `Stock` API** (`subtract_disk` /
`subtract_capsule` / `subtract_arc_sweep` all take a free `radius`) and through any audit whose tool
radius differs from the radius that cut the stock — not through the current `audit_toolpath_engagement`
wiring. **No end-to-end false `cap_certified = True` was constructed**; a 15-configuration sweep
through the real `certify_segment_tea` produced none, and adaptive refinement absorbed every case.

At ρ = r exactly there is additionally a genuine **discontinuity**: TEA reads 0° at δ = 0 (measure-
zero grazing, `found_overlap` discards it) and 180.00° at δ = 1e-9, matching
`2π − 2·arccos(δ/2R)` to six figures. This is an isolated *downward* spike (a station under-reads at
one point while its neighbours over-read), so refinement tends to catch it, but no modulus-of-
continuity bound exists there and the code does not say so.

**Conclusion.** The certificate is not *demonstrated* unsound on the shipped path. It holds by an
accident of the current wiring, its written proof does not support it, and no test defends the
accident. `TEA_GUARD_SAFETY_FACTOR = 2` is documented as burying "any looseness in the lemma and the
~1e-15 relative error of asin/acos" (`src/engagement_2.cpp:476`); it is not slack, it *is* the bound.

**Required contract:** a property test that falsifies the bound against exact ground truth across a
parameterised family; a curvature-aware bound covering the concave case; documentation of the
regime where no bound exists.

### C3 — `requires-python = ">=3.9"` but the new modules need ≥ 3.10

`toolpath.py`, `stock.py`, `engagement.py`, `isolines.py` all use PEP 604 `X | Y` in **evaluated**
annotations with no `from __future__ import annotations` (22 `BitOr` nodes in `toolpath.py` alone).
On 3.9 the module fails at import with `TypeError: unsupported operand type(s) for |`. CI builds
3.10 only (`.github/workflows/build.yml:18`), so nothing catches it. `pyproject.toml` declares
`requires-python = ">=3.9"`.

---

## Important

| # | Finding | Anchor | Evidence |
| --- | --- | --- | --- |
| I1 | The Python arc certifier cannot issue a positive certificate in the generator's operating regime | `engagement.py:40` `AUDIT_ARC_STEP_FRACTION = 0.05` | Positive verdict possible only for trochoid radius < 0.357·R at cap π/2 (< 1.136·R at cap π); trochoid radii are ≥ R. On a 12×8 pocket at cap 90°: **0 of 126** circular cut ops certified. All three tests exercising it use `radius=0.05` against `tool_radius=0.5` — 10× smaller than the tool. |
| I2 | `EngagementReport.max_tea` can report 720° | `engagement_2.cpp:267-282` | Op 190 of the 12×8 audit, station `(0.501, 3.066667)` → `total_tea = max_run_tea = 4π` while stations 0.017 away read 0.036 rad. Duplicate rim sub-arcs are merged twice; run assembly has no `≤ 2π` invariant. Decision direction stays conservative (a doubled run is a superset). |
| I3 | Safety-critical constant duplicated across the language boundary, no contract test | `engagement_2.cpp:476` + `engagement.py:58`, `:119` | `TEA_GUARD_SAFETY_FACTOR` and `tea_growth_bound` exist twice; the comment says "mirroring"; nothing enforces it. |
| I4 | Geometry parameters unvalidated at the seam where angle parameters are meticulous | `engagement_2.cpp:558`, `:589` | `engagement_at(tool_radius=-1.0)` → `(2π, 2π, True)`. `certify_segment_tea(tool_radius=0)` → 25 stations, no error. Non-finite input leaks `RuntimeError: Cannot convert a non-finite number to an integer`. Every `Stock2::subtract_*` *does* check `radius > 0`. Also unvalidated: negative `radial_clearance`, `min > max_trochoid_radius`, `clearance_z ≤ cut_z` (silently degrades to flat linking, `toolpath.cpp:996`). |
| I5 | `link_paths=False` silently emits the traverse `link_paths=True` refuses | `toolpath.cpp:1071` | With linking on and no clearance plane a gouging flat link throws. With linking off no link primitive is emitted, yet `tessellate_operations` concatenates every op into one continuous polyline — the same through-material traverse appears in `ToolpathResult.polyline`, unmarked and uncertified. |
| I6 | The load-bearing shared-root precondition is compiled out of every shipped wheel | `engagement_2.cpp:89` | `as_radpoint` collapses a point's two coordinates to one root; violation silently corrupts every downstream orientation and chord sign. Guarded by `CGAL_assertion`; `pyproject.toml` sets `cmake.build-type = "Release"` → `NDEBUG`. `sign_mixed_radical` already handles two distinct roots. |
| I7 | Two meanings of "exact" in one codebase | `toolpath.h:14,36`, `toolpath.cpp:637`, `toolpath.py:15` | `Boundary` runs on **Epick**, where `K::FT` is `double` and `squared_distance` is an inexact *construction*; `segment_clear` compares two rounded doubles to take a gouge **decision**. Under `docs/exactness.md`'s own deciding/reporting split that is a decision on inexact quantities. The intended sense is "not interpolated". |
| I8 | `certify_bridges` can exit with an uncertified bridge | `toolpath.cpp:465` | The loop breaks when nothing was inserted; if the 8th (last) round *does* insert, the refined run is never re-tested. Structurally masked while `mat_scale ≤ 1`, but the loop's own contract does not hold. |
| I9 | Hole validation is incomplete for the straight-skeleton precondition | `toolpath.cpp:~560` | Only hole *vertices* are tested against the outer boundary, so a hole edge crossing a concave region passes; nested holes pass (pairwise test is edge-edge only). Cost O(H²E²). A violated precondition into `create_interior_straight_skeleton_2` is UB. |
| I10 | `isolines.py` — 169 lines of new public API, zero tests, different standard | `src/compas_cgal/isolines.py` | No `tests/test_isolines.py`. `from typing import List` beside PEP-604 modules; bare `raise ValueError("provide isovalues or n")`; unexplained `threshold: float = 2.0`. It *is* in the nav and API reference. |
| I11 | `compas_cgal.stock` and `compas_cgal.engagement` are absent from the docs | `mkdocs.yml:175-192`, `docs/api/` | 658 lines carrying the branch's thesis, with no API page and not in the nav. Conversely `site/plans/` publishes internal planning and session-state documents; there is no `not_in_nav` / `validation` config. |

---

## Notes

- 20-parameter positional C++ signature returning a 9-tuple of parallel matrices; semantics encoded
  as float64 columns (`meta(i,1) ∈ {0,1,2}`, `meta(i,3) ∈ 0..5`), with Python re-deriving the enum
  from `int(meta[i,3])`.
- `Station{run[i].center, -1.0, -1.0}` as an in-band split marker (`toolpath.cpp:461`).
- `external_tangent` decomposes kernel objects into raw doubles and hand-rolls the perpendicular
  (`vx=-uy, vy=ux`) — the exact pattern `CLAUDE.md` lists under "anti-patterns to REJECT in review";
  `Vector_2::perpendicular` exists.
- `_polygon_to_ccw_vertices` duplicated verbatim in `toolpath.py:107` and `stock.py:19`; the stated
  reason ("the toolpath helper is private") is undercut by `stock.py` importing
  `InvalidPolygonError` from `toolpath.py`, which also inverts the layering.
- No `__all__` in `stock.py` / `engagement.py`; `ToolpathOperation` / `ToolpathResult` mutable while
  `OperationEngagement` / `EngagementReport` are `frozen=True`.
- `Stock2::set()` returns a non-const `Gps&` and `Stock.raw` exposes it — the exact invariants are
  externally mutable.
- `_unmeasured()` returns `cap_certified=True`, so "certified" and "not applicable" are the same
  value and `cap_violations` under-counts by design.
- Bare `assert isinstance(...)` for invariants (stripped under `-O`).
- Undocumented `std::max(4, …)` in `subtract_arc_sweep` (`stock_2.cpp:170`).
- `CHAIN_SLACK_FRACTION`'s rationale cites `radial_clearance`'s *default from another module*
  (`stock_2.cpp:88`), a coupling nothing enforces — `radial_clearance=0` is accepted.
- The disk-chain under-covering certificate omits the double-precision placement of chain centres
  (`cx + guide_r*cos(a)`); the perpendicular direction has zero margin by construction, so the error
  sits in the over-covering (unsafe) direction at ~1e-16 relative.
- `ruff format --check` fails on `toolpath.py:146` (implicit string concatenation).
- `--doctest-glob=*.rst` survives in `pyproject.toml` on a repo whose doctrine bans reStructuredText.

---

## Non-goals

- Rewriting the exact predicate machinery (`sign_mixed_radical`, `run_exceeds_cap`,
  `pessimistic_runs`, `finish_engagement`'s decision half). It is correct; leave it alone.
- Changing the straight-skeleton-vs-MAT choice. That is branch B's subject.
- Performance work. The zone-query result stands.
