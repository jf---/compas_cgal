# ETH Audit Remediation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Close every Critical and Important finding of the 2026-08-19 ETH-grade audit of
`jf/toolpath-redesign` without damaging the exact-kernel machinery the audit found correct.

**Architecture:** Four phases in dependency order. Phase 0 makes the branch executable under the
pixi-exclusive build rule. Phase 1 fixes defects that produce wrong output *today*. Phase 2 repairs
the certifier's integrity — one source of truth for the growth guard, a falsification harness that
pins reality, then a curvature-aware bound. Phase 3 removes the duplicated Python certifier by
moving arc certification into C++. Phase 4 closes platform, contract, and documentation debt. Every
change to a certified path is **add-alongside → validate → swap**, never in-place surgery.

**Tech Stack:** C++20, CGAL 6.0.1 (Epick for `toolpath`, Epeck + `Gps_circle_segment_traits_2` for
`stock_2`/`engagement_2`), nanobind, Eigen, scikit-build-core, pixi, pytest + Hypothesis, ruff.

**Spec:** `docs/superpowers/specs/2026-08-19-eth-audit-findings.md`

**Audit head:** `33cbcb6`. **Plan base:** `ebf0dcb`. The four commits in between add the
`benchmarks/` corpus and `stock_2` instrumentation only; **no anchor cited in the spec moved**
(verified: `toolpath.cpp`, `engagement_2.cpp`, `toolpath.py`, `engagement.py` untouched).

---

## Global Constraints

Every task's requirements implicitly include this section.

- **pixi exclusively.** Never `pip install`, never `conda`, never a bare `python` or `pytest`. Use
  `pixi run <task>`. This branch has no manifest yet — **Task 0 is a hard prerequisite for every
  other task.**
- **Do not work in the branch-B worktree.** `/Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-exact-certified-adaptive-phase1-t9-zero-guide`
  is on `codex/exact-certified-adaptive-phase1-t9-zero-guide` and its pixi tasks import
  `_containment_2` / `_coverage_2` / `_medial_axis_2`, which do not exist here. All work in this
  plan happens in the main checkout on `jf/toolpath-redesign`.
- **Certified paths are add-alongside → validate → swap.** New entry point beside the untouched
  old one; a comparison test proving equivalence or a documented intended difference; only then
  swap. Never edit a certificate in place.
- **Commits require explicit user permission.** Each "Commit" step means: stage the exact files,
  show `git diff --cached --stat`, and ask. Set author *and* committer to
  `Jelle Feringa <jelleferinga@gmail.com>`. No `Co-Authored-By`, no Claude attribution.
- **Exact predicates for decisions.** In `stock_2`/`engagement_2` no epsilon, deflation, or
  `to_double` may reach a branch that changes control flow, output topology, or a certificate.
  Reporting doubles are fine and must never feed back into a decision.
- **No magic numbers.** Every new constant is a named module-level constant with its unit and
  derivation in a comment. Prefer `compas.tolerance.TOL` on the Python side.
- **Google-style docstrings** (`Args:` / `Returns:` / `Raises:`). Markdown, never reStructuredText.
- **Named exceptions**, one per failure mode. Never a bare `raise ValueError("...")` in new code.
- **`ruff` before every commit:** `pixi run lint`, and `ruff format` on files you touched.
- **`git mv` for every file move.** Never delete-and-recreate.
- **Never modify a reference test to make it pass.** Never `xfail`, never `skip`.
- Python floor for new code: the repo declares `>=3.9`; Task 12 decides and enforces the real floor.
  Until Task 12 lands, add `from __future__ import annotations` to any Python file you create.
- **THE CONFLICT-ZONE RULE (added 2026-08-23, see Re-plan below).** `src/engagement_2.{cpp,h}` and
  `src/stock_2.{cpp,h}` are the only files where this branch and
  `codex/exact-certified-adaptive-phase1-t9-zero-guide` both made substantial changes, and a rebase
  onto codex is pending. **Every additional commit to those four files widens a merge that must be
  resolved by hand.** Tasks 8 and 10 finish because they are already inside that zone; after them, no
  task may touch those files until the rebase lands. Every other file this plan edits is untouched on
  codex and costs the same before or after.

---

## Re-plan — 2026-08-23

Fourteen tasks are complete (0, 1, 2, 3, 22, 4, 5, 6, 7a, 7b, 8, 23, 9, and the spec/plan commit).
Execution then paused on a scope challenge: this branch is **not** the project's trunk.

`codex/exact-certified-adaptive-phase1-t9-zero-guide` is a sibling — the branches diverged at
`1860167` with **33 commits here and 244 there** — and it carries three toolpath generators that
machine whole pockets, against this branch's none. The decision recorded in `INTEGRATECODEX.md` is to
**rebase this branch onto codex**, because codex still ships every defect this plan repaired: the
`mat_scale` gouge, the `Stock2` seam, and the unsound engagement-audit certifier that
`benchmarks/`, `scripts/engagement_baseline.py` and `docs/exactness.md` all exercise.

Three facts from that investigation change this plan's ordering:

1. **The rib witness does not falsify codex's newer `continuous_tea_2` certifier** — its event
   partition closes the two-station blind spot by construction. But that certifier sits *beside* the
   audit path rather than replacing it, so the repairs here remain the only fix for what ships.
2. **Codex has no independent falsifier for circle certification.** Segments have one; circles verify
   only the implementation's own output. Task 10 builds exactly that, so its value rose.
3. **Codex never touched `toolpath.cpp`, `toolpath.py`, `engagement.py` or `tests/test_toolpath.py`.**
   Every task against those files ports cleanly and can run before or after the rebase.

The re-ordering and the five new tasks below follow from those. Task 11 in particular rose sharply:
codex's own `docs/benchmarks.md` records that conflating `uncertified` with `truly_exceeding` once
made a generator that had improved 2.4× look like a regression — which is precisely the distinction
Task 11 puts into the type system.

**Revised execution order.** Tasks 29–35 were added independently against codex trunk; none touches
the conflict zone, so all are post-rebase by construction. Two of them reorder everything else:

- **32 goes first.** The full suite has never been observed green — it timed out at ten minutes — and
  four tests fail on trunk. Nothing downstream is trustworthy until that is resolved, including any
  claim this plan's own tasks make about not breaking anything.
- **34 gates 35, and is larger than it looks.** `SegmentSiteMedialAxis.build` accepts **exactly one
  polygon** — a hole-free six-vertex ring whose vertices must be exactly
  `(0,0) (6,0) (6,2) (2,2) (2,6) (0,6)`; every other polygon raises
  `UnsupportedCanonicalMatLShapeGraphError`. Held's examples are arbitrary contours, so 35 cannot
  start until 34 lands.

| when | tasks | rationale |
| --- | --- | --- |
| Before the rebase (conflict zone) | **8** (fix round 2), **10** | already inside the two conflicting files; finishing costs less than leaving them half-done |
| Post-rebase — foundation | **32** | observe the suite green before trusting any later measurement |
| Post-rebase — soundness | **24** | the circle falsifier, using the witnesses the rebase carries |
| Post-rebase — gate | **30**, **29**, **25** | register the generators, add the second cap so a defect can be attributed, then pin the four quoted numbers |
| Post-rebase — reporting | **11** | the `uncertified` / `truly_exceeding` distinction the gate depends on, in the type system |
| Post-rebase — generator hygiene | **14**, **13**, **16**, **17**, **19**, **26** | all on `toolpath.cpp` / the engagement generators; untouched on codex, so they port free |
| Post-rebase — capability | **34**, **31**, **33** | generalise the MAT beyond one polygon; then corner clearing for the residue; then the missing physics KPIs |
| Post-rebase — platform and docs | **12**, **15**, **18**, **27**, **21**, **20**, **28** | lower stakes; **12** is worse on codex, where `adaptive/` imports `typing.Self` (3.11+) against `requires-python = ">=3.9"` |
| Capstone | **35** | reproduce Held's published examples with a per-example KPI delta — the strongest external validation available, and blocked on **34** |
| Re-scope, not rewrite | **23** | `docs/engagement_certificate.md` says the arc certifier is unrepaired, which Task 10 fixes, and it must be reconciled with codex's `docs/continuous_engagement.md` |

!!! warning "Two Held-parity items still have no task"

    `INTEGRATECODEX.md` ranks **breaking the ~141.8° cap saturation** first and **helical or ramped
    entry** second on the critical path. Neither is a task in this plan.

    Saturation is the item that decides whether this is engagement control or a guardrail — asking
    for 20° yields the same measured maximum as asking for 140°. Its cause is named in codex's own
    `docs/benchmarks.md`: the generator regulates advance but not trochoid radius, so once
    medial-axis clearance fixes the radius no advance reaches a tight cap. `radius_regulated_toolpath`
    was built for exactly this and is inert at a 120° cap, which points at the guide's station set —
    and therefore at **Task 34**, since a MAT restricted to one polygon cannot offer a different
    station set to choose from.

    Entry is cheaper and clears two of the six failing gate criteria on its own: every cap exceedance
    on all three gate pockets is a chain-entry cut, verified by set equality.

    Both need briefs before the capability phase, or 35 will reproduce Held's contours at a cap that
    does not bind.
| Machining-quality track (Tasks 29-35) | **32** (nothing is trustworthy without a green run), **29**, **30**, **31**, **33**, **35**, **34** |

The machining-quality track runs beside the certifier track and shares no files with it. Its ordering
is forced in two places only: **29 before 30**, so the new generators land on a gate axis that can
attribute; and **30 before 31**, so the uncut regression is measured by the gate before it is
repaired. **34** is ranked last deliberately — it is the largest task in the plan and its usual
justification is falsified; read its warning box before starting it.

### Landed on codex trunk — 2026-08-23

Four commits, `d9d0e8a` `7af7a24` `8d9adc6` `b62c00a`. Two findings worth carrying forward, both of
them defects in the *instrument* rather than in a generator:

1. **`max_loop_radius_step` compared radii across a retract** — a jump the cutter performs at
   clearance height. Scoped to runs bounded by a rapid *or* a `path_index` change; both boundaries are
   needed, because a chain-ordering generator links chains at cutting depth and has no rapid to split
   on. `rect_20x12` went `3.988` -> `1.024`. The criterion had never been failed by a real defect.
2. **The engagement angle is the RIM ARC, twice the textbook milling angle** — measured against the
   kernel at ratio exactly `2.000` over every radial depth. `a_e = r(1 - cos(theta/2))` was correct
   all along; a confident "correction" to the textbook relation would have reported a plunge as
   **zero** immersion and halved every load figure on the page. The conversion is now explicit in
   `textbook_engagement_deg` and anchored by
   `test_the_kernel_anchors_the_rim_arc_to_radial_depth_relation`, which reads the live kernel — so it
   cannot be "corrected" back on algebra alone, because it already was once.

Both were caught the same way, and neither by reading code: **name the physical event the comparison
corresponds to, then measure which quantity the inputs actually carry.** Anchoring a formula at its
endpoints is the weakest available evidence — it is exactly what a wrong-but-plausible formula is
most likely to satisfy. `docs/machining_metric_validity.md` records all three cases.

---

## File Structure

| File | Responsibility | Tasks |
| --- | --- | --- |
| `pyproject.toml` | pixi workspace, dependencies, tasks; pytest config | 0, 12, 20 |
| `src/engagement_2.h` / `.cpp` | Exact TEA query + motion certificates | 3, 5, 7, 8, 9, 10 |
| `src/stock_2.h` / `.cpp` | Exact stock model; **new:** min-subtraction-radius tracking | 7 |
| `src/toolpath.h` / `.cpp` | Epick generator, leads/links, tessellation | 1, 2, 4, 13, 14, 16, 17, 19 |
| `src/compas_cgal/toolpath.py` | Typed generator façade | 1, 2, 4, 12, 18, 19, 20 |
| `src/compas_cgal/stock.py` | Typed stock façade | 12, 18, 20 |
| `src/compas_cgal/engagement.py` | Audit replay; **shrinks** as certification moves to C++ | 5, 10, 11, 12, 20 |
| `src/compas_cgal/_validation.py` | **New.** Shared polygon → CCW vertex conversion | 18 |
| `src/compas_cgal/isolines.py` | Isoline extraction | 12, 15 |
| `tests/test_growth_bound.py` | **New.** Falsification harness for the analytic guard | 6, 7 |
| `tests/test_isolines.py` | **New.** Isoline coverage | 15 |
| `docs/api/compas_cgal.{stock,engagement}.md` | **New.** API reference pages | 20 |
| `mkdocs.yml` | Nav + `not_in_nav` | 20 |
| `.github/workflows/build.yml` | Python matrix | 12 |

---

# Phase 0 — Make the branch executable

### Task 0: pixi manifest for `jf/toolpath-redesign`

`CLAUDE.md` mandates pixi exclusively, but this checkout has no manifest, so `pixi run` fails and
every other task in this plan is unrunnable. The branch-B manifest cannot be reused: its
`_editable-rebuild` task imports modules that do not exist here.

**Files:**
- Modify: `pyproject.toml` (append pixi tables after the `[tool.cibuildwheel]` section)

**Interfaces:**
- Consumes: nothing.
- Produces: `pixi run baseline`, `pixi run pytest`, `pixi run affected`, `pixi run lint`,
  `pixi run docs`. Every later task uses these and only these.

- [ ] **Step 1: Write the failing check**

Create `tests/test_environment.py`:

```python
from __future__ import annotations

import tomllib
from pathlib import Path


def test_pixi_manifest_declares_the_tasks_this_repo_documents():
    """CLAUDE.md mandates `pixi run baseline|pytest|affected|lint`; the manifest must provide them."""
    manifest = tomllib.loads((Path(__file__).resolve().parents[1] / "pyproject.toml").read_text())
    tasks = manifest["tool"]["pixi"]["tasks"]
    assert {"baseline", "pytest", "affected", "lint"} <= set(tasks)


def test_editable_rebuild_only_imports_modules_this_branch_builds():
    """The rebuild probe must not reference branch-B extensions (_containment_2 etc.)."""
    manifest = tomllib.loads((Path(__file__).resolve().parents[1] / "pyproject.toml").read_text())
    probe = manifest["tool"]["pixi"]["tasks"]["_editable-rebuild"]
    for absent in ("_containment_2", "_coverage_2", "_medial_axis_2"):
        assert absent not in probe
    assert "_stock_2" in probe and "_toolpath" in probe
```

- [ ] **Step 2: Run it to verify it fails**

Run: `python3 -m pytest tests/test_environment.py -v`
(Bare python is permitted for this one step only — pixi does not exist yet.)
Expected: FAIL with `KeyError: 'pixi'`.

- [ ] **Step 3: Append the manifest to `pyproject.toml`**

```toml
# ============================================================================
# pixi
# ============================================================================

[tool.pixi.workspace]
channels = ["conda-forge"]
platforms = ["osx-arm64", "linux-64"]

[tool.pixi.dependencies]
python = "3.12.*"
cmake = ">=3.15"
ninja = "*"
compas = "*"
numpy = "*"
scipy = "*"

[tool.pixi.pypi-dependencies]
compas-cgal = { path = ".", editable = true }
nanobind = ">=1.3.2"
scikit-build-core = ">=0.10"
build = "*"
pytest = ">=7"
pytest-xdist = "*"
pytest-testmon = "*"
hypothesis = "*"
ruff = "*"

[tool.pixi.pypi-options]
no-build-isolation = ["compas-cgal"]

[tool.pixi.feature.docs.dependencies]
python = "3.12.*"
compas = "*"

[tool.pixi.feature.docs.pypi-dependencies]
markdown-callouts = ">=0.4"
mkdocs-material = ">=9.5"
mkdocstrings = { version = "*", extras = ["python"] }

[tool.pixi.feature.docs.tasks]
docs = "mkdocs build --strict"

[tool.pixi.environments]
docs = { features = ["docs"], no-default-feature = true }

[tool.pixi.tasks]
# Import probe: forces scikit-build-core's editable hook to rebuild the extensions
# before any test runs. Lists only the modules THIS branch builds.
_editable-rebuild = '''python -c "from compas_cgal import _stock_2, _toolpath"'''
pytest = { cmd = '''editable_build_dir="$(python -c 'import sys; print(next(f.path for f in sys.meta_path if hasattr(f, "known_wheel_files") and "compas_cgal._stock_2" in f.known_wheel_files))')" && SKBUILD_EDITABLE_SKIP="$editable_build_dir" pytest''', depends-on = ["_editable-rebuild"] }
baseline = { cmd = '''editable_build_dir="$(python -c 'import sys; print(next(f.path for f in sys.meta_path if hasattr(f, "known_wheel_files") and "compas_cgal._stock_2" in f.known_wheel_files))')" && SKBUILD_EDITABLE_SKIP="$editable_build_dir" pytest tests -n auto -q''', depends-on = ["_editable-rebuild"] }
affected = { cmd = '''editable_build_dir="$(python -c 'import sys; print(next(f.path for f in sys.meta_path if hasattr(f, "known_wheel_files") and "compas_cgal._stock_2" in f.known_wheel_files))')" && SKBUILD_EDITABLE_SKIP="$editable_build_dir" pytest tests -n auto --testmon -q''', depends-on = ["_editable-rebuild"] }
lint = "ruff check src/compas_cgal tests benchmarks"
```

- [ ] **Step 4: Install and verify**

Run: `pixi install && pixi run pytest tests/test_environment.py -v`
Expected: PASS, 2 tests.

- [ ] **Step 5: Establish the green baseline**

Run: `pixi run baseline`
Expected: the full suite passes. **Record the exact pass count in the commit message** — every
later task compares against it.

- [ ] **Step 6: Commit** (stage, show, ask)

```bash
git add pyproject.toml tests/test_environment.py
git commit -m "build: pixi manifest for jf/toolpath-redesign"
```

---

# Phase 1 — Defects that produce wrong output today

### Task 1: `mat_scale` contract — close the gouge (spec C1)

`mat_scale > 1` scales the trochoid radius past the available clearance. The trochoid circles are
never certified, so nothing catches it: measured penetration is +2.25 units past the wall at
`mat_scale=1.5` with a ⌀1.0 tool, and the tool *centre* leaves the pocket by 0.40.

**Files:**
- Modify: `src/toolpath.cpp:602-615` (`validate_toolpath_params`), `:772`, `:931` (both call sites)
- Modify: `src/compas_cgal/toolpath.py:282`, `:379` (docstrings for `mat_scale`)
- Test: `tests/test_toolpath.py`

**Interfaces:**
- Consumes: nothing.
- Produces: `validate_toolpath_params(double tool_diameter, double stepover, double pitch,
  double min_trochoid_radius, double max_trochoid_radius, double mat_scale,
  double radial_clearance, int samples_per_cycle, int max_passes)` — two new parameters, in that
  position. Task 4 extends the same function; do not reorder afterwards.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_toolpath.py`:

```python
@pytest.mark.parametrize("bad_scale", [1.0001, 1.5, 2.0, 10.0, 0.0, -1.0])
def test_mat_scale_outside_unit_interval_is_rejected(bad_scale):
    """mat_scale > 1 scales the trochoid radius past the available clearance and gouges.

    The trochoid circles are not certified (only bridges and leads are), so the
    gouge-free guarantee holds *by construction* and only for mat_scale <= 1.
    That precondition must be enforced at the seam, not assumed.
    """
    with pytest.raises(ValueError, match="mat_scale"):
        trochoidal_mat_toolpath(SQUARE, tool_diameter=1.0, mat_scale=bad_scale)


@pytest.mark.parametrize("mat_scale", [0.25, 0.5, 1.0])
def test_no_gouge_circular_primitives_across_admissible_mat_scale(mat_scale):
    """clearance(center) >= trochoid radius + tool radius for every admissible mat_scale."""
    tool_radius = 0.5
    polygon = _dumbbell(2.4)
    poly_xy = [list(pt[:2]) for pt in polygon.points]
    result = trochoidal_mat_toolpath_circular(
        polygon, tool_diameter=1.0, pitch=0.75, clearance_z=3.0, mat_scale=mat_scale
    )
    cut_arcs = [op for op in result.operations if op.operation == "cut" and isinstance(op.geometry, (Arc, Circle))]
    assert len(cut_arcs) > 0
    for op in cut_arcs:
        c = op.geometry.frame.point
        clearance = _distance_to_polygon_boundary_xy([float(c[0]), float(c[1])], poly_xy)
        assert clearance + GOUGE_TOL >= op.geometry.radius + tool_radius


@pytest.mark.parametrize("mat_scale", [0.25, 0.5, 1.0])
def test_no_gouge_tool_centre_across_admissible_mat_scale(mat_scale):
    """Every tessellated tool-centre point keeps exact wall clearance >= tool radius."""
    tool_radius = 0.5
    polygon = _dumbbell(2.4)
    poly_xy = [list(pt[:2]) for pt in polygon.points]
    paths = trochoidal_mat_toolpath(
        polygon, tool_diameter=1.0, pitch=0.75, samples_per_cycle=64, mat_scale=mat_scale
    )
    assert len(paths) > 0
    worst = max(tool_radius - _distance_to_polygon_boundary_xy(pt[:2].tolist(), poly_xy) for path in paths for pt in path)
    assert worst <= GOUGE_TOL
```

- [ ] **Step 2: Run to verify they fail**

Run: `pixi run pytest tests/test_toolpath.py -k mat_scale -v`
Expected: `test_mat_scale_outside_unit_interval_is_rejected` FAILS — `DID NOT RAISE ValueError`.

- [ ] **Step 3: Enforce the contract in C++**

In `src/toolpath.cpp`, extend `validate_toolpath_params`:

```cpp
void
validate_toolpath_params(
    double tool_diameter, double stepover, double pitch,
    double min_trochoid_radius, double max_trochoid_radius,
    double mat_scale, double radial_clearance,
    int samples_per_cycle, int max_passes)
{
    if (tool_diameter <= 0.0) throw std::invalid_argument("tool_diameter should be positive.");
    if (stepover <= 0.0) throw std::invalid_argument("stepover should be positive.");
    if (pitch <= 0.0) throw std::invalid_argument("pitch should be positive.");
    if (min_trochoid_radius < 0.0) throw std::invalid_argument("min_trochoid_radius should be >= 0.");
    if (max_trochoid_radius < 0.0) throw std::invalid_argument("max_trochoid_radius should be >= 0.");
    // GOUGE-FREEDOM PRECONDITION. radius_from_clearance yields
    // r = mat_scale * (clearance - R - radial_clearance); a trochoid circle reaches
    // r + R from its centre. The circles are NOT certified against the boundary (only
    // bridges and leads are), so gouge-freedom holds by construction and ONLY while
    // mat_scale <= 1. Enforced here, at the one seam, rather than assumed.
    if (!(mat_scale > 0.0 && mat_scale <= 1.0))
        throw std::invalid_argument("mat_scale should be in (0, 1]: values above 1 scale the trochoid radius past the certified clearance.");
    if (!(samples_per_cycle >= 4)) throw std::invalid_argument("samples_per_cycle should be at least 4.");
    if (max_passes <= 0) throw std::invalid_argument("max_passes should be positive.");
}
```

Update both call sites (`:772`, `:931`) to pass `mat_scale, radial_clearance`. Leave
`radial_clearance` unvalidated here — Task 4 owns it.

- [ ] **Step 4: Document the bound where the caller reads it**

In `src/toolpath.h`, replace the `@param mat_scale` line:

```
 * @param mat_scale Scale factor for clearance-derived radius; must be in (0, 1].
 *        Values above 1 push the trochoid circle past the certified clearance.
```

In `src/compas_cgal/toolpath.py`, in **both** docstrings (`:312-313` and `:418-419`):

```
    mat_scale
        Scale factor applied to the clearance-derived available radius, in
        ``(0, 1]``. The trochoid circles are gouge-free by construction rather
        than by certification, and that construction requires
        ``mat_scale <= 1``; larger values raise `ValueError`.
```

- [ ] **Step 5: Run to verify they pass**

Run: `pixi run pytest tests/test_toolpath.py -k mat_scale -v`
Expected: PASS, 9 tests.

- [ ] **Step 6: Full suite**

Run: `pixi run baseline`
Expected: the Task 0 baseline count + 9.

- [ ] **Step 7: Commit** (stage, show, ask)

```bash
git add src/toolpath.cpp src/toolpath.h src/compas_cgal/toolpath.py tests/test_toolpath.py
git commit -m "fix(toolpath): mat_scale bounded to (0,1] — gouge-free holds by construction only there"
```

---

### Task 2: certify the inter-path traverse regardless of `link_paths` (spec I5)

With `link_paths=True` and no clearance plane, a gouging flat link raises. With `link_paths=False`
no link primitive is emitted at all — yet `tessellate_operations` concatenates every operation into
one continuous polyline, so the identical through-material traverse appears in
`ToolpathResult.polyline`, unmarked and uncertified. `link_paths` must govern whether a link
*operation is recorded*, never whether the traverse is *checked*.

**Files:**
- Modify: `src/toolpath.cpp:1041-1071` (the `link_paths` branch of the connect block)
- Test: `tests/test_toolpath.py`

**Interfaces:**
- Consumes: `Boundary::segment_clear(const Segment_2&, double) const` (existing).
- Produces: no signature change. Behaviour change: the flat-traverse gouge check now runs on both
  branches.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_toolpath.py`:

```python
def test_unlinked_paths_still_certify_the_implied_traverse():
    """link_paths=False must not smuggle an uncertified through-material traverse.

    The returned polyline concatenates every operation, so consecutive paths are
    joined by an implied straight move whether or not a link primitive is recorded.
    That move is a cutting-height motion and must meet the same wall-clearance bar
    as an explicit link.
    """
    polygon = _dumbbell(1.2)  # pinched waist: cross-path traverses hug the notches
    with pytest.raises(ValueError, match="gouge"):
        trochoidal_mat_toolpath_circular(
            polygon, tool_diameter=1.0, pitch=0.75, link_paths=False, optimize_order=False
        )


def test_unlinked_paths_with_clearance_plane_do_not_raise():
    """A clearance plane lifts the traverse out of the material; no certification needed."""
    result = trochoidal_mat_toolpath_circular(
        _dumbbell(1.2), tool_diameter=1.0, pitch=0.75, link_paths=False, clearance_z=3.0
    )
    assert len(result.operations) > 0
```

- [ ] **Step 2: Run to verify it fails**

Run: `pixi run pytest tests/test_toolpath.py -k unlinked -v`
Expected: `test_unlinked_paths_still_certify_the_implied_traverse` FAILS — `DID NOT RAISE`.

- [ ] **Step 3: Certify on both branches**

Replace the `else` branch at `src/toolpath.cpp:1071` (currently
`} else { cur_xy = lead_in_pt; cur_z = cut_z; }`):

```cpp
        } else {
            // link_paths == false: no LINK primitive is recorded, but the emitted
            // polyline still runs continuously from the previous path's end to this
            // path's entry. That implied traverse is a cut-height motion and gets the
            // SAME certification as an explicit flat link -- link_paths governs what is
            // RECORDED, never what is CHECKED. Without a clearance plane there is no
            // safe fallback, so fail loud (CLAUDE.md: no silent degradation).
            if (cur_xy != lead_in_pt && !use_clearance &&
                !boundary.segment_clear(Segment_2(cur_xy, lead_in_pt), tool_radius)) {
                throw std::invalid_argument(
                    "Implied traverse between unlinked paths would gouge the boundary; "
                    "provide clearance_z for safe Z-linking.");
            }
            cur_xy = lead_in_pt; cur_z = cut_z;
        }
```

- [ ] **Step 4: Run to verify they pass**

Run: `pixi run pytest tests/test_toolpath.py -k unlinked -v`
Expected: PASS, 2 tests.

- [ ] **Step 5: Full suite**

Run: `pixi run baseline`
Expected: previous count + 2. If `test_circular_with_leads_and_links` now raises, that is a **real**
finding — report it rather than relaxing the new check.

- [ ] **Step 6: Commit** (stage, show, ask)

```bash
git add src/toolpath.cpp tests/test_toolpath.py
git commit -m "fix(toolpath): certify the implied traverse when link_paths is off"
```

---

### Task 3: geometry-parameter validation at the exact-kernel seam (spec I4, C++ half)

`engagement_at` validates `cap_chord_ratio` and `gap_close_ratio` to the ulp and never validates
`tool_radius`. A negative radius returns `(2π, 2π, cap_exceeded=True)` — a full-immersion answer for
an impossible tool. Non-finite input leaks nanobind's internal
`RuntimeError: Cannot convert a non-finite number to an integer`.

**Files:**
- Modify: `src/engagement_2.cpp:558-586` (`engagement_at`), `:589-610` (`certify_segment_tea`)
- Test: `tests/test_stock.py`

**Interfaces:**
- Consumes: nothing.
- Produces: `engagement_at` and `certify_segment_tea` raise `std::invalid_argument` (surfacing in
  Python as `ValueError`) for non-positive or non-finite `tool_radius` and non-finite coordinates.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_stock.py`:

```python
import math

import pytest
from compas.geometry import Polygon

from compas_cgal import _stock_2
from compas_cgal.stock import Stock

_CAP_RATIO_90 = 4.0 * math.sin(math.pi / 4.0) ** 2
_UNIT_SQUARE_STOCK = Polygon([(0, 0, 0), (10, 0, 0), (10, 10, 0), (0, 10, 0)])


@pytest.mark.parametrize("bad_radius", [0.0, -1.0, float("nan"), float("inf")])
def test_engagement_at_rejects_non_physical_tool_radius(bad_radius):
    """A tool radius must be positive and finite; -1.0 previously returned a 2*pi answer."""
    stock = Stock(_UNIT_SQUARE_STOCK)
    with pytest.raises(ValueError, match="tool_radius"):
        _stock_2.engagement_at(stock.raw, 5.0, 5.0, bad_radius, _CAP_RATIO_90, 0.0)


@pytest.mark.parametrize("bad_coord", [float("nan"), float("inf"), float("-inf")])
def test_engagement_at_rejects_non_finite_centre(bad_coord):
    """Non-finite station coordinates must be a named domain error, not a nanobind cast error."""
    stock = Stock(_UNIT_SQUARE_STOCK)
    with pytest.raises(ValueError, match="finite"):
        _stock_2.engagement_at(stock.raw, bad_coord, 5.0, 0.5, _CAP_RATIO_90, 0.0)


@pytest.mark.parametrize("bad_radius", [0.0, -1.0, float("nan")])
def test_certify_segment_tea_rejects_non_physical_tool_radius(bad_radius):
    stock = Stock(_UNIT_SQUARE_STOCK)
    with pytest.raises(ValueError, match="tool_radius"):
        _stock_2.certify_segment_tea(stock.raw, 1.0, 1.0, 2.0, 2.0, bad_radius, math.pi / 2)


def test_certify_segment_tea_rejects_non_finite_endpoint():
    stock = Stock(_UNIT_SQUARE_STOCK)
    with pytest.raises(ValueError, match="finite"):
        _stock_2.certify_segment_tea(stock.raw, 1.0, 1.0, float("nan"), 2.0, 0.5, math.pi / 2)
```

- [ ] **Step 2: Run to verify they fail**

Run: `pixi run pytest tests/test_stock.py -k "non_physical or non_finite" -v`
Expected: FAIL — the `tool_radius=-1.0` case returns a tuple instead of raising; the non-finite
cases raise `RuntimeError`, not `ValueError`.

- [ ] **Step 3: Validate at the boundary**

Add near the top of `src/engagement_2.cpp` (after the `using` block), and `#include <cmath>` is
already present:

```cpp
// BOUNDARY GUARD (docs/exactness.md, boundary doctrine). Doubles enter exact-land
// ONCE, here, by exact injection -- which presupposes they ARE rationals. NaN and
// +/-Inf are not, and a non-positive radius is not a cutter. Rejecting both at the
// seam keeps every downstream Epeck::FT construction total, and turns nanobind's
// internal "cannot convert a non-finite number" into a named domain error.
void require_finite(double value, const char* name)
{
    if (!std::isfinite(value))
        throw std::invalid_argument(std::string(name) + " must be finite.");
}

void require_positive_tool_radius(double tool_radius)
{
    require_finite(tool_radius, "tool_radius");
    if (!(tool_radius > 0.0))
        throw std::invalid_argument("tool_radius must be strictly positive.");
}
```

Add `#include <string>` to the include block. Then at the head of `engagement_at`, **before** the
existing `cap_chord_ratio` check:

```cpp
    require_finite(cx, "cx");
    require_finite(cy, "cy");
    require_positive_tool_radius(tool_radius);
```

and at the head of `certify_segment_tea`, before the `cap_radians` check:

```cpp
    require_finite(x0, "x0");
    require_finite(y0, "y0");
    require_finite(x1, "x1");
    require_finite(y1, "y1");
    require_positive_tool_radius(tool_radius);
    require_finite(cap_radians, "cap_radians");
```

- [ ] **Step 4: Run to verify they pass**

Run: `pixi run pytest tests/test_stock.py -k "non_physical or non_finite" -v`
Expected: PASS, 9 tests.

- [ ] **Step 5: Full suite**

Run: `pixi run baseline`

- [ ] **Step 6: Commit** (stage, show, ask)

```bash
git add src/engagement_2.cpp tests/test_stock.py
git commit -m "fix(engagement): validate tool radius and coordinate finiteness at the exact seam"
```

---

### Task 4: generator parameter contracts (spec I4, generator half)

Three parameters accept values that silently change meaning: `radial_clearance` may be negative
(shrinking the safety margin below zero), `min_trochoid_radius` may exceed `max_trochoid_radius`,
and `clearance_z <= cut_z` silently degrades to flat linking because
`use_clearance = has_clearance_z && (clearance_z > cut_z)` (`src/toolpath.cpp:996`).

**Files:**
- Modify: `src/toolpath.cpp:602` (`validate_toolpath_params`), `:996` region
- Modify: `src/compas_cgal/toolpath.py:466-476` (defaults block of `trochoidal_mat_toolpath_circular`)
- Test: `tests/test_toolpath.py`

**Interfaces:**
- Consumes: `validate_toolpath_params(..., mat_scale, radial_clearance, ...)` from Task 1.
- Produces: no signature change beyond Task 1's.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_toolpath.py`:

```python
def test_negative_radial_clearance_is_rejected():
    """radial_clearance is a safety margin; a negative one silently eats into the tool radius."""
    with pytest.raises(ValueError, match="radial_clearance"):
        trochoidal_mat_toolpath(SQUARE, tool_diameter=1.0, radial_clearance=-0.01)


def test_inverted_trochoid_radius_bounds_are_rejected():
    with pytest.raises(ValueError, match="min_trochoid_radius"):
        trochoidal_mat_toolpath(SQUARE, tool_diameter=1.0, min_trochoid_radius=2.0, max_trochoid_radius=1.0)


@pytest.mark.parametrize("clearance_z", [0.0, -1.0])
def test_clearance_plane_at_or_below_cut_plane_is_rejected(clearance_z):
    """A clearance plane that is not above the cut plane cannot lift a traverse.

    Previously this silently fell back to flat linking, which is a different and
    less safe machining strategy than the caller asked for.
    """
    with pytest.raises(ValueError, match="clearance_z"):
        trochoidal_mat_toolpath_circular(SQUARE, tool_diameter=1.0, cut_z=0.0, clearance_z=clearance_z)
```

- [ ] **Step 2: Run to verify they fail**

Run: `pixi run pytest tests/test_toolpath.py -k "radial_clearance or radius_bounds or clearance_plane" -v`
Expected: three FAILs, `DID NOT RAISE`.

- [ ] **Step 3: Enforce in C++**

In `validate_toolpath_params`, after the `max_trochoid_radius` check:

```cpp
    if (radial_clearance < 0.0)
        throw std::invalid_argument("radial_clearance should be >= 0: it is a safety margin subtracted from the available radius.");
    if (max_trochoid_radius > 0.0 && min_trochoid_radius > max_trochoid_radius)
        throw std::invalid_argument("min_trochoid_radius should not exceed max_trochoid_radius.");
```

At `src/toolpath.cpp:996`, replace the silent degradation:

```cpp
    // A clearance plane must be ABOVE the cutting plane to lift a traverse out of
    // material. Accepting one at or below it and silently reverting to flat linking
    // would give the caller a different, less safe strategy than they asked for.
    if (has_clearance_z && !(clearance_z > cut_z))
        throw std::invalid_argument("clearance_z should be strictly greater than cut_z.");
    const bool use_clearance = has_clearance_z;
    const double safe_z = has_clearance_z ? clearance_z : cut_z;
```

- [ ] **Step 4: Run to verify they pass**

Run: `pixi run pytest tests/test_toolpath.py -k "radial_clearance or radius_bounds or clearance_plane" -v`
Expected: PASS, 4 tests.

- [ ] **Step 5: Full suite**

Run: `pixi run baseline`

- [ ] **Step 6: Commit** (stage, show, ask)

```bash
git add src/toolpath.cpp tests/test_toolpath.py
git commit -m "fix(toolpath): contracts for radial_clearance, radius bounds, clearance plane"
```

---

# Phase 2 — Certifier integrity

### Task 5: one source of truth for the growth guard (spec I3)

`TEA_GUARD_SAFETY_FACTOR` and `tea_growth_bound` exist twice — `src/engagement_2.cpp:462,476` and
`src/compas_cgal/engagement.py:58,119` — with a comment saying "mirroring" and nothing enforcing it.
A change to one silently unsounds the other. Bind the C++ functions; delete the mirror.

**Files:**
- Modify: `src/engagement_2.h` (declare), `src/engagement_2.cpp:455-485` (move out of the anonymous
  namespace, bind in `register_engagement`)
- Modify: `src/compas_cgal/engagement.py:56-150`
- Test: `tests/test_engagement_audit.py`

**Interfaces:**
- Produces: `_stock_2.tea_growth_bound(d: float, r: float) -> float` and
  `_stock_2.tea_guard(d: float, r: float) -> float`; Python `_tea_growth_bound` / `_tea_guard`
  become thin forwarders keeping their existing signatures, so Task 6 and Task 10 can import either.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_engagement_audit.py`:

```python
def test_python_guard_is_the_compiled_guard_not_a_copy():
    """The audit's guard must BE the certifier's guard, not a hand-kept mirror of it.

    A mirrored safety constant that drifts turns a conservative certifier into an
    unsound one silently, so the equality is asserted rather than commented.
    """
    from compas_cgal import _stock_2
    from compas_cgal.engagement import _tea_growth_bound, _tea_guard

    for d, r in [(1e-1, 0.5), (1e-2, 0.5), (1e-3, 2.0), (1.0, 0.5), (5.0, 0.5)]:
        assert _tea_growth_bound(d, r) == _stock_2.tea_growth_bound(d, r)
        assert _tea_guard(d, r) == _stock_2.tea_guard(d, r)
```

- [ ] **Step 2: Run to verify it fails**

Run: `pixi run pytest tests/test_engagement_audit.py -k compiled_guard -v`
Expected: FAIL — `AttributeError: module 'compas_cgal._stock_2' has no attribute 'tea_growth_bound'`.

- [ ] **Step 3: Expose the compiled guard**

In `src/engagement_2.h`, after the `CertifiedTea` struct:

```cpp
// Analytic TEA-growth bound and the guard derived from it, exposed so the Python
// audit layer CALLS the certifier's guard instead of mirroring it. A mirrored
// safety constant that drifts turns a conservative certifier unsound in silence.
// REFINEMENT bound only, never a geometric decision (docs/exactness.md).
double tea_growth_bound(double d, double r);
double tea_guard(double d, double r);
```

In `src/engagement_2.cpp`, move `tea_growth_bound`, `TEA_GUARD_SAFETY_FACTOR` and `tea_guard` out of
the anonymous namespace (keep the derivation comments verbatim — they are the proof record), placing
them just above `EngagementSample engagement_at(...)` at `:558`. Add to `register_engagement`:

```cpp
    // Refinement-bound accessors: the Python audit calls these rather than mirroring
    // them, so the guard has exactly one definition (docs/exactness.md).
    m.def("tea_growth_bound", &tea_growth_bound, "d"_a, "r"_a);
    m.def("tea_guard", &tea_guard, "d"_a, "r"_a);
```

- [ ] **Step 4: Delete the Python mirror**

In `src/compas_cgal/engagement.py`, delete the `TEA_GUARD_SAFETY_FACTOR` constant (`:56-58`) and the
bodies of `_tea_growth_bound` (`:119-142`) and `_tea_guard` (`:144-146`), replacing them with:

```python
def _tea_growth_bound(d: float, r: float) -> float:
    """Conservative bound on how far a run's TEA can grow over center travel *d*.

    Forwards to the compiled certifier's own bound (`src/engagement_2.cpp`) so the
    audit and the certificate cannot drift apart. See that function's derivation
    comment for the lemma, its safe failure direction, and its known limits.

    Args:
        d: Euclidean center-travel distance (an upper bound suffices).
        r: Tool radius.

    Returns:
        Upper bound on the run's angular growth (radians).
    """
    return float(_stock_2.tea_growth_bound(d, r))


def _tea_guard(d: float, r: float) -> float:
    """Safety-scaled TEA-growth guard, as the compiled certifier computes it."""
    return float(_stock_2.tea_guard(d, r))
```

- [ ] **Step 5: Run to verify it passes**

Run: `pixi run pytest tests/test_engagement_audit.py -k compiled_guard -v`
Expected: PASS.

- [ ] **Step 6: Full suite**

Run: `pixi run baseline`
Expected: previous count + 1, and the three existing arc-certifier tests still pass unchanged
(their verdicts are bit-identical — the guard values are the same numbers, now computed once).

- [ ] **Step 7: Commit** (stage, show, ask)

```bash
git add src/engagement_2.h src/engagement_2.cpp src/compas_cgal/engagement.py tests/test_engagement_audit.py
git commit -m "refactor(engagement): single compiled source of truth for the TEA growth guard"
```

---

### Task 6: falsification harness for the growth bound (spec C2, evidence)

The certificate's soundness rests on `tea_growth_bound` being an upper bound. The audit proved it is
not, and **no test in the suite would have noticed**. This task builds the instrument. It is
expected to go RED on the concave regime — that RED is the deliverable, and Task 7 turns it green.

**Files:**
- Create: `tests/test_growth_bound.py`

**Interfaces:**
- Consumes: `_stock_2.tea_growth_bound(d, r)`, `_stock_2.engagement_at(...)` from Task 5.
- Produces: `_max_run_tea(stock, x, y, r) -> float` and `_worst_growth_over(...)` helpers, reused by
  Task 7's regression test.

- [ ] **Step 1: Write the harness**

Create `tests/test_growth_bound.py`:

```python
"""Falsification harness for the analytic TEA-growth bound.

`tea_growth_bound(d, r)` is the ONLY thing standing between the guarded-station
method and an unsound certificate: `certify_segment_tea` subtracts
`2 * tea_growth_bound(half_spacing, r)` from the cap and concludes that two
passing stations certify everything between them. If the bound under-estimates
true growth anywhere, that conclusion is false.

These tests measure true growth with the EXACT kernel and compare it against the
bound across the boundary shapes the stock model actually produces. They are
deliberately adversarial: a green run here is the evidence the certificate needs.
"""

from __future__ import annotations

import math

import pytest
from compas.geometry import Polygon

from compas_cgal import _stock_2
from compas_cgal.stock import Stock

TOOL_RADIUS = 0.5
CAP_RATIO_90 = 4.0 * math.sin(math.pi / 4.0) ** 2

# Ambient block, large enough that its own straight walls never reach the probes.
_BLOCK = Polygon([(-6, -6, 0), (6, -6, 0), (6, 6, 0), (-6, 6, 0)])


def _max_run_tea(stock: Stock, x: float, y: float, r: float = TOOL_RADIUS) -> float:
    """Exact largest contiguous engaged run at a station, in radians (reporting)."""
    _total, max_run, _exceeded = _stock_2.engagement_at(stock.raw, x, y, r, CAP_RATIO_90, 0.0)
    return float(max_run)


def _worst_growth_over(stock: Stock, x0: float, y0: float, d: float, samples: int = 361) -> float:
    """Largest TEA increase from (x0, y0) to any point exactly *d* away."""
    base = _max_run_tea(stock, x0, y0)
    worst = 0.0
    for i in range(samples):
        a = 2.0 * math.pi * i / samples
        worst = max(worst, _max_run_tea(stock, x0 + d * math.cos(a), y0 + d * math.sin(a)) - base)
    return worst


def _void_stock(rho: float) -> Stock:
    """Block with one circular void of radius *rho* centred at the origin."""
    stock = Stock(_BLOCK)
    stock.subtract_disk(0.0, 0.0, rho)
    return stock


@pytest.mark.parametrize("d", [1e-1, 1e-2, 1e-3])
def test_bound_holds_against_a_straight_wall(d):
    """The half-plane case: clause (b) is asymptotically tight here, so this is the floor."""
    stock = Stock(Polygon([(-6, -6, 0), (6, -6, 0), (6, 0, 0), (-6, 0, 0)]))
    bound = _stock_2.tea_growth_bound(d, TOOL_RADIUS)
    worst = max(_worst_growth_over(stock, 0.0, y, d) for y in (0.0, 0.1, 0.25, 0.5, 0.75))
    assert worst <= bound, f"straight wall: growth {worst:.6f} exceeds bound {bound:.6f} at d={d}"


@pytest.mark.parametrize("d", [1e-1, 1e-2, 1e-3])
@pytest.mark.parametrize("rho", [0.1, 0.25, 0.5])
def test_bound_holds_against_a_void_no_larger_than_the_tool(d, rho):
    """rho <= r: the rim can never sit inside the void, so the feature behaves convexly."""
    stock = _void_stock(rho)
    bound = _stock_2.tea_growth_bound(d, TOOL_RADIUS)
    worst = max(_worst_growth_over(stock, s, 0.0, d) for s in (0.0, rho, rho + TOOL_RADIUS, rho + TOOL_RADIUS + 0.1))
    assert worst <= bound, f"rho={rho}: growth {worst:.6f} exceeds bound {bound:.6f} at d={d}"


@pytest.mark.parametrize("d", [1e-2, 1e-3, 1e-4])
@pytest.mark.parametrize("rho_over_r", [1.01, 1.1, 1.3, 2.0, 4.0])
def test_bound_holds_against_a_concave_void_larger_than_the_tool(d, rho_over_r):
    """rho > r: the rim can sit INSIDE the void and emerge.

    This is the regime the audit falsified. At internal tangency s0 = rho - r the
    emerged run half-angle satisfies

        psi**2 = (2*d*rho + d**2) / (r * (rho - r + d))

    so growth ~ 2*sqrt(2*d*rho / (r*(rho - r))), which exceeds the current
    4*asin(d/2r) + 2*acos(1 - d/r) without limit as rho -> r+.
    """
    rho = rho_over_r * TOOL_RADIUS
    stock = _void_stock(rho)
    bound = _stock_2.tea_growth_bound(d, TOOL_RADIUS)
    worst = _worst_growth_over(stock, rho - TOOL_RADIUS, 0.0, d)
    assert worst <= bound, (
        f"concave void rho/r={rho_over_r}: growth {worst:.6f} exceeds bound {bound:.6f} "
        f"at d={d} (ratio {worst / bound:.2f}x)"
    )


@pytest.mark.parametrize("d", [1e-2, 1e-3, 1e-4])
def test_bound_holds_on_stock_the_shipped_audit_actually_produces(d):
    """Voids cut by the SAME radius as the query: the audit's own regime.

    Every `Stock2::subtract_*` removes a union of disks of the subtraction radius,
    so within `audit_toolpath_engagement` the void arcs always have rho == r. This
    test pins that the accident holds, so a future change to the audit's wiring
    that breaks it fails here rather than in a certificate.
    """
    stock = Stock(_BLOCK)
    for k in range(6):
        stock.subtract_capsule(-2.0 + k * 0.8, -1.0, -2.0 + k * 0.8, 1.0, TOOL_RADIUS)
    bound = _stock_2.tea_growth_bound(d, TOOL_RADIUS)
    worst = 0.0
    for k in range(6):
        for y in (-1.0, -0.5, 0.0, 0.5, 1.0):
            worst = max(worst, _worst_growth_over(stock, -2.0 + k * 0.8 + TOOL_RADIUS, y, d, samples=121))
    assert worst <= bound, f"depleted stock: growth {worst:.6f} exceeds bound {bound:.6f} at d={d}"
```

- [ ] **Step 2: Run and record exactly which parametrisations are RED**

Run: `pixi run pytest tests/test_growth_bound.py -v`
Expected: the straight-wall, small-void and shipped-audit cases PASS; the
`concave_void_larger_than_the_tool` cases FAIL for `rho_over_r` in `{1.01, 1.1, 1.3}` and pass for
`{2.0, 4.0}` at small `d`. **Paste the failure table into the commit message** — it is the evidence
record and Task 7's acceptance criterion.

- [ ] **Step 3: Do not fix anything yet**

The RED is the deliverable. Do **not** weaken an assertion, add a tolerance, or mark anything
`xfail` — CLAUDE.md forbids all three. Task 7 makes it green by fixing the bound.

- [ ] **Step 4: Commit the harness with its RED recorded** (stage, show, ask)

```bash
git add tests/test_growth_bound.py
git commit -m "test(engagement): falsification harness for the TEA growth bound (RED on concave voids)"
```

---

### Task 7: curvature-aware growth bound (spec C2, fix)

The existing bound is correct for straight walls and convex material and wrong for concave voids
larger than the tool. `Stock2` knows every radius it ever subtracted, so the bound can be
parameterised by the smallest one.

**Files:**
- Modify: `src/stock_2.h:44-76` (add `min_subtraction_radius()`), `src/stock_2.cpp` (track it in
  `subtract_disk` and `subtract_point_chain`)
- Modify: `src/engagement_2.h`, `src/engagement_2.cpp` (new overload + `certify_recursive` swap)
- Test: `tests/test_growth_bound.py`, `tests/test_stock.py`

**Interfaces:**
- Consumes: `Stock2::min_subtraction_radius() const -> double` (this task creates it; returns
  `std::numeric_limits<double>::infinity()` when nothing has been subtracted).
- Produces: `double tea_growth_bound(double d, double r, double rho)` — a **new overload** beside
  the untouched two-argument one, which keeps forwarding with `rho = infinity`. Task 10 uses the
  three-argument form.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_growth_bound.py`:

```python
def test_stock_reports_the_smallest_radius_it_ever_subtracted():
    """The curvature-aware bound needs the tightest concave feature the stock contains."""
    stock = Stock(_BLOCK)
    assert math.isinf(stock.min_subtraction_radius())
    stock.subtract_disk(0.0, 0.0, 0.9)
    assert stock.min_subtraction_radius() == pytest.approx(0.9)
    stock.subtract_capsule(1.0, 1.0, 2.0, 2.0, 0.4)
    assert stock.min_subtraction_radius() == pytest.approx(0.4)
    stock.subtract_disk(3.0, 3.0, 2.0)
    assert stock.min_subtraction_radius() == pytest.approx(0.4)


@pytest.mark.parametrize("d", [1e-2, 1e-3, 1e-4])
@pytest.mark.parametrize("rho_over_r", [1.01, 1.1, 1.3, 2.0, 4.0])
def test_curvature_aware_bound_covers_the_concave_regime(d, rho_over_r):
    """The three-argument bound must dominate true growth where the two-argument one fails."""
    rho = rho_over_r * TOOL_RADIUS
    stock = _void_stock(rho)
    bound = _stock_2.tea_growth_bound_curved(d, TOOL_RADIUS, stock.min_subtraction_radius())
    worst = _worst_growth_over(stock, rho - TOOL_RADIUS, 0.0, d)
    assert worst <= bound, f"rho/r={rho_over_r}: growth {worst:.6f} exceeds curved bound {bound:.6f} at d={d}"


def test_curvature_aware_bound_reduces_to_the_flat_bound_without_concave_features():
    """rho = infinity (nothing subtracted) must reproduce the straight-wall bound exactly."""
    for d, r in [(1e-1, 0.5), (1e-3, 0.5), (1e-2, 2.0)]:
        assert _stock_2.tea_growth_bound_curved(d, r, float("inf")) == _stock_2.tea_growth_bound(d, r)


def test_curvature_aware_bound_is_monotone_in_travel():
    """certify_recursive relies on monotonicity: growth over the nearest-station distance
    must be bounded by growth over the half-spacing."""
    prev = -1.0
    for d in [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.5, 1.0, 2.0]:
        cur = _stock_2.tea_growth_bound_curved(d, TOOL_RADIUS, 0.6)
        assert cur >= prev
        prev = cur
```

Then change `test_bound_holds_against_a_concave_void_larger_than_the_tool` to call
`_stock_2.tea_growth_bound_curved(d, TOOL_RADIUS, rho)` instead of the two-argument form, and rename
it to `test_two_argument_bound_is_documented_as_flat_only` with the assertion inverted for
`rho_over_r in (1.01, 1.1, 1.3)` — the flat bound is *expected* to be exceeded there, and pinning
that keeps the limitation visible:

```python
@pytest.mark.parametrize("d", [1e-2, 1e-3, 1e-4])
@pytest.mark.parametrize("rho_over_r", [1.01, 1.1, 1.3])
def test_two_argument_bound_is_flat_only_and_says_so(d, rho_over_r):
    """The two-argument bound is NOT valid for concave voids larger than the tool.

    Pinned deliberately: the flat bound survives in the codebase for straight-wall
    and convex geometry, and this test is the standing reminder that callers facing
    concave features must use `tea_growth_bound_curved`.
    """
    rho = rho_over_r * TOOL_RADIUS
    worst = _worst_growth_over(_void_stock(rho), rho - TOOL_RADIUS, 0.0, d)
    assert worst > _stock_2.tea_growth_bound(d, TOOL_RADIUS)
```

- [ ] **Step 2: Run to verify they fail**

Run: `pixi run pytest tests/test_growth_bound.py -v`
Expected: `AttributeError: ... has no attribute 'tea_growth_bound_curved'`, and
`Stock` has no `min_subtraction_radius`.

- [ ] **Step 3: Track the smallest subtraction radius**

In `src/stock_2.h`, inside `Stock2`'s public section:

```cpp
    // Smallest radius ever removed from this stock. Every subtract_* removes a union
    // of disks of its own radius, so the resulting CONCAVE boundary arcs all have that
    // radius -- and the tightest one governs how fast the rim's engaged run can grow
    // when the cutter emerges from a void (see tea_growth_bound_curved). Infinity when
    // nothing has been subtracted: the boundary is then straight, the flat-bound case.
    double min_subtraction_radius() const { return min_subtraction_radius_; }
```

and in the private section, beside `Gps set_;`:

```cpp
    double min_subtraction_radius_ = std::numeric_limits<double>::infinity();
```

Add `#include <limits>` to `src/stock_2.h`. In `src/stock_2.cpp`, record the radius in the two
places that actually construct disks — `subtract_disk` (after its `radius <= 0.0` guard) and
`subtract_point_chain` (at the top):

```cpp
    min_subtraction_radius_ = std::min(min_subtraction_radius_, radius);
```

Expose it in the `NB_MODULE` block beside the other `Stock2` methods:

```cpp
        .def("min_subtraction_radius", &Stock2::min_subtraction_radius)
```

and forward it from `src/compas_cgal/stock.py`, beside `arrangement_stats`:

```python
    def min_subtraction_radius(self) -> float:
        """Smallest radius ever removed from this stock.

        Every ``subtract_*`` removes a union of disks of its own radius, so the
        concave boundary arcs the operation leaves all carry that radius. The
        tightest of them governs the TEA-growth bound used to certify motions
        against this stock.

        Returns:
            The smallest subtraction radius, or ``inf`` if nothing was subtracted.
        """
        return float(self._raw.min_subtraction_radius())
```

- [ ] **Step 4: Add the curvature-aware bound alongside the flat one**

In `src/engagement_2.h`, beside the Task 5 declarations:

```cpp
// Curvature-aware TEA-growth bound. `rho` is the smallest CONCAVE boundary-arc
// radius present in the stock (Stock2::min_subtraction_radius(); infinity when the
// boundary is straight).
//
// WHY THE FLAT BOUND IS NOT ENOUGH. tea_growth_bound's newborn-contact term
// 2*acos(1 - d/r) assumes the biting feature is a HALF-PLANE. For a void disk of
// radius rho > r the rim can sit INSIDE the void and emerge; at internal tangency
// s0 = rho - r, travel `delta` produces an emerged run of exactly
//
//     2*psi,  psi^2 = (2*delta*rho + delta^2) / (r * (rho - r + delta))
//
// which tends to 2*sqrt(2*delta*rho / (r*(rho-r))) and grows WITHOUT BOUND as
// rho -> r+. Measured at d = 1e-4, r = 0.5: the flat bound is exceeded 1.03x at
// rho-r = 0.30r, 1.65x at 0.10r and 4.97x at 0.01r -- past the factor-2 guard.
// For a CONVEX feature of radius rho the same derivation gives
// 2*sqrt(2*delta*rho/(r*(rho+r))), which is BELOW the flat term, so convex
// geometry needs no correction. As rho -> infinity the concave term reduces to
// the flat one, so this is a strict generalisation.
//
// REFINEMENT bound only, never a geometric decision (docs/exactness.md). Safe
// failure direction: too large a bound forces extra refinement or a conservative
// "uncertified" verdict, never a false pass.
double tea_growth_bound_curved(double d, double r, double rho);
```

In `src/engagement_2.cpp`, beside `tea_growth_bound` (leave that function **untouched**):

```cpp
double tea_growth_bound_curved(double d, double r, double rho)
{
    const double drift = 4.0 * std::asin(std::min(1.0, d / (2.0 * r)));
    const double flat_newborn = 2.0 * std::acos(std::max(-1.0, 1.0 - d / r));
    // rho <= r: the rim cannot sit inside such a void, so no concave emergence is
    // possible and the flat newborn term governs.
    if (!(rho > r)) return drift + flat_newborn;
    const double psi_sq = (2.0 * d * rho + d * d) / (r * (rho - r + d));
    const double concave_newborn = 2.0 * std::sqrt(psi_sq);
    return drift + std::max(flat_newborn, concave_newborn);
}
```

Bind it in `register_engagement`:

```cpp
    m.def("tea_growth_bound_curved", &tea_growth_bound_curved, "d"_a, "r"_a, "rho"_a);
```

- [ ] **Step 5: Run the new tests**

Run: `pixi run pytest tests/test_growth_bound.py -v`
Expected: all PASS. If `test_curvature_aware_bound_covers_the_concave_regime` still fails at
`rho_over_r = 1.01`, the emergence is not the only mechanism at that curvature — **stop and report**;
do not add a fudge factor.

- [ ] **Step 6: Swap `certify_recursive` onto the curved bound**

In `src/engagement_2.cpp:500-560`, replace the two `tea_guard(half_spacing, r)` uses with the
curvature-aware guard read from the stock:

```cpp
    const double rho = stock.min_subtraction_radius();
    const double guard = TEA_GUARD_SAFETY_FACTOR * tea_growth_bound_curved(half_spacing, r, rho);
```

Leave the rest of `certify_recursive` byte-identical: the gap-closure wiring already derives
`gap_close_ratio` from `guard`, so it follows the corrected value automatically.

- [ ] **Step 7: Verify the swap changes verdicts only in the safe direction**

Run: `pixi run baseline`
Expected: all green. `stations` counts may rise (a larger guard forces more refinement) and some
previously-certified motions on concave stock may now report uncertified — both are the safe
direction. If any test asserting `cap_certified is True` now fails, record which and why before
changing anything.

- [ ] **Step 8: Commit** (stage, show, ask)

```bash
git add src/stock_2.h src/stock_2.cpp src/engagement_2.h src/engagement_2.cpp \
        src/compas_cgal/stock.py tests/test_growth_bound.py
git commit -m "fix(engagement): curvature-aware TEA growth bound closes the concave-void hole"
```

---

### Task 8: the harvest cannot report more than a full turn (spec I2)

`engagement_at` reported `total_tea = max_run_tea = 4π` (720°) at station `(0.501, 3.066667)` while
neighbours 0.017 away read 0.036 rad. Duplicate rim sub-arcs are harvested and merged twice. The
decision stays conservative (a doubled run is a superset), but the reported number is the audit's
headline metric.

**Files:**
- Modify: `src/engagement_2.cpp:390-401` (tail of `engaged_arcs_zone`), `:267-283` (run assembly in
  `finish_engagement`)
- Test: `tests/test_engagement_audit.py`

**Interfaces:**
- Consumes: `Arc { GpsPoint ccw_start; GpsPoint ccw_end; double span; }` (existing).
- Produces: no signature change.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_engagement_audit.py`:

```python
def test_reported_tea_never_exceeds_a_full_turn():
    """TEA is an angle on the cutter rim: no run can exceed 2*pi.

    Regression for a 4*pi (720 deg) report produced by duplicate rim sub-arcs
    being merged twice during run assembly.
    """
    import math

    from compas.geometry import Polygon

    from compas_cgal.engagement import audit_toolpath_engagement
    from compas_cgal.toolpath import trochoidal_mat_toolpath_circular

    polygon = Polygon([(0, 0, 0), (12, 0, 0), (12, 8, 0), (0, 8, 0)])
    result = trochoidal_mat_toolpath_circular(polygon, tool_diameter=1.0, pitch=0.75, clearance_z=3.0)
    report = audit_toolpath_engagement(polygon, result, tool_diameter=1.0, tea_cap=math.pi / 2)

    over = [(e.op_index, e.max_tea) for e in report.operations if e.max_tea > 2 * math.pi + 1e-9]
    assert not over, f"operations reporting more than a full turn of engagement: {over}"
    assert report.max_tea <= 2 * math.pi + 1e-9
```

- [ ] **Step 2: Run to verify it fails**

Run: `pixi run pytest tests/test_engagement_audit.py -k full_turn -v`
Expected: FAIL listing `(190, 12.566370614359172)`.

- [ ] **Step 3: Drop duplicate rim sub-arcs at the harvest**

At the tail of `engaged_arcs_zone` (`src/engagement_2.cpp:390-401`), replace the push with a
duplicate-rejecting insert:

```cpp
    // The zone can report the SAME rim sub-arc more than once (a curve that reaches a
    // face through more than one traversal, or an arrangement whose redundant edges
    // split it). A duplicate is harmless for the DECISION -- merging it twice only
    // enlarges a run, which is conservative -- but it double-counts the REPORTED span
    // and has produced a 4*pi TEA. Reject by EXACT endpoint equality: two rim sub-arcs
    // with identical exact endpoints ARE the same sub-arc (docs/exactness.md: adjacency
    // is exact point equality, never a tolerance).
    for (const GpsXCurve& xc : vis.engaged) {
        GpsPoint s = xc.source();
        GpsPoint t = xc.target();
        if (s == t) continue;   // tangent-touch degeneracy: zero-measure contact
        if (xc.orientation() == CGAL::CLOCKWISE) std::swap(s, t);
        const bool duplicate = std::any_of(arcs.begin(), arcs.end(), [&](const Arc& a) {
            return a.ccw_start == s && a.ccw_end == t;
        });
        if (duplicate) continue;
        const double sx = CGAL::to_double(s.x()), sy = CGAL::to_double(s.y());
        const double tx = CGAL::to_double(t.x()), ty = CGAL::to_double(t.y());
        double span = std::atan2(ty - cy, tx - cx) - std::atan2(sy - cy, sx - cx);
        if (span <= 0.0) span += 2.0 * std::numbers::pi;
        arcs.push_back({s, t, span});
    }
```

- [ ] **Step 4: Add the reporting invariant to run assembly**

In `finish_engagement`, immediately after the wrap-around merge (`src/engagement_2.cpp:283`):

```cpp
    // REPORTING INVARIANT. Spans are reporting doubles summed from atan2, so a small
    // representation slack is expected; a run beyond a full turn is not slack, it is a
    // double-counted sub-arc. Fail loud rather than publish an impossible angle. The
    // slack is one part in 1e-9 of a turn -- far above atan2 accumulation over the
    // handful of sub-arcs a rim can carry, far below any real double count (>= 1 turn).
    constexpr double FULL_TURN_REPORT_SLACK = 1e-9 * 2.0 * std::numbers::pi;
    for (const Arc& run : runs) {
        if (run.span > 2.0 * std::numbers::pi + FULL_TURN_REPORT_SLACK)
            throw std::logic_error("engaged run exceeds a full turn: duplicate rim sub-arc in the harvest.");
    }
```

- [ ] **Step 5: Run to verify it passes**

Run: `pixi run pytest tests/test_engagement_audit.py -k full_turn -v`
Expected: PASS.

- [ ] **Step 6: Full suite**

Run: `pixi run baseline`
Expected: all green. The differential-oracle tests in `tests/test_engagement_oracle.py` must be
unchanged — dedup removes double counting only, so any oracle delta is a real regression.

- [ ] **Step 7: Commit** (stage, show, ask)

```bash
git add src/engagement_2.cpp tests/test_engagement_audit.py
git commit -m "fix(engagement): reject duplicate rim sub-arcs; reported TEA cannot exceed a turn"
```

---

### Task 9: the shared-root precondition must survive `NDEBUG` (spec I6)

`as_radpoint` (`src/engagement_2.cpp:89`) collapses a point's two coordinates onto one root. If the
coordinates ever carried different roots, every downstream orientation and chord sign is silently
wrong. It is guarded by `CGAL_assertion`, and `pyproject.toml` builds `Release` — so the check does
not exist in any shipped wheel.

**Files:**
- Modify: `src/engagement_2.cpp:80-95`
- Test: `tests/test_stock.py`

**Interfaces:**
- Consumes: `RadPoint { FT x0, x1, y0, y1, root; }` (existing).
- Produces: no signature change; `as_radpoint` may now throw `std::logic_error`.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_stock.py`:

```python
def test_release_build_still_enforces_the_shared_root_precondition():
    """The exact cap predicate assumes both coordinates of a point share one root.

    That assumption is load-bearing for every certificate, and CGAL_assertion is
    compiled out of the Release wheels this project ships. The check must be a real
    runtime check. Exercised indirectly: a heavily depleted stock generates every
    circle/circle and circle/line crossing kind the harvest can meet, and none of
    them may trip the guard.
    """
    import math

    from compas.geometry import Polygon

    from compas_cgal import _stock_2
    from compas_cgal.stock import Stock

    stock = Stock(Polygon([(0, 0, 0), (8, 0, 0), (8, 8, 0), (0, 8, 0)]))
    for k in range(12):
        a = 2.0 * math.pi * k / 12.0
        stock.subtract_capsule(4.0, 4.0, 4.0 + 3.0 * math.cos(a), 4.0 + 3.0 * math.sin(a), 0.4)
    ratio = 4.0 * math.sin(math.pi / 4.0) ** 2
    for i in range(60):
        for j in range(60):
            _stock_2.engagement_at(stock.raw, 0.5 + i * 0.12, 0.5 + j * 0.12, 0.5, ratio, 0.0)
```

- [ ] **Step 2: Run it — it should pass already**

Run: `pixi run pytest tests/test_stock.py -k shared_root -v`
Expected: PASS. This is a **guard test**, not a RED-first test: it proves the check does not
misfire on real geometry once Step 3 makes it live. Record the runtime; if it exceeds ~60 s, halve
the grid.

- [ ] **Step 3: Make the precondition real**

Replace `as_radpoint`'s body:

```cpp
RadPoint as_radpoint(const GpsPoint& p)
{
    const CoordNT& X = p.x();
    const CoordNT& Y = p.y();
    FT root(0), x1(0), y1(0);
    if (X.is_extended()) { root = X.root(); x1 = X.a1(); }
    if (Y.is_extended()) {
        // LOAD-BEARING PRECONDITION, checked at runtime rather than asserted. CGAL
        // builds every circle/line and circle/circle intersection point from a single
        // shared discriminant, so both coordinates of one point carry one root. If that
        // ever failed, collapsing to a single `root` here would silently corrupt every
        // orientation and chord sign the certificate rests on -- and CGAL_assertion is
        // compiled out of this project's Release wheels (pyproject.toml build-type).
        if (X.is_extended() && Y.root() != root)
            throw std::logic_error("as_radpoint: point coordinates carry different roots; the shared-root precondition of the exact cap predicate is violated.");
        root = Y.root();
        y1 = Y.a1();
    }
    return { X.a0(), x1, Y.a0(), y1, root };
}
```

- [ ] **Step 4: Run to verify it still passes**

Run: `pixi run pytest tests/test_stock.py -k shared_root -v`
Expected: PASS — the guard is live and does not misfire.

- [ ] **Step 5: Full suite**

Run: `pixi run baseline`

- [ ] **Step 6: Commit** (stage, show, ask)

```bash
git add src/engagement_2.cpp tests/test_stock.py
git commit -m "fix(engagement): enforce the shared-root precondition in release builds"
```

---

# Phase 3 — Audit fidelity

### Task 10: adaptive arc certification in C++ (spec I1, remainder of I3)

`_certify_arc_engagement` is a second certifier written in Python at a fixed 20 stations per turn.
At that density a positive verdict is possible only for trochoid radius < 0.357·R at a 90° cap —
smaller than the tool. Measured end-to-end: **0 of 126** circular cut operations certified. The fix
is not to tune the constant; it is to have one certifier.

**Files:**
- Modify: `src/engagement_2.h`, `src/engagement_2.cpp` (new `certify_arc_tea` beside the untouched
  `certify_segment_tea`)
- Modify: `src/compas_cgal/engagement.py:227-320` (swap `_certify_arc_engagement` onto the binding)
- Test: `tests/test_engagement_audit.py`

**Interfaces:**
- Consumes: `engagement_at`, `tea_growth_bound_curved` (Task 7), `Stock2::min_subtraction_radius()`.
- Produces: `CertifiedTea certify_arc_tea(const Stock2&, double cx, double cy, double radius,
  double start_angle, double sweep, bool cw, double tool_radius, double cap_radians)` bound as
  `_stock_2.certify_arc_tea(stock, cx, cy, radius, start_angle, sweep, cw, tool_radius, cap_radians)
  -> (max_tea, cap_certified, stations)`. `sweep` is the signed-magnitude swept angle in `(0, 2π]`;
  `cw` gives the direction. Angles in radians.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_engagement_audit.py`:

```python
def test_arc_certifier_certifies_ordinary_trochoid_circles():
    """A trochoid circle of tool-scale radius in cleared material must be certifiable.

    At the previous FIXED 20-stations-per-turn density the guard exceeded any cap
    for trochoid radii above ~0.36 * tool radius, so every circular cut operation
    the generator emits was reported uncertified regardless of its true engagement.
    """
    import math

    from compas.geometry import Polygon

    from compas_cgal import _stock_2
    from compas_cgal.stock import Stock

    stock = Stock(Polygon([(0, 0, 0), (10, 0, 0), (10, 10, 0), (0, 10, 0)]))
    # Clear a generous pocket, then trace a tool-scale circle inside the void:
    # engagement is zero everywhere, so the ONLY thing that can refuse is the guard.
    stock.subtract_disk(5.0, 5.0, 3.0)
    max_tea, certified, stations = _stock_2.certify_arc_tea(
        stock.raw, 5.0, 5.0, 1.0, 0.0, 2 * math.pi, True, 0.5, math.pi / 2
    )
    assert certified is True
    assert max_tea == pytest.approx(0.0, abs=1e-9)
    assert stations >= 4


def test_arc_certifier_refuses_a_genuinely_over_cap_circle():
    """A circle cutting virgin material at full immersion must be refused."""
    import math

    from compas.geometry import Polygon

    from compas_cgal import _stock_2
    from compas_cgal.stock import Stock

    stock = Stock(Polygon([(0, 0, 0), (10, 0, 0), (10, 10, 0), (0, 10, 0)]))
    _max_tea, certified, _stations = _stock_2.certify_arc_tea(
        stock.raw, 5.0, 5.0, 1.0, 0.0, 2 * math.pi, True, 0.5, math.pi / 2
    )
    assert certified is False


def test_audit_certifies_a_nonzero_share_of_circular_operations():
    """End-to-end: the audit's circular verdicts must reflect the toolpath, not its own density."""
    import math

    from compas.geometry import Arc, Circle, Polygon

    from compas_cgal.engagement import audit_toolpath_engagement
    from compas_cgal.toolpath import trochoidal_mat_toolpath_circular

    polygon = Polygon([(0, 0, 0), (12, 0, 0), (12, 8, 0), (0, 8, 0)])
    result = trochoidal_mat_toolpath_circular(polygon, tool_diameter=1.0, pitch=0.75, clearance_z=3.0)
    report = audit_toolpath_engagement(polygon, result, tool_diameter=1.0, tea_cap=math.pi)

    circular = [e for e in report.operations if isinstance(result.operations[e.op_index].geometry, (Arc, Circle))]
    assert len(circular) > 0
    certified = sum(1 for e in circular if e.cap_certified)
    assert certified > 0, f"0 of {len(circular)} circular operations certified — density-limited, not toolpath-limited"
```

- [ ] **Step 2: Run to verify they fail**

Run: `pixi run pytest tests/test_engagement_audit.py -k "arc_certifier_certifies_ordinary or over_cap_circle or nonzero_share" -v`
Expected: FAIL — `no attribute 'certify_arc_tea'`, and the end-to-end assertion reports `0 of 126`.

- [ ] **Step 3: Add `certify_arc_tea` beside the untouched segment certifier**

In `src/engagement_2.h`, after `certify_segment_tea`:

```cpp
// Certify TEA(P) <= cap_radians for EVERY cutter center P on the circular motion of
// radius `radius` about (cx, cy), from `start_angle` through `sweep` radians in the
// `cw` direction, against the frozen stock.
//
// Identical method to certify_segment_tea -- adaptive station bisection with the
// guarded exact test -- transposed from segment parameter to arc parameter. The
// station spacing is the CHORD between adjacent stations, which upper-bounds nothing
// and is upper-BOUNDED by the arc length, so the guard is evaluated on the arc length
// (conservative: arc >= chord, and tea_growth_bound_curved is monotone).
//
// This replaces the fixed-density Python mirror in compas_cgal/engagement.py, whose
// 20-stations-per-turn guard exceeded any cap for trochoid radii above ~0.36 * tool
// radius -- so it structurally could not certify the circles the generator emits.
//
// Raises std::invalid_argument if cap_radians is outside (0, pi], if sweep is outside
// (0, 2*pi], if radius <= 0, or if any argument is non-finite.
CertifiedTea certify_arc_tea(const Stock2& stock, double cx, double cy, double radius,
                             double start_angle, double sweep, bool cw,
                             double tool_radius, double cap_radians);
```

In `src/engagement_2.cpp`, beside `certify_recursive` (leave that function untouched), add the arc
analogue in the anonymous namespace:

```cpp
// Arc-parameter analogue of certify_recursive. `a0`/`a1` are absolute angles on the
// guide circle; the recursion bisects the ANGLE, and the guard is evaluated on the
// sub-arc's ARC LENGTH (>= the chord every interior center is measured against, so
// conservative). Structure is otherwise identical -- exact station verdicts plus the
// analytic guard, refine on failure, uncertified at the floor.
void certify_arc_recursive(const Stock2& stock, double cx, double cy, double radius,
                           double a0, double a1, double r, double cap,
                           double cap_chord_ratio, CertifiedTea& acc, int depth)
{
    const double half_arc = 0.5 * radius * std::abs(a1 - a0);
    const double rho = stock.min_subtraction_radius();
    const double guard = TEA_GUARD_SAFETY_FACTOR * tea_growth_bound_curved(half_arc, r, rho);
    const double cap_guarded = cap - guard;
    acc.stations += 1;

    const double x0 = cx + radius * std::cos(a0), y0 = cy + radius * std::sin(a0);
    const double x1 = cx + radius * std::cos(a1), y1 = cy + radius * std::sin(a1);

    if (cap_guarded > 0.0) {
        const double sg = std::sin(0.5 * cap_guarded);
        const double guarded_ratio = 4.0 * sg * sg;
        const double gg = std::sin(0.5 * std::min(guard, std::numbers::pi));
        const double gap_close_ratio = 4.0 * gg * gg;
        const EngagementSample e0 = engagement_at(stock, x0, y0, r, guarded_ratio, gap_close_ratio);
        const EngagementSample e1 = engagement_at(stock, x1, y1, r, guarded_ratio, gap_close_ratio);
        acc.max_tea = std::max({acc.max_tea, e0.max_run_tea, e1.max_run_tea});
        if (!e0.cap_exceeded && !e1.cap_exceeded) return;
    } else {
        const EngagementSample e0 = engagement_at(stock, x0, y0, r, cap_chord_ratio);
        const EngagementSample e1 = engagement_at(stock, x1, y1, r, cap_chord_ratio);
        acc.max_tea = std::max({acc.max_tea, e0.max_run_tea, e1.max_run_tea});
    }

    if (2.0 * half_arc < STATION_FLOOR_FRACTION * r || depth >= CERTIFY_MAX_DEPTH) {
        acc.cap_certified = false;
        return;
    }
    const double am = 0.5 * (a0 + a1);
    certify_arc_recursive(stock, cx, cy, radius, a0, am, r, cap, cap_chord_ratio, acc, depth + 1);
    if (!acc.cap_certified) return;
    certify_arc_recursive(stock, cx, cy, radius, am, a1, r, cap, cap_chord_ratio, acc, depth + 1);
}
```

and the public entry point beside `certify_segment_tea`:

```cpp
CertifiedTea certify_arc_tea(const Stock2& stock, double cx, double cy, double radius,
                             double start_angle, double sweep, bool cw,
                             double tool_radius, double cap_radians)
{
    require_finite(cx, "cx");
    require_finite(cy, "cy");
    require_finite(radius, "radius");
    require_finite(start_angle, "start_angle");
    require_positive_tool_radius(tool_radius);
    if (!(radius > 0.0)) throw std::invalid_argument("radius must be strictly positive.");
    if (!(sweep > 0.0 && sweep <= 2.0 * std::numbers::pi))
        throw std::invalid_argument("sweep must be in (0, 2*pi].");
    if (!(cap_radians > 0.0 && cap_radians <= std::numbers::pi))
        throw std::invalid_argument("cap_radians must be in (0, pi].");

    const double sc = std::sin(0.5 * cap_radians);
    const double cap_chord_ratio = 4.0 * sc * sc;
    const double a1 = cw ? start_angle - sweep : start_angle + sweep;

    CertifiedTea acc{0.0, true, 0};
    certify_arc_recursive(stock, cx, cy, radius, start_angle, a1, tool_radius,
                          cap_radians, cap_chord_ratio, acc, 0);
    return acc;
}
```

Bind it in `register_engagement`:

```cpp
    m.def("certify_arc_tea",
          [](const Stock2& stock, double cx, double cy, double radius, double start_angle,
             double sweep, bool cw, double tool_radius, double cap_radians) {
              CertifiedTea c = certify_arc_tea(stock, cx, cy, radius, start_angle, sweep,
                                               cw, tool_radius, cap_radians);
              return std::make_tuple(c.max_tea, c.cap_certified, c.stations);
          },
          "stock"_a, "cx"_a, "cy"_a, "radius"_a, "start_angle"_a, "sweep"_a, "cw"_a,
          "tool_radius"_a, "cap_radians"_a);
```

- [ ] **Step 4: Run the C++-level tests**

Run: `pixi run pytest tests/test_engagement_audit.py -k "arc_certifier_certifies_ordinary or over_cap_circle" -v`
Expected: PASS, 2 tests.

- [ ] **Step 5: Swap the Python mirror onto the binding**

Replace the body of `_certify_arc_engagement` in `src/compas_cgal/engagement.py` (keep the
signature and the return tuple exactly as they are, so `_replay_arc` is untouched):

```python
def _certify_arc_engagement(stock: Stock, geometry: Arc | Circle, tool_radius: float, tea_cap: float) -> tuple[float, bool, int]:
    """Certify the engagement cap along one circular cut motion.

    Delegates to the compiled adaptive certifier (`certify_arc_tea`,
    `src/engagement_2.cpp`), which bisects the arc parameter under the same guarded
    exact station test `certify_segment_tea` uses for linear motions. The previous
    Python implementation sampled at a FIXED density and could not certify a
    trochoid circle wider than ~0.36 tool radii at a 90 degree cap, so it reported
    the generator's own circles uncertified regardless of their true engagement.

    Args:
        stock: The current (frozen for this measurement) stock.
        geometry: The `Arc` or `Circle` traced by the tool center.
        tool_radius: Tool radius.
        tea_cap: Engagement-angle cap in ``(0, pi]``.

    Returns:
        ``(max_tea, cap_certified, stations)``.
    """
    center = geometry.frame.point
    start = geometry.point_at(0.0)
    start_angle = math.atan2(float(start[1]) - float(center[1]), float(start[0]) - float(center[0]))
    sweep = 2.0 * math.pi if isinstance(geometry, Circle) else abs(float(geometry.angle))
    max_tea, cap_certified, stations = _stock_2.certify_arc_tea(
        stock.raw,
        float(center[0]),
        float(center[1]),
        float(geometry.radius),
        start_angle,
        sweep,
        bool(getattr(geometry, "_clockwise", False)),
        tool_radius,
        tea_cap,
    )
    return float(max_tea), bool(cap_certified), int(stations)
```

The `cw` flag belongs to the operation, not the geometry, so change `_replay_arc`
(`src/compas_cgal/engagement.py:402-411`) to pass it through explicitly rather than reading a
private attribute — add `clockwise: bool` as the last parameter of `_certify_arc_engagement`,
forward `op.clockwise` from `_replay_arc`, and replace the `getattr` above with that parameter.
Then delete `AUDIT_ARC_STEP_FRACTION` (`:36-40`) — it has no remaining reader.

- [ ] **Step 6: Run the end-to-end test and the three legacy arc tests**

Run: `pixi run pytest tests/test_engagement_audit.py -v`
Expected: PASS including `test_audit_certifies_a_nonzero_share_of_circular_operations`. The three
pre-existing arc tests (`refuses_circular_merge_over_cap`, `preserves_benign_single_run_immersion`,
`certifies_non_engaged_arc`) must still pass; adaptive refinement is strictly more capable than the
fixed density, so a previously-refused benign case may now certify — if
`refuses_circular_merge_over_cap` flips to `True`, that is a **real** loss of the gap-closure wiring
guard: stop and report rather than adjusting the test.

- [ ] **Step 7: Full suite**

Run: `pixi run baseline`

- [ ] **Step 8: Commit** (stage, show, ask)

```bash
git add src/engagement_2.h src/engagement_2.cpp src/compas_cgal/engagement.py tests/test_engagement_audit.py
git commit -m "feat(engagement): adaptive arc certification in C++, retire the fixed-density mirror"
```

---

### Task 11: three-valued engagement verdict

`_unmeasured()` returns `cap_certified=True` for retracts and plunges, so "certified" and "not
applicable" are the same value and `cap_violations` under-counts by design. A plunge is a
full-immersion bore — the move most likely to break a tool — and it is currently reported as
certified.

**Files:**
- Modify: `src/compas_cgal/engagement.py:60-115` (add the enum), `:321-323` (`_unmeasured`),
  `:390-411`, `:470-480` (aggregates)
- Test: `tests/test_engagement_audit.py`

**Interfaces:**
- Produces: `class EngagementVerdict(str, Enum)` with members `CERTIFIED = "certified"`,
  `VIOLATED = "violated"`, `NOT_APPLICABLE = "not_applicable"`. `OperationEngagement` gains
  `verdict: EngagementVerdict` and **keeps** `cap_certified: bool` as a derived property
  (`verdict is not VIOLATED`) so no existing consumer breaks. `EngagementReport` gains
  `unmeasured_ops: int`.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_engagement_audit.py`:

```python
def test_unmeasured_moves_are_not_reported_as_certified():
    """A plunge is a full-immersion bore, not a certified cut.

    Collapsing "certified" and "not applicable" into one boolean makes
    cap_violations under-count by construction and flatters plunges.
    """
    import math

    from compas.geometry import Polygon

    from compas_cgal.engagement import EngagementVerdict, audit_toolpath_engagement
    from compas_cgal.toolpath import OperationType, trochoidal_mat_toolpath_circular

    polygon = Polygon([(0, 0, 0), (12, 0, 0), (12, 8, 0), (0, 8, 0)])
    result = trochoidal_mat_toolpath_circular(polygon, tool_diameter=1.0, pitch=0.75, clearance_z=3.0)
    report = audit_toolpath_engagement(polygon, result, tool_diameter=1.0, tea_cap=math.pi / 2)

    plunges = [e for e in report.operations if e.operation == OperationType.PLUNGE]
    assert plunges, "the fixture must contain plunges"
    for e in plunges:
        assert e.verdict is EngagementVerdict.NOT_APPLICABLE
    assert report.unmeasured_ops == sum(1 for e in report.operations if e.verdict is EngagementVerdict.NOT_APPLICABLE)
    assert report.cap_violations == sum(1 for e in report.operations if e.verdict is EngagementVerdict.VIOLATED)
```

- [ ] **Step 2: Run to verify it fails**

Run: `pixi run pytest tests/test_engagement_audit.py -k unmeasured_moves -v`
Expected: FAIL — `ImportError: cannot import name 'EngagementVerdict'`.

- [ ] **Step 3: Add the verdict type**

In `src/compas_cgal/engagement.py`, after the exception classes:

```python
class EngagementVerdict(str, Enum):
    """Outcome of certifying one motion against the engagement cap.

    Three-valued because a boolean conflates two different statements: "the cap was
    checked and holds" and "the cap does not apply to this move". Collapsing them
    makes an unmeasured plunge indistinguishable from a certified cut.
    """

    CERTIFIED = "certified"
    VIOLATED = "violated"
    NOT_APPLICABLE = "not_applicable"
```

Add `from enum import Enum` to the imports. Extend `OperationEngagement`:

```python
    verdict: EngagementVerdict

    @property
    def cap_certified(self) -> bool:
        """``True`` unless the cap was checked and could not be established.

        Retained so existing consumers keep working; `verdict` carries the
        distinction between a certified cut and an unmeasured move.
        """
        return self.verdict is not EngagementVerdict.VIOLATED
```

and remove `cap_certified` from the field list (it becomes derived). Update `_unmeasured` to pass
`verdict=EngagementVerdict.NOT_APPLICABLE`, and both `OperationEngagement(...)` constructions in
`_replay_line` / `_replay_arc` to pass
`verdict=EngagementVerdict.CERTIFIED if cap_certified else EngagementVerdict.VIOLATED`. Add
`unmeasured_ops: int` to `EngagementReport` with a docstring line, and compute the aggregates:

```python
    cap_violations = sum(1 for e in operations if e.verdict is EngagementVerdict.VIOLATED)
    unmeasured_ops = sum(1 for e in operations if e.verdict is EngagementVerdict.NOT_APPLICABLE)
```

Export `EngagementVerdict` — Task 20 adds the `__all__` this joins.

- [ ] **Step 4: Run to verify it passes**

Run: `pixi run pytest tests/test_engagement_audit.py -k unmeasured_moves -v`
Expected: PASS.

- [ ] **Step 5: Full suite**

Run: `pixi run baseline`

- [ ] **Step 6: Commit** (stage, show, ask)

```bash
git add src/compas_cgal/engagement.py tests/test_engagement_audit.py
git commit -m "feat(engagement): three-valued verdict separates certified from not-applicable"
```

---

# Phase 4 — Platform, contracts, documentation

### Task 12: settle the Python floor (spec C3)

`pyproject.toml` declares `requires-python = ">=3.9"`, but `toolpath.py`, `stock.py`,
`engagement.py` and `isolines.py` use PEP 604 `X | Y` in evaluated annotations with no
`from __future__ import annotations`, so they fail at import on 3.9. CI builds 3.10 only.

**Files:**
- Modify: `src/compas_cgal/{toolpath,stock,engagement,isolines}.py` (add the future import)
- Modify: `.github/workflows/build.yml:18`
- Test: `tests/test_environment.py`

**Interfaces:**
- Consumes: `tests/test_environment.py` from Task 0.
- Produces: nothing.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_environment.py`:

```python
import ast
import tomllib
from pathlib import Path

import pytest

_PACKAGE = Path(__file__).resolve().parents[1] / "src" / "compas_cgal"


def _declared_python_floor() -> tuple[int, int]:
    manifest = tomllib.loads((Path(__file__).resolve().parents[1] / "pyproject.toml").read_text())
    spec = manifest["project"]["requires-python"]
    major, minor = spec.lstrip(">=~^ ").split(".")[:2]
    return int(major), int(minor)


@pytest.mark.parametrize("module", sorted(p.name for p in _PACKAGE.glob("*.py")))
def test_pep604_annotations_are_postponed_below_python_310(module):
    """`X | Y` in an evaluated annotation raises TypeError below 3.10.

    Either the declared floor is >= 3.10, or every module using the syntax must
    postpone evaluation with `from __future__ import annotations`.
    """
    if _declared_python_floor() >= (3, 10):
        return  # the floor itself permits the syntax; nothing to postpone
    source = (_PACKAGE / module).read_text()
    tree = ast.parse(source)
    uses_pep604 = any(isinstance(n, ast.BinOp) and isinstance(n.op, ast.BitOr) for n in ast.walk(tree))
    if not uses_pep604:
        return
    postponed = any(
        isinstance(n, ast.ImportFrom) and n.module == "__future__" and any(a.name == "annotations" for a in n.names)
        for n in tree.body
    )
    assert postponed, f"{module} uses PEP 604 unions but the declared floor is < 3.10"
```

Note the early `return` rather than `pytest.skip`: CLAUDE.md forbids skips, and once the declared
floor reaches 3.10 the check is genuinely satisfied rather than unrunnable.

- [ ] **Step 2: Run to verify it fails**

Run: `pixi run pytest tests/test_environment.py -k pep604 -v`
Expected: FAIL for `engagement.py`, `isolines.py`, `stock.py`, `toolpath.py`.

- [ ] **Step 3: Postpone annotation evaluation**

Add as the first statement after the module docstring in each of
`src/compas_cgal/{toolpath,stock,engagement,isolines}.py`:

```python
from __future__ import annotations
```

- [ ] **Step 4: Widen the CI matrix**

In `.github/workflows/build.yml:18`:

```yaml
        python: ["3.9", "3.10", "3.12", "3.13"]
```

- [ ] **Step 5: Run to verify it passes**

Run: `pixi run pytest tests/test_environment.py -v`
Expected: PASS.

- [ ] **Step 6: Full suite**

Run: `pixi run baseline`

- [ ] **Step 7: Commit** (stage, show, ask)

```bash
git add src/compas_cgal/toolpath.py src/compas_cgal/stock.py src/compas_cgal/engagement.py \
        src/compas_cgal/isolines.py .github/workflows/build.yml tests/test_environment.py
git commit -m "fix: honour the declared 3.9 floor — postponed annotations, wider CI matrix"
```

---

### Task 13: `certify_bridges` must not exit on an unverified run (spec I8)

The loop at `src/toolpath.cpp:437-466` breaks when a round inserts nothing. If the **last** round
inserts a midpoint, the refined run is never re-tested, so the contract "every bridge tangent is
certified" does not hold on that exit path.

**Files:**
- Modify: `src/toolpath.cpp:437-466`
- Test: `tests/test_toolpath.py`

**Interfaces:**
- Consumes: `Boundary::segment_clear`, `external_tangent`, `dedup_stations` (existing).
- Produces: no signature change.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_toolpath.py`:

```python
def test_every_emitted_bridge_is_clear_of_the_wall():
    """The bridge certification loop must not exit on a run it never re-tested.

    Directly checks the property the loop claims: every emitted cut LINE keeps at
    least the tool radius from the boundary, on the tightest pocket in the fixture
    set and at a pitch that forces repeated midpoint refinement.
    """
    tool_radius = 0.5
    polygon = _dumbbell(1.2)
    poly_xy = [list(pt[:2]) for pt in polygon.points]
    result = trochoidal_mat_toolpath_circular(polygon, tool_diameter=1.0, pitch=0.2, clearance_z=3.0)
    bridges = [op for op in result.operations if op.operation == "cut" and isinstance(op.geometry, Line)]
    assert bridges, "the fixture must emit bridge lines"
    for op in bridges:
        for t in [i / 32 for i in range(33)]:
            p = op.geometry.point_at(t)
            d = _distance_to_polygon_boundary_xy([float(p[0]), float(p[1])], poly_xy)
            assert d + GOUGE_TOL >= tool_radius, f"bridge point {t:.3f} at {d:.5f} < R={tool_radius}"
```

- [ ] **Step 2: Run to verify it fails or passes**

Run: `pixi run pytest tests/test_toolpath.py -k emitted_bridge -v`
Expected: with `mat_scale <= 1` enforced (Task 1) this may already PASS — the property holds
structurally. That is fine: the test is the contract, and Step 3 makes the loop honour it
unconditionally rather than by luck.

- [ ] **Step 3: Guarantee a clean final pass**

Restructure the loop so it can only exit after a round that inserted nothing:

```cpp
    bool inserted = true;
    for (int round = 0; round < MAX_BRIDGE_REFINE_ROUNDS && inserted; ++round) {
        dedup_stations(run, model);
        if (run.size() < 2) return {std::move(run)};

        inserted = false;
        std::vector<Station> refined;
        refined.reserve(run.size() * 2);

        for (std::size_t i = 0; i + 1 < run.size(); ++i) {
            // ... body unchanged ...
        }
        refined.push_back(run.back());
        run.swap(refined);
    }
    // A round that still inserted on the LAST iteration leaves bridges that were never
    // re-tested. Split the run at every remaining uncertified gap rather than emitting
    // an unverified bridge: the corridor is then simply not machined by this pass.
    if (inserted) {
        dedup_stations(run, model);
        std::vector<Station> final_pass;
        final_pass.reserve(run.size() * 2);
        for (std::size_t i = 0; i + 1 < run.size(); ++i) {
            final_pass.push_back(run[i]);
            const Circle_2 ci(run[i].center, run[i].radius * run[i].radius);
            const Circle_2 cj(run[i + 1].center, run[i + 1].radius * run[i + 1].radius);
            Segment_2 tangent;
            if (!external_tangent(ci, cj, edge_direction, climb_milling, tangent)) continue;
            if (boundary.segment_clear(tangent, model.tool_radius)) continue;
            final_pass.push_back(Station{run[i].center, -1.0, -1.0});  // split marker
        }
        final_pass.push_back(run.back());
        run.swap(final_pass);
    }
```

(Task 16 replaces the `-1.0` sentinel; keep it here so this task stays a one-idea change.)

- [ ] **Step 4: Run to verify it passes**

Run: `pixi run pytest tests/test_toolpath.py -k emitted_bridge -v`
Expected: PASS.

- [ ] **Step 5: Full suite**

Run: `pixi run baseline`

- [ ] **Step 6: Commit** (stage, show, ask)

```bash
git add src/toolpath.cpp tests/test_toolpath.py
git commit -m "fix(toolpath): bridge certification cannot exit on an unverified run"
```

---

### Task 14: complete the hole precondition check (spec I9)

`assemble_domain` tests only that hole *vertices* lie inside the outer boundary, so a hole edge
crossing a concave region passes, and a hole nested inside another hole passes. Both violate
`create_interior_straight_skeleton_2`'s precondition, which is undefined behaviour.

**Files:**
- Modify: `src/toolpath.cpp:548-585` (`assemble_domain`)
- Test: `tests/test_toolpath.py`

**Interfaces:**
- Consumes: `data_to_polygon` (existing).
- Produces: no signature change.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_toolpath.py`:

```python
def test_hole_edge_crossing_a_concave_outer_region_is_rejected():
    """All hole vertices inside is not containment: an edge can still exit a concave outer."""
    outer = Polygon([(0, 0, 0), (10, 0, 0), (10, 10, 0), (6, 5, 0), (0, 10, 0)])  # reflex at (6,5)
    hole = Polygon([(4, 6, 0), (8, 6, 0), (8, 8, 0), (4, 8, 0)])  # vertices inside, edge crosses the notch
    with pytest.raises(ValueError, match="inside|intersect"):
        trochoidal_mat_toolpath(outer, tool_diameter=0.5, holes=[hole])


def test_nested_hole_is_rejected():
    """A hole inside another hole is not a pocket island; the skeleton precondition forbids it."""
    outer = Polygon([(0, 0, 0), (20, 0, 0), (20, 20, 0), (0, 20, 0)])
    big = Polygon([(4, 4, 0), (16, 4, 0), (16, 16, 0), (4, 16, 0)])
    small = Polygon([(8, 8, 0), (12, 8, 0), (12, 12, 0), (8, 12, 0)])
    with pytest.raises(ValueError, match="nested|disjoint"):
        trochoidal_mat_toolpath(outer, tool_diameter=0.5, holes=[big, small])
```

- [ ] **Step 2: Run to verify they fail**

Run: `pixi run pytest tests/test_toolpath.py -k "concave_outer or nested_hole" -v`
Expected: two FAILs — either `DID NOT RAISE` or a CGAL abort.

- [ ] **Step 3: Check edges, not just vertices**

In `assemble_domain`, replace the vertex-containment loop and extend the pairwise loop:

```cpp
    for (const auto& hole_data : holes) {
        Polygon_2 hole = data_to_polygon(hole_data);  // validates simplicity, makes CCW
        // Containment: every vertex inside AND no edge crossing the outer boundary.
        // Vertices alone are insufficient -- across a reflex outer vertex a hole edge
        // can leave and re-enter with both endpoints inside. Feeding that to
        // create_interior_straight_skeleton_2 is undefined behaviour.
        for (auto vertex_iter = hole.vertices_begin(); vertex_iter != hole.vertices_end(); ++vertex_iter) {
            if (outer.bounded_side(*vertex_iter) != CGAL::ON_BOUNDED_SIDE) {
                throw std::invalid_argument("Hole polygon must lie strictly inside the outer boundary.");
            }
        }
        for (auto eh = hole.edges_begin(); eh != hole.edges_end(); ++eh) {
            for (auto eo = outer.edges_begin(); eo != outer.edges_end(); ++eo) {
                if (CGAL::do_intersect(*eh, *eo)) {
                    throw std::invalid_argument("Hole polygon edge intersects the outer boundary; the hole must lie strictly inside.");
                }
            }
        }
        for (const auto& other : hole_polygons) {
            // Nesting: edge-disjointness alone accepts one hole wholly inside another.
            // Test one vertex of each against the other's interior -- with disjoint
            // edges, one vertex decides containment for the whole polygon.
            if (other.bounded_side(*hole.vertices_begin()) == CGAL::ON_BOUNDED_SIDE ||
                hole.bounded_side(*other.vertices_begin()) == CGAL::ON_BOUNDED_SIDE) {
                throw std::invalid_argument("Hole polygons must not be nested.");
            }
            for (auto ea = hole.edges_begin(); ea != hole.edges_end(); ++ea) {
                for (auto eb = other.edges_begin(); eb != other.edges_end(); ++eb) {
                    if (CGAL::do_intersect(*ea, *eb)) {
                        throw std::invalid_argument("Hole polygons must be pairwise disjoint.");
                    }
                }
            }
        }
        // ... rest unchanged ...
```

- [ ] **Step 4: Run to verify they pass**

Run: `pixi run pytest tests/test_toolpath.py -k "concave_outer or nested_hole" -v`
Expected: PASS, 2 tests.

- [ ] **Step 5: Full suite**

Run: `pixi run baseline`
Expected: `test_holes_island_respected` still passes — its hole is genuinely inside and disjoint.

- [ ] **Step 6: Commit** (stage, show, ask)

```bash
git add src/toolpath.cpp tests/test_toolpath.py
git commit -m "fix(toolpath): reject hole edges crossing the outer boundary and nested holes"
```

---

### Task 15: tests for `isolines` (spec I10)

169 lines of public API, in the nav and the API reference, with no test file.

**Files:**
- Create: `tests/test_isolines.py`
- Modify: `src/compas_cgal/isolines.py:151-153` (named errors)

**Interfaces:**
- Consumes: `compas_cgal.isolines.isolines(mesh, scalars, isovalues=None, n=None, ...)` — read the
  current signature before writing the tests; do not assume it.
- Produces: `class InvalidIsovalueSpecError(ValueError)` in `isolines.py`.

- [ ] **Step 1: Read the signature**

Run: `sed -n '88,169p' src/compas_cgal/isolines.py`
Write the tests against what you see, not against this plan's paraphrase.

- [ ] **Step 2: Write the failing tests**

Create `tests/test_isolines.py`:

```python
"""Coverage for `compas_cgal.isolines` — a public, documented module that shipped untested."""

from __future__ import annotations

import numpy as np
import pytest
from compas.datastructures import Mesh

from compas_cgal.isolines import InvalidIsovalueSpecError, isolines


def _unit_grid(n: int = 8) -> tuple[Mesh, list[float]]:
    """Flat n x n triangulated grid on [0, 1]^2 with the scalar field f(x, y) = x."""
    mesh = Mesh()
    keys = [[mesh.add_vertex(x=i / n, y=j / n, z=0.0) for j in range(n + 1)] for i in range(n + 1)]
    for i in range(n):
        for j in range(n):
            mesh.add_face([keys[i][j], keys[i + 1][j], keys[i + 1][j + 1]])
            mesh.add_face([keys[i][j], keys[i + 1][j + 1], keys[i][j + 1]])
    scalars = [mesh.vertex_attribute(v, "x") for v in mesh.vertices()]
    return mesh, scalars


def test_isoline_of_a_linear_field_is_a_straight_line_at_the_isovalue():
    """f(x, y) = x: the 0.5 isoline is x == 0.5 everywhere along it."""
    mesh, scalars = _unit_grid()
    polylines = isolines(mesh, scalars, isovalues=[0.5])
    assert len(polylines) > 0
    for pts in polylines:
        xs = np.asarray(pts, dtype=np.float64)[:, 0]
        assert np.allclose(xs, 0.5, atol=1e-9), f"isoline strays from x=0.5: {xs.min()}..{xs.max()}"


def test_requesting_n_isolines_returns_that_many_levels():
    mesh, scalars = _unit_grid()
    polylines = isolines(mesh, scalars, n=3)
    levels = {round(float(np.asarray(p)[0, 0]), 9) for p in polylines}
    assert len(levels) == 3


def test_neither_isovalues_nor_n_is_rejected():
    mesh, scalars = _unit_grid()
    with pytest.raises(InvalidIsovalueSpecError, match="isovalues or n"):
        isolines(mesh, scalars)


def test_both_isovalues_and_n_is_rejected():
    mesh, scalars = _unit_grid()
    with pytest.raises(InvalidIsovalueSpecError, match="not both"):
        isolines(mesh, scalars, isovalues=[0.5], n=3)


def test_isovalue_outside_the_field_range_yields_no_polylines():
    mesh, scalars = _unit_grid()
    assert isolines(mesh, scalars, isovalues=[5.0]) == []
```

- [ ] **Step 3: Run to verify they fail**

Run: `pixi run pytest tests/test_isolines.py -v`
Expected: `ImportError: cannot import name 'InvalidIsovalueSpecError'`.

- [ ] **Step 4: Add the named error**

In `src/compas_cgal/isolines.py`, after the imports:

```python
class InvalidIsovalueSpecError(ValueError):
    """The isovalue specification is missing or over-specified (`isovalues` xor `n`)."""
```

Replace the two bare raises at `:151-153` with `raise InvalidIsovalueSpecError(...)`, keeping the
existing message text so the `match=` patterns hold. Add it to `__all__`.

- [ ] **Step 5: Run to verify they pass**

Run: `pixi run pytest tests/test_isolines.py -v`
Expected: PASS, 5 tests. If a behavioural test fails, that is a genuine `isolines` finding — report
it before changing the test.

- [ ] **Step 6: Full suite**

Run: `pixi run baseline`

- [ ] **Step 7: Commit** (stage, show, ask)

```bash
git add tests/test_isolines.py src/compas_cgal/isolines.py
git commit -m "test(isolines): first coverage for the public module; named isovalue error"
```

---

### Task 16: typed run boundary instead of the `Station{-1, -1}` sentinel

`certify_bridges` encodes "split here" as `Station{center, -1.0, -1.0}` and the partition loop reads
it back as `clearance < 0.0`. A `Station` with negative clearance is not a station; the invalid
state is representable and one missed check emits it as geometry.

**Files:**
- Modify: `src/toolpath.cpp:437-481`
- Test: covered by the existing `tests/test_toolpath.py` suite (behaviour-preserving)

**Interfaces:**
- Consumes: `Station { Point_2 center; double clearance; double radius; }` (unchanged).
- Produces: `certify_bridges` returns `std::vector<std::vector<Station>>` exactly as before; the
  sentinel disappears from the intermediate representation only.

- [ ] **Step 1: Record the behaviour to preserve**

Run: `pixi run pytest tests/test_toolpath.py -q`
Expected: all green. Note the count — this refactor must not change it.

- [ ] **Step 2: Replace the sentinel with an explicit split index set**

Change the refinement loop to carry `std::vector<Station> refined` plus
`std::vector<std::size_t> split_after` (indices into `refined` after which the run breaks), and
partition on that:

```cpp
    // Split points are carried as INDICES, not as an in-band Station with negative
    // clearance. A Station whose clearance is -1 is not a station; keeping the
    // invalid state unrepresentable means a missed check cannot emit it as geometry.
    std::vector<std::size_t> split_after;
    // ... inside the per-pair loop, where the sentinel was pushed:
    split_after.push_back(refined.size() - 1);
    // ... after the loop:
    std::vector<std::vector<Station>> parts;
    std::vector<Station> current;
    std::size_t next_split = 0;
    for (std::size_t i = 0; i < run.size(); ++i) {
        current.push_back(run[i]);
        if (next_split < split_after.size() && split_after[next_split] == i) {
            ++next_split;
            if (current.size() >= 2) parts.push_back(std::move(current));
            current.clear();
        }
    }
    if (current.size() >= 2) parts.push_back(std::move(current));
    return parts;
```

Carry `split_after` across refinement rounds by rebuilding it each round (it is derived from the
current `run`, so a stale copy is never used).

- [ ] **Step 3: Run to verify behaviour is unchanged**

Run: `pixi run pytest tests/test_toolpath.py -q`
Expected: identical pass count, identical output. Confirm the geometry is byte-identical:

```bash
pixi run pytest tests/test_toolpath.py -k "no_gouge or continuity or stepover" -v
```

- [ ] **Step 4: Full suite**

Run: `pixi run baseline`

- [ ] **Step 5: Commit** (stage, show, ask)

```bash
git add src/toolpath.cpp
git commit -m "refactor(toolpath): split runs by index, drop the negative-clearance sentinel"
```

---

### Task 17: `external_tangent` in kernel primitives

`external_tangent` (`src/toolpath.cpp:322-360`) decomposes `Point_2`/`Vector_2` into raw doubles and
hand-rolls the perpendicular as `vx = -uy, vy = ux` — the exact pattern `CLAUDE.md` lists under
"anti-patterns to REJECT in review". The `sqrt` is genuinely inexact and stays; the rotation and
point construction do not have to be.

**Files:**
- Modify: `src/toolpath.cpp:322-360`
- Test: covered by the existing suite (behaviour-preserving)

**Interfaces:**
- Consumes: `Vector_2::perpendicular(CGAL::Orientation)`, `CGAL::orientation` (existing CGAL).
- Produces: no signature change.

- [ ] **Step 1: Record the behaviour to preserve**

Run: `pixi run pytest tests/test_toolpath.py -q` — note the count.

- [ ] **Step 2: Rewrite the construction in kernel types**

```cpp
    const Vector_2 d = c1.center() - c0.center();
    const double length = approx_length(d);
    const double r0 = approx_radius(c0);
    const double r1 = approx_radius(c1);
    const double delta = r1 - r0;
    if (!(std::abs(delta) < length)) {
        throw std::logic_error("external_tangent: contained circles reached tangent construction");
    }

    // The external-tangent unit normal is m*u + h*v with u the unit centre direction,
    // v its CCW perpendicular, m = -delta/length and h = sqrt(1 - m^2). The sqrt is
    // inherently inexact (construction guard, not a decision), but the direction
    // algebra stays in kernel types: Vector_2::perpendicular for v, Vector_2 scaling
    // and Point_2 + Vector_2 for the tangent endpoints. No coordinate is unpacked.
    const Vector_2 u = d / length;
    const Vector_2 v = u.perpendicular(CGAL::COUNTERCLOCKWISE);
    const double m = -delta / length;
    const double h = std::sqrt(std::max(0.0, 1.0 - m * m));
    const Vector_2 n1 = m * u + h * v;
    const Vector_2 n2 = m * u - h * v;

    const Segment_2 ta(c0.center() + r0 * n1, c1.center() + r1 * n1);
    const Segment_2 tb(c0.center() + r0 * n2, c1.center() + r1 * n2);

    const Point_2& ci = c0.center();
    const auto orient_a = CGAL::orientation(ci, ci + edge_direction, ta.source());
    tangent = climb_milling
        ? (orient_a == CGAL::LEFT_TURN ? ta : tb)
        : (orient_a == CGAL::RIGHT_TURN ? ta : tb);
    return true;
```

- [ ] **Step 3: Run to verify behaviour is unchanged**

Run: `pixi run pytest tests/test_toolpath.py -v`
Expected: identical pass count. `test_tangent_vectors_match_geometry`, `test_consistent_winding`
and `test_climb_parameter_controls_winding` are the discriminating tests here.

- [ ] **Step 4: Full suite**

Run: `pixi run baseline`

- [ ] **Step 5: Commit** (stage, show, ask)

```bash
git add src/toolpath.cpp
git commit -m "refactor(toolpath): external tangent via kernel primitives, no coordinate unpacking"
```

---

### Task 18: one polygon validator

`_polygon_to_ccw_vertices` is duplicated verbatim in `src/compas_cgal/toolpath.py:107` and
`src/compas_cgal/stock.py:19`. The stated reason ("the toolpath helper is private") is undercut by
`stock.py` importing `InvalidPolygonError` from `toolpath.py` anyway — which also makes the
low-level exact stock model depend on the high-level generator.

**Files:**
- Create: `src/compas_cgal/_validation.py`
- Modify: `src/compas_cgal/toolpath.py:58-60, 107-134`, `src/compas_cgal/stock.py:11-66`
- Test: `tests/test_stock.py`

**Interfaces:**
- Produces: `compas_cgal._validation.InvalidPolygonError` and
  `compas_cgal._validation.polygon_to_ccw_vertices(polygon: Polygon) -> np.ndarray`.
  `compas_cgal.toolpath.InvalidPolygonError` remains importable as a re-export, so no consumer
  breaks.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_stock.py`:

```python
def test_stock_does_not_depend_on_the_generator_module():
    """The exact stock model must not import the high-level toolpath generator.

    Layering: `stock` is the primitive `toolpath` builds on. The shared polygon
    validator lives in `_validation`, which neither owns.
    """
    import ast
    from pathlib import Path

    source = (Path(__file__).resolve().parents[1] / "src" / "compas_cgal" / "stock.py").read_text()
    imported = {
        n.module
        for n in ast.walk(ast.parse(source))
        if isinstance(n, ast.ImportFrom) and n.module
    }
    assert "compas_cgal.toolpath" not in imported


def test_both_facades_share_one_polygon_validator():
    from compas_cgal import stock as stock_module
    from compas_cgal import toolpath as toolpath_module
    from compas_cgal._validation import polygon_to_ccw_vertices

    assert stock_module._polygon_to_ccw_vertices is polygon_to_ccw_vertices
    assert toolpath_module._polygon_to_ccw_vertices is polygon_to_ccw_vertices
```

- [ ] **Step 2: Run to verify it fails**

Run: `pixi run pytest tests/test_stock.py -k "does_not_depend or share_one" -v`
Expected: FAIL — `stock.py` imports `compas_cgal.toolpath`; `_validation` does not exist.

- [ ] **Step 3: Create the shared module**

Create `src/compas_cgal/_validation.py` containing `InvalidPolygonError` and
`polygon_to_ccw_vertices` — move the body verbatim from `src/compas_cgal/toolpath.py:107-134`
(it and the `stock.py` copy are byte-identical apart from the docstring; keep the `toolpath.py`
docstring and add the `Args/Returns/Raises` sections from the `stock.py` one).

- [ ] **Step 4: Point both façades at it**

In `src/compas_cgal/toolpath.py`, replace the class definition and the helper with:

```python
from compas_cgal._validation import InvalidPolygonError
from compas_cgal._validation import polygon_to_ccw_vertices as _polygon_to_ccw_vertices
```

keeping `"InvalidPolygonError"` in `__all__` so the public name is unchanged. In
`src/compas_cgal/stock.py`, delete the duplicated function and the `toolpath` import, and use the
same two imports.

- [ ] **Step 5: Run to verify it passes**

Run: `pixi run pytest tests/test_stock.py -k "does_not_depend or share_one" -v`
Expected: PASS, 2 tests.

- [ ] **Step 6: Full suite**

Run: `pixi run baseline`

- [ ] **Step 7: Commit** (stage, show, ask)

```bash
git add src/compas_cgal/_validation.py src/compas_cgal/toolpath.py src/compas_cgal/stock.py tests/test_stock.py
git commit -m "refactor: one polygon validator, stock no longer depends on toolpath"
```

---

### Task 19: say "exact" only where it means exact-kernel (spec I7)

`toolpath.{h,cpp,py}` describe Epick `squared_distance` results as "EXACT clearance", "exact early
rejection" and "Clearance queries are exact". In this repository "exact" has a precise meaning set
by `docs/exactness.md`, and `Boundary` does not meet it: `K::FT` is `double` and `squared_distance`
is an inexact *construction* whose rounded value then takes a gouge **decision**. The intended
claim — "measured at each station, never interpolated" — is true and worth stating.

**Files:**
- Modify: `src/toolpath.h:14-25, 32-38`; `src/toolpath.cpp:170-176, 637-644, 916-925`;
  `src/compas_cgal/toolpath.py:1-17, 243-266, 289-295`
- Modify: `docs/exactness.md` (one paragraph)
- Test: `tests/test_environment.py`

**Interfaces:** documentation only; no signature or behaviour change.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_environment.py`:

```python
def test_epick_modules_do_not_claim_exactness():
    """"exact" is reserved for the Epeck/circle-segment kernel (docs/exactness.md).

    `toolpath` runs on Epick, where squared_distance is an inexact construction, so
    describing its clearance queries as exact collides with the term the certifier
    modules use for a stronger property.
    """
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    offenders = []
    for path in (root / "src" / "toolpath.h", root / "src" / "toolpath.cpp", root / "src" / "compas_cgal" / "toolpath.py"):
        for lineno, line in enumerate(path.read_text().splitlines(), 1):
            lowered = line.lower()
            if "exact" in lowered and "epeck" not in lowered and "not exact" not in lowered:
                offenders.append(f"{path.name}:{lineno}: {line.strip()}")
    assert not offenders, "Epick module claims exactness:\n" + "\n".join(offenders)
```

- [ ] **Step 2: Run to verify it fails**

Run: `pixi run pytest tests/test_environment.py -k exactness -v`
Expected: FAIL listing each occurrence.

- [ ] **Step 3: Restate the claim precisely**

Replace every listed occurrence with the property that is actually true. Examples — apply the same
substitution to each site the test reports:

- `src/toolpath.h:14` → "Interior straight-skeleton vertices with per-vertex boundary clearance."
- `src/toolpath.h:36-37` → "every station radius is measured at its own center (never interpolated
  between skeleton vertices), and every bridge tangent is checked against the boundary at tool
  radius."
- `src/toolpath.cpp:176` → "Certified gouge-free test: every point of *segment* keeps at least
  *clearance* to the boundary, evaluated in Epick (inexact constructions, ~1e-16 relative)."
- `src/compas_cgal/toolpath.py:3-5` → "every trochoid radius is measured at its own station center
  rather than interpolated, and every bridge, lead, and flat link is checked against the boundary at
  tool radius, so emitted motions are gouge-free by construction."
- `src/compas_cgal/toolpath.py:15` → "Clearance queries run in Epick (inexact constructions) and are
  O(edges) per query without an acceleration structure, which is ample for pocket boundaries of
  hundreds of edges."

- [ ] **Step 4: Record the distinction in the doctrine page**

Append to `docs/exactness.md`, under the deciding/reporting section:

```markdown
### Which modules are exact

`stock_2` and `engagement_2` run on `Exact_predicates_exact_constructions_kernel` with
`Gps_circle_segment_traits_2`; their decisions are exact predicates on exact quantities and the
word *exact* applies in the strong sense this page defines.

`toolpath` runs on `Exact_predicates_inexact_constructions_kernel`. Its clearance queries use
`CGAL::squared_distance`, which is a **construction** — the returned `K::FT` is a rounded `double`,
and `Boundary::segment_clear` compares two such doubles to take a gouge decision. That decision is
therefore accurate to roughly 1e-16 relative, not exact. The generator's real guarantee is a
different and still-strong one: **clearance is measured at every station rather than interpolated
between skeleton vertices**, which is what removed the reflex-vertex gouge. Say that; do not say
"exact".
```

- [ ] **Step 5: Run to verify it passes**

Run: `pixi run pytest tests/test_environment.py -k exactness -v`
Expected: PASS.

- [ ] **Step 6: Full suite**

Run: `pixi run baseline`

- [ ] **Step 7: Commit** (stage, show, ask)

```bash
git add src/toolpath.h src/toolpath.cpp src/compas_cgal/toolpath.py docs/exactness.md tests/test_environment.py
git commit -m "docs: reserve \"exact\" for the exact-constructions kernel"
```

---

### Task 20: close the documentation stage (spec I11 + hygiene)

`compas_cgal.stock` and `compas_cgal.engagement` — 658 lines carrying the branch's thesis — have no
API page and are not in the nav, while `site/plans/` publishes internal planning and session-state
documents. `CLAUDE.md` makes stage-closing documentation mandatory.

**Files:**
- Create: `docs/api/compas_cgal.stock.md`, `docs/api/compas_cgal.engagement.md`
- Modify: `mkdocs.yml:138-192`
- Move: `docs/plans/` → `plans/`, `docs/superpowers/` → `superpowers/` (repo root, out of `docs_dir`)
- Modify: `pyproject.toml` (`addopts`), `src/compas_cgal/{stock,engagement}.py` (`__all__`),
  `src/compas_cgal/toolpath.py:87-105` (frozen dataclasses), `src/compas_cgal/toolpath.py:146`
- Modify: `CHANGELOG.md`
- Test: `tests/test_environment.py`

**Interfaces:**
- Consumes: `EngagementVerdict` (Task 11), `polygon_to_ccw_vertices` (Task 18).
- Produces: nothing importable.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_environment.py`:

```python
def test_every_public_module_has_an_api_page_in_the_nav():
    """CLAUDE.md: a stage is incomplete while its API is absent from the docs."""
    import re
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    public = {
        p.stem
        for p in (root / "src" / "compas_cgal").glob("*.py")
        if not p.stem.startswith("_")
    }
    nav = (root / "mkdocs.yml").read_text()
    missing = [m for m in sorted(public) if f"api/compas_cgal.{m}.md" not in nav]
    assert not missing, f"public modules absent from the docs nav: {missing}"
    missing_pages = [m for m in sorted(public) if not (root / "docs" / "api" / f"compas_cgal.{m}.md").exists()]
    assert not missing_pages, f"public modules without an API page: {missing_pages}"


def test_internal_planning_documents_are_not_published():
    """Plans and session-state files are repo artefacts, not documentation for users."""
    from pathlib import Path

    docs = Path(__file__).resolve().parents[1] / "docs"
    assert not (docs / "plans").exists(), "docs/plans is inside docs_dir and gets published"
    assert not (docs / "superpowers").exists(), "docs/superpowers is inside docs_dir and gets published"


def test_public_dataclasses_are_immutable():
    """Result records are values; two modules disagreeing on that is an inconsistency."""
    import dataclasses

    from compas_cgal.toolpath import ToolpathOperation, ToolpathResult

    for cls in (ToolpathOperation, ToolpathResult):
        assert dataclasses.fields(cls) is not None
        assert cls.__dataclass_params__.frozen, f"{cls.__name__} should be frozen"
```

- [ ] **Step 2: Run to verify they fail**

Run: `pixi run pytest tests/test_environment.py -k "api_page or planning_documents or immutable" -v`
Expected: three FAILs.

- [ ] **Step 3: Write the API pages**

`docs/api/compas_cgal.stock.md` — match the format of `docs/api/compas_cgal.toolpath.md` exactly
(read it first), for example:

```markdown
# compas_cgal.stock

::: compas_cgal.stock
```

`docs/api/compas_cgal.engagement.md` likewise for `compas_cgal.engagement`.

- [ ] **Step 4: Update the nav and stop publishing internal documents**

In `mkdocs.yml`, add to the API Reference block in alphabetical position:

```yaml
      - compas_cgal.engagement: api/compas_cgal.engagement.md
      - compas_cgal.stock: api/compas_cgal.stock.md
```

Move the internal trees out of `docs_dir` with history preserved:

```bash
git mv docs/plans plans
git mv docs/superpowers superpowers
```

Then update every intra-repo reference to those paths (`CLAUDE.md`, `CHANGELOG.md`, the plan and
spec files themselves):

```bash
grep -rln "docs/plans\|docs/superpowers" --include="*.md" --include="*.py" --include="*.yml" . | grep -v '^./site/'
```

- [ ] **Step 5: Hygiene sweep**

- `pyproject.toml`: delete `"--doctest-glob=*.rst"` from `addopts` — this repository bans rST.
- `src/compas_cgal/stock.py`: add
  `__all__ = ["Stock"]`.
- `src/compas_cgal/engagement.py`: add
  `__all__ = ["EngagementReport", "EngagementVerdict", "InvalidEngagementCapError",
  "InvalidToolDiameterError", "OperationEngagement", "UnexpectedToolpathGeometryError",
  "audit_toolpath_engagement"]`.
- `src/compas_cgal/toolpath.py:87,99`: change both `@dataclass` to `@dataclass(frozen=True)` to
  match `engagement.py`'s records.
- `src/compas_cgal/toolpath.py:146`: run `ruff format src/compas_cgal/toolpath.py` to collapse the
  implicit string concatenation.
- `CHANGELOG.md`: add a `### Fixed` block under `## Unreleased` listing, one line each, the
  behaviour changes from Tasks 1, 2, 4, 8, 10, 13 and 14. Do not list refactors.

- [ ] **Step 6: Run to verify they pass**

Run: `pixi run pytest tests/test_environment.py -v && pixi run --environment docs docs`
Expected: tests PASS; `mkdocs build --strict` succeeds with no warnings about files outside the nav.

- [ ] **Step 7: Full suite and lint**

Run: `pixi run baseline && pixi run lint && ruff format --check src/compas_cgal tests`
Expected: all green.

- [ ] **Step 8: Commit** (stage, show, ask)

```bash
git add -A
git commit -m "docs: API pages for stock and engagement, internal plans out of the site"
```

---

# Phase 5 — Post-rebase, on codex trunk

Every task in this phase assumes the rebase in `INTEGRATECODEX.md` has landed, so paths are the
**codex** layout. None of them touches `src/engagement_2.*` or `src/stock_2.*`.

---

### Task 24: an independent falsifier for circle certification

Codex's suite is disciplined — 1,022 test functions, zero `skip`/`xfail`, zero `pytest.approx` — with
one hole: segment certification has a real false-certificate hunter and **circle certification has
none**. The three tests asserting a certified circle verify only the implementation's own output,
including `_continuous_tea_2.verify_event_partition(trace.partition)`, which is the same module
verifying its own record.

This branch already built the missing instrument. Wire it up.

**Files:**
- Modify: `tests/adaptive/test_circle_oracle.py`
- Reference: `tests/adaptive/test_segment_oracle.py:97-140` (the pattern to mirror),
  `tests/test_false_certificate.py` (the rib constructions this branch carries)

**Interfaces:**
- Consumes: `_stock_2.engagement_at(stock, cx, cy, r, cap_chord_ratio, gap_close_ratio)` — the
  independent implementation; the rib and spiral-rib stock builders from
  `tests/test_false_certificate.py`.
- Produces: `_assert_no_confirmed_arc_violation(...)`, the circle analogue of
  `_assert_no_confirmed_dyadic_violation`.

- [ ] **Step 1: Read the segment falsifier and copy its shape, not its code**

`tests/adaptive/test_segment_oracle.py:97-140` walks `2**DYADIC_FALSIFIER_DEPTH` stations along a
certified motion and calls `_stock_2.engagement_at` — a *different* implementation from the one under
test — failing if any station's exact predicate fires. Mirror that over the arc parameter.

- [ ] **Step 2: Write the falsifier**

```python
def _assert_no_confirmed_arc_violation(stock, centre, radius, sweep, cw, tool_radius, cap_chord_ratio):
    """Fail if a certified circular motion has an exact violation at a dyadic station.

    Independence is the whole point: `engagement_at` shares no code with the event
    partition being checked, so agreement is evidence rather than tautology.
    """
    denominator = 2**DYADIC_FALSIFIER_DEPTH
    for numerator in range(denominator + 1):
        theta = start_angle + (-sweep if cw else sweep) * numerator / denominator
        cx = centre[0] + radius * math.cos(theta)
        cy = centre[1] + radius * math.sin(theta)
        exceeded = _stock_2.engagement_at(stock, cx, cy, tool_radius, cap_chord_ratio, 0.0)[2]
        if exceeded:
            pytest.fail(f"certified circle has an exact violation at {numerator}/{denominator}")
```

- [ ] **Step 3: Invoke it from every test that asserts a certified circle**

`test_circle_oracle.py:54`, `:83`, `:113`. Each currently asserts only
`trace.exact_verdict == "certified"` plus structural counts.

- [ ] **Step 4: Add the rib as a fixture, and pin liveness independently of the verdict**

Port the annular and spiral rib builders. Each rib test asserts two independent things: the exact
oracle finds a violating centre on the motion (liveness, established without the certifier), and the
certifier's verdict. A witness that can pass by its own construction degenerating is worthless — this
branch learned that the hard way; see `tests/test_false_certificate.py`'s liveness pins.

- [ ] **Step 5: Verify the falsifier can fail**

Inject a false certificate — force `exact_verdict = "certified"` on a rib motion — and confirm the
falsifier fires. A hunter never observed hunting is decoration.

- [ ] **Step 6: Run and commit**

```bash
pixi run baseline
git add tests/adaptive/test_circle_oracle.py
git commit -m "test(adaptive): independent falsifier for circle certification"
```

---

### Task 25: put the four quotable numbers behind the gate

`docs/engagement_controlled_toolpath.md` carries the project's strongest public claims. Four of them
have **zero occurrences** across `src/`, `tests/` and `benchmarks/` — they exist only in prose:
`136.4°` vs `119.5°`, `914.3`, `224 ms` vs `2.6 ms`, and `125.8°` radial saturation.

The *comparisons* are tested — `COMPARABLE_STEPOVER` derives the textbook stepover and asserts the
regulated path wins by ≥2× on exceedance — but the figures are unpinned, and the timing claim is
contradicted internally because no harness times the regulated generators at all.

**Files:**
- Modify: `benchmarks/quality.py` or `benchmarks/gate.py` (add the metrics)
- Modify: `benchmarks/runner.py` (time the regulated generators)
- Modify: `docs/engagement_controlled_toolpath.md`, `docs/radius_regulated_toolpath.md`

**Interfaces:**
- Consumes: `GATE_GENERATORS`, `MeasurementRecord`, the existing gate metric plumbing.
- Produces: gate metrics for peak measured engagement under a dialled stepover, the unregulated
  cutting length, and wall-clock for both engagement generators.

- [ ] **Step 1: Add the unregulated generator to `GATE_GENERATORS`**

`914.3` has no code behind it because the unregulated generator is not a gate generator, so nothing
recomputes it. Adding it also makes the 0.152% residue figure live.

- [ ] **Step 2: Add a stepover-comparison metric**

Emit the peak measured engagement for the textbook-derived stepover alongside the regulated peak, on
the same pocket, so `136.4` and `119.5` become computed values rather than recalled ones.

- [ ] **Step 3: Time the regulated generators**

`perf_counter` currently appears only in `benchmarks/runner.py` (timing the *unregulated* generator)
and `scripts/engagement_baseline.py` (also unregulated). Time both engagement generators on the gate
pockets and record the result.

- [ ] **Step 4: Replace the prose figures with the computed ones**

Every number in those two docs either cites a gate metric by name or is removed. Fix the stale
`72.4 ms` prior-art row while the file is open.

- [ ] **Step 5: Run and commit**

```bash
pixi run python -m benchmarks.cli --corpus all --out docs/benchmarks
pixi run baseline
git add benchmarks/ docs/engagement_controlled_toolpath.md docs/radius_regulated_toolpath.md
git commit -m "bench: gate the four quoted engagement numbers"
```

---

### Task 26: exercise holes on an engagement generator

`holes` is accepted by the engagement generators, passed through, and never exercised. The only holes
test in the repository is against the legacy `trochoidal_mat_toolpath_circular`. Islands change
medial-axis topology, so an untested path here is a real gap rather than a coverage statistic.

**Files:**
- Modify: `tests/test_engagement_toolpath.py`, `tests/test_engagement_radial_toolpath.py`

**Interfaces:** consumes the existing pocket fixtures and `engagement_controlled_toolpath`; produces
no new API.

- [ ] **Step 1: Add a pocket with one island** and assert the emitted programme keeps tool clearance
  from **both** boundaries, reusing the existing clearance assertion.
- [ ] **Step 2: Add a pocket with two disjoint islands**, which forces more than one guide component.
- [ ] **Step 3: Assert the cap still holds away from chain entries** — the property the no-holes
  fixtures already pin, on the harder topology.
- [ ] **Step 4: Run and commit**

```bash
pixi run baseline
git add tests/test_engagement_toolpath.py tests/test_engagement_radial_toolpath.py
git commit -m "test(engagement): exercise island pockets on both engagement generators"
```

---

### Task 27: one clause per figure caption

Panel labels read *"engagement-controlled, 120° cap"* and `docs/toolpath_drawing.md`'s caption is
*"Two trochoidal tool paths over the same pocket"*. A reader arriving from
`docs/continuous_engagement.md` may take "engagement-controlled" for the exact-certified pipeline. It
is not — the disclaimer currently lives only in a module docstring and one prior-art table row.

**Files:** `benchmarks/figures.py`, `benchmarks/qualityfigures.py`, `docs/toolpath_drawing.md`

- [ ] **Step 1: Add one clause to each caption** naming the guarantee's scope — per evaluated
  position, no certificate produced.
- [ ] **Step 2: Regenerate and confirm byte-stability**

```bash
pixi run python -m benchmarks.cli figures
git diff --stat docs/
```

Only caption text should move; `e0ae116` pinned the SVG hashsalt so redrawing identical paths rewrites
identical bytes, and any geometry change in that diff is a finding.

- [ ] **Step 3: Commit**

```bash
git add benchmarks/figures.py benchmarks/qualityfigures.py docs/toolpath_drawing.md
git commit -m "docs: figure captions state the guarantee's scope"
```

---

### Task 28: close the `adaptive/` producer–verifier divergence

`replay_certificate()` ends in an unconditional `raise` (`src/compas_cgal/adaptive/replay.py:1276`),
and independently its operation grammar rejects `RetraceSegmentOperation` (`replay.py:245-261`) — an
operation the generator emits by design (`generator.py:1448-1458`) and exercises in its own reference
fixture. Downstream, `coverage.require_complete()` (`replay.py:1164`) is unreachable, so the
guarantee *"a shipped programme leaves no reachable material"* has no test.

All 17 tests in `tests/adaptive/test_replay.py` are `pytest.raises`, and the helper discards the
return value — so the certificate of record is never observed being produced.

Ranked last in this phase because nothing ships through this path yet.

**Files:** `src/compas_cgal/adaptive/replay.py`, `tests/adaptive/test_replay.py`

- [ ] **Step 1: Decide the intent first, and record it.** Is `replay_certificate` unfinished, or
  deliberately disabled pending the retrace grammar? The unconditional `raise` reads as the latter.
  Do not remove it before answering — ask if the answer is not in the commit history.
- [ ] **Step 2: Admit `RetraceSegmentOperation` into the grammar**, or state in the docstring why the
  generator emits an operation the verifier will not accept.
- [ ] **Step 3: Write the first positive-path test** — bind the return value and assert on it, so a
  regression that made `replay_certificate` raise unconditionally would fail rather than pass.
- [ ] **Step 4: Reach `coverage.require_complete()`** with a fixture that gets past the terminality
  guard, giving `IncompletePocketCoverageError` its first pipeline-level exercise.
- [ ] **Step 5: Run and commit**

```bash
pixi run baseline
git add src/compas_cgal/adaptive/replay.py tests/adaptive/test_replay.py
git commit -m "fix(adaptive): replay verifier accepts the operations the generator emits"
```

---

### Task 29: a second cap, so the gate can attribute a defect again

`benchmarks/gate.py`'s own docstring states the design: *"Three pockets across two generators is the
smallest product in which a defect can be attributed: a finding on every pocket is the generator's, a
finding on every generator is the pocket's, and a finding on one cell is neither."*

**That attribution is currently void.** At `GATE_CAP_DEG = 120` the two gate generators emit
identical paths on all three pockets — cutting length `411.566` / `1018.665` / `565.311`, cut blocks
`49` / `77` / `90`, and the same failing criteria, to six figures. The wiring is correct
(`GATE_GENERATORS` holds two distinct functions); the radius ladder simply never fires at that cap,
which the earlier saturation sweep already showed — at cap 120 the advance-only and radius-regulated
columns agreed exactly, and diverged at 40, 60, 80 and 100. Six cells are measuring three paths.

The plumbing exists. `gate_pocket(name, tool_diameter=..., tea_cap_deg=...)` already takes the cap;
only the call site hard-codes it.

**Files:**
- Modify: `benchmarks/gate.py`, `tests/benchmarks/test_quality.py`

**Interfaces:**
- Consumes: `gate_pocket`, `GATE_POCKET_NAMES`, `GATE_GENERATOR_NAMES`.
- Produces: `GATE_CAP_DEGS: Tuple[float, ...]`, and a gate parametrised over pockets x generators x
  caps.

- [ ] **Step 1: Prove the divergence before relying on it.** Run both generators at each candidate
  cap and record cutting length and cut count. Pick the second cap from measurement, not from the
  old sweep — that sweep predates the rho work and the run-scoped metric.
- [ ] **Step 2: Add `GATE_CAP_DEGS`** with a comment giving each cap's job: 120 is the commissioned
  load a real roughing pass runs at, the second is the cap at which the ladder demonstrably engages.
- [ ] **Step 3: Parametrise the gate over the third axis** and confirm the cell count triples.
- [ ] **Step 4: Assert the axis is load-bearing.** A test that fails if the two generators produce
  identical output at *every* gate cap — otherwise this regresses silently the next time the ladder
  goes inert.
- [ ] **Step 5: Run and commit**

```bash
pixi run pytest tests/benchmarks -q -n auto
git add benchmarks/gate.py tests/benchmarks/test_quality.py
git commit -m "bench: a second cap, so the gate can attribute a defect to a generator"
```

---

### Task 30: register the three new generators in the gate

`rho_regulated_toolpath`, `chain_ordered_toolpath` and `spiral_entry_toolpath` ship
(`7af7a24`) with 19+ tests each, and **none is in `GATE_GENERATORS`** — so the machining-quality gate
never sees them. Their measured wins are recorded only in commit messages and module docstrings,
which is exactly the unpinned-number problem Task 25 exists to fix.

What they claim, on `rect_12x8` at cap 120: degenerate loops `4 -> 0` and tangent breaks `24 -> 0`
(rho regulation); material entries `5 -> 1` (chain ordering); worst chain entry `243.20` degrees
against a proven floor of `240` (spiral entry). Each is a gate criterion or one line from being one.

Ordered after Task 29 so the new cells land on the corrected axis rather than being re-parametrised
twice.

**Files:**
- Modify: `benchmarks/gate.py`, `tests/benchmarks/test_quality.py`

- [ ] **Step 1: Add the three to `GATE_GENERATORS` and `GATE_GENERATOR_NAMES`.** Expect the gate to
  go redder, not greener — the point is measurement, not a green board.
- [ ] **Step 2: Record the new failure profile per cell** before changing anything, so the trade each
  generator makes is documented rather than discovered later.
- [ ] **Step 3: Reconcile the uncut regression with Task 31.** rho regulation raises uncut fraction
  `0.0080 -> 0.0295` on `rect_12x8`; that is a real defect against a reachable-region denominator and
  must not be normalised as the new baseline.
- [ ] **Step 4: Run and commit**

```bash
pixi run pytest tests/benchmarks -q -n auto
git add benchmarks/gate.py tests/benchmarks/test_quality.py
git commit -m "bench: gate the rho-regulated, chain-ordered and spiral-entry generators"
```

---

### Task 31: consume `declined_regions` — the corner the roughing pass cannot reach

`RhoToolpathResult.declined_regions` and `OrderedToolpathResult.declined_regions` name the guide runs
below the degeneracy floor — on `rect_12x8`, 116 guide stations across 4 runs, whose gouge-free radius
is at or under the tool radius. The generators correctly refuse to bore them and correctly report
them. **Nothing downstream consumes the field.**

The consequence is measured: uncut fraction `0.0080 -> 0.0295` on `rect_12x8`, `0.0028 -> 0.0117` on
`rect_20x12`, `0.0069 -> 0.0352` on `L_shape`. Because `uncut_fraction`'s denominator is *what the
tool can reach* — pinned by
`test_the_uncut_denominator_is_what_the_tool_can_reach_and_not_the_pocket` — this is reachable
material left standing, not an artefact of the metric.

`docs/loop_radius_degeneracy.md` proves why no guide radius rescues it: `rho_max = (d - r) / 2`, so
`rho_max <= r` exactly when `d <= 3r`. A corner narrower than three tool radii cannot carry a
trochoid at any radius. The answer is a different operation, not a different radius.

**Files:**
- Create: a corner-clearing pass, beside the generators, never inside them
- Modify: `benchmarks/quality.py` (distinguish declined from silently-missed)

**Interfaces:**
- Consumes: `declined_regions`, the exact stock, `DEGENERATE_LOOP_RATIO`.
- Produces: motions clearing the declined runs, or a typed refusal naming the tool that could.

- [ ] **Step 1: Decide the contract first.** Either the roughing programme clears these corners with a
  non-trochoidal motion, or it declares them for a second operation with a smaller tool. Both are
  defensible; shipping neither is not.
- [ ] **Step 2: Split the metric.** `uncut_fraction` must separate material that was *declined and
  named* from material *silently missed*. A generator that reports what it left is categorically
  better than one that does not, and one number currently conflates them.
- [ ] **Step 3: Implement the chosen contract** as a new path; do not modify either generator.
- [ ] **Step 4: Show the uncut fraction returns to or below the `engagement_controlled` baseline**, or
  state precisely which part is structurally unreachable by this tool and why.
- [ ] **Step 5: Run and commit**

---

### Task 32: a green full suite, and the four adaptive reds

The full suite has not been observed green. It timed out at ten minutes in the session that landed
`d9d0e8a`, and the commits there went in on the targeted benchmark tests alone. Separately, four
tests fail on trunk, all in the candidate-family work and unrelated to anything above:

```
tests/adaptive/test_generator.py::test_task13f_full_continuation
tests/adaptive/test_generator.py::test_real_active_family_stops_at_unresolved_exact_event
tests/adaptive/test_route_retrace_generator.py::test_continuation_rejects_missing_retrace_commit
tests/adaptive/test_route_retrace_generator.py::test_route_retrace_derivation_rejects_unsupported_source_scope[terminal]
```

Sample divergence: expected `cap=0; gouge=56`, actual `cap=10; gouge=46`. Two of the four are
`DID NOT RAISE` — a guard that no longer fires, which is the more serious shape.

- [ ] **Step 1: Establish whether each red is a stale expectation or a live defect.** Do not touch a
  test until that question is answered per failure; the project rule is to fix the code, never the
  reference test.
- [ ] **Step 2: Fix, and record which of the two it was** in the commit message.
- [ ] **Step 3: Get one observed green full run** and record its wall-clock, so the CI timeout is set
  from measurement.
- [ ] **Step 4: Raise the benchmarks workflow timeout** to that measurement plus headroom, or split
  the suite so no single job exceeds it.

---

### Task 33: the physics KPI families still missing

`d9d0e8a` wired the orphaned `MaterialOutcome` / `MachineOutcome` / `ToolLifeOutcome` and added the
program-feasibility group, immersion stationarity, and load-at-design. Three families a machining
audience would expect are still absent, and all three are computable from the existing survey:

1. **Cutting force.** Engagement is a *proxy* for force. The mechanistic model gives the quantity
   itself — peak and RMS resultant force, its rate of change, and the deflection that force implies
   at the wall. Two paths with identical peak engagement can carry very different peak force.
   Requires one new coefficient (specific cutting pressure) in `MaterialModel`, so it is a
   *calibrated* outcome and belongs in `CalibratedOutcome`, never in `PathQuality`.
2. **Exit impact.** Climb versus conventional, and the count of exits taken at high engagement —
   which is what chips a carbide edge. Pure geometry plus engagement; **no calibration**, so it stays
   in `PathQuality`.
3. **Controller starvation.** The program group reports block-length distribution but not the
   starvation verdict, which needs a feed and an interpolation period. `MachineModel` already has the
   feed; add the interpolation period and emit the fraction of cut length in blocks the control
   cannot consume at commanded feed. A trochoidal path is the shape that provokes this — an early
   radius sweep produced 862 cuts at one cap.

Keep the calibration boundary that `benchmarks/models.py` establishes: a metric needing a material or
a machine takes it as an argument and raises a **named** error when absent, rather than substituting
a coefficient.

- [ ] **Step 1: Exit impact first** — it needs no new calibration and no new model field.
- [ ] **Step 2: Add specific cutting pressure to `MaterialModel`** with a cited default and a range
  check, then the force metrics into `CalibratedOutcome`.
- [ ] **Step 3: Add the interpolation period to `MachineModel`** and the starvation fraction.
- [ ] **Step 4: One anchored test per family**, in the style of
  `test_the_kernel_anchors_the_rim_arc_to_radial_depth_relation` — assert against a physical anchor,
  not against the implementation's own output.

---

### Task 34: generalise the native segment-site MAT

`SegmentSiteMedialAxis.build` accepts **exactly one polygon**. The chain
`src/segment_site_mat_bundle.cpp:39` -> `segment_site_catalog_sampling.cpp:283` ->
`segment_site_catalog_neck.cpp:230` -> `segment_site_catalog_graph.cpp:1220` reaches
`require_l_shape_input`, which demands a hole-free six-vertex ring whose vertices are exactly
`(0,0) (6,0) (6,2) (2,2) (2,6) (0,6)`. Measured: every other polygon raises
`UnsupportedCanonicalMatLShapeGraphError`, and at `max_refinement_depth=6` with `r=0.25` even the
fixture raises `ConicSamplingLimitError`. It is one polygon inside a parameter window.

A `canonical_rectangle_mat_graph` / `require_rectangle_input` pair exists
(`segment_site_catalog_graph.cpp:125,905`) but is reached only from `tests/native/task9_*.cpp` — no
Python path calls it, which is itself the evidence that shape-by-shape does not scale.

**This blocks the middle-curve construction entirely** and nothing else can use the true MAT until it
lands.

!!! warning "Do not justify this with the middle curve's smoothness"

    Measured head to head on the fixture, the middle curve is **not** gentler than the straight
    skeleton: guide-radius CoV `0.5719` against `0.5752`, statistically indistinguishable, because
    `rho = (d - r) / 2` is affine in the same clearance and inherits the same relative swing.
    Degenerate loops are *worse* — 24.5% against 13.0% — for the reason
    `docs/loop_radius_degeneracy.md` derives: the middle curve degenerates at `d <= 3r` where the
    MAT-centred construction degenerates at `d <= 2r + c`, a strictly wider band. The honest case for
    this task is **structural** — two wall-hugging strands per axis edge, parabolic arcs at reflex
    vertices, the tool grazing the wall by construction — and it is *untested*, because the only
    polygon that builds has no wide channel to show it on.

- [ ] **Step 1: Route `SegmentSiteMatBundle2::build` through a real segment-site Voronoi
  construction** rather than a canonical fixture graph.
- [ ] **Step 2: Keep the fixture graph as a regression oracle** — the new construction must reproduce
  it exactly on the L, or the difference must be explained.
- [ ] **Step 3: Fix the refinement floor** so the fixture builds at the depths its own tests use.
- [ ] **Step 4: Only then** revisit whether a middle-curve generator earns its place, on a pocket wide
  enough to test the structural claim.

---

### Task 35: reproduce Held's published examples

The project's algorithmic claims trace to Held & Pfeiffer 2025, and Figure 6 is reproduced
(`docs/assets/images/fig6_*.svg`). **The rest of the paper's examples are not.** Reproducing all of
them, with a per-example KPI delta against the published figures, is the strongest external
validation available and needs no new machinery beyond what the corpus already has.

- [ ] **Step 1: Locate or digitise the contours.** They were reported as already extracted from the
  paper and in the repository; an exhaustive search of all local branches, both checkouts including
  ignored files, and the Zotero storage folder beside the PDF did not find them. **Resolve this
  before anything else** — if they are not extracted, digitise from the figures and label the result
  as digitised, as was done for Figure 6.
- [ ] **Step 2: Add each as a corpus family** under `benchmarks/families/`, with provenance in the
  docstring: paper figure number, and measured-versus-digitised.
- [ ] **Step 3: Generate ours beside the published figure** using the existing `draw_comparison`.
- [ ] **Step 4: Tabulate the KPI delta per example** and, where we differ, say whether it is a
  difference in the algorithm, in the parameters, or a defect.
- [ ] **Step 5: Publish to mkdocs** with the deltas, and state plainly which examples we do worse on.

---

# Phase 6 — Held parity

The two items `INTEGRATECODEX.md` ranks first and second on the critical path. Both are post-rebase,
neither touches the conflict zone, and **35 should not run before them** — reproducing Held's contours
at a cap that does not bind would look like validation and would not be.

---

### Task 36: helical or ramped entry

Entry into virgin stock is currently a full slot at 360°: the tool plunges, then cuts a complete
annulus with no cleared side. It is warned about in both module docstrings and never hidden, which is
honest — but it is also the single largest source of failing gate criteria.

Measured on codex trunk: **every** cap exceedance on all three gate pockets is a chain-entry cut,
verified by set equality (`test_only_chain_entry_loops_are_measured_above_the_cap`), as is the 329°
engagement step. Away from entries the worst motion measures 119.19° against a 120° cap. So the cap
already holds everywhere the entry does not reach, and this task clears **two of the six** failing
gate criteria on its own.

**Files:**
- Modify: `src/compas_cgal/engagement_toolpath.py`, `src/compas_cgal/engagement_radial_toolpath.py`
- Modify: `tests/test_engagement_toolpath.py`, `tests/test_engagement_radial_toolpath.py`
- Modify: `benchmarks/gate.py` (the entry-exemption in the cap criterion, once it is no longer needed)

**Interfaces:**
- Consumes: the operation stream's per-operation `z_start` / `z_end`, which already carries plunge and
  retract as differing-z lines.
- Produces: no signature change to `engagement_controlled_toolpath`. The entry becomes a sequence of
  operations rather than a single plunge, so `ToolpathResult.operations` grows a new shape that
  downstream consumers — the audit's `_replay_line`, the benchmark metrics, the figure drawing — must
  all accept.

- [ ] **Step 1: Choose the entry strategy and record why**

Three candidates, and I do not have enough information to pick for you — read the operation stream
first, then decide and justify in the report:

| strategy | what it needs | cost |
| --- | --- | --- |
| **Helical** — descend on a helix inside the first loop | Z motion during an arc; the op stream carries per-op z, so this is a ramped arc sequence | changes the cut-plane assumption the engagement audit makes (`_infer_cut_height` takes the minimum motion z) |
| **Ramped** — descend along the first bridge segment | same z machinery, simpler geometry | a ramp is a 3D cutting move the audit explicitly refuses (`UnexpectedToolpathGeometryError`, "ramped 3D cutting is outside the engagement audit's cut-plane model") |
| **Pre-bored / entry from cleared stock** — start the chain where material is already gone | no Z machinery at all; needs the guide to admit an entry point adjacent to cleared region | pushes the problem to the caller unless the generator can find such a point itself |

Whichever you choose, **the audit must still be able to replay the result.** Both of the first two
break `_infer_cut_height`'s single-plane assumption, and the audit currently raises on a differing-z
line with nonzero XY travel. Say in your report what the audit needs in order to accept your entry.

- [ ] **Step 2: Write the failing test first**

The property is not "the entry looks helical" — it is that **no motion exceeds the cap**:

```python
def test_no_motion_exceeds_the_cap_including_the_entry():
    """The entry cut is the last place the cap is violated; after helical/ramped entry, nothing is.

    Before this task every cap exceedance on every gate pocket was a chain-entry cut,
    verified by set equality. This asserts the set is now empty rather than merely
    equal to the entry set.
    """
    result = engagement_controlled_toolpath(POCKET_12X8, tool_diameter=2.0, tea_cap_deg=120.0)
    over = _measured_over_cap(result, tea_cap_deg=120.0)
    assert over == set(), f"motions above the cap: {sorted(over)}"
```

Reuse whatever `test_only_chain_entry_loops_are_measured_above_the_cap` uses to build its measured
set, so the two tests are comparable and the older one can be rewritten to assert the stronger
property rather than deleted.

- [ ] **Step 3: Implement**

- [ ] **Step 4: Verify the two gate criteria clear**

Run the gate and report the before/after criteria table. The cap criterion and the engagement-step
criterion are the two expected to move; if a third moves, say which and why.

- [ ] **Step 5: Confirm the audit still replays the result**

`audit_toolpath_engagement` must accept the new operation stream. If it raises
`UnexpectedToolpathGeometryError`, that is this task's problem, not a pre-existing limitation to
route around.

- [ ] **Step 6: Commit**

```bash
pixi run baseline
git add src/compas_cgal/engagement_toolpath.py src/compas_cgal/engagement_radial_toolpath.py \
        tests/test_engagement_toolpath.py tests/test_engagement_radial_toolpath.py benchmarks/gate.py
git commit -m "feat(engagement): non-slotting entry, so the cap holds at chain entry too"
```

---

### Task 37: break the ~141.8° cap saturation

**The item that decides whether this is engagement control or a guardrail.** Asking for a 20° cap
yields the same measured maximum engagement as asking for 140°, with 264 forced advances at 20°. A
cap that does not bind below ~141.8° is not control of the engagement angle; it is a limiter.

The cause is diagnosed in codex's own `docs/benchmarks.md`, and it is not a bug in the regulator:

> It regulates **advance** but not **trochoid radius**, so once the radius is fixed by the
> medial-axis clearance no advance can reach a tight cap. Held regulates spacing **and** takes
> smaller circles from the MAT.

`radius_regulated_toolpath` was built precisely to answer that and **did not move it** — it is inert
at a 120° cap. That failure is the most informative fact available and Step 1 is to understand it
before writing any code.

**Depends on Task 34.** If the guide's station set is the constraint, a medial axis that accepts
exactly one hardcoded polygon cannot offer a different set to choose from.

**Files:**
- Modify: `src/compas_cgal/engagement_radial_toolpath.py` (most likely), possibly the guide
  construction it consumes
- Modify: `tests/test_engagement_radial_toolpath.py`
- Modify: `benchmarks/` — the cap sweep that demonstrates the result

**Interfaces:** consumes the medial-axis guide and its station set; produces no new public API unless
Step 1 concludes the guide itself must change shape, in which case report before implementing.

- [ ] **Step 1: Establish why radius regulation did not move the saturation**

This is the whole task's leverage. Candidate explanations, to be distinguished by measurement rather
than argued:

1. **The radii offered are already minimal** — MAT clearance at each station bounds the loop radius
   from below, and the regulator is choosing the smallest available. Then the station *set* is the
   constraint, not the regulator, and the fix is upstream in the guide.
2. **The regulator's acceptance rule saturates** — it stops reducing once some other predicate is
   satisfied, so smaller radii are available and never taken.
3. **Engagement is not monotone in the loop radius** — this is already known and pinned by
   `test_engagement_is_not_monotone_in_the_loop_radius`, so a regulator that assumes monotonicity
   would stall at a local optimum rather than the global one.

Measure which. For a station where the cap is missed, enumerate the radii actually available and the
engagement each would produce, and say whether a feasible radius existed that the regulator declined.

**If the answer is (1), stop and report** — the fix is in the guide, likely behind Task 34, and it is
a different task from this one.

- [ ] **Step 2: Write the failing test — a cap sweep, not a single cap**

```python
@pytest.mark.parametrize("cap_deg", [40.0, 60.0, 90.0, 120.0])
def test_measured_engagement_tracks_the_requested_cap(cap_deg):
    """A cap that does not bind is not engagement control.

    Before this task the measured maximum saturated near 141.8 deg: asking for 20 deg
    produced the same maximum as asking for 140 deg, with 264 forced advances at 20 deg.
    """
    result = radius_regulated_toolpath(POCKET_12X8, tool_diameter=2.0, tea_cap_deg=cap_deg)
    peak = _measured_peak_engagement_deg(result)
    assert peak <= cap_deg + CAP_MEASUREMENT_SLACK_DEG, f"asked {cap_deg}, measured {peak}"
```

`CAP_MEASUREMENT_SLACK_DEG` is a named constant with a derivation — it covers the sampling density of
`_measured_peak_engagement_deg`, nothing else. It is **not** a tuning knob: if it needs to grow to
make a case pass, the case has not passed.

- [ ] **Step 3: Implement whatever Step 1 identified**

- [ ] **Step 4: Report the sweep, including the cost**

Held's length grows 40–60× from a 160° cap to a 20° cap; this generator's grows 6.8×, and codex's own
`docs/benchmarks.md` correctly reads that flatness as the saturation seen from a second direction. So
**a successful fix should make the length curve steeper.** Report the growth factor before and after;
if it is still flat, the saturation is not broken no matter what the peak-engagement test says.

- [ ] **Step 5: Commit**

```bash
pixi run baseline
pixi run python -m benchmarks.cli --corpus all --out docs/benchmarks
git add src/compas_cgal/engagement_radial_toolpath.py tests/test_engagement_radial_toolpath.py benchmarks/
git commit -m "fix(engagement): the cap binds below 141 deg"
```

---

## Self-Review

**Spec coverage.** C1 → Task 1. C2 → Tasks 6 (harness) and 7 (fix); the spec's "not proven
unsound end-to-end" nuance is carried into Task 6's `test_bound_holds_on_stock_the_shipped_audit_actually_produces`
and Task 7's inverted flat-bound test, so both the hole and the accident that hides it are pinned.
C3 → Task 12. I1 → Task 10. I2 → Task 8. I3 → Tasks 5 and 10. I4 → Tasks 3 and 4. I5 → Task 2.
I6 → Task 9. I7 → Task 19. I8 → Task 13. I9 → Task 14. I10 → Task 15. I11 → Task 20. Notes:
sentinel → Task 16; `external_tangent` → Task 17; duplicated validator → Task 18; `__all__`,
frozen dataclasses, ruff format, `--doctest-glob` → Task 20.

**Deliberately not fixed, with reasons.** The 20-parameter C++ signature and the 9-tuple return
(spec Notes) are a genuine type-design defect, but changing them touches every call site and every
test at once — that is a separate plan, not a task here, and nothing in the audit shows them
producing a wrong answer. `Stock2::set()`'s non-const accessor and `Stock.raw` stay: `engagement_at`
needs the non-const arrangement handle for the zone query, and closing the hole properly means a
borrow type, which again is its own change. `subtract_arc_sweep`'s undocumented `std::max(4, …)`,
`CHAIN_SLACK_FRACTION`'s cross-module rationale, the disk-chain placement-error note, and the bare
`assert isinstance` calls are comment-level and are folded into whichever task next touches those
lines; none is worth a dedicated test cycle.

**Ordering constraints.** Task 0 blocks everything. Task 1 changes
`validate_toolpath_params`' signature and Task 4 extends the same function — do 1 before 4.
Task 5 must precede 6, 7 and 10 (they call the bound through the binding). Task 7 must precede 10
(`certify_arc_recursive` uses `tea_growth_bound_curved`). Task 11 must precede 20 (`__all__` lists
`EngagementVerdict`). Task 18 must precede 20 (the `__all__` sweep). Tasks 13, 14, 15, 16, 17, 19
are independent of each other.

**Type consistency.** `validate_toolpath_params` gains `mat_scale, radial_clearance` in Task 1 and
is not reordered afterwards. `tea_growth_bound(d, r)` keeps its two-argument form throughout;
`tea_growth_bound_curved(d, r, rho)` is the new three-argument one — Task 7 introduces it, Task 10
consumes it, and no task renames either. `CertifiedTea` is the return type of both
`certify_segment_tea` and `certify_arc_tea`. `_certify_arc_engagement` keeps its
`(max_tea, cap_certified, stations)` tuple across Task 10's rewrite, so `_replay_arc` is unaffected
apart from the added `clockwise` argument. `OperationEngagement.cap_certified` survives Task 11 as a
property, so every existing reader keeps working.
