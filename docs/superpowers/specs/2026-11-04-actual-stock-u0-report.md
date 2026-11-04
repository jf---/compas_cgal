# U0 report: actual-stock differential oracle (2026-11-04)

Spec: `2026-11-04-actual-stock-refinement-design.md` §3 (D2, D4), §5 (U0). Worktree HEAD `ea4703bb`, nothing committed.

**Exit gate met.** D4 = 0 on fig5 (every circle), skis (every circle) and Monstera (every 20th circle). D2 holds on
skis for every circle. R6 killed 18 of 18 mutations, and all 4 controls behaved as expected.

## Delivered

| file | change |
|---|---|
| `benchmarks/held_actual_stock_last16.py` | +27 / −6. `SweepView` {ACTUAL, CIRCLES} option on `Last16Oracle` (CIRCLES = the circle view C_t: `candidates` drops the connectors). The stock is now built from the typed `self.material: PartPolygon`. Docstring: Omega is P, the oracle is off the production path. |
| `benchmarks/held_material_domain.py` | new, 112 lines. `PartPolygon` / `CentreDomainPolygon` (distinct nominal types over `Point2[WorldXY]`, `build` validates). `PartPolygon.from_case` reads `case.projection.points`. `require_part_material` is an exact row check. Errors: `DegeneratePolygonError`, `MaterialDomainMismatchError`. |
| `benchmarks/held_actual_stock_oracle.py` | new, 366 lines. The typed interface later units call (below). D4 / D2 harnesses, a spawn pool whose results come back in order, and a CLI. |
| `tools/actual_stock_circle_view_probe.py` | new, 171 lines. The design probe `probe_circles_only.py`, ported onto `SweepView`. Typed `ProbeRecord` with a JSON round trip and schema check. |
| `tests/benchmarks/test_held_actual_stock_oracle.py` | new, 19 fast tests |
| `tests/tools/test_actual_stock_circle_view_probe.py` | new, 4 fast tests |
| `benchmarks/acceptance/test_actual_stock_oracle_slow.py` | new slow tier (not collected): D4 for skis, fig5 and Monstera/20; D2 for skis |
| `build/actual-stock-u0/r6/run_r6.py` | R6 harness |

### Interface for U1–U5

- `Last16Oracle(workload, last=16, sweeps=SweepView.ACTUAL | CIRCLES)` and `.material: PartPolygon`.
- `workload_on_part(path, generation_s)` and `differential_workload(name)` build a `ReplayWorkload` and check exactly
  that it is on P (`MaterialDomainMismatchError` otherwise).
- `view_prefixes(oracle, index, last)` is the D1 reference side: `engagement_circle_prefixes` on the circle-only window on P.
- `circle_view_vs_actual(workload, indices, processes, audit_every)` returns a `D4Report`.
  - `.require()` raises `CircleViewRefusalError`.
  - The harness raises `BoundBreakError` if S_t refuses a circle that C_t certifies.
- `production_implies_certificate(workload, indices, processes, sweeps)` returns a `D2Report`.
  - `.require()` raises `EntryCircleNotRefusedError` or `CertificateViolationError`.
  - With `sweeps=CIRCLES` it checks the stronger circle-view statement that U4's decision makes.

Cost structure: D4 decides C_t with the last-16 → full-window escalation. It decides S_t only where C_t refuses, plus
every `audit_every`-th circle as a bound check. Circle 0 has the same stock in both views, so it never counts.

## Evidence

Red first:

- `build/actual-stock-u0/red.log`: ImportError, fast and slow tiers (D4 on skis).
- `red_probe.log`: the same for the probe.

Fast tier:

`nice pixi run --no-install sh build/pyrun.sh -m pytest tests/tools/test_actual_stock_circle_view_probe.py tests/tools/test_actual_stock_profile.py tests/benchmarks/test_held_actual_stock_last16.py tests/benchmarks/test_held_actual_stock_oracle.py -n 2 -q`
gives **46 passed**.

Slow tier (`build/actual-stock-u0/slow_gate.log`, load 18 → 131, 2 processes, `-p no:xdist`): **4 passed in 136 s**.

| check | workload | circles | C_t refused | S_t decided | violations | circle CPU |
|---|---|---:|---:|---:|---:|---:|
| D4 | skis | 1,621 | 1 (circle 0) | 33 | **0** | 37.0 s |
| D4 | fig5 | 2,004 | 1 (circle 0) | 41 | **0** | 46.9 s |
| D4 | Monstera /20 | 222 | 1 (circle 0) | 45 | **0** | 6.7 s |
| D2 (actual) | skis | 1,621 | entry refused, 1,620 last-16 | — | **0** | wall 61 s |

Probe parity: `tools.actual_stock_circle_view_probe figure5_refined --stride 40` gives 51 rows. The actual and
circle-view max runs match `build/actual-stock-design/probe_fig5.log` on every row (to the log's 2 dp). The gap is
p50 0.00, p90 2.60 and max 21.01°, the same as the design.

Types:

`pixi run --no-install mypy --strict --warn-unused-ignores benchmarks/held_material_domain.py benchmarks/held_actual_stock_oracle.py benchmarks/held_actual_stock_last16.py tools/actual_stock_circle_view_probe.py`
reports **Success: no issues found in 4 source files**. The three `benchmarks/` modules are also inside the bare
`mypy` gate (`files = ["benchmarks", …]`).

Lint: `ruff check` and `ruff format --check` pass on every new or edited file.

R6 (`build/actual-stock-u0/r6/r6.log`): **22/22 as expected**.

Controls:

- 3 copies raise `ImportError` with the r6 marker, from the copy's path.
- 1 unmutated copy is green (35 passed). The import locations of all five modules are printed under `m03/`.

Mutations (18, each red on its named node):

- the circle view keeps the connectors;
- **material domain = offset polygon** (`from_path` reads `path.boundary`), killed by the generated-skis guard test;
- the part check is removed, or done with `allclose`;
- a 2-vertex polygon is admitted;
- the D4 violation ignores the actual verdict, or ignores the circle-view verdict;
- the actual is not decided on a C_t refusal;
- the audit sample is dropped;
- the bound-break check is removed;
- the oracle cache ignores the view;
- the worker-init check is removed;
- D2 counts the entry circle, drops the entry check, does not count full-window certificates, or ignores the view;
- in the probe, the D4 statistic ignores the actual, or the schema check is dropped.

## Bug found during build

My first insertion of `SweepView` left the `@dataclass(frozen=True)` of `WindowView` on top of the enum. Its members
then compared equal and hashed equal. `is` still worked, but the per-view oracle cache merged the two views, so D4
silently decided C_t on S_t. The D4 fixture test caught it.

## Proposed pixi lines (for the single pyproject slot; not applied)

```toml
actual-stock-oracle-slow = { cmd = '''editable_build_dir="$(python -c 'import sys; print(next(f.path for f in sys.meta_path if hasattr(f, "known_wheel_files") and "compas_cgal._stock_2" in f.known_wheel_files))')" && SKBUILD_EDITABLE_SKIP="$editable_build_dir" pytest benchmarks/acceptance/test_actual_stock_oracle_slow.py -q -s -p no:xdist''', depends-on = ["_editable-rebuild"], description = "Slow tier of the actual-stock differential oracle (U0): D4 circle view vs actual on skis, Figure 5 every circle and Monstera every 20th; D2 on skis (~2.5 min CPU); not collected by pytest" }
actual-stock-oracle = { cmd = '''editable_build_dir="$(python -c 'import sys; print(next(f.path for f in sys.meta_path if hasattr(f, "known_wheel_files") and "compas_cgal._stock_2" in f.known_wheel_files))')" && SKBUILD_EDITABLE_SKIP="$editable_build_dir" python -m benchmarks.held_actual_stock_oracle''', depends-on = ["_editable-rebuild"], description = "Actual-stock differential oracle CLI: d4|d2 WORKLOAD [--every K] [--audit-every K] [--circles] -> build/actual-stock-u0/<workload>.<check>.every<K>.{json,md}" }
actual-stock-circle-view-probe = { cmd = '''editable_build_dir="$(python -c 'import sys; print(next(f.path for f in sys.meta_path if hasattr(f, "known_wheel_files") and "compas_cgal._stock_2" in f.known_wheel_files))')" && SKBUILD_EDITABLE_SKIP="$editable_build_dir" python -m tools.actual_stock_circle_view_probe''', depends-on = ["_editable-rebuild"], description = "Circle-view probe: max-run on the actual stock vs the circle view per sampled circle -> build/actual-stock-u0/<workload>.circle_view.stride<K>.{json,md}" }
```

To add `tools/actual_stock_circle_view_probe.py` to the `types-benchmarks` path list, append it to that task's
command. The other three modules are already covered by the bare `mypy` gate.

## Not done

- I did not run D2 on fig5 or Monstera. The brief asked for skis only; fig5 and Monstera are U4's exit gate.
- I did not run D2 in the circle view on the real workloads. The harness supports it (`--circles`).
- U-1's outputs in `build/actual-stock-u1/` were not touched.
- No pyproject edit, no commit.
