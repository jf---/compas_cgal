# Census brief, 2026-10-13: shared by every slice

## Why this exists

The project builds engagement-controlled trochoidal / adaptive pocket-clearing
toolpaths with exact (CGAL) geometry. It has not converged, because every past
review read a sample and called it complete. This census reads ALL of it, with
coverage checked by a script, so the plan that follows can be complete.

**The emphasis is GENERATION.** The question for every file is: does it help or
hurt producing world-class toolpaths, ones a machinist would run on a real
machine? Validation is worth nothing without excellent paths to validate.
Read generation code deeply. Treat validation and test code structurally, and
flag it only where it blocks, distorts or wrongly constrains generation.

The bar for "world class": Held & Pfeiffer 2025 (exact engagement, contour-aware
placement, MATHSM comparison), Elber–Cohen–Drake MATHSM, and commercial adaptive
clearing (Fusion Adaptive, VoluMill, iMachining). Concretely that means:
- constant, capped engagement everywhere, entries and corners included;
- a single helical or ramped entry;
- the tool kept down through cleared stock;
- G1/G2 continuity, smooth corners;
- no air-cutting loops;
- full coverage of reachable material;
- short total path and few retracts.

## What is already known (do not rediscover; DO contradict it if your slice shows otherwise)

- Product gate: `tests/benchmarks/test_quality.py::test_the_generated_path_is_worth_running`. 12 criteria on `rect_12x8`, `rect_20x12` and `L_shape`, at caps of 120° and 40°. All 12 cases fail. The gated generators are `engagement_controlled_toolpath` and `radius_regulated_toolpath`; their behaviour has not changed since 2026-08-22.
- Root causes found so far:
  - `src/toolpath.cpp:656-661` walks every straight-skeleton edge narrower end first, and every chain gets its own plunge, so a rectangle gets 5 entries, 4 of them at corners with ρ≈0.02.
  - Walks finish inside the already-cleared centre, emitting loops that remove nothing.
  - Entering virgin stock with a circle is 360° of engagement.
  - Convex-corner material below the ρ>r floor has no clean motion.
- `rho_regulated_toolpath` (not gated) fixes the degenerate loops and the tangent breaks, but leaves 1–3.5% of the corner material uncut.
- The gate's replay counts a LINK at cut depth that removes nothing as a redundant operation, so stay-down links are currently penalised.
- The Held reproduction (`benchmarks/held_*`) reproduces Held's Figures 5, 6 and 8 on Held's own pockets, but it needs publisher-derived start data. It has never been scored on the gate.
- Draft plan (to be made complete by this census): `docs/superpowers/plans/2026-10-13-product-gate-tree-walk.md`. It walks the straight-skeleton guide as a tree from its widest vertex, emits only cuts that remove material, and repositions above the stock. Read its "Design" and "Stage B" sections and say plainly where your slice shows the plan is wrong, incomplete, or ignores a better existing capability. For example: the exact segment-site medial axis instead of the straight skeleton, or the certified `adaptive/` generator instead of sampled probes.

## Rules

- Work in `/Users/jelle/Code/CADCAM/worktrees/compas_cgal_prs-held-task5-continuation` (branch `codex/held-task5-continuation`). Read-only, except your one output file `build/census-2026-10-13/out/<slice>.md`.
- No commits, no edits elsewhere, no stash or checkout.
- Do not run the test suite. If you must check a claim by running code, use `pixi run python ...` and say so.
- Never state a number you did not read or measure. Mark claims `verified` (cite file:line or command output), `unverified`, or `contradicted` (cite).

## Output: `build/census-2026-10-13/out/<slice>.md`

Start with `## Slice summary`:
- the top findings ranked by impact on toolpath quality;
- every generation capability your slice contains that the gated generators do NOT use;
- every place your slice shows the draft plan is wrong or incomplete.

Then one record per assigned file, in slice-list order, with this exact header line (the coverage checker parses it):

```
### <path> | lines <N> | read <ranges>
```

`<ranges>` lists the line ranges you actually read, e.g. `1-958`, or `1-200,300-958`. Generation, kernel and doc slices must read every line (`1-N`). The validation slice may write `structural` (module docstring, public API and every rule that constrains generation), plus full ranges where it read more.

Then these bullets:

- **role:** a generation sub-role (`guide`, `placement`, `ordering`, `linking`, `entry`, `corner`, `engagement-control`, `spacing`, `motion-geometry`, `kernel`), or `validation`, `figure`, `tool`, `test`, `doc`, `build`.
- **used by:** entry points that reach this file (grep the imports and bindings, pixi tasks, other modules), or `no consumer found`.
- **state:** `live`, `parked`, `dead`, or `superseded-by <x>`.
- **generation:** what toolpath behaviour it provides. Quality: `excellent`, `adequate`, `weak` or `broken`, with the reason, judged against the world-class bar above. For validation and doc files, write `n/a` or the generation idea it records.
- **defects:** `file:line` — what is wrong — effect on generated paths (or `none found`).
- **unused capability / ideas:** what exists here that the product generators should use, or that a doc proposed and nobody built.
- **claims:** each docstring or doc claim that matters, marked `verified` / `unverified` / `contradicted`, with evidence.

Doc slices also keep a findings register at the end:

```
F-<slice>-<n> | source file:line | finding or recommendation | open/closed (evidence: commit or test) | generation-relevant y/n
```

## Final message to the orchestrator (≤ 250 words)

- First line: `model: <the exact model name and ID your system prompt names>`.
- Then: files recorded / files assigned, and any file you could not fully read and why.
- Then: your top 5 generation findings, each with `file:line`.
