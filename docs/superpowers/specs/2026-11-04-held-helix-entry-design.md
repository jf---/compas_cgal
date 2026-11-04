# Held helix entry: design note (2026-11-04)

The helix entry is built, but only in the narrow pipeline lane. The Held published
path (`generate("figure5_refined")`) has no entry operation. Its circle 0 starts
cutting in solid stock. The tri-dexel simulation measures that start at about 179°;
the generator's seam pair reports 180°. Adding the helix to the Held path takes more
than wiring: it needs three user declarations and edits to fenced files. No code was
changed.

## Findings

**Where the helix exists.**

- Core model:
  - `src/compas_cgal/helix.py:83-187` defines `Helix`, with `ramp_turns` at `:61` and the closing turn at `HELIX_CLOSING_TURNS` `:50`.
  - `src/compas_cgal/helix_units.py:84-146` defines `HelixPitchPerRadius` and `DEFAULT_HELIX_PITCH_PER_RADIUS`.
  - `src/compas_cgal/helix_units.py:150-190` defines `EntryRampLength`.
- Lane wiring (narrow lane only):
  - `src/compas_cgal/pipeline/generate.py:370-436` is `enter_first_cut`: helix, then tangent-pencil ramp-out, then root loop.
  - `generate_narrow` calls it at `generate.py:679`.
  - `HelicalEntry` is declared at `generate.py:673`.
- Footprint: `src/compas_cgal/pipeline/pencil.py:226,252` (`TangentPencil.footprint`, `helix_footprint`), with `rho_h <= r`.
- Consumers: replay, `benchmarks/helix_observation.py`, `benchmarks/held_path_snapshot.py:394-500` (`HeldHelixSnapshot`), and the post.

**What the Held path lacks.** `benchmarks/held_toolpath_benchmark.py:291-307` (`generate_figure5_refined`) stops at refine, NOTHING removal, connectors and the connector certificate. Its result type, `GeneratedToolpath` (`:193-213`), holds only `circles` and `connectors`. It has no operation stream, no entry and no Z.

The obstacles beyond wiring:

1. **No cut plane.** `HeldReferenceCase` (`benchmarks/held_reference_cases.py:325-344`) declares tool radius, cap and start marker, but no depth or clearance. `Helix.build` needs a `CutPlane`, so a Held helix has no defined descent.
2. **No exact footprint.** `Helix` requires an `IncidentLoop`: an exact grid centre and an entry exactly on the circle (`helix.py:113`, `incident_loop.py:94-103`).
   - Held circles are `Figure5CounterclockwiseCircle`, with binary64 `center`, `radius` and `contact_point` (`benchmarks/held_figure5_boundary_path.py:63-73`).
   - The exported circle-0 start point misses its circle by -3.3e-18 mm² in exact arithmetic (checked on `build/tridexel-cut/path_figure5_refined.json`).
   - So a Held circle cannot be a helix footprint as it stands. Either circle 0 is re-authored as an `IncidentLoop`, which changes the published path and its digests, or a Held-side footprint is certified separately.
3. **Fenced files.** Any wiring edits the files below. The brief fences them for the fig5-levers agent.
   - `benchmarks/held_toolpath_benchmark.py:193-213`: `GeneratedToolpath` gains an entry field.
   - `benchmarks/held_toolpath_benchmark.py:291-307` and `:310-331` (`_emit_closed_polygon_path`): emit the entry.
   - The Held baseline `benchmarks/results/held_toolpath_baseline.json`, whose counts and lengths change.

**Measured geometry** (`figure5_refined` @ea4703bb, from the export):

- Circle 0: centre (52.389, 2.174), radius ρ0 = 0.445 mm, which is ≤ r = 1.0 mm.
- Swept disk: radius 1.445 mm.
- Published start marker: centre (52.804, 1.873), radius 1.001 ≈ r, 0.512 mm from c0. So the swept disk covers the start marker.
- Ramp limit: at the default 2° placeholder ramp, `K = 2π·tan 2° = 0.2194`. At ρ0 the pitch must be ≤ 0.0977 mm per turn.

## Rulings (exact quotes)

- R3, plan `docs/superpowers/plans/2026-10-13-product-pipeline.md:44`: "**R3, entry:** a helical bore as its own operation. Its replay removes its disk at the cut plane, it is scored by ramp angle against the machine model, and the step from the helix to the first loop stays gated."
- R2, plan `:43`: "Held's pockets join the gate." Plan goal (`:9`): "on Held's published pockets: capped engagement everywhere, one entry, ..."
- Plan Task 1.2 (`:176`, `:183`): "Helical entry at the widest point (R3)" ... "the root chosen as the maximum clearance over the guide".
- Ledger `.superpowers/sdd/2026-10-13-product-pipeline/progress.md`:
  - `:411` "USER 2026-10-16 Q4: no plunge-junction gating; declared EntryRampLength = 2r, theta_e = cap*pi/45, decided exactly at gate station 1 on pencil + helix junctions only."
  - `:412` "USER 2026-10-16 Q5: root = 1.1 FIRST_CUT station; MA vertex to Task 2.5 (E10 not adopted)." This supersedes "widest point" for the pipeline lane.
  - `:413` "USER 2026-10-16 Q6: flat closing circle at rho_h after the helix; helix + closing circle = entry, excluded from cap sampling; bore removed at end of closing circle; post maps both."
  - `:420` "Ruling (Q6 shape): the closing turn is part of the HELIX op (ramp turns + one flat turn at bottom_z, same op)".
  - `:334` "rho_h = 1 exact on all 6 ... helix disk == bore disk exactly". This is the pipeline lane: the footprint comes from the root's tangent pencil, with `rho_h <= r`.
  - `:429` "HELIX entry also requires centre-cutting cutter (fail closed until user decides)".
- Ramp angle: `DEFAULT_MAX_RAMP_ANGLE_DEG = 2.0` (`helix_units.py:20`) is a placeholder (D1). No user ruling sets a value.

**What the rulings settle:**

- The helix footprint has `rho_h ≤ r`. The bore is the exact disk `rho_h + r`.
- The helix is one HELIX operation: ramp turns plus one flat closing turn.
- The helix and its closing turn are entry and are not cap-sampled.
- Helix and pencil junctions are gated at θ_e = cap·π/45, which is 5.585° at 80°.

**What they leave open for the Held lane:**

- where the helix sits: Q5 names the pipeline's first-cut station only, and Held has no pencil;
- the Held cut plane;
- the ramp angle, which is still a placeholder;
- how a binary64 Held circle becomes an exact footprint.

The "widest inscribed disk" placement is not current: Q5 moved it to Task 2.5.

## Recommendation

Put the helix on circle 0's own footprint. The helix descends on (c0, ρ0 = 0.445 ≤ r), entering at circle 0's contact point. Circle 0 then *is* the Q6 flat closing turn. The bore equals circle 0's swept disk exactly.

**Why this footprint:**

- The published start (circle 0 from the start marker) is unchanged.
- No ramp-out pencil is needed.
- The helix-to-first-loop step is the generator's existing circle 0 → circle 1 predecessor step, which refinement already holds ≤ cap. Only the Q4 junction check θ_e is new.
- Circle 0's 179° plunge leaves cap sampling (Q6). Whether the seam pair's 180° goes the same way (the seam re-enters a cleared bore) is unverified; the generator's seam model needs a check.

**Needed before code:**

1. User: declare the Held lane's cut plane (depth and clearance) or a per-case default.
2. User: confirm that the Held lane's root is circle 0 (the published start) and not the widest point. Q5 binds the pipeline lane only.
3. User: confirm whether the 2° placeholder stands for Held.
4. Exactness: circle 0 needs an exact footprint. Either author circle 0 as an `IncidentLoop` (centre on the incidence grid, chart-authored entry), which changes the published digests, or inject c0 and ρ0 exactly and certify a re-authored incident entry. Not by tolerance.
5. Fence lift: the fig5-levers agent must release `benchmarks/held_toolpath_benchmark.py:193-213, 291-331` before the entry field and emission land. `held_path_snapshot.py` snapshots pipeline HELIX ops (`:394-500`); whether the Held evidence path consumes one is unverified.

**Acceptance test (red first, once unblocked):**

- The first motion of `generate("figure5_refined")` is HELIX:
  - footprint == circle 0;
  - `ramp_within(DEFAULT_HELIX_PITCH_PER_RADIUS)`;
  - turns == `ramp_turns(descent, ρ0, K)`.
- Circle 0 is excluded as entry.
- The maximum sampled engagement on the remaining motions is ≤ cap, and tri-dexel circle 0 is < 80°.
- The junction from the bore to circle 1 is ≤ θ_e.
- R6 injection in `build/helix-entry/r6/`, three mutations, each of which the test must catch:
  - drop the helix;
  - change the footprint radius to ρ0 + one ulp;
  - set turns − 1.
