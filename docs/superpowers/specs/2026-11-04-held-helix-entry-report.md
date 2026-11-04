# Held helix entry: implementation report (2026-11-04)

Every Held workload now opens with a HELIX entry ahead of the published circle 0. Every junction of the entry is gated at Q4 θ_e (5.585°), up to and including the hand-off into circle 0.

The published circles and connectors are byte-identical: `tools.held_path_digest --check` reports IDENTICAL on all four digest workloads. The digest record gains an additive per-workload `entry` digest. R6 kills 4 of 4 mutants; the control passes and the import proof holds. Nothing is committed.

## What the entry is

- **Placement.** The placement is decided exactly by comparing ρ0 with r (`benchmarks/held_helix_entry.py`):
  - **ρ0 ≤ r** (figure 5, upper lane, upper boundary): the helix descends on circle 0 itself through `HeldCircleFootprint`. This footprint is circle 0's own doubles (not re-authored), so the bore is exactly circle 0's swept disk. There are no ramp-out loops, and the hand-off is circle 0 at its own contact.
  - **ρ0 > r** (crossed skis, Monstera): the pencil root is `inscribed_root`. It is the `IncidentLoop` on circle 0's entry chart with the largest grid radius whose swept disk lies exactly inside circle 0's. Then:
    - The helix descends on `TangentPencil.on_grid(root).footprint` (ρ_h ≤ r).
    - `RampOut.next_member` walks the pencil furthest-first, using the unchanged `enter_first_cut` machinery: coverage bound, then θ_e, then the cap.
    - When the scan admits the root, circle 0 is decided in its place by θ_e at station 1, by the cap, and by an exact core-coverage check. A straight link (gated by the cap) runs from the pencil entry P_R to circle 0's contact.
- **Exact entry P.** Circle 0's contact q misses its circle by −3.3e-18 mm² on figure 5 and +1.9e-18 mm² on upper. The exact entry is therefore defined as P = c0 + ρ0·u/|u| with u = q − c0, which is a one-root point. `HeldCircleFootprint.certifies_incidence` certifies, exactly over the rationals per coordinate, that q is P's nearest double (the neighbour-double interval, repeated squaring). The helix starts and ends at q. No snap.
- **`Helix` reuse via a protocol.**
  - `helix.HelixFootprint(Protocol)` (centre, radius, entry, `certifies_incidence`, `circle`) now types `Helix.footprint`; `IncidentLoop` satisfies it as is.
  - `incident_loop.GuideLoop(Protocol)` widens the `Admissibility` loop queries (`loop_within_cap`, `loop_entry_within_rise`), so circle 0 is decided with the same native predicate.
  - The two consumers that encode the chart parameter (`engagement_audit/operation_identity._helix_bytes` and `tools/pipeline_timing._helix_data`) now narrow by `isinstance` and raise a named error. `held_path_snapshot` already checks `type(...) is IncidentLoop`.
- **Declarations.**
  - `units.CutDepth` (positive, with the stock-top datum and `cut_plane(stock_top)`).
  - `HeldReferenceCase.cut_depth = DEFAULT_HELD_CUT_DEPTH` (1 mm = r).
  - `HeldReferenceCase.centre_cutting = DEFAULT_HELD_CENTRE_CUTTING` (True); `False` raises `post.errors.NonCentreCuttingHelixError`.
  - Stock top `HELD_STOCK_TOP = 0`, cut Z −1.
  - Ramp: `DEFAULT_HELIX_PITCH_PER_RADIUS`, the 2° D1 default, K = 0.2194.
- **Containment.** Every entry sweep lies exactly inside circle 0's swept disk:
  - on circle 0's own footprint, the bore *is* that disk;
  - on the pencil, the root's swept disk is inscribed in it, and every member's disk lies inside the root's.
  - The entry therefore inherits circle 0's own (published) containment. I did not use `TangentPencil.build`'s pocket certificate: `entry_containment` (`ReachableDomain`) did not finish within 170 s on figure 5's 65-vertex polygon, measured today with the call killed, and it is unusable inside a generator stage.
- **Wiring.**
  - `GeneratedToolpath.entry: HeldHelixEntry` is built in `_emit_closed_polygon_path` as a new last stage, `"entry"`, which takes 0.01 to 0.08 s.
  - `benchmarks/tools/held_sequence_cli.py` passes it too.
  - `HeldHelixEntry.operations()` rebuilds the HELIX, the CUT loops and the link with the pipeline emitters.
  - The record itself is frozen dataclasses only, so its digest is canonical.

| workload | footprint | ρ_h (mm) | turns | pitch (mm) vs K·ρ_h | ramp-out loops | hand-off link |
|---|---|---|---|---|---|---|
| figure5_refined | circle 0 | 0.445174 | 11 | 0.0909 ≤ 0.0977 | 0 | none |
| figure8_upper_lane / _boundary | circle 0 | 0.316993 | 15 | 0.0667 ≤ 0.0696 | 0 | none |
| figure8_crossed_skis_boundary | pencil (root 1.2838315, inscribed) | 0.99999999 | 5 | 0.2 ≤ 0.2194 | 1 (ρ 1.13383) | 0.27 µm, CUT |
| figure8_monstera_boundary_qualified_corners | pencil (root 1.0208228, inscribed) | 0.99999985 | 5 | 0.2 ≤ 0.2194 | 0 | 0.22 µm, CUT |

Log: `build/helix-entry/probe_entry*.log`.

## Verification

- **Red first:** `build/helix-entry/red.log` shows a collection error, because the module did not exist yet.
- **New tests:** `tests/benchmarks/test_held_helix_entry.py` covers:
  - the first operation is a HELIX within the rulings;
  - the pitch is ≤ K·ρ and the turn count is minimal;
  - the bore equals circle 0's swept disk (`Stock.exactly_equals` and NOTHING), for ρ0 ≤ r;
  - the pencil entry for ρ0 > r, with an exact inscription check;
  - an independent replay that checks θ_e and the cap at every junction, including the hand-off into circle 0, on all 5 workloads;
  - the non-centre-cutting cutter is refused;
  - the entry digest is deterministic;
  - the existing digests are unchanged and the entry digest is recorded.
- **Test runs:**
  - `pytest tests/benchmarks/test_held_helix_entry.py tests/tools/test_held_path_digest.py tests/benchmarks/test_held_polygon_nothing_stage.py tests/benchmarks/test_held_toolpath_benchmark_connector_gate.py -n 2`: **51 passed** (`build/helix-entry/tests_green2.log`).
  - Related existing tests (helix, admissibility, digest): **154 passed** (`build/helix-entry/unit_related.log`).
  - Regression subset: `tests/engagement_audit tests/post tests/test_incident_loop.py tests/pipeline tests/benchmarks -k "helix or snapshot or held_reference_case or sequence or digest or entry"`, **657 passed** in 315 s (`build/helix-entry/regression.log`). Only the last `-k` applies, to all five paths.
- **Digests.** `python -m tools.held_path_digest` rewrote the record. Its `git diff` shows only the 4 added `"entry"` lines; every existing field is unchanged. `--check` then reported **IDENTICAL** on figure5_refined, figure8_crossed_skis_boundary, figure8_upper_boundary and figure8_monstera_boundary_qualified_corners (`build/helix-entry/digest_check.log`). The entry digests match a separate earlier probe process, so the entry is deterministic across processes.
  - The `WorkloadDigest.entry` field is optional. A record without it still parses, and `to_json` omits it when it is None.
- **mypy --strict:**
  - types-pipeline, types-audit and types-benchmarks all report Success (15, 11 and 86 files).
  - `benchmarks/held_helix_entry.py` and `tools/held_path_digest.py` also report Success.
- **ruff check and format:** clean on all touched files.
- **R6** (`build/helix-entry/r6/run_r6.py`): each mutant's `benchmarks/` is prepended to the package path, and `r6_plugin.py` prints `LOADED <file>`, which is checked against the mutant copy. Diffs and logs sit beside it. Result: R6 OK.

  | mutant | result | killed by |
  |---|---|---|
  | control | PASS | (none) |
  | M1 drop the helix | KILLED | first-op, bore and junction tests |
  | M2 footprint radius ρ0 + 1 ulp | KILLED | bore and digest tests |
  | M3 turns − 1 | KILLED | first-op, pitch-minimal and digest tests |
  | M4 cutter refusal dropped | KILLED | refusal test |

## CLS export and the tri-dexel app

- **Export.**
  - `build/tridexel-cut/export_path.py` now writes `held-tridexel-path/2` with an `"entry"` record; `"moves"` is unchanged.
  - `to_cls.py --helix` writes `build/tridexel-cut/figure5_refined_helix.cls`: 195,388 GOTOs, of which **1,056 change Z** (11 turns × 96 chords, Z 2.0 → 1.0 linear in arc parameter). Then follow the flat closing turn, any ramp-out loops and the hand-off link at Z 1, then the published knots.
  - The plain `figure5_refined.cls` is unchanged in content: 194,236 GOTOs. Only its header commit changed, to 681b61f4, because HEAD moved during the run.
- **What the app needs.** Its exact XY sweep engine refuses Z-changing segments, so it needs a 3D segment sweep for a flat end mill. For the dexel field, a chord removes, per dexel, everything above the lowest tip Z of the tool disk over that dexel, with Z interpolated along the chord.
  - The simplest route is to split each Z-changing GOTO into an XY sweep at the chord's end Z. That over-removes by at most one chord's descent (0.09/96 mm here), so it is not exact.
  - The exact route is a new XY-sweep removal kind with a linear floor.
  - Either way, the removal-history stream (`ordered-xy-sweep-removal-history-v1`) needs a segment kind carrying z0 and z1.
  - The cam repo is untouched.

## Notes and open items

- Editing `held_reference_cases.py` and `units.py` changes the reconstruction-cache identity (`DECIDING_MODULES`). The first load per case after the change reconstructed again; that was minutes under load and is a one-off.
- **Proposed gate line (not applied; pyproject is fenced):** add `benchmarks/held_helix_entry.py tools/held_path_digest.py` to `types-benchmarks`.
- **Not done:** the mkdocs page for the Held entry. `docs/` pages are shared and outside this brief.
- `held_toolpath_baseline.json` metrics (path length) do not include the entry. `_metrics` is unchanged, so the baseline does not move.
- **Superseded:** the stop report from earlier today. The rulings answered both blockers. Circle 0 now follows the entry, so the circle 0 → 1 step stays an ordinary Held step under the cap and figure 5's former 12.38° at circle 1 is no longer an entry junction. The pencil handles ρ0 > r.
