# Held native declared clearance

**Conclusion.** The exact clearance oracle against the declared (Pythagorean)
boundary runs natively in C++/CGAL, and it is the default everywhere the Held
Ω\* lane decides clearance:

- containment of every Ω\* circle during generation;
- the gouge audit of the Ω\* path;
- the gouge audit of the polygon-Ω path.

The Python oracle `benchmarks/held_declared_clearance.py` is kept beside it,
unchanged, as the differential reference. `ClearanceBackend.PYTHON` selects it
at every call site.

On all five Held workloads, the native and Python oracles make **identical
decisions** on the polygon path and on the Ω\* path:

- the same verdict per element;
- the same reported worst primitive, primitive kind and segment;
- reported depths that are bit-identical.

The Ω\* paths audit to **0 gouging circles and 0 gouging connectors**. The
numbers are in [Evidence](#evidence).

## What the oracle decides

A path element is a locus of the tool centre. There are three kinds:

- a binary64 machining circle;
- a segment of an exact ∂Ω\* connector;
- an arc of an exact ∂Ω\* connector (one-root endpoints, rational support).

For an element e and tool radius r, the oracle decides the sign of
dist(e, ∂P) − r exactly, against the declared primitives:

| sign | verdict | meaning |
|---|---|---|
| < 0 | `GOUGE` | the swept tool disk leaves the pocket |
| = 0 | `TOUCH` | the tool runs exactly on the wall (allowed) |
| > 0 | `CLEAR` | |

A distance is the minimum over closed-form critical pairs. These are each
endpoint against the other curve, plus the interior pairs with a common normal.
Crossings are decided exactly, never assumed away. Each candidate comparison is
the sign of a number in one quadratic extension. A sign that mixes two radicals,
u + √β·w, decomposes into `sign(u)`, `sign(w)` and `compare(u², β·w²)`. Nothing
in a decision path is a double.

Which primitive is reported follows the Python oracle exactly: worst sign
first, then the smallest reported distance, then the lowest index. The reported
depth r − dist is a binary64 view that never feeds a decision. The rational is
rounded to its nearest double, ties to even: CGAL's `to_interval` gives a
bracket, and one exact comparison with the bracket midpoint picks the side.
In the subnormal range `to_interval` is not tight, so the enclosure is first
narrowed by exact bisection. Hypothesis found that case.

**Number model.** The native oracle uses the circle-segment traits' own one-root
numbers: `Sqrt_extension` over Epeck's lazy FT, shared with
`_held_centre_domain_2`.

- Same-extension arithmetic is enforced in release builds
  (`CrossRootArithmeticError`).
- A rational-square radicand is folded via `CGAL::is_square`.
- Binary64 lanes are injected exactly at the wire
  (`benchmarks/held_declared_clearance_wire.py` refuses any value that is not
  exactly a double). Rationals travel as text.

**Pruning only prunes.**

- An `AABB_tree` holds the primitive boxes, taken from `CGAL::to_interval`, so
  they are outward-rounded.
- The query box is grown by r in `Interval_nt<true>`.
- Arcs are boxed by their swept extent: the endpoints, plus each axis extreme
  the sweep spans, with the spanning decided exactly.

A pruned primitive is provably farther than r.

## Owner-by-construction TOUCH

Each piece of ∂Ω\* lies on the exact r-offset of the declared feature it offsets.
`CentreDomainPiece.owners` names those features (`BoundaryFeature` LINE i,
ARC i, or JUNCTION k). The audit sends the LINE and ARC owners with each
piece. For an owner, the oracle may decide TOUCH structurally instead of by
candidates, but only where a lemma proves it. Every hypothesis is checked
exactly. If any one fails, the pair is decided in full like any other pair.

- **Lemma L (line owner).** Let L = S + λu, λ ∈ [0, Λ], |u| = 1, n = u⊥.
  - *Hypotheses.* Both segment ends p satisfy n·(p − S) = r and
    u·(p − S) ∈ [0, Λ].
  - *Proof.* Every point x of the segment is a convex combination of the ends.
    So n·(x − S) = r and its foot lies on L. By Pythagoras, dist(x, L) = r.
    The pair is a TOUCH.
- **Lemma A (arc owner).** Let A be the arc of circle (C, ρ) from S to E in
  sense σ, not a full circle.
  - *Hypotheses.* The element lies on a concentric circle of radius R′ with
    |R′ − ρ| = r. Read in A's sense, its sweep runs forward from its first to
    its last direction inside A's span [S, E]. This takes three exact
    `compare_on_arc` calls, and only directions from C count.
  - *Proof.* For every element point x, the ray from C through x meets A, at
    distance | |x − C| − ρ | = r. The pair is a TOUCH.

JUNCTION owners have no lemma. Those pieces are decided in full against every
primitive.

**Soundness against false claims.**

- A non-owner primitive is always decided in full, so an owner claim cannot
  hide a gouge by another primitive.
- A false owner claim fails a hypothesis and falls back to a full decision.
- Owner touches still compute their candidates, so the reported primitive and
  distance are the same as without the lemma.

Two tests cover this:

- `test_an_offset_segment_beyond_its_owner_is_decided_in_full` is the test that
  catches a lemma that skips its hypotheses (defect injection).
- The false-owner test names every primitive as owner and checks that no
  verdict changes.

**Owner order is canonical.** `owners` is an `OwnerSet`
(`benchmarks/held_centre_domain_supports.py`), a frozenset that iterates and
prints in the exact order (kind in declaration order LINE, ARC, JUNCTION, then
index). A plain frozenset iterates in hash-table order, and the enum hashes its
name with the per-process string hash seed. Before this change, identical Ω\*
paths had different `repr` digests across processes, and attribution output
built by iterating owners could differ too.
`tests/benchmarks/test_held_omega_determinism.py` generates the fixtures under
four hash seeds and asserts one digest. The same file pins `OmegaSite`
equality: it is structural, and the lazy point cache takes no part in it.

## How the Ω\* paths use it

`generate_omega(workload, backend=ClearanceBackend.NATIVE)`
(`benchmarks/held_omega_toolpath.py`) builds the native oracle in its `omega`
stage, so the build is timed. The toolpath records the backend in
`generation_backend`. Containment goes through the `CircleClearanceOracle`
protocol (`circles(elements) -> list[ClearanceVerdict]`), which both oracles
satisfy:

- **Initial containment.** `contain_all` (`benchmarks/held_omega_circles.py`)
  decides the first containment of every circle in **one native batch**. A
  circle is contained when its verdict is not GOUGE and one of its points is
  located in Ω\*. A circle that fails is shrunk along its inward normal:
  halving, then a bisection over the bits, one oracle call per probe.
- **Midpoints.** Refinement (`benchmarks/held_omega_refinement.py`) contains one
  midpoint per call. Midpoint k + 1 depends on whether midpoint k, once
  contained, is within the engagement cap of the last emitted circle, so these
  calls cannot be batched.
- **Audit.** `audit_omega_path` and `audit_polygon_path`
  (`benchmarks/held_omega_gouge_audit.py`) both default to NATIVE. The Ω\*
  audit decides each whole ∂Ω\* piece any connector runs on in one call, with
  its owners, then the sub-pieces of gouging pieces in a second call. A
  sub-piece is a subset of its piece, so its clearance is at least the piece's.

## Evidence

Generated by `tools/held_native_clearance_results.py --measure`, which runs
`benchmarks/tools/held_declared_clearance_differential.py` on every workload.
The record is `benchmarks/results/held_native_clearance.json`. Times are wall
seconds on a shared machine, and native times include the native build. Only
same-process ratios mean anything; absolute times move with the load printed
beside them.

--8<-- "benchmarks/results/held_native_clearance.md"

How to read it:

- **Polygon Ω path.** It gouges, by the chord sagitta at concave arcs (see
  [Held exact centre domain](held_exact_centre_domain.md#evidence-the-old-domain)).
  Both oracles find the same gouges, element for element.
- **Ω\* path.** No element gouges. Circles are TOUCH only where the rounded
  centre landed exactly on the normal, and CLEAR otherwise. Connectors run on
  ∂Ω\*, so they TOUCH by construction.
- **Generation.** The oracle decides containment only in the refine stage.
  Initial placement does not call it, and it dominates whole-generation time.
  So the refine column is the one this change moves.

## Tests

- `tests/benchmarks/test_held_declared_clearance_native.py` covers:
  - exact touch and gouge fixtures (2⁻⁴⁰ perturbations);
  - lemma TOUCH on every Ω\* piece of three fixtures;
  - false owner claims, and an offset segment beyond its owner;
  - Hypothesis over binary64 circles and segments, Pythagorean arcs and
    one-root-endpoint elements;
  - batch results equal to singleton results;
  - named errors.
- `tests/benchmarks/test_held_declared_clearance_wire.py`: the exact wire.
- `tests/benchmarks/test_held_native_clearance_default.py`:
  - native is the default in both audits and in `generate_omega`;
  - `contain_all` makes exactly one batch and matches `contain`;
  - the shrink search finds the same circle on both oracles;
  - refinement is `repr`-identical across the two oracles.
- `tests/benchmarks/test_held_omega_determinism.py`: the hash-seed-free digests
  and the canonical owner order, and `OmegaSite` equality.

## Reproducibility

- `pixi run held-native-clearance-results -- --measure` re-measures (minutes).
  Without `--measure`, the task only re-renders
  `benchmarks/results/held_native_clearance.md` from the JSON.
- `pixi run held-native-clearance-differential` is the raw differential alone.
  It writes `build/native-clearance/differential.json`.
- Native sources are in `src/held_declared_clearance_2/`, one responsibility
  per file. The lemmas are stated in `owner_lemma_2.h`.
