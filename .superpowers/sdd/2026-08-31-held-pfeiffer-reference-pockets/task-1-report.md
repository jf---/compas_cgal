# Task 1 report

Status: DONE

Commit: `ee8df08`

Implemented typed PDF/source/world geometry, exact-cycle validation, circle-first
reconstruction, equal-distance G1 biarc fallback, Bernstein polar-monotonicity
certificate, closed radial/control-hull bounds, and bounded polygon projection.

Verification:

- `pixi run --frozen ruff check ...`: passed.
- `pixi run --frozen types-benchmarks`: strict mypy passed.
- focused pytest with `-n auto`: 10 passed.
- `pixi run --frozen affected`: 2288 passed, 16 branch baseline failures; no
  Held-reference failure.

Concern: reconstruction and proof recursion both have a depth ceiling but no
shared node budget. The four publisher cases must record measured primitive
counts before deciding whether a separate work budget is necessary.

## Fix round 1

Status: DONE

Commit: this report's commit

RED evidence:

- centre-local arc regression: 2 failed because malformed translated and
  untranslated sweeps were accepted;
- near-parallel biarc regression: failed with start-tangent residual
  `8.429369702178815e-08` above normalized roundoff;
- emitted-coordinate projection regression: translated circle was accepted at
  the `ulp(2**40)` bound.

Implemented centre-local rotated-radius validation, vector-residual tangent
checks at both endpoints and the biarc join, guarded near-parallel general
solving, and projection measurement from stored chord endpoints in centre-local
coordinates. Static coverage now constructs a valid typed boundary and checks
both `PdfPoint2.build` overloads plus every public return.

Verification:

- focused pytest with `-n auto`: 15 passed;
- Ruff format/check: passed;
- `types-benchmarks`: strict mypy passed;
- `git diff --check`: passed.

Concern: none within Task 1 scope.

## Fix round 2

Status: DONE

Commit: this report's commit

RED evidence:

- stored-coordinate counterexample: both CCW and CW cases failed because the
  midpoint-only projection certificate accepted them;
- collapsed emitted chord: failed with the generic deviation error instead of
  the required named collapse rejection.

Implemented the closed finite-segment maximum for every emitted circular chord.
The candidate set contains interval endpoints, both finite-segment Voronoi
transitions, endpoint-region stationary points, and both interior-region
stationary families. Linear circle equations are solved over the signed sweep
interval without sampling or angular tolerance. The existing ordinary oracle
now checks that reported deviation covers independently measured midpoints
rather than claiming midpoint equality with the continuous maximum.

Verification:

- focused pytest with `-n auto`: 18 passed;
- Ruff format/check: passed;
- `types-benchmarks`: strict mypy passed;
- `git diff --check`: passed.

Concern: none within Task 1 scope.

## Corpus extension

Status: DONE

Commit: this report's commit

RED evidence:

- normalized Monstera source 125 exhausted recursion at depth 24 before the
  chord-local, condition-certified biarc construction;
- proof-carrying reconstruction and path reconstruction imports failed before
  their public API was added;
- the first live run exposed a disconnected stored biarc seam before both
  children were mapped through one shared represented join;
- the completed live oracle then failed only on its deliberately empty census
  until the independently observed counts were hard-coded.

Implemented exact-Fraction `gamma(7)` residual certification, exact-sign root
enclosure, reflection-map and stored-arc G1 error bounds, outward continuous
fidelity bounds, certified line limits, path-level witnessed arc merging, and
explicit PDF Y reflection. Root biarc children remain merge-ineligible because
a whole-cubic witness is not valid for either child circle.

Live publisher results, recorded as `lines + arcs = total` and certified
deviation upper bound:

- Figure 5: `5 + 26 = 31`, `0.07384765937603313`;
- Figure 8 upper: `4 + 66 = 70`, `0.08364126415297131`;
- Figure 8 crossed skis: `2 + 56 = 58`, `0.09649518245132316`;
- Figure 8 Monstera: `103 + 214 = 317`, `0.10229538752363487`.

Verification:

- focused geometry pytest with `-n auto`: 36 passed;
- live publisher reconstruction oracle with `-n auto`: 1 passed;
- Ruff format/check and strict `types-benchmarks`: passed;
- `git diff --check`: passed;
- `affected`: 55 passed and the same 16 branch-baseline failures; no
  Held-reference failure.

Concern: biarc child arcs deliberately do not participate in cross-span merge
until their source cubic is split at a sound arc-length correspondence.

## Certificate fix round 1

Status: DONE; independent certificate re-review remains pending.

Commit: this report's commit

RED evidence:

- the stored-vector root witness certified the old surrogate root although its
  exact bracket signs did not enclose the represented polynomial root;
- mapped G1 bounds admitted infinity and runtime transform construction
  accepted integer/string reflection flags;
- cancellation-sensitive hull, polar, Taylor, stored endpoint, and one-ULP
  boundary regressions exposed non-closed circle/biarc arithmetic;
- Figure 5 source 5 exceeded 2,048 certificate nodes because a genuinely
  invalid root biarc was exhaustively subdivided before source subdivision.

Implemented exact stored-polynomial residual/sign certification, fail-closed
local and mapped G1 perturbation bounds, exact rational hull/polar predicates,
audited trigonometric interval evaluation, stored-endpoint-corrected biarc
correspondence, second-order same-child node bounds, and an exact pointwise
lower witness for early invalid-biarc rejection. Biarc children retain honest
exact source-parameter interval witnesses. Circle recursion retains the stored
candidate radius. A source cubic becomes a line only when its stored controls
are exactly collinear and advance along the authored chord; no live source
cubic takes that branch, so there are no cubic line-hull bounds to report.

Live publisher results (`lines + arcs = total`, certified deviation upper
bound):

- Figure 5: `5 + 44 = 49`, `0.07380041260306111`;
- Figure 8 upper: `4 + 76 = 80`, `0.08300467787418957`;
- Figure 8 crossed skis: `2 + 58 = 60`, `0.09647539567877145`;
- Figure 8 Monstera: `103 + 216 = 319`, `0.10228407943079104`.

Verification:

- focused geometry: 53 passed;
- live publisher oracle: 1 passed in 51.07 seconds;
- strict mypy, Ruff format/check, strict MkDocs, and diff checks: passed;
- affected: 72 passed, same 16 branch-baseline failures; no Held-reference
  failure.

Concern: the conservative corrected proof increases the primitive census; it
is not tuned toward the earlier unsound count. Certificate-review status stays
pending until independent re-review.
