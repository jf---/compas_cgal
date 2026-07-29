# Exact one-root full-circle source design

## Scope

Task 13F reaches a real full-circle event boundary after an exactly certified
link has depleted a diagonal capsule. The post-link stock is valid
`Gps_circle_segment_traits_2` geometry, but some trim vertices have coordinates
of the form

```text
x = x0 + x1 sqrt(alpha)
y = y0 + y1 sqrt(alpha)
```

with exact rational `x0`, `x1`, `y0`, `y1`, and `alpha`. The current
`full-circle-boundary-pullbacks-v3` source admits a vertex only when both
coordinates collapse to rationals. It therefore returns
`full-circle-task5-blocked-v1` before the already implemented same-support and
cross-support cell theorems can run.

This stage extends the single authoritative full-circle source to shared-root
one-root vertices. It does not add a sampled path, a rounded representation, a
second oracle, or permission to search past an unresolved higher-ranked
candidate.

## Existing exact precedent

The segment event substrate already handles the same algebraic shape:

1. decompose a trim vertex into rational and radical parts;
2. form its original passage equation `A(t) + B(t) sqrt(alpha) = 0`;
3. isolate candidate parameters with the integer norm
   `A(t)^2 - alpha B(t)^2`;
4. recheck the original radical equation at every isolated algebraic root; and
5. retain only the physical branch, not the conjugate equation
   `A(t) - B(t) sqrt(alpha) = 0`.

The full-circle extension reuses that doctrine on each rational quarter chart.
It does not introduce a second coefficient domain or a general symbolic
algebra subsystem.

## Exact source grammar

Every trim coordinate, including a rational coordinate, uses one source
grammar:

```text
full-circle-one-root-coordinate-v1(
    base,
    radical_coefficient,
    radicand
)
```

All three fields are reduced exact rationals. A rational coordinate is
represented canonically with `radical_coefficient = 0` and `radicand = 0`.
There is no retained v3 reconstruction path.

One vertex source contains two such coordinates. Construction validates:

- both radicands are non-negative;
- a zero radical coefficient has zero canonical radicand;
- nonzero `x` and `y` radical coefficients use the same radicand; and
- the source vertex ID has one byte-identical coordinate record across every
  incident feature.

Distinct nonzero radicands are outside this bounded theorem and raise
`UnsupportedAlgebraicVertexProjectionError`. The caller keeps the existing
fail-closed unresolved semantics.

## Quarter-chart passage equation

For one center quarter chart, let the exact rational motion center be

```text
Cx(u) = X(u) / D(u)
Cy(u) = Y(u) / D(u)
D(u) = 1 + u^2 > 0
```

and let the cutter radius be the exact rational `rho`. After multiplying by
the known-positive `D(u)^2`, passage through the vertex is

```text
F(u) = A(u) + B(u) sqrt(alpha) = 0
```

where

```text
A = (X - x0 D)^2 + (Y - y0 D)^2
    + alpha (x1^2 + y1^2) D^2 - rho^2 D^2
B = -2 D (x1 (X - x0 D) + y1 (Y - y0 D)).
```

`A` and `B` are jointly cleared to integer polynomials using one strictly
positive common denominator, then normalized only by their shared integer
content and one shared sign convention. Independently making `A` and `B`
primitive would change their relative scale and therefore change the radical
equation. The partition factor is

```text
N(u) = denominator(alpha) A(u)^2
       - numerator(alpha) B(u)^2.
```

The factor is mapped from the local quarter chart into the existing global
parameter. Its roots refine the partition even when the corresponding
conjugate branch is nonphysical.

## Physical-root ownership

Norming is root discovery, not event ownership. For every isolated root `r` of
`N`, the algebraic kernel evaluates the exact signs of `A(r)` and `B(r)`.
The physical equation holds exactly when:

- both evaluations vanish; or
- neither vanishes, their signs oppose, and `N(r) = 0`.

Roots satisfying only `A(r) - B(r) sqrt(alpha) = 0` remain partition roots with
an empty physical event list. Endpoint events and their local witnesses are
attached only to physical roots. Each retained endpoint event records that its
original equation and orientation were rechecked.

This separation is essential: attaching one endpoint event to every norm root
would create false merge/split topology while still producing a perfectly
valid algebraic partition.

## Signed active-sheet predicate

The left/right active sheet cannot use the sign of `N`. A norm has the same
zero set as both conjugate branches and loses the sign of the physical
distance equation.

`ProjectionInput2` and `ProjectionRecord2` therefore gain one optional
`OneRootPredicate2` containing:

- jointly normalized integer coefficients of `A`;
- jointly normalized integer coefficients of `B`; and
- the exact positive rational radicand.

The generic integer partition validates that its norm equals the submitted
projection polynomial up to canonical normalization. At each rational cell
witness, the full-circle fibre classifier evaluates
`sign(A + B sqrt(alpha))` exactly by sign comparison and squared magnitudes.
The one-root predicate, physical-root filtering, and source grammar are all
canonical certificate fields.

The existing integer `signed_predicate_coefficients` remains the representation
for genuinely rational predicates. A projection owns exactly one signed
predicate form.

## Replay and versioning

The authoritative source becomes `full-circle-boundary-pullbacks-v4`.
Reconstruction accepts v4 only. `projection-record-v3` binds the optional
one-root signed predicate. `event-exact-motion-oracle-v5` and
`full-circle-cell-strata-exact-v3` identify the stronger execution graph.

Replay must reject:

- deletion of a one-root signed predicate;
- mutation of its radicand;
- mutation of either original polynomial;
- reassignment of a physical endpoint event to a conjugate root;
- loss of original-equation recheck evidence; and
- a v3 source payload presented under the v4 tag.

Rational-only fixtures run through the same v4 source grammar and must retain
their verdicts, topology, orientation, and exact replay.

## Performance contract

This stage first closes the complete exact path. It may reason away duplicate
source parsing and repeated root reconstruction, but it does not add a
speculative index.

Required structural counts:

- one coordinate decomposition per unique boundary vertex;
- one passage predicate per unique vertex and quarter chart;
- one algebraic-root solve per square-free factor through the existing cache;
- one original-equation recheck per candidate root; and
- no boundary extraction or partition rebuild inside the cell loop.

After Task 13F is terminal, the measured unordered-feature product remains the
next optimization target.

## Acceptance

The stage is complete only when:

1. the diagonal post-capsule positive control no longer returns unresolved;
2. its exact partition reconstructs byte-identically;
3. at least one norm root is retained without a physical endpoint event,
   proving conjugate filtering is exercised;
4. rational full-circle fixtures remain green through v4;
5. every one-root mutation fails replay;
6. the Task 13F fourth link and following circle advance to a proved verdict;
7. the complete adaptive suite, Ruff, strict mypy, strict MkDocs, and diff
   checks pass; and
8. the theorem, limitation, evidence, Held comparison, and performance impact
   are recorded in developer MkDocs before publication.
