# Exact Segment-Site Medial Axis

The exact segment-site medial axis (MAT) is the geometric spine of the
exact-certified adaptive-clearing pipeline. It is not a display skeleton and
it is not a sampled approximation. Its native certificate must preserve the
topology, defining sites, clipping events, and admissible-center boundary
needed by every downstream traversal and engagement decision.

!!! warning "Implementation maturity"

    Task 9 is in progress. The point/point production graph is locally green
    and independently reviewed. The bounded point/segment slice now takes its
    algebraic cell endpoints from live Voronoi halfedge owners, including
    point and segment limiters, and clips rational point/segment parabolas at
    exact true-radius roots. A bounded parallel segment/segment production
    slice now derives its exact feature-overlap interval, binds feature-owned
    live endpoints, binds rational external point and open-segment limiters
    through exact algebraic events, proves open-segment interior ownership, and
    classifies constant clearance. The same path now binds both external
    open-segment endpoints of the rectangle's central S–S cell. Raw
    nonparallel segment/segment branches now have exact source-bound charts,
    strict feature intervals, rational quadratic clearance, and
    polygon-with-holes clipping. The rectangle's lower-left nonparallel cell
    now binds its shared feature endpoint and opposite open-segment limiter
    through the same live adaptor contract. One rectangle production path now
    traverses a single four-segment adaptor build, emits all five interior S–S
    cells, rejects all eight incident exterior P–S rays, and canonicalizes the
    six exact nodes from live degree-three feature incidence. Repeat,
    four-segment reversal, and exact `r² = 0`, `1`, `4`, and `5` states are
    native-gated.
    A bounded normalized nonparallel slice now reconstructs a two-branch
    adaptor cell from its exact raw S–S and P–S chain, preserves separate
    primitive and parent ownership, and survives exact radius clipping.
    Canonical graph nodes now project once into deterministic `int64`
    node-site CSR; the six-node rectangle golden locks offsets and flattened
    stable feature order. A canonical MAT site catalog now consumes Task 3's
    already-normalized rings, emits exact point/open-segment sites with typed
    `(kind, ring, feature)` provenance, retains segment endpoint ownership,
    and maps catalog-bound stable node-site identities to numeric CSR without
    parsing those identities. An additive catalog-fed source path now builds a
    valid segment-Delaunay graph, resolves live generators through one exact
    geometry index, and proves a complete catalog/live-site bijection. A typed
    consuming transform now swaps that validated SDG into the
    degeneracy-removal Voronoi adaptor without a duplicate graph build. The
    same catalog now projects into exact rational point sources and directed
    open-segment supports for the algebraic parameterization layer. Parallel
    and nonparallel S–S endpoint binders now accept the transferred exact
    geometry index directly; a catalog-owned rectangle gate proves that this
    indexed path emits the same algebraic parameters and complete provenance
    as the earlier validated identity table. The same index now binds both
    retained P–S parabolas of an exact concave L-shaped pocket, including all
    four segment-limiter events, with repeat and canonical input-reversal
    invariance. The same L-pocket adaptor is now inventoried at the normalized
    node boundary: eleven finite vertices have three exact incident features
    and one has the four-feature identity `{P2, P4, S2, S3}`. Three-feature
    node bytes remain unchanged; higher-degree nodes use a distinct,
    length-framed identity. A catalog-fed rectangle orchestration now carries
    that owner through the complete bounded graph:
    one adaptor traversal emits five S–S edges and six canonical nodes,
    rejects eight authenticated incident P–S rays, and reproduces the
    independently validated rectangle records at exact `r² = 0`, `1`, `4`,
    and `5` after a structurally validated certificate-namespace projection.
    Repeated construction and canonical input reversal are record-identical
    inside the catalog namespace.
    General arbitrary-pocket traversal, externally limited and composite
    segment/segment cells, endpoint-feature CSR, public numeric-table
    integration, neck evidence, proposal sampling, and the Python binding are
    not complete.

    Do not treat the current native spike APIs as the final public MAT API.
    The maturity table below is the claim boundary.

## Why a separate exact MAT exists

The repository's straight-skeleton trochoidal example is a separate
algorithm. It cannot establish the segment Voronoi facts required by the
certified adaptive design:

- a point/open-segment bisector is a parabola, not a straight-skeleton edge;
- tool-center admissibility is the exact set `C_r = D ⊖ B_r`, not an offset
  inferred from samples;
- graph nodes and adjacency must survive degeneracy without coordinate
  matching;
- necks require exact local clearance minima and a separating graph cut;
- downstream replay must bind exact topology and provenance, not a sampled
  polyline.

The MAT component therefore wraps CGAL's segment-Delaunay/Voronoi machinery
with explicit exact parameterization, clipping, identity, and evidence
contracts.

## Relation to Held and Pfeiffer (2025)

[Held and Pfeiffer's 2025 MATHSM extension][held-pfeiffer-2025] is the direct
algorithmic reference. It is a complete, fast toolpath generator with
published paths and runtime measurements. Phase 1 is not yet a competing
end-to-end result; it is a stricter certification architecture under
construction.

| Dimension | Comparative status | Held–Pfeiffer 2025 | Exact-certified Phase 1 |
| --- | --- | --- | --- |
| Pocket geometry | incomplete | Segments and circular arcs; simply connected; machinability assumed after an `r + ε` transformation | Raw MAT primitives now clip exactly against polygonal domains with holes and exact radius clearance; the rectangle `C_r` graph is unified, but general arbitrary-pocket traversal is incomplete; circular boundaries are not supported |
| MAT backend | stronger exact contract on bounded fixtures; incomplete end-to-end | Vroni/ArcVroni used end-to-end | Exact CGAL point graph plus bounded P–S and parallel/nonparallel S–S cells; the catalog-fed four-segment rectangle path now carries Task 3 identities through exact rational sources, one indexed SDG build, swap-owned adaptor traversal, indexed endpoint binding, exact radius clipping, five edges, six feature-triple nodes, and explicit rejection of eight incident P–S rays; the same indexed owner contract binds two retained concave L-pocket P–S parabolas and their four exact segment-limiter events, while a complete normalized-node inventory proves eleven three-feature vertices and one four-feature `{P2, P4, S2, S3}` vertex without coordinate merging; repeat and canonical input symmetry are record-invariant, but the complete L-pocket graph and arbitrary-pocket traversal remain incomplete |
| Engagement limit | stronger contract; incomplete integration | Analytic circle construction followed by bisection until `θmax − 0.001 <= θ <= θmax` radians | Exact rational chord surrogate and event-exact segment/full-circle certification; integration incomplete |
| Candidate spacing | incomplete | Bisection along the middle curve | Finite candidate lattice; no monotonic-feasibility assumption |
| Machined state | stronger contract | Ordered contour of prior machining disks; transition sweeps omitted from that contour model | Exact stock mutation and exact full-sweep coverage, ordered certify-before-deplete |
| Bottlenecks | incomplete | Graph search, width, and heuristic cap reduction on first passage | Exact neck evidence, separating cut, typed passage state, and bound cap decision; planned |
| Validation | incomplete | Dense engagement sampling for result plots; complete path and runtime experiments | Fresh replay, mutation gates, artifact identity, and exact residual proof; planned |
| End-to-end evidence | weaker | Complete paths, figures, path-length gains, and 3–100 ms path-generation timings excluding Voronoi construction | Tasks 10–16 remain frozen; no complete path, Fig. 5 reproduction, or performance result yet |

### Performance claim boundary

Held and Pfeiffer report complete path-generation time on their Fig. 5 pocket
of about 100 ms at `theta_max = 20 degrees`, 30 ms at `40 degrees`, 9 ms at
`80 degrees`, and 3–6 ms for larger limits. Their scope excludes I/O and
Voronoi construction; the latter took less than 1 ms on their unspecified
"standard desktop PC."

The repository's pre-certification
`trochoidal_mat_toolpath_circular` remains fast. At commit `0365b7a`, a warm
Release build on an Apple M1 Max measured one untimed warm-up followed by seven
in-process generations per pocket:

| Pocket | Operations | Median generation time |
| --- | ---: | ---: |
| square | 271 | 4.586 ms |
| irregular | 402 | 6.763 ms |
| L-shape | 300 | 4.885 ms |
| star | 611 | 16.519 ms |
| kite | 335 | 5.997 ms |
| dumbbell | 732 | 16.363 ms |
| island | 389 | 6.281 ms |

!!! warning "These millisecond ranges do not establish runtime parity"

    The repository timings use different pockets and the measured generator
    does not regulate engagement. It is the separate straight-skeleton
    algorithm described above, not the exact segment-site MAT and adaptive
    traversal under construction. The overlapping numerical range shows that
    fast path assembly survived; it does not show that the certified system is
    as fast as Held and Pfeiffer.

The catalog-fed bounded rectangle path removes three structural costs before
timing: live generators resolve through an exact `O(log S)` geometry index
instead of an `O(k S)` caller-site scan, the validated SDG owner cannot be
copied or moved implicitly, and adaptor construction swaps rather than copies
the complete SDG. The indexed resolver now reaches every parallel and
nonparallel S–S endpoint-owner lookup in the adopted rectangle traversal.
The catalog P–S dispatcher likewise resolves every defining and limiter source
by sorted stable ID and every live owner through the geometry index; it does
not reconstruct a linear generator table. These are architectural performance
safeguards, not a new planner benchmark or evidence for arbitrary-pocket
complexity.

Two other clocks must remain separate. Five warm local Release executions of
the current native Task 9 algebraic fixture suite took 1.79–1.82 s (1.82 s
median), but that is a contract gate rather than one pocket generation. The
existing exact stock replay
measured 87 s for the kite after a 7x end-to-end speedup; that run certifies and
depletes every operation and is not comparable to Held and Pfeiffer's planner
clock. Its cost history and the eliminated incremental-Boolean pathology are
recorded in the
[engagement-audit performance analysis](superpowers/state/sp1-gate-c-analysis.md).

The first defensible performance comparison therefore comes after Tasks 13–15
on the provenance-locked Fig. 5 fixture. It must report, on one machine and
build:

1. segment-site MAT construction;
2. candidate selection and traversal;
3. exact certification plus stock replay;
4. total wall-clock, path length, and emitted operation count.

Until then the accurate status is: **generation mechanics on par in magnitude;
certified end-to-end performance incomplete and currently much slower where it
has been measured.**

The stricter contract changes several seemingly small choices:

- the paper's `0.001`-radian bisection band is a numerical acceptance rule;
  Phase 1 certifies an exact rational surrogate;
- the paper perturbs the machinability transformation by `ε` to avoid
  boundary equality cases; Phase 1 treats equality, tangency, plateaus, and
  cocircularity as explicit producer states;
- the paper's first/second bottleneck passage becomes a typed, replayable
  state transition;
- the contour-aware optimization becomes authoritative stock and coverage
  lineage rather than unbound generator state.

!!! note "What to retain from the reference implementation"

    The paper exploits traversal order: its machined contour is updated by a
    backwards scan, and bottlenecks are found by a graph-linear search.
    Exact certification should retain that structural efficiency. Exactness
    does not justify repeated full arrangement builds or all-pairs rematching.

The paper also identifies an open problem that remains relevant here:
approximating a certified circular/linear path by smoother higher-order
primitives does not automatically preserve its engagement bound. Phase 1
therefore certifies its declared primitive grammar and does not infer a C²
certificate from Hausdorff closeness.

## Proof-boundary dataflow

The target Task 9 dataflow is:

```text
Task 3 canonical rings + exact radius²
                         │
                         ▼
typed point/open-segment site catalog
                         │
                         ▼
exact geometry index + one segment insertion pass
                         │
                         ▼
segment-Delaunay graph with stable caller provenance
                         │
                         ▼
swap-owned degeneracy-removal Voronoi adaptor
                         │
                         ▼
LINE / RAY / SEGMENT / PARABOLA dual ownership
                         │
                         ▼
exact parameterization + exact endpoint binders
                         │
                         ▼
split at D-boundary roots and clearance² = radius² roots
                         │
                         ▼
maximal connected components in C_r
                         │
                         ▼
canonical nodes, edges, feature CSR, and MAT certificate
                  ┌──────┴──────┐
                  ▼             ▼
        exact neck evidence   proposal-only samples
```

Each arrow is a typed transformation. A contract breach raises a named
exception; it does not select a fallback representation.

!!! danger "Do not collapse exact and proposal layers"

    `sample_centers`, `sample_clearance`, chord error, and sagitta exist for
    proposing machining circles and rendering the MAT. They never decide
    topology, clipping, node identity, neck ownership, or certificate
    validity.

## Decision data versus reporting data

The MAT follows the repository-wide
[exact-kernel discipline](exactness.md):

| Concern | Decision representation | Reporting representation |
| --- | --- | --- |
| Site coordinates | exact rational or exact quadratic-field values | `float64` adapters |
| Primitive domain | exact rational/algebraic parameter interval | sampled parameter |
| Domain membership | exact polygon predicates and algebraic signs | rendered XY |
| Tool fit | exact `distance_to_defining_site² >= radius²` | clearance estimate |
| Node identity | stable site/root/curve identity | displayed coordinates |
| Adjacency | exact Voronoi incidence | polyline connectivity |
| Neck class | exact algebraic width comparison | displayed width |

No `double`, `to_double`, handle address, iterator ordinal, or sampled
coordinate may cross from the reporting column into the decision column.

## CGAL topology ownership

### Segment-Delaunay is the substrate, not the output graph

`CGAL::Segment_Delaunay_graph_2` owns the defining sites and the Delaunay
features. Its `primal(edge)` construction supplies the corresponding
Voronoi primitive as a `Line_2`, `Ray_2`, `Segment_2`, or
`Parabola_segment_2`.

Walking `finite_edges_begin()` directly is useful for bounded compile spikes,
but it is not the final topology contract. Raw Delaunay faces and edges still
contain degeneracies that the Voronoi view must normalize.

The production traversal target is:

```cpp
using SegmentSiteAdaptationPolicy2 =
    CGAL::Segment_Delaunay_graph_degeneracy_removal_policy_2<
        SegmentSiteDelaunay2>;

using SegmentSiteVoronoi2 = CGAL::Voronoi_diagram_2<
    SegmentSiteDelaunay2,
    SegmentSiteAdaptationTraits2,
    SegmentSiteAdaptationPolicy2>;
```

The policy rejects degenerate Delaunay edges and faces before the adaptor
exposes Voronoi halfedges and vertices.

### Halfedge access is sided

A Voronoi halfedge exposes the two defining generators through `up()` and
`down()`. Its endpoint limiters are conditional:

- call `left()` only when `has_source()` is true;
- call `right()` only when `has_target()` is true.

An unbounded halfedge can lack one or both endpoints. Accessing an absent
side violates a CGAL precondition and can terminate a Release build. Endpoint
binding therefore proves boundedness before reading the corresponding
limiter.

### Stable site provenance is mandatory

CGAL handles, addresses, and iteration order are process-local artifacts.
They never enter a certificate.

The unified segment-site graph must use
`Segment_Delaunay_graph_storage_traits_with_info_2` with explicit conversion
and merge functors, or an exact post-map whose bijection is proved by native
tests. Stable caller information is recovered through
`storage_site().info()`.

The current point-only slice uses an exact post-map because every point site
is rational and caller-supplied. Segment insertion is more demanding:
derived endpoint point-sites and open-segment sites must preserve their
distinct feature identities through CGAL's storage transformations.

## Exact primitive parameterization

The clipping pipeline does not interrogate rendered coordinates. Each dual
primitive receives a parameter domain and exact coordinate functions.

### Linear primitives

For a point/point bisector:

```text
x(t) = x₀ + x₁ t
y(t) = y₀ + y₁ t
```

The primitive domain distinguishes:

- line: `t ∈ (-∞, +∞)`;
- ray: one exact finite bound;
- segment: two exact finite bounds.

The live CGAL primitive determines which domain applies. Exact Voronoi
vertices orient finite bounds; a limiter predicate orients a ray.

### Point/segment parabolas

A nonincident point/open-segment dual is parameterized from the exact focus
and supporting line. Coordinate functions live in one quadratic extension:

```text
x(t) = x_r(t) + √d x_s(t)
y(t) = y_r(t) + √d y_s(t)
```

The source and target of the open segment bound the feature domain. Live
Voronoi endpoint limiters then bind the actual halfedge endpoints. The
accepted binders distinguish point limiters from segment limiters and reject
ambiguous, unbounded, or mismatched producer states.

#### The supporting-line chart is canonical

A supporting line is geometrically projective: `(a, b, c)` and
`λ(a, b, c)` describe the same line for every nonzero `λ`. The parabola
parameter chart is not scale-free, however. Its tangent basis is `(-b, a)`,
so scaling the line by `λ` maps the same geometric point from `t` to
`t / λ`. Any algebraic-root identity formed in that chart would therefore
depend on endpoint spacing or endpoint order unless the line is canonical.

`canonical_open_segment_source` fixes the chart at the exact-source
boundary:

1. derive `(a, b, c)` from the two exact rational endpoints;
2. reject an empty stable site identity;
3. reject coincident endpoints;
4. clear all rational denominators and divide by the integer gcd;
5. choose the sign whose first nonzero coefficient is positive.

The result is one primitive integer triple for every rational representation
of the same orientation-independent supporting line. For example, both
vertical endpoint orders reduce to `(1, 0, 0)`, while
`(-100, 0) → (100, 0)` reduces from the raw `(0, 200, 0)` to
`(0, 1, 0)`. Production source and limiter construction both pass through
this factory.

This normalization changes coordinates, not geometry. In the bounded
point-limiter fixture, the raw endpoint-length chart produced the factor
`1 + 8t`; the canonical chart produces `1 + t` for the same Voronoi
endpoint. The native golden records the latter as
`AlgebraicRootIdV1({1, 1}, 0)`. The segment-limiter fixture likewise binds
its predeclared canonical root only after both supporting lines pass through
the same boundary.

!!! note "Canonical source is a value-type invariant"

    `MatExactOpenSegmentSource2` is immutable and is not an aggregate. Its
    constructor is private; `canonical_open_segment_source` is the only
    construction boundary used by production and native fixtures. Empty
    identities and zero line normals therefore fail before parameterization,
    clipping, or endpoint binding can observe an invalid chart.

### Parallel segment/segment cells

The first S–S production slice covers two distinct parallel open-segment
supports. Let their canonical supporting lines be:

```text
ℓ₁(x, y) = ax + by + c₁
ℓ₂(x, y) = λ(ax + by) + c₂,  λ > 0
```

Canonical line signs make `λ` positive. After rescaling the second line to the
first normal, the equal-distance locus between the supports is the rational
middle line:

```text
ax + by + (c₁ + c₂/λ)/2 = 0
```

The chart chooses the exact point on this line closest to the coordinate
origin and the tangent `(b, -a)`, with its first nonzero component positive.
Because both supporting lines are canonical, swapping generator arguments or
reversing either segment's endpoints cannot rescale or reverse `t`.

An S–S line is valid only while the perpendicular projection lies in both
open-segment features. Each source endpoint has one exact tangent parameter;
intersecting the two open parameter intervals gives the positive-length
feature cell. An empty or point-only overlap raises
`EmptyParallelSegmentFeatureDomainError`.

The live producer contract is stronger than reconstructing that interval:

1. `up()` and `down()` must be the two expected segment generators;
2. both endpoint sides must exist before reading `left()` or `right()`;
3. `SDG::primal()` must be an exact `Segment_2`;
4. its endpoints must equal the adaptor vertices exactly;
5. each endpoint must equal either a derived feature bound or one exact root
   of its external limiter equation;
6. each live limiter owner must occur in that endpoint's retained provenance.

The feature-only binder remains deliberately strict: a third site without an
exact source record raises `UnsupportedSegmentSegmentLimiterError`. The
production overload now dispatches a live external owner to either its exact
point-source record or its exact open-segment record. Neither owner kind is
treated as a feature bound, and no coordinate-only reconstruction is accepted.

Owner classification precedes feature-bound equality. The adversarial native
fixture places an unrelated point at `q = (5, 3/2)`. Along
`p(t) = (t, 3/2)`, equality with the constant S–S clearance is:

```text
(t - 5)² - 9/4 = (2t - 7)(2t - 13)/4 = 0
```

The canonical source-feature chart remains `t ∈ [0,4]`. Exact root filtering
retains only `t = 7/2`, so the live Voronoi cell is `t ∈ [0,7/2]`; `13/2` is
outside the feature domain. The emitted limiter endpoint carries the
`AlgebraicRootIdV1` for `2t - 7`, `external-limiter`, and the source-factor
multiplicity. The chart domain is never rewritten to impersonate the narrower
live cell. `clip_bounded_linear_clearance_components` independently requires
an increasing, provenance-bearing live interval contained in the feature
domain.

The native failure matrix separately requires:

- a missing point-source record to raise
  `UnsupportedSegmentSegmentLimiterError`;
- duplicate stable limiter identity to raise
  `AmbiguousParallelSegmentPointLimiterError`;
- a non-rational limiter record to raise
  `NonRationalParallelSegmentPointLimiterError` before live-site comparison;
- a rational record at the wrong coordinate or a non-canonical chart to raise
  `MismatchedLiveSegmentSegmentBridgeError`.

An external open-segment limiter uses a different ownership proof. For the
vertical open segment

```text
q(u) = (5, 1) + u(0, 1),  0 < u < 1
```

the squared distance from `p(t)` to its support is `(t - 5)²`, so equality with
the S–S clearance produces the same exact factorization:

```text
(t - 5)² - 9/4 = (2t - 7)(2t - 13)/4 = 0
```

Support equality alone is insufficient: it also describes the infinite line
through the limiter. For every candidate root, the binder proves the
perpendicular foot is strictly inside the live open segment with

```text
0 < ⟨p(t) - q₀, q₁ - q₀⟩ < ‖q₁ - q₀‖².
```

The feature interval again rejects `13/2`, while `t = 7/2` has foot
`(5, 3/2)` and passes both strict inequalities. The endpoint carries the
`AlgebraicRootIdV1` for `2t - 7`, `external-segment-limiter`, and
`parallel-segment-limiter/equation-factor-multiplicity/1`. Reversing all three
source segments leaves the complete graph record unchanged.

The open-segment failure matrix requires a missing source record to raise
`UnsupportedSegmentSegmentLimiterError`, a duplicate stable identity to raise
`AmbiguousParallelSegmentOpenLimiterError`, and a record on a different
support to raise `MismatchedLiveSegmentSegmentBridgeError`.

The rectangle's central cell exercises two independent open-segment limiters
in one live halfedge. Its bottom and top generators give

```text
p(t) = (t, 0),  t ∈ [-4, 4],  δ² = 4.
```

The left support `x = -4` contributes
`(t + 4)² - 4 = (t + 2)(t + 6)` and the right support `x = 4`
contributes `(t - 4)² - 4 = (t - 2)(t - 6)`. Exact feature filtering retains
only `t = -2` and `t = 2`. The production cell is therefore `[-2,2]`, with
separate `left-segment` and `right-segment` event/root provenance. All four
rectangle segments are inserted as real generator sites; no synthetic point
site creates either endpoint.

!!! note "Adaptor topology and endpoint proof are separate contracts"

    The exact rectangle characterization inserts four open boundary segments
    through the single degeneracy-removal adaptor. It yields five bounded S–S
    segments for the interior MAT and eight incident unbounded P–S rays for
    exterior transitions. The adaptor therefore already owns the correct
    interior adjacency; production code must bind and certify its endpoints,
    not rebuild a skeleton by coordinate matching. In particular, the central
    cell's source-feature interval is `[-4,4]`, while its adaptor-owned live
    interval is `[-2,2]`; feature overlap alone overestimates the MAT cell.
    External point or segment records narrow an existing live cell only after
    exact owner identity, equation roots, and open-feature membership all
    agree.

For the native production fixture, the lower segment is
`[-2, 6] × {0}` and the upper segment is `[0, 4] × {3}`. The canonical chart
and feature cell are:

```text
p(t) = (t, 3/2)
t ∈ [0, 4]
```

Clearance is constant:

```text
δ² = distance(p(t), ℓ₁)² = distance(p(t), ℓ₂)² = 9/4
```

`parallel_segment_clearance_boundary` reconstructs the canonical chart,
rejects any rescaled or displaced raw chart, and requires both defining-site
distances to be exactly equal. A raw affine chart therefore cannot bypass
the geometry or root-identity invariant. The production gate proves all three
classifications:

| Exact `r²` | Sign of `δ² - r²` | Production result |
| --- | --- | --- |
| `2` | positive | complete feature cell retained |
| `9/4` | zero | complete closed boundary plateau retained |
| `5/2` | negative | no S–S component emitted |

The external-point and external-open-segment fixtures prove the same three
classifications over their complete live interval `[0,7/2]`.

Complete graph records are identical after reversing both segment endpoint
orders. Horizontal, vertical, and rational diagonal charts have independent
native goldens. Nonparallel S–S support requires a different production
contract because one normalized adaptor halfedge can span more than one raw
exact primitive; the bounded composite slice below exercises that contract.

### A nonparallel rectangle corner has an exact external limiter

The lower-left rectangle cell is defined by the bottom and left open segments.
Their interior signed-distance bisector has canonical chart

```text
p(t) = (-4, -2) + t(1, 1).
```

The two source features permit `t ∈ [0,4]`. The shared point site
`lower-left` owns `t = 0`, while the opposite top segment limits the other
side. Squared-distance equality to its support is

```text
t² = (4 - t)²
16 - 8t = 0
t = 2.
```

The live adaptor cell is therefore `[0,2]`. The target endpoint carries the
simple `AlgebraicRootIdV1` for `t - 2`, `top-segment`, and
`nonparallel-segment-limiter/equation-factor-multiplicity/1`. The binder
separately proves the live up/down generators, exact primal/adaptor equality,
canonical signed branch, exact root membership in `[0,4]`, equality to the
live endpoint, and strict projection into the top open segment.

For a general quadratic-field chart, the segment event is formed from the two
exact squared support distances. Candidate roots come from a rational
zero-set polynomial. When both rational and radical components are present,
that polynomial is the field norm; every candidate is then checked against
the original field equation to reject conjugate artifacts. Mixed equations
carry `norm-factor-multiplicity` provenance; this bounded stage does not
reinterpret that elimination multiplicity as source-equation multiplicity.

!!! warning "Do not norm an equation that is already rational"

    If the radical component is zero, the field norm squares the rational
    equation. The geometric root remains correct, but a simple event factor is
    falsely reported with multiplicity two. `field_zero_set_polynomial`
    therefore preserves an existing rational or pure-radical equation and
    tags its `equation-factor-multiplicity`; it uses and explicitly tags the
    norm only for a genuinely mixed field equation. The two multiplicity
    meanings are never silently conflated.

The bounded clipper requires two increasing, provenance-bearing live
endpoints contained in the source-feature interval. At `r² = 0` it retains
`[0,2]`; at `r² = 1` it retains `[1,2]`; at `r² = 5` it emits no edge.
Negative radius fails before graph construction. Complete records are
unchanged after reversing all four rectangle segments. Missing, duplicate,
and support-mismatched limiter records, reversed live bounds, and ownerless
live bounds each raise distinct named errors.

### One adaptor build yields a unified rectangle graph

The isolated center and lower-left cells are now consumed by one production
traversal over the four rectangle segments. The producer inserts the four
open segments once, proves an eight-feature source bijection for the four
segments and their four endpoint point-sites, and constructs one
degeneracy-removal Voronoi adaptor. Each undirected adaptor edge has one
canonical orientation, selected by the ordered `up()`/`down()` stable IDs.

The traversal authenticates two disjoint classes:

- five bounded S–S halfedges become exact `LINE` graph edges;
- eight unbounded P–S rays are proved incident to their parent segment and
  rejected from the interior MAT.

No coordinate rematching reconstructs adjacency. Parallel versus nonparallel
dispatch is the exact determinant of the two supporting-line normals. Every
S–S endpoint is bound to its live `left()` or `right()` owner, clipped by the
appropriate exact clearance polynomial, and emitted with its original dual,
generator pair, parent pair, root, owner, and multiplicity provenance.

An emitted edge-local endpoint initially has identity
`stable_endpoint_node_identity_v1(dual_id, endpoint)`. At an adaptor vertex,
the producer circulates all incident halfedges and unions their `up()` and
`down()` feature IDs. The rectangle contract requires exactly three distinct
features. Their ordered triple becomes the canonical
`stable_voronoi_node_identity_v1(...)`; point features map back to both
incident parent segments when parent ownership is assembled.

This yields five edges and six nodes:

| Node class | Count | Degree | Exact feature identity |
| --- | ---: | ---: | --- |
| Corner terminal | 4 | 1 | corner point plus its two incident segments |
| Interior junction | 2 | 3 | bottom/top plus left or right segment |

The complete graph record is identical on repeat and after reversing every
source segment. Exact clearance changes the graph without changing the
underlying adaptor:

| Exact `r²` | One-dimensional graph |
| --- | --- |
| `0` | five edges, six feature-triple nodes |
| `1` | five edges, two feature-triple junctions, four exact clearance-root terminals |
| `4` | one central plateau edge and its two junction nodes |
| `5` | empty |

!!! warning "Certified set components are not automatically graph edges"

    The exact clipper returns maximal **closed admissible sets**. At `r² = 4`,
    each corner branch survives at its junction as a valid isolated point.
    Emitting that zero-dimensional component as an edge creates a self-loop
    and corrupts degree. The one-dimensional graph boundary therefore compares
    its two algebraic bounds exactly: increasing intervals become edges,
    equal bounds remain certified-set events but are not emitted as edges, and
    reversed bounds raise `InvalidSegmentSiteGraphComponentError`. The retained
    central plateau edge still owns both shared junction nodes.

!!! note "Canonical node identity comes from incidence, not position"

    A degree-three node belongs to the exact union of features around the live
    adaptor vertex. Edge-local endpoint identities are aliases into that
    feature triple. Matching rendered coordinates would create a second,
    weaker identity system and would not prove which generator sites or parent
    segments own the junction.

### Nonparallel S–S branches have exact source-bound charts

!!! note "Why this slice is unusually difficult"

    A production nonparallel S–S edge must keep four contracts aligned:
    the mathematical angle-bisector branch, the raw SDG primal, the
    degeneracy-normalized Voronoi adjacency, and the stable certificate
    identity. No one CGAL object owns all four. The raw primal owns one exact
    curve primitive; the adaptor owns normalized connectivity and may traverse
    rejected primitives; the caller owns stable generator provenance; and the
    canonical chart owns algebraic parameter and branch identity.

    These layers may agree geometrically while disagreeing structurally.
    Rendered coordinates can therefore look correct even when branch IDs
    collide, endpoint ownership is wrong, or a normalized halfedge spans more
    geometry than its representative primal. The implementation proves each
    correspondence separately rather than treating any one handle, coordinate,
    or iteration order as universal truth.

Let the canonical supporting lines, ordered by stable segment identity, be

```text
ℓ₁(x, y) = a₁x + b₁y + c₁
ℓ₂(x, y) = a₂x + b₂y + c₂
```

and let `nᵢ = aᵢ² + bᵢ²`. Their two Euclidean angle bisectors are represented
without division or floating-point normalization by

```text
n₂ ℓ₁(x, y) + s sqrt(n₁ n₂) ℓ₂(x, y) = 0,  s ∈ {-1, +1}.
```

`nonparallel_segment_bisector_parameterization` constructs both lines in the
exact kernel and requires both endpoints of the live raw `SDG::primal()`
segment to lie on exactly one branch. That exact incidence test binds `s`.
Endpoint order never selects the branch. A zero-length primitive, parallel
supports, or a primitive on neither branch raises a distinct named error.

The chart origin is the rational intersection of `ℓ₁` and `ℓ₂`. Before
orientation normalization, a tangent to branch `s` is

```text
(
  n₂b₁ + s sqrt(n₁n₂)b₂,
 -n₂a₁ - s sqrt(n₁n₂)a₂
).
```

The exact sign of the first nonzero quadratic-field coordinate fixes its
orientation. Stable source ordering fixes the branch equation; exact tangent
sign fixes the parameter direction. These are deliberately separate
invariants. Reversing either source segment, swapping factory arguments, or
reversing the live primitive therefore produces byte-identical chart
coefficients and the same branch identity.

If `n₁n₂` is a rational square, the factory folds every radical coefficient
into its rational coefficient and records canonical radicand `1`. A rational
bisector cannot retain a cosmetically quadratic representation. The native
goldens cover both branches in `Q(sqrt(2))` and a perpendicular rational
fixture.

For a source endpoint `e` on support `ℓᵢ`, let
`vᵢ = (bᵢ, -aᵢ)` be the support tangent and let the bisector chart be
`p(t) = o + ut`. The feature-transition parameter is not obtained by
substituting `e` into the bisector: `e` is generally not on that line.
Instead, the perpendicular foot reaches `e` exactly when

```text
vᵢ · (p(t) - e) = 0,

               vᵢ · (e - o)
t_endpoint = -----------------.
                    vᵢ · u
```

The numerator is rational and the denominator is in `Q(sqrt(d))`.
`nonparallel_segment_tangent_parameter` divides by multiplying the exact
quadratic conjugate, returning the canonical rational/radical coefficient
pair. A foreign segment identity, an endpoint off its declared support, or a
singular projection fails loud.

Each open segment contributes the interval between its two endpoint
parameters. `intersect_nonparallel_segment_feature_domains` sorts both pairs
by exact field sign and takes their strict positive-length intersection.
Equal active bounds union all endpoint provenance. A point-only intersection
is a feature transition, not an S–S cell, and fails loud at this boundary. In
the two-branch native fixture, the bounded raw cells are:

```text
s = +1: t ∈ [-22 - 11sqrt(2), -1 - sqrt(2)]
s = -1: t ∈ [ -8 +  4sqrt(2),  1 - sqrt(2)]
```

Field bounds are compared by the exact sign of their difference. To enter the
shared algebraic endpoint and node-identity substrate, `a + bsqrt(d)` is
embedded as the root selected by the sign of `b` from

```text
(x - a)² - b²d = 0.
```

The primitive integer factor plus its ordered real-root ordinal is identity;
an isolating interval is not. For example, `-22 - 11sqrt(2)` is root `0` of
`x² + 44x + 242`, while `-8 + 4sqrt(2)` is root `1` of
`x² + 16x + 32`. Rational fields collapse before embedding.

#### Quadratic-field coordinates have rational quadratic clearance

The global primitive contract permits a clearance polynomial of degree at
most four. The canonical nonparallel S–S chart has stronger structure. Let
`Δ = a₁b₂ - a₂b₁`. Substituting the unnormalized chart tangent `u` into the
second ordered support gives

```text
ℓ₂(u) = -n₂Δ.
```

The radical terms cancel exactly. Reversing `u` changes only this value's
sign, so both defining-site squared distances are

```text
distance(p(t), ℓ₁)²
  = distance(p(t), ℓ₂)²
  = n₂Δ²t².
```

The true-radius boundary is therefore the rational quadratic

```text
n₂Δ²t² - r² = 0.
```

`nonparallel_segment_clearance_boundary` does not assume the closed form
blindly. It substitutes the stored field chart into both normalized support
distance equations, requires both radical polynomials to vanish, requires
the two rational polynomials to be identical, and checks their quadratic
coefficient against `n₂Δ²` before isolating roots. A malformed source/chart
pair fails before it can emit clearance evidence. The native scaled fixture
proves that this is not an accidental unit-coefficient result:
`n₂Δ² = 125` reduces `125t² - 5` to primitive polynomial
`25t² - 1`.

`maximal_nonparallel_segment_clearance_components` embeds the two feature
bounds into the shared algebraic kernel, inserts only clearance roots inside
that interval, and emits every maximal retained subcell. The native gate
exercises full retention, one-root clipping, and complete rejection on the
upper branch, plus one-root clipping on the lower branch. Physical endpoint
IDs remain the feature-bound provenance; a coincident clearance event is
unioned rather than replacing its owner.

#### Polygon-with-holes clipping stays in one quadratic field

For an affine nonparallel chart

```text
p(t) = o + ut,  o,u ∈ Q(sqrt(d))²
```

and a rational polygon edge `a + es`, the exact intersection solves

```text
ut - es = a - o.
```

`nonparallel_segment_domain_roots` evaluates the two-by-two determinant in
`Q(sqrt(d))`. Both the chart parameter `t` and edge parameter `s` remain
quadratic-field values. Exact field comparisons enforce `s ∈ [0,1]` and the
strict feature interval before `t` is embedded in the shared algebraic root
kernel. Parallel-disjoint edges contribute no root. Collinear supports are
intersected as exact parameter intervals: a disjoint interval contributes
nothing, a one-point contact contributes one provenance-bearing root, and a
positive-length overlap raises `OverlappingDomainBoundaryError` because an
isolated event list cannot represent an interval intersection.

Open cells are classified by a winding predicate whose sidedness and
orientation determinants use the complete quadratic-field coordinates. The
upper `Q(sqrt(2))` fixture clips against an outer rectangle and an interior
rectangular hole. With `r² = 9`, it emits exactly:

```text
t ∈ [-10, -8]  with outer-to-hole endpoint provenance
t ∈ [ -6, -4]  with hole-to-outer endpoint provenance
```

Every boundary endpoint carries both its stable ring-edge ID and its
reconstructible algebraic-root identity. A mismatched chart/feature radicand
raises `MismatchedNonparallelSegmentFeatureDomainError`.

!!! warning "Exact arithmetic does not repair an incomplete expression"

    A winding determinant that retained radical `y` terms but multiplied
    them by only the rational parts of `x` was exactly evaluated and still
    geometrically wrong. Rational fixtures could not expose the omission.
    Domain predicates over `Q(sqrt(d))` must form the complete field
    determinant first and only then take its exact sign.

!!! warning "Zero-set canonicalization can reverse conjugate selection"

    CGAL's bivariate curve cache canonicalizes a polynomial up to a nonzero
    unit. This is valid for root isolation but can reverse the sign returned
    for a source polynomial such as `-4 - t`. Conjugate filtering for
    `R + sqrt(d)S = 0` therefore uses
    `exact_polynomial_sign_at_2`, which restores the original leading-unit
    sign and factor-multiplicity parity. The regression fixture has focus
    `(sqrt(2), 1)`, directrix `y = -1`, and a rectangular hole; it requires
    both exact retained intervals
    `[sqrt(2) - 4, sqrt(2) - 1]` and `[2, 2sqrt(2)]`.

#### Normalized composite adjacency follows parent ownership

!!! warning "A rejected raw edge is not necessarily disposable"

    The exact producer fixture pairs a horizontal support of squared normal
    length `1` with a 45-degree support of squared normal length `2`. The
    degeneracy-removal adaptor exposes two distinct bounded S–S halfedges for
    the same ordered generator pair, one on each
    `Q(sqrt(2))` angle-bisector branch; they meet at the feature transition
    owned by `lower-right`. A generator-pair-only dual ID would collide, and
    the parallel slice's unique-halfedge rule would reject valid topology.
    The producer therefore binds a canonical signed-branch identity to every
    normalized dual before graph emission.

    The same fixture also proves that `halfedge->dual()` is not invariably the
    complete geometric payload after degeneracy removal. One normalized
    halfedge reports adaptor endpoint `(8, 1 + sqrt(2))`, while its underlying
    exact primal segment ends at `(8, 1 - sqrt(2))`; the other branch retains
    adaptor/primal equality. Substituting the adaptor endpoint into that
    representative primal curve would fabricate geometry.

The bounded production slice rebuilds the normalized adjacency from raw exact
SDG primitives. A raw segment feature maps to its own parent segment; a raw
endpoint point maps to every incident parent segment. Each raw feature pair is
then classified before its reject bit is interpreted:

1. a pair whose two features map to the same parent is a self-transition and
   is discarded;
2. a rejected P–S pair that maps to the two distinct target parents is retained
   as part of the normalized branch;
3. ownership that maps simultaneously to self and distinct parents, or to more
   than one distinct parent pair, raises `AmbiguousCompositeSiteOwnerError`.

For the native fixture, the parent sites are:

```text
lower-segment:    (-20, 0) -> (8, 0)
diagonal-segment: (5, -4)  -> (20, 11)
```

The rejected `lower-right`/`lower-segment` LINE is a self-transition and is
discarded. The rejected `lower-right`/`diagonal-segment` PARABOLA belongs to
the distinct parent pair and is retained. The normalized result is therefore
one single-LINE branch and one composite LINE–PARABOLA branch. Exactly one
self-transition is reported as rejected.

The graph keeps three identity layers because they answer different proof
questions:

| Record | Meaning |
| --- | --- |
| `generator_site_ids` | raw features that define this exact conic piece |
| `parent_site_ids` | normalized open segments whose MAT adjacency this piece serves |
| `original_dual_id` | stable signed-bisector branch shared by every piece of one normalized halfedge |

The S–S and P–S pieces initially have different endpoint identities. Their
transition is joined only after their exact CGAL endpoint values compare
equal; the existing certificate node identities are then aliased and their
raw-feature, parent-site, and event provenance is unioned. Coordinates prove
incidence but never become a replacement node ID.

At `r² = 0`, the native gate requires three edges and four nodes: two LINE
pieces, one PARABOLA piece, two degree-two transition nodes, and two terminal
nodes. At `r² = 1`, exact clearance clipping still emits three edges but splits
the graph into five nodes, with one surviving shared transition and four
terminals. Negative `r²` fails at the composite producer boundary before SDG
construction. Repeated construction and reversal of both source segments
must produce byte-for-byte equal graph records, including all three identity
layers. Missing cells, ambiguous ownership, unsupported primitive variants,
or a non-unique chain fail through named exceptions.

!!! note "Current bounded claim"

    This gate proves one two-cell, one-transition normalized nonparallel
    topology with one admissible component per raw piece. The separate
    lower-left rectangle gate proves one externally segment-limited
    nonparallel raw cell. Together they still do not prove unbounded cells,
    externally limited composite chains, multiple retained transitions,
    arbitrary composite-chain length, or the unified pocket traversal.

The following `Parabola_segment_2` APIs are visualization paths and are
forbidden from proof decisions:

- `side_of_parabola`;
- `compute_k`;
- `generate_points`;
- streamed parabola output.

### Algebraic roots stay algebraic

Substituting a primitive into a polygon boundary or clearance equation can
produce roots beyond the coordinate field. A quartic root remains an
`AlgebraicRootIdV1`; it is not coerced into a `MatKernel::Point_2`.

Root identity is:

1. a primitive, square-free integer polynomial with positive leading
   coefficient;
2. the ordered real-root ordinal of that polynomial.

An isolating interval is validation evidence, not identity. Refining it must
not change certificate bytes.

### Coordinate roots and event roots are different evidence

An endpoint can carry more than one valid algebraic-root identity.

- The **coordinate-root identity** answers where the endpoint lies on its
  primitive. Exact graph binding always adds this identity.
- An **event-root identity** answers which exact equation created the split.
  Domain clipping and clearance clipping preserve these identities.

For example, the rational endpoints `t = ±1/4` of a clearance split have
linear coordinate factors, but the event that creates both endpoints is:

```text
16t² - 1 = 0
```

The clearance event therefore retains the quadratic factor plus root ordinal
in addition to each endpoint's coordinate-root identity. These identities
must not be collapsed merely because they reconstruct the same exact
parameter. Native goldens compare the complete canonical provenance vector,
so an omitted or spurious event identity fails the gate.

The original clearance polynomial remains authoritative for sign
classification. Root isolation uses its canonical square-free part, and each
event identity is recorded with its ordinal **before** primitive-domain
filtering. This both satisfies CGAL's square-free solver precondition and
prevents an out-of-domain root from renumbering retained certificate events.

## Exact clipping to the admissible-center domain

Let `D` be the pocket polygon with holes and `r` the tool radius. The MAT
graph is clipped to:

```text
C_r = D ⊖ B_r
```

On a valid Voronoi dual, either defining site realizes the nearest-site
distance. Therefore a center lies in `C_r` exactly when both conditions
hold:

1. the point lies in `D`;
2. `distance_to_defining_site² >= r²`.

The implementation:

1. derives exact roots where the primitive crosses the boundary of `D`;
2. derives exact roots of `clearance² - r²`;
3. orders and coalesces coincident roots;
4. classifies every open parameter cell with exact predicates;
5. emits every maximal admissible connected component.

An unbounded primitive is never discarded merely because its original domain
is unbounded. A line or ray can have one or more bounded components inside a
polygon with holes.

### True-radius point/segment clearance

For a rational point-site focus and canonical open-segment directrix,
`source_parabola_clearance_boundary` substitutes the exact quadratic
parameterization into:

```text
distance(point_site, p(t))² - r²
```

It clears denominators to one primitive integer polynomial of degree at most
four, computes its canonical square-free part, and isolates all real roots
before cell-domain filtering. Root ordinals therefore belong to the complete
canonical curve chart and cannot change when a bounded Voronoi cell retains
only a subset of them.

The native production-path fixture uses focus `(0, 2)`, directrix `y = 0`,
and `r² = 9/4`. Its canonical chart is:

```text
x(t) = -t
y(t) = 1 + t²/4
```

The clearance boundary becomes:

```text
t⁴ + 8t² - 20 = 0
```

The exact roots `-√2` and `+√2` split the bounded source feature into two
retained graph components. Each emitted algebraic endpoint reconstructs the
same `AlgebraicRootIdV1` carried in its clearance provenance.

Negative `r²` raises `NegativeClearanceRadiusSquaredError`. A source focus
with a nonzero quadratic-field coefficient raises
`NonRationalParabolaClearanceError`; normalized pocket point sites are
rational, while general quadratic-field clearance remains outside this
bounded slice. The native graph path accepts exact radius squared, but the
public Python tool-radius binding remains pending.

### Constant-clearance cells

Constant clearance is classified before square-free normalization:

- positive: retain the complete primitive domain;
- negative: reject it;
- zero: retain the complete cell as a closed `C_r`-boundary plateau.

Fabricating roots for an identically zero equation destroys plateau ownership
and is forbidden.

## Endpoint and node provenance

Each emitted edge binds:

- stable original-dual identity;
- connected-component index;
- ordered generator-site identities;
- source and target exact node identities;
- independent endpoint flags for original Voronoi vertex, `D`-boundary
  clip, and clearance root;
- every incident domain feature at tangencies and degenerate vertices.

Coincident events are additive. A point can simultaneously be an original
Voronoi vertex, a polygon-boundary intersection, and a clearance root; one
flag must not overwrite another.

Owned feature endpoints enter domain clipping directly. For the parallel
S–S chart over the rectangle `[0, 4] × [-1, 4]`, the polygon intersections
at `t = 0` and `t = 4` coincide exactly with the feature bounds. The canonical
endpoint records retain the feature owner, algebraic-root identity, and
`D-outer/edge-3` or `D-outer/edge-1`, respectively. A post-clipping endpoint
replacement would discard the polygon evidence and is therefore forbidden;
clipping coalesces provenance by exact parameter equality instead.
`clip_owned_linear_clearance_components` rejects unbounded endpoints,
parameters that differ from the primitive domain, and empty ownership
provenance as separate named failures.

At a degeneracy-normalized node, provenance is the union over every incident
halfedge. A single underlying Delaunay face triple can omit valid incident
generators and is not sufficient evidence.

!!! warning "A halfedge owner is not the complete node provenance"

    In the native parallel S–S fixture, both segment endpoint pairs meet the
    same two transition locations, but the adaptor exposes only the upper
    endpoint point-sites through `left()` and `right()`. Treating those two
    limiter handles as exhaustive would silently discard the coincident lower
    feature evidence. Node CSR must union all incident halfedges and exact
    coincident feature bounds.

### Node-site CSR is a projection of canonical topology

`node_site_csr` consumes the sorted `MatExactGraphNode2` records after adaptor
incidence has already been normalized. It validates that node identities are
strictly ordered and nonempty, and that every node's stable feature IDs are
nonempty, unique, and strictly ordered. It then performs one flattening pass:

```text
node_site_offsets = [0, 3, 6, 9, 12, 15, 18]
node_site_ids     = 18 stable feature identities
```

Those literals are native-gated for the six-node rectangle graph. An empty
graph has the valid CSR sentinel `[0]`. Offsets are `int64`; overflow raises
`MatNodeSiteCsrOverflowError`, while malformed node and site ordering raise
`NonCanonicalMatGraphNodeOrderError` and
`NonCanonicalMatNodeSiteOrderError`.

!!! warning "Do not rebuild incidence during CSR export"

    The difficult generator union belongs to the degeneracy-removal adaptor
    traversal. CSR export must serialize those canonical node records, not
    circulate CGAL halfedges again. Recomputing incidence during export would
    create a second topology algorithm and could make certificate bytes
    disagree with the edge graph they encode.

### Canonical site catalog and numeric identity

`canonical_mat_site_catalog` consumes `CanonicalReachInput2`; it does not
orient, rotate, sort, or otherwise normalize rings again. Task 3 remains the
sole owner of:

- exact binary64 injection;
- outer/hole orientation;
- minimal ring rotation;
- identity sorting of holes;
- canonical feature ordinals.

The catalog enumerates one exact point site and one exact open-segment site per
canonical ring feature. `ReachKernelPoint` and the segment-Delaunay
`MatTraits::Point_2` are compile-time-proven identical exact types, so this
handoff has no coordinate conversion or reporting seam. Every open-segment
record also stores the two canonical point-site identities that own its
endpoints.

The public `site_provenance[S,3]` convention is:

| Column | Encoding |
| --- | --- |
| kind | `0` point site; `1` open-segment site |
| ring | `0` outer ring; `1 + canonical_hole_ordinal` for holes |
| feature | canonical vertex ordinal for a point; canonical directed-edge start ordinal for an open segment |

Catalog rows and versioned stable identities share lexicographic
`(kind, ring, feature)` order. The stable identity is an opaque exact record;
consumers must not parse it. `numeric_node_site_csr` first projects canonical
node incidence through `node_site_csr`, then resolves each stable identity by
catalog lookup. An identity outside the input's catalog raises
`UnknownMatSiteIdentityError`.

!!! warning "One normalization owner, one numeric translation"

    Rebuilding ring orientation or hole order inside Task 9 would make MAT
    site IDs diverge from the `C_r` digest even if the coordinates happened to
    match. Conversely, assigning integers during CGAL traversal would make
    IDs depend on adaptor iteration. Canonical input owns geometry identity;
    the site catalog owns the single stable-to-numeric translation.

Radius is deliberately absent from site provenance because changing tool
radius does not change the pocket's generator features. The Task 3
center-domain digest and final MAT certificate bind the exact radius
independently.

The catalog and numeric mapper are native-gated on rectangle symmetry, exact
site geometry, canonical hole reordering, hand-derived table rows, empty CSR,
and unknown-site rejection. The catalog-fed rectangle graph now consumes
these site records directly. Its node records therefore carry catalog
identities suitable for the numeric mapper; the direct graph-to-public-table
boundary remains part of the pending endpoint-feature CSR slice.

### Exact rational graph sources

The live adaptor and the algebraic clipping layer consume the same catalog in
different exact forms. `CanonicalMatRationalSources2::build` projects each
input point coordinate from `CORE::Expr` to `CORE::BigRat` without decimal
text, `double`, or approximation:

1. obtain a candidate with `BigRatValue()`;
2. inject that candidate back into `CORE::Expr`;
3. require exact equality with the original coordinate.

`BigRatValue()` alone is not sufficient because an algebraic expression can
yield a rational approximation. Failure of candidate extraction or exact
round-trip equality raises `NonRationalCanonicalMatCoordinateError`.

The native non-dyadic oracle locks:

```text
binary64(0.1) = 3602879701896397 / 2^55
binary64(0.2) = 3602879701896397 / 2^54
```

It independently rejects `sqrt(2)`, proving that the seam accepts exact
rational input coordinates rather than arbitrary algebraic values.

Each rational segment record resolves the catalog's ordered source and target
point identities, verifies their exact geometry, and constructs one primitive
integer support line through `canonical_open_segment_source`. Missing endpoint
records and exact endpoint disagreement raise
`MissingCanonicalMatEndpointSourceError` and
`MismatchedCanonicalMatEndpointGeometryError`, respectively.

!!! note "Ordered ownership, unordered live geometry"

    The geometry index treats a segment's endpoint pair as unordered because
    reversing a CGAL segment does not change its generator. Rational source
    projection preserves canonical ring direction because feature intervals,
    limiter ownership, and endpoint provenance do depend on source versus
    target. These contracts are complementary, not interchangeable.

The rectangle gate locks exact point rows, the four support lines `y + 2`,
`x - 4`, `y - 2`, and `x + 4`, endpoint wraparound, non-dyadic binary64
fractions, radical rejection, and complete input-symmetry invariance.

### Catalog-fed segment-Delaunay source

`CanonicalMatDelaunaySource2::build` is the additive input path from the
canonical site catalog into CGAL. It builds two exact sorted indexes:

- point sites keyed by exact `compare_xy`;
- open-segment sites keyed by their exact unordered endpoint pair.

Point and segment key spaces stay distinct because a ring vertex is expected
to coincide with two segment endpoints. Duplicate geometry within either key
space is ambiguous caller ownership and raises
`DuplicateCanonicalMatSiteGeometryError`.

Only the canonical open segments are inserted. CGAL materializes their
endpoint point-sites, so explicit endpoint reinsertion would be redundant.
After insertion, the factory walks every finite Delaunay vertex, resolves it
in `O(log S)` exact lookup, rejects duplicate resolutions, and requires the
matched stable-ID count to equal both the index and catalog cardinalities.
Invalid or incomplete live graphs raise `CanonicalMatDelaunayBijectionError`;
an unknown exact generator raises `UnknownCanonicalMatSiteGeometryError`.

!!! note "Constructed graph ownership is immobile"

    `CanonicalMatDelaunaySource2` is neither copyable nor movable. Older CGAL
    graph types can satisfy an rvalue operation by copying, so a nominal
    move-only wrapper would not prove that a full graph build was transferred.
    C++17 guaranteed copy elision initializes the factory result directly.
    `CanonicalMatVoronoiSource2::build(std::move(source))` is the sole consuming
    boundary and invokes CGAL's explicit swap construction.

The source gate proves exact rectangle cardinality, rotation/reversal
invariance, two-hole insertion, exact catalog/live-site equality, unknown-site
rejection, and duplicate point geometry rejection. It establishes the
catalog-to-SDG boundary independently. The catalog-graph gate described below
now composes that source with the full five-edge rectangle traversal while the
earlier bounded path remains an independent complete-record comparator.

### No-copy Voronoi-adaptor transfer

CGAL's ordinary `SegmentSiteVoronoi2(delaunay)` constructor copies the complete
segment-Delaunay graph. Passing `true` as its `swap_dg` argument transfers the
graph instead. `CanonicalMatVoronoiSource2` makes that easily missed
performance choice structural:

1. reject an already-empty source with
   `ConsumedCanonicalMatDelaunaySourceError`;
2. construct the degeneracy-removal adaptor with `swap_dg=true`;
3. require the source SDG to contain no finite vertices after transfer;
4. revalidate the adaptor and every live generator through the transferred
   exact geometry index;
5. require the post-transfer generator count to equal the index cardinality.

The resulting adaptor owner is also noncopyable and nonmovable. A failed
source-empty, adaptor-validity, or generator-bijection check raises
`InvalidCanonicalMatVoronoiTransferError`.

The rectangle transfer gate independently recovers 13 canonical raw dual
pairs from the new owner: five segment/segment pairs and the eight incident
point/segment pairs that the graph layer rejects. The complete signature is
invariant under ring rotation/reversal, and consuming one source twice fails
loudly.

!!! warning "No-copy ownership still needs graph-level evidence"

    A successful swap proves only the catalog-to-adaptor execution boundary.
    The bounded rectangle adoption is separately gated on complete graph
    records; it was not inferred from source cardinality or raw dual counts.
    General arbitrary-pocket traversal must pass the same graph-level boundary
    before no-copy ownership can support a broader production claim.

### Indexed endpoint binding preserves the exact contract

The parallel and nonparallel S–S endpoint algorithms used to accept only a
`GeneratorSite2` table. Every live up/down and left/right owner lookup scanned
that table by exact geometry. This remained correct, but feeding the canonical
adaptor through that interface would discard the geometry index immediately
before the most lookup-intensive part of traversal.

Both binder families now expose overloads that accept
`CanonicalMatSiteGeometryIndex2`. The original table overloads remain
available. All overloads delegate to one resolver-parameterized implementation,
so the following checks are shared rather than reimplemented:

- exact up/down generator agreement;
- exact primal/adaptor endpoint agreement;
- canonical parallel chart or signed nonparallel branch;
- exact feature interval and live owner;
- external point/open-segment limiter equations;
- algebraic root identity, multiplicity, and endpoint provenance.

!!! warning "Do not rebuild a linear identity table at the binder boundary"

    The transferred geometry index authenticates the relation between every
    live CGAL site and the canonical input catalog. Converting it back into a
    vector would preserve answers on small fixtures but reintroduce
    `O(k S)` rematching and a second identity representation. Pass the index
    through traversal and keep geometry-to-identity resolution in `O(log S)`.

The native equivalence gate runs on the catalog-owned, swap-transferred
rectangle adaptor. It binds the externally segment-limited central parallel
cell and the lower-left nonparallel cell through both overload families, then
requires exact algebraic parameter equality and identical complete provenance
vectors. This establishes the endpoint-binding seam independently; the
catalog-graph gate then composes those indexed binders into the complete
bounded traversal.

#### Concave P–S cells retain the catalog index

The exact L-shaped ring

```text
(0,6) ─ (2,6)
  │       │
  │       (2,2) ─ (6,2)
  │                 │
(0,0) ─────────── (6,0)
```

has one reflex point feature `P3`. The live degeneracy-removal adaptor contains
two bounded nonincident parabolas:

| Defining sites | Live endpoint owners |
| --- | --- |
| `P3–S0` | `S5`, `S2` |
| `P3–S5` | `S3`, `S0` |

`bind_point_segment_cell_endpoints` accepts the canonical rational sources and
`CanonicalMatSiteGeometryIndex2` directly. It derives the directed
open-segment feature interval in the canonical source-parabola chart, proves
the live up/down defining sites, proves the adaptor endpoints belong to the
exact CGAL parabola, and dispatches each left/right owner:

- a source-segment endpoint becomes its exact rational feature bound;
- an external point site delegates to the exact point-limiter equation;
- an external open segment delegates to the exact segment-limiter equation
  with strict feature ownership.

The returned endpoints are strictly ordered algebraic parameters with root,
owner, factor-kind, and multiplicity provenance. The L-shape gate compares
each indexed result with two direct calls to the independently gated
segment-limiter binder, then requires the two-cell signature to survive repeat
and canonical input reversal. Unknown rational sources, an incident
point/own-segment pair, and valid sources that disagree with the live halfedge
raise `UnknownCanonicalMatParabolaSourceError`,
`IncidentCanonicalMatParabolaSourceError`, and
`MismatchedLiveParabolaBridgeError`, respectively.

!!! warning "A site ID does not certify an endpoint event"

    The geometry index proves which canonical site owns a live adaptor handle;
    it does not by itself prove which algebraic event or point on the primitive
    that handle represents. A P–S endpoint must also agree with the canonical
    source chart, the finite feature interval, the adaptor vertex, and the
    exact CGAL parabola endpoint. Keep identity lookup and event certification
    separate, then require both.

This closes the indexed P–S endpoint seam for two retained concave cells. It
does not yet emit or canonicalize the complete L-shaped graph; that traversal
must also handle its P–P ray, bounded S–S cells, incident P–S rejection, domain
and radius clipping, and shared-node provenance.

#### Normalized nodes use the complete feature union

The L-pocket adaptor has 12 finite normalized vertices. Eleven are incident to
three exact canonical features. The remaining vertex is incident to four:

```text
{P2, P4, S2, S3}
```

This is not a coordinate coincidence and it is not an invalid face. It is a
real degeneracy-normalized Voronoi event. A producer that derives node identity
from one raw Delaunay face triple would split or reject this valid topology.
The gate therefore circulates every incident adaptor halfedge, unions the
catalog identities of its `up()` and `down()` generators, and only then derives
the node identity.

`stable_normalized_voronoi_node_identity_v1` requires at least three nonempty,
strictly ordered, unique feature IDs. For exactly three IDs it delegates to
`stable_voronoi_node_identity_v1`, preserving every established node byte. For
four or more IDs it uses the distinct `normalized-voronoi-node/v1` namespace
and length-frames every feature ID. The framing makes delimiter-bearing site
identities unambiguous; the separate namespace prevents a higher-degree record
from masquerading as the older three-feature encoding.

!!! warning "Do not infer normalized topology from a face triple"

    Degeneracy removal changes the producer boundary. The authoritative node
    identity is the sorted union of all exact features incident to the live
    adaptor vertex, not one raw face, one representative primal, or a rendered
    coordinate. Count adaptor degree independently from feature cardinality:
    both happen to be four at this L-pocket event, but they are different
    certified facts.

The native gate requires exactly eleven degree-three and one degree-four node,
unique IDs, repeat identity, canonical ring-reversal identity, byte
compatibility for every three-feature node, strict malformed-input rejection,
and collision-safe framing. This certifies the normalized node-provenance seam;
it does not yet certify the L-pocket's complete edge set or radius-clipped
graph.

### Catalog-fed traversal closes the bounded graph seam

`canonical_rectangle_mat_graph` is the additive orchestration path from one
`CanonicalReachInput2` and exact squared radius to a complete
`MatExactGraph2`. It does not rebuild a friendly generator table or construct a
second SDG. The path:

1. builds the canonical site catalog and exact rational feature sources;
2. validates the exact one-ring rectangle contract;
3. inserts the four directed open segments once and transfers the validated
   SDG into one degeneracy-removal adaptor with `swap_dg=true`;
4. traverses each undirected adaptor dual once, using the exact geometry index
   for all generator and endpoint-owner resolution;
5. dispatches the five S–S cells by exact support determinant, binds their
   parallel or nonparallel endpoints, clips them at exact domain and clearance
   roots, and emits only one-dimensional components;
6. authenticates and rejects the eight incident P–S rays;
7. canonicalizes live adaptor vertices from exact feature incidence, retains
   clearance-created terminals under their exact dual-plus-root identity, and
   validates that no emitted edge collapses.

The acceptance gate first requires complete record identity on repeat and
after canonical ring reversal. It then compares the catalog graph with the
independently validated bounded path at exact `r² = 0`, `1`, `4`, and `5`.
The comparison includes primitive kind, generator and parent sites, original
dual and component identity, algebraic endpoint-root identity, complete
endpoint provenance, node feature/parent/provenance sets and degree, rejected
transition count, and matched-site count.

!!! warning "Project certificate namespaces structurally"

    Opaque catalog site IDs propagate into Voronoi-node IDs, clipped
    endpoint-node IDs, original-dual IDs, component IDs, and domain-boundary
    provenance. Two graphs can therefore certify the same physical features
    while their raw record bytes correctly differ by namespace. Prove the
    source-site mapping first, validate every derived ID against its own
    canonical fields, preserve the dual branch kind, and only then project
    exact identity tokens and their slash-delimited descendants. A live
    Voronoi vertex recomputes from three incident features; a clearance-created
    terminal recomputes from its unique incident dual plus algebraic root.
    Never coerce one node class into the other's identity rule. Coordinate
    matching, count comparison, and unconstrained string replacement are not
    certificate equivalence.

This closes catalog-to-production-graph adoption for the bounded rectangle,
not for arbitrary pockets. The earlier bounded path remains unchanged as an
independent oracle; considering removal or consolidation is a separate,
explicitly authorized step after broader traversal evidence exists.

## Failure anatomy: four cocircular point sites

The square fixture

```text
(-1, -1)  (1, -1)  (1, 1)  (-1, 1)
```

exposes why raw finite-edge traversal is not a topology contract.

All four sites are cocircular. The exact Voronoi graph has one degree-four
node at `(0, 0)`. In the raw point-only walk:

1. CGAL suppresses the zero-length central dual;
2. four unbounded rays remain;
3. the two finite Delaunay faces produce different stable generator triples;
4. both triples have the exact same rational circumcenter;
5. triple-derived IDs therefore split one geometric node into two degree-two
   nodes.

Checking only `source_node_id == target_node_id` cannot detect the defect:
there is no emitted central edge to inspect.

The current point-only boundary registers exact rational Voronoi locations.
Repeated visits with the same stable triple are idempotent; a different
triple at the same exact location raises
`DegeneratePointSiteTopologyError`. This fails loudly instead of returning
malformed topology.

!!! note "Why not silently merge the nodes?"

    A point-only rational circumcenter is enough to detect this case, but it
    is not the canonical cross-curve identity required by the complete
    segment-site graph. Silent coordinate-based merging would establish a
    second, narrower identity system. The unified Voronoi adaptor and exact
    node/provenance contract own the eventual normalization.

## Neck evidence is a graph certificate

Necks are not simply low-clearance samples. The completed graph will carry
four exact `NeckEvidenceV1` variants:

| Variant | Exact event |
| --- | --- |
| `STRICT_EDGE` | relative-interior derivative changes negative to positive |
| `CLEARANCE_ENDPOINT` | one-sided minimum at a `C_r` clip |
| `SHARED_VERTEX` | equal incident minima merged at one exact node |
| `PLATEAU` | maximal connected constant-clearance interval or subgraph |

Every record binds its defining sites, exact squared width, deterministic
owner, and separating-cut partition. Evidence is retained only when removing
the point, vertex, or contracted plateau separates at least two nonempty
traversal sides.

This is why neck extraction belongs after exact graph normalization and before
proposal sampling.

## Proposal-only sampling

Sampling is deliberately downstream:

- lines use exact barycentric parameter proposals;
- parabolas bisect exact parameter intervals;
- reported chord and sagitta errors decide only whether to refine the
  proposal representation;
- exhausting the refinement depth raises `ConicSamplingLimitError`.

Refinement may change the number and position of samples. It must not change:

- node or edge identity;
- connected components or cycle rank;
- endpoint provenance;
- neck ownership;
- MAT certificate bytes.

## File ownership

| File | Responsibility |
| --- | --- |
| `segment_site_catalog.*` | Task 3-derived exact sites, typed numeric provenance, and stable catalog lookup |
| `segment_site_delaunay.*` | exact geometry index, one-pass segment insertion, and catalog/live-site bijection |
| `segment_site_voronoi.*` | explicit no-copy SDG consumption and post-transfer adaptor validation |
| `segment_site_rational_sources.*` | exact rational coordinates, directed endpoint ownership, and normalized support lines |
| `segment_site_parameterization.*` | exact primitive domains and coordinate functions |
| `segment_site_clipping.*` | exact domain and clearance roots; maximal admissible cells |
| `segment_site_provenance.*` | stable site/root provenance, normalized feature-union node identity, and endpoint evidence |
| `segment_site_endpoint_binding.*` | live parabola and bounded parallel/nonparallel S–S endpoint ownership with linear-table and canonical-index overloads sharing one exact algorithm |
| `segment_site_graph_emission.*` | exact one-dimensional component projection and graph record emission |
| `segment_site_graph_csr.*` | validated deterministic projection of canonical node-site provenance |
| `segment_site_catalog_graph.*` | catalog-fed bounded rectangle orchestration, exact adaptor traversal, and canonical node assembly |
| `segment_site_mat.*` | canonical graph records, general graph primitives, and independently validated bounded comparators |
| `segment_site_neck.*` | exact local minima and separating cuts |
| `segment_site_mat_sampling.*` | proposal-only samples |
| `medial_axis_2.cpp` | nanobind adapter only |

File length is not a quality gate. New neck, sampling, binding, or algebraic
responsibilities belong in their named files only when their invariants,
consumers, or evolution differ from graph orchestration.

## Maturity and claim matrix

| Capability | State | Claim |
| --- | --- | --- |
| Point/point linear parameterization | implemented | exact |
| Point/point domain and radius clipping | implemented | exact |
| Point-graph deterministic records | implemented and reviewed | exact |
| Cocircular point topology | guarded | exact collision detection and named fail-loud boundary; normalization pending |
| Point/segment source parameterization | implemented and native-gated | exact canonical supporting-line chart |
| Point/segment endpoint ownership | wired and native-gated | exact for bounded point and segment limiter fixtures |
| Algebraic point/segment cell bounds | implemented | exact |
| True-radius point/segment clearance | implemented and native-gated | exact for rational point-site sources; public binding pending |
| Parallel segment/segment feature domains | implemented and native-gated | exact for positive-length overlap with feature-owned live endpoints |
| External point-limited parallel S–S cell | implemented and native-gated | exact rational quadratic limiter equation, algebraic live endpoint, complete event provenance, and constant-clearance clipping |
| External open-segment-limited parallel S–S cell | implemented and native-gated | exact one- and two-sided rational quadratic support equations, strict open-feature ownership, algebraic live endpoints, complete event provenance, and constant-clearance clipping |
| Nonparallel segment/segment raw cells | implemented and native-gated | exact source-bound branches, feature intervals, algebraic endpoints, true-radius clipping, and polygon-with-holes clipping |
| External open-segment-limited nonparallel S–S cell | implemented and native-gated | exact for the lower-left rectangle branch: mixed feature/segment ownership, field-equation filtering, live interval clipping, and event multiplicity |
| Normalized nonparallel S–S composite fixture | implemented and native-gated | exact two-branch reconstruction; distinct raw-generator, parent-site, and normalized-dual identity; exact transition aliasing and radius clipping |
| Unified rectangle segment-site graph | implemented and native-gated | one adaptor build; five S–S edges, six canonical nodes, eight authenticated rejected incident P–S rays; exact `r² = 0`, `1`, `4`, and `5` topology; complete-record reversal invariance |
| General segment/segment cells | pending | unbounded, arbitrary externally limited, multi-transition, and arbitrary-length composite chains have no production graph claim |
| General degeneracy-removal traversal | pending | the rectangle fixture is production-gated; arbitrary pockets and all primitive combinations have no production claim |
| Degeneracy-normalized node-site CSR | implemented and native-gated | deterministic `int64` offsets and stable feature IDs projected from canonical nodes; six-node rectangle golden and malformed-order failures |
| Canonical typed site catalog | implemented and native-gated | exact point/open-segment sites derived from Task 3 canonical rings; global ring encoding, endpoint ownership, symmetry, and hole-order invariance |
| Catalog-derived rational sources | implemented and native-gated | exact `CORE::Expr` round-trip to `BigRat`, directed endpoint ownership, primitive support lines, non-dyadic goldens, and radical rejection |
| Catalog-fed segment-Delaunay source | implemented and native-gated | exact indexed lookup, one segment insertion pass, complete generator bijection, duplicate-geometry rejection, and immobile graph ownership |
| Catalog-fed Voronoi owner | implemented and native-gated | explicit `swap_dg=true` transfer, source-empty proof, post-transfer indexed bijection, raw rectangle-pair signature, and double-consume rejection |
| Catalog-indexed S–S endpoint binding | implemented and native-gated | parallel and nonparallel catalog-owned rectangle cells match the earlier validated path in exact algebraic parameters and complete provenance; adopted by the bounded catalog traversal |
| Catalog-indexed P–S endpoint binding | implemented and native-gated | both bounded nonincident parabolas of the concave L-shaped ring match direct exact segment-limiter bindings at all four endpoints; repeat and canonical input-reversal invariant; complete L-pocket graph adoption pending |
| Normalized arbitrary-degree node identity | implemented and native-gated | L-pocket adaptor has eleven exact three-feature vertices and one exact four-feature `{P2, P4, S2, S3}` vertex; three-feature bytes remain compatible, higher-degree IDs are length-framed, and repeat/input-reversal identity is exact; complete L-pocket graph adoption pending |
| Catalog-fed unified rectangle graph | implemented and native-gated | one indexed SDG build and swap-owned adaptor; five S–S edges, six canonical nodes, eight authenticated rejected incident P–S rays; repeat/input-symmetry identity and complete-record equivalence at exact `r² = 0`, `1`, `4`, and `5` under a validated certificate-namespace projection |
| Numeric node-site catalog mapping | implemented and native-gated | catalog-bound graph identities map once to `int64` rows; unknown identities fail loud; direct adopted-graph-to-public-table gate pending |
| Endpoint-feature CSR and public numeric table | pending | no public binding or complete certificate claim |
| Exact neck evidence | pending | no production claim |
| Proposal-only sampling | pending | no production claim |
| Python fixed-tuple binding | pending | no public API claim |

## Verification gates

Every native Task 9 slice runs sequentially:

```bash
pixi run _task9-mat-build
pixi run _task9-mat-run
git diff --check
```

The canonical-source slice is bounded by native contracts for horizontal,
vertical, reversed, and rational endpoint pairs; empty identities and
coincident endpoints fail with distinct named exceptions. The same
executable checks production point- and segment-limiter provenance against
exact root goldens, recomputes each golden from the emitted algebraic
parameter, and compares complete graph records after reversing every segment
endpoint order. The parallel S–S gate additionally checks live
`Segment_2`/adaptor equality, horizontal/vertical/diagonal canonical charts,
partial and empty feature overlap, owner-evidence and rescaled-chart
mutations, exact positive/zero/negative clearance, and complete-record
reversal invariance. The external point-limiter gate additionally checks its
quadratic equality, feature-domain root filtering, live-cell provenance,
constant positive/zero/negative clearance, repeatability, reversal, and every
named malformed-record failure. The external open-segment-limiter gate checks
its support-distance equation, strict interior projection, live-cell
provenance, constant positive/zero/negative clearance, repeatability,
reversal, and named missing, duplicate, and mismatched-support failures. Its
two-sided rectangle variant additionally proves independent left/right event
roots, exact narrowing from feature interval `[-4,4]` to live interval
`[-2,2]`, eight-site bijection, and reversal of all four source segments. The
rectangle characterization independently locks five bounded interior S–S
segments and eight incident exterior P–S rays. The nonparallel raw-cell gate
checks both `Q(sqrt(2))` branches, rational-field collapse, exact endpoint
projection, strict feature-domain intersection, algebraic root identity,
rational quadratic clearance, maximal feature/clearance components, and exact
polygon-with-holes clipping. The lower-left rectangle gate additionally checks
mixed feature/segment live ownership, exact narrowing from `[0,4]` to `[0,2]`,
simple event multiplicity, strict segment projection, `r² = 1` clipping to
`[1,2]`, four-segment reversal, and named malformed-source and live-domain
failures. The same executable locks the full
quadratic-field winding determinant and sign-preserving conjugate selection
with analytic P–S and S–S endpoint goldens. The normalized composite gate
additionally proves parent-owned rejection classification, one retained raw
P–S bridge, exact S–S/P–S node aliasing, stable signed-branch identity,
positive-radius clipping, repeatability, and source-endpoint reversal
invariance. The unified rectangle gate additionally proves one-build
enumeration of all five S–S cells, exact rejection of all eight incident P–S
rays, six feature-triple node identities, the two-degree-three/four-degree-one
distribution, exact dimensional filtering at the `r² = 4` plateau, complete
record repeatability, and reversal of all four segments. These gates establish
bounded P–S, feature-owned parallel S–S, point- and open-segment-limited
parallel S–S, raw nonparallel S–S algebra, one externally
open-segment-limited nonparallel rectangle cell, the bounded two-branch
normalized fixture, and one unified rectangle graph. They do not establish
general arbitrary-pocket traversal, externally limited composite S–S cells,
endpoint-feature CSR, or public numeric-table integration. A separate
graph-emission contract retains an increasing algebraic interval, omits an
equal-bound singleton, and raises `InvalidSegmentSiteGraphComponentError` for
unbounded or reversed components. The node-site CSR gate separately locks the
empty sentinel, `int64` offsets, hand-derived rectangle feature order, and
distinct named failures for malformed node or feature ordering. The canonical
site-catalog gate independently locks Task 3 input reuse, point/open-segment
kind codes, global outer/hole ring rows, exact endpoint ownership, rectangle
rotation/reversal invariance, canonical hole reordering, numeric CSR lookup,
the empty sentinel, and named unknown-site rejection. The catalog-fed
segment-Delaunay gate separately locks one-pass segment insertion, exact
catalog/live-generator cardinality and identity, indexed lookup across two
holes, input-symmetry invariance, unknown exact geometry, duplicate exact
point geometry, and noncopyable/nonmovable graph ownership. The catalog-fed
Voronoi gate separately locks source-empty swap transfer, adaptor validity,
post-transfer indexed identity, the five-S–S/eight-P–S raw rectangle
signature, reversal invariance, owner immobility, and named double-consume
rejection. The rational-source gate separately locks exact rectangle point
coordinates, four primitive support lines, directed endpoint wraparound,
binary64 `0.1`/`0.2` fractions, algebraic-radical rejection, and full source
record invariance under input rotation/reversal. The indexed endpoint gate
then uses the catalog-owned swap-transferred adaptor to bind one parallel and
one nonparallel S–S cell through both identity interfaces. It requires exact
equality of both algebraic endpoint parameters and their complete provenance
vectors. The same indexed endpoint gate inventories the concave L-shaped
adaptor by canonical feature identity, binds its two retained P–S parabolas,
and requires each of the four endpoints to equal a direct exact
segment-limiter binding. The complete two-cell signature is repeat- and
input-reversal-invariant; unknown, incident, and live-mismatched source records
fail through distinct named errors. The catalog-graph gate separately
circulates every finite L-pocket adaptor vertex and locks the complete
normalized-node inventory: eleven degree-three feature unions and one
degree-four `{P2, P4, S2, S3}` union, all ID-unique and record-identical after
repeat and canonical input reversal. It also proves old/new three-feature byte
compatibility, length-framing collision resistance, and named rejection of
short, empty, duplicate, or unsorted feature sequences. The same gate composes
the catalog owner and indexed S–S binders through all five rectangle cells,
authenticates all eight rejected incident P–S transitions, validates every
node, dual, component, and root identity, requires complete repeat and
input-symmetry identity, and proves full semantic record equivalence with the
independent rectangle path at exact `r² = 0`, `1`, `4`, and `5`. Its namespace
projector preserves parallel versus signed nonparallel dual kind and translates
only proven identities or their exact slash-delimited descendants. A
three-edge ring and a four-edge non-orthogonal trapezoid independently raise
`UnsupportedCanonicalMatRectangleGraphError`.

The final Task 9 boundary must add the still-absent Python MAT binding test and
then requires:

```bash
pixi run format-adaptive
pixi run lint
pixi run pytest tests/adaptive/test_medial_axis.py -n auto --testmon -q
pixi run pytest tests/test_toolpath.py -n auto -q
```

The native executable currently aggregates many boolean predicates into one
process exit status. When a new gate fails, isolate the new predicate before
changing geometry or production code. Exit 1 proves a failed conjunction; it
does not identify the failed contract.

## Review checklist

Before accepting a MAT change, verify:

1. Does every topology decision remain in exact types?
2. Is the Voronoi adaptor, rather than raw Delaunay iteration, the production
   topology owner?
3. Are unbounded primitives clipped instead of discarded?
4. Are point/incident-segment transition duals rejected exactly while
   nonincident parabolas remain?
5. Are halfedge limiters accessed only on present sides?
6. Is stable caller provenance preserved through site transformations?
7. Are algebraic roots identified without approximate coordinates?
8. Are all coincident endpoint events retained?
9. Is node provenance unioned over all incident halfedges?
10. Can sampling refinement change only proposal fields?
11. Does every unsupported producer state fail with a named exception?
12. Do the native gate, focused Python tests, existing toolpath regression tests, and
    strict MkDocs build pass?

## Authoritative references

- [Exact-Kernel Discipline](exactness.md)
- [Idiomatic Exact CGAL](dev/exact-cgal-idioms.md)
- `CGAL/Segment_Delaunay_graph_2.h` from the locked CGAL 6.0.1 source
- `CGAL/Segment_Delaunay_graph_storage_traits_with_info_2.h`
- `CGAL/Segment_Delaunay_graph_adaptation_policies_2.h`
- `CGAL/Voronoi_diagram_2.h`
- `CGAL/Voronoi_diagram_2/Halfedge.h`

The locked headers are normative for the backend contract. Draw helpers and
streamed geometry are not substitutes for those APIs.

[held-pfeiffer-2025]: https://doi.org/10.14733/cadaps.2025.731-747
