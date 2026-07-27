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
    live endpoints, and classifies constant clearance. Nonparallel and
    externally limited segment/segment cells, the unified segment-site
    traversal, neck evidence, proposal sampling, and the Python binding are
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

| Dimension | Held–Pfeiffer 2025 | Exact-certified Phase 1 |
| --- | --- | --- |
| Pocket geometry | Segments and circular arcs; simply connected; machinability assumed after an `r + ε` transformation | Target: polygonal pockets with holes and exact `C_r`/`M_r`; current MAT subset is incomplete; circular boundaries are not supported |
| MAT backend | Vroni/ArcVroni used end-to-end | Exact CGAL point graph plus bounded P–S and parallel S–S production slices with stable provenance; unified graph incomplete |
| Engagement limit | Analytic circle construction followed by bisection until `θmax − 0.001 <= θ <= θmax` radians | Exact rational chord surrogate and event-exact segment/full-circle certification; integration incomplete |
| Candidate spacing | Bisection along the middle curve | Finite candidate lattice; no monotonic-feasibility assumption |
| Machined state | Ordered contour of prior machining disks; transition sweeps omitted from that contour model | Exact stock mutation and exact full-sweep coverage, ordered certify-before-deplete |
| Bottlenecks | Graph search, width, and heuristic cap reduction on first passage | Exact neck evidence, separating cut, typed passage state, and bound cap decision; planned |
| Validation | Dense engagement sampling for result plots; complete path and runtime experiments | Fresh replay, mutation gates, artifact identity, and exact residual proof; planned |
| End-to-end evidence | Complete paths, figures, path-length gains, and 3–100 ms path-generation timings excluding Voronoi construction | Tasks 10–16 remain frozen; no complete path, Fig. 5 reproduction, or performance result yet |

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
normalized ring vertices + open segments + exact radius²
                         │
                         ▼
segment-Delaunay graph with stable caller provenance
                         │
                         ▼
degeneracy-removal Voronoi adaptor
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
4. its endpoints must equal the adaptor vertices and the derived feature
   bounds exactly;
5. each live limiter owner must occur in that endpoint's feature provenance.

The last rule is deliberately bounded. A third point or segment site that
truncates a parallel S–S cell before a source feature endpoint raises
`UnsupportedSegmentSegmentLimiterError`; its algebraic limiter equation is
not approximated or silently treated as a feature bound.

Owner classification must precede feature-bound point equality. The
adversarial native fixture places an unrelated point at `q = (5, 3/2)`.
Along `p(t) = (t, 3/2)`, equality with the S–S clearance is
`(t - 5)² = 9/4`, so the live cell is truncated at `t = 7/2` before the
upper feature endpoint at `t = 4`. This fixture must reach
`UnsupportedSegmentSegmentLimiterError`; reporting only that the live point
differs from the feature bound loses the actionable failure classification.

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

Complete graph records are identical after reversing both segment endpoint
orders. Horizontal, vertical, and rational diagonal charts have independent
native goldens. Nonparallel S–S supports remain outside this completed slice:
their raw exact primitives now have a canonical quadratic-field
angle-bisector chart, while algebraic feature and limiter binding and graph
emission remain incomplete.

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

!!! warning "A nonparallel generator pair has two branch identities"

    The exact producer fixture pairs a horizontal support of squared normal
    length `1` with a 45-degree support of squared normal length `2`. The
    degeneracy-removal adaptor exposes two distinct bounded S–S halfedges for
    the same ordered generator pair, one on each
    `Q(sqrt(2))` angle-bisector branch; they meet at the feature transition
    owned by `lower-right`. A generator-pair-only dual ID would collide, and
    the parallel slice's unique-halfedge rule would reject valid topology.
    The nonparallel implementation must bind a canonical branch identity to
    every original dual before graph emission.

    The same fixture also proves that `halfedge->dual()` is not invariably the
    complete geometric payload after degeneracy removal. One normalized
    halfedge reports adaptor endpoint `(8, 1 + sqrt(2))`, while its underlying
    exact primal segment ends at `(8, 1 - sqrt(2))`; the other branch retains
    adaptor/primal equality. Nonparallel emission must resolve the adaptor's
    underlying edge-chain semantics or traverse raw SDG primitives and rebuild
    normalized adjacency. Substituting adaptor endpoints into the representative
    primal curve is forbidden.

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
| `segment_site_parameterization.*` | exact primitive domains and coordinate functions |
| `segment_site_clipping.*` | exact domain and clearance roots; maximal admissible cells |
| `segment_site_provenance.*` | stable site/root provenance and endpoint evidence |
| `segment_site_endpoint_binding.*` | live parabola and bounded parallel S–S endpoint ownership |
| `segment_site_mat.*` | graph orchestration and canonical records |
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
| General segment/segment cells | pending | nonparallel and externally limited cells fail loud; no production claim |
| Unified degeneracy-removal traversal | pending | no production claim |
| Degeneracy-normalized feature CSR | pending | no production claim |
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
reversal invariance. These gates establish the bounded P–S and feature-owned
parallel S–S slices. They do not establish nonparallel or externally limited
S–S cells, degeneracy-normalized feature CSR, or the final unified traversal.

The final Task 9 boundary also requires:

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
