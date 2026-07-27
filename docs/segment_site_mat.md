# Exact Segment-Site Medial Axis

The exact segment-site medial axis (MAT) is the geometric spine of the
exact-certified adaptive-clearing pipeline. It is not a display skeleton and
it is not a sampled approximation. Its native certificate must preserve the
topology, defining sites, clipping events, and admissible-center boundary
needed by every downstream traversal and engagement decision.

!!! warning "Implementation maturity"

    Task 9 is in progress. The point/point linear production slice is locally
    green at `b3c0f94`, but its publication review found open stable-identity,
    contract-test, and file-responsibility defects. Point/segment and
    segment/segment helpers exist at different maturity levels, but the
    unified segment-site adaptor, neck evidence, proposal sampling, and Python
    binding are not complete.

    Do not treat the current native spike APIs as the final public MAT API.
    The maturity table below is the claim boundary.

## Why a separate exact MAT exists

The older trochoidal example derives centerlines from a straight skeleton.
That path remains useful and unchanged, but it cannot establish the segment
Voronoi facts required by the certified adaptive design:

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
| Pocket geometry | Segments and circular arcs; simply connected; machinability assumed after an `r + ε` transformation | Polygonal pockets with holes; exact `C_r` and `M_r`; circular boundaries not yet supported |
| MAT backend | Vroni/ArcVroni used end-to-end | Exact CGAL segment-site graph with stable provenance; incomplete |
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

The completed Task 9 dataflow is:

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

At a degeneracy-normalized node, provenance is the union over every incident
halfedge. A single underlying Delaunay face triple can omit valid incident
generators and is not sufficient evidence.

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
| `segment_site_endpoint_binding.*` | live parabola limiter ownership |
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
| Point-graph deterministic records | implemented | publication blockers open |
| Cocircular point topology | guarded | named fail-loud boundary |
| Point/segment source parameterization | implemented helper | not fully wired |
| Point/segment endpoint binders | implemented and reviewed | not production graph |
| True-radius parabola clearance | pending | no production claim |
| Segment/segment feature domains | pending | no production claim |
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
12. Do the native gate, focused Python tests, legacy toolpath tests, and
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
