# Exact Geodesic Distances and Point Sources

![Exact Geodesic Distances and Point Sources](../assets/images/example_geodesics_exact.png)

The exact backend returns the polyhedral geodesic distance itself, and it will take a source
anywhere on the surface. The heat method returns a smoothed approximation of that distance, and
can only be seeded at mesh vertices. This example puts both differences on screen using the same
mesh, the same source sets and the same layout as
[Geodesic Distances and Isolines](example_geodesics.md), so that the backend is the only variable.

The visualization shows two rows:

* **Row 1 (Vertex Sources)**: one source vertex, solved by both backends
* **Row 2 (Point Sources)**: eight sources at face centroids, which only the exact backend accepts

Row 1 displays three columns:

* **Exact Field**: mesh colored by exact geodesic distance (blue → red gradient), source point in black
* **Exact vs Heat**: both fields contoured at the *same* isovalues — exact in red, heat in blue. Where
  the blue curve leaves the red one, the heat method's diffusion has smoothed the field.
* **Discrepancy**: mesh colored by `|heat - exact|` (white → red), which is where that smoothing lands

Row 2 displays three columns:

* **Point Field**: mesh colored by exact geodesic distance from eight face-centroid sources, source
  points in black
* **Nearest Source**: every vertex colored by which source is closest — the geodesic Voronoi partition,
  returned by the same call via `return_sources=True`
* **Snapping Cost**: near-field contours of the same eight sources, taken at their true locations (red,
  black points) and at the nearest mesh vertex (blue, blue points)

Key Features:

* Exact polyhedral geodesic distances from vertex sources via ``exact_geodesic_distances``
* Sources at arbitrary surface locations via ``exact_geodesic_distances_from_points``
* Nearest-source ordinals in the same call via ``return_sources=True``
* Reusable solver for either source kind via ``ExactGeodesicSolver``
* Contouring of any vertex scalar field via [``isolines``](example_isolines.md)

## What "exact" means here

Exact refers to the algorithm, not to the arithmetic. The shortest path across the polyhedron is
computed without algorithmic approximation — no time step, no diffusion, no dependence on triangle
quality — but it is evaluated in double precision under CGAL's
`Exact_predicates_inexact_constructions_kernel`. Geodesic distance is a construction, built from
unfoldings and square roots, so the returned numbers carry ordinary floating-point error. They are
exact in the sense that no algorithmic parameter moves them: what remains is the arithmetic, and the
polyhedron's own departure from whatever surface it stands for.

## Why the source location matters

A vertex-only backend has to move a source that is not a mesh vertex onto the nearest one. That
displacement is bounded by the local edge length, and it is an additive error on the entire field —
every distance is off by up to that much, everywhere, not just near the source. The last panel shows
it directly: the two families of contours are the same computation on the same mesh, differing only
in where the eight sources were allowed to sit. The example prints the measured displacement in
multiples of the mean edge length.

The displacement shrinks only with edge length, so a vertex-only backend buys accuracy by refining
the mesh — which changes the geometry being measured and multiplies the vertex count. A point source
removes the term instead of paying for it.

## Choosing a backend

The two entry points are one identifier apart, but they are not interchangeable:

* `heat_geodesic_distances` factorizes the mesh once and then answers each source set with two sparse
  solves. `HeatGeodesicSolver` amortizes that factorization over many queries, so many source sets on
  one mesh is the case it is built for.
* `exact_geodesic_distances` propagates a window sequence tree per source set, O(n^2 log n) in the worst
  case, and holds that tree in memory while it does. Its cost is dominated by the mesh, not by the
  number of sources.

Reach for the exact backend when the number itself has to be trusted — as a reference field, as the
input to a tolerance, or wherever the sources are not mesh vertices. Reach for the heat method when
a smooth field over many source sets is what is wanted.

```python
---8<--- "docs/examples/example_geodesics_exact.py"
```
