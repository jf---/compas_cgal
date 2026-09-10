from pathlib import Path

import numpy as np
from compas.colors import Color
from compas.colors import ColorMap
from compas.datastructures import Mesh
from compas.geometry import Box
from compas.geometry import Point
from compas.geometry import Polyline
from compas.geometry import Translation
from compas_viewer import Viewer
from compas_viewer.config import Config

from compas_cgal.geodesics import ExactGeodesicSolver
from compas_cgal.geodesics import exact_geodesic_distances
from compas_cgal.geodesics import heat_geodesic_distances
from compas_cgal.isolines import isolines


def make_mesh(V, F, offset):
    m = Mesh.from_vertices_and_faces(V, F)
    m.transform(Translation.from_vector([offset[0], offset[1], 0]))
    return m


def make_vertex_colors(values, colormap):
    return {i: colormap(v, minval=float(values.min()), maxval=float(values.max())) for i, v in enumerate(values)}


def make_partition_colors(ordinals):
    return {i: COLORS[int(s) % len(COLORS)] for i, s in enumerate(ordinals)}


def make_isolines(distances, isovalues, offset):
    m = mesh.copy()
    for key, d in zip(m.vertices(), distances):
        m.vertex_attribute(key, "distance", float(d))
    return [Polyline([[p[0] + offset[0], p[1] + offset[1], p[2]] for p in pts]) for pts in isolines(m, "distance", isovalues=isovalues, resample=False)]


def make_points(points, offset):
    return [Point(p[0] + offset[0], p[1] + offset[1], p[2]) for p in points]


def face_centroid_at(vertex):
    """Centroid of a face incident to a vertex: a source location that is not a mesh vertex."""
    return np.asarray(mesh.face_centroid(mesh.vertex_faces(vertex)[0]))


# =============================================================================
# Load mesh
# =============================================================================

FILE = Path(__file__).parent.parent.parent / "data" / "elephant.off"
mesh = Mesh.from_off(FILE)
mesh.quads_to_triangles()
V, F = mesh.to_vertices_and_faces()
V_np = np.array(V)

# =============================================================================
# Config
# =============================================================================

X_OFF, Y_OFF = 0.75, 1.0
N_ISOVALUES = 24
COLORS = [
    Color.red(),
    Color.orange(),
    Color.yellow(),
    Color.green(),
    Color.cyan(),
    Color.blue(),
    Color.purple(),
    Color.magenta(),
]
distance_cmap = ColorMap.from_two_colors(Color.blue(), Color.red())
error_cmap = ColorMap.from_two_colors(Color.white(), Color.red())

# Snapping a source to the nearest mesh vertex displaces it by a fraction of an
# edge, so the near-field isovalues of the last panel are expressed in multiples
# of the mean edge length rather than in absolute model units.
EDGE = float(np.mean([mesh.edge_length(edge) for edge in mesh.edges()]))
NEAR_ISOVALUES = [k * EDGE for k in (1, 2, 3, 4, 6, 8)]

# =============================================================================
# Vertex sources: the exact backend against the heat method
# =============================================================================

SOURCES = [0]

d_exact = exact_geodesic_distances((V, F), SOURCES)
d_heat = heat_geodesic_distances((V, F), SOURCES)
d_error = np.abs(d_heat - d_exact)

# Both fields are contoured at the same isovalues, which is what makes the two
# sets of curves comparable at all.
ISOVALUES = np.linspace(0.0, float(d_exact.max()), N_ISOVALUES + 2)[1:-1].tolist()

# =============================================================================
# Point sources: locations the heat method cannot express
# =============================================================================

# One source per bounding box corner, as in the heat method example, but taken
# at the centroid of an incident face instead of at the vertex itself.
bbox = Box.from_points(V_np)
corner_vertices = list(dict.fromkeys([int(np.argmin(np.linalg.norm(V_np - c, axis=1))) for c in bbox.vertices]))
POINTS = np.array([face_centroid_at(v) for v in corner_vertices])

# The nearest mesh vertex to each source point: where a vertex-only backend has
# to put the source instead.
SNAPPED = [int(np.argmin(np.linalg.norm(V_np - p, axis=1))) for p in POINTS]

solver = ExactGeodesicSolver((V, F))
d_points, nearest = solver.solve_from_points(POINTS, return_sources=True)  # first point query builds the AABB tree
d_snapped = solver.solve(SNAPPED)  # same solver, vertex sources, no tree needed

snapping = np.linalg.norm(POINTS - V_np[SNAPPED], axis=1)

print(f"mean edge length:            {EDGE:.6f}")
print(f"max |heat - exact|:          {d_error.max():.6f} ({d_error.max() / d_exact.max():.2%} of the field range)")
print(f"max snapping displacement:   {snapping.max():.6f} ({snapping.max() / EDGE:.2f} x mean edge)")

# =============================================================================
# Viz
# =============================================================================

config = Config()
config.camera.target = [X_OFF, -Y_OFF / 2, 0]
config.camera.position = [X_OFF, -2.0, 0.8]

viewer = Viewer(config=config)

# Row 1: Vertex Sources

g1 = viewer.scene.add_group("Vertex Sources")

g1.add(
    make_mesh(V, F, (0, 0)),
    use_vertexcolors=True,
    vertexcolor=make_vertex_colors(d_exact, distance_cmap),
    show_lines=False,
)

for pt in make_points(V_np[SOURCES], (0, 0)):
    g1.add(pt, pointcolor=Color.black(), pointsize=20)

for pl in make_isolines(d_exact, ISOVALUES, (X_OFF, 0)):
    g1.add(pl, linecolor=Color.red(), lineswidth=5)

for pl in make_isolines(d_heat, ISOVALUES, (X_OFF, 0)):
    g1.add(pl, linecolor=Color.blue(), lineswidth=5)

g1.add(
    make_mesh(V, F, (2 * X_OFF, 0)),
    use_vertexcolors=True,
    vertexcolor=make_vertex_colors(d_error, error_cmap),
    show_lines=False,
)

# Row 2: Point Sources

g2 = viewer.scene.add_group("Point Sources")

g2.add(
    make_mesh(V, F, (0, -Y_OFF)),
    use_vertexcolors=True,
    vertexcolor=make_vertex_colors(d_points, distance_cmap),
    show_lines=False,
)

for pt in make_points(POINTS, (0, -Y_OFF)):
    g2.add(pt, pointcolor=Color.black(), pointsize=20)

g2.add(
    make_mesh(V, F, (X_OFF, -Y_OFF)),
    use_vertexcolors=True,
    vertexcolor=make_partition_colors(nearest),
    show_lines=False,
)

for pl in make_isolines(d_points, NEAR_ISOVALUES, (2 * X_OFF, -Y_OFF)):
    g2.add(pl, linecolor=Color.red(), lineswidth=5)

for pl in make_isolines(d_snapped, NEAR_ISOVALUES, (2 * X_OFF, -Y_OFF)):
    g2.add(pl, linecolor=Color.blue(), lineswidth=5)

for pt in make_points(POINTS, (2 * X_OFF, -Y_OFF)):
    g2.add(pt, pointcolor=Color.black(), pointsize=20)

for pt in make_points(V_np[SNAPPED], (2 * X_OFF, -Y_OFF)):
    g2.add(pt, pointcolor=Color.blue(), pointsize=20)

viewer.show()
