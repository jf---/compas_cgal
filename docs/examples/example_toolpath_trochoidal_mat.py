from collections import defaultdict

from compas.colors import Color
from compas.geometry import Arc
from compas.geometry import Circle
from compas.geometry import Line
from compas.geometry import Polygon
from compas.geometry import Polyline
from compas.scene import get_sceneobject_cls
from compas.scene import register
from compas.scene.exceptions import SceneObjectNotRegisteredError
from compas_viewer import Viewer
from compas_viewer.scene.geometryobject import GeometryObject

from compas_cgal.straight_skeleton_2 import interior_straight_skeleton
from compas_cgal.toolpath import trochoidal_mat_toolpath_circular

polygon = Polygon(
    [
        (-1.91, 3.59, 0.0),
        (-5.53, -5.22, 0.0),
        (-0.39, -1.98, 0.0),
        (2.98, -5.51, 0.0),
        (4.83, -2.02, 0.0),
        (9.70, -3.63, 0.0),
        (12.23, 1.25, 0.0),
        (3.42, 0.66, 0.0),
        (2.92, 4.03, 0.0),
    ]
)

operations = trochoidal_mat_toolpath_circular(
    polygon,
    tool_diameter=1.0,
    stepover=0.5,
    pitch=0.75,
    mat_scale=1.0,
    lead_in=0.15,
    lead_out=0.15,
    link_paths=True,
    optimize_order=True,
    cut_z=-0.2,
    clearance_z=2.0,
)


class ArcObject(GeometryObject):
    geometry: Arc

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.show_lines = True

    @property
    def lines(self):
        return self.geometry.to_polyline(n=self.u).lines


def ensure_viewer_arc_registered():
    try:
        get_sceneobject_cls(Arc, context="Viewer")
    except SceneObjectNotRegisteredError:
        register(Arc, ArcObject, context="Viewer")


def add_curve(scene, curve, **kwargs):
    """Add a curve to the scene, falling back to polyline tessellation."""
    if isinstance(curve, Line):
        scene.add(Polyline([curve.start, curve.end]), **kwargs)
    elif isinstance(curve, Circle):
        scene.add(curve.to_polyline(n=64), **kwargs)
    elif isinstance(curve, Arc):
        try:
            scene.add(curve, **kwargs)
        except SceneObjectNotRegisteredError:
            scene.add(curve.to_polyline(n=64), **kwargs)


# ==============================================================================
# Color palette: 10 distinct hues for MAT edge groups
# ==============================================================================

MAT_PALETTE = [
    Color.from_rgb255(230, 60, 50),
    Color.from_rgb255(50, 140, 230),
    Color.from_rgb255(40, 180, 70),
    Color.from_rgb255(220, 160, 30),
    Color.from_rgb255(160, 50, 210),
    Color.from_rgb255(230, 120, 50),
    Color.from_rgb255(50, 190, 190),
    Color.from_rgb255(200, 60, 140),
    Color.from_rgb255(100, 180, 50),
    Color.from_rgb255(100, 80, 200),
]

# Linewidth per operation type
OPERATION_WIDTH = {
    "cut": 2.5,
    "lead_in": 1.8,
    "lead_out": 1.8,
    "link": 1.2,
    "retract": 1.2,
    "plunge": 1.2,
}

# Lighten percentage per operation type (0 = full color, higher = more faded)
OPERATION_LIGHTEN = {
    "cut": 0,
    "lead_in": 30,
    "lead_out": 30,
    "link": 60,
    "retract": 60,
    "plunge": 60,
}


def skeleton_bisector_edges(polygon):
    """Return inner bisector + bisector edges of the straight skeleton as Lines."""
    graph = interior_straight_skeleton(polygon.points)
    edges = []
    for u, v in graph.edges():
        if graph.edge_attribute((u, v), "inner_bisector") or graph.edge_attribute((u, v), "bisector"):
            p0 = graph.node_coordinates(u)
            p1 = graph.node_coordinates(v)
            edges.append(Line(p0, p1))
    return edges


def _point_to_segment_sq(px, py, ax, ay, bx, by):
    """Squared distance from point (px,py) to segment (a,b)."""
    dx, dy = bx - ax, by - ay
    len_sq = dx * dx + dy * dy
    if len_sq < 1e-24:
        return (px - ax) ** 2 + (py - ay) ** 2
    t = max(0.0, min(1.0, ((px - ax) * dx + (py - ay) * dy) / len_sq))
    return (px - (ax + t * dx)) ** 2 + (py - (ay + t * dy)) ** 2


def group_cut_centers(ops):
    """Return list of (x, y) centers of all circle/arc cut operations in a group."""
    centers = []
    for op in ops:
        if op.operation != "cut":
            continue
        g = op.geometry
        if isinstance(g, (Circle, Arc)):
            c = g.frame.point
            centers.append((float(c[0]), float(c[1])))
    return centers


def match_skeleton_edge(edge, groups_centers):
    """Match a skeleton edge to the toolpath group whose cut arc centers lie closest.

    For each group, compute the mean squared distance from its arc centers to the
    skeleton segment. The group with the smallest mean distance wins.
    """
    ax, ay = float(edge.start[0]), float(edge.start[1])
    bx, by = float(edge.end[0]), float(edge.end[1])
    best_idx, best_cost = None, float("inf")
    for pidx, centers in groups_centers.items():
        if not centers:
            continue
        cost = sum(_point_to_segment_sq(cx, cy, ax, ay, bx, by) for cx, cy in centers) / len(centers)
        if cost < best_cost:
            best_cost = cost
            best_idx = pidx
    return best_idx


def fade_color(color, pct):
    """Lighten color by percentage (0 = unchanged, 100 = white)."""
    return color.lightened(int(round(pct)))


# ==============================================================================
# Visualize
# ==============================================================================

ensure_viewer_arc_registered()

viewer = Viewer()
viewer.config.renderer.show_grid = False

viewer.scene.add(polygon, show_faces=False, linecolor=Color.black())

# Group operations by path_index
groups = defaultdict(list)
for op in operations:
    groups[op.path_index].append(op)

# Arc centers per group for skeleton-to-group matching via point-to-segment distance
groups_centers = {pidx: group_cut_centers(ops) for pidx, ops in groups.items()}

skeleton_group = viewer.scene.add_group(name="Skeleton")
unmatched_color = Color.from_rgb255(120, 120, 120)
for edge in skeleton_bisector_edges(polygon):
    pidx = match_skeleton_edge(edge, groups_centers)
    color = fade_color(MAT_PALETTE[pidx % len(MAT_PALETTE)], 50) if pidx is not None else unmatched_color
    add_curve(viewer.scene, edge, parent=skeleton_group, linecolor=color, linewidth=6.0)

for path_index in sorted(groups):
    ops = groups[path_index]
    base_color = MAT_PALETTE[path_index % len(MAT_PALETTE)]

    # Toolpath operations → per-edge group
    toolpath_group = viewer.scene.add_group(name=f"MAT {path_index}")
    for op in ops:
        width = OPERATION_WIDTH.get(op.operation, 1.5)
        lighten = OPERATION_LIGHTEN.get(op.operation, 40)
        color = fade_color(base_color, lighten)
        add_curve(viewer.scene, op.geometry, parent=toolpath_group, linecolor=color, linewidth=width)

viewer.show()
