from typing import Tuple
from typing import Union

import numpy as np
from compas.datastructures import Graph
from compas.geometry import Polygon
from compas.geometry import normal_polygon
from compas.tolerance import TOL

from compas_cgal._cgal import straight_skeleton_2

from .types import IntNx1
from .types import IntNx2
from .types import VerticesNumpy


def graph_from_skeleton_data(points: VerticesNumpy, indices: IntNx1, edges: IntNx2, edge_types: IntNx1) -> Graph:
    """Create a graph from the skeleton data.

    Parameters
    ----------
    points : :class:`numpy.ndarray`
        The vertices of the skeleton, each vertex defined by 3 spatial coordinates.
    indices : :class:`numpy.ndarray`
        The vertex indices of the skeleton, corresponding to the points.
    edges : :class:`numpy.ndarray`
        The edges of the skeleton, each edge defined by 2 vertex indices.
    edge_types : :class:`numpy.ndarray`
        The type per edge, `0` for inner bisector, `1` for bisector, and `2` for boundary.

    Returns
    -------
    :class:`compas.datastructures.Graph`
        The skeleton as a graph.
    """
    graph = Graph()
    for pt, i in zip(points, indices):
        graph.add_node(key=i, x=pt[0], y=pt[1], z=pt[2])

    for edge, etype in zip(edges, edge_types):
        edge = graph.add_edge(*edge)
        if etype == 0:
            graph.edge_attribute(edge, "inner_bisector", True)
        elif etype == 1:
            graph.edge_attribute(edge, "bisector", True)
        else:
            graph.edge_attribute(edge, "boundary", True)
    return graph


def create_interior_straight_skeleton(points, as_graph=True) -> Union[Graph, Tuple[VerticesNumpy, IntNx1, IntNx2, IntNx1]]:
    """Compute the skeleton of a polygon.

    Parameters
    ----------
    points : list of point coordinates or :class:`compas.geometry.Polygon`
        The points of the polygon.
    as_graph : bool, optional
        Whether the skeleton should be returned as a graph, defaults to `True`.

    Returns
    -------
    :attr:`compas.datastructures.Graph` or tuple of (vertices, indices, edges, edge_types)
        The skeleton of the polygon.

    Raises
    ------
    ValueError
        If the normal of the polygon is not [0, 0, 1].
    """
    points = list(points)
    normal = normal_polygon(points, True)
    if not TOL.is_allclose(normal, [0, 0, 1]):
        raise ValueError("The normal of the polygon should be [0, 0, 1]. The normal of the provided polygon is {}".format(normal))
    V = np.asarray(points, dtype=np.float64)
    points, indices, edges, edge_types = straight_skeleton_2.create_interior_straight_skeleton(V)
    if as_graph:
        return graph_from_skeleton_data(points, indices, edges, edge_types)
    return points, indices, edges, edge_types


def create_interior_straight_skeleton_with_holes(points, holes, as_graph=True) -> Union[Graph, Tuple[VerticesNumpy, IntNx1, IntNx2, IntNx1]]:
    """Compute the skeleton of a polygon with holes.

    Parameters
    ----------
    points : list of point coordinates or :class:`compas.geometry.Polygon`
        The points of the polygon.
    holes : list of list of point coordinates or list of :class:`compas.geometry.Polygon`
        The holes of the polygon.
    as_graph : bool, optional
        Whether the skeleton should be returned as a graph, defaults to `True`.

    Returns
    -------
    :attr:`compas.datastructures.Graph` or tuple of (vertices, indices, edges, edge_types)
        The skeleton of the polygon.

    Raises
    ------
    ValueError
        If the normal of the polygon is not [0, 0, 1].
        If the normal of a hole is not [0, 0, -1].
    """
    points = list(points)
    normal = normal_polygon(points, True)
    if not TOL.is_allclose(normal, [0, 0, 1]):
        raise ValueError("The normal of the polygon should be [0, 0, 1]. The normal of the provided polygon is {}".format(normal))
    V = np.asarray(points, dtype=np.float64)

    H = []
    for i, hole in enumerate(holes):
        points = list(hole)
        normal_hole = normal_polygon(points, True)
        if not TOL.is_allclose(normal_hole, [0, 0, -1]):
            raise ValueError("The normal of the hole should be [0, 0, -1]. The normal of the provided {}-th hole is {}".format(i, normal_hole))
        hole = np.asarray(points, dtype=np.float64)
        H.append(hole)
    points, indices, edges, edge_types = straight_skeleton_2.create_interior_straight_skeleton_with_holes(V, H)
    if as_graph:
        return graph_from_skeleton_data(points, indices, edges, edge_types)
    return points, indices, edges, edge_types


def create_offset_polygons_2(points, offset) -> list[Polygon]:
    """Compute the polygon offset.

    Parameters
    ----------
    points : list of point coordinates or :class:`compas.geometry.Polygon`
        The points of the polygon.
    offset : float
        The offset distance. If negative, the offset is outside the polygon, otherwise inside.

    Returns
    -------
    list[:class:`Polygon`]
        The offset polygon(s).

    Raises
    ------
    ValueError
        If the normal of the polygon is not [0, 0, 1].
    """
    points = list(points)
    normal = normal_polygon(points, True)
    if not TOL.is_allclose(normal, [0, 0, 1]):
        raise ValueError("The normal of the polygon should be [0, 0, 1]. The normal of the provided polygon is {}".format(normal))
    V = np.asarray(points, dtype=np.float64)
    offset = float(offset)
    if offset < 0:  # outside
        offset_polygons = straight_skeleton_2.create_offset_polygons_2_outer(V, abs(offset))[1:]  # first one is box
    else:  # inside
        offset_polygons = straight_skeleton_2.create_offset_polygons_2_inner(V, offset)
    return [Polygon(points.tolist()) for points in offset_polygons]


def create_weighted_offset_polygons_2(points, offset, weights) -> list[Polygon]:
    """Compute the polygon offset with weights.

    Parameters
    ----------
    points : list of point coordinates or :class:`compas.geometry.Polygon`
        The points of the polygon.
    offset : float
        The offset distance. If negative, the offset is outside the polygon, otherwise inside.
    weights : list of float
        The weights for each edge, starting with the edge between the last and the first point.

    Returns
    -------
    list[:class:`Polygon`]
        The offset polygon(s).

    Raises
    ------
    ValueError
        If the normal of the polygon is not [0, 0, 1].
    ValueError
        If the number of weights does not match the number of points.
    """
    points = list(points)
    normal = normal_polygon(points, True)
    if not TOL.is_allclose(normal, [0, 0, 1]):
        raise ValueError("The normal of the polygon should be [0, 0, 1]. The normal of the provided polygon is {}".format(normal))

    V = np.asarray(points, dtype=np.float64)
    offset = float(offset)
    W = np.asarray(weights, dtype=np.float64)
    if W.shape[0] != V.shape[0]:
        raise ValueError("The number of weights should be equal to the number of points %d != %d." % (W.shape[0], V.shape[0]))
    if offset < 0:
        offset_polygons = straight_skeleton_2.create_weighted_offset_polygons_2_outer(V, abs(offset), W)[1:]
    else:
        offset_polygons = straight_skeleton_2.create_weighted_offset_polygons_2_inner(V, offset, W)
    return [Polygon(points.tolist()) for points in offset_polygons]
