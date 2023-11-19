from __future__ import annotations
from .types import VerticesFaces
import numpy as np
from compas_cgal._cgal import measure
from compas.plugins import plugin


@plugin(category="trimesh", pluggable_name="trimesh_volume")
def volume(mesh: VerticesFaces) -> float:
    """Compute the volume of a closed triangle mesh.

    Parameters
    ----------
    mesh : :attr:`compas_cgal.types.VerticesFaces`
        The mesh.

    Returns
    -------
    float
        The volume of the mesh.

    Examples
    --------
    >>> from compas.geometry import Box
    >>> from compas_cgal.measure import volume

    >>> box = Box(1)
    >>> mesh = box.to_vertices_and_faces(triangulated=True)

    >>> volume(mesh)
    1.0

    """
    V, F = mesh
    V = np.asarray(V, dtype=np.float64)
    F = np.asarray(F, dtype=np.int32)
    return measure.volume(V, F)
