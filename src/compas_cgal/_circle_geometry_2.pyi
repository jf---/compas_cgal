"""CGAL decisions on exactly injected binary64 geometry; coordinates in mm."""

class InvalidCircleGeometryError(RuntimeError): ...
class NoCircleIntersectionError(RuntimeError): ...
class InvalidBoundaryProjectionError(RuntimeError): ...
class AmbiguousBoundaryProjectionError(RuntimeError): ...
class BoundaryProjectionDistanceError(RuntimeError): ...

def project_boundary_contact(boundary: list[tuple[float, float]], contact: tuple[float, float], maximum_distance: float) -> tuple[int, float]:
    """Choose nearest offset side exactly; report its approximate parameter."""

def disk_contains_disk(first: tuple[float, float], first_radius: float, second: tuple[float, float], second_radius: float) -> bool:
    """Test closed-disk containment, including exact tangency."""

def swept_disk_intersection(first: tuple[float, float], first_radius: float, second: tuple[float, float], second_radius: float, tool_radius: float) -> tuple[float, float]:
    """Report chord x in mm and squared height in mm² after exact existence check."""

def orientation(a: tuple[float, float], b: tuple[float, float], c: tuple[float, float]) -> int:
    """Return -1 clockwise, 0 collinear, +1 counterclockwise."""

def corrected_engagement_cosine_squared(first: tuple[float, float], first_radius: float, second: tuple[float, float], second_radius: float, tool_radius: float) -> float:
    """Report squared chord cosine after exact corrected-intersection checks."""
