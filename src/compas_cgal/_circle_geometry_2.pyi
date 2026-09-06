"""CGAL decisions on exactly injected binary64 geometry; coordinates in mm."""

class InvalidCircleGeometryError(RuntimeError): ...
class NoCircleIntersectionError(RuntimeError): ...

def disk_contains_disk(first: tuple[float, float], first_radius: float, second: tuple[float, float], second_radius: float) -> bool:
    """Test closed-disk containment, including exact tangency."""

def swept_disk_intersection(first: tuple[float, float], first_radius: float, second: tuple[float, float], second_radius: float, tool_radius: float) -> tuple[float, float]:
    """Report chord x in mm and squared height in mm² after exact existence check."""

def orientation(a: tuple[float, float], b: tuple[float, float], c: tuple[float, float]) -> int:
    """Return -1 clockwise, 0 collinear, +1 counterclockwise."""
