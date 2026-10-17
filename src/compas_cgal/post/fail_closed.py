"""The fail-closed post: map, write, parse the text back, check it, and only then hand it out."""

from collections.abc import Sequence

from compas_cgal.post.checks import check_program
from compas_cgal.post.errors import RoundTripMismatchError
from compas_cgal.post.machine import PostSetup
from compas_cgal.post.mapping import SourcePath
from compas_cgal.post.mapping import map_operations
from compas_cgal.post.parser import parse_program
from compas_cgal.post.writer import write_program
from compas_cgal.toolpath import ToolpathOperation


def post_operations(operations: Sequence[ToolpathOperation], setup: PostSetup, source: SourcePath) -> str:
    """G-code text for *operations*, or a named error; never unchecked text.

    The checks run on the text as written (parsed back), and the parse must
    reproduce the mapped program exactly.

    Args:
        operations: The generated path.
        setup: The setup.
        source: The path's planes and tool.

    Returns:
        The checked program text.

    Raises:
        PostError: Any mapping, parsing or check failure (see each stage).
    """
    posted = map_operations(operations, setup, source)
    text = write_program(posted.program, setup)
    parsed = parse_program(text)
    if parsed != posted.program:
        raise RoundTripMismatchError("the written text does not parse back to the mapped program.")
    check_program(parsed, setup)
    return text
