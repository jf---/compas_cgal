"""Shared Markdown table grammar for measurement-claim consumers."""

from __future__ import annotations

from typing import Sequence
from typing import Tuple
from typing import TypeVar

ErrorT = TypeVar("ErrorT", bound=Exception)
_UNSAFE_CELL_TOKENS = ("|", "\r", "\n", "\u2028", "\u2029")


def require_markdown_safe(value: object, field: str, error: type[ErrorT]) -> str:
    """Require one non-empty physical Markdown-table cell.

    Args:
        value: Candidate cell text.
        field: Human-readable owner used in diagnostics.
        error: Named exception type for the consumer boundary.

    Returns:
        The validated string.

    Raises:
        ErrorT: The value is empty, non-string, multiline, or contains a pipe.
    """
    if type(value) is not str or not value or any(token in value for token in _UNSAFE_CELL_TOKENS):
        raise error(f"{field}: must be Markdown-safe and one physical line")
    return value


def parse_pipe_table(
    lines: Sequence[str],
    *,
    header: str,
    rule: str,
    columns: int,
    field: str,
    error: type[ErrorT],
) -> Tuple[Tuple[str, ...], ...]:
    """Parse one contiguous pipe table and reject later pipe rows.

    Args:
        lines: Lines confined to the table's owning section.
        header: Exact table header.
        rule: Exact separator immediately following the header.
        columns: Required number of body cells.
        field: Human-readable owner used in diagnostics.
        error: Named exception type for the consumer boundary.

    Returns:
        Parsed, stripped body cells.

    Raises:
        ErrorT: The table grammar, width, or section remainder is invalid.
    """
    if lines.count(header) != 1 or lines.count(rule) != 1:
        raise error(f"{field}: requires its exact header and separator once")
    header_index = lines.index(header)
    if header_index + 1 >= len(lines) or lines[header_index + 1] != rule:
        raise error(f"{field}: separator must immediately follow its header")
    rows = []
    remainder_start = header_index + 2
    for index, line in enumerate(lines[header_index + 2 :], start=header_index + 2):
        if not line.startswith("|"):
            remainder_start = index
            break
        if not line.endswith("|"):
            raise error(f"{field}: row must start and end with an outer pipe")
        cells = tuple(cell.strip() for cell in line[1:-1].split("|"))
        if len(cells) != columns:
            raise error(f"{field}: row column count must be exactly {columns}, got {len(cells)}")
        rows.append(cells)
        remainder_start = index + 1
    if rows and any(line.startswith("|") for line in lines[remainder_start:]):
        raise error(f"{field}: section remainder contains another pipe row")
    return tuple(rows)
