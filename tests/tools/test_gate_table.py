"""The docs' gate numbers are rendered from the pinned vectors, never typed."""

from __future__ import annotations

from pathlib import Path

from tools.gate_table import render_gate_table
from tools.gate_vectors import GATE_TABLE_PATH
from tools.gate_vectors import GATE_VECTORS_PATH
from tools.gate_vectors import load_gate_vectors

ROOT = Path(__file__).resolve().parents[2]
DOCS_PAGE = ROOT / "docs" / "machining_quality.md"
SNIPPET = f'--8<-- "{GATE_TABLE_PATH.as_posix()}"'


def test_the_committed_table_is_rendered_from_the_committed_vectors() -> None:
    assert (ROOT / GATE_TABLE_PATH).read_text() == render_gate_table(load_gate_vectors(ROOT / GATE_VECTORS_PATH))


def test_every_gate_cell_appears_in_the_table() -> None:
    vectors = load_gate_vectors(ROOT / GATE_VECTORS_PATH)
    table = render_gate_table(vectors)
    for criterion in vectors.criteria:
        assert f"`{criterion.name}`" in table


def test_the_machining_quality_page_includes_the_generated_table() -> None:
    assert SNIPPET in DOCS_PAGE.read_text()
