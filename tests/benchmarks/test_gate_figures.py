from __future__ import annotations

import os
import struct
import subprocess
import sys
from pathlib import Path

import matplotlib
import pytest
from matplotlib import rc_context
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from PIL import Image

# The module draws, so it must draw without a display here too.
matplotlib.use("Agg")

from benchmarks.gate import GATE_CAP_CASES  # noqa: E402
from benchmarks.gate import GATE_CAP_DEG  # noqa: E402
from benchmarks.gate import GATE_DEFAULT_CAP_ID  # noqa: E402
from benchmarks.gate import GATE_GENERATOR_NAMES  # noqa: E402
from benchmarks.errors import FigureDpiMismatchError  # noqa: E402

from benchmarks.gate_figures import HEAD_VECTORS  # noqa: E402
from benchmarks.gate_figures import REPO_ROOT  # noqa: E402
from benchmarks.gate_figures import StaleGateCellsError  # noqa: E402
from benchmarks.gate_figures import UnknownCapError  # noqa: E402
from benchmarks.gate_figures import UnregisteredGeneratorError  # noqa: E402
from benchmarks.gate_figures import VectorsSource  # noqa: E402
from benchmarks.gate_figures import check_caption_criteria  # noqa: E402
from benchmarks.gate_figures import compose_pocket  # noqa: E402
from benchmarks.gate_figures import resolve_cap  # noqa: E402
from benchmarks.gate_figures import resolve_generators  # noqa: E402
from benchmarks.plotting import FIGURE_DPI  # noqa: E402
from benchmarks.plotting import save_figure  # noqa: E402
from tools.gate_measure import case_key  # noqa: E402
from tools.gate_vectors import load_gate_vectors  # noqa: E402

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"
SMOKE_POCKET = "rect_12x8"

# PNG chunks that carry text (keyword NUL value) or a modification time; none
# may hold anything that is not a function of the drawing.
PNG_TEXT_CHUNKS = (b"tEXt", b"zTXt", b"iTXt")
PNG_TIME_CHUNK = b"tIME"
# Text keywords that carry the writer's version or a wall clock.
NON_DRAWING_PNG_KEYWORDS = frozenset({b"Software", b"Creation Time", b"Date"})
# Bytes per PNG chunk header (length + type) and trailer (CRC).
PNG_CHUNK_HEADER_BYTES = 8
PNG_CHUNK_CRC_BYTES = 4

# A savefig resolution no house figure uses, to show the rc setting has no say.
FOREIGN_SAVEFIG_DPI = 50
# A one-inch square figure: its pixel size is its save resolution.
ONE_INCH = 1.0

# The second compose of the byte-stability test runs in a fresh interpreter
# under a hash seed that differs from this process's, so set and dict-of-str
# ordering that leaks into the drawing shows up as a byte difference.
SUBPROCESS_COMPOSE = """
import sys
from pathlib import Path
from benchmarks.gate import GATE_CAP_DEG, GATE_DEFAULT_CAP_ID, GATE_GENERATOR_NAMES
from benchmarks.gate_figures import HEAD_VECTORS, VectorsSource, compose_pocket
from tools.gate_vectors import load_gate_vectors
head = VectorsSource(load_gate_vectors(HEAD_VECTORS), "pinned")
compose_pocket(sys.argv[1], GATE_DEFAULT_CAP_ID, GATE_CAP_DEG, GATE_GENERATOR_NAMES, head, None, False, Path(sys.argv[2]))
"""


def _png_chunks(data: bytes) -> list[tuple[bytes, bytes]]:
    """Split a PNG into (chunk type, chunk payload) pairs, in file order."""
    assert data.startswith(PNG_SIGNATURE)
    chunks = []
    at = len(PNG_SIGNATURE)
    while at < len(data):
        (length,) = struct.unpack(">I", data[at : at + 4])
        kind = data[at + 4 : at + PNG_CHUNK_HEADER_BYTES]
        payload = data[at + PNG_CHUNK_HEADER_BYTES : at + PNG_CHUNK_HEADER_BYTES + length]
        chunks.append((kind, payload))
        at += PNG_CHUNK_HEADER_BYTES + length + PNG_CHUNK_CRC_BYTES
    return chunks


def _other_hash_seed() -> str:
    """A PYTHONHASHSEED value different from this process's."""
    return "2" if os.environ.get("PYTHONHASHSEED") == "1" else "1"


def test_gate_generator_names_resolve_in_gate_order() -> None:
    assert resolve_generators(GATE_GENERATOR_NAMES) == tuple(GATE_GENERATOR_NAMES)


def test_unregistered_generator_raises() -> None:
    with pytest.raises(UnregisteredGeneratorError):
        resolve_generators(("no_such_generator",))


@pytest.mark.parametrize(("cap_id", "cap"), GATE_CAP_CASES)
def test_every_gate_cap_resolves(cap_id: str, cap: float) -> None:
    assert resolve_cap(f"{cap:g}") == (cap_id, cap)


def test_unknown_cap_raises() -> None:
    with pytest.raises(UnknownCapError):
        resolve_cap("7")


def test_caption_criteria_exist_in_the_pinned_vectors() -> None:
    check_caption_criteria(load_gate_vectors(HEAD_VECTORS), str(HEAD_VECTORS))


def test_composing_one_pocket_twice_writes_the_same_bytes(tmp_path: Path) -> None:
    head = VectorsSource(load_gate_vectors(HEAD_VECTORS), "pinned")
    in_process = tmp_path / "in_process"
    in_process.mkdir()
    first = compose_pocket(SMOKE_POCKET, GATE_DEFAULT_CAP_ID, GATE_CAP_DEG, GATE_GENERATOR_NAMES, head, None, False, in_process)
    other = tmp_path / "subprocess"
    other.mkdir()
    env = {**os.environ, "PYTHONHASHSEED": _other_hash_seed()}
    subprocess.run([sys.executable, "-c", SUBPROCESS_COMPOSE, SMOKE_POCKET, str(other)], cwd=REPO_ROOT, env=env, check=True)
    second = other / first.name
    data = first.read_bytes()
    chunks = _png_chunks(data)
    assert PNG_TIME_CHUNK not in {kind for kind, _ in chunks}
    keywords = {payload.split(b"\0", 1)[0] for kind, payload in chunks if kind in PNG_TEXT_CHUNKS}
    assert not keywords & NON_DRAWING_PNG_KEYWORDS
    assert data == second.read_bytes()


def test_check_fresh_refuses_a_stale_cell_before_drawing(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    vectors = load_gate_vectors(HEAD_VECTORS)
    head = VectorsSource(vectors, "pinned")
    stale_name = vectors.freshness_checked[0].name

    def measure_one_cell_changed(case: tuple[str, float, str, str]) -> tuple[str, dict[str, str]]:
        cap_id, _, generator, pocket = case
        key = case_key(cap_id, generator, pocket)
        cells = dict(vectors.cases[key])
        cells[stale_name] = f"{cells[stale_name]}1"
        return key, cells

    monkeypatch.setattr("benchmarks.gate_figures.measure_case", measure_one_cell_changed)
    with pytest.raises(StaleGateCellsError, match=stale_name):
        compose_pocket(SMOKE_POCKET, GATE_DEFAULT_CAP_ID, GATE_CAP_DEG, GATE_GENERATOR_NAMES, head, None, True, tmp_path)
    assert list(tmp_path.iterdir()) == []


def _one_inch_figure(dpi: float) -> Figure:
    figure = Figure(figsize=(ONE_INCH, ONE_INCH), dpi=dpi)
    FigureCanvasAgg(figure)
    return figure


def test_save_figure_resolution_ignores_the_savefig_rc(tmp_path: Path) -> None:
    with rc_context({"savefig.dpi": FOREIGN_SAVEFIG_DPI}):
        written = save_figure(_one_inch_figure(FIGURE_DPI), tmp_path / "square.png")
    with Image.open(written) as image:
        assert image.size == (round(FIGURE_DPI), round(FIGURE_DPI))


def test_save_figure_refuses_a_figure_at_another_resolution(tmp_path: Path) -> None:
    with pytest.raises(FigureDpiMismatchError):
        save_figure(_one_inch_figure(FOREIGN_SAVEFIG_DPI), tmp_path / "square.png")
    assert list(tmp_path.iterdir()) == []
