from __future__ import annotations

import json
from pathlib import Path

import pytest

from benchmarks.errors import MalformedHeldReferenceCaseError
from benchmarks.held_reference_geometry import PdfPoint2
from benchmarks.held_reference_geometry import PdfPointUnit
from benchmarks.held_reference_geometry import SourceLine
from tools.held_reference_case_generator import CANONICAL_FILENAMES
from tools.held_reference_case_generator import _normalization_transform
from tools.held_reference_case_generator import _source_bounds
from tools.held_reference_case_generator import _source_signed_area_twice
from tools.held_reference_case_generator import _write_or_check_payloads


def _clockwise_square(offset_x: float = 0.0, offset_y: float = 0.0) -> tuple[SourceLine, ...]:
    points = tuple(PdfPoint2.build(x + offset_x, y + offset_y) for x, y in ((0.0, 0.0), (0.0, 4.0), (6.0, 4.0), (6.0, 0.0)))
    return tuple(SourceLine.build(points[index], points[(index + 1) % 4]) for index in range(4))


def test_source_bounds_use_true_dimensions_and_are_translation_invariant() -> None:
    original = _clockwise_square()
    translated = _clockwise_square(120.0, -83.0)

    original_minimum, original_maximum = _source_bounds(original)
    translated_minimum, translated_maximum = _source_bounds(translated)
    assert float(original_maximum.x) - float(original_minimum.x) == 6.0
    assert float(original_maximum.y) - float(original_minimum.y) == 4.0
    assert float(translated_maximum.x) - float(translated_minimum.x) == 6.0
    assert float(translated_maximum.y) - float(translated_minimum.y) == 4.0


def test_source_cycle_requires_clockwise_area_before_reflection() -> None:
    assert _source_signed_area_twice(_clockwise_square()) < 0


def test_normalization_is_translation_invariant_and_idempotent_at_unit_scale() -> None:
    original = _clockwise_square()
    translated = _clockwise_square(120.0, -83.0)
    original_transform = _normalization_transform(original, PdfPointUnit(1.0))
    translated_transform = _normalization_transform(translated, PdfPointUnit(1.0))

    original_points = tuple(original_transform.point(line.start) for line in original)
    translated_points = tuple(translated_transform.point(line.start) for line in translated)
    assert original_points == translated_points
    assert float(original_transform.scale) == 1.0


def test_writer_uses_only_fixed_names_and_semantic_check(tmp_path: Path) -> None:
    payloads = {name: {"name": name, "values": [1, 2, 3]} for name in CANONICAL_FILENAMES}
    _write_or_check_payloads(payloads, directory=tmp_path, check=False)

    assert tuple(sorted(path.name for path in tmp_path.iterdir())) == tuple(sorted(CANONICAL_FILENAMES.values()))
    for name, filename in CANONICAL_FILENAMES.items():
        payload = payloads[name]
        (tmp_path / filename).write_text(json.dumps(payload, separators=(",", ":")))
    _write_or_check_payloads(payloads, directory=tmp_path, check=True)


def test_check_fails_on_semantic_difference(tmp_path: Path) -> None:
    payloads = {name: {"name": name} for name in CANONICAL_FILENAMES}
    _write_or_check_payloads(payloads, directory=tmp_path, check=False)
    first = next(iter(CANONICAL_FILENAMES.values()))
    (tmp_path / first).write_text('{"name":"wrong"}')

    with pytest.raises(MalformedHeldReferenceCaseError):
        _write_or_check_payloads(payloads, directory=tmp_path, check=True)
