from __future__ import annotations

import json

import pytest

from benchmarks.external import EXTERNAL_FAMILY
from benchmarks.external import EmptyExternalCorpusError
from benchmarks.external import ExternalCorpusError
from benchmarks.external import load_profiles
from benchmarks.external import profile_from_points

# A tool small enough that every ring below admits it, so a rejection in these
# tests is always the one the test is about.
TEST_TOOL = 1.0
TEST_CAP_DEG = 120.0

UNIT_SQUARE_10 = [[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0]]


def _write(tmp_path, name: str, payload: object) -> None:
    """Write one profile file, payload verbatim so malformed cases are expressible."""
    (tmp_path / f"{name}.json").write_text(json.dumps(payload))


def _profile(tmp_path, name: str, points: list) -> None:
    (tmp_path / f"{name}.json").write_text(json.dumps({"name": name, "points": points}))


def test_loads_every_usable_profile_sorted_by_name(tmp_path) -> None:
    _profile(tmp_path, "b", [[0.0, 0.0], [20.0, 0.0], [20.0, 8.0], [0.0, 8.0]])
    _profile(tmp_path, "a", UNIT_SQUARE_10)
    corpus = load_profiles(tmp_path, tool_diameter=TEST_TOOL, tea_cap_deg=TEST_CAP_DEG)
    assert [s.name for s in corpus.specs] == ["a", "b"]
    assert corpus.rejected == ()


def test_every_spec_is_filed_under_the_external_family(tmp_path) -> None:
    """Report rows have to be attributable to geometry this project did not author."""
    _profile(tmp_path, "a", UNIT_SQUARE_10)
    corpus = load_profiles(tmp_path, tool_diameter=TEST_TOOL, tea_cap_deg=TEST_CAP_DEG)
    assert all(s.family == EXTERNAL_FAMILY for s in corpus.specs)


def test_a_profile_below_the_vertex_floor_is_rejected_with_its_reason(tmp_path) -> None:
    _profile(tmp_path, "tri", [[0.0, 0.0], [10.0, 0.0], [5.0, 9.0]])
    _profile(tmp_path, "quad", UNIT_SQUARE_10)
    corpus = load_profiles(tmp_path, tool_diameter=TEST_TOOL, tea_cap_deg=TEST_CAP_DEG, min_vertices=4)
    assert [s.name for s in corpus.specs] == ["quad"]
    assert [r.name for r in corpus.rejected] == ["tri"]
    assert "3" in corpus.rejected[0].reason


def test_rejections_are_reported_rather_than_silently_dropped(tmp_path) -> None:
    """Dropping the hard profiles is how a resolution rate flatters itself.

    The corpus exists to measure how often the certifier gives up on geometry it
    did not author. A loader that swallowed the awkward rings would bias exactly
    that number, and the bias would be invisible in the result.
    """
    _profile(tmp_path, "ok", UNIT_SQUARE_10)
    _profile(tmp_path, "bowtie", [[0.0, 0.0], [10.0, 10.0], [10.0, 0.0], [0.0, 10.0]])
    corpus = load_profiles(tmp_path, tool_diameter=TEST_TOOL, tea_cap_deg=TEST_CAP_DEG)
    assert [s.name for s in corpus.specs] == ["ok"]
    assert [r.name for r in corpus.rejected] == ["bowtie"]
    assert "PocketNotSimpleError" in corpus.rejected[0].reason


def test_a_profile_too_small_for_the_tool_is_rejected_with_its_reason(tmp_path) -> None:
    _profile(tmp_path, "speck", [[0.0, 0.0], [0.2, 0.0], [0.2, 0.2], [0.0, 0.2]])
    corpus = load_profiles(tmp_path, tool_diameter=TEST_TOOL, tea_cap_deg=TEST_CAP_DEG)
    assert corpus.specs == ()
    assert "DegeneratePocketError" in corpus.rejected[0].reason


def test_a_repeated_closing_vertex_is_the_same_ring(tmp_path) -> None:
    """Both encodings of a closed ring are common in exported profile data.

    Left in place the duplicate makes a zero-length edge, which the simplicity
    test reads as a crossing and the ring is thrown out as dataset noise.
    """
    _profile(tmp_path, "closed", UNIT_SQUARE_10 + [[0.0, 0.0]])
    corpus = load_profiles(tmp_path, tool_diameter=TEST_TOOL, tea_cap_deg=TEST_CAP_DEG)
    assert [s.name for s in corpus.specs] == ["closed"]
    assert len(corpus.specs[0].polygon.points) == 4


def test_missing_directory_fails_loudly(tmp_path) -> None:
    with pytest.raises(ExternalCorpusError):
        load_profiles(tmp_path / "nope", tool_diameter=TEST_TOOL, tea_cap_deg=TEST_CAP_DEG)


def test_a_directory_holding_no_profiles_fails_loudly(tmp_path) -> None:
    """An empty corpus measuring nothing must not read as a corpus that passed."""
    with pytest.raises(EmptyExternalCorpusError):
        load_profiles(tmp_path, tool_diameter=TEST_TOOL, tea_cap_deg=TEST_CAP_DEG)


def test_unparseable_json_fails_loudly(tmp_path) -> None:
    (tmp_path / "bad.json").write_text("{not json")
    with pytest.raises(ExternalCorpusError):
        load_profiles(tmp_path, tool_diameter=TEST_TOOL, tea_cap_deg=TEST_CAP_DEG)


def test_a_profile_missing_its_points_fails_loudly(tmp_path) -> None:
    _write(tmp_path, "bad", {"name": "bad"})
    with pytest.raises(ExternalCorpusError):
        load_profiles(tmp_path, tool_diameter=TEST_TOOL, tea_cap_deg=TEST_CAP_DEG)


def test_a_profile_whose_vertices_are_not_pairs_fails_loudly(tmp_path) -> None:
    """An operator error in the conversion step, not a property of the geometry."""
    _write(tmp_path, "bad", {"name": "bad", "points": [[0.0, 0.0], [1.0], [2.0, 2.0], [0.0, 3.0]]})
    with pytest.raises(ExternalCorpusError):
        load_profiles(tmp_path, tool_diameter=TEST_TOOL, tea_cap_deg=TEST_CAP_DEG)


def test_a_profile_whose_points_are_not_a_list_fails_loudly(tmp_path) -> None:
    _write(tmp_path, "bad", {"name": "bad", "points": "square"})
    with pytest.raises(ExternalCorpusError):
        load_profiles(tmp_path, tool_diameter=TEST_TOOL, tea_cap_deg=TEST_CAP_DEG)


def test_profile_from_points_builds_a_spec() -> None:
    spec = profile_from_points("x", UNIT_SQUARE_10, tool_diameter=TEST_TOOL, tea_cap_deg=TEST_CAP_DEG)
    assert spec.name == "x"
    assert spec.family == EXTERNAL_FAMILY
    assert len(spec.polygon.points) == 4
    assert spec.params["vertices"] == pytest.approx(4.0)


def test_profile_from_points_rejects_a_ring_the_corpus_cannot_machine() -> None:
    """The factory owns the error model; the loader turns it into a rejection."""
    from benchmarks.errors import DegeneratePocketError

    with pytest.raises(DegeneratePocketError):
        profile_from_points("speck", [[0.0, 0.0], [0.2, 0.0], [0.2, 0.2], [0.0, 0.2]], tool_diameter=TEST_TOOL, tea_cap_deg=TEST_CAP_DEG)
