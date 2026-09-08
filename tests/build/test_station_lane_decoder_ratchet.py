"""Mechanically prevent a text crossing from returning to the station lane.

Stage 2 inverted the station source: `StationEventSource2` carries
`compas_cgal::exact::Rational` and *projects* the `exact-rational-v1` text, so
nothing in the lane decodes a number out of a string any more. These checks are
the design's decoder ratchet, scoped to the files stage 2 owns.

**The counts move in neither direction.** A rise means a consumer started
parsing again — the regression the inversion removed. A fall means a decoder
*definition* was deleted, which is stage 6 and needs explicit user permission.
A one-directional ratchet would silently permit the second.

The remaining occurrences are deliberate and named:

- `station_source.cpp` keeps the definition plus the one call inside
  `ExactRational2::build`, the legacy decoder retained so
  `exact_station_attestation_gate` can keep proving the projection agrees with
  it byte for byte.
- `station_classifier.cpp` keeps its definition only; it is dormant and marked
  `[[maybe_unused]]`.
- `segment_strata.cpp`, `circle_strata.cpp` and `segment_oracle.cpp` keep the
  callers that belong to the segment and full-circle lanes, which stages 3 and 4
  own.

## Why the text-read check is scoped by receiver type, not by file

`tool_radius()` and `cap_chord_ratio()` are accessors of **both**
`StationEventSource2` and `SegmentEventSource2`, and `segment_strata.cpp` binds
the identifier `source` to each of them in different functions. A file-wide ban
on `source.tool_radius().text()` therefore fires on
`segment_branches_at` and `segment_branch_pair_dispositions`, which legitimately
read a `SegmentEventSource2` — stage 3's carrier, not stage 2's. A ratchet that
fires on correct code gets edited away by the next person, which is worse than
no ratchet.

So the reads are scoped: `station_bound_regions` splits a translation unit at
its column-0 closing braces and keeps the definitions that bind a
`StationEventSource2` identifier, and the ban applies to that identifier inside
that definition only. `test_the_scope_finder_separates_the_two_source_bindings`
pins that this separation is real and not an accident of the current text.

For the same reason the scoped rule bans *reading a station value as text*, not
decoding in general: `circle_strata.cpp:546` round-trips the **full-circle**
lane's own cap ratio through `parse_rational(exact_rational_text(...))` inside a
function that also binds a station `source`. That crossing is real and stage 4
owns it; the occurrence count above is what holds it still.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict
from typing import List
from typing import Tuple

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
LANE = REPOSITORY_ROOT / "src" / "continuous_tea_2"

# Measured 2026-09-08 at 003e144c. Before stage 2 these were, in the same order:
# 4, 4, 12, 6, 11 -- a total of 37. Stage 2 removed 13 call sites and no
# definitions.
EXPECTED_DECODER_OCCURRENCES = {
    "src/continuous_tea_2/station_source.cpp": 2,
    "src/continuous_tea_2/station_classifier.cpp": 1,
    "src/continuous_tea_2/segment_strata.cpp": 8,
    "src/continuous_tea_2/circle_strata.cpp": 2,
    "src/continuous_tea_2/segment_oracle.cpp": 11,
}

# Repository-wide `parse_rational` definitions. 15 before stage 2 and 15 after:
# the inversion removed call sites only. Deleting a definition is stage 6.
EXPECTED_DECODER_DEFINITIONS = 15

# `[[maybe_unused]] Rational parse_rational(` is a definition too -- the dormant
# one in station_classifier.cpp carries the attribute.
DECODER_DEFINITION = re.compile(
    r"^[ \t]*(?:\[\[[^\]]*\]\][ \t]*)?Rational[ \t]+parse_rational[ \t]*\(",
    re.MULTILINE,
)

# These files are clang-formatted, so every top-level definition ends at a
# closing brace in column 0.
TOP_LEVEL_CLOSE = re.compile(r"^\}", re.MULTILINE)

# A variable or parameter bound to a station source. The trailing character
# separates a binding from a function *declaration* whose return type is
# `StationEventSource2` -- a declaration is followed by `(`, a binding by one of
# these. `StationEventSource2::build` is excluded by `::` matching nothing here.
STATION_BINDING = re.compile(r"\bStationEventSource2\b\s*&?\s*(?P<name>[A-Za-z_]\w*)\s*[=;,)\{]")

STATION_ACCESSORS = ("center_x", "center_y", "tool_radius", "cap_chord_ratio")

# Every top-level definition in these files that binds a station source.
# Asserted, so the scoped read check below cannot silently cover nothing:
# station_source.cpp is the carrier itself and reaches its members through
# `this`; segment_oracle.cpp builds a station and hands it straight to
# `classify_station_cell` without ever binding it.
EXPECTED_STATION_REGIONS = {
    "src/continuous_tea_2/station_source.cpp": 0,
    "src/continuous_tea_2/station_classifier.cpp": 1,
    "src/continuous_tea_2/segment_strata.cpp": 1,
    "src/continuous_tea_2/circle_strata.cpp": 2,
    "src/continuous_tea_2/segment_oracle.cpp": 0,
}

STATION_CONSUMERS = (
    "src/continuous_tea_2/station_classifier.cpp",
    "src/continuous_tea_2/segment_strata.cpp",
    "src/continuous_tea_2/circle_strata.cpp",
    "src/continuous_tea_2/segment_oracle.cpp",
)


def _read(relative_path: str) -> str:
    return (REPOSITORY_ROOT / relative_path).read_text(encoding="utf-8")


def top_level_definitions(source: str) -> List[str]:
    """Split a translation unit at its column-0 closing braces."""
    definitions: List[str] = []
    start = 0
    for match in TOP_LEVEL_CLOSE.finditer(source):
        definitions.append(source[start : match.end()])
        start = match.end()
    tail = source[start:]
    if tail.strip():
        definitions.append(tail)
    return definitions


def station_bound_regions(source: str) -> List[Tuple[str, str]]:
    """Return `(identifier, definition text)` for each station-bound receiver.

    Args:
        source: One C++ translation unit.

    Returns:
        One entry per (definition, identifier) pair where the identifier is
        bound to a `StationEventSource2`.
    """
    regions: List[Tuple[str, str]] = []
    for definition in top_level_definitions(source):
        for name in sorted({match.group("name") for match in STATION_BINDING.finditer(definition)}):
            regions.append((name, definition))
    return regions


def test_station_lane_decoder_count_moves_in_neither_direction() -> None:
    observed = {path: _read(path).count("parse_rational") for path in EXPECTED_DECODER_OCCURRENCES}
    drift = {path: (EXPECTED_DECODER_OCCURRENCES[path], count) for path, count in observed.items() if count != EXPECTED_DECODER_OCCURRENCES[path]}
    assert not drift, "station-lane decoder count moved (path: expected -> observed): " + ", ".join(
        f"{path}: {expected} -> {actual}" for path, (expected, actual) in sorted(drift.items())
    )
    assert sum(observed.values()) == 24


def test_no_decoder_definition_was_removed() -> None:
    """A fall in the counts above must not be a deleted definition.

    Removing a decoder is stage 6 and requires explicit user permission, so this
    asserts equality rather than a bound.
    """
    definitions: Dict[str, int] = {}
    for path in (
        sorted(REPOSITORY_ROOT.glob("src/**/*.cpp"))
        + sorted(REPOSITORY_ROOT.glob("src/**/*.h"))
        + sorted(REPOSITORY_ROOT.glob("tests/**/*.cpp"))
        + sorted(REPOSITORY_ROOT.glob("tests/**/*.h"))
    ):
        count = len(DECODER_DEFINITION.findall(path.read_text(encoding="utf-8")))
        if count:
            definitions[str(path.relative_to(REPOSITORY_ROOT))] = count
    total = sum(definitions.values())
    assert total == EXPECTED_DECODER_DEFINITIONS, (
        f"parse_rational definitions moved {EXPECTED_DECODER_DEFINITIONS} -> {total}; removal is stage 6 and needs explicit user permission. Observed: {definitions}"
    )


def test_the_scope_finder_still_finds_the_station_code() -> None:
    """The scoped read check must not degrade into covering nothing."""
    observed = {path: len(station_bound_regions(_read(path))) for path in EXPECTED_STATION_REGIONS}
    assert observed == EXPECTED_STATION_REGIONS


def test_the_scope_finder_separates_the_two_source_bindings() -> None:
    """`segment_strata.cpp` binds `source` to both carriers; only one is ours.

    `construct_station_cell_stratum` takes a `StationEventSource2`;
    `segment_branches_at` and `segment_branch_pair_dispositions` take a
    `SegmentEventSource2`, stage 3's carrier. If the scope ever stopped
    separating them, the read check would either miss the station or fire on the
    segment.

    The assertion is about the *binding*, not about how the segment lane reads
    its source: stage 3 removes those text reads, and a ratchet that fails when
    the next stage succeeds is a ratchet someone deletes.
    """
    source = _read("src/continuous_tea_2/segment_strata.cpp")
    regions = station_bound_regions(source)
    assert [name for name, _ in regions] == ["source"]
    station_definition = regions[0][1]
    assert "construct_station_cell_stratum" in station_definition
    assert "const StationEventSource2& source" in station_definition
    # The other carrier binds the same identifier, outside every station region.
    assert "const SegmentEventSource2& source" in source
    assert "const SegmentEventSource2& source" not in station_definition


def test_no_station_receiver_is_read_as_text() -> None:
    """No station value is decoded from, or rendered to, text by a consumer."""
    violations: List[str] = []
    for path in EXPECTED_STATION_REGIONS:
        for name, definition in station_bound_regions(_read(path)):
            for accessor in STATION_ACCESSORS:
                for pattern in (
                    f"{name}.{accessor}().text()",
                    f"{name}.attestation().{accessor}.text()",
                    f"{name}.attestation().{accessor}.numerator()",
                    f"{name}.attestation().{accessor}.denominator()",
                ):
                    if pattern in definition:
                        violations.append(f"{path}: {pattern}")
    assert not violations, "a station value is being read as text again: " + ", ".join(violations)


def test_no_consumer_reaches_through_the_attestation_view() -> None:
    """`StationAttestation2` is a derived view, never a carrier.

    Its four `ExactRational2` fields are the digest's, and reading a number back
    out of them would reinstate the crossing stage 2 removed by another route.
    `attestation()` is therefore called nowhere outside its own definition.
    """
    for path in STATION_CONSUMERS:
        assert ".attestation()" not in _read(path), f"{path} reads the station attestation view"


def test_the_station_source_projects_through_the_single_door() -> None:
    """`to_canonical` is called from exactly one place in the station lane."""
    source = _read("src/continuous_tea_2/station_source.cpp")
    assert source.count("to_canonical(") == 1
    assert "ExactRational2 ExactRational2::project(" in source
