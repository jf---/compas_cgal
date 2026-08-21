from __future__ import annotations

import json

from benchmarks.measurement import MeasurementRecord
from benchmarks.report import render_markdown
from benchmarks.report import write_report

# Banned as literal text: the practice is to ORDER the document conclusion-first,
# not to announce that you have (CLAUDE.md, documentation format).
BANNED_LABELS = ("BLUF", "TL;DR", "Bottom Line", "## Summary")


def _record(name: str, certify: float = 1.0, truly_exceeding: int = 0, uncertified: int = 0, unresolved: int = 0, error: str | None = None) -> MeasurementRecord:
    """A measured record, or a failed one when *error* is given."""
    base = MeasurementRecord.failed(name, "necks", {"pinch": 2.0}, 1.0, 120.0, error or "placeholder")
    if error is not None:
        return base
    return MeasurementRecord.from_dict(
        {
            **base.to_dict(),
            "error": None,
            "generate_seconds": 0.01,
            "certify_seconds": certify,
            "operations": 10,
            "cut_operations": 8,
            "stations": 100,
            "max_tea_deg": 130.0,
            "truly_exceeding": truly_exceeding,
            "uncertified": uncertified,
            "unresolved": unresolved,
        }
    )


def test_conclusion_precedes_every_table() -> None:
    text = render_markdown([_record("a", certify=1.0), _record("b", certify=9.0, truly_exceeding=3)])
    assert text.index("over its engagement cap") < text.index("| instance |")


def test_no_banned_summary_labels_appear() -> None:
    text = render_markdown([_record("a")])
    for banned in BANNED_LABELS:
        assert banned not in text


def test_both_cap_columns_are_present_and_distinguished() -> None:
    """A reader meeting these numbers cold must not be able to confuse them.

    Both columns appear in the table, and the legend that separates them is
    emitted before the table rather than left to a footnote nobody reaches.
    """
    text = render_markdown([_record("a", truly_exceeding=2, uncertified=7)])
    header = next(line for line in text.splitlines() if line.startswith("| instance |"))
    assert "truly exceeding" in header
    assert "uncertified" in header
    assert text.index("LOWER BOUND") < text.index("| instance |")
    assert "could not be *proved*" in text


def test_the_conclusion_counts_exceedance_not_uncertifiability() -> None:
    """The headline sentence must state what exceeds, not what could not be proved.

    The regression this guards: an instance with SEVEN unprovable operations and
    ZERO measured exceedances used to be summarised as seven violations, which
    made a generator that had improved read as a regression.
    """
    text = render_markdown([_record("a", truly_exceeding=0, uncertified=7, unresolved=7)])
    conclusion = text.split("## Per-instance measurements")[0]
    assert "**0 of 1** measured instances put the tool over its engagement cap" in conclusion
    assert "at **0** cut motions in total" in conclusion
    assert "A further **1** could not be proved under the cap" in conclusion


def test_an_exceeding_instance_is_not_also_counted_as_merely_unproved() -> None:
    """The two headline sets are disjoint, so no instance is reported twice.

    An instance with a demonstrated exceedance is necessarily uncertified as well;
    counting both sets whole would report it under both headings and re-create the
    ambiguity the split exists to remove.
    """
    text = render_markdown([_record("a", truly_exceeding=2, uncertified=5)])
    conclusion = text.split("## Per-instance measurements")[0]
    assert "**1 of 1** measured instances put the tool over its engagement cap" in conclusion
    assert "A further **0** could not be proved under the cap" in conclusion


def test_exceeding_motions_are_totalled_across_instances() -> None:
    text = render_markdown([_record("a", truly_exceeding=2), _record("b", truly_exceeding=3)])
    conclusion = text.split("## Per-instance measurements")[0]
    assert "**2 of 2** measured instances put the tool over its engagement cap" in conclusion
    assert "at **5** cut motions in total" in conclusion


def test_failures_are_reported_not_silently_dropped() -> None:
    text = render_markdown([_record("ok"), _record("bad", error="DegeneratePocketError: narrow")])
    assert "bad" in text
    assert "DegeneratePocketError" in text


def test_an_all_failed_corpus_claims_no_conclusion() -> None:
    text = render_markdown([_record("bad", error="DegeneratePocketError: narrow")])
    assert "no conclusion is available" in text


def test_write_report_emits_both_artifacts(tmp_path) -> None:
    md_path, json_path = write_report([_record("a", truly_exceeding=1, uncertified=4)], tmp_path)
    assert md_path.exists() and json_path.exists()
    payload = json.loads(json_path.read_text())
    assert payload[0]["name"] == "a"
    # The JSON sidecar carries both columns under their own names, so a consumer
    # cannot silently re-conflate them downstream.
    assert payload[0]["truly_exceeding"] == 1
    assert payload[0]["uncertified"] == 4
