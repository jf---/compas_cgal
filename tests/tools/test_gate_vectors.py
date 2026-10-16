"""The gate-vectors schema: version 2 carries an append-only admissions list (Ruling 0.3-C2-amend)."""

from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from benchmarks.gate import FROZEN_PAIR
from tools.gate_admissions import Admission
from tools.gate_admissions import MalformedAdmissionError
from tools.gate_admissions import ReferenceReassignedError
from tools.gate_vectors import GATE_VECTORS_PATH
from tools.gate_vectors import SCHEMA
from tools.gate_vectors import GateVectors
from tools.gate_vectors import MalformedGateVectorsError
from tools.gate_vectors import parse_gate_vectors

ROOT = Path(__file__).resolve().parents[2]
COMMITTED = ROOT / GATE_VECTORS_PATH
SIBLING = "cap-120-default/engagement_controlled/rect_12x8"
NEW_CASE = "cap-120-default/pipeline/rect_12x8"
RULING = "0.3-C2-amend"


def _committed() -> dict:
    """The committed vectors as the frozen pair pinned them: its cases only, no admission (a synthetic-column base)."""
    raw = json.loads(COMMITTED.read_text())
    raw["cases"] = {key: cells for key, cells in raw["cases"].items() if key.split("/")[1] in FROZEN_PAIR}
    raw["admissions"] = []
    return raw


def _with_case(raw: dict, case: str) -> dict:
    out = copy.deepcopy(raw)
    out["cases"][case] = copy.deepcopy(raw["cases"][SIBLING])
    return out


def _entry(**overrides: object) -> dict:
    entry: dict[str, object] = {"ruling": RULING, "cases": [NEW_CASE], "reference": {"pipeline": None}, "exempt_rows": {}, "excluded": {}}
    entry.update(overrides)
    return entry


def _admitted(*entries: dict, extra_cases: tuple[str, ...] = (NEW_CASE,)) -> dict:
    out = _committed()
    for case in extra_cases:
        out = _with_case(out, case)
    out["admissions"] = list(entries)
    return out


def test_the_schema_is_version_2_and_the_committed_document_admits_only_task_1_1s_columns() -> None:
    raw = json.loads(COMMITTED.read_text())
    assert SCHEMA == "gate-vectors/2"
    assert raw["schema"] == SCHEMA
    (admission,) = GateVectors.build(raw).admissions
    assert admission.ruling == RULING
    assert {case.split("/")[1] for case in admission.cases} == {"rho_regulated", "pipeline"}
    assert dict(admission.reference) == {"rho_regulated": None, "pipeline": "rho_regulated"}


def test_a_schema_1_six_key_document_reads_as_no_admissions() -> None:
    """The parent commit is parsed with the new code: a /1 document is still valid."""
    raw = _committed()
    del raw["admissions"]
    raw["schema"] = "gate-vectors/1"
    vectors = GateVectors.build(raw)
    assert vectors.admissions == ()
    assert vectors.to_json()["schema"] == SCHEMA


def test_a_schema_1_document_with_a_seventh_key_is_malformed() -> None:
    raw = _committed()
    raw["schema"] = "gate-vectors/1"
    with pytest.raises(MalformedGateVectorsError, match="admissions"):
        GateVectors.build(raw)


def test_a_schema_2_document_without_admissions_is_malformed() -> None:
    raw = _committed()
    del raw["admissions"]
    with pytest.raises(MalformedGateVectorsError, match="admissions"):
        GateVectors.build(raw)


def test_an_admission_round_trips_through_the_written_json() -> None:
    entry = _entry(
        cases=[NEW_CASE, "cap-120-default/rho_regulated/rect_12x8"],
        reference={"pipeline": "rho_regulated", "rho_regulated": None},
        exempt_rows={"cap-120-default/rho_regulated/rect_12x8": ["misplaced_plunges"]},
        excluded={"cap-120-default/rho_regulated/L_shape": "compas_cgal.errors.ConstructionPreconditionError"},
    )
    raw = _admitted(entry, extra_cases=(NEW_CASE, "cap-120-default/rho_regulated/rect_12x8"))
    vectors = GateVectors.build(raw)
    assert vectors.to_json()["admissions"] == [entry]
    assert parse_gate_vectors(vectors.dumps()) == vectors
    (admission,) = vectors.admissions
    assert admission == Admission(
        ruling=RULING,
        cases=(NEW_CASE, "cap-120-default/rho_regulated/rect_12x8"),
        reference={"pipeline": "rho_regulated", "rho_regulated": None},
        exempt_rows={"cap-120-default/rho_regulated/rect_12x8": ("misplaced_plunges",)},
        excluded={"cap-120-default/rho_regulated/L_shape": "compas_cgal.errors.ConstructionPreconditionError"},
    )
    assert vectors.excluded == {"cap-120-default/rho_regulated/L_shape": "compas_cgal.errors.ConstructionPreconditionError"}
    assert vectors.exempt_rows == {"cap-120-default/rho_regulated/rect_12x8": frozenset({"misplaced_plunges"})}


@pytest.mark.parametrize(
    ("overrides", "reason"),
    [
        ({"ruling": ""}, "ruling"),
        ({"ruling": "the user said yes"}, "ruling"),
        ({"cases": []}, "non-empty"),
        ({"cases": [NEW_CASE, NEW_CASE]}, "more than once"),
        ({"cases": ["not-a-case-key"]}, "not-a-case-key"),
        ({"reference": {}}, "pipeline"),
        ({"reference": {"pipeline": None, "other": None}}, "other"),
        ({"reference": {"pipeline": "pipeline"}}, "itself"),
        ({"exempt_rows": {"cap-120-default/other/rect_12x8": ["misplaced_plunges"]}}, "other"),
        ({"exempt_rows": {NEW_CASE: ["tangent_breaks"]}}, "tangent_breaks"),
        ({"exempt_rows": {NEW_CASE: []}}, "non-empty"),
        ({"exempt_rows": {NEW_CASE: ["misplaced_plunges", "misplaced_plunges"]}}, "more than once"),
        ({"excluded": {NEW_CASE: "pkg.errors.SomeError"}}, NEW_CASE),
        ({"excluded": {"cap-120-default/pipeline/L_shape": "not a class"}}, "not a class"),
    ],
    ids=[
        "empty-ruling",
        "ruling-not-an-id",
        "no-cases",
        "duplicate-case",
        "bad-case-key",
        "reference-missing-generator",
        "reference-extra-generator",
        "reference-to-itself",
        "exempt-case-not-admitted",
        "exempt-not-a-stage-0-row",
        "exempt-empty",
        "exempt-duplicate",
        "excluded-and-admitted",
        "excluded-bad-class-name",
    ],
)
def test_a_malformed_admission_is_refused_by_name(overrides: dict, reason: str) -> None:
    with pytest.raises(MalformedAdmissionError, match=reason):
        Admission.build(_entry(**overrides))


def test_an_admission_with_unknown_or_missing_keys_is_refused() -> None:
    with pytest.raises(MalformedAdmissionError, match="exactly"):
        Admission.build({**_entry(), "note": "x"})
    entry = _entry()
    del entry["excluded"]
    with pytest.raises(MalformedAdmissionError, match="exactly"):
        Admission.build(entry)
    with pytest.raises(MalformedAdmissionError, match="object"):
        Admission.build([])


def test_admissions_must_be_a_list() -> None:
    raw = _admitted()
    raw["admissions"] = {}
    with pytest.raises(MalformedAdmissionError, match="list"):
        GateVectors.build(raw)


def test_an_admitted_case_must_be_in_the_document() -> None:
    with pytest.raises(MalformedAdmissionError, match=NEW_CASE):
        GateVectors.build(_admitted(_entry(), extra_cases=()))


def test_an_excluded_case_must_not_be_in_the_document() -> None:
    excluded = "cap-120-default/pipeline/L_shape"
    with pytest.raises(MalformedAdmissionError, match=excluded):
        GateVectors.build(_admitted(_entry(excluded={excluded: "pkg.errors.SomeError"}), extra_cases=(NEW_CASE, excluded)))


def test_a_case_admitted_twice_is_refused() -> None:
    with pytest.raises(MalformedAdmissionError, match="already admitted"):
        GateVectors.build(_admitted(_entry(), _entry(ruling="1.1-closure")))


def test_a_case_excluded_twice_is_refused() -> None:
    excluded = {"cap-120-default/pipeline/L_shape": "pkg.errors.SomeError"}
    other = "cap-120-default/pipeline/rect_20x12"
    raw = _admitted(_entry(excluded=excluded), _entry(ruling="1.1-closure", cases=[other], excluded=excluded), extra_cases=(NEW_CASE, other))
    with pytest.raises(MalformedAdmissionError, match="already excluded"):
        GateVectors.build(raw)


def test_an_exempt_criterion_must_be_a_criterion_of_the_document() -> None:
    raw = _admitted(_entry(exempt_rows={NEW_CASE: ["misplaced_plunges"]}))
    del raw["criteria"]["misplaced_plunges"]
    for cells in raw["cases"].values():
        del cells["misplaced_plunges"]
    with pytest.raises(MalformedGateVectorsError, match="misplaced_plunges"):
        GateVectors.build(raw)


def test_a_later_admission_lifts_an_earlier_exclusion() -> None:
    """Append-only, so the way out of an exclusion is a later entry admitting the case."""
    lifted = "cap-120-default/pipeline/L_shape"
    first = _entry(excluded={lifted: "pkg.errors.SomeError"})
    second = _entry(ruling="1.1-closure", cases=[lifted], reference={"pipeline": None})
    vectors = GateVectors.build(_admitted(first, second, extra_cases=(NEW_CASE, lifted)))
    assert vectors.excluded == {}
    with pytest.raises(MalformedAdmissionError, match="already admitted"):
        GateVectors.build(_admitted(second, _entry(ruling="1.1-x", cases=[NEW_CASE], excluded={lifted: "pkg.errors.SomeError"}), extra_cases=(NEW_CASE, lifted)))


def test_an_excluded_exception_is_named_by_its_qualified_class_name() -> None:
    """A bare class name could be matched by a same-named class of another module."""
    with pytest.raises(MalformedAdmissionError, match="ConstructionPreconditionError"):
        Admission.build(_entry(excluded={"cap-120-default/pipeline/L_shape": "ConstructionPreconditionError"}))


def test_a_generator_reference_is_fixed_at_its_first_admission() -> None:
    """Per-column reference: a later entry cannot re-point a generator, to null included."""
    later = "cap-120-default/pipeline/rect_20x12"
    rho = "cap-120-default/rho_regulated/rect_12x8"
    first = _entry(cases=[NEW_CASE, rho], reference={"pipeline": "rho_regulated", "rho_regulated": None})
    for changed in (None, "engagement_controlled"):
        raw = _admitted(first, _entry(ruling="1.1-closure", cases=[later], reference={"pipeline": changed}), extra_cases=(NEW_CASE, rho, later))
        with pytest.raises(ReferenceReassignedError, match="pipeline"):
            GateVectors.build(raw)
    kept = _admitted(first, _entry(ruling="1.1-closure", cases=[later], reference={"pipeline": "rho_regulated"}), extra_cases=(NEW_CASE, rho, later))
    assert len(GateVectors.build(kept).admissions) == 2


def test_a_document_remembers_the_schema_it_was_read_as() -> None:
    raw = _committed()
    assert GateVectors.build(raw).schema == SCHEMA
    del raw["admissions"]
    raw["schema"] = "gate-vectors/1"
    assert GateVectors.build(raw).schema == "gate-vectors/1"
