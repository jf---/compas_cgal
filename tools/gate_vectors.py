"""The pinned gate vectors: one results file, validated on load.

`benchmarks/results/gate_vectors.json` is the single place the machining-quality
gate's measured numbers live. Every gate cell is a (case, criterion) pair whose
value is a decimal STRING at the precision the gate reports it, so a diff of the
file is a diff of the numbers a reader sees and two files compare without float
churn. This module owns the schema and nothing else: measuring is
`tools.gate_measure`, comparing is `tools.gate_ratchet`, rendering is
`tools.gate_table`; one admission entry's shape is `tools.gate_admissions`.

Schema ``gate-vectors/2`` adds a seventh key, ``admissions``: the append-only
list of ledger rulings that admitted new cases (Ruling 0.3-C2-amend). A
``gate-vectors/1`` document (exactly the six older keys) still loads, with no
admissions, because a parent commit is parsed with the current code; it is
always written back as ``/2``.

Criterion directions:

* ``lower`` -- a gated criterion; lower is better and the ratchet compares it.
* ``reported`` -- never ratcheted; checked for freshness only when it is an
  integer (no decimals), because a reported float's last digit is not
  guaranteed to agree across platforms (Ruling 0.3-CI-fresh).
* ``wallclock`` -- a time in seconds; never ratcheted and never checked for
  freshness, because wall clock is noise at the precision it is stored.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from decimal import Decimal
from decimal import InvalidOperation
from pathlib import Path
from typing import Literal
from typing import Mapping

from typing_extensions import Self

from tools.gate_admissions import Admission
from tools.gate_admissions import MalformedAdmissionError
from tools.gate_admissions import ReferenceReassignedError
from tools.gate_case_keys import is_case_key

SCHEMA = "gate-vectors/2"

# The schema before admissions; read (never written) so a parent commit's
# document still parses. Its documents carry exactly `LEGACY_KEYS`.
LEGACY_SCHEMA = "gate-vectors/1"
LEGACY_KEYS = frozenset({"schema", "rules_version", "instrument", "criteria", "cases", "rebaseline"})
DOCUMENT_KEYS = LEGACY_KEYS | {"admissions"}

# The gate rules the vectors were measured under. A bump re-baselines only the
# criteria its recorded approval lists as changed; every other gated criterion is
# still ratcheted across the bump (Rulings 0.3-C2, 0.3-bump-worsens): never
# silently.
RULES_VERSION = "1.2"

GATE_VECTORS_PATH = Path("benchmarks") / "results" / "gate_vectors.json"
GATE_TABLE_PATH = Path("benchmarks") / "results" / "gate_vectors.md"

Direction = Literal["lower", "reported", "wallclock"]
DIRECTIONS: frozenset[str] = frozenset({"lower", "reported", "wallclock"})
GROUPS: tuple[str, ...] = ("elementary", "cut", "speed", "longevity", "cost")

# A re-baseline approval must cite the ledger ruling that approved it, e.g.
# "Ruling 0.4-C1: ...". The ledger is the authority; this only refuses an
# approval that names none.
LEDGER_RULING_ID = re.compile(r"\bRuling \d+\.\d+-[A-Za-z0-9][A-Za-z0-9-]*\b")

# What a gated cell must satisfy, as the gate states it: a literal decimal, the
# case's cap in degrees, or the pocket's centre-domain component count.
REQUIRED_CAP = "cap"
REQUIRED_COMPONENTS = "centre_domain_components"


class MalformedGateVectorsError(Exception):
    """The gate-vectors document does not satisfy the gate-vectors schema."""


class MalformedRebaselineRecordError(MalformedGateVectorsError):
    """The ``rebaseline`` record of a gate-vectors document is malformed."""


REBASELINE_KEYS = frozenset({"from", "approval", "criteria", "retired"})


@dataclass(frozen=True)
class RebaselineRecord:
    """The recorded approval of a rules-version bump (Ruling 0.3-bump-worsens).

    Attributes:
        from_version: The parent document's rules version.
        approval: The user-approved ledger line; cites a ruling id.
        criteria: The gated criteria whose definition the bump changed, in
            record order. Only these are re-baselined; every other gated
            criterion is still checked for worsening on the bump commit.
        retired: The criteria the parent gated that this document no longer
            gates (removed, or demoted to ``reported``), in record order. Each
            must be acknowledged here; the ratchet refuses an unacknowledged
            drop and a retirement of a criterion the parent did not gate.
    """

    from_version: str
    approval: str
    criteria: tuple[str, ...]
    retired: tuple[str, ...]

    @classmethod
    def build(cls, raw: object, gated: tuple[str, ...]) -> Self | None:
        """Validate a decoded ``rebaseline`` value against the document's gated criteria.

        Args:
            raw: The decoded JSON value; an empty object means no record.
            gated: The document's gated criterion names.

        Returns:
            The validated record, or None for an empty object.

        Raises:
            MalformedRebaselineRecordError: The record violates the schema.
        """
        if not isinstance(raw, dict):
            raise MalformedRebaselineRecordError("rebaseline must be an object.")
        if not raw:
            return None
        if set(raw) != REBASELINE_KEYS:
            raise MalformedRebaselineRecordError(f"a rebaseline record carries exactly {sorted(REBASELINE_KEYS)} (the changed criteria), got {sorted(raw)}.")
        from_version, approval, criteria = raw["from"], raw["approval"], raw["criteria"]
        if not isinstance(from_version, str) or not from_version:
            raise MalformedRebaselineRecordError("rebaseline 'from' must be a non-empty string.")
        if not isinstance(approval, str) or not LEDGER_RULING_ID.search(approval):
            raise MalformedRebaselineRecordError(f"a rebaseline approval must cite a ledger ruling id ('Ruling <stage>.<task>-<id>'), got {approval!r}.")
        changed = _criterion_list(raw["criteria"], "criteria")
        retired = _criterion_list(raw["retired"], "retired")
        if not changed and not retired:
            raise MalformedRebaselineRecordError("a rebaseline record must be non-empty: a bump lists at least one changed ('criteria') or retired ('retired') criterion.")
        unknown = [name for name in changed if name not in gated]
        if unknown:
            raise MalformedRebaselineRecordError(f"rebaseline 'criteria' names {unknown}, which are not gated criteria of this document.")
        still_gated = [name for name in retired if name in gated]
        if still_gated:
            raise MalformedRebaselineRecordError(f"rebaseline 'retired' names {still_gated}, which this document still gates.")
        return cls(from_version=from_version, approval=approval, criteria=changed, retired=retired)

    def to_json(self) -> dict[str, object]:
        """The record as it is written."""
        return {"from": self.from_version, "approval": self.approval, "criteria": list(self.criteria), "retired": list(self.retired)}


def _criterion_list(value: object, key: str) -> tuple[str, ...]:
    """One list of criterion names of a rebaseline record: strings, no duplicates (possibly empty).

    Raises:
        MalformedRebaselineRecordError: *value* is not such a list.
    """
    if not isinstance(value, list) or not all(isinstance(name, str) for name in value):
        raise MalformedRebaselineRecordError(f"rebaseline {key!r} must be a list of criterion names, got {value!r}.")
    duplicates = sorted({name for name in value if value.count(name) > 1})
    if duplicates:
        raise MalformedRebaselineRecordError(f"rebaseline {key!r} lists {duplicates} more than once.")
    return tuple(value)


@dataclass(frozen=True)
class CriterionSpec:
    """How one criterion is stored, compared, and shown.

    Attributes:
        name: The JSON key, also the metric's field name in `benchmarks.quality`.
        group: The quality group the metric belongs to, for the table.
        direction: ``lower``, ``reported`` or ``wallclock``.
        decimals: Digits after the point at which the value is stored; 0 for a
            count.
        required: For a gated criterion, the gate's threshold: a decimal string,
            `REQUIRED_CAP`, or `REQUIRED_COMPONENTS`. None when not gated.
    """

    name: str
    group: str
    direction: Direction
    decimals: int
    required: str | None

    @classmethod
    def build(cls, name: str, raw: object) -> Self:
        """Validate one criterion entry.

        Args:
            name: The criterion's key.
            raw: The decoded JSON value.

        Returns:
            The validated specification.

        Raises:
            MalformedGateVectorsError: The entry violates the schema.
        """
        if not isinstance(raw, dict) or set(raw) != {"group", "direction", "decimals", "required"}:
            raise MalformedGateVectorsError(f"criterion {name!r} must be an object with exactly group, direction, decimals and required.")
        group, direction, decimals, required = raw["group"], raw["direction"], raw["decimals"], raw["required"]
        if group not in GROUPS:
            raise MalformedGateVectorsError(f"criterion {name!r} has unknown group {group!r}.")
        if direction not in DIRECTIONS:
            raise MalformedGateVectorsError(f"criterion {name!r} has unknown direction {direction!r}.")
        if type(decimals) is not int or decimals < 0:
            raise MalformedGateVectorsError(f"criterion {name!r} needs a non-negative integer number of decimals.")
        if (direction == "lower") != (required is not None):
            raise MalformedGateVectorsError(f"criterion {name!r}: exactly the gated criteria carry a required value.")
        if required is not None and required not in (REQUIRED_CAP, REQUIRED_COMPONENTS):
            _decimal(required, where=f"criterion {name!r} required")
        return cls(name=name, group=group, direction=direction, decimals=decimals, required=required)

    def format(self, value: float | int) -> str:
        """The stored string for *value* at this criterion's precision.

        Args:
            value: The measured number.

        Returns:
            The decimal string.
        """
        return f"{value:.{self.decimals}f}"

    def to_json(self) -> dict[str, object]:
        """The criterion entry as it is written."""
        return {"group": self.group, "direction": self.direction, "decimals": self.decimals, "required": self.required}


def _decimal(value: object, *, where: str) -> Decimal:
    if not isinstance(value, str):
        raise MalformedGateVectorsError(f"{where} must be a decimal string, got {type(value).__name__}.")
    try:
        number = Decimal(value)
    except InvalidOperation as error:
        raise MalformedGateVectorsError(f"{where} is not a decimal: {value!r}.") from error
    if not number.is_finite():
        raise MalformedGateVectorsError(f"{where} is not finite: {value!r}.")
    return number


@dataclass(frozen=True)
class GateVectors:
    """The whole validated document.

    Attributes:
        rules_version: The gate rules the cells were measured under.
        instrument: The instrument settings (tool, probes, grid).
        criteria: Every criterion, in table order.
        cases: Case key to criterion name to decimal string.
        rebaseline: For a document whose rules version differs from its
            parent's, the recorded approval of that re-baseline; None
            otherwise (written as ``{}``).
        admissions: The ledger rulings that admitted new cases, in append
            order; empty for a ``gate-vectors/1`` document.
        schema: The schema the document was read as; always written as
            `SCHEMA`. The ratchet accepts `LEGACY_SCHEMA` only for a parent.
    """

    rules_version: str
    instrument: Mapping[str, object]
    criteria: tuple[CriterionSpec, ...]
    cases: Mapping[str, Mapping[str, str]]
    rebaseline: RebaselineRecord | None
    admissions: tuple[Admission, ...] = ()
    schema: str = SCHEMA

    @classmethod
    def build(cls, raw: object) -> Self:
        """Validate a decoded gate-vectors document.

        Args:
            raw: The decoded JSON value.

        Returns:
            The validated vectors.

        Raises:
            MalformedGateVectorsError: The document violates the schema.
            MalformedAdmissionError: An admission entry violates the schema, or
                the admissions disagree with the document's cases.
        """
        if not isinstance(raw, dict):
            raise MalformedGateVectorsError("the gate-vectors document must be a JSON object.")
        schema = raw.get("schema")
        if schema not in (SCHEMA, LEGACY_SCHEMA):
            raise MalformedGateVectorsError(f"schema must be {SCHEMA!r} (or the readable {LEGACY_SCHEMA!r}), got {schema!r}.")
        expected_keys = DOCUMENT_KEYS if schema == SCHEMA else LEGACY_KEYS
        if set(raw) != expected_keys:
            raise MalformedGateVectorsError(f"a {schema} document must have exactly the keys {sorted(expected_keys)}, got {sorted(raw)}.")
        rules_version = raw["rules_version"]
        if not isinstance(rules_version, str) or not rules_version:
            raise MalformedGateVectorsError("rules_version must be a non-empty string.")
        instrument = raw["instrument"]
        if not isinstance(instrument, dict):
            raise MalformedGateVectorsError("instrument must be an object.")
        criteria_raw = raw["criteria"]
        if not isinstance(criteria_raw, dict) or not criteria_raw:
            raise MalformedGateVectorsError("criteria must be a non-empty object.")
        criteria = tuple(CriterionSpec.build(name, value) for name, value in criteria_raw.items())
        if not any(criterion.direction == "lower" for criterion in criteria):
            raise MalformedGateVectorsError("at least one criterion must be gated.")
        rebaseline = RebaselineRecord.build(raw["rebaseline"], tuple(criterion.name for criterion in criteria if criterion.direction == "lower"))
        cases_raw = raw["cases"]
        if not isinstance(cases_raw, dict) or not cases_raw:
            raise MalformedGateVectorsError("cases must be a non-empty object.")
        names = [criterion.name for criterion in criteria]
        cases: dict[str, dict[str, str]] = {}
        for key, cells in cases_raw.items():
            if not is_case_key(key):
                raise MalformedGateVectorsError(f"case key {key!r} must be 'cap-id/generator/pocket'.")
            if not isinstance(cells, dict) or list(cells) != names:
                raise MalformedGateVectorsError(f"case {key!r} must carry exactly the criteria, in criteria order.")
            for criterion in criteria:
                _decimal(cells[criterion.name], where=f"{key}/{criterion.name}")
            cases[key] = dict(cells)
        admissions = _admissions(raw.get("admissions", []), cases, names)
        return cls(rules_version=rules_version, instrument=dict(instrument), criteria=criteria, cases=cases, rebaseline=rebaseline, admissions=admissions, schema=schema)

    @property
    def gated(self) -> tuple[CriterionSpec, ...]:
        """The criteria the ratchet compares."""
        return tuple(criterion for criterion in self.criteria if criterion.direction == "lower")

    @property
    def freshness_checked(self) -> tuple[CriterionSpec, ...]:
        """The criteria a fresh measurement must reproduce exactly: gated cells and reported integers."""
        return tuple(criterion for criterion in self.criteria if criterion.direction == "lower" or (criterion.direction == "reported" and criterion.decimals == 0))

    @property
    def excluded(self) -> dict[str, str]:
        """The cases currently excluded, each with the exception class its generator raises.

        An exclusion holds until a later admission admits the case.
        """
        return replay_admissions(self.admissions)[1]

    @property
    def exempt_rows(self) -> dict[str, frozenset[str]]:
        """Per admitted case, the stage-0 all-case row criteria it is exempt from."""
        exempt: dict[str, frozenset[str]] = {}
        for admission in self.admissions:
            for case, names in admission.exempt_rows.items():
                exempt[case] = frozenset(names)
        return exempt

    def value(self, case: str, criterion: str) -> Decimal:
        """One cell as an exact decimal.

        Args:
            case: The case key.
            criterion: The criterion name.

        Returns:
            The stored value.
        """
        return Decimal(self.cases[case][criterion])

    def to_json(self) -> dict[str, object]:
        """The document as it is written."""
        return {
            "schema": SCHEMA,
            "rules_version": self.rules_version,
            "instrument": dict(self.instrument),
            "criteria": {criterion.name: criterion.to_json() for criterion in self.criteria},
            "rebaseline": {} if self.rebaseline is None else self.rebaseline.to_json(),
            "admissions": [admission.to_json() for admission in self.admissions],
            "cases": {key: dict(cells) for key, cells in self.cases.items()},
        }

    def dumps(self) -> str:
        """The canonical text of the document: two-space indent, trailing newline."""
        return json.dumps(self.to_json(), indent=2) + "\n"


def replay_admissions(admissions: tuple[Admission, ...]) -> tuple[set[str], dict[str, str]]:
    """Replay the admissions in order: the admitted cases and the cases still excluded.

    Args:
        admissions: The entries, in append order.

    Returns:
        The admitted case keys and the excluded case keys with their exception
        class names.

    Raises:
        MalformedAdmissionError: A case is admitted twice, or excluded while it
            is admitted or already excluded.
        ReferenceReassignedError: A later admission gives a generator another
            reference than its first admission did.
    """
    admitted: set[str] = set()
    excluded: dict[str, str] = {}
    references: dict[str, str | None] = {}
    for admission in admissions:
        for generator, reference in admission.reference.items():
            if generator in references and references[generator] != reference:
                raise ReferenceReassignedError(
                    f"admission {admission.ruling} references generator {generator!r} to {reference!r}; its first admission fixed {references[generator]!r}."
                )
            references.setdefault(generator, reference)
        for case in admission.cases:
            if case in admitted:
                raise MalformedAdmissionError(f"admission {admission.ruling} admits {case!r}, which is already admitted.")
            admitted.add(case)
            excluded.pop(case, None)
        for case, name in admission.excluded.items():
            if case in admitted:
                raise MalformedAdmissionError(f"admission {admission.ruling} excludes {case!r}, which is already admitted.")
            if case in excluded:
                raise MalformedAdmissionError(f"admission {admission.ruling} excludes {case!r}, which is already excluded.")
            excluded[case] = name
    return admitted, excluded


def _admissions(raw: object, cases: Mapping[str, object], criteria: list[str]) -> tuple[Admission, ...]:
    """Validate the ``admissions`` list against the document's cases and criteria.

    Args:
        raw: The decoded ``admissions`` value (``[]`` for a /1 document).
        cases: The document's case keys.
        criteria: The document's criterion names.

    Returns:
        The validated entries, in append order.

    Raises:
        MalformedAdmissionError: The list or an entry is malformed, an admitted
            case is not in the document, or an excluded case is.
        MalformedGateVectorsError: An exemption names a criterion the document
            does not carry.
    """
    if not isinstance(raw, list):
        raise MalformedAdmissionError(f"admissions must be a list, got {type(raw).__name__}.")
    admissions = tuple(Admission.build(entry) for entry in raw)
    admitted, excluded = replay_admissions(admissions)
    absent = sorted(case for case in admitted if case not in cases)
    if absent:
        raise MalformedAdmissionError(f"admitted cases {absent} are not in the document's cases.")
    present = sorted(case for case in excluded if case in cases)
    if present:
        raise MalformedAdmissionError(f"excluded cases {present} record cells; an excluded case is not in the document's cases.")
    for admission in admissions:
        for case, names in admission.exempt_rows.items():
            unknown = [name for name in names if name not in criteria]
            if unknown:
                raise MalformedGateVectorsError(f"admission {admission.ruling} exempts {case!r} from {unknown}, which the document does not carry.")
    return admissions


def load_gate_vectors(path: Path) -> GateVectors:
    """Read and validate a gate-vectors file.

    Args:
        path: The JSON file.

    Returns:
        The validated vectors.

    Raises:
        FileNotFoundError: The file does not exist.
        MalformedGateVectorsError: The file is not a valid gate-vectors document.
    """
    return parse_gate_vectors(path.read_text())


def parse_gate_vectors(text: str) -> GateVectors:
    """Validate the text of a gate-vectors document.

    Args:
        text: The JSON text.

    Returns:
        The validated vectors.

    Raises:
        MalformedGateVectorsError: The text is not valid JSON or not a valid
            gate-vectors document.
    """
    try:
        raw = json.loads(text)
    except json.JSONDecodeError as error:
        raise MalformedGateVectorsError(f"gate vectors are not valid JSON: {error}.") from error
    return GateVectors.build(raw)
