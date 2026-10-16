"""The admission record of the gate vectors (Ruling 0.3-C2-amend, route 2).

A new case enters the pinned vectors in one of two ways: it beats the existing
generators' minimum for its (cap, pocket) somewhere and loses nowhere (Ruling
0.3-C2 option (i)), or a ledger ruling admits it by name. An admission is an
entry of the document's append-only ``admissions`` list:

* ``ruling`` -- the ledger ruling id that admits the cases, e.g. ``0.3-C2-amend``;
* ``cases`` -- the admitted case keys, each measured and present in the document;
* ``reference`` -- each admitted generator's reference generator, or null. An
  admitted case whose generator has a reference must be no worse than the
  reference's case at the same (cap, pocket) on every gated cell; with null,
  the admitted cells become the baseline on the ruling's authority;
* ``exempt_rows`` -- per admitted case, the stage-0 all-case acceptance rows
  (`benchmarks.stage_acceptance`) the case is exempt from (G4);
* ``excluded`` -- per case whose generator raises, the exception class name it
  raises; such a case records no cells (G3).

This module owns one entry's shape. The document-level invariants (admitted
cases present, excluded cases absent, nothing admitted twice) are
`tools.gate_vectors`; the commit-level ones (append-only, reference rule, every
new case covered) are `tools.gate_ratchet`.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Mapping

from typing_extensions import Self

from benchmarks.stage_acceptance import ALL_CASES
from benchmarks.stage_acceptance import STAGE_ACCEPTANCE
from benchmarks.stage_acceptance import VectorRow
from tools.gate_case_keys import CASE_KEY_SEPARATOR
from tools.gate_case_keys import case_generator
from tools.gate_case_keys import is_case_key

# A ledger ruling id as the ledger writes it: "<stage>.<task>-<id>", e.g.
# "0.3-C2-amend" or "1.1-closure".
RULING_ID = re.compile(r"\d+\.\d+-[A-Za-z0-9][A-Za-z0-9-]*")

# An exception class by its qualified name, ``module.Class`` as
# `qualified_class_name` prints it; a bare name could match a same-named class
# of another module.
EXCEPTION_CLASS_NAME = re.compile(r"[A-Za-z_][A-Za-z0-9_]*(\.[A-Za-z_][A-Za-z0-9_]*)+")

ADMISSION_KEYS = frozenset({"ruling", "cases", "reference", "exempt_rows", "excluded"})

# The stage whose all-case vector rows an admission may exempt a case from.
EXEMPTABLE_STAGE = "stage-0"


class MalformedAdmissionError(Exception):
    """An entry of the gate-vectors ``admissions`` list violates the schema."""


class ReferenceReassignedError(MalformedAdmissionError):
    """A later admission gives a generator another reference than its first admission did."""


def qualified_class_name(error: BaseException) -> str:
    """The ``module.Class`` name of *error*'s class, as ``excluded`` records it."""
    return f"{type(error).__module__}.{type(error).__qualname__}"


def exemptable_criteria() -> frozenset[str]:
    """The criteria of the stage-0 all-case vector rows: the only rows ``exempt_rows`` may name."""
    return frozenset(row.criterion for row in STAGE_ACCEPTANCE[EXEMPTABLE_STAGE] if isinstance(row, VectorRow) and row.case == ALL_CASES)


@dataclass(frozen=True)
class Admission:
    """One ruling's admission of new cases into the gate vectors.

    Attributes:
        ruling: The ledger ruling id that admits the cases.
        cases: The admitted case keys, in record order; unique, non-empty.
        reference: For each generator of `cases`, its reference generator, or
            None for an unconditional baseline. Fixed at the generator's first
            admission; a later entry must repeat it.
        exempt_rows: Per admitted case, the stage-0 all-case row criteria it is
            exempt from.
        excluded: Per case that is not measured, the qualified class name
            (``module.Class``) of the exception its generator raises; disjoint
            from `cases`.
    """

    ruling: str
    cases: tuple[str, ...]
    reference: Mapping[str, str | None]
    exempt_rows: Mapping[str, tuple[str, ...]]
    excluded: Mapping[str, str]

    @classmethod
    def build(cls, raw: object) -> Self:
        """Validate one decoded ``admissions`` entry.

        Args:
            raw: The decoded JSON value.

        Returns:
            The validated admission.

        Raises:
            MalformedAdmissionError: The entry violates the schema.
        """
        if not isinstance(raw, dict):
            raise MalformedAdmissionError(f"an admission must be an object, got {type(raw).__name__}.")
        if set(raw) != ADMISSION_KEYS:
            raise MalformedAdmissionError(f"an admission carries exactly {sorted(ADMISSION_KEYS)}, got {sorted(raw)}.")
        ruling = raw["ruling"]
        if not isinstance(ruling, str) or not RULING_ID.fullmatch(ruling):
            raise MalformedAdmissionError(f"an admission's ruling must be a ledger ruling id ('<stage>.<task>-<id>'), got {ruling!r}.")
        cases = _unique_strings(raw["cases"], f"admission {ruling} 'cases'")
        if not cases:
            raise MalformedAdmissionError(f"admission {ruling} must admit a non-empty list of cases.")
        for case in cases:
            _require_case_key(case, ruling)
        reference = _reference(raw["reference"], cases, ruling)
        exempt_rows = _exempt_rows(raw["exempt_rows"], cases, ruling)
        excluded = _excluded(raw["excluded"], cases, ruling)
        return cls(ruling=ruling, cases=cases, reference=reference, exempt_rows=exempt_rows, excluded=excluded)

    def to_json(self) -> dict[str, object]:
        """The entry as it is written."""
        return {
            "ruling": self.ruling,
            "cases": list(self.cases),
            "reference": dict(self.reference),
            "exempt_rows": {case: list(names) for case, names in self.exempt_rows.items()},
            "excluded": dict(self.excluded),
        }


def _require_case_key(case: str, ruling: str) -> None:
    """Refuse a string that is not a 'cap-id/generator/pocket' key.

    Raises:
        MalformedAdmissionError: *case* is not a case key.
    """
    if not is_case_key(case):
        raise MalformedAdmissionError(f"admission {ruling} names {case!r}, which is not a 'cap-id/generator/pocket' case key.")


def _unique_strings(value: object, where: str) -> tuple[str, ...]:
    """A JSON list of strings with no duplicates (possibly empty).

    Raises:
        MalformedAdmissionError: *value* is not such a list.
    """
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        raise MalformedAdmissionError(f"{where} must be a list of strings, got {value!r}.")
    duplicates = sorted({item for item in value if value.count(item) > 1})
    if duplicates:
        raise MalformedAdmissionError(f"{where} lists {duplicates} more than once.")
    return tuple(value)


def _reference(value: object, cases: tuple[str, ...], ruling: str) -> dict[str, str | None]:
    """The reference map: exactly the admitted generators, each to another generator or null.

    Raises:
        MalformedAdmissionError: The map is not of that shape.
    """
    if not isinstance(value, dict):
        raise MalformedAdmissionError(f"admission {ruling} 'reference' must be an object, got {value!r}.")
    generators = {case_generator(case) for case in cases}
    missing = sorted(generators - set(value))
    extra = sorted(set(value) - generators)
    if missing or extra:
        raise MalformedAdmissionError(f"admission {ruling} 'reference' must name exactly the admitted generators {sorted(generators)}; missing {missing}, extra {extra}.")
    for generator, target in value.items():
        if target is None:
            continue
        if not isinstance(target, str) or not target or CASE_KEY_SEPARATOR in target:
            raise MalformedAdmissionError(f"admission {ruling} reference of {generator!r} must be a generator name or null, got {target!r}.")
        if target == generator:
            raise MalformedAdmissionError(f"admission {ruling} references generator {generator!r} to itself.")
    return dict(value)


def _exempt_rows(value: object, cases: tuple[str, ...], ruling: str) -> dict[str, tuple[str, ...]]:
    """The exemption map: admitted cases to non-empty lists of stage-0 all-case row criteria.

    Raises:
        MalformedAdmissionError: The map is not of that shape.
    """
    if not isinstance(value, dict):
        raise MalformedAdmissionError(f"admission {ruling} 'exempt_rows' must be an object, got {value!r}.")
    allowed = exemptable_criteria()
    exempt: dict[str, tuple[str, ...]] = {}
    for case, names in value.items():
        if case not in cases:
            raise MalformedAdmissionError(f"admission {ruling} exempts {case!r}, which it does not admit.")
        criteria = _unique_strings(names, f"admission {ruling} exempt_rows of {case!r}")
        if not criteria:
            raise MalformedAdmissionError(f"admission {ruling} exempt_rows of {case!r} must be non-empty.")
        unknown = [name for name in criteria if name not in allowed]
        if unknown:
            raise MalformedAdmissionError(f"admission {ruling} exempts {case!r} from {unknown}, which are not {EXEMPTABLE_STAGE} all-case rows {sorted(allowed)}.")
        exempt[case] = criteria
    return exempt


def _excluded(value: object, cases: tuple[str, ...], ruling: str) -> dict[str, str]:
    """The exclusion map: case keys not admitted here to exception class names.

    Raises:
        MalformedAdmissionError: The map is not of that shape.
    """
    if not isinstance(value, dict):
        raise MalformedAdmissionError(f"admission {ruling} 'excluded' must be an object, got {value!r}.")
    for case, name in value.items():
        _require_case_key(case, ruling)
        if case in cases:
            raise MalformedAdmissionError(f"admission {ruling} both admits and excludes {case!r}.")
        if not isinstance(name, str) or not EXCEPTION_CLASS_NAME.fullmatch(name):
            raise MalformedAdmissionError(f"admission {ruling} excludes {case!r} under {name!r}, which is not a qualified exception class name ('module.Class').")
    return dict(value)
