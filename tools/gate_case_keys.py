"""The case key of the gate vectors: ``cap-id/generator/pocket``."""

from __future__ import annotations

# Separator of the three components of a case key: cap id, generator, pocket.
CASE_KEY_SEPARATOR = "/"
CASE_KEY_PARTS = 3


def is_case_key(key: object) -> bool:
    """Whether *key* is a string of exactly the three case-key components."""
    return isinstance(key, str) and len(key.split(CASE_KEY_SEPARATOR)) == CASE_KEY_PARTS


def case_generator(case: str) -> str:
    """The generator component of a case key."""
    return case.split(CASE_KEY_SEPARATOR)[1]


def cap_and_pocket(case: str) -> tuple[str, str]:
    """The cap id and pocket components of a case key."""
    cap_id, _generator, pocket = case.split(CASE_KEY_SEPARATOR)
    return cap_id, pocket
