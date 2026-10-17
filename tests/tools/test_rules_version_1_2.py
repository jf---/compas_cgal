"""The "1.2" rules bump (Ruling 1.2-D2): the measured vectors carry rules 1.2, re-baselined from 0.3."""

from __future__ import annotations

from tools.gate_vectors import RULES_VERSION

# The bump's parent rules version: the re-baseline record's ``from``
# (`build/task12/g2-rebaseline-draft.json`).
PARENT_RULES_VERSION = "0.3"
BUMPED_RULES_VERSION = "1.2"


def test_the_gate_rules_are_bumped_to_1_2_from_0_3() -> None:
    assert RULES_VERSION == BUMPED_RULES_VERSION
    assert RULES_VERSION != PARENT_RULES_VERSION
