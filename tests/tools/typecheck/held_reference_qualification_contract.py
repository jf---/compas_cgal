from pathlib import Path
from typing import assert_type

from benchmarks.spec import PocketSpec
from tools.held_reference_qualification import QualificationOutcome
from tools.held_reference_qualification import Seconds
from tools.held_reference_qualification import qualify_cases
from tools.held_reference_qualification import render_qualification_markdown
from tools.held_reference_qualification import write_qualification_report


class _GeneratedPath:
    @property
    def operations(self) -> tuple[object, ...]:
        return (object(),)


def _generator(_spec: PocketSpec) -> _GeneratedPath:
    return _GeneratedPath()


outcomes = qualify_cases(_generator)
assert_type(outcomes, tuple[QualificationOutcome, ...])
assert_type(outcomes[0].generation_seconds, Seconds)
assert_type(render_qualification_markdown(outcomes), str)
assert_type(write_qualification_report(outcomes), Path)
