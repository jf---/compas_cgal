from typing import assert_type

from benchmarks.held_reference_cases import Degree
from benchmarks.held_reference_cases import Figure7Observation
from benchmarks.held_reference_cases import HeldReferenceCase
from benchmarks.held_reference_cases import load_all_held_reference_cases
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.spec import PocketSpec

case = assert_type(load_held_reference_case("figure5"), HeldReferenceCase)
assert_type(load_all_held_reference_cases(), tuple[HeldReferenceCase, ...])
assert_type(case.pocket_spec(), PocketSpec)
assert_type(case.tea_cap, Degree)
assert_type(case.figure7_observation, Figure7Observation | None)
