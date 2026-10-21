"""The Held path digest record: canonical serialisation, schema round trip, the --check comparison and the worker protocol."""

from __future__ import annotations

import hashlib
import json
from collections.abc import Callable
from dataclasses import dataclass
from dataclasses import replace
from typing import Any

import pytest

from tools.held_path_digest import DIGEST_WORKLOADS
from tools.held_path_digest import GENERATED_BY
from tools.held_path_digest import WORKER_MARKER
from tools.held_path_digest import DigestMismatchError
from tools.held_path_digest import HeldPathDigests
from tools.held_path_digest import MalformedDigestRecordError
from tools.held_path_digest import MissingWorkerResultError
from tools.held_path_digest import NonCanonicalReprError
from tools.held_path_digest import WorkloadDigest
from tools.held_path_digest import WorkloadMeasurement
from tools.held_path_digest import canonical_digest
from tools.held_path_digest import mismatches
from tools.held_path_digest import parse_worker_output
from tools.held_path_digest import require_identical


@dataclass(frozen=True)
class _Point:
    x: float
    y: float


class _AddressRepr:
    """Default ``object.__repr__``: carries the memory address, so it is not canonical."""


def _digest(seed: str) -> WorkloadDigest:
    sha = hashlib.sha256(seed.encode()).hexdigest()
    return WorkloadDigest(initial=sha, circles=sha, connectors=sha, circle_count=3, connector_count=3, seam="connector")


def _record() -> HeldPathDigests:
    return HeldPathDigests({name: _digest(name) for name in DIGEST_WORKLOADS})


def test_canonical_digest_is_the_sha256_of_the_repr() -> None:
    value = (_Point(0.1, 2.0), _Point(-0.0, 1e-300))
    assert canonical_digest(value) == hashlib.sha256(repr(value).encode()).hexdigest()


def test_canonical_digest_separates_values_that_differ_in_the_last_bit() -> None:
    assert canonical_digest(_Point(0.1, 0.2)) != canonical_digest(_Point(0.1, 0.2 + 2.0**-55))


def test_a_repr_carrying_a_memory_address_is_refused() -> None:
    with pytest.raises(NonCanonicalReprError):
        canonical_digest((_Point(0.0, 0.0), _AddressRepr()))


def test_the_record_round_trips_through_json_loads() -> None:
    record = _record()
    text = record.to_json()
    raw = json.loads(text)
    assert raw["generated_by"] == GENERATED_BY
    assert set(raw["workloads"]) == set(DIGEST_WORKLOADS)
    assert HeldPathDigests.from_json(text) == record


@pytest.mark.parametrize(
    "corrupt",
    [
        lambda raw: raw.update(schema="held-path-digests/0"),
        lambda raw: raw.pop("workloads"),
        lambda raw: raw["workloads"]["figure5_refined"].update(initial="not-hex"),
        lambda raw: raw["workloads"]["figure5_refined"].update(circle_count="3"),
        lambda raw: raw["workloads"]["figure5_refined"].pop("seam"),
        lambda raw: raw["workloads"]["figure5_refined"].update(extra=1),
    ],
    ids=["schema", "no-workloads", "digest-not-hex", "count-not-int", "missing-field", "extra-field"],
)
def test_a_malformed_record_fails_loud(corrupt: Callable[[dict[str, Any]], object]) -> None:
    raw = json.loads(_record().to_json())
    corrupt(raw)
    with pytest.raises(MalformedDigestRecordError):
        HeldPathDigests.from_json(json.dumps(raw))


def test_identical_records_have_no_mismatch() -> None:
    assert mismatches(_record(), _record()) == []
    require_identical(_record(), _record())


def test_one_changed_field_is_named_with_its_workload() -> None:
    recorded = _record()
    current = HeldPathDigests({**recorded.workloads, "figure5_refined": replace(recorded.workloads["figure5_refined"], connectors="0" * 64)})
    found = mismatches(recorded, current)
    assert len(found) == 1
    assert "figure5_refined" in found[0]
    assert "connectors" in found[0]
    with pytest.raises(DigestMismatchError, match="figure5_refined"):
        require_identical(recorded, current)


def test_a_workload_missing_from_the_current_run_is_a_mismatch() -> None:
    recorded = _record()
    current = HeldPathDigests({name: digest for name, digest in recorded.workloads.items() if name != "figure8_upper_boundary"})
    assert any("figure8_upper_boundary" in line for line in mismatches(recorded, current))


def test_checking_a_subset_compares_only_that_subset() -> None:
    recorded = _record()
    current = HeldPathDigests({"figure5_refined": recorded.workloads["figure5_refined"]})
    assert mismatches(recorded, current, workloads=("figure5_refined",)) == []


def test_the_worker_result_is_read_from_its_marker_line_among_other_output() -> None:
    measurement = WorkloadMeasurement("figure5_refined", _digest("a"), 1.25, {"guide": 0.25, "initial": 1.0})
    output = f"some generator print\n{WORKER_MARKER}{measurement.to_json_line()}\ntrailing\n"
    assert parse_worker_output(output) == measurement


def test_worker_output_without_a_marker_line_fails_loud() -> None:
    with pytest.raises(MissingWorkerResultError):
        parse_worker_output("generator crashed before reporting\n")
