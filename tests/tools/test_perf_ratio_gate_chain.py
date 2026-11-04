"""The ratio gate's anchor chain: links, the reset rule, the replay, and their named refusals."""

from __future__ import annotations

import dataclasses
import json
from collections.abc import Callable
from collections.abc import Mapping
from collections.abc import Sequence
from fractions import Fraction
from pathlib import Path
from typing import Any

import pytest
from tools.perf_history_store import REPOSITORY as REPOSITORY_ROOT
from tools.perf_ratio_gate_chain import AnchorChainBrokenError
from tools.perf_ratio_gate_chain import AnchorChainEndError
from tools.perf_ratio_gate_chain import AnchorChainMismatchError
from tools.perf_ratio_gate_chain import AnchorLink
from tools.perf_ratio_gate_chain import AnchorLinkSelfLoopError
from tools.perf_ratio_gate_chain import AnchorRecordDigestError
from tools.perf_ratio_gate_chain import AnchorRecordMissingError
from tools.perf_ratio_gate_chain import AnchorRecordOutsideRootError
from tools.perf_ratio_gate_chain import AnchorRecordUnloadableError
from tools.perf_ratio_gate_chain import AnchorValueMismatchError
from tools.perf_ratio_gate_chain import DuplicateAnchorRecordError
from tools.perf_ratio_gate_chain import ExtensionPairCapError
from tools.perf_ratio_gate_chain import InvalidAnchorSessionError
from tools.perf_ratio_gate_chain import RegressedAnchorLinkError
from tools.perf_ratio_gate_chain import UnimprovedAnchorLinkError
from tools.perf_ratio_gate_chain import UnrequestedExtensionError
from tools.perf_ratio_gate_chain import carry_anchor
from tools.perf_ratio_gate_chain import derive_anchors
from tools.perf_ratio_gate_citations import SAFETY
from tools.perf_ratio_gate_citations import CalibrationRecordMissingError
from tools.perf_ratio_gate_citations import CalibrationRecordOutsideRootError
from tools.perf_ratio_gate_citations import InvalidCalibrationSessionError
from tools.perf_ratio_gate_citations import MarginMismatchError
from tools.perf_ratio_gate_decide import MAX_PAIRS
from tools.perf_ratio_gate_decide import ExactRatio
from tools.perf_ratio_gate_merge import ExtensionMismatchError
from tools.perf_ratio_gate_registry import PROMOTION_ROOT
from tools.perf_ratio_gate_registry import Anchor
from tools.perf_ratio_gate_registry import GateRegistry
from tools.perf_ratio_gate_registry import MalformedRegistryError
from tools.perf_ratio_gate_registry import Margin
from tools.perf_ratio_gate_registry import RegistryDigestMismatchError
from tools.perf_ratio_gate_registry import SupersededRegistrySchemaError
from tools.perf_ratio_gate_registry import build_registry
from tools.perf_ratio_gate_registry import load_registry
from tools.perf_ratio_gate_registry import promotion_dir
from tools.perf_ratio_gate_registry import write_registry
from tools.perf_ratio_record import ArmName
from tools.perf_ratio_record import CpuSeconds
from tools.perf_ratio_record import Injection
from tools.perf_ratio_record import InconsistentRatioRecordError
from tools.perf_ratio_record import MalformedRatioRecordError
from tools.perf_ratio_record import PerfRatioSession
from tools.run_provenance import sha256_file

from tests.tools.perf_ratio_sessions import HARNESS
from tests.tools.perf_ratio_sessions import REF
from tests.tools.perf_ratio_sessions import _dirty
from tests.tools.perf_ratio_sessions import _dump_unchecked
from tests.tools.perf_ratio_sessions import _registry
from tests.tools.perf_ratio_sessions import _rewrite
from tests.tools.perf_ratio_sessions import cite
from tests.tools.perf_ratio_sessions import link_calibration
from tests.tools.perf_ratio_sessions import promotion_link
from tests.tools.perf_ratio_sessions import synthetic_session
from tests.tools.perf_ratio_sessions import write_promotion
from tests.tools.perf_ratio_sessions import write_record


# --- the anchor chain (spec amendment 2026-11-01) ---------------------------------------------------


# Commits of a synthetic promotion chain REF -> REV_B -> REV_C (-> REV_D).
REV_B = "b0" * 20
REV_C = "c0" * 20
REV_D = "d0" * 20
ONE = Fraction(1)
# One link's pairs: "w" drifts to 1.02 (carried under m = 0.03), "v" improves to 0.95 (a certified reset).
CARRY_RESET = {"w": ((100.0, 102.0),), "v": ((100.0, 95.0),)}
# A link after CARRY_RESET: "w" holds at 1.0 (A_w * rho_w = 1.02, carried), "v" improves again to 0.95.
HOLD_RESET = {"w": ((100.0, 100.0),), "v": ((100.0, 95.0),)}
# A gate session whose "w" is red (1.05 > 1 + 0.03): the gate extends it; "v" is the link's certified improvement.
RED_W = {"w": ((100.0, 105.0),), "v": ((100.0, 95.0),)}
# A harness digest other than HARNESS: the harness a chain's earlier links ran before a rebase.
OLD_HARNESS = "g" * 64


def _chain(root: Path, links: Sequence[AnchorLink], workloads: Sequence[str] = ("v", "w"), anchor_ref: str = REF) -> GateRegistry:
    ref_commit = links[-1].to_ref if links else anchor_ref
    return build_registry(ref_commit=ref_commit, harness_sha256=HARNESS, calibrations=(), workloads=workloads, anchor_ref=anchor_ref, anchor_links=links, root=root)


def _round_trip(root: Path, registry: GateRegistry) -> GateRegistry:
    path = root / "perf_ratio_gate.json"
    write_registry(registry, path, root=root)
    return load_registry(path, root=root)


def _raw_chain(links: Sequence[AnchorLink], workloads: tuple[str, ...] = ("v", "w"), *, ref_commit: str = REV_B, anchor_ref: str = REF) -> GateRegistry:
    """An unchecked registry over *links*, every stored anchor 1, as a hand edit would write it."""
    return GateRegistry(
        ref_commit=ref_commit,
        harness_sha256=HARNESS,
        safety=SAFETY,
        workloads=workloads,
        calibrations=(),
        anchor_ref=anchor_ref,
        anchor_links=tuple(links),
        exact_anchors={name: Anchor(ONE) for name in workloads},
    )


def _ratio(pairs: Sequence[tuple[float, float]]) -> Fraction:
    return Fraction(min(head for _, head in pairs) / min(ref for ref, _ in pairs))


def test_an_empty_chain_anchors_every_workload_at_one(tmp_path: Path) -> None:
    registry = _registry(tmp_path)
    assert registry.anchor_ref == REF
    assert registry.anchor_links == ()
    assert registry.anchors_exact() == {"v": ONE, "w": ONE}
    loaded = _round_trip(tmp_path, registry)
    assert loaded.anchors_exact() == {"v": ONE, "w": ONE}
    document = json.loads((tmp_path / "perf_ratio_gate.json").read_text(encoding="utf-8"))
    assert document["schema"] == "perf-ratio-gate/2"
    assert document["anchor_ref"] == REF
    assert document["anchor_links"] == []
    assert document["anchors"] == {"v": 1.0, "w": 1.0}


def test_a_carried_link_multiplies_the_anchor_by_rho(tmp_path: Path) -> None:
    link = promotion_link(tmp_path, REF, REV_B, CARRY_RESET)
    margin = link.calibration.margins["w"]
    rho = _ratio(CARRY_RESET["w"])
    assert margin == pytest.approx(0.03)
    assert 1 - Fraction(margin) <= rho <= 1 + Fraction(margin)
    loaded = _round_trip(tmp_path, _chain(tmp_path, [link]))
    assert loaded.anchors_exact()["w"] == rho == Fraction(102.0 / 100.0)
    assert loaded.anchors["w"] == 1.02


def test_a_certified_improvement_resets_the_anchor_to_one(tmp_path: Path) -> None:
    link = promotion_link(tmp_path, REF, REV_B, CARRY_RESET)
    assert _ratio(CARRY_RESET["v"]) < 1 - Fraction(link.calibration.margins["v"])
    loaded = _round_trip(tmp_path, _chain(tmp_path, [link]))
    assert loaded.anchors_exact()["v"] == ONE


def test_a_carried_drift_is_reset_by_a_later_certified_improvement(tmp_path: Path) -> None:
    first = promotion_link(tmp_path, REF, REV_B, CARRY_RESET)
    second = promotion_link(tmp_path, REV_B, REV_C, {"w": ((100.0, 95.0),), "v": ((100.0, 95.0),)})
    product = Fraction(1.02) * Fraction(0.95)
    assert product < 1 - Fraction(second.calibration.margins["w"])
    loaded = _round_trip(tmp_path, _chain(tmp_path, [first, second]))
    assert loaded.anchors_exact() == {"v": ONE, "w": ONE}


def test_the_rule_decides_on_the_anchored_product_not_the_link_ratio(tmp_path: Path) -> None:
    # rho = 0.96 alone is a certified improvement (< 1 - 0.03); A * rho = 1.02 * 0.96 = 0.9792 is not.
    first = promotion_link(tmp_path, REF, REV_B, CARRY_RESET)
    second = promotion_link(tmp_path, REV_B, REV_C, {"w": ((100.0, 96.0),), "v": ((100.0, 95.0),)})
    margin = Fraction(second.calibration.margins["w"])
    product = Fraction(102.0 / 100.0) * Fraction(96.0 / 100.0)
    assert Fraction(0.96) < 1 - margin <= product
    loaded = _round_trip(tmp_path, _chain(tmp_path, [first, second]))
    assert loaded.anchors_exact()["w"] == product
    # The exact anchor is re-derived from the chain, never read back from its stored reporting float.
    assert Fraction(float(product)) != product
    assert loaded.anchors["w"] == float(product)


def test_a_link_ratio_replays_its_extensions(tmp_path: Path) -> None:
    # The gate session alone is red (1.05 > 1.03), still red after the first extension (min head 104),
    # and the second brings min head down to 101: each extension targets a workload red at that point.
    extensions = [("w", ((100.0, 104.0), (100.0, 106.0))), ("w", ((100.0, 101.0),))]
    link = promotion_link(tmp_path, REF, REV_B, {"w": ((100.0, 105.0),), "v": ((100.0, 95.0),)}, extensions)
    assert len(link.sessions) == 3
    loaded = _round_trip(tmp_path, _chain(tmp_path, [link]))
    assert loaded.anchors_exact()["w"] == Fraction(101.0 / 100.0)


def test_workloads_entering_and_leaving_the_chain(tmp_path: Path) -> None:
    # x is carried at the first link, then left out of the second: dropped, and re-added later at 1.
    # v first appears at the second link (carried from 1); u was never measured by a promotion.
    # y is each link's certified improvement (0.95 < 1 - 0.03), so both links are promotable; w carries
    # 1.02 then 1.02 * 1.005 = 1.0251 <= 1.03.
    first = promotion_link(tmp_path, REF, REV_B, {"w": ((100.0, 102.0),), "x": ((100.0, 102.0),), "y": ((100.0, 95.0),)})
    second = promotion_link(tmp_path, REV_B, REV_C, {"w": ((100.0, 100.5),), "v": ((100.0, 101.0),), "y": ((100.0, 95.0),)})
    anchors = derive_anchors([first, second], ("u", "v", "w", "x"), root=tmp_path)
    assert anchors == {"u": ONE, "v": Fraction(1.01), "w": Fraction(1.02) * Fraction(1.005), "x": ONE}
    assert derive_anchors([first], ("w",), root=tmp_path) == {"w": Fraction(1.02)}
    loaded = _round_trip(tmp_path, _chain(tmp_path, [first, second], workloads=("u", "v", "w", "x")))
    assert loaded.anchors_exact() == anchors


def test_carry_anchor_is_the_reset_rule_in_exact_arithmetic() -> None:
    margin = Margin(0.03)
    assert carry_anchor(Anchor(ONE), ExactRatio(Fraction(1.02)), margin) == Fraction(1.02)
    assert carry_anchor(Anchor(ONE), ExactRatio(Fraction(0.95)), margin) == ONE
    boundary = ExactRatio(Fraction(float(1 - Fraction(margin))))
    # The boundary float sits within half an ulp of 1 - m: the branch follows the exact comparison.
    expected = Fraction(boundary) if Fraction(boundary) >= 1 - Fraction(margin) else ONE
    assert carry_anchor(Anchor(ONE), boundary, margin) == expected


def test_every_chain_refusal_is_an_anchor_chain_mismatch() -> None:
    for error in (
        AnchorRecordMissingError,
        AnchorRecordDigestError,
        AnchorRecordOutsideRootError,
        AnchorRecordUnloadableError,
        AnchorChainBrokenError,
        AnchorChainEndError,
        AnchorValueMismatchError,
        InvalidAnchorSessionError,
        AnchorLinkSelfLoopError,
        DuplicateAnchorRecordError,
        RegressedAnchorLinkError,
        UnimprovedAnchorLinkError,
        UnrequestedExtensionError,
        ExtensionPairCapError,
    ):
        assert issubclass(error, AnchorChainMismatchError)
    assert issubclass(SupersededRegistrySchemaError, MalformedRegistryError)


def test_promotion_dir_names_both_references() -> None:
    assert promotion_dir(REF, REV_B) == PROMOTION_ROOT / "15fbed70-b0b0b0b0"


def _chained_file(root: Path) -> Path:
    links = [promotion_link(root, REF, REV_B, CARRY_RESET), promotion_link(root, REV_B, REV_C, HOLD_RESET)]
    path = root / "perf_ratio_gate.json"
    write_registry(_chain(root, links), path, root=root)
    return path


def test_a_missing_link_record_is_refused(tmp_path: Path) -> None:
    path = _chained_file(tmp_path)
    (tmp_path / promotion_dir(REV_B, REV_C).relative_to(REPOSITORY_ROOT) / "gate" / "perf_ratio_session.json").unlink()
    with pytest.raises(AnchorRecordMissingError):
        load_registry(path, root=tmp_path)


def test_a_rewritten_link_record_is_refused(tmp_path: Path) -> None:
    path = _chained_file(tmp_path)
    record = tmp_path / promotion_dir(REF, REV_B).relative_to(REPOSITORY_ROOT) / "gate" / "perf_ratio_session.json"
    write_record(synthetic_session({"w": ((100.0, 101.0),), "v": ((100.0, 95.0),)}, head_rev=REV_B), record)
    with pytest.raises(AnchorRecordDigestError):
        load_registry(path, root=tmp_path)


def test_a_link_record_outside_the_root_is_refused(tmp_path: Path) -> None:
    path = _chained_file(tmp_path)

    def escape(document: dict[str, Any]) -> None:
        document["anchor_links"][0]["sessions"][0]["path"] = "../escaped/perf_ratio_session.json"

    _rewrite(path, escape, redigest=True)
    with pytest.raises(AnchorRecordOutsideRootError):
        load_registry(path, root=tmp_path)


def _inconsistent(data: bytes) -> bytes:
    # The record's first stored workload ratio doubled: valid JSON, but its derived value no longer recomputes.
    raw = json.loads(data.decode("utf-8"))
    raw["workloads"][0]["ratio"] *= 2
    return json.dumps(raw).encode("utf-8")


# Rewrites of a cited record that keep it citable but not loadable, each with the loader failure it causes.
UNLOADABLE_RECORDS: dict[str, tuple[Callable[[bytes], bytes], type[Exception]]] = {
    "not-json": (lambda data: b"not json", MalformedRatioRecordError),
    "not-utf-8": (lambda data: b"\xff\xfe{}", UnicodeDecodeError),
    "inconsistent": (_inconsistent, InconsistentRatioRecordError),
}


@pytest.mark.parametrize("content", sorted(UNLOADABLE_RECORDS))
def test_a_digest_pinned_link_record_that_does_not_load_is_refused(tmp_path: Path, content: str) -> None:
    # The record's digest is re-pinned (and the registry re-digested): only the load can fail. That is a refusal of
    # the cited evidence, not an infra failure, so it is an AnchorChainMismatchError with the loader's error chained.
    path = _chained_file(tmp_path)
    record = tmp_path / promotion_dir(REF, REV_B).relative_to(REPOSITORY_ROOT) / "gate" / "perf_ratio_session.json"
    rewrite, cause = UNLOADABLE_RECORDS[content]
    record.write_bytes(rewrite(record.read_bytes()))

    def repin(document: dict[str, Any]) -> None:
        document["anchor_links"][0]["sessions"][0]["record_sha256"] = sha256_file(record)

    _rewrite(path, repin, redigest=True)
    with pytest.raises(AnchorRecordUnloadableError) as raised:
        load_registry(path, root=tmp_path)
    assert isinstance(raised.value, AnchorChainMismatchError)
    assert isinstance(raised.value.__cause__, cause)
    assert "has its cited digest but does not load" in str(raised.value)


# The ruling an unloadable cited record's refusal states (CI-N3 follow-up), named independently of the gate modules.
UNLOADABLE_RULING = (
    "The defect is in the harness change that broke loading: that change must keep committed records loadable or ship a migration; the gate has no recovery path of its own."
)


def test_an_unloadable_link_record_names_its_path_its_cause_and_the_ruling(tmp_path: Path) -> None:
    path = _chained_file(tmp_path)
    cited = promotion_dir(REF, REV_B).relative_to(REPOSITORY_ROOT) / "gate" / "perf_ratio_session.json"
    record = tmp_path / cited
    record.write_bytes(_inconsistent(record.read_bytes()))

    def repin(document: dict[str, Any]) -> None:
        document["anchor_links"][0]["sessions"][0]["record_sha256"] = sha256_file(record)

    _rewrite(path, repin, redigest=True)
    with pytest.raises(AnchorRecordUnloadableError) as raised:
        load_registry(path, root=tmp_path)
    message = str(raised.value)
    assert f"cited record {cited.as_posix()} " in message
    assert "InconsistentRatioRecordError: " in message
    assert UNLOADABLE_RULING in message


def test_a_stored_anchor_that_disagrees_is_refused(tmp_path: Path) -> None:
    path = _chained_file(tmp_path)

    def creep(document: dict[str, Any]) -> None:
        document["anchors"]["w"] = 1.0

    _rewrite(path, creep, redigest=True)
    with pytest.raises(AnchorValueMismatchError):
        load_registry(path, root=tmp_path)


def test_a_stored_anchor_set_that_misses_a_workload_is_refused(tmp_path: Path) -> None:
    path = _chained_file(tmp_path)
    _rewrite(path, lambda document: document["anchors"].pop("v"), redigest=True)
    with pytest.raises(AnchorValueMismatchError):
        load_registry(path, root=tmp_path)


def test_write_registry_refuses_a_stored_anchor_that_disagrees(tmp_path: Path) -> None:
    link = promotion_link(tmp_path, REF, REV_B, CARRY_RESET)
    path = tmp_path / "perf_ratio_gate.json"
    with pytest.raises(AnchorValueMismatchError):
        write_registry(_raw_chain([link]), path, root=tmp_path)
    assert not path.exists()


def test_a_link_that_does_not_start_where_the_previous_ended_is_refused(tmp_path: Path) -> None:
    links = [promotion_link(tmp_path, REF, REV_B, CARRY_RESET), promotion_link(tmp_path, REV_C, REV_D, CARRY_RESET)]
    with pytest.raises(AnchorChainBrokenError):
        _chain(tmp_path, links)
    path = _dump_unchecked(_raw_chain(links, ref_commit=REV_D), tmp_path / "perf_ratio_gate.json")
    with pytest.raises(AnchorChainBrokenError):
        load_registry(path, root=tmp_path)


def test_a_first_link_that_does_not_start_at_the_anchor_ref_is_refused(tmp_path: Path) -> None:
    links = [promotion_link(tmp_path, REF, REV_B, CARRY_RESET)]
    with pytest.raises(AnchorChainBrokenError):
        _chain(tmp_path, links, anchor_ref=REV_D)
    path = _dump_unchecked(_raw_chain(links, anchor_ref=REV_D), tmp_path / "perf_ratio_gate.json")
    with pytest.raises(AnchorChainBrokenError):
        load_registry(path, root=tmp_path)


def test_a_chain_that_does_not_end_at_the_reference_is_refused(tmp_path: Path) -> None:
    links = [promotion_link(tmp_path, REF, REV_B, CARRY_RESET)]
    with pytest.raises(AnchorChainEndError):
        build_registry(ref_commit=REV_C, harness_sha256=HARNESS, calibrations=(), workloads=("v", "w"), anchor_ref=REF, anchor_links=links, root=tmp_path)
    path = _dump_unchecked(_raw_chain(links, ref_commit=REV_C), tmp_path / "perf_ratio_gate.json")
    with pytest.raises(AnchorChainEndError):
        load_registry(path, root=tmp_path)


def test_an_empty_chain_anchored_elsewhere_is_refused(tmp_path: Path) -> None:
    with pytest.raises(AnchorChainEndError):
        build_registry(ref_commit=REF, harness_sha256=HARNESS, calibrations=(), workloads=("v", "w"), anchor_ref=REV_B)
    path = _dump_unchecked(_raw_chain([], ref_commit=REF, anchor_ref=REV_B), tmp_path / "perf_ratio_gate.json")
    with pytest.raises(AnchorChainEndError):
        load_registry(path, root=tmp_path)


def _gate_variant(**overrides: object) -> Callable[[], tuple[PerfRatioSession, ...]]:
    def sessions() -> tuple[PerfRatioSession, ...]:
        options: dict[str, object] = {"ref_commit": REF, "head_rev": REV_B}
        options.update(overrides)
        return (synthetic_session(CARRY_RESET, **options),)  # type: ignore[arg-type]

    return sessions


def _dirty_gate() -> tuple[PerfRatioSession, ...]:
    return (_dirty(synthetic_session(CARRY_RESET, head_rev=REV_B)),)


def _extension_of_an_unmeasured_workload() -> tuple[PerfRatioSession, ...]:
    return (synthetic_session(CARRY_RESET, head_rev=REV_B), synthetic_session({"u": ((100.0, 100.0),)}, head_rev=REV_B))


def _extension_of_two_workloads() -> tuple[PerfRatioSession, ...]:
    # Pairs other than the gate session's, so the extension is not a byte copy of it (a duplicate citation).
    return (synthetic_session(CARRY_RESET, head_rev=REV_B), synthetic_session({"w": ((100.0, 101.0),), "v": ((100.0, 96.0),)}, head_rev=REV_B))


# Each way a link's cited sessions fail to be the honest record of promoting REF -> REV_B.
INVALID_LINK_SESSIONS: Mapping[str, Callable[[], tuple[PerfRatioSession, ...]]] = {
    "none": lambda: (),
    "injected": _gate_variant(injected=Injection(arm=ArmName.HEAD, seconds=CpuSeconds(1.0))),
    "ref-elsewhere": _gate_variant(ref_commit=REV_C),
    "head-elsewhere": _gate_variant(head_rev=REV_C),
    "dirty": _dirty_gate,
    "cpu-differs": _gate_variant(cpu_model="Intel Xeon"),
    "extension-of-an-unmeasured-workload": _extension_of_an_unmeasured_workload,
    "extension-of-two-workloads": _extension_of_two_workloads,
}
# Each case's named refusal reason (a `match=` regex), so a fixture change cannot move a case to another check;
# a case missing here fails its test with a KeyError.
INVALID_LINK_REFUSALS: Mapping[str, str] = {
    "none": r"cites no session; a link cites its promoted gate session first",
    "injected": r"the gate session is an audit session \(injected 1\.0 s on head\)",
    "ref-elsewhere": rf"the gate session measured against {REV_C}, not the link's from_ref",
    "head-elsewhere": rf"the gate session's HEAD arm ran at {REV_C}, not the link's to_ref",
    "dirty": r"the gate session's head arm tree was dirty",
    "cpu-differs": r"the gate session ran on 'Intel Xeon', not the link calibration's",
    "extension-of-an-unmeasured-workload": r"measures \['u'\], not one workload of the gate session",
    "extension-of-two-workloads": r"measures \['w', 'v'\], not one workload of the gate session",
}


@pytest.mark.parametrize("reason", sorted(INVALID_LINK_SESSIONS))
def test_a_link_citing_a_dishonest_promotion_is_refused(tmp_path: Path, reason: str) -> None:
    sessions = INVALID_LINK_SESSIONS[reason]()
    cited = write_promotion(tmp_path, REF, REV_B, sessions[0], sessions[1:]) if sessions else ()
    link = AnchorLink(from_ref=REF, to_ref=REV_B, sessions=cited, calibration=link_calibration(tmp_path, REF, ("v", "w")))
    with pytest.raises(InvalidAnchorSessionError, match=INVALID_LINK_REFUSALS[reason]):
        _chain(tmp_path, [link])
    path = _dump_unchecked(_raw_chain([link]), tmp_path / "perf_ratio_gate.json")
    with pytest.raises(InvalidAnchorSessionError, match=INVALID_LINK_REFUSALS[reason]):
        load_registry(path, root=tmp_path)


def test_a_link_extension_measured_elsewhere_is_refused(tmp_path: Path) -> None:
    # w is red (1.05 > 1.03), so its extension is requested; it ran another harness.
    gate = synthetic_session(RED_W, head_rev=REV_B)
    elsewhere = synthetic_session({"w": ((100.0, 101.0),)}, head_rev=REV_B, harness="g" * 64)
    link = AnchorLink(from_ref=REF, to_ref=REV_B, sessions=write_promotion(tmp_path, REF, REV_B, gate, [elsewhere]), calibration=link_calibration(tmp_path, REF, ("v", "w")))
    with pytest.raises(ExtensionMismatchError):
        _chain(tmp_path, [link])


def test_a_link_calibration_with_a_wrong_margin_is_refused(tmp_path: Path) -> None:
    link = promotion_link(tmp_path, REF, REV_B, CARRY_RESET)
    widened = dataclasses.replace(link.calibration, margins={name: Margin(margin * 2.0) for name, margin in link.calibration.margins.items()})
    with pytest.raises(MarginMismatchError):
        _chain(tmp_path, [dataclasses.replace(link, calibration=widened)])


def test_a_link_calibration_at_another_reference_is_refused(tmp_path: Path) -> None:
    link = promotion_link(tmp_path, REF, REV_B, CARRY_RESET)
    elsewhere = link_calibration(tmp_path, REV_C, ("v", "w"))
    with pytest.raises(InvalidCalibrationSessionError):
        _chain(tmp_path, [dataclasses.replace(link, calibration=elsewhere)])


def test_a_link_calibration_over_other_workloads_is_refused(tmp_path: Path) -> None:
    link = promotion_link(tmp_path, REF, REV_B, CARRY_RESET)
    narrower = link_calibration(tmp_path, REF, ("w",))
    with pytest.raises(InvalidCalibrationSessionError):
        _chain(tmp_path, [dataclasses.replace(link, calibration=narrower)])


# --- review fix wave: harness rebase, promotability on replay, self-loops, duplicate citations -----


def test_a_harness_rebase_keeps_the_chain(tmp_path: Path) -> None:
    # The link and its calibration ran under the old harness; the registry was rebased to HARNESS since.
    gate = synthetic_session(CARRY_RESET, head_rev=REV_B, harness=OLD_HARNESS)
    calibration = link_calibration(tmp_path, REF, ("v", "w"), harness=OLD_HARNESS)
    link = AnchorLink(REF, REV_B, write_promotion(tmp_path, REF, REV_B, gate), calibration)
    registry = build_registry(ref_commit=REV_B, harness_sha256=HARNESS, calibrations=(), workloads=("v", "w"), anchor_ref=REF, anchor_links=[link], root=tmp_path)
    loaded = _round_trip(tmp_path, registry)
    assert loaded.harness_sha256 == HARNESS
    assert loaded.anchors_exact() == {"v": ONE, "w": Fraction(1.02)}


def test_a_link_calibration_under_another_harness_than_its_gate_session_is_refused(tmp_path: Path) -> None:
    gate = synthetic_session(CARRY_RESET, head_rev=REV_B)
    calibration = link_calibration(tmp_path, REF, ("v", "w"), harness=OLD_HARNESS)
    link = AnchorLink(REF, REV_B, write_promotion(tmp_path, REF, REV_B, gate), calibration)
    with pytest.raises(InvalidCalibrationSessionError, match=OLD_HARNESS):
        _chain(tmp_path, [link])
    path = _dump_unchecked(_raw_chain([link]), tmp_path / "perf_ratio_gate.json")
    with pytest.raises(InvalidCalibrationSessionError):
        load_registry(path, root=tmp_path)


def test_a_link_calibration_harness_refusal_names_the_gate_sessions_harness(tmp_path: Path) -> None:
    # The registry's harness is HARNESS too, so only the owner the message names tells the two paths apart.
    gate = synthetic_session(CARRY_RESET, head_rev=REV_B)
    calibration = link_calibration(tmp_path, REF, ("v", "w"), harness=OLD_HARNESS)
    link = AnchorLink(REF, REV_B, write_promotion(tmp_path, REF, REV_B, gate), calibration)
    with pytest.raises(InvalidCalibrationSessionError, match=rf"not the gate session's {HARNESS}") as built:
        _chain(tmp_path, [link])
    assert "registry's" not in str(built.value)
    path = _dump_unchecked(_raw_chain([link]), tmp_path / "perf_ratio_gate.json")
    with pytest.raises(InvalidCalibrationSessionError, match=rf"not the gate session's {HARNESS}"):
        load_registry(path, root=tmp_path)


def _refused_both_ways(root: Path, links: Sequence[AnchorLink], error: type[Exception]) -> None:
    """*links* are refused with *error* when built, and when a hand-written registry over them is loaded."""
    with pytest.raises(error):
        _chain(root, links)
    path = _dump_unchecked(_raw_chain(links, ref_commit=links[-1].to_ref), root / "perf_ratio_gate.json")
    with pytest.raises(error):
        load_registry(path, root=root)


def test_a_link_that_left_a_workload_regressed_is_refused(tmp_path: Path) -> None:
    # w at 1.10 > 1 + 0.03 and never extended: promote would have refused this session (P6).
    link = promotion_link(tmp_path, REF, REV_B, {"w": ((100.0, 110.0),), "v": ((100.0, 95.0),)})
    _refused_both_ways(tmp_path, [link], RegressedAnchorLinkError)


def test_a_link_that_regressed_only_under_its_anchor_is_refused(tmp_path: Path) -> None:
    # rho_w = 1.02 alone passes; after the first link A_w = 1.02, and 1.02 * 1.02 > 1.03.
    first = promotion_link(tmp_path, REF, REV_B, CARRY_RESET)
    second = promotion_link(tmp_path, REV_B, REV_C, CARRY_RESET)
    _refused_both_ways(tmp_path, [first, second], RegressedAnchorLinkError)


def test_a_link_without_a_certified_improvement_is_refused(tmp_path: Path) -> None:
    link = promotion_link(tmp_path, REF, REV_B, {"w": ((100.0, 102.0),), "v": ((100.0, 100.0),)})
    _refused_both_ways(tmp_path, [link], UnimprovedAnchorLinkError)


def test_a_link_improved_only_without_its_anchor_is_refused(tmp_path: Path) -> None:
    # rho_w = 0.96 alone is IMPROVED (< 0.97); under A_w = 1.02 the product 0.9792 is not, and v holds at 1.
    first = promotion_link(tmp_path, REF, REV_B, CARRY_RESET)
    second = promotion_link(tmp_path, REV_B, REV_C, {"w": ((100.0, 96.0),), "v": ((100.0, 100.0),)})
    _refused_both_ways(tmp_path, [first, second], UnimprovedAnchorLinkError)


def test_an_extension_of_a_workload_that_was_not_red_is_refused(tmp_path: Path) -> None:
    # v at 0.975 is PASS, so the gate never extends it; a cherry-picked extension at 0.96 would reset A_v (P7).
    link = promotion_link(tmp_path, REF, REV_B, {"w": ((100.0, 95.0),), "v": ((100.0, 97.5),)}, [("v", ((100.0, 96.0),))])
    _refused_both_ways(tmp_path, [link], UnrequestedExtensionError)


def test_an_extension_after_the_workload_turned_green_is_refused(tmp_path: Path) -> None:
    # The first extension brings w to 1.01 (PASS); the gate stops there, so a second one was never requested.
    link = promotion_link(tmp_path, REF, REV_B, RED_W, [("w", ((100.0, 101.0),)), ("w", ((100.0, 99.0),))])
    _refused_both_ways(tmp_path, [link], UnrequestedExtensionError)


def test_an_extension_red_only_under_its_anchor_is_replayed(tmp_path: Path) -> None:
    # rho_w = 1.02 alone is PASS; under A_w = 1.02 it is red (1.0404), so its extension to 1.0 was requested.
    first = promotion_link(tmp_path, REF, REV_B, CARRY_RESET)
    second = promotion_link(tmp_path, REV_B, REV_C, CARRY_RESET, [("w", ((100.0, 100.0),))])
    loaded = _round_trip(tmp_path, _chain(tmp_path, [first, second]))
    assert loaded.anchors_exact() == {"v": ONE, "w": Fraction(1.02)}


def test_a_link_whose_merged_pairs_exceed_max_pairs_is_refused(tmp_path: Path) -> None:
    # w holds MAX_PAIRS red pairs before its extension: the gate stops extending at MAX_PAIRS.
    gate = {"w": ((100.0, 105.0),) * MAX_PAIRS, "v": ((100.0, 95.0),)}
    link = promotion_link(tmp_path, REF, REV_B, gate, [("w", ((100.0, 101.0),))])
    _refused_both_ways(tmp_path, [link], ExtensionPairCapError)


def test_a_link_whose_merged_pairs_reach_max_pairs_is_replayed(tmp_path: Path) -> None:
    gate = {"w": ((100.0, 105.0),) * (MAX_PAIRS - 1), "v": ((100.0, 95.0),)}
    link = promotion_link(tmp_path, REF, REV_B, gate, [("w", ((100.0, 101.0),))])
    assert _round_trip(tmp_path, _chain(tmp_path, [link])).anchors_exact()["w"] == Fraction(1.01)


def test_a_self_loop_link_is_refused(tmp_path: Path) -> None:
    # An A/A session cited as a promotion REV_B -> REV_B would move A_w by noise while the reference stays (P4).
    first = promotion_link(tmp_path, REF, REV_B, CARRY_RESET)
    loop = promotion_link(tmp_path, REV_B, REV_B, {"w": ((100.0, 98.6),), "v": ((100.0, 95.0),)})
    _refused_both_ways(tmp_path, [first, loop], AnchorLinkSelfLoopError)


def test_a_link_citing_one_record_twice_is_refused(tmp_path: Path) -> None:
    link = promotion_link(tmp_path, REF, REV_B, RED_W, [("w", ((100.0, 101.0),))])
    doubled = dataclasses.replace(link, sessions=(*link.sessions, link.sessions[1]))
    _refused_both_ways(tmp_path, [doubled], DuplicateAnchorRecordError)


def test_a_link_citing_a_byte_copy_of_one_record_is_refused(tmp_path: Path) -> None:
    link = promotion_link(tmp_path, REF, REV_B, RED_W, [("w", ((100.0, 101.0),))])
    original = tmp_path / link.sessions[1].path
    copy = original.parent.parent / "ext-copy" / original.name
    copy.parent.mkdir()
    copy.write_bytes(original.read_bytes())
    copied = dataclasses.replace(link, sessions=(*link.sessions, cite(tmp_path, copy)))
    _refused_both_ways(tmp_path, [copied], DuplicateAnchorRecordError)


def test_write_registry_refuses_exact_anchors_that_are_not_the_derived_ones(tmp_path: Path) -> None:
    # Off the derived 1.02 by far less than half an ulp: the reporting floats agree, the decision values do not.
    link = promotion_link(tmp_path, REF, REV_B, CARRY_RESET)
    nudged = Anchor(Fraction(1.02) + Fraction(1, 10**30))
    assert float(nudged) == 1.02
    registry = dataclasses.replace(_raw_chain([link]), exact_anchors={"v": Anchor(ONE), "w": nudged})
    path = tmp_path / "perf_ratio_gate.json"
    with pytest.raises(AnchorValueMismatchError):
        write_registry(registry, path, root=tmp_path)
    assert not path.exists()


# Each chain field the self-digest covers, edited in place without re-digesting.
UNDIGESTED_CHAIN_EDITS: Mapping[str, Callable[[dict[str, Any]], None]] = {
    "link-margin": lambda document: document["anchor_links"][0]["calibration"]["margins"].update(w=0.5),
    "link-to-ref": lambda document: document["anchor_links"][0].update(to_ref=REV_D),
    "anchor-ref": lambda document: document.update(anchor_ref=REV_D),
    "anchor-value": lambda document: document["anchors"].update(w=1.5),
}


@pytest.mark.parametrize("edit", sorted(UNDIGESTED_CHAIN_EDITS))
def test_the_self_digest_covers_the_anchor_chain(tmp_path: Path, edit: str) -> None:
    path = _chained_file(tmp_path)
    _rewrite(path, UNDIGESTED_CHAIN_EDITS[edit], redigest=False)
    with pytest.raises(RegistryDigestMismatchError):
        load_registry(path, root=tmp_path)


def _drop_link_key(document: dict[str, Any]) -> None:
    del document["anchor_links"][0]["calibration"]


MALFORMED_CHAIN_EDITS: Mapping[str, Callable[[dict[str, Any]], None]] = {
    "anchor-a-decimal-string": lambda document: document["anchors"].update(w="1.02"),
    "anchor-ref-not-a-string": lambda document: document.update(anchor_ref=1),
    "anchor-links-not-a-list": lambda document: document.update(anchor_links={}),
    "link-missing-key": _drop_link_key,
    "link-extra-key": lambda document: document["anchor_links"][0].update(extra=1),
    "link-calibration-not-an-object": lambda document: document["anchor_links"][0].update(calibration="x"),
}


@pytest.mark.parametrize("edit", sorted(MALFORMED_CHAIN_EDITS))
def test_a_malformed_anchor_chain_is_refused(tmp_path: Path, edit: str) -> None:
    path = _chained_file(tmp_path)
    _rewrite(path, MALFORMED_CHAIN_EDITS[edit], redigest=True)
    with pytest.raises(MalformedRegistryError):
        load_registry(path, root=tmp_path)


@pytest.mark.parametrize("fault", ["missing", "outside"])
def test_a_link_calibration_refusal_names_its_link(tmp_path: Path, fault: str) -> None:
    link = promotion_link(tmp_path, REF, REV_B, CARRY_RESET)
    record = link.calibration.sessions[0]
    if fault == "missing":
        (tmp_path / record.path).unlink()
        error: type[Exception] = CalibrationRecordMissingError
    else:
        record = dataclasses.replace(record, path="../escaped/perf_ratio_session.json")
        error = CalibrationRecordOutsideRootError
    calibration = dataclasses.replace(link.calibration, sessions=(record, *link.calibration.sessions[1:]))
    with pytest.raises(error, match=rf"anchor link {REF} -> {REV_B}"):
        _chain(tmp_path, [dataclasses.replace(link, calibration=calibration)])
