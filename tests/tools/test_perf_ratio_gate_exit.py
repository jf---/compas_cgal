"""The ratio gate CLI's exit-code classes: 0 pass, 1 regression, 2 usage, 3 refused, 4 infra.

The codes are the contract a CI job branches on, so the tests assert the literal numbers,
never the module constants: a constant edited to another value must turn a test red. Every
session is synthetic and every record lives under ``tmp_path`` (``root=tmp_path``), as in
`tests.tools.test_perf_ratio_gate`, whose fixtures these tests reuse.
"""

from __future__ import annotations

import inspect
import json
import subprocess
from collections.abc import Callable
from pathlib import Path
from types import ModuleType

import pytest
from tools import perf_ratio_gate
from tools import perf_ratio_gate_chain
from tools import perf_ratio_gate_citations
from tools import perf_ratio_gate_decide
from tools import perf_ratio_gate_errors
from tools import perf_ratio_gate_merge
from tools import perf_ratio_gate_registry
from tools import perf_ratio_gate_verdict
from tools.perf_ratio import HarnessSourceMissingError
from tools.perf_ratio_arms import HeadCheckoutFailedError
from tools.perf_ratio_arms import HeadRevUnknownError
from tools.perf_ratio_gate import main
from tools.perf_ratio_gate_errors import GateRefusal
from tools.perf_ratio_gate_merge import ExtensionMismatchError
from tools.perf_ratio_gate_verdict import ExtensionRequest
from tools.perf_ratio_gate_verdict import RatioRegressionError
from tools.perf_ratio_record import InconsistentRatioRecordError
from tools.perf_ratio_record import MalformedRatioRecordError
from tools.perf_ratio_record import NonPositiveCpuError
from tools.perf_ratio_record import PerfRatioSession
from tools.perf_ratio_record import load_session
from tools.pipeline_timing_1_1_base import MeasuredTreeError
from tools.run_provenance import MachineIdentityError

from tests.tools.perf_ratio_sessions import CPU_MODEL
from tests.tools.perf_ratio_sessions import HARNESS
from tests.tools.perf_ratio_sessions import OTHER_REV
from tests.tools.perf_ratio_sessions import REF
from tests.tools.perf_ratio_sessions import write_record
from tests.tools.test_perf_ratio_gate import EXTENDED_TO
from tests.tools.test_perf_ratio_gate import GT
from tests.tools.test_perf_ratio_gate import HEAD
from tests.tools.test_perf_ratio_gate import HEAD2
from tests.tools.test_perf_ratio_gate import IMPROVED
from tests.tools.test_perf_ratio_gate import NEW_HARNESS
from tests.tools.test_perf_ratio_gate import OMEGA
from tests.tools.test_perf_ratio_gate import PASSING
from tests.tools.test_perf_ratio_gate import RED_FIRST
from tests.tools.test_perf_ratio_gate import REGRESSED
from tests.tools.test_perf_ratio_gate import FakeCheckout
from tests.tools.test_perf_ratio_gate import FakeMeasure
from tests.tools.test_perf_ratio_gate import aa_records
from tests.tools.test_perf_ratio_gate import calibrate_at
from tests.tools.test_perf_ratio_gate import calibrated
from tests.tools.test_perf_ratio_gate import gate_record
from tests.tools.test_perf_ratio_gate import gated
from tests.tools.test_perf_ratio_gate import no_checkout
from tests.tools.test_perf_ratio_gate import no_measure
from tests.tools.test_perf_ratio_gate import promoted_once
from tests.tools.test_perf_ratio_gate import promotion_target

# The worktree root, whose tools/ holds the gate modules.
ROOT = Path(__file__).resolve().parents[2]

# Every module of the gate. A new gate module joins this tuple, so its exceptions are classified too;
# `test_gate_modules_are_every_perf_ratio_gate_module` fails until it does.
GATE_MODULES: tuple[ModuleType, ...] = (
    perf_ratio_gate,
    perf_ratio_gate_registry,
    perf_ratio_gate_verdict,
    perf_ratio_gate_citations,
    perf_ratio_gate_chain,
    perf_ratio_gate_decide,
    perf_ratio_gate_merge,
    perf_ratio_gate_errors,
)

# Gate exceptions that are not refusals, each with its class. Nothing else may stay outside `GateRefusal`.
NOT_REFUSALS: frozenset[type[BaseException]] = frozenset(
    {
        RatioRegressionError,  # regression, exit 1: the gate decided, and HEAD is slower.
    }
)

# The exception classes the gate modules define (48 at 438cfe33, plus `GateRefusal`, plus the CI-N3 fix wave's
# `CalibrationRecordUnloadableError` and `AnchorRecordUnloadableError`): a floor, so an enumeration that silently
# finds nothing (a renamed module attribute, a broken filter) fails instead of passing vacuously.
MIN_GATE_EXCEPTIONS = 51

# The stderr prefix of every classified failure.
TAG = "ratio gate"


def gate_exceptions() -> list[type[BaseException]]:
    """Every exception class defined (not merely re-exported) in a gate module."""
    found: dict[str, type[BaseException]] = {}
    for module in GATE_MODULES:
        for _name, value in inspect.getmembers(module, inspect.isclass):
            if issubclass(value, BaseException) and value.__module__ == module.__name__:
                found[f"{value.__module__}.{value.__qualname__}"] = value
    return [found[key] for key in sorted(found)]


def run_main(argv: list[str], tmp_path: Path, **seams: object) -> int:
    options: dict[str, object] = {"root": tmp_path, "measure": no_measure, "checkout": no_checkout}
    options.update(seams)
    return main(argv, **options)  # type: ignore[arg-type]


def stderr_lines(capsys: pytest.CaptureFixture[str]) -> list[str]:
    return [line for line in capsys.readouterr().err.splitlines() if line.strip()]


class RaisingMeasure:
    """A `Measure` that raises what *make* builds, on its first request."""

    def __init__(self, make: Callable[[], BaseException]) -> None:
        self.make = make

    def __call__(self, request: ExtensionRequest) -> PerfRatioSession:
        raise self.make()


def red_session(tmp_path: Path) -> Path:
    """A calibrated registry's regressed session record: gating it requests an extension."""
    return write_record(gated({GT: PASSING, OMEGA: REGRESSED}), tmp_path / "gate" / "perf_ratio_session.json")


# --- completeness ------------------------------------------------------------------------------------


def test_the_enumeration_sees_every_gate_module_and_the_merge_refusal() -> None:
    found = gate_exceptions()
    assert len(found) >= MIN_GATE_EXCEPTIONS
    assert ExtensionMismatchError in found
    assert GateRefusal in found
    assert {cls.__module__ for cls in found} == {module.__name__ for module in GATE_MODULES if module is not perf_ratio_gate_decide}


def test_gate_modules_are_every_perf_ratio_gate_module() -> None:
    # A new ``tools/perf_ratio_gate*.py`` module must join GATE_MODULES, or its exceptions escape the completeness
    # check. A gate module without that prefix (analysis §7's ``tools/perf_ratio_ci.py``) is not globbed: add it by
    # hand, or name it ``perf_ratio_gate_*``.
    globbed = {f"tools.{path.stem}" for path in (ROOT / "tools").glob("perf_ratio_gate*.py")}
    assert {module.__name__ for module in GATE_MODULES} == globbed


@pytest.mark.parametrize("cls", gate_exceptions(), ids=lambda cls: cls.__qualname__)
def test_every_gate_exception_is_a_refusal_or_allow_listed(cls: type[BaseException]) -> None:
    if cls in NOT_REFUSALS:
        assert not issubclass(cls, GateRefusal), f"{cls.__qualname__} is allow-listed as a non-refusal but is a GateRefusal"
    else:
        assert issubclass(cls, GateRefusal), f"{cls.__module__}.{cls.__qualname__} is neither a GateRefusal nor allow-listed"


def test_gate_refusal_is_an_exception_but_no_harness_error() -> None:
    assert issubclass(GateRefusal, Exception)
    for harness in (MalformedRatioRecordError, OSError, subprocess.CalledProcessError):
        assert not issubclass(GateRefusal, harness)
        assert not issubclass(harness, GateRefusal)


# --- one CLI test per class ----------------------------------------------------------------------------


def test_pass_exits_0(tmp_path: Path) -> None:
    registry = calibrated(tmp_path)
    record = write_record(gated({GT: PASSING, OMEGA: PASSING}), tmp_path / "gate" / "perf_ratio_session.json")
    assert run_main(["check", str(record), "--registry", str(registry)], tmp_path) == 0


def test_a_regression_exits_1_and_tags_stderr(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    registry = calibrated(tmp_path)
    record = red_session(tmp_path)
    capsys.readouterr()
    assert run_main(["check", str(record), "--registry", str(registry)], tmp_path, measure=FakeMeasure(REGRESSED)) == 1
    (line,) = stderr_lines(capsys)
    assert line.startswith(f"{TAG}: regression: RatioRegressionError: ratio regression: {OMEGA}: ")


def test_an_absent_registry_exits_3(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    record = write_record(gated({GT: PASSING, OMEGA: PASSING}), tmp_path / "gate" / "perf_ratio_session.json")
    assert run_main(["check", str(record), "--registry", str(tmp_path / "absent.json")], tmp_path) == 3
    (line,) = stderr_lines(capsys)
    assert line.startswith(f"{TAG}: refused: MissingRegistryError: no ratio gate registry at ")


def test_a_schema_one_registry_exits_3(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    registry = calibrated(tmp_path)
    document = json.loads(registry.read_text(encoding="utf-8"))
    document["schema"] = "perf-ratio-gate/1"
    registry.write_text(json.dumps(document), encoding="utf-8")
    record = write_record(gated({GT: PASSING, OMEGA: PASSING}), tmp_path / "gate" / "perf_ratio_session.json")
    capsys.readouterr()
    assert run_main(["check", str(record), "--registry", str(registry)], tmp_path) == 3
    (line,) = stderr_lines(capsys)
    assert line.startswith(f"{TAG}: refused: SupersededRegistrySchemaError: ")


def test_an_extension_that_fails_to_run_exits_4(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    registry = calibrated(tmp_path)
    record = red_session(tmp_path)
    failing = RaisingMeasure(lambda: subprocess.CalledProcessError(1, ("perf_ratio",)))
    capsys.readouterr()
    assert run_main(["check", str(record), "--registry", str(registry)], tmp_path, measure=failing) == 4
    (line,) = stderr_lines(capsys)
    assert line.startswith(f"{TAG}: infra: CalledProcessError: ")


def test_a_session_record_that_fails_to_load_exits_4(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    registry = calibrated(tmp_path)
    record = tmp_path / "gate" / "perf_ratio_session.json"
    record.parent.mkdir(parents=True)
    record.write_text("not json", encoding="utf-8")
    capsys.readouterr()
    assert run_main(["check", str(record), "--registry", str(registry)], tmp_path) == 4
    (line,) = stderr_lines(capsys)
    assert line.startswith(f"{TAG}: infra: MalformedRatioRecordError: ")


# The harness errors `main` classifies as infra, named here independently of `perf_ratio_gate.INFRA_ERRORS`.
HARNESS_INFRA: dict[str, Callable[[], BaseException]] = {
    "CalledProcessError": lambda: subprocess.CalledProcessError(1, ("perf_ratio",)),
    "OSError": lambda: OSError("disk"),
    "FileNotFoundError": lambda: FileNotFoundError("ext-1.json"),
    "MalformedRatioRecordError": lambda: MalformedRatioRecordError("bad"),
    "InconsistentRatioRecordError": lambda: InconsistentRatioRecordError("bad"),
    "NonPositiveCpuError": lambda: NonPositiveCpuError("bad"),
    "HeadRevUnknownError": lambda: HeadRevUnknownError("bad"),
    "HeadCheckoutFailedError": lambda: HeadCheckoutFailedError("bad"),
    "MachineIdentityError": lambda: MachineIdentityError("bad"),
    "MeasuredTreeError": lambda: MeasuredTreeError("bad"),
    "HarnessSourceMissingError": lambda: HarnessSourceMissingError("bad"),
}


@pytest.mark.parametrize("name", sorted(HARNESS_INFRA))
def test_each_harness_error_is_classified_infra_without_a_traceback(name: str, tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    registry = calibrated(tmp_path)
    record = red_session(tmp_path)
    capsys.readouterr()
    assert run_main(["check", str(record), "--registry", str(registry)], tmp_path, measure=RaisingMeasure(HARNESS_INFRA[name])) == 4
    err = capsys.readouterr().err
    assert "Traceback" not in err
    (line,) = [line for line in err.splitlines() if line.strip()]
    assert line.startswith(f"{TAG}: infra: {name}: ")


def test_an_unknown_exception_exits_4_with_its_traceback(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    registry = calibrated(tmp_path)
    record = red_session(tmp_path)
    capsys.readouterr()
    assert run_main(["check", str(record), "--registry", str(registry)], tmp_path, measure=RaisingMeasure(lambda: ZeroDivisionError("a bug"))) == 4
    err = capsys.readouterr().err
    assert "Traceback (most recent call last)" in err
    assert "ZeroDivisionError: a bug" in err
    assert [line for line in err.splitlines() if line.strip()][-1] == f"{TAG}: infra (unclassified): ZeroDivisionError: a bug"


def test_a_refusal_raised_inside_an_extension_is_refused_not_infra(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    # The measured extension ran on another machine: the merge refuses it (ExtensionMismatchError, a GateRefusal).
    registry = calibrated(tmp_path)
    record = red_session(tmp_path)
    capsys.readouterr()
    assert run_main(["check", str(record), "--registry", str(registry)], tmp_path, measure=FakeMeasure(REGRESSED, cpu_model="Intel Xeon W-3245")) == 3
    (line,) = stderr_lines(capsys)
    assert line.startswith(f"{TAG}: refused: ExtensionMismatchError: ")


def test_usage_exits_2(tmp_path: Path) -> None:
    with pytest.raises(SystemExit) as raised:
        run_main(["--bogus"], tmp_path)
    assert raised.value.code == 2


@pytest.mark.parametrize(
    "head_cpus",
    [
        {GT: IMPROVED, OMEGA: REGRESSED},
        {GT: PASSING, OMEGA: PASSING},
    ],
    ids=["regression", "no-improvement"],
)
def test_promote_refusals_exit_3(head_cpus: dict[str, float], tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    registry = calibrated(tmp_path)
    before = registry.read_bytes()
    record = write_record(gated(head_cpus), tmp_path / "gate" / "perf_ratio_session.json")
    capsys.readouterr()
    assert run_main(["promote", str(record), "--registry", str(registry)], tmp_path, measure=FakeMeasure(REGRESSED), checkout=FakeCheckout()) == 3
    (line,) = stderr_lines(capsys)
    assert line.startswith(f"{TAG}: refused: PromotionRefusedError: ")
    assert registry.read_bytes() == before


def test_calibrate_refusals_exit_3(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    registry = calibrated(tmp_path)
    argv = ["calibrate", "--sessions", "1", "--registry", str(registry)]
    assert run_main(argv, tmp_path, cpu_model=lambda: CPU_MODEL, tree_harness=lambda: HARNESS) == 3
    (line,) = stderr_lines(capsys)
    assert line.startswith(f"{TAG}: refused: TooFewCalibrationSessionsError: ")


# --- cited evidence that no longer loads: refused, not infra (CI-N3 review M1, m3, m4) --------------------

# The harness record loader's failures (its module is a HARNESS_SOURCE): what a parser change raises on cited
# records whose bytes, and so whose digests, did not change. UnicodeDecodeError: bytes that are not UTF-8.
RECORD_LOAD_FAILURES: dict[str, Callable[[Path], Exception]] = {
    "MalformedRatioRecordError": lambda path: MalformedRatioRecordError(f"{path}: perf-ratio/2 expected"),
    "InconsistentRatioRecordError": lambda path: InconsistentRatioRecordError(f"{path}: ratio recomputes otherwise"),
    "NonPositiveCpuError": lambda path: NonPositiveCpuError(f"{path}: a run's CPU is 0"),
    "UnicodeDecodeError": lambda path: UnicodeDecodeError("utf-8", b"\xff", 0, 1, "invalid start byte"),
}


def tightened_loader(fail: Callable[[Path], Exception], under: Path) -> Callable[[Path], PerfRatioSession]:
    """The harness loader after a parser change: a record under *under* fails with *fail*, every other one loads."""

    def load(path: Path) -> PerfRatioSession:
        if path.resolve().is_relative_to(under.resolve()):
            raise fail(path)
        return load_session(path)

    return load


@pytest.mark.parametrize("failure", sorted(RECORD_LOAD_FAILURES))
def test_a_cited_calibration_record_that_no_longer_loads_exits_3(failure: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]) -> None:
    registry = calibrated(tmp_path)
    record = write_record(gated({GT: PASSING, OMEGA: PASSING}), gate_record(tmp_path))
    monkeypatch.setattr(perf_ratio_gate_citations, "load_session", tightened_loader(RECORD_LOAD_FAILURES[failure], tmp_path / "records"))
    capsys.readouterr()
    assert run_main(["check", str(record), "--registry", str(registry)], tmp_path) == 3
    (line,) = stderr_lines(capsys)
    assert line.startswith(f"{TAG}: refused: CalibrationRecordUnloadableError: ")
    assert f"does not load: {failure}: " in line


def test_a_cited_anchor_record_that_no_longer_loads_exits_3(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]) -> None:
    # One link REF -> HEAD, then calibrated at HEAD: only the link's promotion records fail to load.
    registry = promoted_once(tmp_path)
    calibrate_at(tmp_path, registry, HEAD)
    record = write_record(gated({GT: PASSING, OMEGA: PASSING}, ref_commit=HEAD, head_rev=HEAD2), gate_record(tmp_path))
    fail = RECORD_LOAD_FAILURES["MalformedRatioRecordError"]
    monkeypatch.setattr(perf_ratio_gate_citations, "load_session", tightened_loader(fail, promotion_target(tmp_path, REF, HEAD)))
    capsys.readouterr()
    assert run_main(["check", str(record), "--registry", str(registry)], tmp_path) == 3
    (line,) = stderr_lines(capsys)
    assert line.startswith(f"{TAG}: refused: AnchorRecordUnloadableError: ")


def test_rebase_harness_on_a_registry_whose_cited_records_no_longer_load_exits_3(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]) -> None:
    # The documented recovery after a harness edit loads the registry first: its unloadable citations refuse it.
    registry = calibrated(tmp_path)
    rebased = aa_records(tmp_path / "rebased", harness=NEW_HARNESS)
    fail = RECORD_LOAD_FAILURES["MalformedRatioRecordError"]
    monkeypatch.setattr(perf_ratio_gate_citations, "load_session", tightened_loader(fail, tmp_path / "records"))
    capsys.readouterr()
    argv = ["calibrate", "--records", *map(str, rebased), "--registry", str(registry), "--rebase-harness"]
    assert run_main(argv, tmp_path) == 3
    (line,) = stderr_lines(capsys)
    assert line.startswith(f"{TAG}: refused: CalibrationRecordUnloadableError: ")


def test_a_persisted_extension_that_no_longer_loads_exits_3_on_a_promote_rerun(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    registry = calibrated(tmp_path)
    record = write_record(gated({GT: IMPROVED, OMEGA: RED_FIRST}), gate_record(tmp_path))
    argv = ["promote", str(record), "--registry", str(registry)]
    # The first run persists ext-1, then its checkout lands elsewhere (a refusal); ext-1 is then made unloadable.
    assert run_main(argv, tmp_path, measure=FakeMeasure(EXTENDED_TO), checkout=FakeCheckout(answer=OTHER_REV)) == 3
    (promotion_target(tmp_path, REF, HEAD) / "ext-1" / "perf_ratio_session.json").write_text("not json", encoding="utf-8")
    capsys.readouterr()
    assert run_main(argv, tmp_path) == 3
    (line,) = stderr_lines(capsys)
    assert line.startswith(f"{TAG}: refused: PromotionRecordsExistError: ")


def test_a_registry_that_is_not_utf_8_exits_3(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    registry = calibrated(tmp_path)
    registry.write_bytes(b"\xff\xfe{}")
    record = write_record(gated({GT: PASSING, OMEGA: PASSING}), gate_record(tmp_path))
    capsys.readouterr()
    assert run_main(["check", str(record), "--registry", str(registry)], tmp_path) == 3
    (line,) = stderr_lines(capsys)
    assert line.startswith(f"{TAG}: refused: MalformedRegistryError: ")


def test_a_promoted_record_that_is_not_utf_8_is_classified_infra(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    # A session, so infra; but a named one, not "infra (unclassified)" with a traceback.
    registry = calibrated(tmp_path)
    record = gate_record(tmp_path)
    record.parent.mkdir(parents=True, exist_ok=True)
    record.write_bytes(b"\xff\xfe{}")
    capsys.readouterr()
    assert run_main(["promote", str(record), "--registry", str(registry)], tmp_path) == 4
    err = capsys.readouterr().err
    assert "Traceback" not in err
    (line,) = [line for line in err.splitlines() if line.strip()]
    assert line.startswith(f"{TAG}: infra: MalformedRatioRecordError: ")


def test_a_checked_record_that_is_not_utf_8_is_classified_infra(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    # A session, so infra; but a named one, not "infra (unclassified)" with a traceback.
    registry = calibrated(tmp_path)
    record = gate_record(tmp_path)
    record.parent.mkdir(parents=True, exist_ok=True)
    record.write_bytes(b"\xff\xfe{}")
    capsys.readouterr()
    assert run_main(["check", str(record), "--registry", str(registry)], tmp_path) == 4
    err = capsys.readouterr().err
    assert "Traceback" not in err
    (line,) = [line for line in err.splitlines() if line.strip()]
    assert line.startswith(f"{TAG}: infra: MalformedRatioRecordError: ")
    assert f"{record} is not UTF-8: " in line
