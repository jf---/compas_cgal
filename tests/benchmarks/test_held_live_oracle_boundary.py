from datetime import datetime
from datetime import timezone
from typing import cast

import pytest

import tests.benchmarks.held_figure5_path_live_oracle as oracle_module
from benchmarks.held_path_characterize import CharacterizationPhase


class MissingReport:
    def exists(self) -> bool:
        return False


def test_live_oracle_records_terminal_failure_and_bare_reraises(monkeypatch: pytest.MonkeyPatch) -> None:
    failure = RuntimeError("sentinel quality failure")
    ledgers: list[dict[str, object]] = []
    phases = cast(
        tuple[CharacterizationPhase, ...],
        ("generation", "guarded_audit", "survey", "quality_reduction"),
    )

    class FixedDateTime:
        @staticmethod
        def now(zone: timezone) -> datetime:
            return datetime(2026, 9, 5, 10, 0, tzinfo=zone)

    def writer(*, pixi_command: str, phase_observer: object, path: object) -> object:
        assert pixi_command == oracle_module.LIVE_COMMAND
        assert path is oracle_module.REPORT_PATH
        for phase in phases:
            cast(object, phase_observer)(phase)  # type: ignore[operator]
        raise failure

    monotonic = iter((10.0, 11.0, 12.0, 13.0, 788.885679, 790.83))
    monkeypatch.setattr(oracle_module, "REPORT_PATH", MissingReport())
    monkeypatch.setattr(oracle_module, "DEFAULT_REPORT_PATH", oracle_module.REPORT_PATH)
    monkeypatch.setattr(oracle_module, "datetime", FixedDateTime)
    monkeypatch.setattr(oracle_module.time, "monotonic", lambda: next(monotonic))
    monkeypatch.setattr(oracle_module, "write_held_figure5_report", writer)
    monkeypatch.setattr(oracle_module, "_write_live_run_ledger", lambda **values: ledgers.append(values))

    with pytest.raises(RuntimeError) as caught:
        oracle_module.test_held_figure5_path_live_oracle()

    assert caught.value is failure
    terminal = ledgers[-1]
    assert terminal["outcome"] == "failed"
    assert terminal["observed_phases"] == phases
    assert terminal["elapsed_seconds"] == 778.885679
    assert terminal["report_generated_at"] is None
    assert terminal["qualification_verdict"] == "not evaluated"
    assert terminal["qualification_failures"] is None
    assert terminal["command_result"] == ("pytest failed after 780.830000 seconds: RuntimeError: sentinel quality failure")
