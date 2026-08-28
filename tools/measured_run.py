"""Produce one authenticated all-corpus reporting bundle from clean source."""

from __future__ import annotations

import argparse
import datetime
import pathlib
import subprocess
import sys
from typing import Optional
from typing import Sequence

from benchmarks.cli import AGGREGATE_CORPUS_NAMES
from benchmarks.cli import DEFAULT_CAP_DEG
from benchmarks.cli import DEFAULT_TOOL_DIAMETER
from benchmarks.report import JSON_NAME
from benchmarks.report import MARKDOWN_NAME
from tools import corpus_result
from tools import measurement_artifact
from tools.corpus_result import InvalidMeasuredResultError  # noqa: F401
from tools.corpus_result import NoMeasuredResultError  # noqa: F401
from tools.corpus_result import latest_result  # noqa: F401
from tools.corpus_result import validate_result  # noqa: F401
from tools.measurement_artifact import ENVELOPE_KEYS as STAMP_KEYS  # noqa: F401
from tools.measurement_artifact import DirtyMeasurementTreeError as DirtyMeasuredRunError  # noqa: F401
from tools.measurement_artifact import MeasurementArtifactCollisionError as MeasuredResultCollisionError
from tools.measurement_artifact import MeasurementArtifactError as MeasuredRunError
from tools.measurement_artifact import MeasurementInputChangedError as MeasuredRunInputChangedError  # noqa: F401


class MeasuredRunChildError(MeasuredRunError):
    """The corpus child could not start or terminated by signal."""


def _repository_root() -> pathlib.Path:
    return measurement_artifact._resolve_repository(pathlib.Path(__file__).resolve().parents[1])


def _utc_now() -> datetime.datetime:
    return datetime.datetime.now(datetime.timezone.utc)


def _run_child(command: tuple[str, ...], cwd: pathlib.Path) -> int:
    try:
        completed = subprocess.run(command, cwd=cwd, check=False)
    except OSError as exc:
        raise MeasuredRunChildError(f"could not spawn corpus child in {cwd}") from exc
    if completed.returncode < 0:
        raise MeasuredRunChildError(f"corpus child terminated by signal {-completed.returncode}")
    return completed.returncode


def _input_payload() -> dict[str, object]:
    return {
        "name": "all",
        "tool_diameter": DEFAULT_TOOL_DIAMETER,
        "cap_deg": DEFAULT_CAP_DEG,
        "collect_digits": True,
        "aggregate_corpora": list(AGGREGATE_CORPUS_NAMES),
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Measure the fixed authored aggregate and publish its authenticated bundle."""
    parser = argparse.ArgumentParser(description="Run the all-corpus benchmark from clean committed source.")
    parser.parse_args(argv)
    repository = _repository_root()
    started = _utc_now()
    source = measurement_artifact.capture_clean_source(repository)
    logical_name = f"{started.date().isoformat()}-{str(source.commit)[:12]}"
    results = repository / "benchmarks" / "results"
    final = results / logical_name
    if final.exists() or final.is_symlink():
        raise MeasuredResultCollisionError(f"measurement result already exists: {final}")
    with measurement_artifact.publication_stage(results, logical_name) as stage:
        stage_argument = stage.relative_to(repository).as_posix()
        command = (
            sys.executable,
            "-m",
            "benchmarks.cli",
            "corpus",
            "--name",
            "all",
            "--out",
            stage_argument,
        )
        status = _run_child(command, repository)
        if status != 0:
            return status
        corpus_result._validate_corpus_payloads(stage)
        payloads = {
            MARKDOWN_NAME: (stage / MARKDOWN_NAME).read_bytes(),
            JSON_NAME: (stage / JSON_NAME).read_bytes(),
        }
        envelope = measurement_artifact.build_envelope(
            artifact_kind=corpus_result.ARTIFACT_KIND,
            source=source,
            started=started,
            finished=_utc_now(),
            argv=command,
            input_version=measurement_artifact.IdentityVersion(corpus_result.INPUT_VERSION),
            input_payload=_input_payload(),
            result_version=measurement_artifact.IdentityVersion(corpus_result.RESULT_VERSION),
            payloads=payloads,
        )
        measurement_artifact.write_envelope(stage, envelope)
        corpus_result._validate_result(stage, repository=repository, logical_name=logical_name)
        published = measurement_artifact.publish_stage(source=source, stage=stage, final=final)
    print(published.relative_to(repository).as_posix())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
