from __future__ import annotations

import json

import pytest

from benchmarks.cli import AGGREGATE_CORPUS_NAMES
from benchmarks.cli import CORPUS_NAMES
from benchmarks.cli import EXIT_CORPUS_ERROR
from benchmarks.cli import EXIT_OK
from benchmarks.cli import build_corpus
from benchmarks.cli import build_parser
from benchmarks.cli import main
from benchmarks.errors import MissingExternalDirectoryError
from benchmarks.errors import UnknownCorpusError
from benchmarks.figures import DEFAULT_FIGURE_FORMATS
from benchmarks.figures import DEFAULT_FIGURES_OUT
from benchmarks.palette import Theme

TEST_TOOL_DIAMETER = 1.0
TEST_CAP_DEG = 120.0

# Corpora whose instances are authored in this repository, so they can be built
# without an operator supplying anything.
AUTHORED_CORPUS_NAMES = tuple(n for n in CORPUS_NAMES if n != "external")

# One cap and one spacing: the cheapest arguments that still drive both
# generators, the audit, the exceedance replay, and a successful baseline
# selection end to end. The published sweep is minutes; this is seconds. The
# spacing is one whose measured engagement is below the cap, so the selection
# path is exercised rather than the no-compliant-spacing path.
SMOKE_FIGURE6_CAP = "160"
SMOKE_FIGURE6_SPACINGS = ("0.3",)


def test_every_authored_corpus_builds_at_least_one_instance() -> None:
    for name in AUTHORED_CORPUS_NAMES:
        build = build_corpus(name, tool_diameter=TEST_TOOL_DIAMETER, tea_cap_deg=TEST_CAP_DEG)
        assert len(build.specs) > 0, name
        assert build.rejected == (), name


def test_the_aggregate_is_exactly_its_parts() -> None:
    """`all` must not quietly drop or duplicate a family as families are added."""
    expected = [spec.name for sub in AGGREGATE_CORPUS_NAMES for spec in build_corpus(sub, TEST_TOOL_DIAMETER, TEST_CAP_DEG).specs]
    assert [s.name for s in build_corpus("all", TEST_TOOL_DIAMETER, TEST_CAP_DEG).specs] == expected


def test_instance_names_are_unique_across_the_aggregate() -> None:
    """The name is the row key in every report; a collision silently merges rows."""
    names = [s.name for s in build_corpus("all", TEST_TOOL_DIAMETER, TEST_CAP_DEG).specs]
    assert len(set(names)) == len(names)


def test_an_unknown_corpus_name_is_rejected_by_name() -> None:
    with pytest.raises(UnknownCorpusError):
        build_corpus("does-not-exist", TEST_TOOL_DIAMETER, TEST_CAP_DEG)


def test_the_external_corpus_without_a_directory_is_rejected_by_name() -> None:
    with pytest.raises(MissingExternalDirectoryError):
        build_corpus("external", TEST_TOOL_DIAMETER, TEST_CAP_DEG)


def test_an_unknown_corpus_name_on_the_command_line_exits() -> None:
    with pytest.raises(SystemExit):
        main(["corpus", "--name", "does-not-exist"])


def test_no_command_exits_rather_than_doing_something_arbitrary() -> None:
    with pytest.raises(SystemExit):
        main([])


def test_the_external_corpus_without_a_directory_exits_non_zero(capsys) -> None:
    """A usage error is reported and exits; it must not surface as a traceback."""
    assert main(["corpus", "--name", "external"]) == EXIT_CORPUS_ERROR
    assert "MissingExternalDirectoryError" in capsys.readouterr().err


def test_a_smoke_corpus_run_writes_both_report_artifacts(tmp_path) -> None:
    assert main(["corpus", "--name", "smoke", "--out", str(tmp_path), "--no-digits"]) == EXIT_OK
    report = tmp_path / "benchmark_report.md"
    assert report.exists()
    assert (tmp_path / "benchmark_report.json").exists()
    assert "| instance |" in report.read_text()


def test_an_external_corpus_run_reports_its_rejections(tmp_path, capsys) -> None:
    """A thinner corpus must never read as a corpus that had nothing to reject."""
    profiles = tmp_path / "profiles"
    profiles.mkdir()
    (profiles / "ok.json").write_text(json.dumps({"name": "ok", "points": [[0, 0], [10, 0], [10, 10], [0, 10]]}))
    (profiles / "speck.json").write_text(json.dumps({"name": "speck", "points": [[0, 0], [0.2, 0], [0.2, 0.2], [0, 0.2]]}))
    out = tmp_path / "out"
    assert main(["corpus", "--name", "external", "--external-dir", str(profiles), "--out", str(out), "--no-digits"]) == EXIT_OK
    assert "rejected speck" in capsys.readouterr().err


def test_a_figure6_run_writes_both_figure_artifacts(tmp_path) -> None:
    exit_code = main(["figure6", "--out", str(tmp_path), "--caps", SMOKE_FIGURE6_CAP, "--spacings", *SMOKE_FIGURE6_SPACINGS])
    assert exit_code == EXIT_OK
    markdown = (tmp_path / "figure6.md").read_text()
    assert "| cap (deg) |" in markdown
    payload = json.loads((tmp_path / "figure6.json").read_text())
    assert [p["cap_deg"] for p in payload["points"]] == [float(SMOKE_FIGURE6_CAP)]
    assert len(payload["spacing_trials"]) == len(SMOKE_FIGURE6_SPACINGS)
    # The baseline was actually selected, so the whole comparison path ran rather
    # than short-circuiting on an empty trial set.
    assert payload["points"][0]["mathsm_length"] is not None
    assert payload["points"][0]["length_ratio"] is not None


def test_every_command_is_wired_to_a_handler() -> None:
    """A subcommand that parses but dispatches nowhere fails only in production."""
    for command in ("corpus", "figure6", "figures"):
        assert getattr(build_parser().parse_args([command]), "handler", None) is not None, command


def test_the_figures_command_defaults_to_the_published_location_and_format() -> None:
    """The docs link the figure by path, so its default output is part of the interface."""
    args = build_parser().parse_args(["figures"])
    assert args.out == DEFAULT_FIGURES_OUT
    assert tuple(args.formats) == DEFAULT_FIGURE_FORMATS
    assert args.theme == Theme.LIGHT.value
