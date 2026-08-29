from __future__ import annotations

import ast
import datetime
import importlib
import json
import pathlib
import subprocess
import sys
from typing import Any

import pytest


SOURCE = "a" * 40
HISTORY = "70049dd991e7d5c7b93d512785393a49a7f03564"
MC013_HISTORY = "1e7d48e3d6b115d1ab4cb61c43ca21d2bc9bb6fd"
CAPS = [20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0]
SPACINGS = [0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6]
PROJECT_ROOT = pathlib.Path(__file__).parents[2]


def _module() -> Any:
    return importlib.import_module("tools.measurement_claim_result")


def _raw() -> dict[str, object]:
    return {
        "pocket": {"name": "rect_20x12", "family": "analytic", "tool_diameter": 2.0, "params": {"width": 20.0, "height": 12.0}, "holes": []},
        "engagement_measured_at_cap_deg": 180.0,
        "points": [
            {
                "cap_deg": cap,
                "controlled_length": 200.0,
                "controlled_cut_motions": 10,
                "controlled_entry_cuts": 1,
                "controlled_max_tea_deg": 180.0,
                "controlled_max_tea_after_entry_deg": 120.0,
                "controlled_exceedances_after_entry": 0,
                "controlled_meets_cap": cap >= 120.0,
                "mathsm_spacing_tool_diameters": 0.1,
                "mathsm_length": 300.0,
                "mathsm_max_tea_after_entry_deg": 100.0,
                "length_ratio": 2.0 / 3.0,
            }
            for cap in CAPS
        ],
        "spacing_trials": [
            {
                "spacing_tool_diameters": spacing,
                "length": 200.0,
                "cut_motions": 10,
                "entry_cuts": 1,
                "max_tea_deg": 180.0,
                "max_tea_after_entry_deg": 131.14 if spacing == 0.025 else (100.0 if spacing == 0.1 else 120.0),
            }
            for spacing in SPACINGS
        ],
    }


def _markdown() -> bytes:
    comparison = "\n".join("| a | b | c | d | e | f | g | h | i |" for _ in CAPS)
    trials = "\n".join("| a | b | c | d | e |" for _ in SPACINGS)
    return (
        "# Figure 6 reproduction — path length against engagement cap (rect_20x12)\n"
        "## Per-cap comparison\n"
        "| cap (deg) | controlled length | controlled max TEA after entry (deg) | controlled meets cap | "
        "controlled exceedances after entry | spacing (tool diam.) | constant-spacing length | "
        "constant-spacing max TEA after entry (deg) | length ratio |\n"
        "| ---: | ---: | ---: | :--- | ---: | ---: | ---: | ---: | ---: |\n"
        f"{comparison}\n## Constant-spacing trials\n"
        "| spacing (tool diam.) | length | cut motions | max TEA after entry (deg) | raw max TEA (deg) |\n"
        "| ---: | ---: | ---: | ---: | ---: |\n"
        f"{trials}\n"
    ).encode("utf-8")


def _payload(module: Any) -> dict[str, object]:
    return module.compose_benchmark_payload(source_commit=SOURCE, figure6_payload=_raw(), figure6_markdown=_markdown())


def test_benchmark_payload_maps_exact_four_claims_and_json_values() -> None:
    module = _module()
    payload = _payload(module)
    validated = module.validate_benchmark_payload(payload, figure6_payload=_raw(), figure6_markdown=_markdown())
    assert [claim["claim_id"] for claim in validated["claims"]] == ["MC-011", "MC-012", "MC-013", "MC-014"]
    assert [claim["disposition"] for claim in validated["claims"]] == ["not-a-claim", "deleted", "corrected", "not-a-claim"]
    assert [claim["history_commit"] for claim in validated["claims"]] == [HISTORY, HISTORY, MC013_HISTORY, HISTORY]
    assert [claim["missing_inputs"] for claim in validated["claims"]] == [[], ["generator", "circle-selection", "entry-treatment", "operation-enumeration"], [], []]
    assert validated["claims"][2]["selected_values"]["fine_max_tea_after_entry"] == 131.14


@pytest.mark.parametrize(
    "damage",
    [
        "holes",
        "spacing-order",
        "selected-equality",
        "selected-inversion",
        "selected-absent",
        "selected-duplicate",
        "derived-selected-value",
        "claim-source",
    ],
)
def test_benchmark_validation_rejects_raw_or_derived_contract_damage(damage: str) -> None:
    module = _module()
    raw = _raw()
    payload = _payload(module)
    if damage == "holes":
        raw["pocket"]["holes"] = [[[0.0, 0.0, 0.0]]]
    elif damage == "spacing-order":
        raw["spacing_trials"].reverse()
    elif damage == "selected-equality":
        raw["spacing_trials"][0]["max_tea_after_entry_deg"] = 100.0
    elif damage == "selected-inversion":
        raw["spacing_trials"][0]["max_tea_after_entry_deg"] = 90.0
    elif damage == "selected-absent":
        raw["spacing_trials"][0]["spacing_tool_diameters"] = 0.03
    elif damage == "selected-duplicate":
        raw["spacing_trials"][1]["spacing_tool_diameters"] = 0.025
    elif damage == "derived-selected-value":
        payload["claims"][2]["selected_values"]["fine_max_tea_after_entry"] = 999.0
    else:
        payload["claims"][0]["source_commit"] = "b" * 40
    with pytest.raises(module.InvalidMeasurementClaimPayloadError):
        module.validate_benchmark_payload(payload, figure6_payload=raw, figure6_markdown=_markdown())


@pytest.mark.parametrize(
    "damage",
    [
        "invalid-utf8",
        "empty",
        "duplicate-h1",
        "duplicate-h2",
        "reordered",
        "comparison-heading",
        "comparison-header",
        "comparison-rule",
        "comparison-table-moved",
        "comparison-seven-rows",
        "comparison-nine-rows",
        "comparison-width",
        "comparison-missing-closing-pipe",
        "trial-heading",
        "trial-header",
        "trial-rule",
        "trial-eleven-rows",
        "trial-thirteen-rows",
        "trial-missing-closing-pipe",
        "width",
    ],
)
def test_markdown_structure_is_strict_but_numbers_are_not_authority(damage: str) -> None:
    module = _module()
    raw = _raw()
    markdown = _markdown()
    if damage == "invalid-utf8":
        markdown = b"\xff"
    elif damage == "empty":
        markdown = b"  \n"
    elif damage == "duplicate-h1":
        markdown += b"# duplicate\n"
    elif damage == "duplicate-h2":
        markdown += b"## Constant-spacing trials\n"
    elif damage == "reordered":
        markdown = (
            markdown.replace(b"## Per-cap comparison", b"## Later")
            .replace(b"## Constant-spacing trials", b"## Per-cap comparison")
            .replace(b"## Later", b"## Constant-spacing trials")
        )
    elif damage == "comparison-heading":
        markdown = markdown.replace(b"## Per-cap comparison", b"## Comparison")
    elif damage == "comparison-header":
        markdown = markdown.replace(b"| cap (deg) |", b"| cap |", 1)
    elif damage == "comparison-rule":
        markdown = markdown.replace(b"| ---: | ---: | ---: | :--- |", b"| --- | ---: | ---: | :--- |", 1)
    elif damage == "comparison-table-moved":
        h1, rest = markdown.split(b"## Per-cap comparison\n", 1)
        comparison, trials = rest.split(b"## Constant-spacing trials\n", 1)
        markdown = h1 + b"## Per-cap comparison\n" + trials + b"## Constant-spacing trials\n" + comparison
    elif damage == "comparison-seven-rows":
        markdown = markdown.replace(b"| a | b | c | d | e | f | g | h | i |\n", b"", 1)
    elif damage == "comparison-nine-rows":
        markdown = markdown.replace(
            b"## Constant-spacing trials",
            b"| a | b | c | d | e | f | g | h | i |\n## Constant-spacing trials",
        )
    elif damage == "comparison-width":
        markdown = markdown.replace(b"| a | b | c | d | e | f | g | h | i |", b"| a | b | c | d | e | f | g | h |", 1)
    elif damage == "comparison-missing-closing-pipe":
        markdown = markdown.replace(b"| a | b | c | d | e | f | g | h | i |", b"| a | b | c | d | e | f | g | h | i", 1)
    elif damage == "trial-heading":
        markdown = markdown.replace(b"## Constant-spacing trials", b"## Spacing trials")
    elif damage == "trial-header":
        markdown = markdown.replace(b"| spacing (tool diam.) | length |", b"| spacing | length |", 1)
    elif damage == "trial-rule":
        markdown = markdown.replace(b"| ---: | ---: | ---: | ---: | ---: |", b"| --- | ---: | ---: | ---: | ---: |", 1)
    elif damage == "trial-eleven-rows":
        markdown = markdown.replace(b"| a | b | c | d | e |\n", b"", 1)
    elif damage == "trial-thirteen-rows":
        markdown += b"| a | b | c | d | e |\n"
    elif damage == "trial-missing-closing-pipe":
        before, after = markdown.rsplit(b"| a | b | c | d | e |", 1)
        markdown = before + b"| a | b | c | d | e" + after
    elif damage == "width":
        markdown = markdown.replace(b"| a | b | c | d | e |", b"| a | b | c | d |", 1)
    with pytest.raises(module.InvalidMeasurementClaimPayloadError):
        module.validate_benchmark_payload(_payload(module), figure6_payload=raw, figure6_markdown=markdown)

    numeric_mutation = _markdown().replace(b"| a | b | c |", b"| 999 | b | c |", 1)
    assert module.validate_benchmark_payload(_payload(module), figure6_payload=raw, figure6_markdown=numeric_mutation)


@pytest.mark.parametrize("table", ["comparison", "trials"])
@pytest.mark.parametrize("gap", [b"\n", b"prose before rows\n\n"], ids=["blank", "prose-then-blank"])
def test_markdown_rejects_gap_after_rule_but_allows_prose_after_completed_rows(table: str, gap: bytes) -> None:
    module = _module()
    rules = {
        "comparison": b"| ---: | ---: | ---: | :--- | ---: | ---: | ---: | ---: | ---: |\n",
        "trials": b"| ---: | ---: | ---: | ---: | ---: |\n",
    }
    markdown = _markdown()
    if table == "comparison":
        broken = markdown.replace(rules[table], rules[table] + gap, 1)
    else:
        comparison, trials = markdown.split(b"## Constant-spacing trials\n", 1)
        broken = comparison + b"## Constant-spacing trials\n" + trials.replace(rules[table], rules[table] + gap, 1)
    with pytest.raises(module.InvalidMeasurementClaimPayloadError, match=rf"figure6\.md {table}: requires exactly"):
        module.validate_benchmark_payload(_payload(module), figure6_payload=_raw(), figure6_markdown=broken)

    valid = _markdown()
    if table == "comparison":
        valid = valid.replace(b"\n## Constant-spacing trials", b"\nprose after comparison rows\n\n## Constant-spacing trials", 1)
    else:
        valid += b"prose after trial rows\n"
    assert module.validate_benchmark_payload(_payload(module), figure6_payload=_raw(), figure6_markdown=valid)


def test_result_facade_contains_reexports_only() -> None:
    path = pathlib.Path(importlib.import_module("tools.measurement_claim_result").__file__)
    tree = ast.parse(path.read_text(encoding="utf-8"))
    assert not [node for node in tree.body if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))]


def test_benchmark_artifact_round_trips_common_and_family_authentication(tmp_path: pathlib.Path) -> None:
    module = _module()
    artifact = importlib.import_module("tools.measurement_artifact")
    repository = tmp_path / "repository"
    repository.mkdir()
    subprocess.run(["git", "-C", str(repository), "init", "-q"], check=True)
    (repository / "pixi.lock").write_bytes(b"benchmark lock\n")
    subprocess.run(["git", "-C", str(repository), "add", "pixi.lock"], check=True)
    subprocess.run(
        ["git", "-C", str(repository), "-c", "user.name=Jelle Feringa", "-c", "user.email=jelleferinga@gmail.com", "commit", "-qm", "lock"],
        check=True,
    )
    source = artifact.capture_clean_source(repository.resolve())
    raw_bytes = (json.dumps(_raw(), indent=2, allow_nan=False) + "\n").encode("utf-8")
    markdown = _markdown()
    payload = module.compose_benchmark_payload(source_commit=source.commit, figure6_payload=_raw(), figure6_markdown=markdown)
    payload_bytes = (json.dumps(payload, indent=2, allow_nan=False) + "\n").encode("utf-8")
    semantic = module.benchmark_semantic_input(source.commit)
    digest = artifact.input_identity_sha256(version=artifact.IdentityVersion(module.BENCHMARK_INPUT_VERSION), payload=semantic)
    started = datetime.datetime(2026, 8, 29, 12, 0, 0, 123456, tzinfo=datetime.timezone.utc)
    logical = f"2026-08-29-{str(source.commit)[:12]}-benchmark-{str(digest)[:12]}"
    result = repository / "benchmarks" / "measurement_claim_results" / logical
    result.mkdir(parents=True)
    payloads = {"figure6.md": markdown, "figure6.json": raw_bytes, "benchmark-claims.json": payload_bytes}
    for name, data in payloads.items():
        (result / name).write_bytes(data)
    argv = ("python", "-m", "benchmarks.cli", "figure6", "--out", str(result), *module.FIGURE6_SEMANTIC_COMMAND[3:])
    envelope = artifact.build_envelope(
        artifact_kind=artifact.ArtifactKind(module.BENCHMARK_ARTIFACT_KIND),
        source=source,
        started=started,
        finished=started + datetime.timedelta(seconds=1),
        argv=argv,
        input_version=artifact.IdentityVersion(module.BENCHMARK_INPUT_VERSION),
        input_payload=semantic,
        result_version=artifact.IdentityVersion(module.BENCHMARK_RESULT_VERSION),
        payloads=payloads,
    )
    artifact.write_envelope(result, envelope)
    validated = module.validate_benchmark_claim_artifact(result)
    assert validated[0] == payload
    assert validated[1].input_sha256 == digest

    (result / "figure6.md").write_bytes(markdown.replace(b"| a | b | c |", b"| 999 | b | c |", 1))
    with pytest.raises(artifact.InvalidMeasurementEnvelopeError, match="digest"):
        module.validate_benchmark_claim_artifact(result)


@pytest.mark.parametrize(
    ("holes", "message"),
    [
        ((), r"pocket\.holes: must be an exact JSON array"),
        ([()], r"pocket\.holes\[0\]: must be an exact JSON array"),
        ([[(0.0, 1.0, 2.0)]], r"pocket\.holes\[0\]\[0\]: must be an exact JSON array"),
        ([[[0.0, 1.0]]], r"pocket\.holes\[0\]\[0\]: must contain exactly three"),
        ([[[0.0, 1.0, 2.0, 3.0]]], r"pocket\.holes\[0\]\[0\]: must contain exactly three"),
        ([[[0.0, 1.0, float("nan")]]], r"pocket\.holes\[0\]\[0\]\[2\]: must be finite numeric"),
        ([[[0.0, 1.0, True]]], r"pocket\.holes\[0\]\[0\]\[2\]: must be finite numeric"),
        ([[[0.0, 1.0, "2.0"]]], r"pocket\.holes\[0\]\[0\]\[2\]: must be finite numeric"),
    ],
)
def test_runtime_hole_points_reject_wrong_wire_shape_or_value(holes: object, message: str) -> None:
    module = _module()
    raw = _raw()
    raw["pocket"]["holes"] = holes
    with pytest.raises(module.InvalidMeasurementClaimPayloadError, match=message):
        module.validate_benchmark_payload(_payload(module), figure6_payload=raw, figure6_markdown=_markdown())


def test_well_formed_nonempty_holes_reach_fixed_no_holes_contract() -> None:
    module = _module()
    raw = _raw()
    raw["pocket"]["holes"] = [[[0, 1.0, 2]]]
    with pytest.raises(module.InvalidMeasurementClaimPayloadError, match=r"pocket\.holes: fixed Figure-6 config requires exactly \[\]"):
        module.validate_benchmark_payload(_payload(module), figure6_payload=raw, figure6_markdown=_markdown())


def test_commit1_python_modules_parse_at_the_python39_grammar_floor() -> None:
    for relative in (
        "tools/measurement_claim_artifact_validation.py",
        "tools/measurement_claim_benchmark_identity.py",
        "tools/measurement_claim_benchmark_schema.py",
        "tools/measurement_claim_benchmark_semantic_input.py",
        "tools/measurement_claim_benchmark_validation.py",
        "tools/measurement_claim_figure6.py",
        "tools/measurement_claim_ledger.py",
        "tools/measurement_claim_probes.py",
        "tools/measurement_claim_result.py",
        "tools/measurement_claim_task6_ledger.py",
    ):
        path = PROJECT_ROOT / relative
        ast.parse(path.read_text(encoding="utf-8"), filename=str(path), feature_version=(3, 9))


def test_benchmark_artifact_rejects_extra_authenticated_payload(tmp_path: pathlib.Path) -> None:
    module = _module()
    artifact = importlib.import_module("tools.measurement_artifact")
    repository = tmp_path / "repository"
    repository.mkdir()
    subprocess.run(["git", "-C", str(repository), "init", "-q"], check=True)
    (repository / "pixi.lock").write_bytes(b"benchmark lock\n")
    subprocess.run(["git", "-C", str(repository), "add", "pixi.lock"], check=True)
    subprocess.run(
        ["git", "-C", str(repository), "-c", "user.name=Jelle Feringa", "-c", "user.email=jelleferinga@gmail.com", "commit", "-qm", "lock"],
        check=True,
    )
    source = artifact.capture_clean_source(repository.resolve())
    raw_bytes = (json.dumps(_raw(), indent=2, allow_nan=False) + "\n").encode("utf-8")
    markdown = _markdown()
    payload = module.compose_benchmark_payload(source_commit=source.commit, figure6_payload=_raw(), figure6_markdown=markdown)
    payload_bytes = (json.dumps(payload, indent=2, allow_nan=False) + "\n").encode("utf-8")
    semantic = module.benchmark_semantic_input(source.commit)
    digest = artifact.input_identity_sha256(version=artifact.IdentityVersion(module.BENCHMARK_INPUT_VERSION), payload=semantic)
    started = datetime.datetime(2026, 8, 29, 12, 0, 0, 123456, tzinfo=datetime.timezone.utc)
    logical = f"2026-08-29-{str(source.commit)[:12]}-benchmark-{str(digest)[:12]}"
    result = repository / "benchmarks" / "measurement_claim_results" / logical
    result.mkdir(parents=True)
    payloads = {
        "figure6.md": markdown,
        "figure6.json": raw_bytes,
        "benchmark-claims.json": payload_bytes,
        "extra.json": b"{}\n",
    }
    for name, data in payloads.items():
        (result / name).write_bytes(data)
    argv = ("python", "-m", "benchmarks.cli", "figure6", "--out", str(result), *module.FIGURE6_SEMANTIC_COMMAND[3:])
    envelope = artifact.build_envelope(
        artifact_kind=artifact.ArtifactKind(module.BENCHMARK_ARTIFACT_KIND),
        source=source,
        started=started,
        finished=started + datetime.timedelta(seconds=1),
        argv=argv,
        input_version=artifact.IdentityVersion(module.BENCHMARK_INPUT_VERSION),
        input_payload=semantic,
        result_version=artifact.IdentityVersion(module.BENCHMARK_RESULT_VERSION),
        payloads=payloads,
    )
    artifact.write_envelope(result, envelope)

    with pytest.raises(module.InvalidMeasurementClaimPayloadError, match="result_identity.payloads"):
        module.validate_benchmark_claim_artifact(result)


@pytest.mark.parametrize("field", ["semantic_command", "config", "claim_sources"])
def test_benchmark_semantic_projection_uses_returned_payload_fields(field: str) -> None:
    module = _module()
    payload = _payload(module)
    if field == "semantic_command":
        payload["semantic_command"] = ["mutated"]
    elif field == "config":
        payload["config"]["width"] = 99.0
    else:
        payload["claims"][0]["source"] = "mutated.py:1"

    projected = module.benchmark_payload_semantic_input(payload)

    if field == "claim_sources":
        assert projected[field][0]["source"] == "mutated.py:1"
    else:
        assert projected[field] == payload[field]


def test_benchmark_ledger_evidence_projects_authenticated_command_and_config() -> None:
    module = _module()
    validation = importlib.import_module("tools.measurement_claim_benchmark_validation")
    artifact = importlib.import_module("tools.measurement_artifact")
    payload = _payload(module)
    envelope = artifact.ValidatedEnvelope.build(
        finished=datetime.datetime(2026, 8, 29, 12, 0, 1, tzinfo=datetime.timezone.utc),
        commit=SOURCE,
        input_sha256="b" * 64,
        result_sha256="c" * 64,
        payload_sha256={"figure6.md": "d" * 64, "figure6.json": "e" * 64, "benchmark-claims.json": "f" * 64},
    )
    cells = validation.render_benchmark_ledger_evidence(
        payload,
        envelope,
        started=datetime.datetime(2026, 8, 29, 12, 0, 0, tzinfo=datetime.timezone.utc),
        artifact_directory=pathlib.PurePosixPath("benchmarks/measurement_claim_results/2026-08-29-aaaaaaaaaaaa-benchmark-bbbbbbbbbbbb"),
    )
    expected_command = json.dumps(payload["semantic_command"], sort_keys=True, separators=(",", ":"), allow_nan=False)
    expected_config = json.dumps(payload["config"], sort_keys=True, separators=(",", ":"), allow_nan=False)
    assert all(f"semantic_command={expected_command}" in cell for cell in cells.values())
    assert all(f"config={expected_config}" in cell for cell in cells.values())
    assert all("\n" not in cell and "\r" not in cell and "|" not in cell for cell in cells.values())


def test_figure6_world_point_contract_is_strict_mypy_checked(tmp_path: pathlib.Path) -> None:
    positive = tmp_path / "positive.py"
    negative = tmp_path / "negative.py"
    positive.write_text(
        """from tools.measurement_claim_benchmark_schema import Figure6ConfigPayload
from tools.measurement_claim_schema import Degrees, Millimetres, ToolDiameters, WorldMillimetres, WorldPointMillimetres

point: WorldPointMillimetres = (
    WorldMillimetres(0.0),
    WorldMillimetres(1.0),
    WorldMillimetres(2.0),
)
config: Figure6ConfigPayload = {
    "width": Millimetres(20.0),
    "height": Millimetres(12.0),
    "tool_diameter": Millimetres(2.0),
    "holes": [[point]],
    "reporting_cap": Degrees(180.0),
    "caps": [Degrees(20.0)],
    "spacings": [ToolDiameters(0.1)],
    "length_unit": "mm",
    "angle_unit": "degree",
    "spacing_unit": "tool-diameter",
}
""",
        encoding="utf-8",
    )
    negative.write_text(
        """from tools.measurement_claim_schema import Millimetres, WorldMillimetres, WorldPointMillimetres

too_short: WorldPointMillimetres = (WorldMillimetres(0.0), WorldMillimetres(1.0))
too_long: WorldPointMillimetres = (
    WorldMillimetres(0.0),
    WorldMillimetres(1.0),
    WorldMillimetres(2.0),
    WorldMillimetres(3.0),
)
wrong_unit: WorldPointMillimetres = (
    WorldMillimetres(0.0),
    WorldMillimetres(1.0),
    Millimetres(2.0),
)
""",
        encoding="utf-8",
    )
    positive_result = subprocess.run(
        [sys.executable, "-m", "mypy", "--strict", "--warn-unused-ignores", str(positive)],
        cwd=PROJECT_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    negative_result = subprocess.run(
        [sys.executable, "-m", "mypy", "--strict", "--warn-unused-ignores", str(negative)],
        cwd=PROJECT_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    assert positive_result.returncode == 0, positive_result.stdout + positive_result.stderr
    assert negative_result.returncode != 0
    diagnostics = negative_result.stdout + negative_result.stderr
    assert diagnostics.count("Incompatible types in assignment") == 3
