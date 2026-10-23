"""The four coverage-certificate records load from their committed JSON and re-render byte for byte.

Each tool (acceptance, circle core, cost model, raster) writes its record with
``to_json`` and renders it with ``render_markdown``. ``load_record`` is the
inverse of ``to_json``: it refuses a record of another schema with the tool's
named stale-record error (re-measure, never render), and on the committed
record it reproduces both the JSON and the docs Markdown exactly.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType

import pytest

from tools import coverage_certificate_acceptance
from tools import coverage_certificate_circle_core
from tools import coverage_certificate_cost_model
from tools import coverage_certificate_raster


@dataclass(frozen=True)
class RecordTool:
    """One tool under test and the name of its stale-record error."""

    module: ModuleType
    stale_error: str


TOOLS = (
    RecordTool(coverage_certificate_acceptance, "StaleAcceptanceRecordError"),
    RecordTool(coverage_certificate_circle_core, "StaleCircleCoreRecordError"),
    RecordTool(coverage_certificate_cost_model, "StaleCostModelRecordError"),
    RecordTool(coverage_certificate_raster, "StaleRasterRecordError"),
)
IDS = [tool.module.__name__.rsplit(".", 1)[1] for tool in TOOLS]


@pytest.mark.parametrize("tool", TOOLS, ids=IDS)
def test_committed_record_is_at_the_current_schema(tool: RecordTool) -> None:
    raw = json.loads(tool.module.JSON_PATH.read_text(encoding="utf-8"))
    assert raw["schema"] == tool.module.SCHEMA


@pytest.mark.parametrize("tool", TOOLS, ids=IDS)
def test_load_record_inverts_to_json(tool: RecordTool) -> None:
    text = tool.module.JSON_PATH.read_text(encoding="utf-8")
    assert tool.module.load_record(tool.module.JSON_PATH).to_json() == text


@pytest.mark.parametrize("tool", TOOLS, ids=IDS)
def test_committed_markdown_is_the_render_of_the_loaded_record(tool: RecordTool) -> None:
    record = tool.module.load_record(tool.module.JSON_PATH)
    assert tool.module.render_markdown(record).encode("utf-8") == tool.module.MD_PATH.read_bytes()


@pytest.mark.parametrize("tool", TOOLS, ids=IDS)
def test_a_record_of_another_schema_is_refused(tool: RecordTool, tmp_path: Path) -> None:
    raw = json.loads(tool.module.JSON_PATH.read_text(encoding="utf-8"))
    raw["schema"] = f"{tool.module.SCHEMA}-older"
    stale = tmp_path / "stale.json"
    stale.write_text(json.dumps(raw), encoding="utf-8")
    error = getattr(tool.module, tool.stale_error)
    with pytest.raises(error, match="re-measure"):
        tool.module.load_record(stale)
