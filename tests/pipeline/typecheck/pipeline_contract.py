"""Static contract for the pipeline's public API; checked by `pixi run types-pipeline`, never executed."""

from typing import assert_type

from compas.geometry import Polygon

from compas_cgal.incident_loop import IncidentLoop
from compas_cgal.pipeline.admissibility import Admissibility
from compas_cgal.pipeline.admissibility import CapChordRatio
from compas_cgal.pipeline.admissibility import cap_ratio
from compas_cgal.pipeline.guide import Guide
from compas_cgal.pipeline.guide import GuideChain
from compas_cgal.pipeline.guide import GuideStation
from compas_cgal.pipeline.guide import pipeline_spacing
from compas_cgal.pipeline.ladder import NoAdmissibleCandidate
from compas_cgal.pipeline.ladder import furthest_first
from compas_cgal.pipeline.narrow import Advance
from compas_cgal.pipeline.narrow import ChainSplit
from compas_cgal.pipeline.narrow import OrientedChain
from compas_cgal.pipeline.narrow import RadiusStepCondition
from compas_cgal.pipeline.narrow import advance_ladder
from compas_cgal.pipeline.narrow import station_loop
from compas_cgal.pipeline.generate import generate_narrow
from compas_cgal.pipeline.linking import BeltLink
from compas_cgal.pipeline.linking import NoLink
from compas_cgal.pipeline.linking import SameCircleArc
from compas_cgal.pipeline.linking import StraightLink
from compas_cgal.pipeline.linking import link_between
from compas_cgal.pipeline.result import DeclinedRun
from compas_cgal.pipeline.result import ForcedDecisionLedger
from compas_cgal.pipeline.result import ForcedEntry
from compas_cgal.pipeline.result import ForcedEntryReason
from compas_cgal.pipeline.result import LinkRecord
from compas_cgal.pipeline.result import OmittedLoop
from compas_cgal.pipeline.result import PipelineResult
from compas_cgal.pipeline.result import ForcedAdvance
from compas_cgal.stock import Stock
from compas_cgal.units import ClearanceZ
from compas_cgal.units import GuideRadius
from compas_cgal.units import Radian
from compas_cgal.units import Spacing
from compas_cgal.units import ToolRadius


def guide_contract(guide: Guide) -> None:
    chain = guide.chains[0]
    assert_type(chain, GuideChain)
    assert_type(chain.reversed(), GuideChain)
    station = chain.stations[0]
    assert_type(station, GuideStation)
    assert_type(station.radius, GuideRadius)
    assert_type(pipeline_spacing(2.0, Radian(1.0)), Spacing)


def ladder_contract(stock: Stock, chain: OrientedChain, tool: ToolRadius, admissibility: Admissibility) -> None:
    def ints(_stock: Stock, _origin: int, candidate: int, /) -> bool:
        return candidate > 0

    assert_type(furthest_first(stock, 0, (2, 1), (ints,)), int | NoAdmissibleCandidate)
    advance = advance_ladder(stock, chain, chain.first, chain.last, (RadiusStepCondition(tool),))
    assert_type(advance, Advance | ChainSplit)
    if isinstance(advance, Advance):
        assert_type(advance.forced, ForcedAdvance | None)
    else:
        assert_type(advance.gap, DeclinedRun)
    assert_type(chain.declined, tuple[DeclinedRun, ...])
    assert_type(station_loop(chain.entry_station, tool), IncidentLoop)
    loop = station_loop(chain.entry_station, tool)
    assert_type(admissibility.loop_within_cap(stock, loop, True), bool)
    assert_type(admissibility.cap_ratio, CapChordRatio)
    assert_type(cap_ratio(Radian(1.0)), CapChordRatio)


def link_contract(a: IncidentLoop, b: IncidentLoop) -> None:
    link = link_between(a, b, True, True)
    assert_type(link, SameCircleArc | BeltLink | StraightLink | NoLink)
    if not isinstance(link, NoLink):
        assert_type(link.next_loop, IncidentLoop)


def generate_contract(polygon: Polygon, tool: ToolRadius, clearance: ClearanceZ) -> None:
    result = generate_narrow(polygon, [], tool, Radian(1.0), clearance)
    assert_type(result, PipelineResult)
    assert_type(result.ledger, ForcedDecisionLedger)
    assert_type(result.ledger.entries, tuple[ForcedEntry, ...])
    assert_type(result.omitted, tuple[OmittedLoop, ...])
    assert_type(result.links, tuple[LinkRecord, ...])
    assert_type(result.ledger.entries[0].reason, ForcedEntryReason)
