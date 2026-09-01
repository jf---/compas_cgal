from typing import assert_type

from benchmarks.held_reference_cases import Degree
from benchmarks.held_reference_cases import Figure7Observation
from benchmarks.held_reference_cases import HeldReferenceCase
from benchmarks.held_reference_cases import load_all_held_reference_cases
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.spec import PocketSpec
from compas_cgal.adaptive.units import Point2
from compas_cgal.adaptive.units import WorldXY

case = assert_type(load_held_reference_case("figure5"), HeldReferenceCase)
assert_type(load_all_held_reference_cases(), tuple[HeldReferenceCase, ...])
assert_type(case.pocket_spec(), PocketSpec)
assert_type(case.tea_cap, Degree)
assert_type(case.start_marker, Point2[WorldXY] | None)
assert_type(case.figure7_observation, Figure7Observation | None)
assert_type(
    Figure7Observation.build(
        publication_page=744,
        pdf_page=14,
        panels=("a", "b", "c"),
        role="shape_only_tool_centre_observation",
        boundary_authority=False,
        numeric_fidelity_gate=False,
    ),
    Figure7Observation,
)
assert_type(
    HeldReferenceCase.build(
        name=case.name,
        authors=case.authors,
        title=case.title,
        doi=case.doi,
        publication_page=case.publication_page,
        pdf_page=case.pdf_page,
        figure=case.figure,
        subfigure=case.subfigure,
        source_crop=case.source_crop,
        source_tool_centre=case.source_tool_centre,
        source_tool_radius=case.source_tool_radius,
        boundary=case.boundary,
        reconstruction=case.reconstruction,
        projection=case.projection,
        tool_radius=case.tool_radius,
        tea_cap=case.tea_cap,
        start_marker=case.start_marker,
        start_marker_radius=case.start_marker_radius,
        figure7_observation=case.figure7_observation,
    ),
    HeldReferenceCase,
)
