from __future__ import annotations

import math
from dataclasses import dataclass
from io import BytesIO
from typing import Any
from typing import List
from typing import Optional
from typing import Tuple

import matplotlib
import pytest
from compas.geometry import Arc
from compas.geometry import Circle
from compas.geometry import Frame
from compas.geometry import Line
from compas.geometry import Polygon

# No display anywhere in this suite: the drawing module builds its own Agg canvas
# and never imports pyplot, and this pins the global backend as well so an
# accidental pyplot import in a future test cannot start a window server.
matplotlib.use("Agg")

from matplotlib.backends.backend_agg import FigureCanvasAgg  # noqa: E402
from matplotlib.colors import to_rgba  # noqa: E402

from benchmarks.errors import AmbiguousPanelError  # noqa: E402
from benchmarks.errors import EmptyComparisonError  # noqa: E402
from benchmarks.errors import EmptyToolpathError  # noqa: E402
from benchmarks.errors import EngagementLengthMismatchError  # noqa: E402
from benchmarks.errors import InvalidBandCountError  # noqa: E402
from benchmarks.errors import MissingEngagementDataError  # noqa: E402
from benchmarks.errors import MissingToolDiameterError  # noqa: E402
from benchmarks.errors import UnknownOperationClassError  # noqa: E402
from benchmarks.errors import UnplottableBoundaryError  # noqa: E402
from benchmarks.errors import UnplottableGeometryError  # noqa: E402
from benchmarks.palette import DARK  # noqa: E402
from benchmarks.palette import LIGHT  # noqa: E402
from benchmarks.palette import Theme  # noqa: E402
from benchmarks.palette import band_edges  # noqa: E402
from benchmarks.palette import band_index  # noqa: E402
from benchmarks.marks import CUT_WIDTH_PT  # noqa: E402
from benchmarks.pathgeometry import ARC_DEGREES_PER_SAMPLE  # noqa: E402
from benchmarks.plotting import POINTS_PER_INCH  # noqa: E402
from benchmarks.plotting import ColourBy  # noqa: E402
from benchmarks.plotting import draw_comparison  # noqa: E402
from benchmarks.plotting import draw_toolpath  # noqa: E402

# The pocket every drawing in this module is made over: ten tool diameters by six,
# the same shape `benchmarks.figure6` reports on, so a figure drawn here has the
# aspect ratio the real ones do.
POCKET_WIDTH = 20.0
POCKET_HEIGHT = 12.0
TOOL_DIAMETER = 2.0

# Machining circles are drawn at this radius; small enough that several fit along
# the pocket, large enough that a chord would be visible if one were drawn.
CIRCLE_RADIUS = 1.5

# Height a rapid retracts to. Only its existence matters here, never its value.
CLEARANCE_Z = 4.0


@dataclass(frozen=True)
class FakeOperation:
    """A toolpath operation with nothing but the three fields a drawing reads."""

    geometry: Any
    operation: str
    path_index: int


@dataclass(frozen=True)
class FakeResult:
    """A toolpath result carrying only its operation stream."""

    operations: List[FakeOperation]


class ExplodingPolylineResult:
    """A result whose tessellated polyline cannot be read without failing the test."""

    def __init__(self, operations: List[FakeOperation]) -> None:
        self.operations = operations

    @property
    def polyline(self) -> Any:
        raise AssertionError("The drawing read `result.polyline`; arcs must be sampled from their own parametrisation.")


def _pocket() -> Polygon:
    """The rectangular pocket every drawing is made over."""
    return Polygon([[0.0, 0.0, 0.0], [POCKET_WIDTH, 0.0, 0.0], [POCKET_WIDTH, POCKET_HEIGHT, 0.0], [0.0, POCKET_HEIGHT, 0.0]])


def _circle(cx: float, cy: float, radius: float = CIRCLE_RADIUS) -> Circle:
    """A machining circle centred at *cx*, *cy* in the cutting plane."""
    return Circle(radius=radius, frame=Frame([cx, cy, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]))


def _arc(cx: float, cy: float, start: float, end: float) -> Arc:
    """An open arc centred at *cx*, *cy*, sweeping from *start* to *end* radians."""
    return Arc(radius=CIRCLE_RADIUS, start_angle=start, end_angle=end, frame=Frame([cx, cy, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]))


def _line(x0: float, y0: float, x1: float, y1: float, z0: float = 0.0, z1: float = 0.0) -> Line:
    """A straight motion between two stations."""
    return Line([x0, y0, z0], [x1, y1, z1])


def _chain(path_index: int, x: float, y: float) -> List[FakeOperation]:
    """One machining chain: plunge, two circles joined by a cut, retract."""
    return [
        FakeOperation(_line(x, y, x, y, CLEARANCE_Z, 0.0), "plunge", path_index),
        FakeOperation(_circle(x, y), "cut", path_index),
        FakeOperation(_line(x, y, x + 2.0, y), "cut", path_index),
        FakeOperation(_circle(x + 2.0, y), "cut", path_index),
        FakeOperation(_line(x + 2.0, y, x + 2.0, y, 0.0, CLEARANCE_Z), "retract", path_index),
    ]


def _path(chains: int = 2) -> FakeResult:
    """A path of *chains* chains, linked at clearance height."""
    operations: List[FakeOperation] = []
    exit_point: Optional[Tuple[float, float]] = None
    for index in range(chains):
        x = 2.0 + 4.0 * (index % 4)
        y = 3.0 + 3.0 * (index // 4)
        if exit_point is not None:
            operations.append(FakeOperation(_line(exit_point[0], exit_point[1], x, y, CLEARANCE_Z, CLEARANCE_Z), "link", index))
        operations.extend(_chain(index, x, y))
        exit_point = (x + 2.0, y)
    return FakeResult(operations=operations)


def _engagement(result: FakeResult) -> List[Optional[float]]:
    """One plausible measurement per operation, None for everything that never cuts."""
    return [float(30 + 10 * index) if operation.operation == "cut" else None for index, operation in enumerate(result.operations)]


def _strokes(axis: Any, colour: str) -> List[Any]:
    """Every drawn line of one colour."""
    return [line for line in axis.lines if line.get_color() == colour]


def _stroke_colours(axis: Any) -> Tuple[str, ...]:
    """Every colour drawn as a stroke, deduplicated."""
    return tuple({line.get_color() for line in axis.lines if line.get_linestyle() != "None"})


def test_every_colour_mode_draws_a_legend_and_saves_without_a_display() -> None:
    result = _path()
    for mode in ColourBy:
        drawing = draw_toolpath(
            result,
            boundary=_pocket(),
            colour_by=mode,
            engagement_deg=_engagement(result) if mode is ColourBy.ENGAGEMENT else None,
            title="pocket",
        )
        assert len(drawing.axes) == 1, mode
        assert len(drawing.legend_labels) >= 2, mode
        assert isinstance(drawing.figure.canvas, FigureCanvasAgg), mode
        buffer = BytesIO()
        drawing.figure.savefig(buffer, format="svg")
        assert buffer.getvalue().startswith(b"<?xml"), mode


def test_engagement_mode_without_measurements_raises_its_own_error() -> None:
    with pytest.raises(MissingEngagementDataError):
        draw_toolpath(_path(), boundary=_pocket(), colour_by=ColourBy.ENGAGEMENT)


def test_engagement_mode_rejects_a_measurement_list_of_the_wrong_length() -> None:
    result = _path()
    with pytest.raises(EngagementLengthMismatchError):
        draw_toolpath(result, boundary=_pocket(), colour_by=ColourBy.ENGAGEMENT, engagement_deg=[10.0, 20.0])


def test_engagement_bands_span_the_supplied_range() -> None:
    result = _path()
    values = _engagement(result)
    drawing = draw_toolpath(result, boundary=_pocket(), colour_by=ColourBy.ENGAGEMENT, engagement_deg=values)
    measured = [value for value in values if value is not None]
    assert drawing.legend_labels[0].startswith(f"{min(measured):.0f}")
    assert any(label.endswith(f"–{max(measured):.0f}°") for label in drawing.legend_labels)
    assert "not measured" in drawing.legend_labels


def test_a_ninth_traversal_folds_instead_of_cycling_hues() -> None:
    result = _path(chains=12)
    drawing = draw_toolpath(result, boundary=_pocket(), colour_by=ColourBy.TRAVERSAL)
    drawn = set(_stroke_colours(drawing.axis))
    assert drawn <= {*LIGHT.ramp, LIGHT.muted, LIGHT.ink}
    assert len([colour for colour in LIGHT.ramp if colour in drawn]) == len(LIGHT.ramp)
    folded = [label for label in drawing.legend_labels if label.startswith("other (")]
    assert folded == [f"other ({12 - len(LIGHT.ramp)} chains)"]


def test_a_path_within_the_ramp_capacity_never_folds() -> None:
    drawing = draw_toolpath(_path(chains=3), boundary=_pocket(), colour_by=ColourBy.TRAVERSAL)
    assert [label for label in drawing.legend_labels if label.startswith("other (")] == []
    assert "chain 2" in drawing.legend_labels


def test_a_circle_is_drawn_as_an_arc_not_as_one_segment_per_operation() -> None:
    operations = [FakeOperation(_circle(POCKET_WIDTH / 2.0, POCKET_HEIGHT / 2.0), "cut", 0)]
    result = FakeResult(operations=operations)
    drawing = draw_toolpath(result, boundary=_pocket())
    drawn = _strokes(drawing.axis, LIGHT.cut)
    assert len(drawn) == 1
    vertices = len(drawn[0].get_xydata())
    assert vertices > len(result.operations)
    assert vertices == int(360.0 / ARC_DEGREES_PER_SAMPLE) + 1


def test_an_arc_is_sampled_in_proportion_to_its_sweep() -> None:
    quarter = FakeResult(operations=[FakeOperation(_arc(10.0, 6.0, 0.0, math.pi / 2.0), "cut", 0)])
    drawing = draw_toolpath(quarter, boundary=_pocket())
    assert len(_strokes(drawing.axis, LIGHT.cut)[0].get_xydata()) == int(90.0 / ARC_DEGREES_PER_SAMPLE) + 1


def test_the_tessellated_polyline_is_never_read() -> None:
    result = ExplodingPolylineResult(operations=list(_path().operations))
    drawing = draw_toolpath(result, boundary=_pocket())
    assert len(drawing.axes) == 1


def test_plunges_and_retracts_are_shaped_point_markers_not_strokes() -> None:
    drawing = draw_toolpath(_path(chains=1), boundary=_pocket())
    markers = {line.get_marker() for line in drawing.axis.lines if line.get_linestyle() == "None"}
    assert markers == {"v", "^"}
    assert "plunge" in drawing.legend_labels
    assert "retract" in drawing.legend_labels
    for line in drawing.axis.lines:
        if line.get_marker() in {"v", "^"}:
            assert line.get_color() == LIGHT.secondary


def test_the_tool_envelope_needs_a_diameter_to_draw() -> None:
    with pytest.raises(MissingToolDiameterError):
        draw_toolpath(_path(), boundary=_pocket(), show_tool_envelope=True)


def test_the_tool_envelope_is_a_width_in_data_units_and_survives_a_resize() -> None:
    drawing = draw_toolpath(_path(), boundary=_pocket(), tool_diameter=TOOL_DIAMETER, show_tool_envelope=True)
    axis, figure = drawing.axis, drawing.figure
    envelope = [line for line in _strokes(axis, LIGHT.grid) if line.get_linestyle() != "None"]
    assert envelope

    def expected() -> float:
        origin = axis.transData.transform((0.0, 0.0))
        offset = axis.transData.transform((TOOL_DIAMETER, 0.0))
        return math.hypot(float(offset[0] - origin[0]), float(offset[1] - origin[1])) * POINTS_PER_INCH / float(figure.dpi)

    figure.canvas.draw()
    before = expected()
    assert envelope[0].get_linewidth() == pytest.approx(before)
    assert before > CUT_WIDTH_PT

    figure.set_size_inches(figure.get_figwidth() / 2.0, figure.get_figheight() / 2.0)
    figure.canvas.draw()
    assert envelope[0].get_linewidth() == pytest.approx(expected())
    assert envelope[0].get_linewidth() < before


def test_a_comparison_stacks_panels_that_share_one_view_and_one_legend() -> None:
    panels = {"unregulated": _path(chains=2), "controlled": _path(chains=3)}
    drawing = draw_comparison(
        panels,
        boundary=_pocket(),
        colour_by=ColourBy.TRAVERSAL,
        panel_subtitles={"unregulated": "10 ops"},
        title="two paths",
        subtitle="one pocket",
    )
    assert len(drawing.axes) == 2
    first, second = drawing.axes
    assert first.get_xlim() == second.get_xlim()
    assert first.get_ylim() == second.get_ylim()
    assert first.get_position().y0 > second.get_position().y0
    assert len(drawing.figure.legends) == 1
    assert tuple(text.get_text() for text in drawing.figure.legends[0].get_texts()) == drawing.legend_labels
    assert "chain 2" in drawing.legend_labels


def test_a_comparison_needs_at_least_one_panel() -> None:
    with pytest.raises(EmptyComparisonError):
        draw_comparison({}, boundary=_pocket())


def test_a_comparison_in_engagement_mode_needs_every_panel_measured() -> None:
    result = _path()
    with pytest.raises(MissingEngagementDataError):
        draw_comparison(
            {"a": result, "b": result},
            boundary=_pocket(),
            colour_by=ColourBy.ENGAGEMENT,
            engagement_deg={"a": _engagement(result)},
        )


def test_the_single_panel_shortcut_refuses_a_multi_panel_drawing() -> None:
    drawing = draw_comparison({"a": _path(), "b": _path()}, boundary=_pocket())
    with pytest.raises(AmbiguousPanelError):
        drawing.axis


def test_the_dark_variant_is_the_dark_steps_and_not_an_inversion() -> None:
    drawing = draw_toolpath(_path(), boundary=_pocket(), theme=Theme.DARK)
    assert drawing.figure.get_facecolor() == to_rgba(DARK.surface)
    assert _strokes(drawing.axis, DARK.cut)
    assert _strokes(drawing.axis, DARK.ink)
    assert not _strokes(drawing.axis, LIGHT.cut)
    for text in drawing.figure.legends[0].get_texts():
        assert text.get_color() == DARK.ink


def test_a_geometry_panel_carries_no_axis_furniture_and_an_equal_aspect() -> None:
    drawing = draw_toolpath(_path(), boundary=_pocket())
    axis = drawing.axis
    assert axis.axison is False
    assert axis.get_aspect() == 1.0
    assert axis.get_facecolor() == to_rgba(LIGHT.surface)


def test_the_boundary_is_drawn_in_ink_and_closed() -> None:
    holes = [Polygon([[8.0, 5.0, 0.0], [12.0, 5.0, 0.0], [12.0, 7.0, 0.0], [8.0, 7.0, 0.0]])]
    drawing = draw_toolpath(_path(), boundary=_pocket(), holes=holes)
    rings = _strokes(drawing.axis, LIGHT.ink)
    assert len(rings) == 2
    for ring in rings:
        points = ring.get_xydata()
        assert tuple(points[0]) == tuple(points[-1])


def test_traversal_labels_land_on_the_first_cutting_move_of_every_chain() -> None:
    result = _path(chains=3)
    drawing = draw_toolpath(result, boundary=_pocket(), colour_by=ColourBy.TRAVERSAL, annotate_traversals=True)
    labels = {text.get_text(): text.get_position() for text in drawing.axis.texts}
    assert set(labels) == {"0", "1", "2"}
    first_cut = next(operation for operation in result.operations if operation.operation == "cut")
    expected = first_cut.geometry.point_at(0.0)
    assert labels["0"] == pytest.approx((expected[0], expected[1]))


def test_an_empty_result_is_refused_rather_than_drawn_blank() -> None:
    with pytest.raises(EmptyToolpathError):
        draw_toolpath(FakeResult(operations=[]), boundary=_pocket())


def test_an_unknown_operation_class_is_refused_rather_than_drawn_as_a_cut() -> None:
    result = FakeResult(operations=[FakeOperation(_line(0.0, 0.0, 1.0, 1.0), "dwell", 0)])
    with pytest.raises(UnknownOperationClassError):
        draw_toolpath(result, boundary=_pocket())


def test_geometry_with_no_path_is_refused() -> None:
    result = FakeResult(operations=[FakeOperation(object(), "cut", 0)])
    with pytest.raises(UnplottableGeometryError):
        draw_toolpath(result, boundary=_pocket())


def test_a_boundary_of_two_points_is_refused() -> None:
    with pytest.raises(UnplottableBoundaryError):
        draw_toolpath(_path(), boundary=[[0.0, 0.0], [1.0, 1.0]])


def test_a_caller_supplied_axes_is_drawn_into_and_keeps_its_figure() -> None:
    from matplotlib.figure import Figure

    figure = Figure(figsize=(3.0, 2.0))
    axis = figure.add_subplot()
    drawing = draw_toolpath(_path(), boundary=_pocket(), ax=axis, title="borrowed")
    assert drawing.figure is figure
    assert drawing.axes == (axis,)
    assert axis.get_legend() is not None
    buffer = BytesIO()
    drawing.figure.savefig(buffer, format="svg")
    assert buffer.getvalue()


def test_the_palette_is_the_documented_instance_and_not_an_eyeballed_one() -> None:
    assert (LIGHT.surface, LIGHT.ink, LIGHT.secondary, LIGHT.muted, LIGHT.grid) == ("#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9")
    assert (LIGHT.cut, LIGHT.link, LIGHT.lead) == ("#2a78d6", "#eb6834", "#1baf7a")
    assert LIGHT.ramp == ("#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281")
    assert (DARK.surface, DARK.ink, DARK.secondary, DARK.muted, DARK.grid) == ("#1a1a19", "#ffffff", "#c3c2b7", "#898781", "#2c2c2a")
    assert (DARK.cut, DARK.link, DARK.lead) == ("#3987e5", "#d95926", "#199e70")
    assert DARK.ramp == ("#184f95", "#256abf", "#3987e5", "#6da7ec", "#9ec5f4")
    assert len(LIGHT.ramp) == len(DARK.ramp)


def test_a_degenerate_range_puts_every_value_in_one_band() -> None:
    assert band_edges(5.0, 5.0, 5) == ((5.0, 5.0),)
    assert band_index(5.0, 5.0, 5.0, 5) == 0


def test_a_band_index_is_clamped_to_the_ramp() -> None:
    assert band_index(-100.0, 0.0, 10.0, 5) == 0
    assert band_index(100.0, 0.0, 10.0, 5) == 4
    assert band_index(10.0, 0.0, 10.0, 5) == 4


def test_a_ramp_needs_at_least_one_band() -> None:
    with pytest.raises(InvalidBandCountError):
        band_edges(0.0, 1.0, 0)
    with pytest.raises(InvalidBandCountError):
        band_index(0.5, 0.0, 1.0, 0)
