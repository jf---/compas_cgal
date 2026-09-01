from __future__ import annotations

from pathlib import Path

from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.held_reference_figures import RasterPoint2
from benchmarks.held_reference_figures import measure_figure7_observation
from benchmarks.held_reference_figures import measure_figure7_panels


def test_live_figure7_panels_close_the_approved_600_dpi_page_grammar(held_publisher_pdf: Path) -> None:
    figure5 = load_held_reference_case("figure5")

    panels = measure_figure7_panels(held_publisher_pdf, figure5)

    assert tuple(panel.panel for panel in panels) == ("a", "b", "c")
    assert tuple((panel.rgb.shape[1], panel.rgb.shape[0]) for panel in panels) == (
        (1567, 1092),
        (1567, 1092),
        (1568, 1092),
    )
    assert tuple((panel.colour_support.minimum, panel.colour_support.maximum) for panel in panels) == (
        (RasterPoint2.build(476, 740), RasterPoint2.build(1970, 1685)),
        (RasterPoint2.build(2421, 740), RasterPoint2.build(3915, 1685)),
        (RasterPoint2.build(1449, 2087), RasterPoint2.build(2943, 3032)),
    )
    assert all(panel.registration.x_scale != panel.registration.y_scale for panel in panels)
    assert measure_figure7_observation(held_publisher_pdf, figure5) == figure5.figure7_observation
