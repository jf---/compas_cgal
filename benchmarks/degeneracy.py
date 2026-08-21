"""Hand-built degenerate instances — the expected sources of `unresolved`.

Exactness makes degeneracies RELIABLE, not absent. Each instance here targets one
documented ordinary-not-exceptional case from the engagement kernel's inventory,
so an `unresolved` regression is attributable to a named cause instead of being
noticed as a statistic. Every instance is machinable: an instance the generator
declines to enter produces no operations and therefore no evidence, which is the
one thing a degeneracy corpus must never contain.

Three constructions were rejected on measurement rather than on taste
(2026-08-21, tool 1.0, cap 120 deg):

* **A repeated vertex cannot be authored here.** `compas.geometry.Polygon`
  collapses a repeated vertex on construction -- and collapses a pair one ULP
  apart too, so it is a proximity rule, not exact equality. A `duplicate_vertex`
  instance would be a bit-identical copy of the plain box wearing a name that
  claims otherwise. `cocircular_square` takes its place; the omission is pinned
  by a test so that a future compas which stops collapsing is noticed.
* **A standalone slot exactly 2r wide yields ZERO operations.** It survives
  `PocketSpec.build` (area 20 against a 1.0 floor), so nothing rejects it -- it
  simply produces an empty toolpath, and an empty instance contributes nothing to
  any statistic while diluting every average. `half_turn_arm` keeps the exact
  half-turn geometry and attaches it to a pocket the cutter can reach it from:
  405 operations for the bare box, 525 with the arm, so the arm is machined.
* **Zero-measure contact BETWEEN rings is outside the stock model.** Two islands
  touching at a corner raise "Hole polygons must be pairwise disjoint", and an
  island edge lying on the outer wall raises "Hole polygon must lie strictly
  inside the outer boundary". Those are `Stock`'s own documented invariants and
  the corpus respects them rather than probing behind them.

One error in the source plan is corrected here: it placed the tangent island two
tool DIAMETERS from the wall while claiming simultaneous tangency, which needs
one. A cutter of radius r centred at x = r touches x = 0 and x = 2r together.

MEASURED (2026-08-21, tool 1.0, cap 120 deg, via `benchmarks.runner.run_spec`)
-- every instance machines, and none of them errors:

    instance             operations   cut   uncertified   truly_exceeding
    tangent_island              487   424           235                89
    collinear_run               411   395           205                68
    half_turn_arm               525   435           261                43
    pinch_exactly_tool          705   635           354                94
    cocircular_square           335   324           168                63

Certification runs 15-25 s per instance at this tool size, so the corpus is a
report-time cost, not a test-time one; the tests in `tests/benchmarks` only build
the specs and assert their exact defining quantities.
"""

from __future__ import annotations

from typing import Mapping

from compas.geometry import Polygon

from benchmarks.spec import PocketSpec

# Outer pocket used by most degenerate instances.
BOX_W = 20.0
BOX_H = 12.0

# Island extent, and where the collinear run and the neck sit along the box.
ISLAND_SIDE = 4.0
ISLAND_LOW_Y = 4.0
COLLINEAR_MID_X = 10.0
NECK_HALF_SPAN = 1.0

# Length of the exactly-inscribed arm. Long enough to carry several cutter
# positions at the half-turn engagement rather than a single tangential touch.
ARM_LENGTH = 6.0

# The corpus's instance names, in emission order. Exported so a report or a test
# can assert the corpus's membership without rebuilding its geometry.
DEGENERACY_NAMES: tuple[str, ...] = (
    "tangent_island",
    "collinear_run",
    "half_turn_arm",
    "pinch_exactly_tool",
    "cocircular_square",
)


def degeneracy_corpus(tool_diameter: float, tea_cap_deg: float) -> list[PocketSpec]:
    """Every degenerate instance, each targeting one named kernel case.

    Args:
        tool_diameter: Cutter diameter. Every instance's defining feature is
            expressed in it, so the corpus stays exact at any tool size.
        tea_cap_deg: Engagement cap in degrees.

    Returns:
        The corpus, in `DEGENERACY_NAMES` order.

    Raises:
        DegeneratePocketError: A tool large enough to swallow the fixed box was
            requested, so an instance no longer admits it.
        NonPositiveToolError: The diameter is NaN, zero, or negative.
        InvalidCapError: The cap is NaN or outside (0, 180].
    """
    r = 0.5 * tool_diameter
    mid_y = 0.5 * BOX_H
    box = Polygon([[0.0, 0.0, 0.0], [BOX_W, 0.0, 0.0], [BOX_W, BOX_H, 0.0], [0.0, BOX_H, 0.0]])

    def spec(name: str, polygon: Polygon, params: Mapping[str, float], holes: tuple[Polygon, ...] = ()) -> PocketSpec:
        return PocketSpec.build(name=name, family="degeneracy", polygon=polygon, tool_diameter=tool_diameter, tea_cap_deg=tea_cap_deg, holes=holes, params=params)

    # The gap from the wall to the island is EXACTLY the tool diameter, so the
    # cutter rim meets both at once: a zero-measure contact on each side, decided
    # by an exact predicate rather than by which one a tolerance sees first.
    island_x = tool_diameter
    tangent_island = Polygon(
        [
            [island_x, ISLAND_LOW_Y, 0.0],
            [island_x, ISLAND_LOW_Y + ISLAND_SIDE, 0.0],
            [island_x + ISLAND_SIDE, ISLAND_LOW_Y + ISLAND_SIDE, 0.0],
            [island_x + ISLAND_SIDE, ISLAND_LOW_Y, 0.0],
        ]
    )

    # Three exactly collinear boundary vertices: CGAL::orientation returns
    # COLLINEAR, which is a verdict and not a near-miss.
    collinear = Polygon([[0.0, 0.0, 0.0], [COLLINEAR_MID_X, 0.0, 0.0], [BOX_W, 0.0, 0.0], [BOX_W, BOX_H, 0.0], [0.0, BOX_H, 0.0]])

    # An arm exactly 2r wide, opening off the box: the cutter is exactly inscribed
    # there, so its engaged run is an exact half turn -- the pi case the kernel
    # decides by orientation rather than by chord comparison. Attached to the box
    # because the same arm standing alone is never entered at all.
    arm = Polygon(
        [
            [0.0, 0.0, 0.0],
            [BOX_W, 0.0, 0.0],
            [BOX_W, mid_y - r, 0.0],
            [BOX_W + ARM_LENGTH, mid_y - r, 0.0],
            [BOX_W + ARM_LENGTH, mid_y + r, 0.0],
            [BOX_W, mid_y + r, 0.0],
            [BOX_W, BOX_H, 0.0],
            [0.0, BOX_H, 0.0],
        ]
    )

    # Neck exactly equal to the tool diameter: traversable with zero clearance.
    left_x, right_x = COLLINEAR_MID_X - NECK_HALF_SPAN, COLLINEAR_MID_X + NECK_HALF_SPAN
    pinch = Polygon(
        [
            [0.0, 0.0, 0.0],
            [left_x, 0.0, 0.0],
            [COLLINEAR_MID_X, mid_y - r, 0.0],
            [right_x, 0.0, 0.0],
            [BOX_W, 0.0, 0.0],
            [BOX_W, BOX_H, 0.0],
            [right_x, BOX_H, 0.0],
            [COLLINEAR_MID_X, mid_y + r, 0.0],
            [left_x, BOX_H, 0.0],
            [0.0, BOX_H, 0.0],
        ]
    )

    # A square: all four walls are equidistant from the centre, so the medial axis
    # carries ONE degree-4 node where a rectangle carries two degree-3 nodes. The
    # coincidence is exact, and the kernel must decide it as a coincidence rather
    # than cluster two nearly-equal circumradii.
    square = Polygon([[0.0, 0.0, 0.0], [BOX_H, 0.0, 0.0], [BOX_H, BOX_H, 0.0], [0.0, BOX_H, 0.0]])

    return [
        spec("tangent_island", box, {"wall_gap": tool_diameter}, (tangent_island,)),
        spec("collinear_run", collinear, {"collinear_cross": 0.0}),
        spec("half_turn_arm", arm, {"arm_width": tool_diameter}),
        spec("pinch_exactly_tool", pinch, {"neck_width": tool_diameter}),
        spec("cocircular_square", square, {"inradius": 0.5 * BOX_H}),
    ]
