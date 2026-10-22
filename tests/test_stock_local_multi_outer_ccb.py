"""Local depletion after a Gps boolean that misfiled a hole as an outer CCB.

In a bounded planar arrangement every bounded face has exactly one outer CCB.
CGAL's `Gps_on_surface_base_2::_remove_redundant_edges` breaks that: its last
tagging loop files every still-undecided cycle of a merged face as OUTER, so
after the global circle sweep plus `subtract_line_sweep_over` below, an island
of material inside the removed void is filed as a second outer CCB of the void
face (src/gps_ccb_repair_2.h).

The local update consumed that state: its flood walked only
`Face::outer_ccb()` (precondition compiled out of release builds) and
`Arrangement_2::remove_edge` then called `erase_inner_ccb(nullptr)` on the
misfiled cycle -- SIGSEGV, exit 139. The tangency of the last circle sweep with
the line sweep's end disk at one exact point only decides that the repair
reaches the misfiled cycle.

Each case runs in a child process so a hard fault is one failed test rather
than a dead pytest worker, and asserts the ANSWER -- exact point-set equality
with the global `Gps::difference` path, the Gps representation invariant and
the arrangement size -- not just clean termination.
"""

import subprocess
import sys
from fractions import Fraction

import pytest

from compas_cgal import _held_centre_domain_2

_SETUP = """
import numpy as np
from compas_cgal import _stock_2

POLY = np.array([(-10.0, 0.0, 0.0), (25.0, 0.0, 0.0), (25.0, 20.0, 0.0), (-10.0, 20.0, 0.0)])
PARAMS = [k / 46 for k in range(47)]


def replay(ops, ratio, local):
    stock = _stock_2.Stock2(POLY, [])
    for op in ops:
        if op[0] == "c":
            sweep = stock.subtract_circle_sweep_local if local else stock.subtract_circle_sweep
            sweep(op[1], op[2], op[3], 1.0)
        else:
            if ratio is not None:
                _stock_2.engagement_line_prefixes(stock, *op[1:], PARAMS, 1.0, ratio, 0.0)
            stock.subtract_line_sweep_over(*op[1:], 1.0)
    return stock


def check(ops, ratio):
    reference = replay(ops, ratio, local=False)
    local = replay(ops, ratio, local=True)
    assert local.exactly_equals(reference)
    assert reference.exactly_equals(local)
    assert local.representation_is_valid()
    assert local.arrangement_stats() == reference.arrangement_stats()
    print("ok")
"""

# Minimal: one global circle sweep, one non-Pythagorean line sweep whose end
# disk touches y = 0 at x = 0.6371..., then a local circle sweep whose outer
# circle touches y = 0 at that same exact point.
_MINIMAL = """
check([
    ("c", 0.0, 2.436736686061722, 1.436736686061722),
    ("p", 0.31855087196269893, 1.7, 0.6371017439253979, 1.0),
    ("c", 0.6371017439253979, 1.612007087599138, 0.6120070875991379),
], None)
"""

# The connector-engagement sequences that surfaced the defect, engagement
# queries included, every circle sweep on the local path.
_CONNECTOR_A = """
check([
    ("c", 0.0, 2.436736686061722, 1.436736686061722),
    ("p", 0.0, 1.0, 0.31855087196269893, 1.7),
    ("p", 0.31855087196269893, 1.7, 0.6371017439253979, 1.0),
    ("c", 0.6371017439253979, 1.612007087599138, 0.6120070875991379),
], 0.9690393558190852)
"""

_CONNECTOR_B = """
check([
    ("c", 0.0, 2.8545278605748545, 1.8545278605748545),
    ("p", 0.0, 1.0, 0.8771326995040981, 1.7),
    ("p", 0.8771326995040981, 1.7, 1.7542653990081962, 1.0),
    ("c", 1.7542653990081962, 1.492131634513695, 0.4921316345136949),
    ("p", 1.7542653990081962, 1.0, 3.5173154843182495, 1.0),
], 3.4812474752758726)
"""

# Child-process wall clock bound, seconds: each case is four exact sweeps on a
# rectangle (well under a second unloaded); the bound only turns a hang into a
# failure on a machine loaded to several hundred.
_CHILD_TIMEOUT_S = 600


@pytest.mark.parametrize(
    "body",
    [_MINIMAL, _CONNECTOR_A, _CONNECTOR_B],
    ids=["minimal_tangent_triple", "connector_A", "connector_B"],
)
def test_local_circle_sweep_after_line_sweep_matches_global(body):
    completed = subprocess.run(
        [sys.executable, "-c", _SETUP + body],
        capture_output=True,
        text=True,
        timeout=_CHILD_TIMEOUT_S,
    )
    assert completed.returncode == 0, f"child exited {completed.returncode}\nstdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    assert completed.stdout.strip().endswith("ok")


_REFILE_COUNTS = """
import numpy as np
from compas_cgal import _stock_2

POLY = np.array([(-10.0, 0.0, 0.0), (25.0, 0.0, 0.0), (25.0, 20.0, 0.0), (-10.0, 20.0, 0.0)])
stock = _stock_2.Stock2(POLY, [])
stock.subtract_circle_sweep(0.0, 2.436736686061722, 1.436736686061722, 1.0)
planar = stock.refile_misfiled_outer_ccbs()
stock.subtract_line_sweep_over(0.31855087196269893, 1.7, 0.6371017439253979, 1.0, 1.0)
reference = _stock_2.Stock2(POLY, [])
reference.subtract_circle_sweep(0.0, 2.436736686061722, 1.436736686061722, 1.0)
reference.subtract_line_sweep_over(0.31855087196269893, 1.7, 0.6371017439253979, 1.0, 1.0)
misfiled = stock.refile_misfiled_outer_ccbs()
again = stock.refile_misfiled_outer_ccbs()
assert stock.exactly_equals(reference) and reference.exactly_equals(stock)
assert stock.arrangement_stats() == reference.arrangement_stats()
print(planar, misfiled, again)
"""


def test_the_gps_misfile_is_present_and_refiled_exactly_once():
    """The circle sweep leaves a planar arrangement; the line sweep misfiles the island once.

    Pins that the repair is not dead code: a CGAL upgrade that fixes the misfile
    turns this red, and the re-filed set is the same point set and the same size.
    """
    completed = subprocess.run([sys.executable, "-c", _REFILE_COUNTS], capture_output=True, text=True, timeout=_CHILD_TIMEOUT_S)
    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.split() == ["0", "1", "0"]


def _wire(value):
    value = Fraction(value)
    return str(value.numerator), str(value.denominator)


def _line(start, end, normal):
    return (0, 0, [_wire(start[0]), _wire(start[1]), _wire(end[0]), _wire(end[1]), _wire(normal[0]), _wire(normal[1])])


def test_a_two_component_centre_domain_builds_with_a_planar_arrangement():
    """A dumbbell whose neck is narrower than 2r erodes into two Ω* components.

    The centre domain re-files misfiled outer CCBs right after its Gps
    difference and then requires exactly one outer CCB per contained face: two
    components must build, each with one outer cycle and no holes.
    """
    # Two 10 x 10 rooms joined by a 1-wide, 4-long neck; r = 1 closes the neck.
    corners = [(0, 0), (10, 0), (10, 4), (14, 4), (14, 0), (24, 0), (24, 10), (14, 10), (14, 5), (10, 5), (10, 10), (0, 10)]
    primitives = []
    for (ax, ay), (bx, by) in zip(corners, corners[1:] + corners[:1]):
        dx, dy = bx - ax, by - ay
        length = abs(dx) + abs(dy)
        primitives.append(_line((ax, ay), (bx, by), (-dy // length, dx // length)))
    domain = _held_centre_domain_2.CentreDomain2.build(primitives, _wire(1))
    components = domain.components()
    assert len(components) == 2
    assert all(len(component) == 1 for component in components)
    assert domain.refiled_outer_ccbs == 0
