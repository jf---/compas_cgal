from __future__ import annotations

import benchmarks.exceedance as exceedance
from benchmarks.exceedance import count_truly_exceeding
from benchmarks.exceedance import exceedance_positions
from benchmarks.families.analytic import rectangle
from benchmarks.runner import generate_toolpath
from benchmarks.spec import PocketSpec
from compas_cgal.engagement import audit_toolpath_engagement

# A denser sampling than the module default, used to show the count has converged
# rather than being an artifact of how finely each motion is walked.
REFINED_SAMPLES_PER_MOTION = 2 * exceedance.EXCEEDANCE_SAMPLES_PER_MOTION

# Caps bracketing the unregulated generator's behaviour: it slots into virgin
# stock, so the tight cap is exceeded and the loose one must be exceeded no more
# often. 180 degrees is the kernel's contractual maximum.
TIGHT_CAP_DEG = 90.0
LOOSE_CAP_DEG = 180.0

# The cap crosses into the exact kernel as the rational surrogate 4*sin^2(cap/2),
# so a run the exact predicate decides is over the cap can be REPORTED a sub-ulp
# below the typed angle. This bounds that documented API gap (docs/exactness.md,
# boundary doctrine). It is not a decision tolerance -- the decision was already
# made, exactly, before this number was produced.
CAP_SURROGATE_SLACK_DEG = 1e-9


def _pocket(tea_cap_deg: float = 120.0) -> PocketSpec:
    """The smallest instance whose unregulated toolpath still exceeds a 120 deg cap."""
    return rectangle(width=8.0, height=6.0, tool_diameter=2.0, tea_cap_deg=tea_cap_deg)


def test_unregulated_generator_is_measurably_over_the_cap() -> None:
    """The unregulated generator slots into virgin stock, so exceedance is demonstrable.

    A zero here would mean the measurement never fires, which no later assertion
    about its value could distinguish from a correct zero.
    """
    spec = _pocket()
    assert count_truly_exceeding(spec, generate_toolpath(spec)) > 0


def test_count_is_stable_under_a_finer_walk(monkeypatch) -> None:
    """Doubling the sample density does not change the count.

    The count is a sampled lower bound, so it can only rise with density. That it
    does not is the evidence behind `EXCEEDANCE_SAMPLES_PER_MOTION`: the reported
    number is the converged one, not a function of the walk.
    """
    spec = _pocket()
    result = generate_toolpath(spec)
    coarse = count_truly_exceeding(spec, result)
    monkeypatch.setattr(exceedance, "EXCEEDANCE_SAMPLES_PER_MOTION", REFINED_SAMPLES_PER_MOTION)
    assert count_truly_exceeding(spec, result) == coarse


def test_a_looser_cap_is_never_exceeded_more_often() -> None:
    """Monotonicity in the cap, which any correct exceedance measure must have."""
    result = generate_toolpath(_pocket())
    tight = count_truly_exceeding(_pocket(TIGHT_CAP_DEG), result)
    loose = count_truly_exceeding(_pocket(LOOSE_CAP_DEG), result)
    assert loose <= tight


def test_exceedance_never_outruns_uncertifiability() -> None:
    """A motion demonstrably over the cap cannot also have been certified under it.

    The soundness contract between the corpus's two cap columns, and the reason
    they are not interchangeable. A sampled position over the TRUE cap is also
    over the certifier's guarded (strictly smaller) cap, so a sound certificate
    must refuse that motion: `truly_exceeding <= uncertified` always. The converse
    fails freely -- `uncertified` counts operations that were never measured at
    all -- which is precisely why reporting it as a violation count misleads.
    """
    spec = _pocket()
    result = generate_toolpath(spec)
    report = audit_toolpath_engagement(spec.polygon, result, spec.tool_diameter, spec.tea_cap_rad, holes=list(spec.holes))
    truly_exceeding = count_truly_exceeding(spec, result)
    assert truly_exceeding <= report.cap_violations


def test_positions_attribute_the_count_to_motions() -> None:
    """Every counted motion is locatable, and no other motion appears in the rows."""
    spec = _pocket()
    result = generate_toolpath(spec)
    rows = exceedance_positions(spec, result)
    assert rows, "no over-cap position was reported for a toolpath known to exceed"
    assert len({index for index, _x, _y, _tea in rows}) == count_truly_exceeding(spec, result)
    for index, _x, _y, tea_deg in rows:
        assert 0 <= index < len(result.operations)
        # A position is reported only when the exact predicate fires, so its
        # reported run must sit at or above the cap it was tested against.
        assert tea_deg >= spec.tea_cap_deg - CAP_SURROGATE_SLACK_DEG
