"""Test support for the post: a hypothetical envelope and a re-validating setup rebuild."""

from dataclasses import fields
from fractions import Fraction

from compas_cgal.post.machine import PostSetup
from compas_cgal.post.machine import TravelEnvelope
from compas_cgal.post.quantities import ProgramMillimetre

# A hypothetical machine travel box in the G54 work frame, mm. Q10 has no
# default: this box exists only so the tests can build a setup; it is generous
# around the 36 x 24 mm pocket placed at the work origin.
TEST_ENVELOPE = TravelEnvelope.build(
    x=(Fraction(-100), Fraction(400)),
    y=(Fraction(-100), Fraction(300)),
    z=(Fraction(-50), Fraction(100)),
)

# The pocket's single forced plunge (FIRST_CUT, into material) is accepted.
ACCEPTED_PLUNGES = 1


def rebuild(setup: PostSetup, **changes: object) -> PostSetup:
    """Re-run `PostSetup.build` on *setup*'s fields with *changes* applied (validation included)."""
    arguments: dict[str, object] = {field.name: getattr(setup, field.name) for field in fields(setup)}
    arguments.update(changes)
    return PostSetup.build(**arguments)  # type: ignore[arg-type]  # field-name kwargs, checked at run time


def mm(value: int | Fraction) -> ProgramMillimetre:
    """An exact work-frame length, mm."""
    return ProgramMillimetre(Fraction(value))
