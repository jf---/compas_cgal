"""The Gps CCB repair errors are named Python exceptions in every module that runs the repair.

`src/gps_ccb_repair_2.h` raises three corrupt-arrangement errors. Each module
registers its own instantiation, so a throw from one module can never surface
as the other module's class.
"""

import pytest

from compas_cgal import _held_centre_domain_2
from compas_cgal import _stock_2

REPAIR_ERROR_NAMES = (
    "AmbiguousOuterCcbError",
    "BoundedFaceWithoutOuterCcbError",
    "UnboundedFaceOuterCcbError",
)


@pytest.mark.parametrize("name", REPAIR_ERROR_NAMES)
def test_repair_errors_are_named_runtime_errors_per_module(name):
    stock_error = getattr(_stock_2, name)
    domain_error = getattr(_held_centre_domain_2, name)
    assert issubclass(stock_error, RuntimeError)
    assert issubclass(domain_error, RuntimeError)
    assert stock_error is not domain_error
    assert stock_error.__module__ == _stock_2.__name__
    assert domain_error.__module__ == _held_centre_domain_2.__name__
