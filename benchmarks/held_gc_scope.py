"""A scoped pause of Python's cyclic garbage collector around a generator run.

The Held generators allocate about 10^6 short-lived exact objects (``Fraction``,
frozen dataclasses) that reference counting frees; the cyclic collector's
generation-2 passes over that live set take 13-26 % of generation and reclaim
almost nothing (8 objects after a full Monstera run, lever-2 plan). The pause
decides nothing: reference counting still frees every acyclic object, and the
one collection on exit reclaims any cycles created inside the scope.
"""

from __future__ import annotations

import contextlib
import gc
from collections.abc import Iterator


@contextlib.contextmanager
def paused_cyclic_gc() -> Iterator[None]:
    """Disable the cyclic collector for the scope; restore the prior state on every exit.

    Only the scope that found the collector enabled re-enables it, then runs one
    ``gc.collect()``. A nested scope, or a scope entered with the collector
    already disabled, leaves it disabled and does not collect. The restore runs
    on normal exit and on exceptions; exceptions propagate unchanged.

    Yields:
        None.
    """
    was_enabled = gc.isenabled()
    gc.disable()
    try:
        yield
    finally:
        if was_enabled:
            gc.enable()
            gc.collect()
