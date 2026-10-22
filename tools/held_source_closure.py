"""Content identity of measured code: the repository import closure of a measuring process, hashed.

A record's ``source_sha256`` must change whenever any code that decides its
numbers changes. A hand-picked file list cannot promise that (a transitively
imported module goes unlisted), so the identity is taken from ``sys.modules``:
every loaded Python module whose file lies under ``CLOSURE_ROOTS`` of the
repository, repo-relative and sorted. Native extensions are not hashed here;
their identity is the build identity's.
"""

from __future__ import annotations

import hashlib
import pathlib
import sys
from collections.abc import Iterable
from collections.abc import Mapping

from benchmarks.performance_record import REPOSITORY

# Repository subtrees whose modules can decide a measured number (repo-relative, POSIX).
CLOSURE_ROOTS = ("benchmarks/", "tools/", "src/compas_cgal/")
PYTHON_SUFFIX = ".py"
# Width of the big-endian length prefix hashed before each path and each file's bytes (lengths < 2**64).
LENGTH_PREFIX_BYTES = 8


def source_closure() -> tuple[str, ...]:
    """Repo-relative paths of every loaded Python module under ``CLOSURE_ROOTS``, sorted.

    Call it after the measured code has run: modules imported lazily during
    the run are only in the closure once they are loaded.
    """
    paths = set()
    for module in list(sys.modules.values()):
        filename = getattr(module, "__file__", None)
        if filename is None:
            continue
        path = pathlib.Path(filename).resolve()
        if path.suffix != PYTHON_SUFFIX or not path.is_relative_to(REPOSITORY):
            continue
        relative = path.relative_to(REPOSITORY).as_posix()
        if relative.startswith(CLOSURE_ROOTS):
            paths.add(relative)
    return tuple(sorted(paths))


def capture(paths: Iterable[str]) -> dict[str, bytes]:
    """The current bytes of each repo-relative path."""
    return {path: (REPOSITORY / path).read_bytes() for path in paths}


def digest_contents(contents: Mapping[str, bytes]) -> str:
    """SHA-256 over length-prefixed (path, bytes) pairs in sorted path order: changes with any path or any byte."""
    digest = hashlib.sha256()
    for path in sorted(contents):
        for field in (path.encode("utf-8"), contents[path]):
            digest.update(len(field).to_bytes(LENGTH_PREFIX_BYTES, "big"))
            digest.update(field)
    return digest.hexdigest()
