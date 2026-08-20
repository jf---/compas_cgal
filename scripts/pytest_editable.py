"""Run pytest with scikit-build-core's editable rebuild suppressed for this project.

scikit-build-core's editable import hook re-runs CMake whenever a compiled extension is
imported. Under ``pytest -n auto`` every xdist worker imports the extensions, so every worker
would launch its own concurrent CMake run. Setting ``SKBUILD_EDITABLE_SKIP`` to the finder's
build directory suppresses that; the ``_editable-rebuild`` pixi task performs the single real
rebuild beforehand, so the extensions are already current.

The marker is an ``os.pathsep``-joined list of build directories, not a boolean, so the exact
path matters — see ``scikit_build_core/_editable_redirect.py`` (``rebuild(path=self.path)``
and the ``if path in env.get(MARKER, "").split(os.pathsep)`` guard).

This runs as a script rather than inline shell because pixi's task shell (deno_task_shell)
does not propagate a command substitution's exit status to the enclosing assignment: in
``build_dir="$(failing-command)" && pytest`` the ``&&`` does NOT short-circuit, so pytest
would run with an empty marker and silently fall back to a rebuild storm. Resolving the path
in Python instead makes a failed lookup abort loudly, before pytest starts.
"""

from __future__ import annotations

import os
import sys

MARKER = "SKBUILD_EDITABLE_SKIP"
PROBE_MODULE = "compas_cgal._stock_2"
"""Compiled extension used to identify this project's editable finder on ``sys.meta_path``."""


class EditableBuildDirNotFound(RuntimeError):
    """No scikit-build-core editable finder on ``sys.meta_path`` provides the probe module."""


def find_editable_build_dir(module: str = PROBE_MODULE) -> str:
    """Return the build directory of the editable finder that provides ``module``.

    Args:
        module: Fully qualified name of a compiled extension built by this project.

    Returns:
        The build directory that scikit-build-core rebuilds into, exactly as it appears in
        the ``SKBUILD_EDITABLE_SKIP`` marker.

    Raises:
        EditableBuildDirNotFound: If no editable finder on ``sys.meta_path`` claims ``module``.
    """
    for finder in sys.meta_path:
        known_wheel_files = getattr(finder, "known_wheel_files", None)
        if known_wheel_files is not None and module in known_wheel_files:
            return finder.path
    raise EditableBuildDirNotFound(
        f"no scikit-build-core editable finder on sys.meta_path provides {module!r}, "
        f"so {MARKER} cannot be set and every pytest-xdist worker would re-run CMake "
        "concurrently. Run `pixi install` to recreate the editable install; if the "
        f"extension was renamed, update PROBE_MODULE in {__file__}."
    )


def main(argv: list[str]) -> None:
    """Replace this process with pytest, with the editable-skip marker set.

    Args:
        argv: Arguments to forward to pytest.

    Raises:
        EditableBuildDirNotFound: If the editable build directory cannot be resolved, in
            which case pytest is never started.
    """
    os.environ[MARKER] = find_editable_build_dir()
    os.execv(sys.executable, [sys.executable, "-m", "pytest", *argv])


if __name__ == "__main__":
    main(sys.argv[1:])
