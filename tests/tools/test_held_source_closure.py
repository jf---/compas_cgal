"""Contract of tools/held_source_closure.py: the measured code's identity is its whole repository import closure."""

from __future__ import annotations

import tools.held_connector_repair_budget_results  # noqa: F401  (loads both tools' closures)
import tools.held_splice_ab_results  # noqa: F401
from tools.held_source_closure import CLOSURE_ROOTS
from tools.held_source_closure import capture
from tools.held_source_closure import digest_contents
from tools.held_source_closure import source_closure

# Decision-relevant modules both tools import only transitively; the hand-picked SOURCES lists missed them.
TRANSITIVE = ("benchmarks/held_figure5_path.py", "benchmarks/held_figure5_boundary_path.py", "benchmarks/held_reference_cases.py", "benchmarks/spec.py")


def test_the_closure_holds_every_imported_repository_module_sorted_and_repo_relative() -> None:
    closure = source_closure()
    assert list(closure) == sorted(set(closure))
    assert all(path.startswith(CLOSURE_ROOTS) and path.endswith(".py") for path in closure)
    for path in ("tools/held_connector_repair_budget_results.py", "tools/held_splice_ab_results.py", *TRANSITIVE):
        assert path in closure
    assert any(path.startswith("src/compas_cgal/") for path in closure)
    assert not any(path.startswith("tests/") for path in closure)


def test_a_change_in_a_transitively_imported_module_changes_the_digest() -> None:
    contents = capture(source_closure())
    changed = dict(contents)
    changed["benchmarks/held_figure5_path.py"] = contents["benchmarks/held_figure5_path.py"] + b"\n"
    assert digest_contents(contents) == digest_contents(capture(source_closure()))
    assert digest_contents(changed) != digest_contents(contents)


def test_the_digest_names_the_paths_too() -> None:
    contents = {"a.py": b"x", "b.py": b""}
    assert digest_contents(contents) != digest_contents({"a.py": b"", "b.py": b"x"})
    assert digest_contents(contents) != digest_contents({"a.py": b"x"})
