"""Assign every ledger file to exactly one census slice; fail on any file left unassigned or assigned twice."""

from __future__ import annotations

import re
import sys
from pathlib import Path

HERE = Path(__file__).parent
LEDGER = HERE / "ledger.tsv"

# Slice id -> ordered list of regex rules (first matching slice wins, in SLICE_ORDER).
RULES: dict[str, list[str]] = {
    "g1_python_generators_and_guide": [
        r"^src/compas_cgal/[^/]+\.pyi?$",
        r"^src/toolpath\.(cpp|h)$",
        r"^src/straight_skeleton_2\.(cpp|h)$",
        r"^src/medial_axis_2\.(cpp|h)$",
        r"^src/boundary_normal_circle_2\.(cpp|h)$",
        r"^src/circle_geometry_2.*$",
    ],
    "g2_adaptive": [r"^src/compas_cgal/adaptive/"],
    "g3_segment_site_mat": [r"^src/segment_site"],
    "g5_tea_and_exact": [r"^src/continuous_tea_2", r"^src/exact/", r"^src/exact_algebraic"],
    "g4_cpp_kernels": [r"^src/[^/]+\.(cpp|h|hpp)$", r"^src/[^/]+/"],
    "g6_held_generators": [r"^benchmarks/held_", r"^benchmarks/tools/"],
    "d2_plans": [r"^docs/superpowers/plans/"],
    "d1_docs_specs_reviews": [r"^docs/", r"^[^/]+\.md$", r"^\.superpowers/"],
    "v_validation_tests_tools_build": [r".*"],
}
SLICE_ORDER = list(RULES)

ledger = [line.split("\t") for line in LEDGER.read_text().splitlines()]
assignment: dict[str, list[tuple[int, str]]] = {name: [] for name in SLICE_ORDER}
for count, path in ledger:
    for name in SLICE_ORDER:
        if any(re.search(rule, path) for rule in RULES[name]):
            assignment[name].append((int(count), path))
            break
    else:
        sys.exit(f"unassigned: {path}")

assigned = [path for rows in assignment.values() for _, path in rows]
assert len(assigned) == len(set(assigned)) == len(ledger), "partition is not a bijection onto the ledger"

slices = HERE / "slices"
slices.mkdir(exist_ok=True)
for name, rows in assignment.items():
    (slices / f"{name}.tsv").write_text("".join(f"{count}\t{path}\n" for count, path in rows))
    print(f"{sum(c for c, _ in rows):7d} lines {len(rows):4d} files  {name}")
