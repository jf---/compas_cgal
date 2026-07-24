# CGAL 6.0.1 adaptive-package license gate

Status: `CGAL_LICENSE_GATE_ERROR`

The repository is distributed under LGPL-3.0, as recorded in `LICENSE`,
`docs/license.md`, and `pyproject.toml`. No repository evidence records either
a GPL-compatible distribution decision for the adaptive extension or a
commercial CGAL entitlement covering CGAL 6.0.1. The README's general
compatibility statement does not resolve the package-specific GPL terms below.

The authoritative package classifications are the
[CGAL 6.0.1 package overview](https://doc.cgal.org/6.0.1/Manual/packages.html)
and the SPDX declarations in the locked CGAL 6.0.1 headers. The exact archive
and extracted-source digest are fixed in
`docs/build/native-dependency-lock.json`.

| Package | Planned headers or types | Tasks | CGAL 6.0.1 license | Decision |
|---|---|---:|---|---|
| [2D Arrangements](https://doc.cgal.org/6.0.1/Arrangement_on_surface_2/index.html) | `Arr_circle_segment_traits_2.h`, `Arr_algebraic_segment_traits_2.h` | Existing `_stock_2`; 3–8 | GPL-3.0-or-later or commercial | Blocked |
| [Algebraic Kernel](https://doc.cgal.org/6.0.1/Algebraic_kernel_d/index.html) with bundled CORE | `Algebraic_kernel_d_1.h`, `Algebraic_kernel_d_2.h`, `CORE::BigInt` | 5–10 | LGPL-3.0-or-later | Compatible in isolation |
| [2D Regularized Boolean Set Operations](https://doc.cgal.org/6.0.1/Boolean_set_operations_2/index.html) | `Boolean_set_operations_2.h`, `General_polygon_set_2.h`, `Gps_circle_segment_traits_2.h` | Existing `_stock_2`; 3, 7, 11 | GPL-3.0-or-later or commercial | Blocked |
| [2D Segment Delaunay Graphs](https://doc.cgal.org/6.0.1/Segment_Delaunay_graph_2/index.html) | `Segment_Delaunay_graph_2.h`, `Segment_Delaunay_graph_traits_2.h`, `Segment_Delaunay_graph_storage_traits_with_info_2.h` | 9–10 | GPL-3.0-or-later or commercial | Blocked |
| [2D Voronoi Diagram Adaptor](https://doc.cgal.org/6.0.1/Voronoi_diagram_2/index.html) | `Voronoi_diagram_2.h`, `Segment_Delaunay_graph_adaptation_traits_2.h`, `Segment_Delaunay_graph_adaptation_policies_2.h` | 9–10 | GPL-3.0-or-later or commercial | Blocked |
| [2D Apollonius Graphs](https://doc.cgal.org/6.0.1/Apollonius_graph_2/index.html) | `Parabola_segment_2.h` | 9–10 | GPL-3.0-or-later or commercial | Blocked |

## Morphology selection

Phase 1 selects the exact line/circle-arrangement construction already defined
by Task 3: offset-edge lines and vertex-radius circles are arranged, cells are
selected by exact disk containment, and regularized Boolean operations form
`C_r`, `M_r`, and `D \ M_r`. The
[2D Minkowski Sums](https://doc.cgal.org/6.0.1/Minkowski_sum_2/index.html)
package is therefore not selected and its headers must not be instantiated.
Selecting it later would add another GPL-3.0-or-later-or-commercial gate.

## Gate effect

The existing `_stock_2` target already instantiates both 2D Arrangements and
2D Regularized Boolean Set Operations through `stock_2.h`. Tasks 1 and 2
modify and rebuild that target; Task 1A consumes Task 1's types and identities.
Consequently, Tasks 1–16 are blocked until the repository records one of:

1. a GPL-compatible distribution decision applying to the adaptive extension;
2. a commercial CGAL entitlement covering CGAL 6.0.1 and every instantiated
   package.

No sampled geometry, different kernel, package substitution, or silent
fallback satisfies this gate.
