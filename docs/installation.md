# Installation

## Stable

Stable releases of `compas_cgal` are available from conda-forge through
Pixi.

```bash
pixi init
pixi add -c conda-forge compas_cgal
```

Several examples use the COMPAS Viewer for visualisation.
Add `compas_viewer` to the same workspace when needed:

```bash
pixi add -c conda-forge compas_viewer
```

## Dev Install

Clone the repository, create its locked environment, and run the parallel
test gate:

```bash
git clone https://github.com/compas-dev/compas_cgal.git
cd compas_cgal
pixi install --locked
pixi run pytest -- tests -n auto
```

### Clean native builds

The first editable build downloads and verifies the locked CGAL, Boost, and
Eigen source trees before compiling any target that consumes their headers.
That ordering is part of the build contract, not a manual bootstrap step.
`coverage_exact_core` and every nanobind module therefore depend explicitly on
the `external_downloads` target. A clean worktree must be able to run the test
gate directly; requiring a prior `native-lock-check` would conceal a broken
dependency graph.

Use the standalone lock gate when auditing the native inputs themselves:

```bash
pixi run native-lock-check
```
