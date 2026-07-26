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
