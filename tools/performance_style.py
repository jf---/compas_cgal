"""Shared look of the performance heatmaps: palette, surface, typography, cell geometry."""

from __future__ import annotations

from matplotlib.colors import LinearSegmentedColormap

# Sequential blue ramp (dataviz reference palette, steps 100 -> 700).
RAMP = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
SURFACE, INK, MUTED, NA = "#fbfaf8", "#1f1f1d", "#6b6b66", "#efeeea"
CMAP = LinearSegmentedColormap.from_list("seqblue", RAMP)
FIGURE_DPI = 200
# Cells are inset by this fraction of a unit cell on each side.
CELL_INSET = 0.02
CELL_SIDE = 1.0 - 2 * CELL_INSET
# Share above which a cell's label is bold / drawn in white, percent.
BOLD_FROM_PERCENT = 20.0
WHITE_TEXT_FROM_PERCENT = 45.0
PERCENT = 100.0
VALUE_FONT_SIZE = 11
NA_FONT_SIZE = 9
TITLE_FONT_SIZE = 13
FOOTER_FONT_SIZE = 7.5
