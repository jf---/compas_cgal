# Task 0.1 gate analysis (2026-10-14)

**Missed row:** "every engagement criterion (cap exceedances, max step, slotting) moves on at least one pocket". Cap exceedances and max step move. Slotting is identical on all 12 cases. Every other row passes.

**Cause:** this is a mis-predicted acceptance row, not a defect. Slotting counts straight motions whose engagement exceeds `SLOT_ENGAGEMENT_FRACTION × cap` (`benchmarks/survey.py:100`, 0.5: 60° at cap 120, 20° at cap 40). The in-process model lowers a straight slot from 360° to 180°. Both values sit far above either threshold, so no straight motion crosses it. Every counted motion engages at least half the cap, and that is physically true. Only a definition change moves the count. The plan puts that change in Task 0.3 rule (c): slotting = a straight motion over the cap, the same rule as loops.

## Evidence (orchestrator re-measurement on the uncommitted 0.1 tree, `build/task01/gate_vectors.py`)

| case | cap exceedances, old → new | max step (deg), old → new | slotting, old → new |
| --- | --- | --- | --- |
| 120 ec rect_12x8 | 5 → 1 | 329.49 → 157.15 | 0 → 0 |
| 120 rr rect_12x8 | 5 → 1 | 329.49 → 157.15 | 0 → 0 |
| 120 ec rect_20x12 | 9 → 9 | 337.54 → 164.20 | 4 → 4 |
| 120 rr rect_20x12 | 9 → 9 | 337.54 → 164.20 | 4 → 4 |
| 120 ec L_shape | 6 → 1 | 325.42 → 152.77 | 0 → 0 |
| 120 rr L_shape | 6 → 1 | 325.42 → 152.77 | 0 → 0 |
| 40 ec rect_12x8 | 145 → 21 | 356.18 → 177.13 | 40 → 40 |
| 40 rr rect_12x8 | 54 → 21 | 349.89 → 169.89 | 41 → 41 |
| 40 ec rect_20x12 | 185 → 105 | 357.19 → 178.02 | 108 → 108 |
| 40 rr rect_20x12 | 125 → 65 | 357.19 → 178.02 | 116 → 116 |
| 40 ec L_shape | 221 → 6 | 355.04 → 176.08 | 15 → 15 |
| 40 rr L_shape | 229 → 105 | 347.52 → 167.52 | 62 → 62 |

Uncut, degenerate, redundant and tangent breaks are identical in all 12 cases, and gouge-free holds. Straight-slot probe: `[0.0, 180.0, 180.0] … [180.0, 180.0, 180.0]`, where it was `[0.0, 192.8, 205.7] … [360.0, …]` before.

## Smallest amendment

Move "slotting moves" from Task 0.1's acceptance to Task 0.3, where rule (c) redefines slotting against the cap. Task 0.1's acceptance becomes: "cap exceedances and max engagement step move on at least one pocket; slotting, uncut, gouge, tangent and continuity stay unchanged."
