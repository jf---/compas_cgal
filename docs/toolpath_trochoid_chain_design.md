# Trochoid Chain: Circle Stitching at Varying Radii

When generating a trochoidal toolpath along a MAT edge, the trochoid radius varies
with distance to the boundary.  Adjacent circles in the chain have different radii,
which means the external tangent arrival point on circle *i* does not coincide with
the tangent departure point for the next pair.  How these are stitched together
affects both geometric correctness and output cleanliness.

## The geometry

Each circle in the chain connects to the next via an external tangent line.
The tangent touches circle *i* at a **departure point** and circle *i+1* at
an **arrival point**.  When all radii are equal, the arrival from the previous
tangent and the departure for the next tangent land on the same point — no gap.

When radii vary, these two points separate by a small angular offset:

```
            arrival from prev tangent
                 ↓
                 A·────────────── tangent to circle i+1
                / \
               /   \
     circle i /  C  \   ← center
               \   /
                \ /
                 D·────────────── tangent from circle i-1
                 ↑
            departure for next tangent
```

The arc from **A** (arrival) to **D** (departure) is the gap that must be covered to
keep the chain continuous.

## Three approaches

### 1. Bridge arc + full circle (overlapping)

```cpp
// Bridge: short arc A → D
chain.push_back(make_arc(center, r, A, D, cw));
// Full circle at D
chain.push_back(make_circle(center, r, D, ori));
```

Correct for CNC (full material removal), but produces **two overlapping primitives**:
the bridge covers an angular range that the full circle also covers.  This causes
visual duplication in viewers (a smooth circle + faceted arc on top of each other)
and inflates toolpath length calculations.

### 2. Single long-way arc (crude, no overlap)

```cpp
// One arc: A → D taking the long way around (~355°)
chain.push_back(make_arc(center, r, A, D, cw));
```

Eliminates overlap but misses a small wedge (~5°) of material.  For dense pitch
relative to tool radius, adjacent circle overlap covers the wedge in practice.
But geometrically imprecise — the output does not describe a full revolution.

### 3. Fixed-winding long-way arc (broken)

```cpp
// Always use nominal winding direction
chain.push_back(make_arc(center, r, A, D, cw));
```

Fails because the gap between A and D can fall on **either side** of the circle
depending on the radius gradient direction.  When the gap is on the winding side,
`sweep()` returns the short way (~5-90°) instead of the long way (~355°).
This produces visible gaps and uneven spacing in the toolpath.

### 4. Adaptive-winding long-way arc (current implementation)

```cpp
if (has_prev && prev_arrival != tangent.source()) {
    // Try nominal winding; flip if it takes the short way
    auto arc = make_arc(center, r, prev_arrival, tangent.source(), cw);
    if (arc.sweep() < pi) {
        arc = make_arc(center, r, prev_arrival, tangent.source(), !cw);
    }
    chain.push_back(arc);
} else {
    // First circle or constant radius: full 360°.
    chain.push_back(make_circle(center, r,
        has_prev ? prev_arrival : tangent.source(), ori));
}
```

Always takes the long way (sweep ≥ π, typically ~355°).  Uses exact `Point_2`
equality (free with Epick) for the constant-radius case.  One primitive per
circle, no overlap.

The occasional winding flip affects only the direction of the last ~5° of the arc.
For CNC this is negligible — the tool still makes a near-full revolution at each
station and cutting forces are dominated by the other ~355°.

## Why exact equality works

With `Exact_predicates_inexact_constructions_kernel`, point comparison is exact.
When all radii are truly equal, the tangent geometry is symmetric and `prev_arrival`
and `tangent.source()` are constructed from the same arithmetic — yielding identical
`Point_2` values.  When radii differ even slightly, the points separate and the
varying-radius path is taken.  No epsilon threshold is needed.

## Impact on downstream

| Consumer | Bridge + circle | Fixed winding | Adaptive winding |
|----------|----------------|---------------|------------------|
| G-code post-processor | Two moves, small overlap | One move, sometimes short | One move, always long |
| Toolpath length | Over-counts overlap | Under-counts when short | ≤1.5% short for varying radii |
| Viewer rendering | Doubled geometry | Gaps visible | Clean |
| Material coverage | Full (with overlap) | Misses wedges | ≥π at every station |
