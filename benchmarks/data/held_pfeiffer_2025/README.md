# Held–Pfeiffer 2025 reference pockets

The four pocket cases reconstruct the green outlines in Figures 5 and 8 of Martin Held and Josef Pfeiffer, “Trochoidal Tool Paths for Pocket Machining with Full Control of the Tool Engagement Angle,” *Computer-Aided Design & Applications* 22(4), 731–747, DOI `10.14733/cadaps.2025.731-747`. A fifth artifact, `figure5_publisher_turns.json`, records Figure 5(a)'s ordered publisher stream: fitted turn geometry and diagnostics, stream endpoints, and every connector primitive in reconstructed world XY.

The coordinates come from the publisher PDF vectors. They are normalized so the depicted red tool radius equals 1 mm; the paper does not publish a physical drawing scale. Publisher-page Y is reflected into a CCW world-XY boundary. Figure 7 contributes only a shape-level tool-centre observation for Figure 5 and is neither boundary authority nor a numeric fidelity gate.

Check the fifth artifact separately against an official page-12 SVG with `pixi run held-figure5-publisher-check -- path/to/publisher-page-12.svg`.

Each pocket-case JSON retains the selected PDF primitives, the normalization transform, a certified line/arc reconstruction, and its deterministic polygon projection. Regenerate and verify the four pocket cases from the publisher PDF with:

```bash
pixi run held-reference-generate -- path/to/paper.pdf
pixi run held-reference-generate -- --check path/to/paper.pdf
```

The pocket-case generator writes only its four canonical filenames. Their loader rejects duplicate keys, non-finite values, schema additions, unsupported versions, and evidence that does not reproduce from stored source primitives. The separate Figure 5 source-check CLI re-extracts the official SVG and compares every tracked publisher-stream field.
