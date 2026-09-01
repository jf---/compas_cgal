# Held–Pfeiffer 2025 reference pockets

These four cases reconstruct the green pocket outlines in Figures 5 and 8 of Martin Held and Josef Pfeiffer, “Trochoidal Tool Paths for Pocket Machining with Full Control of the Tool Engagement Angle,” *Computer-Aided Design & Applications* 22(4), 731–747, DOI `10.14733/cadaps.2025.731-747`.

The coordinates come from the publisher PDF vectors. They are normalized so the depicted red tool radius equals 1 mm; the paper does not publish a physical drawing scale. Publisher-page Y is reflected into a CCW world-XY boundary. Figure 7 contributes only a shape-level tool-centre observation for Figure 5 and is neither boundary authority nor a numeric fidelity gate.

Each JSON document retains the selected PDF primitives, the normalization transform, a certified line/arc reconstruction, and its deterministic polygon projection. Regenerate and verify them from the publisher PDF with:

```bash
pixi run held-reference-generate -- path/to/paper.pdf
pixi run held-reference-generate -- --check path/to/paper.pdf
```

The generator writes only the four canonical filenames in this directory. The loader rejects duplicate keys, non-finite values, schema additions, unsupported versions, and evidence that does not reproduce from the stored source primitives.
