# Measured constants

Every numeric constant in this repository that was *measured* rather than derived
carries its provenance in a comment at its definition — the pocket, the cap, the
tool diameter, and the observed value. Those comments are the authority. Search
for `Measured on` or `MEASURED, NOT DERIVED` to find them.

The generator/benchmark claim-ledger that used to mirror those comments into
signed artifacts was removed on 2026-09-08: it had no consumer outside its own
tests, and its live-scan gate had begun failing legitimate work by asserting that
the repository contained exactly one measurement comment.

## Engagement is not monotone in spacing

This ruling is not recoverable from the source comments, and
`benchmarks/mathsm.py` depends on it.

Over twelve sampled spacings (0.025 … 0.6 tool diameters) on `rect_20x12` with a
2 mm tool, maximum engagement after entry is **non-monotone**: spacing 0.025
measures 131.14°, while the coarser 0.1 measures 98.73°. Tightening the spacing
stops helping and then hurts.

What that establishes, and what it does not:

- It **is** a counterexample to assuming a generic spacing order. A selector may
  not bisect on spacing or stop at the first compliant sample.
- It is **not** a universal spacing law, a common peak station, or a causal
  geometric mechanism. It is one configuration.

`benchmarks.mathsm` therefore measures the complete requested sweep and selects
the shortest compliant trial, because its protocol carries no monotonicity
contract. The same non-monotonicity holds in guide radius; see
`benchmarks/mathsm.py` for the sweep definition.
