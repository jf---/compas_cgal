# Rejected generator measurement bundle

> **status: rejected** — preserved byte-for-byte outside the accepted-result
> root on 2026-08-29.

## Identity

- Input identity SHA-256:
  `8a73a14731d21986efba86c8354245dc960cf6d6f0279d8ff706fa300d3eec09`
- Result identity SHA-256:
  `38a0775a821532a33196a8e9429551c07ec9846bf6fe933a6b81364e7554a446`
- Payload file SHA-256 (`generator-claims.json`):
  `e1d8eeb6386f159dee98deed7302531dfafc1b0cec31f2ebecd101acf6881826`
- Stamp file SHA-256 (`stamp.json`):
  `da16c47308aefe8172381c83af240457dba22849af92c97b8c49b717cc3edffc`

The payload is 107,494 bytes and the stamp is 30,259 bytes. Rehashing both
files after the move produced the same file hashes recorded above.

## Semantic rejection

This bundle declares `artifact_kind: generator-measurement-claims/v1`. The
sole current ledger consumer requires the canonical v2 envelope and rejects
this bundle with `artifact_kind does not match consumer`; v1 does not provide
the required distinct execution and source-correction identities. Its numbers
remain evidence worth preserving, but the bundle cannot satisfy Task-6 ledger
acceptance and must never reside under
`benchmarks/measurement_claim_results/`.
