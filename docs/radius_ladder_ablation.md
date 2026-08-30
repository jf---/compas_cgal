# Radius-Ladder Ablation Record

*neither half of the radius-ladder repair is safe alone; the gate exists
because both together UNGATED is worse than ranking alone.*

This page preserves a four-way development ablation as historical evidence.
The configuration labels describe internal variants used for that comparison;
they are not recipes or supported product settings.

!!! warning "Provenance"

    These four configurations were measured during development (2026-08-22)
    by toggling internals that shipped as a single gated design; they are **not
    reconstructible from shipped knobs** and are recorded as historical
    evidence. The one shipped knob, `RADIUS_LADDER_SUBDIVISIONS`, has a
    reproducible sweep committed beside it
    (`src/compas_cgal/engagement_radial_toolpath.py`, commit `29050b0`).

    The `20x12` baseline row below records 8 over-cap circles. A distinct,
    removed source assertion authenticated during Task 6 records 12. The
    available configuration descriptions do not establish that both assertions
    describe the same run. Both therefore remain separate historical assertions;
    neither is promoted as reconstructed truth.

## Recorded four-way comparison

Each pocket reports three measurements in the same order: worst tool-engagement
angle (TEA), number of over-cap machining circles, and cutting length. The
requested cap was 40 degrees for the `6x4` pocket and 60 degrees for the `20x12`
pocket.

| Configuration | `6x4`: worst TEA (deg) | `6x4`: over-cap circles | `6x4`: cutting length | `20x12`: worst TEA (deg) | `20x12`: over-cap circles | `20x12`: cutting length |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Baseline | 86.4 | 34 | 295 | 126.1 | 8 | 3236 |
| Ranking only | 54.5 | 41 | 389 | 88.6 | 12 | 3382 |
| Refinement only | **213.6** | 1 | 357 | 126.1 | 4 | 3388 |
| Both, ungated | 121.0 | 11 | 830 | 90.6 | 12 | 3379 |

The halves pull in different directions. Ranking alone lowers the recorded peak
but increases the over-cap count on both pockets. Refinement alone lowers that
count, yet the `6x4` run contains the recorded **213.6-degree slotting cut**.
Running both halves without the gate also gives a worse peak than ranking alone
on both pockets. The gate exists to keep those two mechanisms coupled only in
the measured regime selected by the shipped design.

## Forced-station constraint

In the recorded design, at a forced station every candidate is over-cap by
definition and a sub-maximal circle does not finish the station, so lowering the
peak always costs a circle. A one-step lookahead verified that result at 32/32
stations in the recorded configuration. This is a structural opposition,
not a tuning trade-off.

The 32-of-32 result is a verified historical constraint for that configuration,
not a universal mathematical proof for other geometries, caps, or generator
revisions.
