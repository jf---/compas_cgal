# Refutation Soundness Cross-Check

The exact event partition in `continuous_tea_2` is **incomplete**, not unsound,
on the lateral link segments examined here. It never certifies a motion the
station refutation probe refutes. On those motions it either agrees
(`cap_exceeded`) or fails to build a partition at all
(`IncompleteSegmentPartitionError`, and on one family a native process abort).
The probe's contribution is therefore **completeness plus cost**: it decides
motions the partition cannot decide, and it decides the rest about `160x`
cheaper.

Measured over `125` full audits dispatched on the branch family with the probe
disabled (the sweep was stopped there; the distribution had been stable for
tens of audits):

| probe | full audit | count |
| --- | --- | ---: |
| silent | `certified` | `86` |
| `REFUTED` | `cap_exceeded` | `24` |
| `REFUTED` | raised `IncompleteSegmentPartitionError` | `15` |
| **`REFUTED`** | **`certified`** | **`0`** |

This page is the evidence, and it carries two corrections rather than hiding
them. A first draft compared the wrong `engagement_at` field. A second draft
then claimed the partition was *certifying* violating motions; that claim was
inferred from a stale recorded test expectation and is **retracted** — see
[Retracted: the partition was never shown to certify a violation](#retracted-the-partition-was-never-shown-to-certify-a-violation).

## What was measured, and against what

`_stock_2.engagement_at(stock, cx, cy, tool_radius, cap_chord_ratio, gap_close_ratio)`
returns `(total_tea, max_run_tea, cap_exceeded)`.

!!! danger "The cap is a per-run bound, not a total"

    `engagement_2.cpp` decides the cap on each maximal engaged run. Two
    disjoint runs of `1.70 rad` give `total_tea = 3.40` while neither run
    exceeds `pi`, so comparing `total_tea` against `pi` proves nothing at all.

    The first draft of this finding compared `total_tea`. That comparison was
    invalid as written. It is retracted and replaced by the table below.
    The reported numbers did not move, because every witness station in the
    corpus carries exactly one engaged run and `total_tea == max_run_tea` in
    all `31` rows — but that coincidence was not knowable in advance and is not
    what the argument may rest on.

The authoritative field is the third. `engagement_2.h` states that
`cap_exceeded` is "the only DECISION carried out here and it is decided EXACTLY
on the exact arrangement, never from these reported doubles", while `total_tea`
and `max_run_tea` are reporting doubles. Measurement settings:

| parameter | value | why |
| --- | --- | --- |
| `cap_chord_ratio` | `4.0` | exact rational surrogate `4*sin^2(pi/2)` for cap `= pi`, so `cap_exceeded` is directly comparable to the probe's own cap |
| `gap_close_ratio` | `0.0` | no gap-closure pessimism; pessimistic runs equal true runs and every field is the pre-pessimism value bit-for-bit |
| `certify_segment_tea` cap | `math.pi` | branch-A's independent whole-motion certifier, exact station predicates plus the analytic growth guard |

`engagement_2.cpp` and `continuous_tea_2/` share no geometry code. The station
probe decides through `classify_station_cell`; `engagement_at` decides through
`run_exceeds_cap` and `sign_mixed_radical`; `certify_segment_tea` decides by
adaptive station sampling under a guarded cap. Three implementations, one
question.

## Task 13F route 2, operation 7

All `10` links that survive gouge containment. `pi = 3.141593`,
`2*pi = 6.283185`.

| station | exact binary64 point | `total_tea` | `max_run_tea` | `cap_exceeded` | `certify_segment_tea` `max_tea` | `cap_certified` | stations |
| --- | --- | ---: | ---: | --- | ---: | --- | ---: |
| `1/1` | yes | `3.285247` | `3.285247` | `True` | `3.285247` | `False` | `160` |
| `7/8` | no | `3.304183` | `3.304183` | `True` | `3.457686` | `False` | `140` |
| `1/1` | yes | `3.285247` | `3.285247` | `True` | `3.285247` | `False` | `160` |
| `1/1` | yes | `3.290448` | `3.290448` | `True` | `3.290448` | `False` | `159` |
| `1/1` | yes | `3.290448` | `3.290448` | `True` | `3.290448` | `False` | `159` |
| `7/8` | no | `3.306509` | `3.306509` | `True` | `3.459960` | `False` | `139` |
| `1/1` | yes | `3.437048` | `3.437048` | `True` | `3.437048` | `False` | `141` |
| `1/1` | yes | `3.439415` | `3.439415` | `True` | `3.439415` | `False` | `142` |
| `7/8` | no | `3.304183` | `3.304183` | `True` | `3.457686` | `False` | `140` |
| `7/8` | no | `3.306509` | `3.306509` | `True` | `3.459960` | `False` | `139` |

## Branch family, operation 3

All `21` refutations raised while advancing the active family. `certify_segment_tea`
reports a larger `max_tea` than the witness station on the `7/8` rows because it
samples `255` stations across the whole motion and finds a worse one elsewhere.

| station | exact binary64 point | `total_tea` | `max_run_tea` | `cap_exceeded` | `certify_segment_tea` `max_tea` | `cap_certified` | stations |
| --- | --- | ---: | ---: | --- | ---: | --- | ---: |
| `1/1` | yes | `5.150175` | `5.150175` | `True` | `5.150175` | `False` | `250` |
| `1/1` | yes | `3.873024` | `3.873024` | `True` | `3.873024` | `False` | `266` |
| `7/8` | no | `3.341484` | `3.341484` | `True` | `4.613387` | `False` | `255` |
| `1/1` | yes | `5.072280` | `5.072280` | `True` | `5.072280` | `False` | `252` |
| `1/1` | yes | `3.767857` | `3.767857` | `True` | `3.767857` | `False` | `264` |
| `7/8` | no | `3.341484` | `3.341484` | `True` | `4.613387` | `False` | `255` |
| `1/1` | yes | `5.002183` | `5.002183` | `True` | `5.002183` | `False` | `251` |
| `1/1` | yes | `3.665474` | `3.665474` | `True` | `3.665474` | `False` | `267` |
| `7/8` | no | `3.341484` | `3.341484` | `True` | `4.613387` | `False` | `255` |
| `1/1` | yes | `4.938825` | `4.938825` | `True` | `4.938825` | `False` | `250` |
| `1/1` | yes | `3.565592` | `3.565592` | `True` | `3.565592` | `False` | `269` |
| `7/8` | no | `3.341484` | `3.341484` | `True` | `4.613387` | `False` | `255` |
| `1/1` | yes | `4.881447` | `4.881447` | `True` | `4.881447` | `False` | `255` |
| `1/1` | yes | `3.467988` | `3.467988` | `True` | `3.467988` | `False` | `271` |
| `7/8` | no | `3.341484` | `3.341484` | `True` | `4.613387` | `False` | `255` |
| `1/1` | yes | `4.829485` | `4.829485` | `True` | `4.829485` | `False` | `251` |
| `1/1` | yes | `3.372486` | `3.372486` | `True` | `3.372486` | `False` | `271` |
| `7/8` | no | `3.341484` | `3.341484` | `True` | `4.613387` | `False` | `255` |
| `1/1` | yes | `4.782513` | `4.782513` | `True` | `4.782513` | `False` | `254` |
| `1/1` | yes | `3.278950` | `3.278950` | `True` | `3.278950` | `False` | `274` |
| `7/8` | no | `3.341484` | `3.341484` | `True` | `4.613387` | `False` | `255` |

## Totals

| | count |
| --- | ---: |
| refutations examined | `31` |
| `cap_exceeded` `True` | `31` |
| `cap_exceeded` `False` | `0` |
| `cap_certified` `True` (whole motion) | `0` |
| `max_run_tea > pi` | `31` |
| `max_run_tea <= pi` | `0` |
| `max_run_tea` at `2*pi` (`MATERIAL`) | `0` |

`max_run_tea` spans `3.278950` to `5.150175 rad`: the largest single engaged run
exceeds the `pi` cap by `8` to `115` degrees. None is marginal.

## The probe's notion of violation is the sampled oracle's, not an approximation

The two fixture tables above are real geometry but a small corpus. A separate
synthetic differential, pinned as
`test_no_refutation_contradicts_the_independent_exact_cap_flag`, compares the
station probe against `engagement_at`'s exact `cap_exceeded` flag over a `6x4`
pocket, five disk configurations, six segments, caps `pi` and `120 degrees`, two
tool radii, and the whole nine-station ladder. Only stations landing on an exact
binary64 point are compared, so both implementations are asked about the
identical point.

| | count |
| --- | ---: |
| comparisons | `894` |
| probe `REFUTED` | `543` |
| `cap_exceeded` `True` | `579` |
| **probe `REFUTED` while `cap_exceeded` `False`** | **`0`** |
| `cap_exceeded` `True` while probe silent | `36` |

The relation is one-way by design and measured to be one-way in fact. A
refutation is always supported by the independent flag; the flag fires `36`
times where the probe stays silent, which costs only a full audit.

### The start station is the one that could flip

All `36` gaps sit at station `0/1`, and that is not a coincidence:
`segment_oracle.cpp` carries a dedicated `start_disk_has_no_material_interior`
predicate because at `t = 0` the cutter is already at the start. At `0/1` the
corpus records `68` comparisons, `28` probe refutations and `32` flags — the
same `4`-case safe gap an independently written differential found.

!!! warning "Polarity enforced, not assumed"

    Today's polarity is safe: the probe declines to refute where the flag
    fires. The opposite polarity at the same station would be an **unsound
    refutation**, and the station ladder includes `0/1`.
    `test_start_station_refutation_never_contradicts_the_exact_cap_flag` pins
    the implication and asserts the corpus actually refutes at `0/1`, so the
    property is enforced rather than incidental.

## `MATERIAL` versus `CAP_EXCEEDED` (reporting-only)

`segment_station_cap_exceeded_exact` returns true for a station cell classified
`MATERIAL` **or** `CAP_EXCEEDED`. `MATERIAL` means the whole rim is buried, a
full `2*pi`, which violates any cap `<= pi` — so either classification is a
sound refutation. The native predicate returns only a bool, so the split is not
available exactly at the Python boundary.

As a reporting-only discriminator, no row has `max_run_tea` at `2*pi`: the
largest is `5.150175` against `6.283185`. Independently, an earlier run of the
same corpus reported `stock.contains(cx, cy)` `False` at all `10` route-2
witness stations, which rules out a buried rim there. On this evidence all `31`
refutations are `CAP_EXCEEDED`, none `MATERIAL`. If the exact split is wanted as
a hard fact rather than an inference, it needs a native binding returning the
`StationCellDecision` rather than a bool.

## What this does and does not establish

**Established.** At `31` exact rational stations spanning two independent
fixtures, three independent implementations agree the cap is violated. The
refutation probe has produced no false positive.

**Established.** The event partition is incomplete on this family. Of `39` refuted
motions dispatched, `24` produce a matching `cap_exceeded` verdict and `15`
raise `IncompleteSegmentPartitionError`, so the partition cannot decide roughly
a third of them. It is unhealthy in a second way on route 2: a direct
`audit_segment_tea_event_exact` on the first refuted link terminates the process
with `boost float_next<double>: Argument must be finite, but got inf`, and a
probe-disabled continuation run reaches `6.3 GB` RSS before aborting. Because
that abort kills the process, the remaining nine route-2 links are unreachable
and their partition verdicts are unknown.

**Not established — and specifically retracted.** That the partition ever
*certified* a violating motion. See the retraction below.

**Not established.** The root cause inside exact event discovery. This page
shows where the partition fails to decide; it does not show which event class is
missed. That work is open.

**Not established.** Any claim about motions outside these two fixtures. `31`
refutations is the whole corpus the probe has produced so far, not a sample from
a larger validated population.

## Retracted: the partition was never shown to certify a violation

An earlier draft of this page and of the commit that introduced it stated that
`MotionCertifier.certify` returned a `MotionWitness` for these motions, and
therefore that the certifier was **unsound**. That is withdrawn. It was never
measured.

The bad inference: `test_task13f_full_continuation` records an expected
`attempts=56; cap=0; gouge=56`. Reading that as "these ten links were counted as
gouges" implies the link certified and the following circle failed containment.
The recorded expectation is stale — the current tree aborts the process at that
point, so the counts cannot be reproduced and cannot support any inference about
what `certify` returned.

Directly measured instead, `125` full audits with the probe disabled: `0` cases
of `REFUTED` with a `certified` verdict. The partition agrees or gives up.

The distinction is the whole severity of the finding:

| | consequence |
| --- | --- |
| what was claimed — silently certifies violations | unsound; a bad toolpath could be emitted |
| what is measured — cannot decide violations | incomplete and fragile, but **fails safe** |

The second is still worth fixing: an exact certifier that cannot decide motions
which three independent implementations find cleanly over cap by `8` to `115`
degrees, and whose failure mode on one family is a process abort, is not
finished. But it did not emit a wrong certificate, and this page previously said
it did.

## Consequences already visible

Two `tests/adaptive/test_generator.py` cases changed. Neither was modified, and
the analysis lives alongside them in
[What the probe exposed in the full partition](continuous_engagement.md#what-the-probe-exposed-in-the-full-partition).

- `test_task13f_full_continuation` was already failing before the probe, as a
  `worker crashed` native abort. With the probe the run completes and reports
  `attempts=56; cap=10; gouge=46`; its recorded `cap=0; gouge=56` predates the
  crash and is not reproducible on the current tree.
- `test_real_active_family_stops_at_unresolved_exact_event` newly fails. Its
  guarded proof gap is now a proof of violation — the single run of
  `5.150175 rad` in the first branch-family row above — so the search skips that
  candidate soundly and accepts a later one. This changes the accepted toolpath
  for that fixture. It is sound, and it is not observationally equivalent.

Both changes come from the motions in the incomplete bucket. Where the
partition decides, the probe agrees with it and the only difference is cost;
where the partition gives up, the probe supplies a verdict the partition could
not. No certificate was ever wrong, and no incorrect toolpath was emitted.
