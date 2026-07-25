# Segment event substrate design

## Scope

Task 5 will expose a source-bound, exact event partition for one nonzero
segment motion. It consumes the exact Task 4 `Stock2` boundary records and
supplies the complete event strata required by Task 6. Task 6 continues to own
material seeding at `t = 0`, cyclic-run assembly, verdict selection, and
`EventTrace`.

Full-circle motion is a separate follow-on slice.

## Exact source boundary

`SegmentEventSource2::from_binary64` accepts public binary64
`x0`, `y0`, `x1`, `y1`, `tool_radius`, and `cap_chord_ratio`. It validates all
values are finite, the motion is nonzero, the radius is positive, and the cap
ratio is in its closed physical domain. Each accepted binary64 value is lifted
once to a canonical exact rational numerator/denominator pair and stored in the
source. No later stage reads the original floating-point values.

Invalid inputs raise a distinct named error for nonfinite input, zero-length
motion, nonpositive radius, or invalid cap domain. A future rational-input
factory is out of scope.

## Records

`SegmentEventPartition2` binds:

- the exact-lifted source;
- the canonical Task 4 boundary feature IDs used to construct it;
- all supported line/circle pullback records;
- trim-accepted active branches;
- globally ordered cells and event fibres;
- simultaneous events and endpoint-order events;
- every ordered active branch pair's orientation and cap disposition;
- one canonical byte string and SHA-256 digest.

`SegmentEventStratum2` represents either an open cell or an event fibre.
`BranchPairDisposition2` identifies an ordered branch pair and records its
orientation and cap disposition. Algebraic multiplicity remains separate from
geometric contact disposition.

## Construction

`construct_segment_event_partition(stock, source)` extracts Task 4 records,
orders them by canonical feature identity, derives each line/circle pullback
from the stored exact source, applies exact trim acceptance, partitions every
projection, globally coalesces equal exact roots, and constructs ordered cell
and fibre strata. Event fibres retain all simultaneous events. Each stratum
contains the complete active-branch set and ordered pair dispositions.

The convenience binary64 overload delegates immediately to
`SegmentEventSource2::from_binary64`.

Unsupported algebraic degree, unsupported trim geometry, or an incomplete
partition raises its named substrate error. There is no sampling, fallback, or
soft unresolved construction path.

## Verification

`verify_segment_event_partition(stock, source, candidate)` reconstructs the
entire partition from the bound stock and exact source. It verifies the claimed
canonical bytes and digest against both the candidate records and the
independently reconstructed records. A mismatch returns
`UNRESOLVED_DEGENERACY`; construction errors remain loud named failures.

Verification never trusts diagnostic isolating intervals and never decides
from sampled values.

## Tests

Literal native-facing tests prove:

- exact binary64 lifting and all named validation failures;
- deterministic construction independent of Task 4 record insertion order;
- line and circle support pullbacks originate from the bound stock/source;
- exact global root coalescing retains simultaneous events;
- endpoint-order events and ordered pair dispositions are present;
- every open cell and event fibre carries complete active branch IDs;
- canonical bytes/digest bind the exact source and Task 4 feature identities;
- deleting or altering a source, branch, event, pair disposition, cell, fibre,
  byte string, or digest makes verification unresolved.

The focused Task 5 suite, affected tests, Ruff, and strict mypy must pass before
the segment substrate commit.
