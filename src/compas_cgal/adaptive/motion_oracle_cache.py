"""Content-addressed memo for the exact native lateral-motion TEA oracle.

`MotionCertifier.certify` spends nearly all of its time inside one native
call. Route retrace and candidate retries re-issue that call for byte-identical
stock content and byte-identical motions, from *different* certifier instances,
so the memo lives at module scope and is keyed on content identity rather than
on object identity of the certifier.

The memo is deliberately narrow. It stores only the `(verdict, trace)` pair the
native oracle returns. Everything a caller adds around that pair -- the
operation ordinal, the user cap, the witness -- is rebuilt on every call, so a
memo hit can never lend one operation's ordinal to another.

Three properties keep the memo semantically invisible:

- The key carries every argument the native oracle consumes, encoded with the
  frozen CCAN grammar. It is the complete canonical record, not a hash of it,
  so distinct inputs cannot collide at all.
- An entry belongs to the exact oracle object that produced it. Replacing the
  native entry point (a substitute oracle, a fault injection) yields a fresh
  key, so a substitute is always consulted rather than answered from a result
  the real oracle returned earlier.
- Raised native exceptions are never stored. Only a returned pair is memoized.
"""

from __future__ import annotations

import hashlib
from collections import OrderedDict
from collections.abc import Callable
from typing import Final
from typing import Literal
from typing import TypeAlias

from compas_cgal import _continuous_tea_2
from compas_cgal import _stock_2
from compas_cgal.adaptive.canonical import canonical_task1_bytes
from compas_cgal.adaptive.canonical import encode_binary64
from compas_cgal.adaptive.canonical import encode_bytes
from compas_cgal.adaptive.canonical import encode_component_map
from compas_cgal.adaptive.canonical import encode_tagged_union
from compas_cgal.adaptive.errors import InvalidMotionOracleCacheKeyError
from compas_cgal.adaptive.motion import EngagementCap
from compas_cgal.adaptive.motion import ExactCircleMotion
from compas_cgal.adaptive.motion import ExactSegmentMotion
from compas_cgal.adaptive.units import ToolRadius

NativeMotionVerdict: TypeAlias = Literal["certified", "cap_exceeded", "unresolved"]
NativeMotionAudit: TypeAlias = tuple[NativeMotionVerdict, "_continuous_tea_2.EventTrace2"]
NativeMotionOracle: TypeAlias = Callable[..., NativeMotionAudit]
_MemoKey: TypeAlias = tuple[NativeMotionOracle, bytes]

MOTION_ORACLE_CACHE_KEY_VERSION: Final[bytes] = b"motion-oracle-audit-key-v1"
SEGMENT_ORACLE_TAG: Final[bytes] = b"audit-segment-tea-event-exact-v1"
CIRCLE_ORACLE_TAG: Final[bytes] = b"audit-full-circle-tea-event-exact-v1"
# One operation certifies a bounded candidate family, and route retrace revisits
# only the recent tail of that history; 512 entries hold several operations of
# certification history while capping the native traces one process retains.
MOTION_ORACLE_CACHE_CAPACITY: Final[int] = 512

_DIGEST_SIZE: Final[int] = hashlib.sha256().digest_size
_ORACLE_TAGS: Final[frozenset[bytes]] = frozenset((SEGMENT_ORACLE_TAG, CIRCLE_ORACLE_TAG))
_NATIVE_MOTION_AUDIT_CACHE: Final[OrderedDict[_MemoKey, NativeMotionAudit]] = OrderedDict()


def _require_content_digest(value: bytes, name: str) -> None:
    if type(value) is not bytes or len(value) != _DIGEST_SIZE:
        raise InvalidMotionOracleCacheKeyError(f"{name} must be one exact SHA-256 digest.")


def _require_key_inputs(
    motion: ExactSegmentMotion | ExactCircleMotion,
    tool_radius: ToolRadius,
    effective_cap: EngagementCap,
    stock_lineage_digest: bytes,
    canonical_boundary_digest: bytes,
) -> None:
    _require_content_digest(stock_lineage_digest, "memo stock lineage digest")
    _require_content_digest(canonical_boundary_digest, "memo stock boundary digest")
    if type(motion) not in (ExactSegmentMotion, ExactCircleMotion):
        raise InvalidMotionOracleCacheKeyError("memo key requires one exact motion type.")
    if type(tool_radius) is not ToolRadius or type(effective_cap) is not EngagementCap:
        raise InvalidMotionOracleCacheKeyError("memo key requires an exact tool radius and engagement cap.")


def canonical_motion_audit_key_bytes(
    *,
    stock_lineage_digest: bytes,
    canonical_boundary_digest: bytes,
    motion: ExactSegmentMotion | ExactCircleMotion,
    tool_radius: ToolRadius,
    effective_cap: EngagementCap,
    oracle_tag: bytes,
) -> bytes:
    """Encode every input the native motion oracle consumes.

    The two stock digests identify the exact native stock content, the motion
    record carries its exact binary64 coordinates, and the cap carries the
    native binary64 chord-ratio surrogate. The caller's operation ordinal and
    user cap are absent on purpose: neither reaches the native oracle.

    Args:
        stock_lineage_digest: SHA-256 identity of the observed depletion lineage.
        canonical_boundary_digest: SHA-256 identity of the exact stock boundary.
        motion: Exact segment or full-circle motion.
        tool_radius: Typed cutter radius passed to the oracle.
        effective_cap: Exact policy-derived cap passed to the oracle.
        oracle_tag: Closed tag naming the native entry point.

    Returns:
        Complete canonical CCAN record binding every native argument.

    Raises:
        InvalidMotionOracleCacheKeyError: If any key input is not one exact
            owned value.
    """
    _require_key_inputs(
        motion,
        tool_radius,
        effective_cap,
        stock_lineage_digest,
        canonical_boundary_digest,
    )
    if type(oracle_tag) is not bytes or oracle_tag not in _ORACLE_TAGS:
        raise InvalidMotionOracleCacheKeyError("memo key requires one closed native entry-point tag.")
    return encode_tagged_union(
        MOTION_ORACLE_CACHE_KEY_VERSION,
        encode_component_map(
            {
                b"canonical-boundary-digest": encode_bytes(canonical_boundary_digest),
                b"effective-cap": canonical_task1_bytes(effective_cap),
                b"motion": canonical_task1_bytes(motion),
                b"native-oracle": encode_tagged_union(oracle_tag, b""),
                b"stock-lineage-digest": encode_bytes(stock_lineage_digest),
                b"tool-radius": encode_tagged_union(
                    b"tool-radius-mm-v1",
                    encode_binary64(float(tool_radius.value)),
                ),
            }
        ),
    )


def _segment_dispatch(
    stock: _stock_2.Stock2,
    motion: ExactSegmentMotion,
    tool_radius: ToolRadius,
    effective_cap: EngagementCap,
) -> tuple[NativeMotionOracle, bytes, tuple[object, ...]]:
    return (
        _continuous_tea_2.audit_segment_tea_event_exact,
        SEGMENT_ORACLE_TAG,
        (
            stock,
            motion.start.x,
            motion.start.y,
            motion.end.x,
            motion.end.y,
            tool_radius.value,
            effective_cap.chord_ratio,
        ),
    )


def _circle_dispatch(
    stock: _stock_2.Stock2,
    motion: ExactCircleMotion,
    tool_radius: ToolRadius,
    effective_cap: EngagementCap,
) -> tuple[NativeMotionOracle, bytes, tuple[object, ...]]:
    return (
        _continuous_tea_2.audit_full_circle_tea_event_exact,
        CIRCLE_ORACLE_TAG,
        (
            stock,
            motion.center.x,
            motion.center.y,
            motion.phase_vector.x,
            motion.phase_vector.y,
            motion.clockwise,
            tool_radius.value,
            effective_cap.chord_ratio,
        ),
    )


def audit_motion_tea_event_exact(
    *,
    stock: _stock_2.Stock2,
    motion: ExactSegmentMotion | ExactCircleMotion,
    tool_radius: ToolRadius,
    effective_cap: EngagementCap,
    stock_lineage_digest: bytes,
    canonical_boundary_digest: bytes,
) -> NativeMotionAudit:
    """Return the exact native TEA audit for one lateral motion.

    A memoized pair is returned unchanged; a miss dispatches to the native
    oracle and stores its returned pair. Native exceptions propagate to the
    caller and are never stored, so every failing request re-runs the oracle
    and keeps its own traceback and cause chain.

    Args:
        stock: Immutable native stock snapshot the oracle reads.
        motion: Exact segment or full-circle motion.
        tool_radius: Typed cutter radius.
        effective_cap: Exact policy-derived engagement cap.
        stock_lineage_digest: SHA-256 identity of the observed depletion lineage.
        canonical_boundary_digest: SHA-256 identity of the exact stock boundary.

    Returns:
        The native `(verdict, trace)` pair, memoized or freshly computed.

    Raises:
        InvalidMotionOracleCacheKeyError: If a key input or the declared memo
            capacity is not one exact owned value.
    """
    capacity = MOTION_ORACLE_CACHE_CAPACITY
    if type(capacity) is not int or capacity < 1:
        raise InvalidMotionOracleCacheKeyError("motion-oracle memo capacity must be one positive exact integer.")
    if type(stock) is not _stock_2.Stock2:
        raise InvalidMotionOracleCacheKeyError("motion audit requires one owned native Stock2 snapshot.")
    _require_key_inputs(
        motion,
        tool_radius,
        effective_cap,
        stock_lineage_digest,
        canonical_boundary_digest,
    )
    if type(motion) is ExactSegmentMotion:
        oracle, oracle_tag, arguments = _segment_dispatch(
            stock,
            motion,
            tool_radius,
            effective_cap,
        )
    elif type(motion) is ExactCircleMotion:
        oracle, oracle_tag, arguments = _circle_dispatch(
            stock,
            motion,
            tool_radius,
            effective_cap,
        )
    else:
        raise InvalidMotionOracleCacheKeyError("motion audit requires one exact motion type.")
    key: _MemoKey = (
        oracle,
        canonical_motion_audit_key_bytes(
            stock_lineage_digest=stock_lineage_digest,
            canonical_boundary_digest=canonical_boundary_digest,
            motion=motion,
            tool_radius=tool_radius,
            effective_cap=effective_cap,
            oracle_tag=oracle_tag,
        ),
    )
    memoized = _NATIVE_MOTION_AUDIT_CACHE.get(key)
    if memoized is not None:
        _NATIVE_MOTION_AUDIT_CACHE.move_to_end(key)
        return memoized
    audit = oracle(*arguments)
    _NATIVE_MOTION_AUDIT_CACHE[key] = audit
    while len(_NATIVE_MOTION_AUDIT_CACHE) > capacity:
        _NATIVE_MOTION_AUDIT_CACHE.popitem(last=False)
    return audit


def clear_native_motion_audit_cache() -> None:
    """Drop every memoized native audit.

    Discarding memoized results can never change a verdict, only the work spent
    re-deriving it. Callers use this to measure oracle dispatch deterministically
    and to release retained native traces between independent runs.
    """
    _NATIVE_MOTION_AUDIT_CACHE.clear()


def native_motion_audit_cache_size() -> int:
    """Return the number of memoized native audits.

    Returns:
        Current entry count, never above `MOTION_ORACLE_CACHE_CAPACITY`.
    """
    return len(_NATIVE_MOTION_AUDIT_CACHE)
