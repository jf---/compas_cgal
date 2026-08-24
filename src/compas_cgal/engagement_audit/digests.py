"""Typed SHA-256 domains for the engagement-audit boundary."""

from typing import NewType

AuditInputDigest = NewType("AuditInputDigest", bytes)
AuditNativeRequestDigest = NewType("AuditNativeRequestDigest", bytes)
AuthenticatedOperationDigest = NewType("AuthenticatedOperationDigest", bytes)
NativeMotionDigest = NewType("NativeMotionDigest", bytes)
AuditPolicyDigest = NewType("AuditPolicyDigest", bytes)
