"""Content-addressed identity of the executable audit proof stack."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from typing import Final
from typing import NewType
from typing import Self

from compas_cgal.adaptive.canonical import encode_component_map
from compas_cgal.adaptive.canonical import encode_sequence
from compas_cgal.adaptive.canonical import encode_tagged_union
from compas_cgal.adaptive.identity import ComponentIdentity
from compas_cgal.adaptive.identity import IdentityDigest
from compas_cgal.adaptive.identity import NativeSourceTreeDigest
from compas_cgal.engagement_audit.errors import DuplicateBuildComponentError
from compas_cgal.engagement_audit.errors import InvalidBuildIdentityError

BUILD_IDENTITY_VERSION: Final[bytes] = b"build-identity-v1"
PythonSourceTreeDigest = NewType("PythonSourceTreeDigest", bytes)
PixiLockDigest = NewType("PixiLockDigest", bytes)


def _require_sha256_digest(value: object, name: str) -> bytes:
    if type(value) is not bytes or len(value) != hashlib.sha256().digest_size:
        raise InvalidBuildIdentityError(f"{name} must be exactly one 32-byte SHA-256 digest.")
    return value


def _validated_components(
    components: object,
    *,
    require_canonical_order: bool,
) -> tuple[ComponentIdentity, ...]:
    if type(components) is not tuple:
        raise InvalidBuildIdentityError("build components must be an exact tuple.")
    if not components:
        raise InvalidBuildIdentityError("build identity requires at least one component.")
    if any(type(component) is not ComponentIdentity for component in components):
        raise InvalidBuildIdentityError("every build component must be an exact ComponentIdentity.")

    typed_components = components
    ordered = tuple(sorted(typed_components, key=lambda component: component.component_domain))
    domains = tuple(component.component_domain for component in ordered)
    for left, right in zip(domains, domains[1:]):
        if left == right:
            domain = bytes(left).decode("utf-8", errors="backslashreplace")
            raise DuplicateBuildComponentError(f"build component domain {domain!r} is duplicated.")
    if require_canonical_order and typed_components != ordered:
        raise InvalidBuildIdentityError("raw build components must already use canonical domain order.")
    return ordered


@dataclass(frozen=True)
class BuildIdentity:
    """Complete source, lock, and component identity for one audit build."""

    components: tuple[ComponentIdentity, ...]
    native_source_tree_digest: NativeSourceTreeDigest
    python_source_tree_digest: PythonSourceTreeDigest
    pixi_lock_digest: PixiLockDigest

    def __post_init__(self) -> None:
        _validated_components(self.components, require_canonical_order=True)
        _require_sha256_digest(self.native_source_tree_digest, "native source-tree digest")
        _require_sha256_digest(self.python_source_tree_digest, "Python source-tree digest")
        _require_sha256_digest(self.pixi_lock_digest, "pixi.lock digest")

    @classmethod
    def build(
        cls,
        *,
        components: tuple[ComponentIdentity, ...],
        native_source_tree_digest: NativeSourceTreeDigest,
        python_source_tree_digest: PythonSourceTreeDigest,
        pixi_lock_digest: PixiLockDigest,
    ) -> Self:
        ordered = _validated_components(components, require_canonical_order=False)
        return cls(
            components=ordered,
            native_source_tree_digest=native_source_tree_digest,
            python_source_tree_digest=python_source_tree_digest,
            pixi_lock_digest=pixi_lock_digest,
        )

    @property
    def canonical_bytes(self) -> bytes:
        if type(self) is not BuildIdentity:
            raise InvalidBuildIdentityError("build identity must be exact BuildIdentity, not a subclass.")
        return encode_tagged_union(
            BUILD_IDENTITY_VERSION,
            encode_component_map(
                {
                    b"components": encode_sequence(tuple(component.canonical_bytes for component in self.components)),
                    b"native-source-tree-digest": bytes(self.native_source_tree_digest),
                    b"pixi-lock-digest": bytes(self.pixi_lock_digest),
                    b"python-source-tree-digest": bytes(self.python_source_tree_digest),
                }
            ),
        )

    @property
    def digest(self) -> IdentityDigest:
        return IdentityDigest(hashlib.sha256(self.canonical_bytes).digest())
