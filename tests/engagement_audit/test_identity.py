from __future__ import annotations

import hashlib
from dataclasses import replace

import pytest

from compas_cgal.adaptive.canonical import encode_bytes
from compas_cgal.adaptive.canonical import encode_tagged_union
from compas_cgal.adaptive.identity import ComponentDomainTag
from compas_cgal.adaptive.identity import ComponentIdentity
from compas_cgal.adaptive.identity import NativeSourceTreeDigest
from compas_cgal.adaptive.identity import SourceRevision
from compas_cgal.adaptive.identity import StrategyVersion
from compas_cgal.engagement_audit.errors import DuplicateBuildComponentError
from compas_cgal.engagement_audit.errors import InvalidBuildIdentityError
from compas_cgal.engagement_audit.identity import BuildIdentity
from compas_cgal.engagement_audit.identity import PixiLockDigest
from compas_cgal.engagement_audit.identity import PythonSourceTreeDigest


def _digest(seed: bytes) -> bytes:
    return hashlib.sha256(seed).digest()


def _component(domain: bytes, *, version: bytes = b"v1") -> ComponentIdentity:
    return ComponentIdentity.build(
        component_domain=ComponentDomainTag(domain),
        strategy_version=StrategyVersion(version),
        source_revision=SourceRevision(b"source-revision"),
        native_source_tree_digest=NativeSourceTreeDigest(_digest(domain + b"-native")),
        canonical_parameter_bytes=encode_tagged_union(b"test-parameters-v1", encode_bytes(domain)),
    )


def _build_identity(
    *,
    components: tuple[ComponentIdentity, ...] | None = None,
) -> BuildIdentity:
    return BuildIdentity.build(
        components=components or (_component(b"stock"), _component(b"engagement")),
        native_source_tree_digest=NativeSourceTreeDigest(_digest(b"native-tree")),
        python_source_tree_digest=PythonSourceTreeDigest(_digest(b"python-tree")),
        pixi_lock_digest=PixiLockDigest(_digest(b"pixi-lock")),
    )


def test_build_identity_binds_every_component_and_source_digest() -> None:
    identity = _build_identity()

    assert identity.components == tuple(sorted(identity.components, key=lambda item: item.component_domain))
    assert bytes(identity.digest) == hashlib.sha256(identity.canonical_bytes).digest()
    assert bytes(identity.native_source_tree_digest) in identity.canonical_bytes
    assert bytes(identity.python_source_tree_digest) in identity.canonical_bytes
    assert bytes(identity.pixi_lock_digest) in identity.canonical_bytes
    for component in identity.components:
        assert component.canonical_bytes in identity.canonical_bytes


def test_component_order_does_not_change_build_identity() -> None:
    stock = _component(b"stock")
    engagement = _component(b"engagement")

    forward = _build_identity(components=(stock, engagement))
    reverse = _build_identity(components=(engagement, stock))

    assert forward.canonical_bytes == reverse.canonical_bytes
    assert forward.digest == reverse.digest


def test_component_version_changes_build_identity() -> None:
    original = _build_identity(components=(_component(b"stock"),))
    changed = _build_identity(components=(_component(b"stock", version=b"v2"),))

    assert original.digest != changed.digest


def test_build_identity_rejects_duplicate_component_domains() -> None:
    component = _component(b"stock")

    with pytest.raises(DuplicateBuildComponentError, match="stock"):
        _build_identity(components=(component, component))


def test_build_identity_rejects_empty_component_set() -> None:
    with pytest.raises(InvalidBuildIdentityError, match="at least one component"):
        BuildIdentity.build(
            components=(),
            native_source_tree_digest=NativeSourceTreeDigest(_digest(b"native-tree")),
            python_source_tree_digest=PythonSourceTreeDigest(_digest(b"python-tree")),
            pixi_lock_digest=PixiLockDigest(_digest(b"pixi-lock")),
        )


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("native_source_tree_digest", b"short"),
        ("python_source_tree_digest", b"short"),
        ("pixi_lock_digest", b"short"),
    ],
)
def test_build_identity_rejects_non_sha256_digest(field: str, value: bytes) -> None:
    identity = _build_identity()

    with pytest.raises(InvalidBuildIdentityError, match="32-byte SHA-256"):
        replace(identity, **{field: value})  # type: ignore[arg-type]


def test_raw_construction_cannot_bypass_component_validation() -> None:
    with pytest.raises(InvalidBuildIdentityError, match="exact tuple"):
        BuildIdentity(
            components=[_component(b"stock")],  # type: ignore[arg-type]
            native_source_tree_digest=NativeSourceTreeDigest(_digest(b"native-tree")),
            python_source_tree_digest=PythonSourceTreeDigest(_digest(b"python-tree")),
            pixi_lock_digest=PixiLockDigest(_digest(b"pixi-lock")),
        )
