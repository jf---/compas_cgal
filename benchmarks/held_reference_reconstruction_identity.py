"""Build identity of one cached Held reconstruction: inputs hash plus version components.

The cached value is ``reconstruct_source_path`` applied to inputs that the
case loader derives from the case document. It is a pure function of

* the case document bytes (``document_digest``),
* the Python source that parses the document and decides the result
  (``code_digest``, a hash over :data:`DECIDING_SOURCES`),
* the interpreter (``fractions``/``math`` semantics) and the platform
  (``libm`` results for ``sin``/``cos``/``atan2``),

so the cache key is a hash over exactly those components. The code version is
derived from source bytes rather than a hand-bumped constant: a hand-bumped
version needs a guard test that pins the very same source hash, so it adds a
manual step without adding a guarantee. The cost is that a cosmetic edit to a
deciding module invalidates the cache (one cold load per case). Which modules
decide the result is enforced by an import-closure test.

The digest describes the code as IMPORTED, not as on disk at first use: the
source bytes of every deciding module are snapshotted when this module is
imported, which imports every deciding module first, and nothing re-reads
them afterwards. An edit between import and the first case load therefore
cannot relabel the running code. The residual window is the import statement
itself (milliseconds, the same window Python's own pyc mtime check has).
``importlib.reload`` of a deciding module is not tracked: reload, then
re-import this module, or start a fresh interpreter.
"""

from __future__ import annotations

import hashlib
import importlib
import json
import platform
import sys
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import NewType
from typing import Self

from benchmarks.held_reference_reconstruction_cache_errors import InvalidReconstructionBuildIdentityError

RECONSTRUCTION_CACHE_FORMAT = "held-reference-reconstruction-cache/v1"

CaseDocumentDigest = NewType("CaseDocumentDigest", bytes)
ReconstructionCodeDigest = NewType("ReconstructionCodeDigest", bytes)

_SHA256_SIZE = hashlib.sha256().digest_size
_PYTHON_SOURCE_SUFFIX = ".py"


@dataclass(frozen=True)
class DecidingSource:
    """One module whose source bytes, as imported, are a version component of the cached result."""

    module: str
    path: Path
    source_digest: bytes

    @classmethod
    def of_imported(cls, module: ModuleType) -> Self:
        """Snapshot the source of an imported module; call it right after the import.

        Raises:
            InvalidReconstructionBuildIdentityError: The module has no readable ``.py`` source file.
        """
        file = getattr(module, "__file__", None)
        if not isinstance(file, str):
            raise InvalidReconstructionBuildIdentityError(f"deciding module {module.__name__!r} has no source file.")
        path = Path(file).resolve()
        if path.suffix != _PYTHON_SOURCE_SUFFIX:
            raise InvalidReconstructionBuildIdentityError(f"deciding module {module.__name__!r} is not Python source: {path}.")
        try:
            data = path.read_bytes()
        except OSError as error:
            raise InvalidReconstructionBuildIdentityError(f"cannot read deciding module {module.__name__!r}: {error}") from error
        return cls(module.__name__, path, hashlib.sha256(data).digest())


DECIDING_MODULES: tuple[str, ...] = (
    "benchmarks.held_reference_cases",
    "benchmarks.held_reference_certification",
    "benchmarks.held_reference_geometry",
    "benchmarks.held_reference_reconstruction_codec",
    "benchmarks.held_reference_reconstruction_identity",
    "compas_cgal.units",
)
"""Modules that parse the document into reconstruction inputs or compute and encode the result."""

NON_DECIDING_FIRST_PARTY_IMPORTS: dict[str, str] = {
    "benchmarks.errors": "exception types only",
    "benchmarks.held_reference_reconstruction_cache_errors": "exception types only",
    "benchmarks.held_reference_reconstruction_cache": "moves entry bytes; the format constant and integrity digest version it",
    "benchmarks.spec": "PocketSpec, built from a validated case after reconstruction",
    "compas_cgal.errors": "exception types only",
}
"""First-party imports of deciding modules that cannot change the cached value, with the reason."""


def code_digest_of(sources: Iterable[DecidingSource]) -> ReconstructionCodeDigest:
    """SHA-256 over each source's module name and snapshotted source digest, in module order."""
    digest = hashlib.sha256()
    for source in sorted(sources, key=lambda item: item.module):
        digest.update(source.module.encode())
        digest.update(b"\0")
        digest.update(source.source_digest)
    return ReconstructionCodeDigest(digest.digest())


def reconstruction_code_digest() -> ReconstructionCodeDigest:
    """Digest of the deciding sources as imported into this process (see the module docstring)."""
    return _RUNNING_CODE_DIGEST


def current_interpreter() -> str:
    """Interpreter version component: implementation plus full version string."""
    return f"{sys.implementation.name}:{sys.version}"


def current_platform() -> str:
    """Platform version component, including the OS release that ships ``libm``."""
    return platform.platform()


def _require_digest(value: object, name: str) -> bytes:
    if type(value) is not bytes or len(value) != _SHA256_SIZE:
        raise InvalidReconstructionBuildIdentityError(f"{name} must be one {_SHA256_SIZE}-byte SHA-256 digest.")
    return value


def _require_label(value: object, name: str) -> str:
    if type(value) is not str or not value:
        raise InvalidReconstructionBuildIdentityError(f"{name} must be a non-empty string.")
    return value


@dataclass(frozen=True)
class ReconstructionBuildIdentity:
    """Content address of one cached reconstruction."""

    document_digest: CaseDocumentDigest
    code_digest: ReconstructionCodeDigest
    interpreter: str
    platform: str

    @classmethod
    def build(
        cls,
        *,
        document_digest: CaseDocumentDigest,
        code_digest: ReconstructionCodeDigest,
        interpreter: str,
        platform: str,
    ) -> Self:
        """Validate every component.

        Raises:
            InvalidReconstructionBuildIdentityError: A digest is not 32 bytes or a label is empty.
        """
        return cls(
            CaseDocumentDigest(_require_digest(document_digest, "document digest")),
            ReconstructionCodeDigest(_require_digest(code_digest, "code digest")),
            _require_label(interpreter, "interpreter"),
            _require_label(platform, "platform"),
        )

    @classmethod
    def for_document(cls, document: bytes) -> Self:
        """Identity of ``document`` under the running code, interpreter and platform."""
        return cls.build(
            document_digest=CaseDocumentDigest(hashlib.sha256(document).digest()),
            code_digest=reconstruction_code_digest(),
            interpreter=current_interpreter(),
            platform=current_platform(),
        )

    def to_json(self) -> dict[str, str]:
        """Canonical JSON form, embedded in each entry and hashed into :attr:`key`."""
        return {
            "format": RECONSTRUCTION_CACHE_FORMAT,
            "document_sha256": self.document_digest.hex(),
            "code_sha256": self.code_digest.hex(),
            "interpreter": self.interpreter,
            "platform": self.platform,
        }

    @property
    def key(self) -> str:
        """Hex SHA-256 of the canonical identity: the entry's file stem."""
        return hashlib.sha256(canonical_json_bytes(self.to_json())).hexdigest()


def canonical_json_bytes(value: object) -> bytes:
    """Deterministic JSON bytes: sorted keys, no whitespace, ASCII only."""
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False).encode("ascii")


# Snapshot last: importing the deciding modules may re-enter this module (the
# case loader imports it), which then needs every definition above.
DECIDING_SOURCES: tuple[DecidingSource, ...] = tuple(DecidingSource.of_imported(importlib.import_module(name)) for name in DECIDING_MODULES)
"""The deciding modules' sources, snapshotted at import."""

_RUNNING_CODE_DIGEST = code_digest_of(DECIDING_SOURCES)
