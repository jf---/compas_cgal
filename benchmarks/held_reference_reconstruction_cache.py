"""Content-addressed on-disk store of certified Held reconstructions.

One entry per :class:`ReconstructionBuildIdentity`, at ``<directory>/<key>.json``.
An entry is used only if its format, embedded identity and integrity digest
all match; anything else raises and is never used. Writes are atomic
(temporary file plus rename), so concurrent writers of the same identity are
safe: they write identical bytes.

Entries of other code or another cache format are never touched by a load;
they are removed only by an explicit :func:`prune_stale_entries`, also run as
``python -m benchmarks.held_reference_reconstruction_cache prune``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import tempfile
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Self

from benchmarks.held_reference_geometry import ReferenceReconstruction
from benchmarks.held_reference_reconstruction_cache_errors import CorruptedReconstructionCacheEntryError
from benchmarks.held_reference_reconstruction_cache_errors import MalformedReconstructionEncodingError
from benchmarks.held_reference_reconstruction_cache_errors import ReconstructionCacheIdentityMismatchError
from benchmarks.held_reference_reconstruction_cache_errors import ReconstructionCacheWriteError
from benchmarks.held_reference_reconstruction_cache_errors import UnsupportedReconstructionCacheFormatError
from benchmarks.held_reference_reconstruction_codec import decode_reconstruction
from benchmarks.held_reference_reconstruction_codec import encode_reconstruction
from benchmarks.held_reference_reconstruction_identity import RECONSTRUCTION_CACHE_FORMAT
from benchmarks.held_reference_reconstruction_identity import ReconstructionBuildIdentity
from benchmarks.held_reference_reconstruction_identity import canonical_json_bytes
from benchmarks.held_reference_reconstruction_identity import reconstruction_code_digest

_ENTRY_KEYS = frozenset({"format", "identity", "reconstruction_sha256", "reconstruction"})
_ENTRY_SUFFIX = ".json"

DEFAULT_CACHE_DIRECTORY = Path(__file__).resolve().parents[1] / "build" / "held-case-cache"
"""The gitignored default cache directory, under the repository's build tree."""


@dataclass(frozen=True)
class ReconstructionCache:
    """A directory of reconstruction entries; nothing touches disk until used."""

    directory: Path

    @classmethod
    def build(cls, directory: Path) -> Self:
        """Anchor a cache at ``directory`` (resolved; created on first store)."""
        return cls(Path(directory).resolve())

    def entry_path(self, identity: ReconstructionBuildIdentity) -> Path:
        """Where the entry for ``identity`` lives."""
        return self.directory / f"{identity.key}.json"

    def load(self, identity: ReconstructionBuildIdentity) -> ReferenceReconstruction | None:
        """Return the cached reconstruction, or ``None`` when no entry exists.

        Raises:
            CorruptedReconstructionCacheEntryError: Unreadable, not JSON, wrong shape, digest mismatch, or undecodable.
            UnsupportedReconstructionCacheFormatError: The entry was written in another cache format.
            ReconstructionCacheIdentityMismatchError: The entry's embedded identity is not ``identity``.
        """
        path = self.entry_path(identity)
        try:
            data = path.read_bytes()
        except FileNotFoundError:
            return None
        except OSError as error:
            raise CorruptedReconstructionCacheEntryError(f"{path}: unreadable: {error}") from error
        try:
            entry = json.loads(data)
        except (UnicodeError, json.JSONDecodeError) as error:
            raise CorruptedReconstructionCacheEntryError(f"{path}: not JSON: {error}") from error
        if not isinstance(entry, dict) or set(entry) != _ENTRY_KEYS:
            raise CorruptedReconstructionCacheEntryError(f"{path}: entry must have exactly the keys {sorted(_ENTRY_KEYS)}.")
        if entry["format"] != RECONSTRUCTION_CACHE_FORMAT:
            raise UnsupportedReconstructionCacheFormatError(f"{path}: format {entry['format']!r}, expected {RECONSTRUCTION_CACHE_FORMAT!r}.")
        if entry["identity"] != identity.to_json():
            raise ReconstructionCacheIdentityMismatchError(f"{path}: embedded identity {entry['identity']!r} is not {identity.to_json()!r}.")
        encoded = entry["reconstruction"]
        if hashlib.sha256(canonical_json_bytes(encoded)).hexdigest() != entry["reconstruction_sha256"]:
            raise CorruptedReconstructionCacheEntryError(f"{path}: reconstruction fails its integrity digest.")
        try:
            return decode_reconstruction(encoded)
        except MalformedReconstructionEncodingError as error:
            raise CorruptedReconstructionCacheEntryError(f"{path}: {error}") from error

    def store(self, identity: ReconstructionBuildIdentity, reconstruction: ReferenceReconstruction) -> Path:
        """Write the entry for ``identity`` atomically and return its path.

        Raises:
            ReconstructionCacheWriteError: The directory or the entry cannot be written.
        """
        encoded = encode_reconstruction(reconstruction)
        entry = {
            "format": RECONSTRUCTION_CACHE_FORMAT,
            "identity": identity.to_json(),
            "reconstruction_sha256": hashlib.sha256(canonical_json_bytes(encoded)).hexdigest(),
            "reconstruction": encoded,
        }
        path = self.entry_path(identity)
        temporary: Path | None = None
        try:
            self.directory.mkdir(parents=True, exist_ok=True)
            with tempfile.NamedTemporaryFile(dir=self.directory, prefix=f".{identity.key}.", suffix=".tmp", delete=False) as handle:
                temporary = Path(handle.name)
                handle.write(canonical_json_bytes(entry))
            temporary.replace(path)
        except OSError as error:
            if temporary is not None:
                temporary.unlink(missing_ok=True)
            raise ReconstructionCacheWriteError(f"{path}: cannot write entry: {error}") from error
        return path


@dataclass(frozen=True)
class ReconstructionCachePruneReport:
    """Outcome of one prune: entries of the running code and format kept, all others removed."""

    directory: Path
    kept: int
    removed: int

    def __str__(self) -> str:
        return f"{self.directory}: kept {self.kept}, removed {self.removed}"


def _entry_is_stale(path: Path, current_code_sha256: str) -> bool:
    """Whether the entry at ``path`` belongs to other code or another cache format.

    Raises:
        CorruptedReconstructionCacheEntryError: Unreadable, not JSON, or the wrong shape.
        ReconstructionCacheIdentityMismatchError: A current-format entry is filed under a key that is not its identity's.
    """
    try:
        entry = json.loads(path.read_bytes())
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise CorruptedReconstructionCacheEntryError(f"{path}: unreadable or not JSON: {error}") from error
    if not isinstance(entry, dict) or "format" not in entry:
        raise CorruptedReconstructionCacheEntryError(f"{path}: entry is not an object with a format.")
    if entry["format"] != RECONSTRUCTION_CACHE_FORMAT:
        return True
    identity = entry["identity"] if set(entry) == _ENTRY_KEYS else None
    code = identity.get("code_sha256") if isinstance(identity, dict) else None
    if not isinstance(code, str):
        raise CorruptedReconstructionCacheEntryError(f"{path}: entry must have exactly the keys {sorted(_ENTRY_KEYS)} and a code digest.")
    if hashlib.sha256(canonical_json_bytes(identity)).hexdigest() != path.stem:
        raise ReconstructionCacheIdentityMismatchError(f"{path}: embedded identity {identity!r} is not filed under its key.")
    return code != current_code_sha256


def prune_stale_entries(cache_dir: Path) -> ReconstructionCachePruneReport:
    """Remove every entry not written by the running code in the current cache format.

    Every entry is classified before anything is removed, so a corrupted or
    misfiled entry aborts the prune with the directory untouched. In-flight
    temporaries of concurrent writers are not entries and are left alone. A
    missing directory holds no entries.

    Args:
        cache_dir: The cache directory.

    Returns:
        How many entries were kept and how many removed.

    Raises:
        CorruptedReconstructionCacheEntryError: An entry is unreadable, not JSON, or the wrong shape.
        ReconstructionCacheIdentityMismatchError: A current-format entry is misfiled.
        ReconstructionCacheWriteError: A stale entry cannot be removed.
    """
    directory = Path(cache_dir).resolve()
    if not directory.exists():
        return ReconstructionCachePruneReport(directory=directory, kept=0, removed=0)
    current = reconstruction_code_digest().hex()
    entries = sorted(path for path in directory.iterdir() if path.suffix == _ENTRY_SUFFIX and not path.name.startswith("."))
    stale = [path for path in entries if _entry_is_stale(path, current)]
    for path in stale:
        try:
            path.unlink(missing_ok=True)
        except OSError as error:
            raise ReconstructionCacheWriteError(f"{path}: cannot remove stale entry: {error}") from error
    return ReconstructionCachePruneReport(directory=directory, kept=len(entries) - len(stale), removed=len(stale))


def main(argv: Sequence[str] | None = None) -> int:
    """Command line: ``prune [--cache-dir DIR]`` prints the prune report."""
    parser = argparse.ArgumentParser(prog="python -m benchmarks.held_reference_reconstruction_cache", description="Maintain the Held reconstruction cache.")
    commands = parser.add_subparsers(dest="command", required=True)
    prune = commands.add_parser("prune", help="remove entries of other code or another cache format")
    prune.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE_DIRECTORY)
    arguments = parser.parse_args(argv)
    print(prune_stale_entries(arguments.cache_dir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
