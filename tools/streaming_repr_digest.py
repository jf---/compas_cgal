"""SHA-256 of ``repr(value)``, streamed: the digest is byte-identical, the string is never built.

The byte-identity digests (``tools/held_path_digest.py``, the declared-clearance
differential) hash ``repr`` of a generated path. On Monstera that text is
722 M characters, and ``repr`` keeps the string, the per-element parts its
``join`` builds and the ``encode`` copy alive at once (about 2.9 GB of peak RSS).

``repr_sha256`` feeds SHA-256 the same characters in order, chunk by chunk:

- ``tuple`` and ``list`` (exact types) print their brackets, the ``", "``
  separator and the one-tuple ``",)"`` themselves, as CPython does, including
  ``(...)`` / ``[...]`` for a container already being printed;
- a dataclass whose ``__repr__`` is the one ``dataclasses`` generated prints
  ``Qualname(f=..., ...)`` over its ``repr=True`` fields, and ``...`` on
  re-entry, as ``reprlib.recursive_repr`` does;
- every other value is a leaf and is printed by its own ``repr``, and so is
  a container or dataclass whose children are all leaves (its ``repr`` is
  exactly the streamed text, and only one level deep).

The memory-address guard of the callers (a ``repr`` carrying ``" at 0x...>"``
differs run to run) runs on the streamed text, across chunk boundaries.

Not reproduced: a reference cycle that passes through a leaf (a ``dict``, a
``namedtuple``, a custom ``__repr__``) back into a streamed container prints
that container in full once more inside the leaf, where ``repr`` would print
the re-entry marker. Generated paths are acyclic frozen values.
"""

from __future__ import annotations

import dataclasses
import functools
import hashlib
import re

# ``object.__repr__`` and friends: "<module.Class object at 0x10a2b3c4d>".
MEMORY_ADDRESS = re.compile(r" at 0x[0-9a-fA-F]+>")
# A suffix of the streamed text that a later chunk could complete into a
# ``MEMORY_ADDRESS`` match: a non-empty prefix of " at 0x<hex>*>" without the ">".
_ADDRESS_PREFIX_AT_END = re.compile(r" (?:a(?:t(?: (?:0(?:x[0-9a-fA-F]*)?)?)?)?)?\Z")
# Characters of text before an address match quoted in the error, for the reader.
CONTEXT_CHARS = 80
# Text parts (brackets, separators, leaf reprs) collected before one hash update:
# 4 Ki parts of a few to a few hundred characters keep a chunk near a megabyte at
# most, and the per-update overhead negligible.
DEFAULT_CHUNK_PARTS = 1 << 12
# The source file name ``dataclasses`` compiles its generated methods under (``exec`` of a string).
_GENERATED_SOURCE = "<string>"


class ReprCarriesMemoryAddressError(Exception):
    """The streamed ``repr`` carries a memory address, so it differs run to run.

    Attributes:
        context: The match with up to ``CONTEXT_CHARS`` characters of text before it.
    """

    def __init__(self, context: str) -> None:
        super().__init__(f"repr is not canonical, it carries a memory address: ...{context}")
        self.context = context


class AddressGuard:
    """Refuse a ``MEMORY_ADDRESS`` match in text that arrives in pieces.

    The tail kept between pieces is the longer of the quoted context and the
    suffix a later piece could complete into a match, so a match is found
    exactly when it occurs in the concatenated text.
    """

    def __init__(self) -> None:
        self._tail = ""

    def scan(self, text: str) -> None:
        """Scan the next piece of text.

        Raises:
            ReprCarriesMemoryAddressError: the text so far contains a match.
        """
        window = self._tail + text
        found = MEMORY_ADDRESS.search(window)
        if found is not None:
            raise ReprCarriesMemoryAddressError(window[max(0, found.start() - CONTEXT_CHARS) : found.end()])
        partial = _ADDRESS_PREFIX_AT_END.search(window)
        keep = max(CONTEXT_CHARS, 0 if partial is None else len(window) - partial.start())
        self._tail = window[-keep:]


class _ChunkedSha256:
    """Collect text parts, then hash their UTF-8 encoding and guard them one chunk at a time."""

    def __init__(self, chunk_parts: int) -> None:
        self.chunk_parts = chunk_parts
        self.parts: list[str] = []
        self._hash = hashlib.sha256()
        self._guard = AddressGuard()

    def flush(self) -> None:
        text = "".join(self.parts)
        self.parts.clear()
        self._guard.scan(text)
        self._hash.update(text.encode())

    def hexdigest(self) -> str:
        self.flush()
        return self._hash.hexdigest()


@functools.lru_cache(maxsize=None)
def _generated_repr_fields(cls: type) -> tuple[str, ...] | None:
    """The fields a ``dataclasses``-generated ``__repr__`` of ``cls`` prints, or ``None`` if it has none.

    The ``__repr__`` ``cls`` resolves to must be defined by a dataclass that
    asked for one (``repr=True``) and be the ``reprlib.recursive_repr`` wrapper
    of a function ``dataclasses`` compiled from a string; the names are that
    dataclass's ``repr=True`` fields, in order, as the generated code bakes them.
    """
    owner = next((base for base in cls.__mro__ if "__repr__" in base.__dict__), None)
    if owner is None:
        return None
    params = owner.__dict__.get("__dataclass_params__")
    if params is None or not params.repr:
        return None
    inner = getattr(owner.__dict__["__repr__"], "__wrapped__", None)
    code = getattr(inner, "__code__", None)
    if code is None or code.co_filename != _GENERATED_SOURCE or code.co_name != "__repr__":
        return None
    return tuple(item.name for item in dataclasses.fields(owner) if item.repr)


@functools.lru_cache(maxsize=None)
def _is_leaf_type(kind: type) -> bool:
    """Whether values of ``kind`` are printed by their own ``repr`` rather than streamed."""
    return kind is not tuple and kind is not list and _generated_repr_fields(kind) is None


def _emit(value: object, sink: _ChunkedSha256, active: set[int]) -> None:
    """Append the ``repr`` text of ``value`` to ``sink``, flushing it after a child once a chunk is full."""
    write = sink.parts.append
    kind = type(value)
    if kind is tuple or kind is list:
        items: tuple[object, ...] | list[object] = value  # type: ignore[assignment]
        if all(_is_leaf_type(type(item)) for item in items):
            # Only leaves below: ``repr`` prints exactly the streamed text, one level deep.
            write(repr(value))
            return
        key = id(value)
        if key in active:
            write("(...)" if kind is tuple else "[...]")
            return
        active.add(key)
        write("(" if kind is tuple else "[")
        for index, item in enumerate(items):
            if index:
                write(", ")
            if _is_leaf_type(type(item)):
                write(repr(item))
            else:
                _emit(item, sink, active)
            if len(sink.parts) >= sink.chunk_parts:
                sink.flush()
        write((",)" if len(items) == 1 else ")") if kind is tuple else "]")
        active.discard(key)
        return
    names = _generated_repr_fields(kind)
    if names is None:
        write(repr(value))
        return
    values = [getattr(value, name) for name in names]
    if all(_is_leaf_type(type(field_value)) for field_value in values):
        write(repr(value))
        return
    key = id(value)
    if key in active:
        write("...")
        return
    active.add(key)
    write(value.__class__.__qualname__ + "(")
    for index, (name, field_value) in enumerate(zip(names, values)):
        write(f", {name}=" if index else f"{name}=")
        if _is_leaf_type(type(field_value)):
            write(repr(field_value))
        else:
            _emit(field_value, sink, active)
        if len(sink.parts) >= sink.chunk_parts:
            sink.flush()
    write(")")
    active.discard(key)


def repr_sha256(value: object, *, chunk_parts: int = DEFAULT_CHUNK_PARTS) -> str:
    """SHA-256 hex digest of ``repr(value).encode()``, computed without building the string.

    Args:
        value: The value to digest.
        chunk_parts: Text parts collected before one hash update; never changes the digest.

    Raises:
        ReprCarriesMemoryAddressError: the ``repr`` carries a memory address.
    """
    sink = _ChunkedSha256(chunk_parts)
    _emit(value, sink, set())
    return sink.hexdigest()
