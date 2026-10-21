from __future__ import annotations

import ast
import hashlib
import importlib
import json
import math
import os
import subprocess
import sys
import uuid
from collections.abc import Sequence
from pathlib import Path

import pytest
from hypothesis import assume
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st

import benchmarks.held_reference_cases as cases_module
from benchmarks.errors import InvalidPublishedPrimitiveError
from benchmarks.errors import MalformedHeldReferenceCaseError
from benchmarks.held_reference_cases import load_held_reference_case
from benchmarks.held_reference_geometry import ReferenceArc
from benchmarks.held_reference_geometry import ReferenceLine
from benchmarks.held_reference_geometry import ReferencePrimitive
from benchmarks.held_reference_geometry import ReferenceReconstruction
from benchmarks.held_reference_geometry import SourcePrimitive
from benchmarks.held_reference_geometry import SourceToWorld
from benchmarks.held_reference_geometry import reconstruct_source_path
from benchmarks.held_reference_reconstruction_cache import DEFAULT_CACHE_DIRECTORY
from benchmarks.held_reference_reconstruction_cache import ReconstructionCache
from benchmarks.held_reference_reconstruction_cache import ReconstructionCachePruneReport
from benchmarks.held_reference_reconstruction_cache import main as cache_main
from benchmarks.held_reference_reconstruction_cache import prune_stale_entries
from benchmarks.held_reference_reconstruction_cache_errors import CorruptedReconstructionCacheEntryError
from benchmarks.held_reference_reconstruction_cache_errors import InvalidReconstructionBuildIdentityError
from benchmarks.held_reference_reconstruction_cache_errors import MalformedReconstructionEncodingError
from benchmarks.held_reference_reconstruction_cache_errors import ReconstructionCacheIdentityMismatchError
from benchmarks.held_reference_reconstruction_cache_errors import UnsupportedReconstructionCacheFormatError
from benchmarks.held_reference_reconstruction_codec import decode_reconstruction
from benchmarks.held_reference_reconstruction_codec import encode_reconstruction
from benchmarks.held_reference_reconstruction_identity import DECIDING_SOURCES
from benchmarks.held_reference_reconstruction_identity import DecidingSource
from benchmarks.held_reference_reconstruction_identity import code_digest_of
from benchmarks.held_reference_reconstruction_identity import NON_DECIDING_FIRST_PARTY_IMPORTS
from benchmarks.held_reference_reconstruction_identity import CaseDocumentDigest
from benchmarks.held_reference_reconstruction_identity import ReconstructionBuildIdentity
from benchmarks.held_reference_reconstruction_identity import ReconstructionCodeDigest
from benchmarks.held_reference_reconstruction_identity import reconstruction_code_digest
from compas_cgal.units import Millimetre
from compas_cgal.units import Point2
from compas_cgal.units import Radian
from compas_cgal.units import WorldXY

FIGURE5_DOCUMENT = (cases_module.DATA_DIRECTORY / "figure5.json").read_bytes()
FIRST_PARTY_PACKAGES = ("benchmarks", "compas_cgal", "tools")


def _float_fingerprint(reconstruction: ReferenceReconstruction) -> tuple[object, ...]:
    """Bitwise identity: equality plus the sign of zero, which ``==`` ignores."""

    def point(value: Point2[WorldXY]) -> tuple[str, str]:
        return (float(value.x).hex(), float(value.y).hex())

    def primitive(value: ReferencePrimitive) -> tuple[object, ...]:
        if isinstance(value, ReferenceLine):
            return ("line", point(value.start), point(value.end))
        return ("arc", point(value.start), point(value.end), point(value.centre), float(value.sweep).hex())

    return (tuple(primitive(value) for value in reconstruction.primitives), float(reconstruction.deviation_upper_bound).hex())


class _RecordingReconstructor:
    def __init__(self) -> None:
        self.results: list[ReferenceReconstruction] = []

    def __call__(self, sources: Sequence[SourcePrimitive], transform: SourceToWorld, limit: Millimetre) -> ReferenceReconstruction:
        result = reconstruct_source_path(sources, transform, limit)
        self.results.append(result)
        return result


def _refuse_reconstruction(sources: Sequence[SourcePrimitive], transform: SourceToWorld, limit: Millimetre) -> ReferenceReconstruction:
    raise AssertionError("a cache hit must not recompute the reconstruction")


@pytest.fixture(scope="module")
def figure5_fresh() -> ReferenceReconstruction:
    recorder = _RecordingReconstructor()
    cases_module.validate_case_payload(json.loads(FIGURE5_DOCUMENT), expected_name="figure5", reconstruct=recorder)
    assert len(recorder.results) == 1
    return recorder.results[0]


@pytest.fixture(scope="module")
def figure5_cache(tmp_path_factory: pytest.TempPathFactory) -> ReconstructionCache:
    cache = ReconstructionCache.build(tmp_path_factory.mktemp("held-case-cache"))
    cases_module.decode_case_document(FIGURE5_DOCUMENT, expected_name="figure5", reconstruction_cache=cache)
    return cache


def _figure5_identity() -> ReconstructionBuildIdentity:
    return ReconstructionBuildIdentity.for_document(FIGURE5_DOCUMENT)


def _copy_cache(source: ReconstructionCache, tmp_path: Path) -> ReconstructionCache:
    target = ReconstructionCache.build(tmp_path / "cache")
    target.directory.mkdir()
    for entry in source.directory.iterdir():
        (target.directory / entry.name).write_bytes(entry.read_bytes())
    return target


def test_cache_hit_is_exactly_the_fresh_figure5_reconstruction(
    figure5_cache: ReconstructionCache,
    figure5_fresh: ReferenceReconstruction,
) -> None:
    cached = figure5_cache.load(_figure5_identity())

    assert cached is not None
    assert cached == figure5_fresh
    assert _float_fingerprint(cached) == _float_fingerprint(figure5_fresh)


def test_warm_load_uses_the_entry_and_yields_the_same_case(
    figure5_cache: ReconstructionCache,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fresh_case = cases_module.decode_case_document(FIGURE5_DOCUMENT, expected_name="figure5", reconstruction_cache=None)
    monkeypatch.setattr(cases_module, "reconstruct_source_path", _refuse_reconstruction)

    warm_case = cases_module.decode_case_document(FIGURE5_DOCUMENT, expected_name="figure5", reconstruction_cache=figure5_cache)

    assert warm_case == fresh_case


def test_load_held_reference_case_threads_the_cache(figure5_cache: ReconstructionCache, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(cases_module, "reconstruct_source_path", _refuse_reconstruction)

    case = load_held_reference_case("figure5", reconstruction_cache=figure5_cache)

    assert case.name == "figure5"


def test_miss_computes_and_stores_one_entry(tmp_path: Path) -> None:
    cache = ReconstructionCache.build(tmp_path / "cache")

    cases_module.decode_case_document(FIGURE5_DOCUMENT, expected_name="figure5", reconstruction_cache=cache)

    assert [entry.name for entry in cache.directory.iterdir()] == [f"{_figure5_identity().key}.json"]


def test_a_document_that_fails_validation_is_never_stored(tmp_path: Path) -> None:
    payload = json.loads(FIGURE5_DOCUMENT)
    payload["machining"]["tea_cap_deg"] = 79.0
    cache = ReconstructionCache.build(tmp_path / "cache")

    with pytest.raises(MalformedHeldReferenceCaseError):
        cases_module.decode_case_document(json.dumps(payload).encode(), expected_name="figure5", reconstruction_cache=cache)

    assert not cache.directory.exists() or not any(cache.directory.iterdir())


def test_key_changes_when_document_bytes_change() -> None:
    original = ReconstructionBuildIdentity.for_document(FIGURE5_DOCUMENT)
    whitespace = ReconstructionBuildIdentity.for_document(FIGURE5_DOCUMENT + b"\n")
    edited = ReconstructionBuildIdentity.for_document(FIGURE5_DOCUMENT.replace(b"figure5", b"figure6", 1))

    assert len({original.key, whitespace.key, edited.key}) == 3
    assert original.key == ReconstructionBuildIdentity.for_document(bytes(FIGURE5_DOCUMENT)).key


def test_key_changes_when_any_version_component_changes() -> None:
    original = _figure5_identity()
    other_code = ReconstructionBuildIdentity.build(
        document_digest=original.document_digest,
        code_digest=ReconstructionCodeDigest(bytes(32)),
        interpreter=original.interpreter,
        platform=original.platform,
    )
    other_interpreter = ReconstructionBuildIdentity.build(
        document_digest=original.document_digest,
        code_digest=original.code_digest,
        interpreter=original.interpreter + "+other",
        platform=original.platform,
    )
    other_platform = ReconstructionBuildIdentity.build(
        document_digest=original.document_digest,
        code_digest=original.code_digest,
        interpreter=original.interpreter,
        platform=original.platform + "+other",
    )

    assert len({original.key, other_code.key, other_interpreter.key, other_platform.key}) == 4


def test_identity_factory_rejects_a_non_digest() -> None:
    with pytest.raises(InvalidReconstructionBuildIdentityError):
        ReconstructionBuildIdentity.build(
            document_digest=CaseDocumentDigest(b"short"),
            code_digest=reconstruction_code_digest(),
            interpreter="cpython",
            platform="any",
        )


def test_code_digest_covers_every_first_party_module_that_decides_the_result() -> None:
    covered = {source.module for source in DECIDING_SOURCES}
    uncovered: list[str] = []
    for source in DECIDING_SOURCES:
        tree = ast.parse(source.path.read_text())
        for node in ast.walk(tree):
            names: list[str] = []
            if isinstance(node, ast.ImportFrom) and node.module is not None and node.level == 0:
                names = [node.module]
            elif isinstance(node, ast.Import):
                names = [alias.name for alias in node.names]
            for name in names:
                top = name.split(".")[0]
                if top not in FIRST_PARTY_PACKAGES:
                    continue
                if name not in covered and name not in NON_DECIDING_FIRST_PARTY_IMPORTS:
                    uncovered.append(f"{source.module} imports {name}")

    assert uncovered == []


def test_reconstruction_kernels_import_only_the_standard_library_outside_first_party() -> None:
    kernels = [source for source in DECIDING_SOURCES if source.module != "benchmarks.held_reference_cases"]
    foreign: list[str] = []
    for source in kernels:
        for node in ast.walk(ast.parse(source.path.read_text())):
            names: list[str] = []
            if isinstance(node, ast.ImportFrom) and node.module is not None and node.level == 0:
                names = [node.module]
            elif isinstance(node, ast.Import):
                names = [alias.name for alias in node.names]
            for name in names:
                top = name.split(".")[0]
                if top not in FIRST_PARTY_PACKAGES and top not in sys.stdlib_module_names and top != "__future__":
                    foreign.append(f"{source.module} imports {name}")

    assert foreign == []


def test_an_entry_under_another_identity_raises(figure5_cache: ReconstructionCache, tmp_path: Path) -> None:
    cache = _copy_cache(figure5_cache, tmp_path)
    original = _figure5_identity()
    stale_code = ReconstructionBuildIdentity.build(
        document_digest=original.document_digest,
        code_digest=ReconstructionCodeDigest(bytes(32)),
        interpreter=original.interpreter,
        platform=original.platform,
    )
    cache.entry_path(original).replace(cache.entry_path(stale_code))

    with pytest.raises(ReconstructionCacheIdentityMismatchError):
        cache.load(stale_code)


def test_an_entry_in_another_format_raises(figure5_cache: ReconstructionCache, tmp_path: Path) -> None:
    cache = _copy_cache(figure5_cache, tmp_path)
    path = cache.entry_path(_figure5_identity())
    entry = json.loads(path.read_bytes())
    entry["format"] = "held-reference-reconstruction-cache/v0"
    path.write_text(json.dumps(entry))

    with pytest.raises(UnsupportedReconstructionCacheFormatError):
        cache.load(_figure5_identity())


def test_a_truncated_entry_raises(figure5_cache: ReconstructionCache, tmp_path: Path) -> None:
    cache = _copy_cache(figure5_cache, tmp_path)
    path = cache.entry_path(_figure5_identity())
    path.write_bytes(path.read_bytes()[:-17])

    with pytest.raises(CorruptedReconstructionCacheEntryError):
        cache.load(_figure5_identity())


def test_a_tampered_but_still_decodable_value_raises(figure5_cache: ReconstructionCache, tmp_path: Path) -> None:
    cache = _copy_cache(figure5_cache, tmp_path)
    path = cache.entry_path(_figure5_identity())
    entry = json.loads(path.read_bytes())
    bound = float.fromhex(entry["reconstruction"]["deviation_upper_bound"])
    entry["reconstruction"]["deviation_upper_bound"] = math.nextafter(bound, 0.0).hex()
    path.write_text(json.dumps(entry))
    decode_reconstruction(entry["reconstruction"])

    with pytest.raises(CorruptedReconstructionCacheEntryError):
        cache.load(_figure5_identity())


def test_a_corrupted_entry_is_not_reported_as_a_malformed_case(figure5_cache: ReconstructionCache, tmp_path: Path) -> None:
    cache = _copy_cache(figure5_cache, tmp_path)
    cache.entry_path(_figure5_identity()).write_bytes(b"{")

    with pytest.raises(CorruptedReconstructionCacheEntryError):
        cases_module.decode_case_document(FIGURE5_DOCUMENT, expected_name="figure5", reconstruction_cache=cache)


def test_the_default_cache_lives_in_the_gitignored_build_tree() -> None:
    repository = Path(cases_module.__file__).resolve().parents[1]
    directory = cases_module.DEFAULT_RECONSTRUCTION_CACHE.directory

    assert directory.is_relative_to(repository / "build")
    assert directory == DEFAULT_CACHE_DIRECTORY


def _stale_code_identity() -> ReconstructionBuildIdentity:
    original = _figure5_identity()
    return ReconstructionBuildIdentity.build(
        document_digest=original.document_digest,
        code_digest=ReconstructionCodeDigest(bytes(32)),
        interpreter=original.interpreter,
        platform=original.platform,
    )


def _populated_cache(figure5_cache: ReconstructionCache, tmp_path: Path) -> tuple[ReconstructionCache, set[str]]:
    """A tmp cache with the current figure5 entry, a stale-code entry and an old-format entry."""
    cache = _copy_cache(figure5_cache, tmp_path)
    current = cache.entry_path(_figure5_identity())
    reconstruction = cache.load(_figure5_identity())
    assert reconstruction is not None
    stale_code = cache.store(_stale_code_identity(), reconstruction)
    old_format = cache.directory / f"{'0' * 64}.json"
    entry = json.loads(current.read_bytes())
    entry["format"] = "held-reference-reconstruction-cache/v0"
    old_format.write_text(json.dumps(entry))
    return cache, {stale_code.name, old_format.name}


def test_prune_removes_entries_of_other_code_or_format_and_keeps_current(figure5_cache: ReconstructionCache, tmp_path: Path) -> None:
    cache, stale = _populated_cache(figure5_cache, tmp_path)

    report = prune_stale_entries(cache.directory)

    assert report == ReconstructionCachePruneReport(directory=cache.directory, kept=1, removed=2)
    assert [entry.name for entry in cache.directory.iterdir()] == [cache.entry_path(_figure5_identity()).name]
    assert not any((cache.directory / name).exists() for name in stale)


def test_prune_is_idempotent(figure5_cache: ReconstructionCache, tmp_path: Path) -> None:
    cache, _ = _populated_cache(figure5_cache, tmp_path)
    prune_stale_entries(cache.directory)

    assert prune_stale_entries(cache.directory) == ReconstructionCachePruneReport(directory=cache.directory, kept=1, removed=0)


def test_a_load_never_prunes(figure5_cache: ReconstructionCache, tmp_path: Path) -> None:
    cache, stale = _populated_cache(figure5_cache, tmp_path)

    cases_module.decode_case_document(FIGURE5_DOCUMENT, expected_name="figure5", reconstruction_cache=cache)

    assert all((cache.directory / name).exists() for name in stale)


def test_prune_leaves_in_flight_temporaries_alone(figure5_cache: ReconstructionCache, tmp_path: Path) -> None:
    cache = _copy_cache(figure5_cache, tmp_path)
    temporary = cache.directory / f".{_figure5_identity().key}.abc.tmp"
    temporary.write_bytes(b"{")

    assert prune_stale_entries(cache.directory) == ReconstructionCachePruneReport(directory=cache.directory, kept=1, removed=0)
    assert temporary.exists()


def test_prune_of_a_missing_directory_reports_nothing(tmp_path: Path) -> None:
    directory = tmp_path / "never-created"

    assert prune_stale_entries(directory) == ReconstructionCachePruneReport(directory=directory.resolve(), kept=0, removed=0)


def test_prune_raises_on_a_corrupted_entry_and_removes_nothing(figure5_cache: ReconstructionCache, tmp_path: Path) -> None:
    cache, stale = _populated_cache(figure5_cache, tmp_path)
    (cache.directory / f"{'f' * 64}.json").write_bytes(b"{")
    before = sorted(entry.name for entry in cache.directory.iterdir())

    with pytest.raises(CorruptedReconstructionCacheEntryError):
        prune_stale_entries(cache.directory)

    assert sorted(entry.name for entry in cache.directory.iterdir()) == before


def test_prune_raises_on_a_misfiled_current_format_entry(figure5_cache: ReconstructionCache, tmp_path: Path) -> None:
    cache = _copy_cache(figure5_cache, tmp_path)
    cache.entry_path(_figure5_identity()).replace(cache.directory / f"{'e' * 64}.json")

    with pytest.raises(ReconstructionCacheIdentityMismatchError):
        prune_stale_entries(cache.directory)

    assert [entry.name for entry in cache.directory.iterdir()] == [f"{'e' * 64}.json"]


def test_prune_cli_reports_the_outcome(figure5_cache: ReconstructionCache, tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    cache, _ = _populated_cache(figure5_cache, tmp_path)

    status = cache_main(["prune", "--cache-dir", str(cache.directory)])

    assert status == 0
    assert capsys.readouterr().out.strip() == f"{cache.directory}: kept 1, removed 2"


def test_prune_cli_runs_as_a_module(figure5_cache: ReconstructionCache, tmp_path: Path) -> None:
    cache, _ = _populated_cache(figure5_cache, tmp_path)
    repository = Path(cases_module.__file__).resolve().parents[1]
    environment = {**os.environ, "PYTHONPATH": os.pathsep.join([str(repository), os.environ.get("PYTHONPATH", "")])}

    completed = subprocess.run(
        [sys.executable, "-m", "benchmarks.held_reference_reconstruction_cache", "prune", "--cache-dir", str(cache.directory)],
        cwd=repository,
        env=environment,
        check=True,
        capture_output=True,
        text=True,
    )

    assert completed.stdout.strip() == f"{cache.directory}: kept 1, removed 2"


def test_code_digest_is_the_snapshot_of_the_deciding_modules_as_imported() -> None:
    assert reconstruction_code_digest() == code_digest_of(DECIDING_SOURCES)
    for source in DECIDING_SOURCES:
        module_file = sys.modules[source.module].__file__
        assert module_file is not None
        assert source.path == Path(module_file).resolve()


def test_the_identity_never_rereads_deciding_sources_after_import(monkeypatch: pytest.MonkeyPatch) -> None:
    before = _figure5_identity()

    def _refuse(self: Path) -> bytes:
        raise AssertionError(f"read {self} after import")

    monkeypatch.setattr(Path, "read_bytes", _refuse)

    assert ReconstructionBuildIdentity.for_document(FIGURE5_DOCUMENT) == before


def test_an_edit_after_import_does_not_relabel_the_running_code(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    deciding = Path(sys.modules["benchmarks.held_reference_reconstruction_codec"].__file__ or "")
    name = f"scratch_codec_{uuid.uuid4().hex}"
    scratch = tmp_path / f"{name}.py"
    original = deciding.read_bytes()
    scratch.write_bytes(original)
    monkeypatch.syspath_prepend(str(tmp_path))
    monkeypatch.setattr(sys, "dont_write_bytecode", True)
    module = importlib.import_module(name)
    monkeypatch.setitem(sys.modules, name, module)  # removed again at teardown
    imported = DecidingSource.of_imported(module)

    edited = original + b"\nSCRATCH_EDIT = 1\n"
    scratch.write_bytes(edited)

    assert not hasattr(module, "SCRATCH_EDIT")
    assert imported.source_digest == hashlib.sha256(original).digest()
    assert code_digest_of((imported,)) != code_digest_of((DecidingSource(name, scratch.resolve(), hashlib.sha256(edited).digest()),))

    importlib.reload(module)
    reloaded = DecidingSource.of_imported(module)

    assert module.SCRATCH_EDIT == 1
    assert reloaded.source_digest == hashlib.sha256(edited).digest()
    assert code_digest_of((reloaded,)) != code_digest_of((imported,))


_FINITE = st.floats(allow_nan=False, allow_infinity=False, width=64)
_BOUNDED = st.floats(min_value=-1.0e6, max_value=1.0e6, allow_nan=False, allow_infinity=False, width=64)
_RADIUS = st.floats(min_value=1.0e-3, max_value=1.0e3, allow_nan=False, allow_infinity=False, width=64)
_SWEEP = st.floats(min_value=-math.tau, max_value=math.tau, allow_nan=False, allow_infinity=False, width=64).filter(lambda value: value != 0.0)


@st.composite
def _reconstructions(draw: st.DrawFn) -> ReferenceReconstruction:
    cursor = Point2[WorldXY].build(draw(_FINITE), draw(_FINITE))
    primitives: list[ReferencePrimitive] = []
    for _ in range(draw(st.integers(min_value=1, max_value=6))):
        if draw(st.booleans()):
            end = Point2[WorldXY].build(draw(_FINITE), draw(_FINITE))
            assume(end != cursor)
            primitives.append(ReferenceLine.build(cursor, end))
            cursor = end
            continue
        start = Point2[WorldXY].build(draw(_BOUNDED), draw(_BOUNDED))
        if primitives:
            start = cursor
        assume(abs(float(start.x)) <= 1.0e6 and abs(float(start.y)) <= 1.0e6)
        radius = draw(_RADIUS)
        bearing = draw(st.floats(min_value=-math.pi, max_value=math.pi, allow_nan=False, width=64))
        sweep = draw(_SWEEP)
        centre = Point2[WorldXY].build(float(start.x) - radius * math.cos(bearing), float(start.y) - radius * math.sin(bearing))
        end = Point2[WorldXY].build(float(centre.x) + radius * math.cos(bearing + sweep), float(centre.y) + radius * math.sin(bearing + sweep))
        try:
            primitives.append(ReferenceArc.build(start, end, centre, Radian(sweep)))
        except InvalidPublishedPrimitiveError:
            assume(False)
        cursor = end
    return ReferenceReconstruction.build(primitives, Millimetre(abs(draw(_FINITE))))


@settings(max_examples=300, deadline=None)
@given(_reconstructions())
def test_codec_round_trips_any_reconstruction_bitwise(reconstruction: ReferenceReconstruction) -> None:
    encoded = encode_reconstruction(reconstruction)
    decoded = decode_reconstruction(json.loads(json.dumps(encoded)))

    assert decoded == reconstruction
    assert _float_fingerprint(decoded) == _float_fingerprint(reconstruction)
    assert encode_reconstruction(decoded) == encoded


@settings(max_examples=100, deadline=None)
@given(st.sampled_from([-0.0, 0.0, 5e-324, -5e-324, sys.float_info.max, -sys.float_info.max, 0.1]))
def test_codec_preserves_signed_zero_and_extremes(value: float) -> None:
    other = 1.0 if value != 1.0 else 2.0
    reconstruction = ReferenceReconstruction.build(
        (ReferenceLine.build(Point2[WorldXY].build(value, value), Point2[WorldXY].build(other, other)),),
        Millimetre(0.0),
    )

    decoded = decode_reconstruction(encode_reconstruction(reconstruction))

    assert _float_fingerprint(decoded) == _float_fingerprint(reconstruction)


@pytest.mark.parametrize(
    "encoded",
    [
        {},
        {"primitives": [], "deviation_upper_bound": "0x0.0p+0"},
        {"primitives": [{"kind": "line", "start": ["0x1.0p+0", "0x0.0p+0"], "end": ["inf", "0x0.0p+0"]}], "deviation_upper_bound": "0x0.0p+0"},
        {"primitives": [{"kind": "spline", "start": ["0x1.0p+0", "0x0.0p+0"]}], "deviation_upper_bound": "0x0.0p+0"},
        {"primitives": [{"kind": "line", "start": [1.0, 0.0], "end": [2.0, 0.0]}], "deviation_upper_bound": "0x0.0p+0"},
    ],
)
def test_codec_rejects_malformed_encodings(encoded: dict[str, object]) -> None:
    with pytest.raises(MalformedReconstructionEncodingError):
        decode_reconstruction(encoded)
