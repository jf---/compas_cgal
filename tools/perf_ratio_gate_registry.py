"""The ratio gate registry: ``benchmarks/results/perf_ratio_gate.json``.

The registry holds the reference commit the gate decides against, the harness
digest every gated session must carry, the gated workload set, and per CPU model
one A/A calibration: the margins ``m_w`` and the cited ``perf-ratio/1`` records
they come from (spec §Margins: ``m_w = SAFETY * max_k |rho_{w,k} - 1|``).

The registry also stores the anchor chain and refuses any stored anchor ``A_w``
that disagrees with the chain's replay (`tools.perf_ratio_gate_chain`).

The registry is content-addressed. Its ``registry_sha256`` is the SHA-256 of its
canonical JSON content, and each cited record carries its file's SHA-256. On load
nothing stored is trusted: the self-digest is re-derived, every cited record is
re-hashed, re-read and re-validated as an honest A/A session of this registry,
and every margin is re-derived from those records and compared with ``==``
(the same float operations on the same inputs give the same bits, and JSON
round-trips a float exactly, so no tolerance applies). `write_registry` runs
the same verification before it writes, so it never writes a registry that
would not load. A ``perf-ratio-gate/1`` registry (no anchor chain) is refused
with its migration: rebuild it with ``perf-ratio-calibrate --records``.

The registry is the top of three layers, imports pointing downward only:
`tools.perf_ratio_gate_citations` (cited A/A records, margins, calibration
verification), `tools.perf_ratio_gate_chain` (anchor links, the reset rule, the
replay), and this module (the registry file: build, load, write, its schema and
self-digest). The names of the lower layers are re-exported here.
"""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Mapping
from collections.abc import Sequence
from dataclasses import dataclass
from dataclasses import field as dataclass_field
from dataclasses import replace
from fractions import Fraction
from pathlib import Path
from typing import Any
from typing import Optional

from tools.perf_history_json import MalformedTimingRecordError
from tools.perf_history_json import field
from tools.perf_history_json import items
from tools.perf_history_json import number
from tools.perf_history_json import obj
from tools.perf_history_json import text
from tools.perf_history_store import REPOSITORY
from tools.perf_ratio_gate_chain import AnchorChainBrokenError as AnchorChainBrokenError
from tools.perf_ratio_gate_chain import AnchorChainEndError as AnchorChainEndError
from tools.perf_ratio_gate_chain import AnchorChainMismatchError as AnchorChainMismatchError
from tools.perf_ratio_gate_chain import AnchorLink as AnchorLink
from tools.perf_ratio_gate_chain import AnchorLinkSelfLoopError as AnchorLinkSelfLoopError
from tools.perf_ratio_gate_chain import AnchorRecordDigestError as AnchorRecordDigestError
from tools.perf_ratio_gate_chain import AnchorRecordMissingError as AnchorRecordMissingError
from tools.perf_ratio_gate_chain import AnchorRecordOutsideRootError as AnchorRecordOutsideRootError
from tools.perf_ratio_gate_chain import AnchorRecordUnloadableError as AnchorRecordUnloadableError
from tools.perf_ratio_gate_chain import AnchorValueMismatchError as AnchorValueMismatchError
from tools.perf_ratio_gate_chain import DuplicateAnchorRecordError as DuplicateAnchorRecordError
from tools.perf_ratio_gate_chain import ExtensionPairCapError as ExtensionPairCapError
from tools.perf_ratio_gate_chain import InvalidAnchorSessionError as InvalidAnchorSessionError
from tools.perf_ratio_gate_chain import RegressedAnchorLinkError as RegressedAnchorLinkError
from tools.perf_ratio_gate_chain import UnimprovedAnchorLinkError as UnimprovedAnchorLinkError
from tools.perf_ratio_gate_chain import UnrequestedExtensionError as UnrequestedExtensionError
from tools.perf_ratio_gate_chain import carry_anchor as carry_anchor
from tools.perf_ratio_gate_chain import derive_anchors as derive_anchors
from tools.perf_ratio_gate_chain import link_ratios as link_ratios
from tools.perf_ratio_gate_citations import CALIBRATION_SESSIONS as CALIBRATION_SESSIONS
from tools.perf_ratio_gate_citations import SAFETY as SAFETY
from tools.perf_ratio_gate_citations import Calibration as Calibration
from tools.perf_ratio_gate_citations import CalibrationRecordDigestError as CalibrationRecordDigestError
from tools.perf_ratio_gate_citations import CalibrationRecordMissingError as CalibrationRecordMissingError
from tools.perf_ratio_gate_citations import CalibrationRecordOutsideRootError as CalibrationRecordOutsideRootError
from tools.perf_ratio_gate_citations import CalibrationRecordUnloadableError as CalibrationRecordUnloadableError
from tools.perf_ratio_gate_citations import CitedRecord as CitedRecord
from tools.perf_ratio_gate_citations import DuplicateCalibrationSessionError as DuplicateCalibrationSessionError
from tools.perf_ratio_gate_citations import InvalidCalibrationSessionError as InvalidCalibrationSessionError
from tools.perf_ratio_gate_citations import MarginMismatchError as MarginMismatchError
from tools.perf_ratio_gate_citations import TooFewCalibrationSessionsError as TooFewCalibrationSessionsError
from tools.perf_ratio_gate_citations import _verify_calibration
from tools.perf_ratio_gate_citations import calibration_of as calibration_of
from tools.perf_ratio_gate_citations import derive_margins as derive_margins
from tools.perf_ratio_gate_decide import UNIT_ANCHOR as UNIT_ANCHOR
from tools.perf_ratio_gate_decide import Anchor as Anchor
from tools.perf_ratio_gate_decide import Margin as Margin
from tools.perf_ratio_gate_errors import GateRefusal
from tools.perf_ratio_record import REF_SHORT_CHARS

REGISTRY_SCHEMA = "perf-ratio-gate/2"
# Registry schemas this gate no longer reads, each refused with its migration.
SUPERSEDED_SCHEMAS = frozenset({"perf-ratio-gate/1"})
# The calibrate CLI's pixi task; the one source of the command every refusal names.
CALIBRATE_TASK = "pixi run perf-ratio-calibrate"
# How a superseded registry is rebuilt: over its A/A records, at REGISTRY_SCHEMA with an empty chain.
MIGRATION_COMMAND = f"{CALIBRATE_TASK} --records"
REGISTRY_PATH = REPOSITORY / "benchmarks" / "results" / "perf_ratio_gate.json"
CALIBRATION_ROOT = REPOSITORY / "benchmarks" / "results" / "perf_ratio_calibration"
PROMOTION_ROOT = REPOSITORY / "benchmarks" / "results" / "perf_ratio_promotions"
# Runs of characters a CPU model slug replaces with one hyphen.
_SLUG_SEPARATORS = re.compile(r"[^a-z0-9]+")


class MalformedRegistryError(GateRefusal):
    """The registry file is not JSON, or has the wrong schema, keys or field types."""


class SupersededRegistrySchemaError(MalformedRegistryError):
    """The registry has a schema this gate no longer reads (`SUPERSEDED_SCHEMAS`); its message names the migration."""


class RegistryDigestMismatchError(GateRefusal):
    """The stored ``registry_sha256`` is not the digest of the registry's content."""


class SafetyMismatchError(GateRefusal):
    """The registry's ``safety`` is not this module's `SAFETY`."""


class DuplicateCalibrationError(GateRefusal):
    """A registry holds two calibrations for one ``cpu_model``."""


class CalibrationWorkloadMismatchError(GateRefusal):
    """A calibration's workload set is not the registry's."""


class MissingWorkloadsError(GateRefusal):
    """A registry's workload set is empty (a calibration-less registry is built without one)."""


class DuplicateRegistryWorkloadError(GateRefusal):
    """A registry's workload set names one workload more than once."""


class UnsortedRegistryWorkloadsError(GateRefusal):
    """A registry's workload set is not in sorted order."""


class UncalibratedMachineError(GateRefusal):
    """The registry holds no calibration for a ``cpu_model``."""


class UnsluggableCpuModelError(GateRefusal):
    """A CPU model has no letter or digit, so it has no calibration directory name."""


@dataclass(frozen=True)
class GateRegistry:
    """The gate's reference, harness, workload set and per-machine margins.

    Raw construction is unchecked; `build_registry` and `load_registry` build checked ones.

    Attributes:
        ref_commit: Full commit id of the reference the gate decides against.
        harness_sha256: The harness digest every calibration and gated session carries.
        safety: The margin factor the margins were derived with (`SAFETY`).
        workloads: The gated workload names, sorted.
        calibrations: One calibration per CPU model, sorted by ``cpu_model``.
        anchor_ref: The commit the anchor chain starts at; None (raw construction only): ``ref_commit``.
        anchor_links: The promotions since ``anchor_ref``, in order.
        exact_anchors: Workload -> its exact ``A_w``, as `derive_anchors` gives it from the chain.
    """

    ref_commit: str
    harness_sha256: str
    safety: float
    workloads: tuple[str, ...]
    calibrations: tuple[Calibration, ...]
    anchor_ref: Optional[str] = None
    anchor_links: tuple[AnchorLink, ...] = ()
    exact_anchors: Mapping[str, Anchor] = dataclass_field(default_factory=dict)

    @property
    def chain_start(self) -> str:
        """The commit the anchor chain starts at: ``anchor_ref``, or ``ref_commit`` when it is None."""
        return self.ref_commit if self.anchor_ref is None else self.anchor_ref

    def anchors_exact(self) -> Mapping[str, Anchor]:
        """Workload -> its exact anchor ``A_w`` (the decision value)."""
        return self.exact_anchors

    @property
    def anchors(self) -> dict[str, float]:
        """Workload -> its anchor as the nearest double: a reporting copy, never decided on."""
        return {name: float(anchor) for name, anchor in self.exact_anchors.items()}

    def margins_for(self, cpu_model: str) -> Mapping[str, Margin]:
        """The margins calibrated on *cpu_model*.

        Raises:
            UncalibratedMachineError: No calibration names *cpu_model*.
        """
        for calibration in self.calibrations:
            if calibration.cpu_model == cpu_model:
                return calibration.margins
        calibrated = ", ".join(repr(calibration.cpu_model) for calibration in self.calibrations) or "none"
        raise UncalibratedMachineError(f"no A/A calibration for cpu_model {cpu_model!r} (calibrated: {calibrated}); run an A/A calibration on that machine first.")

    def content(self) -> dict[str, object]:
        """The registry's JSON content, without its digest."""
        return {
            "schema": REGISTRY_SCHEMA,
            "ref_commit": self.ref_commit,
            "harness_sha256": self.harness_sha256,
            "safety": self.safety,
            "workloads": list(self.workloads),
            "calibrations": [calibration.to_json() for calibration in self.calibrations],
            "anchor_ref": self.chain_start,
            "anchor_links": [link.to_json() for link in self.anchor_links],
            "anchors": self.anchors,
        }

    @property
    def registry_sha256(self) -> str:
        """SHA-256 of the canonical JSON content (sorted keys, ``(",", ":")`` separators)."""
        return _digest(self.content())

    def dumps(self) -> str:
        """The registry as deterministic JSON text (sorted keys, two-space indent, trailing newline), digest included."""
        document = self.content()
        document["registry_sha256"] = self.registry_sha256
        return json.dumps(document, indent=2, sort_keys=True) + "\n"


def calibrated_slug(cpu_model: str) -> str:
    """The directory name of *cpu_model*: lowercase, runs of non-``[a-z0-9]`` as one ``-``, stripped of ``-``.

    Raises:
        UnsluggableCpuModelError: *cpu_model* has no ASCII letter or digit.
    """
    slug = _SLUG_SEPARATORS.sub("-", cpu_model.lower()).strip("-")
    if not slug:
        raise UnsluggableCpuModelError(f"cpu_model {cpu_model!r} has no letter or digit to name a calibration directory.")
    return slug


def calibration_dir(ref_commit: str, cpu_model: str) -> Path:
    """Where the A/A records of *ref_commit* on *cpu_model* live: ``CALIBRATION_ROOT/<ref short>/<cpu slug>``."""
    return CALIBRATION_ROOT / ref_commit[:REF_SHORT_CHARS] / calibrated_slug(cpu_model)


def promotion_dir(from_ref: str, to_ref: str) -> Path:
    """Where the records of promoting *from_ref* to *to_ref* live: ``PROMOTION_ROOT/<from short>-<to short>``."""
    return PROMOTION_ROOT / f"{from_ref[:REF_SHORT_CHARS]}-{to_ref[:REF_SHORT_CHARS]}"


def build_registry(
    *,
    ref_commit: str,
    harness_sha256: str,
    calibrations: Sequence[Calibration],
    workloads: Optional[Sequence[str]] = None,
    anchor_ref: Optional[str] = None,
    anchor_links: Sequence[AnchorLink] = (),
    root: Path = REPOSITORY,
) -> GateRegistry:
    """A registry over *calibrations*, its workload set theirs, its anchors derived from *anchor_links*.

    Args:
        ref_commit: Full commit id of the reference.
        harness_sha256: The harness digest.
        calibrations: At most one calibration per CPU model; may be empty (a freshly promoted reference).
        workloads: The gated workload names; required when *calibrations* is empty, and when given
            beside calibrations it must be their workload set.
        anchor_ref: The commit the chain starts at; None: *ref_commit* (an empty chain).
        anchor_links: The promotions since *anchor_ref*, in order; their records are verified.
        root: The directory the links' cited paths are relative to (`REPOSITORY`).

    Returns:
        The registry, workloads sorted and calibrations sorted by ``cpu_model``, `SAFETY` as its safety,
        its exact anchors derived by `derive_anchors`.

    Raises:
        DuplicateCalibrationError: Two calibrations name one ``cpu_model``.
        MissingWorkloadsError: No calibrations and no (or an empty) workload set.
        DuplicateRegistryWorkloadError: *workloads* names one workload twice.
        CalibrationWorkloadMismatchError: A calibration's workload set is not the registry's.
        AnchorChainBrokenError: A link does not start where the chain stands.
        AnchorChainEndError: The chain does not end at *ref_commit*.
        InvalidAnchorSessionError: As `derive_anchors`.
        AnchorRecordOutsideRootError: As `derive_anchors`.
        AnchorRecordMissingError: As `derive_anchors`.
        AnchorRecordDigestError: As `derive_anchors`.
        ExtensionMismatchError: As `derive_anchors`.
        InvalidCalibrationSessionError: As `derive_anchors`.
        MarginMismatchError: As `derive_anchors`.
    """
    _require_one_calibration_per_cpu(calibrations)
    if workloads is not None:
        names = tuple(sorted(workloads))
    elif calibrations:
        names = tuple(sorted(calibrations[0].margins))
    else:
        names = ()
    _require_workload_shape(names)
    for calibration in calibrations:
        _require_calibration_workloads(calibration, names)
    start = ref_commit if anchor_ref is None else anchor_ref
    _require_chain_continuity(start, anchor_links, ref_commit)
    return GateRegistry(
        ref_commit=ref_commit,
        harness_sha256=harness_sha256,
        safety=SAFETY,
        workloads=names,
        calibrations=tuple(sorted(calibrations, key=lambda calibration: calibration.cpu_model)),
        anchor_ref=start,
        anchor_links=tuple(anchor_links),
        exact_anchors=derive_anchors(anchor_links, names, root=root),
    )


def load_registry(path: Path = REGISTRY_PATH, *, root: Path = REPOSITORY) -> GateRegistry:
    """The registry at *path*, every stored value re-derived.

    Checks, in order: schema and shape; the self-digest; safety; the workload set (non-empty, sorted,
    no repeats); one calibration per CPU; then per calibration: K, each cited record's digest, each session's validity, the margin re-derivation;
    then the anchor chain: continuity from ``anchor_ref`` to ``ref_commit``, each link's records and calibration,
    and the stored anchors against the re-derived ones. The exact anchors returned are the re-derived ones,
    never the stored reporting floats.

    Args:
        path: The registry file.
        root: The directory cited record paths are relative to (`REPOSITORY`).

    Returns:
        The verified registry.

    Raises:
        MalformedRegistryError: Not UTF-8, not JSON, or the wrong schema, keys or field types.
        SupersededRegistrySchemaError: A ``perf-ratio-gate/1`` registry; the message names the migration.
        RegistryDigestMismatchError: The stored ``registry_sha256`` is not the content's.
        SafetyMismatchError: The stored safety is not `SAFETY`.
        MissingWorkloadsError: The workload set is empty.
        UnsortedRegistryWorkloadsError: The workload set is not sorted.
        DuplicateRegistryWorkloadError: The workload set names one workload twice.
        DuplicateCalibrationError: Two calibrations name one ``cpu_model``.
        TooFewCalibrationSessionsError: A calibration cites fewer than `CALIBRATION_SESSIONS` records.
        DuplicateCalibrationSessionError: A calibration cites one record twice.
        CalibrationRecordOutsideRootError: A cited path leaves *root*.
        CalibrationRecordMissingError: A cited record does not exist.
        CalibrationRecordDigestError: A cited record's bytes changed.
        CalibrationRecordUnloadableError: A cited record has its digest but does not load (a harness
            parser change, or a hand edit).
        InvalidCalibrationSessionError: A cited session is not an honest A/A session of this registry.
        MarginMismatchError: A stored margin set differs from the re-derived one.
        AnchorChainBrokenError: A link does not start where the chain stands.
        AnchorChainEndError: The chain does not end at the reference.
        InvalidAnchorSessionError: As `derive_anchors`.
        AnchorRecordOutsideRootError: As `derive_anchors`.
        AnchorRecordMissingError: As `derive_anchors`.
        AnchorRecordDigestError: As `derive_anchors`.
        AnchorRecordUnloadableError: As `derive_anchors`.
        ExtensionMismatchError: As `derive_anchors`.
        AnchorValueMismatchError: The stored anchors are not the re-derived ones.
    """
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except UnicodeDecodeError as error:
        raise MalformedRegistryError(f"{path} is not UTF-8: {error}") from error
    except json.JSONDecodeError as error:
        raise MalformedRegistryError(f"{path} is not JSON: {error}") from error
    try:
        registry, stored_digest = _parse(raw, str(path))
    except MalformedTimingRecordError as error:
        raise MalformedRegistryError(str(error)) from error
    if stored_digest != registry.registry_sha256:
        raise RegistryDigestMismatchError(f"{path}: stored registry_sha256 {stored_digest} is not its content's {registry.registry_sha256}.")
    return replace(registry, exact_anchors=_verify(registry, root))


def write_registry(registry: GateRegistry, path: Path = REGISTRY_PATH, *, root: Path = REPOSITORY) -> None:
    """Write *registry* to *path* after verifying it as `load_registry` would.

    The text is written to a sibling file and moved over *path*, so *path* is
    never left half-written.

    Args:
        registry: The registry to write.
        path: The registry file.
        root: The directory cited record paths are relative to (`REPOSITORY`).

    Raises:
        SafetyMismatchError: As `load_registry`; *path* is then untouched.
        MissingWorkloadsError: As `load_registry`.
        UnsortedRegistryWorkloadsError: As `load_registry`.
        DuplicateRegistryWorkloadError: As `load_registry`.
        DuplicateCalibrationError: As `load_registry`.
        TooFewCalibrationSessionsError: As `load_registry`.
        DuplicateCalibrationSessionError: As `load_registry`.
        CalibrationRecordOutsideRootError: As `load_registry`.
        CalibrationRecordMissingError: As `load_registry`.
        CalibrationRecordDigestError: As `load_registry`.
        InvalidCalibrationSessionError: As `load_registry`.
        MarginMismatchError: As `load_registry`.
        AnchorChainBrokenError: As `load_registry`.
        AnchorChainEndError: As `load_registry`.
        InvalidAnchorSessionError: As `load_registry`.
        AnchorRecordOutsideRootError: As `load_registry`.
        AnchorRecordMissingError: As `load_registry`.
        AnchorRecordDigestError: As `load_registry`.
        ExtensionMismatchError: As `load_registry`.
        AnchorValueMismatchError: As `load_registry`; also when the exact anchors are not the
            re-derived ones, even where their reporting floats agree.
    """
    derived = _verify(registry, root)
    if dict(registry.exact_anchors) != derived:
        raise AnchorValueMismatchError(f"exact anchors {dict(registry.exact_anchors)} are not the ones the anchor chain re-derives, {derived}.")
    staging = path.with_name(f"{path.name}.writing")
    staging.write_text(registry.dumps(), encoding="utf-8")
    staging.replace(path)


def _verify(registry: GateRegistry, root: Path) -> dict[str, Anchor]:
    # Every check load and write share; returns the anchors the chain re-derives.
    if registry.safety != SAFETY:
        raise SafetyMismatchError(f"registry safety is {registry.safety!r}, not this gate's SAFETY {SAFETY!r}.")
    _require_workload_shape(registry.workloads)
    _require_one_calibration_per_cpu(registry.calibrations)
    workloads = frozenset(registry.workloads)
    for calibration in registry.calibrations:
        _verify_calibration(
            calibration,
            ref_commit=registry.ref_commit,
            harness_sha256=registry.harness_sha256,
            harness_owner="registry's",
            workloads=workloads,
            root=root,
            where=f"calibration {calibration.cpu_model!r}",
        )
    _require_chain_continuity(registry.chain_start, registry.anchor_links, registry.ref_commit)
    derived = derive_anchors(registry.anchor_links, registry.workloads, root=root)
    reported = {name: float(anchor) for name, anchor in derived.items()}
    if registry.anchors != reported:
        raise AnchorValueMismatchError(f"stored anchors {registry.anchors} are not the ones the anchor chain re-derives, {reported}.")
    return derived


def _require_chain_continuity(anchor_ref: str, links: Sequence[AnchorLink], ref_commit: str) -> None:
    standing = anchor_ref
    for index, link in enumerate(links):
        if link.from_ref == link.to_ref:
            raise AnchorLinkSelfLoopError(f"anchor link {index} promotes {link.from_ref} to itself; a promotion moves the reference.")
        if link.from_ref != standing:
            raise AnchorChainBrokenError(f"anchor link {index} starts at {link.from_ref}, not where the chain stands, {standing}.")
        standing = link.to_ref
    if standing != ref_commit:
        raise AnchorChainEndError(f"the anchor chain from {anchor_ref} ends at {standing}, not the registry's reference {ref_commit}.")


def _require_workload_shape(workloads: Sequence[str]) -> None:
    # The one shape check build, write and load share: non-empty, sorted, no repeats.
    if not workloads:
        raise MissingWorkloadsError("a registry names a non-empty workload set; a calibration-less registry is built with workloads.")
    if len(set(workloads)) != len(workloads):
        raise DuplicateRegistryWorkloadError(f"a registry names each workload once, not {list(workloads)}.")
    if list(workloads) != sorted(workloads):
        raise UnsortedRegistryWorkloadsError(f"a registry's workloads are sorted, not {list(workloads)}.")


def _require_one_calibration_per_cpu(calibrations: Sequence[Calibration]) -> None:
    models = [calibration.cpu_model for calibration in calibrations]
    repeated = sorted({model for model in models if models.count(model) > 1})
    if repeated:
        raise DuplicateCalibrationError(f"a registry holds one calibration per cpu_model; more than one for {', '.join(repr(model) for model in repeated)}.")


def _require_calibration_workloads(calibration: Calibration, workloads: Sequence[str]) -> None:
    if frozenset(calibration.margins) != frozenset(workloads):
        raise CalibrationWorkloadMismatchError(f"calibration {calibration.cpu_model!r} covers {sorted(calibration.margins)}, not the registry's workloads {list(workloads)}.")


_REGISTRY_KEYS = frozenset({"schema", "ref_commit", "harness_sha256", "safety", "workloads", "calibrations", "anchor_ref", "anchor_links", "anchors", "registry_sha256"})
_CALIBRATION_KEYS = frozenset({"cpu_model", "sessions", "margins"})
_LINK_KEYS = frozenset({"from_ref", "to_ref", "sessions", "calibration"})
_CITED_KEYS = frozenset({"path", "record_sha256"})


def _parse(raw: object, where: str) -> tuple[GateRegistry, str]:
    record = obj(raw, where)
    schema = text(field(record, "schema", where), f"{where}.schema")
    if schema in SUPERSEDED_SCHEMAS:
        raise SupersededRegistrySchemaError(
            f"{where} is a {schema} registry, without the anchor chain; this gate reads {REGISTRY_SCHEMA}. "
            f"Migrate it: move it aside, then run one {MIGRATION_COMMAND} <A/A records> per machine's A/A records "
            f"(a calibration never mixes cpu_models); this rebuilds it at {REGISTRY_SCHEMA} with an empty chain."
        )
    if schema != REGISTRY_SCHEMA:
        raise MalformedRegistryError(f"{where}.schema is {REGISTRY_SCHEMA!r}, not {schema!r}.")
    _require_keys(record, _REGISTRY_KEYS, where)
    workloads = tuple(text(item, f"{where}.workloads[{index}]") for index, item in enumerate(items(field(record, "workloads", where), f"{where}.workloads")))
    calibration_items = items(field(record, "calibrations", where), f"{where}.calibrations")
    link_items = items(field(record, "anchor_links", where), f"{where}.anchor_links")
    anchors = obj(field(record, "anchors", where), f"{where}.anchors")
    registry = GateRegistry(
        ref_commit=text(field(record, "ref_commit", where), f"{where}.ref_commit"),
        harness_sha256=text(field(record, "harness_sha256", where), f"{where}.harness_sha256"),
        safety=_json_number(field(record, "safety", where), f"{where}.safety"),
        workloads=workloads,
        calibrations=tuple(_parse_calibration(item, f"{where}.calibrations[{index}]") for index, item in enumerate(calibration_items)),
        anchor_ref=text(field(record, "anchor_ref", where), f"{where}.anchor_ref"),
        anchor_links=tuple(_parse_link(item, f"{where}.anchor_links[{index}]") for index, item in enumerate(link_items)),
        # The stored reporting floats, as exact rationals so the digest re-derives; load replaces them with the chain's.
        exact_anchors={name: Anchor(Fraction(_json_number(value, f"{where}.anchors.{name}"))) for name, value in anchors.items()},
    )
    return registry, text(field(record, "registry_sha256", where), f"{where}.registry_sha256")


def _parse_calibration(raw: object, where: str) -> Calibration:
    record = obj(raw, where)
    _require_keys(record, _CALIBRATION_KEYS, where)
    margins = obj(field(record, "margins", where), f"{where}.margins")
    session_items = items(field(record, "sessions", where), f"{where}.sessions")
    return Calibration(
        cpu_model=text(field(record, "cpu_model", where), f"{where}.cpu_model"),
        sessions=tuple(_parse_cited(item, f"{where}.sessions[{index}]") for index, item in enumerate(session_items)),
        margins={name: Margin(_json_number(value, f"{where}.margins.{name}")) for name, value in margins.items()},
    )


def _parse_link(raw: object, where: str) -> AnchorLink:
    record = obj(raw, where)
    _require_keys(record, _LINK_KEYS, where)
    session_items = items(field(record, "sessions", where), f"{where}.sessions")
    return AnchorLink(
        from_ref=text(field(record, "from_ref", where), f"{where}.from_ref"),
        to_ref=text(field(record, "to_ref", where), f"{where}.to_ref"),
        sessions=tuple(_parse_cited(item, f"{where}.sessions[{index}]") for index, item in enumerate(session_items)),
        calibration=_parse_calibration(field(record, "calibration", where), f"{where}.calibration"),
    )


def _parse_cited(raw: object, where: str) -> CitedRecord:
    record = obj(raw, where)
    _require_keys(record, _CITED_KEYS, where)
    return CitedRecord(
        path=text(field(record, "path", where), f"{where}.path"),
        record_sha256=text(field(record, "record_sha256", where), f"{where}.record_sha256"),
    )


def _json_number(value: object, where: str) -> float:
    # The shared accessor also accepts decimal strings; this registry writes only JSON numbers.
    if isinstance(value, str):
        raise MalformedRegistryError(f"{where} is a JSON number, not the string {value!r}.")
    return number(value, where)


def _require_keys(record: Mapping[str, Any], expected: frozenset[str], where: str) -> None:
    keys = frozenset(record)
    if keys != expected:
        raise MalformedRegistryError(f"{where} keys differ: missing {sorted(expected - keys)}, unexpected {sorted(keys - expected)}.")


def _canonical(document: Mapping[str, object]) -> str:
    return json.dumps(document, sort_keys=True, separators=(",", ":"))


def _digest(content: Mapping[str, object]) -> str:
    return hashlib.sha256(_canonical(content).encode("utf-8")).hexdigest()
