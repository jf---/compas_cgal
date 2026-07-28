import hashlib
from dataclasses import dataclass
from typing import Self

import numpy as np

from compas_cgal.adaptive.canonical import ExactRationalV1
from compas_cgal.adaptive.errors import InvalidNeckEvidenceError
from compas_cgal.adaptive.errors import InvalidNeckPassageTransitionError
from compas_cgal.adaptive.errors import TerminalNeckPassageError
from compas_cgal.adaptive.identity import IdentityDigest
from compas_cgal.adaptive.medial_axis import MatEdgeId
from compas_cgal.adaptive.medial_axis import MatSiteId
from compas_cgal.adaptive.medial_axis import MedialAxis
from compas_cgal.adaptive.operation import NeckCapDecision
from compas_cgal.adaptive.operation import NeckOwnerId
from compas_cgal.adaptive.operation import NeckTraversalOrientation
from compas_cgal.adaptive.operation import OrientedNeckScope
from compas_cgal.adaptive.operation import WidthClassId
from compas_cgal.adaptive.policy import NeckPolicy
from compas_cgal.adaptive.policy import PassageState

_NEXT_PASSAGE_STATE = {
    PassageState.UNVISITED: PassageState.FIRST_PASS_COMPLETE,
    PassageState.FIRST_PASS_COMPLETE: PassageState.SECOND_PASS_COMPLETE,
    PassageState.SECOND_PASS_COMPLETE: PassageState.TERMINAL,
}
_SHA256_BYTES = hashlib.sha256().digest_size


def _exact_bytes(value: object, name: str, *, digest: bool = False) -> bytes:
    if type(value) is not bytes or not value:
        raise InvalidNeckEvidenceError(f"{name} must be nonempty exact bytes.")
    if digest and len(value) != _SHA256_BYTES:
        raise InvalidNeckEvidenceError(f"{name} must be one SHA-256 digest.")
    return value


@dataclass(frozen=True)
class ClassifiedNeck:
    owner_id: NeckOwnerId
    defining_site_ids: tuple[MatSiteId, ...]
    evidence_bytes: bytes
    evidence_digest: IdentityDigest
    width_class_id: WidthClassId
    comparison_certificate: bytes
    separating_cut_edge_ids: tuple[MatEdgeId, ...]

    def __post_init__(self) -> None:
        _exact_bytes(self.owner_id, "neck owner identity")
        if type(self.defining_site_ids) is not tuple or not self.defining_site_ids or self.defining_site_ids != tuple(sorted(set(self.defining_site_ids))):
            raise InvalidNeckEvidenceError("neck defining-site identities must be nonempty and canonical.")
        _exact_bytes(self.evidence_bytes, "neck evidence")
        _exact_bytes(self.evidence_digest, "neck evidence digest", digest=True)
        if hashlib.sha256(self.evidence_bytes).digest() != self.evidence_digest:
            raise InvalidNeckEvidenceError("neck evidence digest contradicts its exact bytes.")
        if type(self.width_class_id) is not WidthClassId:
            raise InvalidNeckEvidenceError("neck width class must be one exact typed identity.")
        _exact_bytes(self.comparison_certificate, "neck comparison certificate")
        if bytes(self.evidence_digest) not in self.comparison_certificate:
            raise InvalidNeckEvidenceError("neck comparison certificate omits its evidence digest.")
        if type(self.separating_cut_edge_ids) is not tuple or not self.separating_cut_edge_ids or self.separating_cut_edge_ids != tuple(sorted(set(self.separating_cut_edge_ids))):
            raise InvalidNeckEvidenceError("neck separating cut must contain canonical edge identities.")

    @property
    def forward(self) -> "NeckPassage":
        return NeckPassage.build(
            neck=self,
            orientation=NeckTraversalOrientation.FORWARD,
            state=PassageState.UNVISITED,
        )

    @property
    def reverse(self) -> "NeckPassage":
        return NeckPassage.build(
            neck=self,
            orientation=NeckTraversalOrientation.REVERSE,
            state=PassageState.UNVISITED,
        )


@dataclass(frozen=True)
class NeckPassage:
    neck: ClassifiedNeck
    orientation: NeckTraversalOrientation
    state: PassageState

    def __post_init__(self) -> None:
        if type(self.neck) is not ClassifiedNeck:
            raise InvalidNeckPassageTransitionError("neck passage requires one exact classified neck.")
        if type(self.orientation) is not NeckTraversalOrientation:
            raise InvalidNeckPassageTransitionError("neck passage orientation must use the closed exact type.")
        if type(self.state) is not PassageState:
            raise InvalidNeckPassageTransitionError("neck passage state must use the closed exact type.")

    @classmethod
    def build(
        cls,
        *,
        neck: ClassifiedNeck,
        orientation: NeckTraversalOrientation,
        state: PassageState,
    ) -> Self:
        return cls(neck, orientation, state)

    @property
    def scope(self) -> OrientedNeckScope:
        return OrientedNeckScope.build(
            neck_owner_id=self.neck.owner_id,
            orientation=self.orientation,
        )

    def propose_cap_decision(self, policy: NeckPolicy) -> NeckCapDecision:
        if type(policy) is not NeckPolicy:
            raise InvalidNeckPassageTransitionError("neck cap proposal requires one exact neck policy.")
        passage_after = _NEXT_PASSAGE_STATE.get(self.state)
        if passage_after is None:
            raise TerminalNeckPassageError("terminal neck passage cannot propose another cap decision.")
        effective_cap = policy.effective_cap(
            self.neck.width_class_id.value,
            self.state,
        )
        return NeckCapDecision.build(
            neck_evidence_digest=self.neck.evidence_digest,
            width_class_id=self.neck.width_class_id,
            passage_before=self.state,
            passage_after=passage_after,
            user_cap=policy.user_cap,
            effective_cap=effective_cap,
        )

    def advance(self, decision: NeckCapDecision) -> Self:
        if type(decision) is not NeckCapDecision:
            raise InvalidNeckPassageTransitionError("neck passage advance requires one exact cap decision.")
        expected_after = _NEXT_PASSAGE_STATE.get(self.state)
        if expected_after is None:
            raise TerminalNeckPassageError("terminal neck passage cannot advance.")
        if (
            decision.neck_evidence_digest != self.neck.evidence_digest
            or decision.width_class_id != self.neck.width_class_id
            or decision.passage_before is not self.state
            or decision.passage_after is not expected_after
        ):
            raise InvalidNeckPassageTransitionError("neck cap decision contradicts its owned passage transition.")
        return type(self).build(
            neck=self.neck,
            orientation=self.orientation,
            state=expected_after,
        )


@dataclass(frozen=True)
class NeckInventory:
    axis: MedialAxis
    policy: NeckPolicy
    necks: tuple[ClassifiedNeck, ...]

    def __post_init__(self) -> None:
        if type(self.axis) is not MedialAxis:
            raise InvalidNeckEvidenceError("neck inventory requires one exact typed MAT.")
        if type(self.policy) is not NeckPolicy:
            raise InvalidNeckEvidenceError("neck inventory requires one exact neck policy.")
        if type(self.necks) is not tuple or any(type(neck) is not ClassifiedNeck for neck in self.necks):
            raise InvalidNeckEvidenceError("neck inventory must own one immutable typed tuple.")
        owner_ids = tuple(neck.owner_id for neck in self.necks)
        if owner_ids != tuple(sorted(set(owner_ids))):
            raise InvalidNeckEvidenceError("neck inventory owner identities must be unique and canonical.")
        if owner_ids != self.axis.native_owner.neck_owner_ids:
            raise InvalidNeckEvidenceError("neck inventory owner identities contradict its exact MAT owner.")
        if any(site_id not in self.axis.site_by_id for neck in self.necks for site_id in neck.defining_site_ids):
            raise InvalidNeckEvidenceError("neck inventory references an unknown MAT site.")
        if any(edge_id not in self.axis.edge_by_id for neck in self.necks for edge_id in neck.separating_cut_edge_ids):
            raise InvalidNeckEvidenceError("neck inventory references an unknown MAT edge.")

    @classmethod
    def build(
        cls,
        *,
        axis: MedialAxis,
        policy: NeckPolicy,
    ) -> Self:
        if type(axis) is not MedialAxis:
            raise InvalidNeckEvidenceError("neck classification requires one exact typed MAT.")
        if type(policy) is not NeckPolicy:
            raise InvalidNeckEvidenceError("neck classification requires one exact neck policy.")

        projection = axis.native_owner.projection
        evidence = projection[15]
        boundaries = tuple(ExactRationalV1.from_fraction(boundary).canonical_bytes for boundary in policy.squared_width_boundaries)
        classes, comparison_certificates = axis.native_owner.validate_and_classify_necks(
            axis.mat_certificate,
            evidence,
            boundaries,
        )
        if not isinstance(classes, np.ndarray) or classes.dtype != np.dtype(np.int64) or classes.ndim != 1 or classes.shape != (len(evidence),):
            raise InvalidNeckEvidenceError("native neck classes must be one complete int64 vector.")
        if type(comparison_certificates) is not tuple or len(comparison_certificates) != len(evidence):
            raise InvalidNeckEvidenceError("native neck comparison certificates have wrong cardinality.")
        owner_ids = axis.native_owner.neck_owner_ids
        defining_site_ids = axis.native_owner.neck_defining_site_ids
        if len(owner_ids) != len(evidence) or len(defining_site_ids) != len(evidence):
            raise InvalidNeckEvidenceError("native neck identities contradict evidence cardinality.")

        cut_offsets = projection[16]
        cut_edge_indices = projection[17]
        if (
            not isinstance(cut_offsets, np.ndarray)
            or cut_offsets.dtype != np.dtype(np.int64)
            or cut_offsets.shape != (len(evidence) + 1,)
            or not isinstance(cut_edge_indices, np.ndarray)
            or cut_edge_indices.dtype != np.dtype(np.int64)
            or cut_edge_indices.ndim != 1
        ):
            raise InvalidNeckEvidenceError("native neck separating-cut CSR has an invalid shape.")

        necks: list[ClassifiedNeck] = []
        class_count = len(policy.squared_width_boundaries) + 1
        for index, evidence_bytes in enumerate(evidence):
            class_id = int(classes[index])
            if class_id < 0 or class_id >= class_count:
                raise InvalidNeckEvidenceError("native neck width class is outside the policy partition.")
            begin = int(cut_offsets[index])
            end = int(cut_offsets[index + 1])
            cut_indices = cut_edge_indices[begin:end]
            if np.any(cut_indices < 0) or np.any(cut_indices >= len(axis.edges)):
                raise InvalidNeckEvidenceError("native neck separating cut references an unknown MAT edge.")
            necks.append(
                ClassifiedNeck(
                    NeckOwnerId(_exact_bytes(owner_ids[index], "neck owner identity")),
                    tuple(MatSiteId(_exact_bytes(site_id, "neck defining-site identity")) for site_id in defining_site_ids[index]),
                    _exact_bytes(evidence_bytes, "neck evidence"),
                    IdentityDigest(hashlib.sha256(evidence_bytes).digest()),
                    WidthClassId.build(class_id),
                    _exact_bytes(
                        comparison_certificates[index],
                        "neck comparison certificate",
                    ),
                    tuple(axis.edges[int(edge_index)].identity for edge_index in cut_indices),
                )
            )
        return cls(axis, policy, tuple(necks))
