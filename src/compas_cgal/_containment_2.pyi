from compas_cgal._coverage_2 import ExactRegion2

class ContainmentConstructionError(RuntimeError): ...

class ContainmentRecord2:
    @property
    def contained(self) -> bool: ...
    @property
    def guide_anchor_in_center_domain(self) -> bool: ...
    @property
    def outer_disk_contained(self) -> bool: ...
    @property
    def disk_sweep(self) -> bool: ...
    @property
    def strategy_version(self) -> bytes: ...
    @property
    def structural_record(self) -> bytes: ...
    def matches_exact_segment(
        self,
        authority_digest: bytes,
        x0: float,
        y0: float,
        x1: float,
        y1: float,
        tool_radius: float,
    ) -> bool: ...
    def matches_exact_full_circle(
        self,
        authority_digest: bytes,
        center_x: float,
        center_y: float,
        phase_x: float,
        phase_y: float,
        tool_radius: float,
    ) -> bool: ...
    def matches_exact_entry_disk(
        self,
        authority_digest: bytes,
        center_x: float,
        center_y: float,
        entry_radius: float,
        tool_radius: float,
    ) -> bool: ...
    def matches_exact_full_circle_in_disk(
        self,
        authority_digest: bytes,
        entry_center_x: float,
        entry_center_y: float,
        entry_radius: float,
        circle_center_x: float,
        circle_center_y: float,
        phase_x: float,
        phase_y: float,
        tool_radius: float,
    ) -> bool: ...

def evaluate_exact_segment_containment(
    design: ExactRegion2,
    center_domain: ExactRegion2,
    authority_digest: bytes,
    x0: float,
    y0: float,
    x1: float,
    y1: float,
    tool_radius: float,
) -> ContainmentRecord2: ...
def evaluate_exact_full_circle_containment(
    design: ExactRegion2,
    center_domain: ExactRegion2,
    authority_digest: bytes,
    center_x: float,
    center_y: float,
    phase_x: float,
    phase_y: float,
    tool_radius: float,
) -> ContainmentRecord2: ...
def evaluate_exact_entry_disk_containment(
    design: ExactRegion2,
    center_domain: ExactRegion2,
    authority_digest: bytes,
    center_x: float,
    center_y: float,
    entry_radius: float,
    tool_radius: float,
) -> ContainmentRecord2: ...
def evaluate_exact_full_circle_in_disk(
    authority_digest: bytes,
    entry_center_x: float,
    entry_center_y: float,
    entry_radius: float,
    circle_center_x: float,
    circle_center_y: float,
    phase_x: float,
    phase_y: float,
    tool_radius: float,
) -> ContainmentRecord2: ...
def containment_strategy_version() -> bytes: ...
