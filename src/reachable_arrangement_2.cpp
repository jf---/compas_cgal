#include "reachable_arrangement_2.h"

#include "exact_sweep_2.h"
#include "reachable_errors_2.h"
#include "reachable_input_2.h"

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#include <CGAL/Boolean_set_operations_2/Gps_polygon_validation.h>
#include <CGAL/enum.h>

namespace {

template <typename T>
std::vector<T> sorted_unique(std::vector<T> values)
{
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());
    return values;
}

ReachPolygon ring_polygon(const CanonicalReachRing2& ring)
{
    ReachPolygon polygon;
    for (std::size_t index = 0; index < ring.points.size(); ++index) {
        polygon.push_back(ReachXCurve(
            ring.points[index],
            ring.points[(index + 1) % ring.points.size()]));
    }
    return polygon;
}

ReachPolygonWithHoles design_polygon(
    const CanonicalReachInput2& input)
{
    ReachPolygon outer = ring_polygon(input.outer);
    std::vector<ReachPolygon> holes;
    holes.reserve(input.holes.size());
    for (const CanonicalReachRing2& hole : input.holes) {
        holes.push_back(ring_polygon(hole));
    }
    ReachPolygonWithHoles design(
        std::move(outer),
        holes.begin(),
        holes.end());
    const ReachTraits traits;
    if (!CGAL::is_valid_polygon_with_holes<ReachTraits>(
            design,
            traits)) {
        throw InvalidReachableDomainInputError(
            "reachable-domain polygon with holes is invalid");
    }
    return design;
}

std::string ring_id(const CanonicalReachRing2& ring)
{
    if (ring.outer) {
        return "outer";
    }
    return reach_tagged_record(
        "canonical-hole-id-v1",
        {reach_u64_record(ring.canonical_ordinal)});
}

std::string primitive_id(
    const CanonicalReachRing2& ring,
    std::size_t ordinal,
    std::string_view role,
    double binary64_radius)
{
    return reach_tagged_record(
        "reach-primitive-v1",
        {
            ring_id(ring),
            reach_u64_record(ordinal),
            std::string(role),
            reach_binary64_record(binary64_radius),
        });
}

std::string design_primitive_id(const CanonicalReachRing2& ring)
{
    return reach_tagged_record(
        "reach-design-primitive-v1",
        {ring_id(ring)});
}

std::string source_record(
    const CanonicalReachRing2& ring,
    std::size_t ordinal,
    std::string_view role,
    double binary64_radius)
{
    return reach_tagged_record(
        "source-curve-v1",
        {
            ring_id(ring),
            reach_u64_record(ordinal),
            std::string(role),
            reach_binary64_record(binary64_radius),
        });
}

std::string source_piece_id(
    const std::string& source,
    std::size_t ordinal)
{
    return reach_tagged_record(
        "source-piece-v1",
        {source, reach_u64_record(ordinal)});
}

ReachDataTraits2::X_monotone_curve_2 labelled_curve(
    const ReachXCurve& curve,
    std::vector<std::string> source_piece_ids,
    std::vector<std::string> primitive_ids)
{
    return {
        curve,
        ReachCurveLabels2{
            sorted_unique(std::move(source_piece_ids)),
            sorted_unique(std::move(primitive_ids)),
        },
    };
}

void register_primitive(
    ReachPrimitiveKinds2& kinds,
    const std::string& id,
    ReachPrimitiveKind2 kind)
{
    if (!kinds.emplace(id, kind).second) {
        throw ReachableArrangementTopologyError(
            "reachable primitive identity collision");
    }
}

void append_ring_primitives(
    const CanonicalReachRing2& ring,
    const ReachFT& radius,
    double binary64_radius,
    std::vector<ReachDataTraits2::X_monotone_curve_2>& curves,
    ReachPrimitiveKinds2& primitive_kinds,
    std::vector<std::string>& source_records)
{
    const std::string design_id = design_primitive_id(ring);
    register_primitive(
        primitive_kinds,
        design_id,
        ring.outer
            ? ReachPrimitiveKind2::Outer
            : ReachPrimitiveKind2::Hole);
    for (std::size_t index = 0; index < ring.points.size(); ++index) {
        curves.push_back(labelled_curve(
            ReachXCurve(
                ring.points[index],
                ring.points[(index + 1) % ring.points.size()]),
            {},
            {design_id}));
    }

    for (std::size_t index = 0; index < ring.points.size(); ++index) {
        const ReachKernelPoint& start = ring.points[index];
        const ReachKernelPoint& end =
            ring.points[(index + 1) % ring.points.size()];
        const ReachKernelVector direction = end - start;
        if (direction == ReachKernelVector(CGAL::NULL_VECTOR)) {
            throw InvalidReachableDomainInputError(
                "reachable-domain edges must be nondegenerate");
        }
        const ReachFT length = CGAL::sqrt(direction.squared_length());
        const ReachKernelVector normal =
            direction.perpendicular(CGAL::COUNTERCLOCKWISE)
            * (radius / length);
        const ReachKernelPoint minus_start = start - normal;
        const ReachKernelPoint minus_end = end - normal;
        const ReachKernelPoint plus_start = start + normal;
        const ReachKernelPoint plus_end = end + normal;
        const std::string strip_id = primitive_id(
            ring,
            index,
            "edge-strip",
            binary64_radius);
        register_primitive(
            primitive_kinds,
            strip_id,
            ReachPrimitiveKind2::Forbidden);

        const std::string minus_source = source_record(
            ring,
            index,
            "offset-minus",
            binary64_radius);
        const std::string plus_source = source_record(
            ring,
            index,
            "offset-plus",
            binary64_radius);
        source_records.push_back(minus_source);
        source_records.push_back(plus_source);
        curves.push_back(labelled_curve(
            ReachXCurve(minus_start, minus_end),
            {source_piece_id(minus_source, 0)},
            {strip_id}));
        curves.push_back(labelled_curve(
            ReachXCurve(minus_end, plus_end),
            {},
            {strip_id}));
        curves.push_back(labelled_curve(
            ReachXCurve(plus_end, plus_start),
            {source_piece_id(plus_source, 0)},
            {strip_id}));
        curves.push_back(labelled_curve(
            ReachXCurve(plus_start, minus_start),
            {},
            {strip_id}));
    }

    for (std::size_t index = 0; index < ring.points.size(); ++index) {
        const std::string disk_id = primitive_id(
            ring,
            index,
            "vertex-disk",
            binary64_radius);
        register_primitive(
            primitive_kinds,
            disk_id,
            ReachPrimitiveKind2::Forbidden);
        const std::string circle_source = source_record(
            ring,
            index,
            "vertex-circle",
            binary64_radius);
        source_records.push_back(circle_source);
        const ReachPolygon disk =
            reach_disk_polygon(ring.points[index], radius);
        std::size_t piece_ordinal = 0;
        for (auto piece = disk.curves_begin();
             piece != disk.curves_end();
             ++piece) {
            curves.push_back(labelled_curve(
                *piece,
                {source_piece_id(circle_source, piece_ordinal)},
                {disk_id}));
            ++piece_ordinal;
        }
        if (piece_ordinal != 2) {
            throw ReachableArrangementTopologyError(
                "vertex disk did not emit two x-monotone source pieces");
        }
    }
}

using FaceHandle2 = ReachArrangement2::Face_handle;
using HalfedgeHandle2 = ReachArrangement2::Halfedge_handle;

template <typename Handle>
const void* handle_address(const Handle& handle)
{
    return static_cast<const void*>(std::addressof(*handle));
}

struct DenseHalfedge2 {
    HalfedgeHandle2 handle;
    std::size_t face_slot;
    std::size_t neighbor_slot;
    std::size_t edge_slot;
};

struct DenseArrangementIndex2 {
    std::vector<FaceHandle2> faces;
    std::vector<DenseHalfedge2> halfedges;
    std::vector<std::vector<std::size_t>> halfedge_slots_by_face;
    std::vector<std::vector<std::size_t>> primitive_slots_by_edge;
    std::vector<ReachPrimitiveKind2> primitive_kinds;
    std::unordered_map<const void*, std::size_t>
        halfedge_slot_by_address;
    std::size_t unbounded_face_slot = 0;
};

DenseArrangementIndex2 compile_dense_arrangement(
    ReachArrangement2& arrangement,
    const ReachPrimitiveKinds2& primitive_kinds,
    ReachableDomainBuildAudit2& audit)
{
    DenseArrangementIndex2 dense;
    const std::size_t face_count = arrangement.number_of_faces();
    const std::size_t halfedge_count =
        arrangement.number_of_halfedges();
    dense.faces.reserve(face_count);
    dense.halfedges.reserve(halfedge_count);
    dense.halfedge_slots_by_face.resize(face_count);
    dense.primitive_slots_by_edge.reserve(
        arrangement.number_of_edges());
    dense.primitive_kinds.reserve(primitive_kinds.size());
    dense.halfedge_slot_by_address.reserve(halfedge_count);

    std::unordered_map<const void*, std::size_t>
        face_slot_by_address;
    face_slot_by_address.reserve(face_count);
    for (auto face = arrangement.faces_begin();
         face != arrangement.faces_end();
         ++face) {
        const std::size_t slot = dense.faces.size();
        if (!face_slot_by_address.emplace(
                handle_address(face),
                slot).second) {
            throw ReachableArrangementTopologyError(
                "dense arrangement contains a duplicate face address");
        }
        dense.faces.push_back(face);
        ++audit.dense_face_visits;
    }

    std::unordered_map<std::string, std::size_t>
        primitive_slot_by_id;
    primitive_slot_by_id.reserve(primitive_kinds.size());
    std::size_t outer_count = 0;
    for (const auto& [id, kind] : primitive_kinds) {
        const std::size_t slot = dense.primitive_kinds.size();
        primitive_slot_by_id.emplace(id, slot);
        dense.primitive_kinds.push_back(kind);
        outer_count += kind == ReachPrimitiveKind2::Outer ? 1 : 0;
    }
    if (outer_count != 1) {
        throw ReachableArrangementTopologyError(
            "reachable parity requires exactly one outer primitive");
    }

    const auto face_slot = [&face_slot_by_address](
                               FaceHandle2 face) {
        const auto found =
            face_slot_by_address.find(handle_address(face));
        if (found == face_slot_by_address.end()) {
            throw ReachableArrangementTopologyError(
                "dense arrangement references an unknown face");
        }
        return found->second;
    };
    for (auto edge = arrangement.edges_begin();
         edge != arrangement.edges_end();
         ++edge) {
        const std::vector<std::string>& primitive_ids =
            edge->curve().data().primitive_ids;
        if (!std::is_sorted(
                primitive_ids.begin(),
                primitive_ids.end())
            || std::adjacent_find(
                   primitive_ids.begin(),
                   primitive_ids.end())
                != primitive_ids.end()) {
            throw ReachableArrangementTopologyError(
                "arrangement primitive labels must be sorted and unique");
        }
        std::vector<std::size_t> primitive_slots;
        primitive_slots.reserve(primitive_ids.size());
        for (const std::string& id : primitive_ids) {
            ++audit.primitive_label_resolutions;
            const auto found = primitive_slot_by_id.find(id);
            if (found == primitive_slot_by_id.end()) {
                throw ReachableArrangementTopologyError(
                    "arrangement edge references an unknown primitive");
            }
            primitive_slots.push_back(found->second);
        }
        const std::size_t edge_slot =
            dense.primitive_slots_by_edge.size();
        dense.primitive_slots_by_edge.push_back(
            std::move(primitive_slots));
        for (const HalfedgeHandle2 halfedge :
             {HalfedgeHandle2(edge), edge->twin()}) {
            const std::size_t halfedge_slot =
                dense.halfedges.size();
            const std::size_t owner_slot =
                face_slot(halfedge->face());
            if (!dense.halfedge_slot_by_address.emplace(
                    handle_address(halfedge),
                    halfedge_slot).second) {
                throw ReachableArrangementTopologyError(
                    "dense arrangement contains a duplicate halfedge address");
            }
            dense.halfedges.push_back(
                {
                    halfedge,
                    owner_slot,
                    face_slot(halfedge->twin()->face()),
                    edge_slot,
                });
            dense.halfedge_slots_by_face[owner_slot].push_back(
                halfedge_slot);
            ++audit.dense_halfedge_visits;
        }
    }
    if (dense.faces.size() != face_count
        || dense.halfedges.size() != halfedge_count) {
        throw ReachableArrangementTopologyError(
            "dense arrangement indexing did not cover the complete DCEL");
    }
    dense.unbounded_face_slot =
        face_slot(arrangement.unbounded_face());
    return dense;
}

class ActivePrimitiveParity2 {
public:
    explicit ActivePrimitiveParity2(
        const std::vector<ReachPrimitiveKind2>& primitive_kinds)
        : primitive_kinds_(primitive_kinds),
          active_(primitive_kinds.size(), false)
    {
    }

    void toggle(const std::vector<std::size_t>& primitive_slots)
    {
        for (const std::size_t slot : primitive_slots) {
            if (slot >= active_.size()) {
                throw ReachableArrangementTopologyError(
                    "primitive parity references an invalid dense slot");
            }
            const bool activated = active_[slot] == 0;
            const std::size_t increment = activated ? 1 : 0;
            const std::size_t decrement = activated ? 0 : 1;
            active_[slot] = activated;
            active_count_ += increment;
            active_count_ -= decrement;
            std::size_t& count =
                count_for(primitive_kinds_[slot]);
            if (count < decrement) {
                throw ReachableArrangementTopologyError(
                    "primitive parity count underflow");
            }
            count += increment;
            count -= decrement;
        }
    }

    ReachFaceState2 face_state() const
    {
        if (active_outer_ > 1) {
            throw ReachableArrangementTopologyError(
                "multiple outer primitives are active");
        }
        ReachFaceState2 state;
        state.classified = true;
        state.outer_active = active_outer_ == 1;
        state.active_holes = active_holes_;
        state.active_forbidden = active_forbidden_;
        state.selected =
            state.outer_active
            && state.active_holes == 0
            && state.active_forbidden == 0;
        return state;
    }

    bool empty() const
    {
        return active_count_ == 0
            && active_outer_ == 0
            && active_holes_ == 0
            && active_forbidden_ == 0;
    }

private:
    std::size_t& count_for(ReachPrimitiveKind2 kind)
    {
        switch (kind) {
        case ReachPrimitiveKind2::Outer:
            return active_outer_;
        case ReachPrimitiveKind2::Hole:
            return active_holes_;
        case ReachPrimitiveKind2::Forbidden:
            return active_forbidden_;
        }
        throw ReachableArrangementTopologyError(
            "unknown reachable primitive kind");
    }

    const std::vector<ReachPrimitiveKind2>& primitive_kinds_;
    std::vector<unsigned char> active_;
    std::size_t active_count_ = 0;
    std::size_t active_outer_ = 0;
    std::size_t active_holes_ = 0;
    std::size_t active_forbidden_ = 0;
};

bool same_face_state(
    const ReachFaceState2& left,
    const ReachFaceState2& right)
{
    return left.selected == right.selected
        && left.outer_active == right.outer_active
        && left.active_holes == right.active_holes
        && left.active_forbidden == right.active_forbidden;
}

struct ParityFrame2 {
    std::size_t face_slot;
    std::size_t next_halfedge = 0;
    std::optional<std::size_t> entering_edge_slot;
};

void classify_dense_faces(
    const DenseArrangementIndex2& dense,
    ReachableDomainBuildAudit2& audit)
{
    for (const FaceHandle2 face : dense.faces) {
        face->set_data(ReachFaceState2{});
    }
    ActivePrimitiveParity2 active(dense.primitive_kinds);
    dense.faces[dense.unbounded_face_slot]->set_data(
        active.face_state());
    std::vector<ParityFrame2> stack{
        {dense.unbounded_face_slot, 0, std::nullopt},
    };
    while (!stack.empty()) {
        ParityFrame2& frame = stack.back();
        const std::vector<std::size_t>& face_halfedges =
            dense.halfedge_slots_by_face[frame.face_slot];
        if (frame.next_halfedge == face_halfedges.size()) {
            const std::optional<std::size_t> entering =
                frame.entering_edge_slot;
            stack.pop_back();
            if (entering.has_value()) {
                active.toggle(
                    dense.primitive_slots_by_edge[*entering]);
            }
            continue;
        }
        const DenseHalfedge2& crossed =
            dense.halfedges[face_halfedges[frame.next_halfedge++]];
        ++audit.parity_halfedge_visits;
        active.toggle(
            dense.primitive_slots_by_edge[crossed.edge_slot]);
        const ReachFaceState2 predicted = active.face_state();
        if (!dense.faces[crossed.neighbor_slot]
                 ->data()
                 .classified) {
            dense.faces[crossed.neighbor_slot]->set_data(
                predicted);
            stack.push_back(
                {crossed.neighbor_slot, 0, crossed.edge_slot});
            continue;
        }
        if (!same_face_state(
                dense.faces[crossed.neighbor_slot]->data(),
                predicted)) {
            throw ReachableArrangementTopologyError(
                "dual routes disagree on aggregate face parity");
        }
        active.toggle(
            dense.primitive_slots_by_edge[crossed.edge_slot]);
    }
    if (!active.empty()) {
        throw ReachableArrangementTopologyError(
            "dual parity traversal did not restore the unbounded state");
    }
    for (const FaceHandle2 face : dense.faces) {
        if (!face->data().classified) {
            throw ReachableArrangementTopologyError(
                "dual parity traversal left a face unclassified");
        }
    }
}

std::size_t selected_component_count(
    const DenseArrangementIndex2& dense,
    ReachableDomainBuildAudit2& audit)
{
    std::vector<unsigned char> visited(dense.faces.size(), false);
    std::size_t component_count = 0;
    for (std::size_t face_slot = 0;
         face_slot < dense.faces.size();
         ++face_slot) {
        if (!dense.faces[face_slot]->data().selected
            || visited[face_slot] != 0) {
            continue;
        }
        ++component_count;
        std::vector<std::size_t> pending{face_slot};
        visited[face_slot] = true;
        while (!pending.empty()) {
            const std::size_t current = pending.back();
            pending.pop_back();
            for (const std::size_t halfedge_slot :
                 dense.halfedge_slots_by_face[current]) {
                ++audit.component_halfedge_visits;
                const std::size_t neighbor =
                    dense.halfedges[halfedge_slot].neighbor_slot;
                if (dense.faces[neighbor]->data().selected
                    && visited[neighbor] == 0) {
                    visited[neighbor] = true;
                    pending.push_back(neighbor);
                }
            }
        }
    }
    return component_count;
}

std::size_t halfedge_slot(
    const DenseArrangementIndex2& dense,
    HalfedgeHandle2 halfedge)
{
    const auto found = dense.halfedge_slot_by_address.find(
        handle_address(halfedge));
    if (found == dense.halfedge_slot_by_address.end()) {
        throw ReachableArrangementTopologyError(
            "boundary traversal references an unknown dense halfedge");
    }
    return found->second;
}

std::size_t next_selected_boundary(
    const DenseArrangementIndex2& dense,
    std::size_t current_slot,
    ReachableDomainBuildAudit2& audit)
{
    HalfedgeHandle2 candidate =
        dense.halfedges[current_slot].handle->next();
    for (std::size_t step = 0;
         step < dense.halfedges.size();
         ++step) {
        ++audit.boundary_rotation_halfedge_visits;
        if (!candidate->face()->data().selected) {
            throw ReachableArrangementTopologyError(
                "selected boundary traversal left the selected face set");
        }
        if (!candidate->twin()->face()->data().selected) {
            return halfedge_slot(dense, candidate);
        }
        candidate = candidate->twin()->next();
    }
    throw ReachableArrangementTopologyError(
        "selected boundary traversal did not close around a vertex");
}

ReachXCurve directed_halfedge_curve(HalfedgeHandle2 halfedge)
{
    const ReachTraits traits;
    const auto compare = traits.compare_xy_2_object();
    ReachXCurve curve =
        static_cast<const ReachXCurve&>(halfedge->curve());
    if (compare(curve.source(), halfedge->source()->point())
        != CGAL::EQUAL) {
        curve = traits.construct_opposite_2_object()(curve);
    }
    if (compare(curve.source(), halfedge->source()->point())
            != CGAL::EQUAL
        || compare(curve.target(), halfedge->target()->point())
            != CGAL::EQUAL) {
        throw ReachableArrangementTopologyError(
            "selected halfedge curve does not match its exact DCEL endpoints");
    }
    return curve;
}

ReachPolygonWithHoles extract_center_polygon(
    const DenseArrangementIndex2& dense,
    ReachableDomainBuildAudit2& audit)
{
    std::vector<unsigned char> visited(
        dense.halfedges.size(),
        false);
    std::vector<ReachPolygon> outer_cycles;
    std::vector<ReachPolygon> hole_cycles;
    for (std::size_t halfedge_slot_value = 0;
         halfedge_slot_value < dense.halfedges.size();
         ++halfedge_slot_value) {
        ++audit.boundary_scan_halfedge_visits;
        const HalfedgeHandle2 halfedge =
            dense.halfedges[halfedge_slot_value].handle;
        const bool selected_boundary =
            halfedge->face()->data().selected
            && !halfedge->twin()->face()->data().selected;
        if (!selected_boundary
            || visited[halfedge_slot_value] != 0) {
            continue;
        }
        ReachPolygon cycle;
        const std::size_t first = halfedge_slot_value;
        std::size_t current = first;
        for (std::size_t step = 0;
             step < dense.halfedges.size();
             ++step) {
            if (visited[current] != 0) {
                throw ReachableArrangementTopologyError(
                    "selected boundary cycle revisited a halfedge before closure");
            }
            visited[current] = true;
            ++audit.boundary_cycle_halfedge_visits;
            const HalfedgeHandle2 current_halfedge =
                dense.halfedges[current].handle;
            if (current_halfedge->curve().data().source_piece_ids.empty()) {
                throw ReachableArrangementTopologyError(
                    "selected boundary edge has no propagated source provenance");
            }
            cycle.push_back(
                directed_halfedge_curve(current_halfedge));
            current = next_selected_boundary(
                dense,
                current,
                audit);
            if (current == first) {
                break;
            }
        }
        if (current != first || cycle.is_empty()) {
            throw ReachableArrangementTopologyError(
                "selected boundary cycle did not close");
        }
        const CGAL::Orientation orientation = cycle.orientation();
        if (orientation == CGAL::COUNTERCLOCKWISE) {
            outer_cycles.push_back(std::move(cycle));
        }
        else if (orientation == CGAL::CLOCKWISE) {
            hole_cycles.push_back(std::move(cycle));
        }
        else {
            throw ReachableArrangementTopologyError(
                "selected boundary cycle has zero exact area");
        }
    }
    if (outer_cycles.size() != 1) {
        throw ReachableArrangementTopologyError(
            "selected center domain must have exactly one outer cycle");
    }
    ReachPolygonWithHoles center(
        std::move(outer_cycles.front()),
        hole_cycles.begin(),
        hole_cycles.end());
    const ReachTraits traits;
    if (!CGAL::is_valid_polygon_with_holes<ReachTraits>(
            center,
            traits)) {
        throw ReachableArrangementTopologyError(
            "selected center polygon with holes is invalid");
    }
    return center;
}

} // namespace

ReachCurveLabels2 MergeReachCurveLabels2::operator()(
    const ReachCurveLabels2& left,
    const ReachCurveLabels2& right) const
{
    std::vector<std::string> source_piece_ids = left.source_piece_ids;
    source_piece_ids.insert(
        source_piece_ids.end(),
        right.source_piece_ids.begin(),
        right.source_piece_ids.end());
    std::vector<std::string> primitive_ids = left.primitive_ids;
    primitive_ids.insert(
        primitive_ids.end(),
        right.primitive_ids.begin(),
        right.primitive_ids.end());
    return {
        sorted_unique(std::move(source_piece_ids)),
        sorted_unique(std::move(primitive_ids)),
    };
}

void classify_faces_by_primitive_parity(
    ReachArrangement2& arrangement,
    const ReachPrimitiveKinds2& primitive_kinds)
{
    ReachableDomainBuildAudit2 audit;
    classify_faces_by_primitive_parity(
        arrangement,
        primitive_kinds,
        audit);
}

void classify_faces_by_primitive_parity(
    ReachArrangement2& arrangement,
    const ReachPrimitiveKinds2& primitive_kinds,
    ReachableDomainBuildAudit2& audit)
{
    const DenseArrangementIndex2 dense =
        compile_dense_arrangement(
            arrangement,
            primitive_kinds,
            audit);
    classify_dense_faces(dense, audit);
}

ReachableArrangement2 build_reachable_arrangement(
    CanonicalReachInput2 input)
{
    validate_canonical_reach_input(input);
    ReachableArrangement2 result{
        ReachArrangement2{},
        std::move(input),
        {},
        ReachPolygonWithHoles{},
        ReachPolygonWithHoles{},
        {},
    };
    result.audit.input_vertex_count = result.input.input_vertex_count_;
    result.audit.ring_rotation_comparisons =
        result.input.ring_rotation_comparisons_;
    result.design_polygon = design_polygon(result.input);

    std::vector<ReachDataTraits2::X_monotone_curve_2> labelled;
    ReachPrimitiveKinds2 primitive_kinds;
    append_ring_primitives(
        result.input.outer,
        result.input.radius,
        result.input.binary64_radius,
        labelled,
        primitive_kinds,
        result.source_records);
    for (const CanonicalReachRing2& hole : result.input.holes) {
        append_ring_primitives(
            hole,
            result.input.radius,
            result.input.binary64_radius,
            labelled,
            primitive_kinds,
            result.source_records);
    }
    std::sort(result.source_records.begin(), result.source_records.end());
    CGAL::insert(
        result.arrangement,
        labelled.begin(),
        labelled.end());
    ++result.audit.provenance_arrangements;
    const DenseArrangementIndex2 dense =
        compile_dense_arrangement(
            result.arrangement,
            primitive_kinds,
            result.audit);
    classify_dense_faces(dense, result.audit);

    for (const FaceHandle2 face : dense.faces) {
        if (face->data().selected) {
            ++result.audit.selected_face_count;
        }
    }
    for (std::size_t halfedge_slot_value = 0;
         halfedge_slot_value < dense.halfedges.size();
         halfedge_slot_value += 2) {
        const DenseHalfedge2& halfedge =
            dense.halfedges[halfedge_slot_value];
        if (dense.faces[halfedge.face_slot]->data().selected
            && dense.faces[halfedge.neighbor_slot]
                   ->data()
                   .selected) {
            ++result.audit.selected_adjacency_count;
        }
    }
    const std::size_t component_count =
        selected_component_count(dense, result.audit);
    if (component_count == 0) {
        throw PocketNotMachinableError(
            "Phase 1 center domain is empty for the declared tool radius");
    }
    if (component_count != 1) {
        throw PocketNotMachinableError(
            "Phase 1 requires one connected center domain");
    }
    result.center_polygon =
        extract_center_polygon(dense, result.audit);
    ++result.audit.center_extractions;
    return result;
}
