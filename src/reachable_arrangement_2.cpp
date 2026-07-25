#include "reachable_arrangement_2.h"

#include "canonical_rotation_2.h"
#include "exact_sweep_2.h"
#include "reachable_errors_2.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <memory>
#include <set>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <CGAL/Boolean_set_operations_2/Gps_polygon_validation.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/enum.h>
#include <CGAL/number_utils.h>

namespace {

template <typename T>
std::vector<T> sorted_unique(std::vector<T> values)
{
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());
    return values;
}

bool point_equal(
    const ReachKernelPoint& left,
    const ReachKernelPoint& right)
{
    return CGAL::compare_xy(left, right) == CGAL::EQUAL;
}

std::string ring_record(
    bool outer,
    const std::vector<std::array<double, 2>>& points)
{
    std::vector<std::string> fields;
    fields.reserve(points.size() + 1);
    fields.emplace_back(outer ? "outer" : "hole");
    for (const std::array<double, 2>& point : points) {
        fields.push_back(reach_tagged_record(
            "binary64-point-v1",
            {
                reach_binary64_record(point[0]),
                reach_binary64_record(point[1]),
            }));
    }
    return reach_tagged_record("input-ring-v1", fields);
}

void rotate_ring(
    CanonicalReachRing2& ring,
    std::size_t& comparison_count)
{
    const std::size_t first = minimal_rotation_index(
        ring.points,
        [](const ReachKernelPoint& left, const ReachKernelPoint& right) {
            return CGAL::compare_xy(left, right);
        },
        comparison_count);
    std::rotate(
        ring.points.begin(),
        ring.points.begin() + static_cast<std::ptrdiff_t>(first),
        ring.points.end());
    std::rotate(
        ring.binary64_points.begin(),
        ring.binary64_points.begin() + static_cast<std::ptrdiff_t>(first),
        ring.binary64_points.end());
}

CanonicalReachRing2 canonical_ring(
    Eigen::Ref<const compas::RowMatrixXd> vertices,
    bool outer,
    std::size_t& comparison_count)
{
    if (vertices.cols() < 2 || vertices.rows() < 3) {
        throw InvalidReachableDomainInputError(
            "reachable-domain rings require at least three XY vertices");
    }
    CanonicalReachRing2 ring{outer, 0, {}, {}, {}};
    ring.points.reserve(static_cast<std::size_t>(vertices.rows()));
    ring.binary64_points.reserve(
        static_cast<std::size_t>(vertices.rows()));
    for (Eigen::Index row = 0; row < vertices.rows(); ++row) {
        const double x = vertices(row, 0);
        const double y = vertices(row, 1);
        if (!std::isfinite(x) || !std::isfinite(y)) {
            throw InvalidReachableDomainInputError(
                "reachable-domain vertices must be finite binary64 values");
        }
        const double canonical_x = x == 0.0 ? 0.0 : x;
        const double canonical_y = y == 0.0 ? 0.0 : y;
        ring.points.emplace_back(canonical_x, canonical_y);
        ring.binary64_points.push_back({canonical_x, canonical_y});
    }
    if (ring.points.size() > 3
        && point_equal(ring.points.front(), ring.points.back())) {
        ring.points.pop_back();
        ring.binary64_points.pop_back();
    }
    if (ring.points.size() < 3) {
        throw InvalidReachableDomainInputError(
            "reachable-domain rings require three distinct vertices");
    }
    for (std::size_t index = 0; index < ring.points.size(); ++index) {
        if (point_equal(
                ring.points[index],
                ring.points[(index + 1) % ring.points.size()])) {
            throw InvalidReachableDomainInputError(
                "reachable-domain edges must be nondegenerate");
        }
    }

    CGAL::Polygon_2<ReachKernel> polygon(
        ring.points.begin(),
        ring.points.end());
    if (!polygon.is_simple()) {
        throw InvalidReachableDomainInputError(
            "reachable-domain rings must be exact simple polygons");
    }
    const CGAL::Orientation orientation = polygon.orientation();
    if (orientation == CGAL::COLLINEAR) {
        throw InvalidReachableDomainInputError(
            "reachable-domain rings must have nonzero exact area");
    }
    if ((orientation == CGAL::COUNTERCLOCKWISE) != outer) {
        std::reverse(ring.points.begin(), ring.points.end());
        std::reverse(
            ring.binary64_points.begin(),
            ring.binary64_points.end());
    }
    rotate_ring(ring, comparison_count);
    ring.record = ring_record(outer, ring.binary64_points);
    return ring;
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

template <typename CcbIterator>
void append_ccb_halfedges(
    CcbIterator begin,
    CcbIterator end,
    std::vector<HalfedgeHandle2>& halfedges)
{
    for (auto ccb = begin; ccb != end; ++ccb) {
        auto first = *ccb;
        auto current = first;
        do {
            halfedges.push_back(current);
            ++current;
        } while (current != first);
    }
}

std::vector<HalfedgeHandle2> face_halfedges(FaceHandle2 face)
{
    std::vector<HalfedgeHandle2> halfedges;
    append_ccb_halfedges(
        face->outer_ccbs_begin(),
        face->outer_ccbs_end(),
        halfedges);
    append_ccb_halfedges(
        face->inner_ccbs_begin(),
        face->inner_ccbs_end(),
        halfedges);
    return halfedges;
}

class ActivePrimitiveParity2 {
public:
    explicit ActivePrimitiveParity2(
        const ReachPrimitiveKinds2& primitive_kinds)
        : primitive_kinds_(primitive_kinds)
    {
    }

    void toggle(const std::vector<std::string>& primitive_ids)
    {
        for (const std::string& id : primitive_ids) {
            const auto kind = primitive_kinds_.find(id);
            if (kind == primitive_kinds_.end()) {
                throw ReachableArrangementTopologyError(
                    "arrangement edge references an unknown primitive");
            }
            const bool activated = active_ids_.insert(id).second;
            const std::size_t increment = activated ? 1 : 0;
            const std::size_t decrement = activated ? 0 : 1;
            if (!activated && active_ids_.erase(id) != 1) {
                throw ReachableArrangementTopologyError(
                    "primitive parity state could not deactivate an active ID");
            }
            std::size_t& count = count_for(kind->second);
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
        return active_ids_.empty()
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

    const ReachPrimitiveKinds2& primitive_kinds_;
    std::set<std::string> active_ids_;
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
    FaceHandle2 face;
    std::vector<HalfedgeHandle2> halfedges;
    std::size_t next_halfedge = 0;
    std::vector<std::string> entering_primitives;
};

std::size_t selected_component_count(
    ReachArrangement2& arrangement)
{
    std::set<const void*> visited;
    std::size_t component_count = 0;
    for (auto face = arrangement.faces_begin();
         face != arrangement.faces_end();
         ++face) {
        if (!face->data().selected
            || visited.contains(std::addressof(*face))) {
            continue;
        }
        ++component_count;
        std::vector<FaceHandle2> pending{face};
        visited.insert(std::addressof(*face));
        while (!pending.empty()) {
            const FaceHandle2 current = pending.back();
            pending.pop_back();
            for (const HalfedgeHandle2 halfedge :
                 face_halfedges(current)) {
                const FaceHandle2 neighbor =
                    halfedge->twin()->face();
                if (neighbor->data().selected
                    && visited.insert(
                        std::addressof(*neighbor)).second) {
                    pending.push_back(neighbor);
                }
            }
        }
    }
    return component_count;
}

HalfedgeHandle2 next_selected_boundary(
    ReachArrangement2& arrangement,
    HalfedgeHandle2 current)
{
    HalfedgeHandle2 candidate = current->next();
    for (std::size_t step = 0;
         step < arrangement.number_of_halfedges();
         ++step) {
        if (!candidate->face()->data().selected) {
            throw ReachableArrangementTopologyError(
                "selected boundary traversal left the selected face set");
        }
        if (!candidate->twin()->face()->data().selected) {
            return candidate;
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
    ReachArrangement2& arrangement)
{
    std::set<const void*> visited;
    std::vector<ReachPolygon> outer_cycles;
    std::vector<ReachPolygon> hole_cycles;
    for (auto halfedge = arrangement.halfedges_begin();
         halfedge != arrangement.halfedges_end();
         ++halfedge) {
        const bool selected_boundary =
            halfedge->face()->data().selected
            && !halfedge->twin()->face()->data().selected;
        if (!selected_boundary
            || visited.contains(std::addressof(*halfedge))) {
            continue;
        }
        ReachPolygon cycle;
        const HalfedgeHandle2 first = halfedge;
        HalfedgeHandle2 current = first;
        for (std::size_t step = 0;
             step < arrangement.number_of_halfedges();
             ++step) {
            if (!visited.insert(std::addressof(*current)).second) {
                throw ReachableArrangementTopologyError(
                    "selected boundary cycle revisited a halfedge before closure");
            }
            if (current->curve().data().source_piece_ids.empty()) {
                throw ReachableArrangementTopologyError(
                    "selected boundary edge has no propagated source provenance");
            }
            cycle.push_back(directed_halfedge_curve(current));
            current = next_selected_boundary(arrangement, current);
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

void validate_canonical_radius(const CanonicalReachInput2& input)
{
    if (!std::isfinite(input.binary64_radius)
        || !(input.binary64_radius > 0.0)
        || CGAL::compare(input.radius, ReachFT(0)) != CGAL::LARGER
        || CGAL::compare(
               input.radius,
               ReachFT(input.binary64_radius))
            != CGAL::EQUAL) {
        throw InvalidReachableDomainInputError(
            "canonical reach radius must be a positive exact binary64 injection");
    }
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

CanonicalReachInput2 canonical_reach_input(
    Eigen::Ref<const compas::RowMatrixXd> boundary,
    const std::vector<compas::RowMatrixXd>& holes,
    double tool_radius)
{
    if (!std::isfinite(tool_radius) || !(tool_radius > 0.0)) {
        throw InvalidReachableDomainInputError(
            "reachable-domain tool radius must be finite and positive");
    }
    std::size_t comparison_count = 0;
    CanonicalReachInput2 input;
    input.outer = canonical_ring(boundary, true, comparison_count);
    input.radius = ReachFT(tool_radius);
    input.binary64_radius = tool_radius;
    input.holes.reserve(holes.size());
    for (const compas::RowMatrixXd& hole : holes) {
        input.holes.push_back(
            canonical_ring(hole, false, comparison_count));
    }
    std::sort(
        input.holes.begin(),
        input.holes.end(),
        [](const CanonicalReachRing2& left,
           const CanonicalReachRing2& right) {
            return left.record < right.record;
        });
    for (std::size_t index = 0; index < input.holes.size(); ++index) {
        if (index != 0
            && input.holes[index - 1].record
                == input.holes[index].record) {
            throw InvalidReachableDomainInputError(
                "reachable-domain holes must have distinct exact identities");
        }
        input.holes[index].canonical_ordinal = index;
    }
    input.input_vertex_count_ = input.outer.points.size();
    for (const CanonicalReachRing2& hole : input.holes) {
        input.input_vertex_count_ += hole.points.size();
    }
    input.ring_rotation_comparisons_ = comparison_count;
    std::vector<std::string> fields{
        input.outer.record,
        reach_binary64_record(input.binary64_radius),
    };
    for (const CanonicalReachRing2& hole : input.holes) {
        fields.push_back(hole.record);
    }
    input.recipe_record =
        reach_tagged_record("reachable-input-recipe-v1", fields);
    return input;
}

void classify_faces_by_primitive_parity(
    ReachArrangement2& arrangement,
    const ReachPrimitiveKinds2& primitive_kinds)
{
    const std::size_t outer_count = static_cast<std::size_t>(
        std::count_if(
            primitive_kinds.begin(),
            primitive_kinds.end(),
            [](const auto& item) {
                return item.second == ReachPrimitiveKind2::Outer;
            }));
    if (outer_count != 1) {
        throw ReachableArrangementTopologyError(
            "reachable parity requires exactly one outer primitive");
    }
    for (auto face = arrangement.faces_begin();
         face != arrangement.faces_end();
         ++face) {
        face->set_data(ReachFaceState2{});
    }

    ActivePrimitiveParity2 active(primitive_kinds);
    const FaceHandle2 unbounded = arrangement.unbounded_face();
    unbounded->set_data(active.face_state());
    std::vector<ParityFrame2> stack{
        {
            unbounded,
            face_halfedges(unbounded),
            0,
            {},
        },
    };
    while (!stack.empty()) {
        ParityFrame2& frame = stack.back();
        if (frame.next_halfedge == frame.halfedges.size()) {
            const std::vector<std::string> entering =
                std::move(frame.entering_primitives);
            stack.pop_back();
            active.toggle(entering);
            continue;
        }
        const HalfedgeHandle2 crossed =
            frame.halfedges[frame.next_halfedge++];
        const std::vector<std::string>& crossed_primitives =
            crossed->curve().data().primitive_ids;
        active.toggle(crossed_primitives);
        const FaceHandle2 neighbor = crossed->twin()->face();
        const ReachFaceState2 predicted = active.face_state();
        if (!neighbor->data().classified) {
            neighbor->set_data(predicted);
            stack.push_back(
                {
                    neighbor,
                    face_halfedges(neighbor),
                    0,
                    crossed_primitives,
                });
            continue;
        }
        if (!same_face_state(neighbor->data(), predicted)) {
            throw ReachableArrangementTopologyError(
                "dual routes disagree on aggregate face parity");
        }
        active.toggle(crossed_primitives);
    }
    if (!active.empty()) {
        throw ReachableArrangementTopologyError(
            "dual parity traversal did not restore the unbounded state");
    }
    for (auto face = arrangement.faces_begin();
         face != arrangement.faces_end();
         ++face) {
        if (!face->data().classified) {
            throw ReachableArrangementTopologyError(
                "dual parity traversal left a face unclassified");
        }
    }
}

ReachableArrangement2 build_reachable_arrangement(
    CanonicalReachInput2 input)
{
    validate_canonical_radius(input);
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
    classify_faces_by_primitive_parity(
        result.arrangement,
        primitive_kinds);

    for (auto face = result.arrangement.faces_begin();
         face != result.arrangement.faces_end();
         ++face) {
        if (face->data().selected) {
            ++result.audit.selected_face_count;
        }
    }
    for (auto edge = result.arrangement.edges_begin();
         edge != result.arrangement.edges_end();
         ++edge) {
        if (edge->face()->data().selected
            && edge->twin()->face()->data().selected) {
            ++result.audit.selected_adjacency_count;
        }
    }
    const std::size_t component_count =
        selected_component_count(result.arrangement);
    if (component_count == 0) {
        throw PocketNotMachinableError(
            "Phase 1 center domain is empty for the declared tool radius");
    }
    if (component_count != 1) {
        throw PocketNotMachinableError(
            "Phase 1 requires one connected center domain");
    }
    result.center_polygon =
        extract_center_polygon(result.arrangement);
    ++result.audit.center_extractions;
    return result;
}
