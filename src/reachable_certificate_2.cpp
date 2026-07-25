#include "reachable_certificate_2.h"

#include "canonical_rotation_2.h"
#include "reachable_errors_2.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <map>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <CGAL/enum.h>

namespace {

constexpr std::string_view REACHABLE_STRATEGY_VERSION =
    "exact-reachable-arrangement-v2";
constexpr std::size_t MINIMAL_ROTATION_COMPARISON_FACTOR = 3;

using ReachHalfedgeConstHandle2 =
    ReachArrangement2::Halfedge_const_handle;

template <typename Handle>
const void* handle_address(const Handle& handle)
{
    return static_cast<const void*>(std::addressof(*handle));
}

std::size_t decode_u64(
    const std::string& bytes,
    std::size_t& offset,
    std::string_view context)
{
    constexpr std::size_t U64_BYTES = 8;
    if (offset > bytes.size()
        || bytes.size() - offset < U64_BYTES) {
        throw ReachableArrangementTopologyError(
            std::string(context)
            + " has a truncated unsigned field");
    }
    std::uint64_t value = 0;
    for (std::size_t index = 0; index < U64_BYTES; ++index) {
        value = (value << 8)
            | static_cast<unsigned char>(bytes[offset + index]);
    }
    offset += U64_BYTES;
    if (value > std::numeric_limits<std::size_t>::max()) {
        throw ReachableArrangementTopologyError(
            std::string(context)
            + " exceeds the native structural range");
    }
    return static_cast<std::size_t>(value);
}

std::size_t decode_u64_field(
    const std::string& bytes,
    std::string_view context)
{
    std::size_t offset = 0;
    const std::size_t value = decode_u64(
        bytes,
        offset,
        context);
    if (offset != bytes.size()) {
        throw ReachableArrangementTopologyError(
            std::string(context)
            + " has trailing bytes");
    }
    return value;
}

std::vector<std::string> record_fields(
    const std::string& record,
    std::string_view tag)
{
    std::string prefix(tag);
    prefix.push_back('\0');
    if (!record.starts_with(prefix)) {
        throw ReachableArrangementTopologyError(
            std::string(tag)
            + " structural record has the wrong tag");
    }
    std::size_t offset = prefix.size();
    const std::size_t count = decode_u64(
        record,
        offset,
        tag);
    constexpr std::size_t LENGTH_PREFIX_BYTES = 8;
    if (count
        > (record.size() - offset) / LENGTH_PREFIX_BYTES) {
        throw ReachableArrangementTopologyError(
            std::string(tag)
            + " structural record declares impossible fields");
    }
    std::vector<std::string> fields;
    fields.reserve(count);
    for (std::size_t index = 0; index < count; ++index) {
        const std::size_t size = decode_u64(
            record,
            offset,
            tag);
        if (offset > record.size()
            || size > record.size() - offset) {
            throw ReachableArrangementTopologyError(
                std::string(tag)
                + " structural field exceeds its record");
        }
        fields.push_back(record.substr(offset, size));
        offset += size;
    }
    if (offset != record.size()) {
        throw ReachableArrangementTopologyError(
            std::string(tag)
            + " structural record has trailing bytes");
    }
    return fields;
}

struct SourceContract2 {
    std::size_t piece_count;
};

using SourceContracts2 =
    std::unordered_map<std::string, SourceContract2>;

std::size_t expected_source_record_count(
    const ReachableArrangement2& reachable)
{
    std::size_t input_vertex_count =
        reachable.input.outer.points.size();
    for (const CanonicalReachRing2& hole :
         reachable.input.holes) {
        if (input_vertex_count
            > std::numeric_limits<std::size_t>::max()
                - hole.points.size()) {
            throw ReachableArrangementTopologyError(
                "certificate input vertex count overflows");
        }
        input_vertex_count += hole.points.size();
    }
    if (input_vertex_count
        != reachable.audit.input_vertex_count) {
        throw ReachableArrangementTopologyError(
            "certificate input vertex count contradicts the canonical recipe");
    }
    constexpr std::size_t SOURCE_ROLES_PER_VERTEX = 3;
    if (input_vertex_count
        > std::numeric_limits<std::size_t>::max()
            / SOURCE_ROLES_PER_VERTEX) {
        throw ReachableArrangementTopologyError(
            "certificate source-record count overflows");
    }
    return SOURCE_ROLES_PER_VERTEX * input_vertex_count;
}

const CanonicalReachRing2& source_ring(
    const ReachableArrangement2& reachable,
    const std::vector<std::string>& fields)
{
    const std::size_t ring_ordinal =
        decode_u64_field(fields[1], "source ring ordinal");
    if (fields[0] == "outer") {
        if (ring_ordinal != 0) {
            throw ReachableArrangementTopologyError(
                "outer source record has a nonzero ring ordinal");
        }
        return reachable.input.outer;
    }
    if (fields[0] == "hole"
        && ring_ordinal < reachable.input.holes.size()) {
        return reachable.input.holes[ring_ordinal];
    }
    throw ReachableArrangementTopologyError(
        "source record references an unknown canonical ring");
}

SourceContracts2 validate_source_records(
    const ReachableArrangement2& reachable)
{
    if (reachable.source_records.size()
            != expected_source_record_count(reachable)
        || !std::is_sorted(
            reachable.source_records.begin(),
            reachable.source_records.end())
        || std::adjacent_find(
               reachable.source_records.begin(),
               reachable.source_records.end())
            != reachable.source_records.end()) {
        throw ReachableArrangementTopologyError(
            "source catalog cardinality, order, or uniqueness "
            "contradicts the recipe");
    }
    SourceContracts2 contracts;
    contracts.reserve(reachable.source_records.size());
    for (const std::string& source : reachable.source_records) {
        const std::vector<std::string> fields =
            record_fields(source, "source-curve-v2");
        if (fields.size() != 5) {
            throw ReachableArrangementTopologyError(
                "source record does not match the direct v2 schema");
        }
        const CanonicalReachRing2& ring =
            source_ring(reachable, fields);
        const std::size_t feature_ordinal =
            decode_u64_field(
                fields[2],
                "source feature ordinal");
        if (feature_ordinal >= ring.points.size()) {
            throw ReachableArrangementTopologyError(
                "source record references an unknown ring feature");
        }
        if (fields[4]
            != reach_binary64_record(
                reachable.input.binary64_radius)) {
            throw ReachableArrangementTopologyError(
                "source record radius diverges from canonical input");
        }
        std::size_t piece_count = 0;
        if (fields[3] == "offset-minus"
            || fields[3] == "offset-plus") {
            piece_count = 1;
        }
        else if (fields[3] == "vertex-circle") {
            piece_count = 2;
        }
        else {
            throw ReachableArrangementTopologyError(
                "source record has an unknown construction role");
        }
        if (!contracts.emplace(
                source,
                SourceContract2{piece_count})
                 .second) {
            throw ReachableArrangementTopologyError(
                "source record identity collision");
        }
    }
    return contracts;
}

void validate_source_piece(
    const std::string& source_piece,
    const SourceContracts2& contracts,
    std::unordered_set<std::string>& referenced_sources)
{
    const std::vector<std::string> fields =
        record_fields(source_piece, "source-piece-v1");
    if (fields.size() != 2) {
        throw ReachableArrangementTopologyError(
            "source-piece record has the wrong field count");
    }
    const auto source = contracts.find(fields[0]);
    if (source == contracts.end()) {
        throw ReachableArrangementTopologyError(
            "arrangement source piece is absent from the source catalog");
    }
    if (decode_u64_field(
            fields[1],
            "source-piece ordinal")
        >= source->second.piece_count) {
        throw ReachableArrangementTopologyError(
            "arrangement source piece has an invalid ordinal");
    }
    referenced_sources.insert(source->first);
}

bool validate_complete_source_provenance(
    const ReachableArrangement2& reachable,
    const SourceContracts2& contracts)
{
    std::unordered_set<std::string> referenced_sources;
    referenced_sources.reserve(contracts.size());
    for (auto edge = reachable.arrangement.edges_begin();
         edge != reachable.arrangement.edges_end();
         ++edge) {
        const std::vector<std::string>& source_piece_ids =
            edge->curve().data().source_piece_ids;
        if (!std::is_sorted(
                source_piece_ids.begin(),
                source_piece_ids.end())
            || std::adjacent_find(
                   source_piece_ids.begin(),
                   source_piece_ids.end())
                != source_piece_ids.end()) {
            throw ReachableArrangementTopologyError(
                "arrangement source-piece labels must be sorted and unique");
        }
        for (const std::string& source_piece :
             source_piece_ids) {
            validate_source_piece(
                source_piece,
                contracts,
                referenced_sources);
        }
    }
    if (referenced_sources.size() != contracts.size()) {
        throw ReachableArrangementTopologyError(
            "source catalog contains an unbound source curve");
    }
    return true;
}

struct ReachPointLess2 {
    bool operator()(
        const ReachPoint& left,
        const ReachPoint& right) const
    {
        return ReachTraits().compare_xy_2_object()(
                   left,
                   right)
            == CGAL::SMALLER;
    }
};

using ReachVertexRecordMap2 =
    std::map<ReachPoint, std::string, ReachPointLess2>;

ReachVertexRecordMap2 build_vertex_records(
    const ReachArrangement2& arrangement)
{
    std::map<std::string, std::vector<ReachPoint>>
        points_by_incidence;
    for (auto vertex = arrangement.vertices_begin();
         vertex != arrangement.vertices_end();
         ++vertex) {
        std::vector<std::string> incidences;
        auto first = vertex->incident_halfedges();
        if (first != nullptr) {
            auto halfedge = first;
            do {
                const std::vector<std::string>& source_piece_ids =
                    halfedge->curve().data().source_piece_ids;
                incidences.insert(
                    incidences.end(),
                    source_piece_ids.begin(),
                    source_piece_ids.end());
            } while (++halfedge != first);
        }
        std::sort(incidences.begin(), incidences.end());
        points_by_incidence[reach_tagged_record(
            "arrangement-incidence-multiset-v1",
            incidences)]
            .push_back(vertex->point());
    }

    ReachVertexRecordMap2 vertex_records;
    for (auto& [incidence, points] : points_by_incidence) {
        std::sort(
            points.begin(),
            points.end(),
            ReachPointLess2{});
        for (std::size_t ordinal = 0;
             ordinal < points.size();
             ++ordinal) {
            if (!vertex_records.emplace(
                    points[ordinal],
                    reach_tagged_record(
                        "arrangement-vertex-v1",
                        {
                            incidence,
                            reach_u64_record(ordinal),
                        }))
                     .second) {
                throw ReachableArrangementTopologyError(
                    "arrangement contains a duplicate exact vertex");
            }
        }
    }
    if (vertex_records.size()
        != arrangement.number_of_vertices()) {
        throw ReachableArrangementTopologyError(
            "vertex records do not cover the complete arrangement");
    }
    return vertex_records;
}

const std::string& vertex_record(
    const ReachVertexRecordMap2& records,
    const ReachPoint& point)
{
    const auto found = records.find(point);
    if (found == records.end()) {
        throw ReachableArrangementTopologyError(
            "selected cycle references an unknown exact vertex");
    }
    return found->second;
}

ReachXCurve directed_halfedge_curve(
    ReachHalfedgeConstHandle2 halfedge)
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
            "certificate halfedge does not match its exact DCEL endpoints");
    }
    return curve;
}

std::string exact_direction_record(
    const ReachXCurve& curve)
{
    if (curve.is_linear()) {
        return reach_tagged_record(
            "curve-direction-v1",
            {
                curve.is_directed_right()
                    ? "linear-right"
                    : "linear-left",
            });
    }
    return reach_tagged_record(
        "curve-direction-v1",
        {
            curve.orientation() == CGAL::COUNTERCLOCKWISE
                ? "circular-counterclockwise"
                : "circular-clockwise",
        });
}

std::vector<std::string> canonical_cycle_rotation(
    const std::vector<std::string>& elements,
    ReachableDomainBuildAudit2& audit)
{
    if (elements.empty()) {
        throw ReachableArrangementTopologyError(
            "selected boundary cycle must be nonempty");
    }
    if (audit.cycle_element_count
        > std::numeric_limits<std::size_t>::max()
            - elements.size()) {
        throw ReachableArrangementTopologyError(
            "certificate cycle-element count overflows");
    }
    audit.cycle_element_count += elements.size();
    std::size_t comparisons = 0;
    const std::size_t first = minimal_rotation_index(
        elements,
        [](const std::string& left, const std::string& right) {
            if (left < right) {
                return CGAL::SMALLER;
            }
            if (right < left) {
                return CGAL::LARGER;
            }
            return CGAL::EQUAL;
        },
        comparisons);
    if (audit.cycle_rotation_comparisons
        > std::numeric_limits<std::size_t>::max()
            - comparisons) {
        throw ReachableArrangementTopologyError(
            "certificate cycle-comparison count overflows");
    }
    audit.cycle_rotation_comparisons += comparisons;
    std::vector<std::string> rotated;
    rotated.reserve(elements.size());
    rotated.insert(
        rotated.end(),
        elements.begin()
            + static_cast<std::ptrdiff_t>(first),
        elements.end());
    rotated.insert(
        rotated.end(),
        elements.begin(),
        elements.begin()
            + static_cast<std::ptrdiff_t>(first));
    return rotated;
}

std::string cycle_record(
    const std::vector<ReachHalfedgeConstHandle2>& halfedges,
    std::string_view role,
    const ReachVertexRecordMap2& vertex_records,
    ReachableDomainBuildAudit2& audit)
{
    std::vector<std::string> elements;
    elements.reserve(halfedges.size());
    for (const ReachHalfedgeConstHandle2 halfedge : halfedges) {
        const std::vector<std::string>& source_piece_ids =
            halfedge->curve().data().source_piece_ids;
        if (source_piece_ids.empty()) {
            throw ReachableArrangementTopologyError(
                "selected cell edge has no propagated source provenance");
        }
        const ReachXCurve curve =
            directed_halfedge_curve(halfedge);
        elements.push_back(reach_tagged_record(
            "selected-cycle-element-v1",
            {
                reach_tagged_record(
                    "source-endpoint-v1",
                    {vertex_record(
                        vertex_records,
                        halfedge->source()->point())}),
                reach_tagged_record(
                    "target-endpoint-v1",
                    {vertex_record(
                        vertex_records,
                        halfedge->target()->point())}),
                reach_tagged_record(
                    "curve-source-set-v1",
                    source_piece_ids),
                exact_direction_record(curve),
            }));
    }
    return reach_tagged_record(
        "selected-boundary-cycle-v1",
        {
            std::string(role),
            reach_tagged_record(
                "cycle-elements-v1",
                canonical_cycle_rotation(elements, audit)),
        });
}

template <typename Circulator>
std::vector<ReachHalfedgeConstHandle2> ccb_halfedges(
    Circulator first)
{
    std::vector<ReachHalfedgeConstHandle2> halfedges;
    Circulator halfedge = first;
    do {
        halfedges.push_back(halfedge);
    } while (++halfedge != first);
    return halfedges;
}

std::string selected_face_record(
    ReachArrangement2::Face_const_handle face,
    const ReachVertexRecordMap2& vertex_records,
    ReachableDomainBuildAudit2& audit)
{
    if (!face->data().selected || face->is_unbounded()) {
        throw ReachableArrangementTopologyError(
            "selected-cell record requires a bounded selected face");
    }
    const std::string outer = cycle_record(
        ccb_halfedges(face->outer_ccb()),
        "outer",
        vertex_records,
        audit);
    std::vector<std::string> holes;
    holes.reserve(face->number_of_inner_ccbs());
    for (auto inner = face->inner_ccbs_begin();
         inner != face->inner_ccbs_end();
         ++inner) {
        holes.push_back(cycle_record(
            ccb_halfedges(*inner),
            "hole",
            vertex_records,
            audit));
    }
    std::sort(holes.begin(), holes.end());
    return reach_tagged_record(
        "selected-cell-v1",
        {
            outer,
            reach_tagged_record(
                "selected-cell-holes-v1",
                holes),
        });
}

bool exact_cell_selection(
    const ReachableArrangement2& reachable)
{
    std::size_t selected_count = 0;
    for (auto face = reachable.arrangement.faces_begin();
         face != reachable.arrangement.faces_end();
         ++face) {
        const ReachFaceState2& state = face->data();
        if (!state.classified
            || state.selected
                != (state.outer_active
                    && state.active_holes == 0
                    && state.active_forbidden == 0)) {
            throw ReachableArrangementTopologyError(
                "certificate face selection contradicts exact parity");
        }
        selected_count += state.selected ? 1 : 0;
    }
    if (selected_count
        != reachable.audit.selected_face_count) {
        throw ReachableArrangementTopologyError(
            "certificate selected-face count contradicts the build audit");
    }
    return true;
}

ReachHalfedgeConstHandle2 next_selected_boundary(
    ReachHalfedgeConstHandle2 current,
    std::size_t halfedge_count)
{
    ReachHalfedgeConstHandle2 candidate = current->next();
    for (std::size_t step = 0;
         step < halfedge_count;
         ++step) {
        if (!candidate->face()->data().selected) {
            throw ReachableArrangementTopologyError(
                "component boundary traversal left the selected face set");
        }
        if (!candidate->twin()->face()->data().selected) {
            return candidate;
        }
        candidate = candidate->twin()->next();
    }
    throw ReachableArrangementTopologyError(
        "component boundary did not close around a vertex");
}

std::vector<std::string> component_boundary_cycles(
    const ReachArrangement2& arrangement,
    const ReachVertexRecordMap2& vertex_records,
    ReachableDomainBuildAudit2& audit)
{
    const std::size_t halfedge_count =
        arrangement.number_of_halfedges();
    std::unordered_set<const void*> visited;
    visited.reserve(halfedge_count);
    std::vector<std::string> records;
    std::size_t selected_boundary_count = 0;
    for (auto candidate = arrangement.halfedges_begin();
         candidate != arrangement.halfedges_end();
         ++candidate) {
        const bool selected_boundary =
            candidate->face()->data().selected
            && !candidate->twin()->face()->data().selected;
        if (!selected_boundary) {
            continue;
        }
        ++selected_boundary_count;
        const ReachHalfedgeConstHandle2 first = candidate;
        if (visited.contains(handle_address(first))) {
            continue;
        }
        std::vector<ReachHalfedgeConstHandle2> halfedges;
        ReachPolygon polygon;
        ReachHalfedgeConstHandle2 current = first;
        for (std::size_t step = 0;
             step < halfedge_count;
             ++step) {
            if (!visited.insert(handle_address(current)).second) {
                throw ReachableArrangementTopologyError(
                    "component boundary revisited a halfedge before closure");
            }
            halfedges.push_back(current);
            polygon.push_back(
                directed_halfedge_curve(current));
            current = next_selected_boundary(
                current,
                halfedge_count);
            if (current == first) {
                break;
            }
        }
        if (current != first || halfedges.empty()) {
            throw ReachableArrangementTopologyError(
                "component boundary cycle did not close");
        }
        const CGAL::Orientation orientation =
            polygon.orientation();
        std::string_view role;
        if (orientation == CGAL::COUNTERCLOCKWISE) {
            role = "outer";
        }
        else if (orientation == CGAL::CLOCKWISE) {
            role = "hole";
        }
        else {
            throw ReachableArrangementTopologyError(
                "component boundary cycle has zero exact area");
        }
        records.push_back(cycle_record(
            halfedges,
            role,
            vertex_records,
            audit));
    }
    if (visited.size() != selected_boundary_count
        || records.empty()) {
        throw ReachableArrangementTopologyError(
            "component boundary records do not cover the selected boundary");
    }
    std::sort(records.begin(), records.end());
    return records;
}

std::vector<std::string> adjacency_records(
    const ReachArrangement2& arrangement,
    const std::unordered_map<const void*, std::string>&
        cell_id_by_face,
    const std::size_t expected_selected_adjacency_edge_count)
{
    std::vector<std::pair<std::string, std::string>>
        face_pairs;
    std::size_t selected_adjacency_edge_count = 0;
    for (auto edge = arrangement.edges_begin();
         edge != arrangement.edges_end();
         ++edge) {
        if (!edge->face()->data().selected
            || !edge->twin()->face()->data().selected) {
            continue;
        }
        ++selected_adjacency_edge_count;
        if (edge->face() == edge->twin()->face()) {
            throw ReachableArrangementTopologyError(
                "selected-cell adjacency cannot reference one face twice");
        }
        const auto left = cell_id_by_face.find(
            handle_address(edge->face()));
        const auto right = cell_id_by_face.find(
            handle_address(edge->twin()->face()));
        if (left == cell_id_by_face.end()
            || right == cell_id_by_face.end()) {
            throw ReachableArrangementTopologyError(
                "selected adjacency references an unknown cell");
        }
        if (left->second == right->second) {
            throw ReachableArrangementTopologyError(
                "distinct selected faces share one cell identity");
        }
        face_pairs.emplace_back(
            std::min(left->second, right->second),
            std::max(left->second, right->second));
    }
    if (selected_adjacency_edge_count
        != expected_selected_adjacency_edge_count) {
        throw ReachableArrangementTopologyError(
            "selected adjacency edge count contradicts the build audit");
    }
    std::sort(face_pairs.begin(), face_pairs.end());
    face_pairs.erase(
        std::unique(face_pairs.begin(), face_pairs.end()),
        face_pairs.end());
    std::vector<std::string> records;
    records.reserve(face_pairs.size());
    for (const auto& [left, right] : face_pairs) {
        records.push_back(reach_tagged_record(
            "selected-cell-adjacency-v1",
            {left, right}));
    }
    return records;
}

} // namespace

bool ReachableDomainCertificate2::matches_exact_inputs(
    Eigen::Ref<const compas::RowMatrixXd> boundary,
    const std::vector<compas::RowMatrixXd>& holes,
    double tool_radius) const
{
    try {
        return canonical_reach_input(
                   boundary,
                   holes,
                   tool_radius)
                .recipe_record
            == input_recipe_record;
    }
    catch (const InvalidReachableDomainInputError&) {
        return false;
    }
}

ReachableDomainCertificate2 build_reachable_certificate(
    ReachableArrangement2& reachable,
    bool reachable_subset_of_design)
{
    validate_canonical_reach_input(reachable.input);
    const SourceContracts2 source_contracts =
        validate_source_records(reachable);
    const bool complete_provenance =
        validate_complete_source_provenance(
            reachable,
            source_contracts);
    const ReachVertexRecordMap2 vertex_records =
        build_vertex_records(reachable.arrangement);
    std::vector<std::string> arrangement_vertex_records;
    arrangement_vertex_records.reserve(vertex_records.size());
    for (const auto& [point, record] : vertex_records) {
        static_cast<void>(point);
        arrangement_vertex_records.push_back(record);
    }
    std::sort(
        arrangement_vertex_records.begin(),
        arrangement_vertex_records.end());
    if (arrangement_vertex_records.size()
            != reachable.arrangement.number_of_vertices()
        || std::adjacent_find(
               arrangement_vertex_records.begin(),
               arrangement_vertex_records.end())
            != arrangement_vertex_records.end()) {
        throw ReachableArrangementTopologyError(
            "certificate vertex records are incomplete or nonunique");
    }

    std::vector<std::string> selected_cell_records;
    selected_cell_records.reserve(
        reachable.audit.selected_face_count);
    std::unordered_map<const void*, std::string>
        cell_record_by_face;
    cell_record_by_face.reserve(
        reachable.audit.selected_face_count);
    for (auto face = reachable.arrangement.faces_begin();
         face != reachable.arrangement.faces_end();
         ++face) {
        if (!face->data().selected) {
            continue;
        }
        std::string record = selected_face_record(
            face,
            vertex_records,
            reachable.audit);
        if (!cell_record_by_face.emplace(
                handle_address(face),
                record)
                 .second) {
            throw ReachableArrangementTopologyError(
                "selected face has a duplicate ephemeral identity");
        }
        selected_cell_records.push_back(std::move(record));
    }
    std::sort(
        selected_cell_records.begin(),
        selected_cell_records.end());
    if (selected_cell_records.size()
            != reachable.audit.selected_face_count
        || std::adjacent_find(
               selected_cell_records.begin(),
               selected_cell_records.end())
            != selected_cell_records.end()) {
        throw ReachableArrangementTopologyError(
            "certificate cell records are incomplete or nonunique");
    }
    std::unordered_map<std::string, std::string>
        cell_id_by_record;
    cell_id_by_record.reserve(selected_cell_records.size());
    std::vector<std::string> selected_cell_ids;
    selected_cell_ids.reserve(selected_cell_records.size());
    for (std::size_t ordinal = 0;
         ordinal < selected_cell_records.size();
         ++ordinal) {
        std::string id = reach_tagged_record(
            "selected-cell-id-v1",
            {reach_u64_record(ordinal)});
        cell_id_by_record.emplace(
            selected_cell_records[ordinal],
            id);
        selected_cell_ids.push_back(std::move(id));
    }
    for (auto& [address, record] : cell_record_by_face) {
        static_cast<void>(address);
        const auto id = cell_id_by_record.find(record);
        if (id == cell_id_by_record.end()) {
            throw ReachableArrangementTopologyError(
                "selected face record has no canonical cell identity");
        }
        record = id->second;
    }

    const std::vector<std::string> adjacencies =
        adjacency_records(
            reachable.arrangement,
            cell_record_by_face,
            reachable.audit.selected_adjacency_count);
    const std::vector<std::string> boundary_cycles =
        component_boundary_cycles(
            reachable.arrangement,
            vertex_records,
            reachable.audit);
    if (reachable.audit.cycle_element_count
            <= std::numeric_limits<std::size_t>::max()
                / MINIMAL_ROTATION_COMPARISON_FACTOR
        && reachable.audit.cycle_rotation_comparisons
            > MINIMAL_ROTATION_COMPARISON_FACTOR
                * reachable.audit.cycle_element_count) {
        throw ReachableArrangementTopologyError(
            "cycle canonicalization exceeded its linear comparison bound");
    }

    std::vector<std::string> component_records{
        reach_tagged_record(
            "selected-component-v1",
            {
                reach_tagged_record(
                    "selected-component-cells-v1",
                    selected_cell_ids),
                reach_tagged_record(
                    "selected-component-adjacencies-v1",
                    adjacencies),
                reach_tagged_record(
                    "selected-component-boundary-cycles-v1",
                    boundary_cycles),
            }),
    };
    return {
        std::string(REACHABLE_STRATEGY_VERSION),
        reachable.source_records,
        std::move(arrangement_vertex_records),
        std::move(selected_cell_records),
        std::move(component_records),
        exact_cell_selection(reachable),
        complete_provenance,
        reachable_subset_of_design,
        reachable.input.recipe_record,
    };
}
