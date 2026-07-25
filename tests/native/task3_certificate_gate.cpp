#include "task3_certificate_gate.h"

#include "reachable_arrangement_2.h"
#include "reachable_certificate_2.h"
#include "reachable_errors_2.h"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <CGAL/Arrangement_2.h>

namespace {

void require(bool condition, const char* message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

std::size_t decode_u64(
    const std::string& record,
    std::size_t& offset,
    const char* message)
{
    constexpr std::size_t U64_BYTES = 8;
    if (offset > record.size()
        || record.size() - offset < U64_BYTES) {
        throw std::runtime_error(message);
    }
    std::uint64_t value = 0;
    for (std::size_t index = 0; index < U64_BYTES; ++index) {
        value = (value << 8)
            | static_cast<unsigned char>(record[offset + index]);
    }
    offset += U64_BYTES;
    if (value > std::numeric_limits<std::size_t>::max()) {
        throw std::runtime_error(message);
    }
    return static_cast<std::size_t>(value);
}

std::size_t decode_u64_field(
    const std::string& field,
    const char* message)
{
    std::size_t offset = 0;
    const std::size_t value =
        decode_u64(field, offset, message);
    require(offset == field.size(), message);
    return value;
}

std::vector<std::string> record_fields(
    const std::string& record,
    const std::string& tag)
{
    const std::string prefix = tag + '\0';
    if (!record.starts_with(prefix)) {
        throw std::runtime_error(
            "structural record has the wrong tag");
    }
    std::size_t offset = prefix.size();
    const std::size_t count = decode_u64(
        record,
        offset,
        "structural record has no field count");
    constexpr std::size_t LENGTH_PREFIX_BYTES = 8;
    require(
        count <= (record.size() - offset) / LENGTH_PREFIX_BYTES,
        "structural record declares impossible fields");
    std::vector<std::string> fields;
    fields.reserve(count);
    for (std::size_t index = 0; index < count; ++index) {
        const std::size_t size = decode_u64(
            record,
            offset,
            "structural record has no field length");
        if (offset > record.size()
            || size > record.size() - offset) {
            throw std::runtime_error(
                "structural record field exceeds its bytes");
        }
        fields.push_back(record.substr(offset, size));
        offset += size;
    }
    require(
        offset == record.size(),
        "structural record has trailing bytes");
    return fields;
}

compas::RowMatrixXd rectangle_matrix()
{
    compas::RowMatrixXd boundary(4, 3);
    boundary << 0, 0, 0, 10, 0, 0, 10, 8, 0, 0, 8, 0;
    return boundary;
}

CanonicalReachInput2 rectangle_input()
{
    return canonical_reach_input(
        rectangle_matrix(),
        {},
        1.0);
}

std::vector<std::string> source_piece_records(
    const std::vector<std::string>& sources)
{
    std::vector<std::string> pieces;
    pieces.reserve(sources.size());
    for (const std::string& source : sources) {
        pieces.push_back(reach_tagged_record(
            "source-piece-v1",
            {source, reach_u64_record(0)}));
    }
    std::sort(pieces.begin(), pieces.end());
    return pieces;
}

void append_labelled_outer(
    std::vector<ReachDataTraits2::X_monotone_curve_2>& curves,
    const std::vector<std::string>& source_pieces)
{
    const std::vector<ReachKernelPoint> points{
        {0, 0},
        {4, 0},
        {4, 3},
        {0, 3},
    };
    for (std::size_t index = 0; index < points.size(); ++index) {
        curves.push_back({
            ReachXCurve(
                points[index],
                points[(index + 1) % points.size()]),
            ReachCurveLabels2{
                source_pieces,
                {"synthetic-outer"},
            },
        });
    }
}

ReachableArrangement2 classified_synthetic(
    std::vector<std::string> source_records,
    std::vector<ReachDataTraits2::X_monotone_curve_2> curves)
{
    ReachableArrangement2 reachable{
        ReachArrangement2{},
        rectangle_input(),
        std::move(source_records),
        ReachPolygonWithHoles{},
        ReachPolygonWithHoles{},
        {},
    };
    reachable.audit.input_vertex_count = 4;
    CGAL::insert(
        reachable.arrangement,
        curves.begin(),
        curves.end());
    classify_faces_by_primitive_parity(
        reachable.arrangement,
        {{"synthetic-outer", ReachPrimitiveKind2::Outer}});
    for (auto face = reachable.arrangement.faces_begin();
         face != reachable.arrangement.faces_end();
         ++face) {
        reachable.audit.selected_face_count +=
            face->data().selected ? 1 : 0;
    }
    for (auto edge = reachable.arrangement.edges_begin();
         edge != reachable.arrangement.edges_end();
         ++edge) {
        reachable.audit.selected_adjacency_count +=
            edge->face()->data().selected
                && edge->twin()->face()->data().selected
            ? 1
            : 0;
    }
    return reachable;
}

ReachableArrangement2 split_synthetic(
    bool reverse_insertion,
    bool complete_catalog)
{
    const ReachableArrangement2 source_owner =
        build_reachable_arrangement(rectangle_input());
    std::vector<std::string> sources =
        source_owner.source_records;
    if (!complete_catalog) {
        sources.resize(1);
    }
    const std::vector<std::string> pieces =
        source_piece_records(sources);
    std::vector<ReachDataTraits2::X_monotone_curve_2> curves;
    append_labelled_outer(curves, pieces);
    curves.push_back({
        ReachXCurve(
            ReachKernelPoint(2, 0),
            ReachKernelPoint(2, 1)),
        ReachCurveLabels2{pieces, {}},
    });
    curves.push_back({
        ReachXCurve(
            ReachKernelPoint(2, 1),
            ReachKernelPoint(2, 3)),
        ReachCurveLabels2{pieces, {}},
    });
    if (reverse_insertion) {
        std::reverse(curves.begin(), curves.end());
    }
    return classified_synthetic(
        std::move(sources),
        std::move(curves));
}

ReachableArrangement2 self_adjacency_synthetic()
{
    const ReachableArrangement2 source_owner =
        build_reachable_arrangement(rectangle_input());
    const std::vector<std::string> pieces =
        source_piece_records(source_owner.source_records);
    std::vector<ReachDataTraits2::X_monotone_curve_2> curves;
    append_labelled_outer(curves, pieces);
    curves.push_back({
        ReachXCurve(
            ReachKernelPoint(2, 1),
            ReachKernelPoint(2, 2)),
        ReachCurveLabels2{pieces, {}},
    });
    return classified_synthetic(
        source_owner.source_records,
        std::move(curves));
}

ReachableArrangement2 incomplete_catalog_synthetic()
{
    const ReachableArrangement2 source_owner =
        build_reachable_arrangement(rectangle_input());
    std::vector<std::string> sources{
        source_owner.source_records.front(),
    };
    const std::vector<std::string> pieces =
        source_piece_records(sources);
    std::vector<ReachDataTraits2::X_monotone_curve_2> curves;
    append_labelled_outer(curves, pieces);
    return classified_synthetic(
        std::move(sources),
        std::move(curves));
}

template <typename Build>
void require_topology_failure(
    Build build,
    const char* message)
{
    try {
        build();
    }
    catch (const ReachableArrangementTopologyError&) {
        return;
    }
    throw std::runtime_error(message);
}

void incomplete_catalog_gate()
{
    ReachableArrangement2 incomplete =
        incomplete_catalog_synthetic();
    require(
        incomplete.source_records.size() == 1,
        "incomplete-catalog fixture is not incomplete");
    require_topology_failure(
        [&incomplete]() {
            static_cast<void>(
                build_reachable_certificate(
                    incomplete,
                    true));
        },
        "certificate accepted a recipe-incomplete source catalog");
}

void require_vertex_record(const std::string& record)
{
    const std::vector<std::string> fields =
        record_fields(record, "arrangement-vertex-v1");
    require(
        fields.size() == 2,
        "vertex record lost incidence or ordinal");
    const std::vector<std::string> incidences =
        record_fields(
            fields[0],
            "arrangement-incidence-multiset-v1");
    require(
        !incidences.empty(),
        "selected vertex lost source incidence");
    static_cast<void>(
        decode_u64_field(
            fields[1],
            "vertex ordinal is not u64"));
}

void require_cycle_record(
    const std::string& record,
    const std::string& expected_role)
{
    const std::vector<std::string> fields =
        record_fields(
            record,
            "selected-boundary-cycle-v1");
    require(
        fields.size() == 2
            && fields[0] == expected_role,
        "selected cycle lost its boundary role");
    const std::vector<std::string> elements =
        record_fields(fields[1], "cycle-elements-v1");
    require(
        !elements.empty(),
        "selected cycle lost its elements");
    for (const std::string& element : elements) {
        const std::vector<std::string> element_fields =
            record_fields(
                element,
                "selected-cycle-element-v1");
        require(
            element_fields.size() == 4,
            "cycle element lost endpoint, source, or direction evidence");
        const std::vector<std::string> source_endpoint =
            record_fields(
                element_fields[0],
                "source-endpoint-v1");
        const std::vector<std::string> target_endpoint =
            record_fields(
                element_fields[1],
                "target-endpoint-v1");
        require(
            source_endpoint.size() == 1
                && target_endpoint.size() == 1
                && source_endpoint[0] != target_endpoint[0],
            "cycle element has invalid endpoint roles");
        require_vertex_record(source_endpoint[0]);
        require_vertex_record(target_endpoint[0]);
        require(
            !record_fields(
                 element_fields[2],
                 "curve-source-set-v1")
                 .empty(),
            "cycle element lost its source set");
        const std::vector<std::string> direction =
            record_fields(
                element_fields[3],
                "curve-direction-v1");
        require(
            direction.size() == 1
                && (direction[0] == "linear-right"
                    || direction[0] == "linear-left"
                    || direction[0]
                        == "circular-counterclockwise"
                    || direction[0]
                        == "circular-clockwise"),
            "cycle element lost exact direction role");
    }
}

void require_cell_schema(
    const ReachableDomainCertificate2& certificate,
    std::size_t expected_inner_ccbs)
{
    require(
        certificate.selected_cell_records.size() == 1,
        "schema fixture requires one selected cell");
    const std::vector<std::string> fields =
        record_fields(
            certificate.selected_cell_records.front(),
            "selected-cell-v1");
    require(
        fields.size() == 2,
        "selected cell lost outer or inner CCB records");
    require_cycle_record(fields[0], "outer");
    const std::vector<std::string> holes =
        record_fields(fields[1], "selected-cell-holes-v1");
    require(
        holes.size() == expected_inner_ccbs,
        "selected cell omitted an inner CCB");
    for (const std::string& hole : holes) {
        require_cycle_record(hole, "hole");
    }
}

void cell_schema_gate()
{
    ReachableArrangement2 rectangle =
        build_reachable_arrangement(rectangle_input());
    const ReachableDomainCertificate2 rectangle_certificate =
        build_reachable_certificate(rectangle, true);
    require_cell_schema(rectangle_certificate, 0);

    compas::RowMatrixXd outer(4, 3);
    outer << 0, 0, 0, 20, 0, 0, 20, 20, 0, 0, 20, 0;
    compas::RowMatrixXd island(4, 3);
    island << 8, 8, 0, 8, 12, 0, 12, 12, 0, 12, 8, 0;
    ReachableArrangement2 island_reachable =
        build_reachable_arrangement(
            canonical_reach_input(
                outer,
                {island},
                1.0));
    const ReachableDomainCertificate2 island_certificate =
        build_reachable_certificate(
            island_reachable,
            true);
    require_cell_schema(island_certificate, 1);
}

void adjacency_gate()
{
    ReachableArrangement2 reachable =
        split_synthetic(false, true);
    require(
        reachable.audit.selected_face_count == 2
            && reachable.audit.selected_adjacency_count == 2,
        "parallel-adjacency fixture does not contain two shared DCEL edges");
    const ReachableDomainCertificate2 certificate =
        build_reachable_certificate(reachable, true);
    const std::vector<std::string> component =
        record_fields(
            certificate.component_records.front(),
            "selected-component-v1");
    require(
        component.size() == 3,
        "component record lost adjacency structure");
    const std::vector<std::string> cell_ids =
        record_fields(
            component[0],
            "selected-component-cells-v1");
    const std::vector<std::string> adjacencies =
        record_fields(
            component[1],
            "selected-component-adjacencies-v1");
    require(
        certificate.selected_cell_records.size() == 2
            && cell_ids.size() == 2
            && cell_ids[0] != cell_ids[1]
            && adjacencies.size() == 1,
        "parallel edges did not produce two cells and one canonical pair");
    const std::set<std::string> members(
        cell_ids.begin(),
        cell_ids.end());
    const std::vector<std::string> pair =
        record_fields(
            adjacencies.front(),
            "selected-cell-adjacency-v1");
    require(
        pair.size() == 2
            && pair[0] != pair[1]
            && members.contains(pair[0])
            && members.contains(pair[1]),
        "adjacency is self-referential or outside component membership");
    for (std::size_t ordinal = 0;
         ordinal < cell_ids.size();
         ++ordinal) {
        const std::vector<std::string> id =
            record_fields(
                cell_ids[ordinal],
                "selected-cell-id-v1");
        require(
            id.size() == 1
                && decode_u64_field(
                       id[0],
                       "cell ordinal is not u64")
                    == ordinal
                && ordinal
                    < certificate.selected_cell_records.size(),
            "component cell ID does not bind a selected-cell record");
    }
}

void self_adjacency_gate()
{
    ReachableArrangement2 reachable =
        self_adjacency_synthetic();
    require(
        reachable.audit.selected_adjacency_count == 1,
        "self-adjacency fixture has no same-face DCEL edge");
    require_topology_failure(
        [&reachable]() {
            static_cast<void>(
                build_reachable_certificate(
                    reachable,
                    true));
        },
        "certificate accepted a selected self adjacency");
}

void insertion_order_gate()
{
    ReachableArrangement2 forward =
        split_synthetic(false, true);
    ReachableArrangement2 reverse =
        split_synthetic(true, true);
    const ReachableDomainCertificate2 forward_certificate =
        build_reachable_certificate(forward, true);
    const ReachableDomainCertificate2 reverse_certificate =
        build_reachable_certificate(reverse, true);
    require(
        forward_certificate.source_curve_records
                == reverse_certificate.source_curve_records
            && forward_certificate.selected_cell_records
                == reverse_certificate.selected_cell_records
            && forward_certificate.component_records
                == reverse_certificate.component_records,
        "insertion order changed multi-incidence certificate bytes");
}

void source_and_input_gate()
{
    constexpr std::size_t RECTANGLE_SOURCE_COUNT = 12;
    constexpr std::size_t MAX_DIRECT_SOURCE_RECORD_BYTES = 128;
    ReachableArrangement2 rectangle =
        build_reachable_arrangement(rectangle_input());
    const ReachableDomainCertificate2 certificate =
        build_reachable_certificate(rectangle, true);
    require(
        certificate.strategy_version
                == "exact-reachable-arrangement-v2"
            && certificate.source_curve_records.size()
                == RECTANGLE_SOURCE_COUNT
            && certificate.exact_cell_selection
            && certificate.complete_source_provenance
            && certificate.reachable_subset_of_design,
        "certificate lost strategy, source, or exact structural claims");
    std::size_t max_source_record_bytes = 0;
    for (const std::string& source :
         certificate.source_curve_records) {
        const std::vector<std::string> fields =
            record_fields(source, "source-curve-v2");
        require(
            fields.size() == 5
                && fields[0] == "outer"
                && fields[1] == reach_u64_record(0)
                && (fields[3] == "offset-minus"
                    || fields[3] == "offset-plus"
                    || fields[3] == "vertex-circle")
                && source.size()
                    <= MAX_DIRECT_SOURCE_RECORD_BYTES,
            "source record does not use the bounded direct v2 schema");
        max_source_record_bytes =
            std::max(max_source_record_bytes, source.size());
    }
    compas::RowMatrixXd changed = rectangle_matrix();
    changed(1, 0) = 11.0;
    require(
        certificate.matches_exact_inputs(
            rectangle_matrix(), {}, 1.0)
            && !certificate.matches_exact_inputs(
                changed, {}, 1.0)
            && !certificate.matches_exact_inputs(
                rectangle_matrix(), {}, 1.5),
        "certificate exact-input binding accepted a mutation");
    require(
        rectangle.audit.cycle_rotation_comparisons
            <= 3 * rectangle.audit.cycle_element_count,
        "cycle canonicalization exceeded linear comparison bound");
    std::cout
        << "task3-certificate"
        << " sources=" << certificate.source_curve_records.size()
        << " vertices=" << rectangle.arrangement.number_of_vertices()
        << " cells=" << certificate.selected_cell_records.size()
        << " components=" << certificate.component_records.size()
        << " cycle-elements=" << rectangle.audit.cycle_element_count
        << " cycle-comparisons="
        << rectangle.audit.cycle_rotation_comparisons
        << " max-source-bytes=" << max_source_record_bytes
        << '\n';
}

} // namespace

void task3_certificate_gate()
{
    incomplete_catalog_gate();
    cell_schema_gate();
    adjacency_gate();
    self_adjacency_gate();
    insertion_order_gate();
    source_and_input_gate();
}
