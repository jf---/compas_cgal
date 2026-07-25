#include "exact_region_2.h"
#include "exact_sweep_2.h"
#include "reachable_arrangement_2.h"
#include "reachable_certificate_2.h"
#include "reachable_errors_2.h"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <CGAL/Arrangement_2.h>
#include <CGAL/Boolean_set_operations_2/Gps_polygon_validation.h>
#include <CGAL/iterator.h>
#include <CGAL/number_utils.h>

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

std::vector<std::string> record_fields(
    const std::string& record,
    const std::string& tag)
{
    const std::string prefix = tag + '\0';
    if (!record.starts_with(prefix)) {
        throw std::runtime_error("structural record has the wrong tag");
    }
    std::size_t offset = prefix.size();
    const std::size_t count = decode_u64(
        record,
        offset,
        "structural record has no field count");
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
    if (offset != record.size()) {
        throw std::runtime_error(
            "structural record has trailing bytes");
    }
    return fields;
}
void exact_region_storage_gate()
{
    ReachSet disk;
    disk.insert(reach_disk_polygon(ReachKernelPoint(0, 0), ReachFT(2)));
    const ExactRegion2 original = ExactRegion2::build(
        std::move(disk),
        ExactRegionRole2::ReachableMaterial,
        "native-gate-region");
    const ExactRegion2 clone = original.clone();
    require(
        original.shares_storage_with_for_audit(clone),
        "read-only clone deep-copied exact storage");
    require(clone.contains(2.0, 0.0), "closed disk lost exact boundary");
}

void require_counterclockwise_parts(
    const std::vector<ReachPolygon>& parts,
    const char* message)
{
    for (const ReachPolygon& part : parts) {
        require(part.orientation() == CGAL::COUNTERCLOCKWISE, message);
    }
}

void irrational_capsule_gate()
{
    const ReachFT root_thirteen = CGAL::sqrt(ReachFT(13));
    const std::vector<ReachPolygon> capsule = reach_capsule_parts(
        ReachKernelPoint(0, 0),
        ReachKernelPoint(2, 3),
        ReachFT(5));
    require_counterclockwise_parts(
        capsule,
        "irrational capsule emitted non-CCW standalone polygon");
    const ReachSet capsule_set = reach_join_parts(capsule, {});
    const ReachPoint irrational_side_boundary(
        ReachFT(1) - ReachFT(15) / root_thirteen,
        ReachFT(3) / ReachFT(2) + ReachFT(10) / root_thirteen);
    require(
        capsule_set.oriented_side(irrational_side_boundary)
            == CGAL::ON_ORIENTED_BOUNDARY,
        "sqrt(13) capsule lost exact irrational side boundary");
}

ReachXCurve quarter_guide_arc(
    CGAL::Orientation orientation,
    const ReachFT& radius)
{
    const ReachKernelPoint center(0, 0);
    const ReachKernelPoint start(radius, ReachFT(0));
    const ReachKernelPoint end = orientation == CGAL::COUNTERCLOCKWISE
        ? ReachKernelPoint(ReachFT(0), radius)
        : ReachKernelPoint(ReachFT(0), -radius);
    const ReachKernel::Circle_2 circle(
        center,
        radius * radius,
        orientation);
    const ReachCurve curve(
        circle,
        ReachPoint(start.x(), start.y()),
        ReachPoint(end.x(), end.y()));
    std::vector<ReachXCurve> pieces;
    ReachTraits().make_x_monotone_2_object()(
        curve,
        CGAL::dispatch_or_drop_output<ReachXCurve>(
            std::back_inserter(pieces)));
    require(
        pieces.size() == 1,
        "quarter guide arc did not remain x-monotone");
    return pieces.front();
}

void arc_sweep_case(
    CGAL::Orientation orientation,
    const ReachFT& guide_radius,
    const ReachFT& tool_radius,
    CGAL::Oriented_side expected_center_side)
{
    const ReachXCurve guide =
        quarter_guide_arc(orientation, guide_radius);
    const std::vector<ReachPolygon> parts =
        reach_arc_sweep_parts(guide, tool_radius);
    require(parts.size() == 3, "arc sweep lost primitive parts");
    require_counterclockwise_parts(
        parts,
        "arc sweep emitted non-CCW standalone polygon");
    const ReachSet sweep = reach_join_parts(parts, {});
    require(
        reach_component_count(sweep) == 1,
        "arc sweep did not produce one connected component");
    require(
        sweep.oriented_side(ReachPoint(ReachFT(0), ReachFT(0)))
            == expected_center_side,
        "arc sweep center topology contradicts radius branch");
    require(
        sweep.oriented_side(
            ReachPoint(guide_radius + tool_radius, ReachFT(0)))
            == CGAL::ON_ORIENTED_BOUNDARY,
        "arc sweep lost exact outer endpoint boundary");
}

void arc_sweep_orientation_and_radius_gate()
{
    arc_sweep_case(
        CGAL::COUNTERCLOCKWISE,
        ReachFT(2),
        ReachFT(1),
        CGAL::ON_NEGATIVE_SIDE);
    arc_sweep_case(
        CGAL::CLOCKWISE,
        ReachFT(2),
        ReachFT(1),
        CGAL::ON_NEGATIVE_SIDE);
    arc_sweep_case(
        CGAL::COUNTERCLOCKWISE,
        ReachFT(1),
        ReachFT(1),
        CGAL::ON_ORIENTED_BOUNDARY);
    arc_sweep_case(
        CGAL::CLOCKWISE,
        ReachFT(1),
        ReachFT(1),
        CGAL::ON_ORIENTED_BOUNDARY);
    arc_sweep_case(
        CGAL::COUNTERCLOCKWISE,
        ReachFT(1),
        ReachFT(2),
        CGAL::ON_POSITIVE_SIDE);
    arc_sweep_case(
        CGAL::CLOCKWISE,
        ReachFT(1),
        ReachFT(2),
        CGAL::ON_POSITIVE_SIDE);
}

void full_circle_case(
    const ReachKernelVector& phase_vector,
    const ReachFT& tool_radius,
    const ReachFT& outer_radius,
    std::size_t expected_hole_count,
    CGAL::Oriented_side expected_center_side)
{
    const ReachSet sweep = reach_full_circle_sweep(
        ReachKernelPoint(0, 0),
        phase_vector,
        tool_radius);
    std::vector<ReachPolygonWithHoles> components;
    sweep.polygons_with_holes(std::back_inserter(components));
    require(
        components.size() == 1,
        "full-circle sweep did not produce one component");
    require(
        static_cast<std::size_t>(std::distance(
            components.front().holes_begin(),
            components.front().holes_end()))
            == expected_hole_count,
        "full-circle sweep hole count contradicts radius branch");
    require(
        sweep.oriented_side(ReachPoint(ReachFT(0), ReachFT(0)))
            == expected_center_side,
        "full-circle sweep center topology contradicts radius branch");
    require(
        sweep.oriented_side(ReachPoint(outer_radius, ReachFT(0)))
            == CGAL::ON_ORIENTED_BOUNDARY,
        "full-circle sweep lost exact outer boundary");
}

void full_circle_radius_gate()
{
    const ReachSet annulus = reach_full_circle_sweep(
        ReachKernelPoint(0, 0),
        ReachKernelVector(2, 0),
        ReachFT(1));
    require(
        annulus.oriented_side(ReachPoint(ReachFT(1), ReachFT(0)))
            == CGAL::ON_ORIENTED_BOUNDARY,
        "full-circle annulus lost exact inner boundary");
    full_circle_case(
        ReachKernelVector(2, 0),
        ReachFT(1),
        ReachFT(3),
        1,
        CGAL::ON_NEGATIVE_SIDE);
    full_circle_case(
        ReachKernelVector(1, 0),
        ReachFT(1),
        ReachFT(2),
        0,
        CGAL::ON_POSITIVE_SIDE);
    full_circle_case(
        ReachKernelVector(1, 0),
        ReachFT(2),
        ReachFT(3),
        0,
        CGAL::ON_POSITIVE_SIDE);
}

compas::RowMatrixXd rectangle_matrix()
{
    compas::RowMatrixXd boundary(4, 3);
    boundary << 0, 0, 0, 10, 0, 0, 10, 8, 0, 0, 8, 0;
    return boundary;
}

CanonicalReachInput2 rectangle_input()
{
    return canonical_reach_input(rectangle_matrix(), {}, 1.0);
}

void require_all_faces_classified(const ReachableArrangement2& reachable)
{
    for (auto face = reachable.arrangement.faces_begin();
         face != reachable.arrangement.faces_end();
         ++face) {
        require(face->data().classified, "reachable face remained unclassified");
    }
}

void require_selected_boundary_provenance(
    const ReachableArrangement2& reachable)
{
    std::size_t boundary_count = 0;
    for (auto halfedge = reachable.arrangement.halfedges_begin();
         halfedge != reachable.arrangement.halfedges_end();
         ++halfedge) {
        const bool selected_boundary =
            halfedge->face()->data().selected
            && !halfedge->twin()->face()->data().selected;
        if (!selected_boundary) {
            continue;
        }
        ++boundary_count;
        require(
            !halfedge->curve().data().source_piece_ids.empty(),
            "selected boundary lost source-piece provenance");
    }
    require(boundary_count != 0, "reachable domain exposed no selected boundary");
}

void require_reachable_audit(
    const ReachableArrangement2& reachable,
    std::size_t source_record_count)
{
    require(
        reachable.audit.provenance_arrangements == 1,
        "reachable build constructed multiple arrangements");
    require(
        reachable.audit.center_extractions == 1,
        "reachable build extracted the center domain multiple times");
    require(
        reachable.audit.source_geometric_rematches == 0,
        "reachable build rematched sources geometrically");
    require(
        reachable.source_records.size() == source_record_count,
        "reachable build emitted the wrong source-record count");
    require(
        reachable.audit.ring_rotation_comparisons
            <= 3 * reachable.audit.input_vertex_count,
        "ring canonicalization exceeded linear comparison bound");
    std::size_t propagated_primitive_labels = 0;
    for (auto edge = reachable.arrangement.edges_begin();
         edge != reachable.arrangement.edges_end();
         ++edge) {
        propagated_primitive_labels +=
            edge->curve().data().primitive_ids.size();
    }
    const std::size_t face_count =
        reachable.arrangement.number_of_faces();
    const std::size_t halfedge_count =
        reachable.arrangement.number_of_halfedges();
    require(
        reachable.audit.dense_face_visits == face_count,
        "dense index revisited arrangement faces");
    require(
        reachable.audit.dense_halfedge_visits == halfedge_count,
        "dense index revisited arrangement halfedges");
    require(
        reachable.audit.primitive_label_resolutions
            == propagated_primitive_labels,
        "dense index did not resolve each primitive label exactly once");
    require(
        reachable.audit.parity_halfedge_visits == halfedge_count,
        "parity traversal did not perform exactly one directed halfedge pass");
    require(
        reachable.audit.component_halfedge_visits <= halfedge_count,
        "component traversal exceeded one directed halfedge pass");
    require(
        reachable.audit.boundary_scan_halfedge_visits == halfedge_count,
        "boundary extraction rescanned arrangement halfedges");
    require(
        reachable.audit.boundary_cycle_halfedge_visits <= halfedge_count,
        "boundary extraction revisited selected boundary halfedges");
    require(
        reachable.audit.boundary_rotation_halfedge_visits <= halfedge_count,
        "boundary successor search exceeded one directed halfedge pass");
    require_all_faces_classified(reachable);
    require_selected_boundary_provenance(reachable);
    const ReachTraits traits;
    require(
        CGAL::is_valid_polygon_with_holes<ReachTraits>(
            reachable.design_polygon,
            traits),
        "canonical design polygon is invalid");
    require(
        CGAL::is_valid_polygon_with_holes<ReachTraits>(
            reachable.center_polygon,
            traits),
        "selected center polygon is invalid");
}

void rectangle_reachable_gate()
{
    const ReachableArrangement2 rectangle =
        build_reachable_arrangement(rectangle_input());
    require_reachable_audit(rectangle, 12);

    ReachSet center;
    center.insert(rectangle.center_polygon);
    require(
        center.oriented_side(ReachPoint(ReachFT(5), ReachFT(4)))
            == CGAL::ON_POSITIVE_SIDE,
        "rectangle center domain lost its exact interior");
    require(
        center.oriented_side(ReachPoint(ReachFT(1), ReachFT(1)))
            == CGAL::ON_ORIENTED_BOUNDARY,
        "rectangle center domain lost exact tangency");
    require(
        center.oriented_side(ReachPoint(ReachFT(0), ReachFT(0)))
            == CGAL::ON_NEGATIVE_SIDE,
        "rectangle center domain retained a forbidden corner");
}

void canonical_input_invariance_gate()
{
    compas::RowMatrixXd rotated(4, 3);
    rotated << 10, 8, 0, 0, 8, 0, 0, 0, 0, 10, 0, 0;
    compas::RowMatrixXd reversed(4, 3);
    reversed << 0, 0, 0, 0, 8, 0, 10, 8, 0, 10, 0, 0;
    const CanonicalReachInput2 reference = rectangle_input();
    require(
        canonical_reach_input(rotated, {}, 1.0).recipe_record
            == reference.recipe_record,
        "rotated input changed canonical reach identity");
    require(
        canonical_reach_input(reversed, {}, 1.0).recipe_record
            == reference.recipe_record,
        "reversed input changed canonical reach identity");
}

void acute_corner_gate()
{
    compas::RowMatrixXd boundary(3, 3);
    boundary << 0, 0, 0, 12, 0, 0, 1, 7, 0;
    ReachableArrangement2 acute = build_reachable_arrangement(
        canonical_reach_input(boundary, {}, 0.5));
    require_reachable_audit(acute, 9);
    const ReachableDomainCertificate2 certificate =
        build_reachable_certificate(acute, true);
    require(
        certificate.source_curve_records.size() == 9
            && certificate.selected_cell_records.size()
                == acute.audit.selected_face_count
            && certificate.component_records.size() == 1,
        "acute certificate lost structural coverage");
}

void island_gate()
{
    compas::RowMatrixXd boundary(4, 3);
    boundary << 0, 0, 0, 20, 0, 0, 20, 20, 0, 0, 20, 0;
    compas::RowMatrixXd island(4, 3);
    island << 8, 8, 0, 8, 12, 0, 12, 12, 0, 12, 8, 0;
    ReachableArrangement2 reachable = build_reachable_arrangement(
        canonical_reach_input(boundary, {island}, 1.0));
    require_reachable_audit(reachable, 24);
    require(
        static_cast<std::size_t>(std::distance(
            reachable.center_polygon.holes_begin(),
            reachable.center_polygon.holes_end()))
            == 1,
        "island reachable domain lost its exact hole");
    const ReachableDomainCertificate2 certificate =
        build_reachable_certificate(reachable, true);
    const std::vector<std::string> component =
        record_fields(
            certificate.component_records.front(),
            "selected-component-v1");
    const std::vector<std::string> boundary_cycles =
        record_fields(
            component[2],
            "selected-component-boundary-cycles-v1");
    require(
        certificate.source_curve_records.size() == 24
            && certificate.selected_cell_records.size()
                == reachable.audit.selected_face_count
            && certificate.component_records.size() == 1
            && boundary_cycles.size() == 2,
        "island certificate lost cells, sources, or boundary cycles");
}

void overlap_label_gate()
{
    ReachArrangement2 arrangement;
    std::vector<ReachDataTraits2::X_monotone_curve_2> curves{
        {
            ReachXCurve(
                ReachKernelPoint(0, 0),
                ReachKernelPoint(4, 0)),
            ReachCurveLabels2{
                {"source-b", "source-a", "source-a"},
                {"primitive-b", "primitive-a"},
            },
        },
        {
            ReachXCurve(
                ReachKernelPoint(1, 0),
                ReachKernelPoint(3, 0)),
            ReachCurveLabels2{
                {"source-c", "source-b"},
                {"primitive-c", "primitive-a"},
            },
        },
    };
    CGAL::insert(arrangement, curves.begin(), curves.end());
    bool found_overlap = false;
    for (auto edge = arrangement.edges_begin();
         edge != arrangement.edges_end();
         ++edge) {
        const ReachCurveLabels2& labels = edge->curve().data();
        if (labels.source_piece_ids
            == std::vector<std::string>{
                "source-a",
                "source-b",
                "source-c",
            }) {
            found_overlap = true;
            require(
                labels.primitive_ids
                    == std::vector<std::string>{
                        "primitive-a",
                        "primitive-b",
                        "primitive-c",
                    },
                "overlap merge lost sorted unique primitive labels");
        }
    }
    require(
        found_overlap,
        "overlap merge lost sorted unique source-piece labels");
}

void append_labelled_rectangle(
    std::vector<ReachDataTraits2::X_monotone_curve_2>& curves,
    const std::vector<std::string>& primitive_ids,
    const std::vector<std::string>& source_piece_ids = {})
{
    const std::vector<ReachKernelPoint> points{
        {0, 0},
        {4, 0},
        {4, 3},
        {0, 3},
    };
    for (std::size_t index = 0; index < points.size(); ++index) {
        curves.push_back(
            {
                ReachXCurve(
                    points[index],
                    points[(index + 1) % points.size()]),
                ReachCurveLabels2{
                    source_piece_ids,
                    primitive_ids,
                },
            });
    }
}

void require_identical_boundary_parity(
    bool add_forbidden,
    bool expected_bounded_selected,
    std::size_t expected_bounded_forbidden)
{
    ReachArrangement2 arrangement;
    std::vector<ReachDataTraits2::X_monotone_curve_2> curves;
    append_labelled_rectangle(curves, {"outer"});
    append_labelled_rectangle(curves, {"outer"});
    ReachPrimitiveKinds2 primitive_kinds{
        {"outer", ReachPrimitiveKind2::Outer},
    };
    if (add_forbidden) {
        append_labelled_rectangle(curves, {"forbidden"});
        append_labelled_rectangle(curves, {"forbidden"});
        primitive_kinds.emplace(
            "forbidden",
            ReachPrimitiveKind2::Forbidden);
    }
    CGAL::insert(arrangement, curves.begin(), curves.end());
    classify_faces_by_primitive_parity(
        arrangement,
        primitive_kinds);

    std::size_t bounded_faces = 0;
    for (auto face = arrangement.faces_begin();
         face != arrangement.faces_end();
         ++face) {
        if (face == arrangement.unbounded_face()) {
            require(
                !face->data().outer_active
                    && face->data().active_forbidden == 0
                    && !face->data().selected,
                "identical boundary changed unbounded geometric parity");
            continue;
        }
        ++bounded_faces;
        require(
            face->data().outer_active,
            "identical outer boundaries cancelled geometric interior");
        require(
            face->data().active_holes == 0,
            "identical boundary invented hole parity");
        require(
            face->data().active_forbidden
                == expected_bounded_forbidden,
            "identical forbidden boundaries changed geometric parity");
        require(
            face->data().selected == expected_bounded_selected,
            "identical boundary selected the wrong geometric face");
    }
    require(
        bounded_faces == 1,
        "identical rectangle boundaries changed exact face count");
}

void identical_boundary_parity_gate()
{
    require_identical_boundary_parity(true, false, 1);
    require_identical_boundary_parity(false, true, 0);
}

void disconnected_domain_gate()
{
    compas::RowMatrixXd boundary(12, 3);
    boundary
        << 0, 0, 0,
        4, 0, 0,
        4, 1.5, 0,
        8, 1.5, 0,
        8, 0, 0,
        12, 0, 0,
        12, 4, 0,
        8, 4, 0,
        8, 2.5, 0,
        4, 2.5, 0,
        4, 4, 0,
        0, 4, 0;
    bool raised = false;
    try {
        static_cast<void>(build_reachable_arrangement(
            canonical_reach_input(boundary, {}, 0.75)));
    }
    catch (const PocketNotMachinableError&) {
        raised = true;
    }
    require(
        raised,
        "disconnected center domain did not raise PocketNotMachinableError");
}

template <typename Mutate>
void require_invalid_canonical_mutation(
    Mutate mutate,
    const char* message)
{
    CanonicalReachInput2 input = rectangle_input();
    mutate(input);
    try {
        static_cast<void>(
            build_reachable_arrangement(std::move(input)));
    }
    catch (const InvalidReachableDomainInputError&) {
        return;
    }
    catch (...) {
        throw std::runtime_error(message);
    }
    throw std::runtime_error(message);
}

void canonical_input_bypass_gate()
{
    try {
        static_cast<void>(build_reachable_arrangement(
            CanonicalReachInput2{}));
        throw std::runtime_error(
            "reachable build accepted a raw default input");
    }
    catch (const InvalidReachableDomainInputError&) {
    }
    require_invalid_canonical_mutation(
        [](CanonicalReachInput2& input) {
            input.radius = ReachFT(0);
        },
        "reachable build leaked invalid radius past its input error boundary");
    require_invalid_canonical_mutation(
        [](CanonicalReachInput2& input) {
            input.recipe_record += "forged";
        },
        "reachable build accepted a forged input recipe");
    require_invalid_canonical_mutation(
        [](CanonicalReachInput2& input) {
            input.outer.record += "forged";
        },
        "reachable build accepted a forged ring record");
    require_invalid_canonical_mutation(
        [](CanonicalReachInput2& input) {
            input.outer.points[1] = ReachKernelPoint(11, 0);
        },
        "reachable build accepted exact-point/binary64 divergence");
    require_invalid_canonical_mutation(
        [](CanonicalReachInput2& input) {
            input.outer.canonical_ordinal = 1;
        },
        "reachable build accepted an invalid outer ordinal");
    require_invalid_canonical_mutation(
        [](CanonicalReachInput2& input) {
            input.outer.outer = false;
        },
        "reachable build leaked an invalid ring role to topology");
    require_invalid_canonical_mutation(
        [](CanonicalReachInput2& input) {
            input.outer.points[1] = ReachKernelPoint(11, 0);
            input.outer.points[2] = ReachKernelPoint(11, 8);
            input.outer.binary64_points[1][0] = 11.0;
            input.outer.binary64_points[2][0] = 11.0;
        },
        "reachable build accepted valid geometry detached from identity");
}

void selected_adjacency_certificate_gate()
{
    const ReachableArrangement2 source_owner =
        build_reachable_arrangement(rectangle_input());
    const std::string source = source_owner.source_records.front();
    const std::string piece = reach_tagged_record(
        "source-piece-v1", {source, reach_u64_record(0)});
    ReachableArrangement2 adjacent{
        ReachArrangement2{},
        rectangle_input(),
        {source},
        ReachPolygonWithHoles{}, ReachPolygonWithHoles{}, {},
    };
    std::vector<ReachDataTraits2::X_monotone_curve_2> curves;
    append_labelled_rectangle(curves, {"synthetic-outer"}, {piece});
    curves.push_back({
        ReachXCurve(ReachKernelPoint(2, 0), ReachKernelPoint(2, 3)),
        ReachCurveLabels2{{piece}, {}},
    });
    CGAL::insert(adjacent.arrangement, curves.begin(), curves.end());
    classify_faces_by_primitive_parity(
        adjacent.arrangement,
        {{"synthetic-outer", ReachPrimitiveKind2::Outer}});
    for (auto face = adjacent.arrangement.faces_begin();
         face != adjacent.arrangement.faces_end();
         ++face) {
        adjacent.audit.selected_face_count +=
            face->data().selected ? 1 : 0;
    }
    for (auto edge = adjacent.arrangement.edges_begin();
         edge != adjacent.arrangement.edges_end();
         ++edge) {
        adjacent.audit.selected_adjacency_count +=
            edge->face()->data().selected
                && edge->twin()->face()->data().selected
            ? 1
            : 0;
    }
    const ReachableDomainCertificate2 certificate =
        build_reachable_certificate(adjacent, true);
    const std::vector<std::string> component = record_fields(
        certificate.component_records.front(), "selected-component-v1");
    const std::vector<std::string> adjacencies = record_fields(
        component[1], "selected-component-adjacencies-v1");
    require(
        certificate.selected_cell_records.size() == 2
            && adjacent.audit.selected_adjacency_count == 1
            && adjacencies.size() == 1,
        "certificate did not bind one exact selected-cell adjacency");
}

void certificate_gate()
{
    constexpr std::size_t RECTANGLE_SOURCE_COUNT = 12;
    constexpr std::size_t MAX_DIRECT_SOURCE_RECORD_BYTES = 128;
    ReachableArrangement2 rectangle =
        build_reachable_arrangement(rectangle_input());
    const ReachableDomainCertificate2 certificate =
        build_reachable_certificate(rectangle, true);

    require(
        certificate.strategy_version
            == "exact-reachable-arrangement-v2",
        "certificate strategy version is not structural v2");
    require(
        certificate.source_curve_records.size()
            == RECTANGLE_SOURCE_COUNT,
        "certificate omitted rectangle source curves");
    require(
        certificate.selected_cell_records.size()
            == rectangle.audit.selected_face_count,
        "certificate omitted selected arrangement cells");
    require(
        certificate.component_records.size() == 1,
        "certificate did not bind the selected component");
    require(
        certificate.exact_cell_selection,
        "certificate lost exact face selection");
    require(
        certificate.complete_source_provenance,
        "certificate lost propagated source provenance");
    require(
        certificate.reachable_subset_of_design,
        "certificate lost its supplied exact subset decision");
    require(
        certificate.input_recipe_record
            == rectangle.input.recipe_record,
        "certificate lost canonical input identity");
    require(
        std::is_sorted(
            certificate.source_curve_records.begin(),
            certificate.source_curve_records.end())
            && std::is_sorted(
                certificate.selected_cell_records.begin(),
                certificate.selected_cell_records.end())
            && std::is_sorted(
                certificate.component_records.begin(),
                certificate.component_records.end()),
        "certificate record collections are not canonical");

    std::size_t max_source_record_bytes = 0;
    for (const std::string& source :
         certificate.source_curve_records) {
        const std::vector<std::string> fields =
            record_fields(source, "source-curve-v2");
        require(
            fields.size() == 5,
            "source record does not use the direct v2 schema");
        require(
            fields[0] == "outer"
                && fields[1] == reach_u64_record(0),
            "rectangle source record lost canonical ring identity");
        require(
            fields[3] == "offset-minus"
                || fields[3] == "offset-plus"
                || fields[3] == "vertex-circle",
            "source record has an unknown construction role");
        require(
            source.size() <= MAX_DIRECT_SOURCE_RECORD_BYTES,
            "source record embeds nonlocal ring data");
        max_source_record_bytes =
            std::max(max_source_record_bytes, source.size());
    }

    const std::vector<std::string> component_fields =
        record_fields(
            certificate.component_records.front(),
            "selected-component-v1");
    require(
        component_fields.size() == 3,
        "component record lost cells, adjacency, or boundary cycles");
    const std::vector<std::string> cells =
        record_fields(
            component_fields[0],
            "selected-component-cells-v1");
    const std::vector<std::string> adjacencies =
        record_fields(
            component_fields[1],
            "selected-component-adjacencies-v1");
    const std::vector<std::string> boundary_cycles =
        record_fields(
            component_fields[2],
            "selected-component-boundary-cycles-v1");
    require(
        cells.size() == certificate.selected_cell_records.size(),
        "component does not bind every selected cell exactly once");
    for (std::size_t ordinal = 0;
         ordinal < cells.size();
         ++ordinal) {
        const std::vector<std::string> id_fields =
            record_fields(cells[ordinal], "selected-cell-id-v1");
        require(
            id_fields.size() == 1
                && id_fields[0] == reach_u64_record(ordinal),
            "component cell identity is not canonical");
    }
    require(
        adjacencies.size()
            == rectangle.audit.selected_adjacency_count,
        "component omitted selected-cell adjacency");
    require(
        std::is_sorted(adjacencies.begin(), adjacencies.end())
            && std::adjacent_find(
                adjacencies.begin(),
                adjacencies.end())
                == adjacencies.end(),
        "component adjacency is repeated or noncanonical");
    require(
        boundary_cycles.size() == 1,
        "rectangle component lost its regularized boundary");
    require(
        rectangle.audit.cycle_rotation_comparisons
            <= 3 * rectangle.audit.cycle_element_count,
        "cycle canonicalization exceeded linear comparison bound");

    compas::RowMatrixXd changed_boundary = rectangle_matrix();
    changed_boundary(1, 0) = 11.0;
    require(
        certificate.matches_exact_inputs(
            rectangle_matrix(),
            {},
            1.0),
        "certificate rejected its exact inputs");
    require(
        !certificate.matches_exact_inputs(
            changed_boundary,
            {},
            1.0),
        "certificate accepted changed exact geometry");
    require(
        !certificate.matches_exact_inputs(
            rectangle_matrix(),
            {},
            1.5),
        "certificate accepted a changed exact radius");

    compas::RowMatrixXd rotated(4, 3);
    rotated << 10, 8, 0, 0, 8, 0, 0, 0, 0, 10, 0, 0;
    compas::RowMatrixXd reversed(4, 3);
    reversed << 0, 0, 0, 0, 8, 0, 10, 8, 0, 10, 0, 0;
    require(
        certificate.matches_exact_inputs(rotated, {}, 1.0)
            && certificate.matches_exact_inputs(
                reversed,
                {},
                1.0),
        "certificate exact-input match is not canonical");
    ReachableArrangement2 rotated_reachable =
        build_reachable_arrangement(
            canonical_reach_input(rotated, {}, 1.0));
    ReachableArrangement2 reversed_reachable =
        build_reachable_arrangement(
            canonical_reach_input(reversed, {}, 1.0));
    const ReachableDomainCertificate2 rotated_certificate =
        build_reachable_certificate(rotated_reachable, true);
    const ReachableDomainCertificate2 reversed_certificate =
        build_reachable_certificate(reversed_reachable, true);
    require(
        certificate.source_curve_records
                == rotated_certificate.source_curve_records
            && certificate.source_curve_records
                == reversed_certificate.source_curve_records
            && certificate.selected_cell_records
                == rotated_certificate.selected_cell_records
            && certificate.selected_cell_records
                == reversed_certificate.selected_cell_records
            && certificate.component_records
                == rotated_certificate.component_records
            && certificate.component_records
                == reversed_certificate.component_records,
        "canonical ring variants changed structural certificate");

    ReachableArrangement2 missing_source =
        build_reachable_arrangement(rectangle_input());
    missing_source.source_records.pop_back();
    try {
        static_cast<void>(
            build_reachable_certificate(missing_source, true));
        throw std::runtime_error(
            "certificate accepted an unbound source piece");
    }
    catch (const ReachableArrangementTopologyError&) {
    }

    std::cout
        << "task3-certificate"
        << " sources=" << certificate.source_curve_records.size()
        << " vertices=" << rectangle.arrangement.number_of_vertices()
        << " cells=" << certificate.selected_cell_records.size()
        << " adjacencies=" << adjacencies.size()
        << " components=" << certificate.component_records.size()
        << " boundary-cycles=" << boundary_cycles.size()
        << " cycle-elements=" << rectangle.audit.cycle_element_count
        << " cycle-comparisons="
        << rectangle.audit.cycle_rotation_comparisons
        << " max-source-bytes=" << max_source_record_bytes
        << '\n';
}
} // namespace

int main()
{
    exact_region_storage_gate();
    irrational_capsule_gate();
    arc_sweep_orientation_and_radius_gate();
    full_circle_radius_gate();
    rectangle_reachable_gate();
    canonical_input_invariance_gate();
    acute_corner_gate();
    island_gate();
    overlap_label_gate();
    identical_boundary_parity_gate();
    disconnected_domain_gate();
    canonical_input_bypass_gate();
    selected_adjacency_certificate_gate();
    certificate_gate();
}
