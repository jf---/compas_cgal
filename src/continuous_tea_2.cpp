#include "compas.h"
#include "continuous_tea_2/boundary_events.h"
#include "continuous_tea_2/result.h"

#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

namespace {

nb::bytes as_bytes(const std::string& value)
{
    return nb::bytes(value.data(), value.size());
}

nb::list coefficient_bytes(
    const std::vector<std::string>& values)
{
    nb::list result;
    for (const std::string& value : values) {
        result.append(as_bytes(value));
    }
    return result;
}

} // namespace

NB_MODULE(_continuous_tea_2, m)
{
    nb::exception<BoundaryExtractionError> extraction_error(
        m,
        "BoundaryExtractionError");
    nb::exception<DegenerateBoundarySupportError>(
        m,
        "DegenerateBoundarySupportError",
        extraction_error.ptr());
    nb::exception<MissingBoundaryEndpointError>(
        m,
        "MissingBoundaryEndpointError",
        extraction_error.ptr());
    nb::exception<MissingBoundaryIntersectionError>(
        m,
        "MissingBoundaryIntersectionError",
        extraction_error.ptr());
    nb::exception<BoundaryFeatureIndexError>(
        m,
        "BoundaryFeatureIndexError",
        extraction_error.ptr());

    nb::enum_<ContinuousTeaVerdict>(
        m,
        "ContinuousTeaVerdict")
        .value(
            "CERTIFIED",
            ContinuousTeaVerdict::CERTIFIED)
        .value(
            "CAP_EXCEEDED",
            ContinuousTeaVerdict::CAP_EXCEEDED)
        .value(
            "UNRESOLVED_DEGENERACY",
            ContinuousTeaVerdict::UNRESOLVED_DEGENERACY);

    m.attr("BOUNDARY_EVENT_KINDS") = nb::make_tuple(
        "transverse",
        "tangent",
        "vertex",
        "overlap",
        "seam");

    nb::class_<BoundaryFeatureRecord2>(
        m,
        "BoundaryFeatureRecord2")
        .def_prop_ro(
            "support_kind",
            [](const BoundaryFeatureRecord2& record) {
                return record.support_kind;
            })
        .def_prop_ro(
            "support_coefficients",
            [](const BoundaryFeatureRecord2& record) {
                return coefficient_bytes(
                    record.support_coefficients);
            })
        .def_ro(
            "primitive_coefficients",
            &BoundaryFeatureRecord2::primitive_coefficients)
        .def_prop_ro(
            "support_id",
            [](const BoundaryFeatureRecord2& record) {
                return as_bytes(record.support_id);
            })
        .def_prop_ro(
            "source_exact",
            [](const BoundaryFeatureRecord2& record) {
                return as_bytes(record.source_exact);
            })
        .def_prop_ro(
            "target_exact",
            [](const BoundaryFeatureRecord2& record) {
                return as_bytes(record.target_exact);
            })
        .def_prop_ro(
            "source_vertex_id",
            [](const BoundaryFeatureRecord2& record) {
                return as_bytes(record.source_vertex_id);
            })
        .def_prop_ro(
            "target_vertex_id",
            [](const BoundaryFeatureRecord2& record) {
                return as_bytes(record.target_vertex_id);
            })
        .def_prop_ro(
            "material_side",
            [](const BoundaryFeatureRecord2& record) {
                return record.material_side;
            })
        .def_prop_ro(
            "trim_predicate",
            [](const BoundaryFeatureRecord2& record) {
                return as_bytes(record.trim_predicate);
            })
        .def_prop_ro(
            "feature_id",
            [](const BoundaryFeatureRecord2& record) {
                return as_bytes(record.feature_id);
            })
        .def_ro(
            "overlap_multiplicity",
            &BoundaryFeatureRecord2::overlap_multiplicity);

    nb::class_<BoundaryEvent2>(m, "BoundaryEvent2")
        .def_ro("kind", &BoundaryEvent2::kind)
        .def_prop_ro(
            "first_feature_id",
            [](const BoundaryEvent2& event) {
                return as_bytes(event.first_feature_id);
            })
        .def_prop_ro(
            "second_feature_id",
            [](const BoundaryEvent2& event) {
                return as_bytes(event.second_feature_id);
            })
        .def_prop_ro(
            "vertex_id",
            [](const BoundaryEvent2& event) {
                return as_bytes(event.vertex_id);
            })
        .def_prop_ro(
            "overlap_id",
            [](const BoundaryEvent2& event) {
                return as_bytes(event.overlap_id);
            })
        .def_prop_ro(
            "exact_overlap_record",
            [](const BoundaryEvent2& event) {
                return as_bytes(
                    event.exact_overlap_record);
            })
        .def_ro(
            "multiplicity",
            &BoundaryEvent2::multiplicity);

    m.def(
        "extract_boundary_records",
        &extract_boundary_records,
        "stock"_a);
    m.def(
        "classify_boundary_pair",
        &classify_boundary_pair,
        "stock"_a,
        "first_index"_a,
        "second_index"_a);
    m.def(
        "extract_boundary_events",
        &extract_boundary_events,
        "stock"_a);
}
