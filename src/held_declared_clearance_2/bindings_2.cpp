// nanobind surface of the native declared-boundary clearance oracle
// (declared_clearance_2.h). One call per batch: every element of the batch is
// decided in C++ with the GIL released; per element the result is
// (clearance, primitive, depth) with clearance -1 GOUGE / 0 TOUCH / 1 CLEAR,
// primitive -1 when CLEAR, depth a reported binary64.

#include "held_declared_clearance_2/declared_clearance_2.h"
#include "held_declared_clearance_2/errors_2.h"
#include "held_declared_clearance_2/wire_2.h"

#include <nanobind/nanobind.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/unique_ptr.h>
#include <nanobind/stl/vector.h>

#include <cstddef>
#include <memory>
#include <tuple>
#include <vector>

namespace nb = nanobind;
using namespace nb::literals;

namespace {

namespace hdc = compas_cgal::held_declared_clearance;

using VerdictWire = std::tuple<int, std::ptrdiff_t, double>;

VerdictWire verdict_wire(const hdc::ClearanceVerdict& verdict)
{
    return {static_cast<int>(verdict.clearance), verdict.primitive, verdict.depth};
}

// Reported counters of the last batch (never a decision input).
hdc::ClearanceCounters last_counters;

void require_equal_sizes(std::size_t first, std::size_t second)
{
    if (first != second) throw hdc::InvalidClearanceWireError("coordinate arrays differ in length");
}

}  // namespace

NB_MODULE(_held_declared_clearance_2, m)
{
    m.doc() = "Native exact clearance of Held path elements against the declared (Pythagorean) boundary.";
    auto base = nb::exception<hdc::DeclaredClearanceError>(m, "DeclaredClearanceError", PyExc_ValueError);
    nb::exception<hdc::InvalidDeclaredPrimitiveError>(m, "InvalidDeclaredPrimitiveError", base.ptr());
    nb::exception<hdc::NonPositiveClearanceRadiusError>(m, "NonPositiveClearanceRadiusError", base.ptr());
    nb::exception<hdc::DegenerateClearanceElementError>(m, "DegenerateClearanceElementError", base.ptr());
    nb::exception<hdc::NegativeRadicandError>(m, "NegativeRadicandError", base.ptr());
    nb::exception<hdc::CrossRootArithmeticError>(m, "CrossRootArithmeticError", base.ptr());
    nb::exception<hdc::MixedRootPointError>(m, "MixedRootPointError", base.ptr());
    nb::exception<hdc::InvalidClearanceWireError>(m, "InvalidClearanceWireError", base.ptr());
    nb::exception<hdc::ReportingOverflowError>(m, "ReportingOverflowError", base.ptr());

    m.def("last_counters", [] {
        return std::make_tuple(last_counters.elements, last_counters.pairs, last_counters.owner_touches,
                               last_counters.crossings);
    });

    nb::class_<hdc::DeclaredClearance2>(m, "DeclaredClearance2")
        .def_static(
            "build",
            [](const std::vector<hdc::PrimitiveWire>& primitives, const hdc::RationalWire& tool_radius) {
                std::vector<hdc::Primitive> parsed;
                parsed.reserve(primitives.size());
                for (std::size_t index = 0; index < primitives.size(); ++index) {
                    parsed.push_back(hdc::primitive_from_wire(index, primitives[index]));
                }
                const hdc::FT radius = hdc::rational_from_wire(tool_radius);
                nb::gil_scoped_release release;
                return hdc::DeclaredClearance2::build(std::move(parsed), radius);
            },
            "primitives"_a, "tool_radius"_a)
        .def_prop_ro("primitive_count", &hdc::DeclaredClearance2::primitive_count)
        .def(
            "circles_binary64",
            [](const hdc::DeclaredClearance2& oracle, const std::vector<double>& xs, const std::vector<double>& ys,
               const std::vector<double>& radii) {
                require_equal_sizes(xs.size(), ys.size());
                require_equal_sizes(xs.size(), radii.size());
                std::vector<hdc::CircleElement> circles;
                circles.reserve(xs.size());
                for (std::size_t i = 0; i < xs.size(); ++i) {
                    circles.push_back({hdc::RationalPoint(hdc::binary64_to_exact(xs[i]), hdc::binary64_to_exact(ys[i])),
                                       hdc::binary64_to_exact(radii[i])});
                }
                std::vector<VerdictWire> out;
                out.reserve(circles.size());
                hdc::ClearanceWorkspace workspace;
                {
                    nb::gil_scoped_release release;
                    for (const hdc::CircleElement& circle : circles) out.push_back(verdict_wire(oracle.circle(circle, workspace)));
                }
                last_counters = workspace.counters;
                return out;
            },
            "xs"_a, "ys"_a, "radii"_a)
        .def(
            "segments_binary64",
            [](const hdc::DeclaredClearance2& oracle, const std::vector<double>& x0, const std::vector<double>& y0,
               const std::vector<double>& x1, const std::vector<double>& y1) {
                require_equal_sizes(x0.size(), y0.size());
                require_equal_sizes(x0.size(), x1.size());
                require_equal_sizes(x0.size(), y1.size());
                std::vector<hdc::SegmentElement> segments;
                segments.reserve(x0.size());
                for (std::size_t i = 0; i < x0.size(); ++i) {
                    segments.push_back(hdc::SegmentElement::between(
                        hdc::RationalPoint(hdc::binary64_to_exact(x0[i]), hdc::binary64_to_exact(y0[i])),
                        hdc::RationalPoint(hdc::binary64_to_exact(x1[i]), hdc::binary64_to_exact(y1[i]))));
                }
                const std::vector<std::size_t> no_owners;
                std::vector<VerdictWire> out;
                out.reserve(segments.size());
                hdc::ClearanceWorkspace workspace;
                {
                    nb::gil_scoped_release release;
                    for (const hdc::SegmentElement& segment : segments) {
                        out.push_back(verdict_wire(oracle.segment(segment, no_owners, workspace)));
                    }
                }
                last_counters = workspace.counters;
                return out;
            },
            "x0"_a, "y0"_a, "x1"_a, "y1"_a)
        .def(
            "elements",
            [](const hdc::DeclaredClearance2& oracle, const std::vector<hdc::ElementWire>& wires) {
                std::vector<hdc::DecodedElement> elements;
                elements.reserve(wires.size());
                for (const hdc::ElementWire& wire : wires) {
                    elements.push_back(hdc::element_from_wire(wire, oracle.primitive_count()));
                }
                std::vector<VerdictWire> out;
                out.reserve(elements.size());
                hdc::ClearanceWorkspace workspace;
                {
                    nb::gil_scoped_release release;
                    for (const hdc::DecodedElement& element : elements) {
                        const hdc::ClearanceVerdict verdict = element.is_arc
                                                                  ? oracle.arc(element.arc, element.owners, workspace)
                                                                  : oracle.segment(element.segment, element.owners, workspace);
                        out.push_back(verdict_wire(verdict));
                    }
                }
                last_counters = workspace.counters;
                return out;
            },
            "elements"_a);
}
