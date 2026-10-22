// nanobind surface of the exact Held centre domain (held_centre_domain_2.h).
//
// Rationals cross as (numerator, denominator) decimal text pairs, one-root
// numbers as (rational, coefficient, radicand) triples of those pairs. The
// binary64 batches inject their doubles exactly; `place_binary64` returns the
// correctly rounded or exactly rounded-down doubles its placement decided
// (held_centre_domain_placement_2.h), never an approximation that decides.

#include "gps_ccb_repair_2.h"
#include "held_centre_domain_2.h"

#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/unique_ptr.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace {
namespace hcd = compas_cgal::held_centre_domain;

using RationalWire = std::pair<std::string, std::string>;
// (kind, sense, rationals): a line is (0, 0, [sx, sy, ex, ey, nx, ny]); an arc is
// (1, sense, [sx, sy, ex, ey, cx, cy, radius]).
using PrimitiveWire = std::tuple<int, int, std::vector<RationalWire>>;

constexpr std::size_t LINE_RATIONALS = 6;
constexpr std::size_t ARC_RATIONALS = 7;

hcd::FT rational(const RationalWire& wire) { return hcd::rational_from_text(wire.first, wire.second); }

hcd::Primitive primitive_from_wire(const PrimitiveWire& wire)
{
    const auto& [kind, sense, values] = wire;
    hcd::Primitive primitive{hcd::PrimitiveKind::Line, {}, {}, {}, {}, hcd::FT(0), 0};
    if (kind == static_cast<int>(hcd::PrimitiveKind::Line)) {
        if (values.size() != LINE_RATIONALS) throw hcd::InvalidCentreDomainInputError("a line needs 6 rationals");
        primitive.normal = hcd::Kernel::Vector_2(rational(values[4]), rational(values[5]));
    } else if (kind == static_cast<int>(hcd::PrimitiveKind::Arc)) {
        if (values.size() != ARC_RATIONALS) throw hcd::InvalidCentreDomainInputError("an arc needs 7 rationals");
        primitive.kind = hcd::PrimitiveKind::Arc;
        primitive.centre = hcd::RationalPoint(rational(values[4]), rational(values[5]));
        primitive.radius = rational(values[6]);
        primitive.sense = sense;
    } else {
        throw hcd::InvalidCentreDomainInputError("unknown primitive kind");
    }
    primitive.start = hcd::RationalPoint(rational(values[0]), rational(values[1]));
    primitive.end = hcd::RationalPoint(rational(values[2]), rational(values[3]));
    return primitive;
}

nb::tuple rational_wire(const hcd::RationalText& text) { return nb::make_tuple(text.numerator, text.denominator); }

nb::tuple one_root_wire(const hcd::OneRootText& text)
{
    return nb::make_tuple(rational_wire(text.rational), rational_wire(text.coefficient), rational_wire(text.radicand));
}

nb::tuple point_wire(const hcd::PointText& text) { return nb::make_tuple(one_root_wire(text.x), one_root_wire(text.y)); }

nb::tuple piece_wire(const hcd::Piece& piece)
{
    return nb::make_tuple(piece.circular,
                          nb::make_tuple(rational_wire(piece.support[0]), rational_wire(piece.support[1]),
                                         rational_wire(piece.support[2])),
                          piece.sense, point_wire(piece.source), point_wire(piece.target));
}

nb::tuple one_root_value_wire(const hcd::OneRoot& value)
{
    return nb::make_tuple(rational_wire(hcd::rational_to_text(value.a)), rational_wire(hcd::rational_to_text(value.b)),
                          rational_wire(hcd::rational_to_text(value.r)));
}

nb::tuple nearest_wire(const hcd::NearestDecision& decision, const hcd::NearestValues& values)
{
    return nb::make_tuple(decision.ambiguous ? 1 : 0, static_cast<int>(decision.kind), decision.piece,
                          one_root_value_wire(values.squared_distance), one_root_value_wire(values.x),
                          one_root_value_wire(values.y));
}

// Status codes of a batch projection (held_centre_domain_native.py mirrors them).
constexpr int PROJECTED = 0;
constexpr int AMBIGUOUS = 1;
constexpr int BEYOND_GAP = 2;

// Queries of the last batch the interval stage could not decide (reporting only).
std::size_t last_exact_fallbacks = 0;

hcd::FT binary64(double value)
{
    if (!std::isfinite(value)) throw hcd::InvalidCentreDomainInputError("query coordinates must be finite binary64");
    return hcd::FT(value);  // exact: every finite binary64 is a dyadic rational
}

nb::list components_wire(const hcd::CentreDomain2& domain)
{
    nb::list components;
    for (const hcd::Component& component : domain.components()) {
        nb::list cycles;
        for (const hcd::Cycle& cycle : component.cycles) {
            nb::list pieces;
            for (const hcd::Piece& piece : cycle) pieces.append(piece_wire(piece));
            cycles.append(pieces);
        }
        components.append(cycles);
    }
    return components;
}

}  // namespace

NB_MODULE(_held_centre_domain_2, m)
{
    m.doc() = "Exact Held centre domain Omega* (docs/held_exact_centre_domain.md); CORE-free.";
    nb::exception<hcd::InvalidCentreDomainInputError>(m, "InvalidCentreDomainInputError", PyExc_ValueError);
    nb::exception<hcd::BrokenCentreDomainChainError>(m, "BrokenCentreDomainChainError", PyExc_RuntimeError);
    nb::exception<hcd::NoRationalUnitNormalError>(m, "NoRationalUnitNormalError", PyExc_RuntimeError);
    // The Gps CCB repair found the domain arrangement corrupt beyond the misfile
    // it repairs (gps_ccb_repair_2.h): broken invariants, never argument faults.
    nb::exception<AmbiguousOuterCcbError<hcd::CentreDomainCcbRepairDomain>>(
        m, "AmbiguousOuterCcbError", PyExc_RuntimeError);
    nb::exception<BoundedFaceWithoutOuterCcbError<hcd::CentreDomainCcbRepairDomain>>(
        m, "BoundedFaceWithoutOuterCcbError", PyExc_RuntimeError);
    nb::exception<UnboundedFaceOuterCcbError<hcd::CentreDomainCcbRepairDomain>>(
        m, "UnboundedFaceOuterCcbError", PyExc_RuntimeError);

    m.def("last_exact_fallbacks", [] { return last_exact_fallbacks; });
    nb::class_<hcd::CentreDomain2>(m, "CentreDomain2")
        .def_static(
            "build",
            [](const std::vector<PrimitiveWire>& primitives, const RationalWire& tool_radius) {
                std::vector<hcd::Primitive> parsed;
                parsed.reserve(primitives.size());
                for (const PrimitiveWire& wire : primitives) parsed.push_back(primitive_from_wire(wire));
                const hcd::FT radius = rational(tool_radius);
                nb::gil_scoped_release release;
                return hcd::CentreDomain2::build(parsed, radius);
            },
            "primitives"_a, "tool_radius"_a)
        .def("components", &components_wire)
        .def_prop_ro("member_count", &hcd::CentreDomain2::member_count)
        .def_prop_ro("construction_seconds", &hcd::CentreDomain2::construction_seconds)
        .def_prop_ro("refiled_outer_ccbs", &hcd::CentreDomain2::refiled_outer_ccbs)
        .def(
            "locate",
            [](const hcd::CentreDomain2& domain, const RationalWire& x, const RationalWire& y) {
                const hcd::RationalPoint point(rational(x), rational(y));
                nb::gil_scoped_release release;
                return static_cast<int>(domain.locate(point));
            },
            "x"_a, "y"_a)
        .def(
            "nearest",
            [](const hcd::CentreDomain2& domain, const RationalWire& x, const RationalWire& y) {
                const hcd::RationalPoint point(rational(x), rational(y));
                hcd::NearestDecision decision{};
                hcd::NearestValues values{};
                {
                    nb::gil_scoped_release release;
                    hcd::QueryWorkspace workspace;
                    decision = domain.decide_nearest(point, nullptr, workspace);
                    values = domain.nearest_values(point, decision);
                }
                return nearest_wire(decision, values);
            },
            "x"_a, "y"_a)
        .def(
            "project_binary64",
            // Batch contact decisions for binary64 points (injected exactly):
            // per point (status, kind, piece). Coordinates are materialised on
            // demand through `nearest`.
            [](const hcd::CentreDomain2& domain, const std::vector<double>& xs, const std::vector<double>& ys,
               const RationalWire& gap) {
                if (xs.size() != ys.size()) throw hcd::InvalidCentreDomainInputError("x and y counts differ");
                const hcd::FT gap_value = rational(gap);
                const hcd::FT gap_squared = gap_value * gap_value;
                std::vector<hcd::RationalPoint> points;
                points.reserve(xs.size());
                for (std::size_t i = 0; i < xs.size(); ++i) points.emplace_back(binary64(xs[i]), binary64(ys[i]));
                std::vector<hcd::NearestDecision> decisions;
                decisions.reserve(points.size());
                hcd::QueryWorkspace workspace;
                {
                    nb::gil_scoped_release release;
                    for (const hcd::RationalPoint& point : points) {
                        decisions.push_back(domain.decide_nearest(point, &gap_squared, workspace));
                    }
                }
                nb::list out;
                last_exact_fallbacks = workspace.exact_fallbacks;
                for (const hcd::NearestDecision& decision : decisions) {
                    const int status = decision.beyond_gap ? BEYOND_GAP : decision.ambiguous ? AMBIGUOUS : PROJECTED;
                    out.append(nb::make_tuple(status, static_cast<int>(decision.kind), decision.piece));
                }
                return out;
            },
            "xs"_a, "ys"_a, "gap"_a)
        .def(
            "place_binary64",
            // Batch placement at decided contacts (binary64 probes injected
            // exactly; `pieces` and `vertices` from `project_binary64`): per
            // contact (status, placement, half, cx, cy, radius, limited,
            // normal_x, normal_y).
            [](const hcd::CentreDomain2& domain, const std::vector<double>& xs, const std::vector<double>& ys,
               const std::vector<std::size_t>& pieces, const std::vector<int>& vertices,
               const std::vector<double>& guide_radii) {
                const std::size_t count = xs.size();
                if (ys.size() != count || pieces.size() != count || vertices.size() != count || guide_radii.size() != count) {
                    throw hcd::InvalidCentreDomainInputError("placement batch columns differ in length");
                }
                for (const int vertex : vertices) {
                    if (vertex != 0 && vertex != 1) throw hcd::InvalidCentreDomainInputError("a vertex flag must be 0 or 1");
                }
                std::vector<hcd::PlacedCircle> placed;
                placed.reserve(count);
                {
                    nb::gil_scoped_release release;
                    for (std::size_t i = 0; i < count; ++i) {
                        placed.push_back(domain.place_on_normal(pieces[i], xs[i], ys[i], vertices[i] == 1, guide_radii[i]));
                    }
                }
                nb::list out;
                for (const hcd::PlacedCircle& circle : placed) {
                    out.append(nb::make_tuple(static_cast<int>(circle.status), static_cast<int>(circle.placement), circle.half,
                                              circle.cx, circle.cy, circle.radius, circle.limited, circle.normal_x,
                                              circle.normal_y));
                }
                return out;
            },
            "xs"_a, "ys"_a, "pieces"_a, "vertices"_a, "guide_radii"_a);
}
