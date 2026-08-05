#include "boundary_events.h"
#include "result.h"

#include <algorithm>
#include <cstdint>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <tuple>
#include <utility>
#include <variant>

#include <boost/multiprecision/cpp_int.hpp>
#include <CGAL/number_utils.h>

namespace {

using IntersectionPoint = std::pair<GpsPoint, unsigned int>;
using IntersectionResult = std::variant<IntersectionPoint, GpsXCurve>;
using BigInt = boost::multiprecision::cpp_int;

struct VertexState {
    GpsPoint point;
    std::string exact_record;
    std::vector<std::pair<std::string, std::string>> incidents;
    std::string vertex_id;
};

std::string u64_record(std::size_t value)
{
    std::string result;
    for (int shift = 56; shift >= 0; shift -= 8) {
        result.push_back(
            static_cast<char>(
                (static_cast<std::uint64_t>(value) >> shift)
                & 0xffU));
    }
    return result;
}

std::string tagged_record(
    std::string_view tag,
    const std::vector<std::string>& fields)
{
    std::string result(tag);
    result.push_back('\0');
    result += u64_record(fields.size());
    for (const std::string& field : fields) {
        result += u64_record(field.size());
        result += field;
    }
    return result;
}

std::string ccan_node(char kind, const std::string& payload)
{
    std::string result("CCAN\0\1", 6);
    result.push_back(kind);
    result += u64_record(payload.size());
    result += payload;
    return result;
}

std::string ccan_bytes(const std::string& value)
{
    return ccan_node('B', value);
}

std::string ccan_integer(const BigInt& value)
{
    BigInt magnitude = value < 0 ? -value : value;
    std::vector<unsigned char> bytes;
    if (magnitude != 0) {
        export_bits(
            magnitude,
            std::back_inserter(bytes),
            8,
            true);
    }
    std::string payload(
        1,
        value < 0
            ? static_cast<char>(1)
            : static_cast<char>(0));
    payload.append(
        reinterpret_cast<const char*>(bytes.data()),
        bytes.size());
    return ccan_node('I', payload);
}

std::string ccan_sequence(
    const std::vector<std::string>& values)
{
    std::string payload = u64_record(values.size());
    for (const std::string& value : values) {
        payload += ccan_bytes(value);
    }
    return ccan_node('S', payload);
}

std::string ccan_tagged(
    const std::string& tag,
    const std::string& payload)
{
    return ccan_node(
        'T',
        ccan_bytes(tag) + ccan_bytes(payload));
}

std::string ccan_map(
    std::vector<std::pair<std::string, std::string>> fields)
{
    std::vector<std::pair<std::string, std::string>> encoded;
    encoded.reserve(fields.size());
    for (const auto& [key, value] : fields) {
        encoded.emplace_back(
            ccan_bytes(key),
            ccan_bytes(value));
    }
    std::sort(encoded.begin(), encoded.end());
    std::string payload = u64_record(encoded.size());
    for (const auto& [key, value] : encoded) {
        payload += key;
        payload += value;
    }
    return ccan_node('M', payload);
}

BigInt greatest_common_divisor(BigInt left, BigInt right)
{
    left = left < 0 ? -left : left;
    right = right < 0 ? -right : right;
    while (right != 0) {
        const BigInt remainder = left % right;
        left = right;
        right = remainder;
    }
    return left;
}

std::pair<BigInt, BigInt> exact_fraction(
    const Epeck::FT& value)
{
    const auto& exact = CGAL::exact(value);
    return {
        boost::multiprecision::numerator(exact),
        boost::multiprecision::denominator(exact),
    };
}

std::string exact_rational_record(const Epeck::FT& value)
{
    std::ostringstream stream;
    stream << CGAL::exact(value);
    return tagged_record("exact-rational-v1", {stream.str()});
}

std::string exact_coordinate_record(const GpsPoint::CoordNT& value)
{
    return tagged_record(
        "one-root-coordinate-v1",
        {
            value.is_extended() ? "extended" : "rational",
            exact_rational_record(value.a0()),
            exact_rational_record(value.a1()),
            exact_rational_record(value.root()),
        });
}

std::string exact_point_record(const GpsPoint& point)
{
    return tagged_record(
        "exact-one-root-point-v1",
        {
            exact_coordinate_record(point.x()),
            exact_coordinate_record(point.y()),
        });
}

bool point_equal(const GpsPoint& left, const GpsPoint& right)
{
    return GpsTraits().compare_xy_2_object()(left, right)
        == CGAL::EQUAL;
}

struct SupportState {
    std::string kind;
    std::vector<std::string> exact_coefficients;
    std::vector<std::string> primitive_coefficients;
    std::string support_id;
};

SupportState support_state(const GpsXCurve& curve)
{
    std::vector<Epeck::FT> coefficients;
    const std::string kind =
        curve.is_linear() ? "line" : "circle";
    if (curve.is_linear()) {
        const Epeck::Line_2 line = curve.supporting_line();
        coefficients = {line.a(), line.b(), line.c()};
    } else {
        const Epeck::Circle_2 circle = curve.supporting_circle();
        coefficients = {
            Epeck::FT(1),
            Epeck::FT(-2) * circle.center().x(),
            Epeck::FT(-2) * circle.center().y(),
            circle.center().x() * circle.center().x()
                + circle.center().y() * circle.center().y()
                - circle.squared_radius(),
        };
    }
    std::vector<std::string> exact_records;
    exact_records.reserve(coefficients.size());
    std::vector<std::pair<BigInt, BigInt>> fractions;
    fractions.reserve(coefficients.size());
    BigInt denominator_lcm = 1;
    for (const Epeck::FT& coefficient : coefficients) {
        exact_records.push_back(
            exact_rational_record(coefficient));
        const auto fraction = exact_fraction(coefficient);
        fractions.push_back(fraction);
        denominator_lcm =
            denominator_lcm
            / greatest_common_divisor(
                denominator_lcm,
                fraction.second)
            * fraction.second;
    }
    std::vector<BigInt> primitive;
    primitive.reserve(fractions.size());
    BigInt common_factor = 0;
    for (const auto& [numerator, denominator] : fractions) {
        primitive.push_back(
            numerator * (denominator_lcm / denominator));
        common_factor = greatest_common_divisor(
            common_factor,
            primitive.back());
    }
    if (common_factor == 0) {
        throw DegenerateBoundarySupportError(
            "boundary support has zero coefficients");
    }
    for (BigInt& coefficient : primitive) {
        coefficient /= common_factor;
    }
    const auto first_nonzero = std::find_if(
        primitive.begin(),
        primitive.end(),
        [](const BigInt& value) {
            return value != 0;
        });
    if (*first_nonzero < 0) {
        for (BigInt& coefficient : primitive) {
            coefficient = -coefficient;
        }
    }
    std::vector<std::string> primitive_text;
    std::vector<std::string> encoded_integers;
    primitive_text.reserve(primitive.size());
    encoded_integers.reserve(primitive.size());
    for (const BigInt& coefficient : primitive) {
        primitive_text.push_back(
            coefficient.convert_to<std::string>());
        encoded_integers.push_back(
            ccan_integer(coefficient));
    }
    const std::string support_id = ccan_tagged(
        "incident-support-id-v1",
        ccan_map(
            {
                {
                    "coefficients",
                    ccan_sequence(encoded_integers),
                },
                {
                    "kind",
                    ccan_tagged(
                        curve.is_linear()
                            ? "line-support-v1"
                            : "circle-support-v1",
                        {}),
                },
            }));
    return {
        kind,
        std::move(exact_records),
        std::move(primitive_text),
        support_id,
    };
}

void append_polygon_curves(
    const GpsPolygon& polygon,
    std::vector<GpsXCurve>& curves)
{
    curves.insert(
        curves.end(),
        polygon.curves_begin(),
        polygon.curves_end());
}

std::vector<GpsXCurve> stock_curves(const Stock2& stock)
{
    std::vector<GpsPolygonWithHoles> components;
    stock.set().polygons_with_holes(
        std::back_inserter(components));
    std::vector<GpsXCurve> curves;
    for (const GpsPolygonWithHoles& component : components) {
        append_polygon_curves(
            component.outer_boundary(),
            curves);
        for (auto hole = component.holes_begin();
             hole != component.holes_end();
             ++hole) {
            append_polygon_curves(*hole, curves);
        }
    }
    return curves;
}

std::size_t find_vertex(
    const std::vector<VertexState>& vertices,
    const GpsPoint& point)
{
    const auto found = std::find_if(
        vertices.begin(),
        vertices.end(),
        [&point](const VertexState& vertex) {
            return point_equal(vertex.point, point);
        });
    if (found == vertices.end()) {
        throw MissingBoundaryEndpointError(
            "boundary endpoint is absent from vertex state");
    }
    return static_cast<std::size_t>(
        std::distance(vertices.begin(), found));
}

std::string incident_record(
    const std::string& support_id,
    std::string_view orientation_tag)
{
    return ccan_tagged(
        "incident-support-v1",
        ccan_map(
            {
                {
                    "orientation",
                    ccan_tagged(
                        std::string(orientation_tag),
                        {}),
                },
                {"support-id", support_id},
            }));
}

std::vector<VertexState> build_vertices(
    const std::vector<GpsXCurve>& curves,
    const std::vector<std::string>& support_ids)
{
    std::vector<VertexState> vertices;
    for (std::size_t index = 0; index < curves.size(); ++index) {
        for (const auto& [point, orientation] :
             {
                 std::pair{
                     curves[index].source(),
                     std::string_view("trim-leaving-v1")},
                 std::pair{
                     curves[index].target(),
                     std::string_view("trim-entering-v1")},
             }) {
            auto found = std::find_if(
                vertices.begin(),
                vertices.end(),
                [&point](const VertexState& vertex) {
                    return point_equal(vertex.point, point);
                });
            if (found == vertices.end()) {
                vertices.push_back(
                    {point, exact_point_record(point), {}, {}});
                found = std::prev(vertices.end());
            }
            found->incidents.emplace_back(
                support_ids[index],
                incident_record(
                    support_ids[index],
                    orientation));
        }
    }
    for (VertexState& vertex : vertices) {
        std::sort(
            vertex.incidents.begin(),
            vertex.incidents.end());
    }
    std::sort(
        vertices.begin(),
        vertices.end(),
        [](const VertexState& left, const VertexState& right) {
            return GpsTraits().compare_xy_2_object()(
                       left.point,
                       right.point)
                == CGAL::SMALLER;
        });
    for (std::size_t index = 0; index < vertices.size(); ++index) {
        std::size_t ordinal = 0;
        std::vector<std::string> support_ids_at_vertex;
        std::vector<std::string> incident_records;
        std::vector<std::string> all_incident_records;
        for (const auto& [support_id, incident] :
             vertices[index].incidents) {
            all_incident_records.push_back(incident);
            if (support_ids_at_vertex.empty()
                || support_ids_at_vertex.back() != support_id) {
                support_ids_at_vertex.push_back(support_id);
                incident_records.push_back(incident);
            }
        }
        for (std::size_t previous = 0;
             previous < index;
             ++previous) {
            std::vector<std::string> previous_supports;
            for (const auto& [support_id, incident] :
                 vertices[previous].incidents) {
                static_cast<void>(incident);
                if (previous_supports.empty()
                    || previous_supports.back() != support_id) {
                    previous_supports.push_back(support_id);
                }
            }
            if (previous_supports == support_ids_at_vertex) {
                ++ordinal;
            }
        }
        if (support_ids_at_vertex.size() >= 2) {
            const bool has_repeated_support =
                all_incident_records.size()
                != support_ids_at_vertex.size();
            vertices[index].vertex_id = ccan_tagged(
                "boundary-vertex-id-v1",
                ccan_tagged(
                    has_repeated_support
                        ? "multi-incidence-intersection-v1"
                        : "support-intersection-v1",
                    ccan_map(
                        {
                            {
                                "incident-supports",
                                ccan_sequence(
                                    has_repeated_support
                                        ? all_incident_records
                                        : incident_records),
                            },
                            {
                                "intersection-ordinal",
                                ccan_integer(ordinal),
                            },
                        })));
        } else {
            vertices[index].vertex_id = ccan_tagged(
                "boundary-vertex-id-v1",
                ccan_tagged(
                    "parameter-seam-v1",
                    ccan_map(
                        {
                            {
                                "incident-supports",
                                ccan_sequence(
                                    all_incident_records),
                            },
                            {
                                "seam-ordinal",
                                ccan_integer(ordinal),
                            },
                        })));
        }
    }
    return vertices;
}

std::vector<BoundaryFeatureRecord2> records_from_curves(
    const std::vector<GpsXCurve>& curves)
{
    std::vector<SupportState> supports;
    std::vector<std::string> support_ids;
    supports.reserve(curves.size());
    support_ids.reserve(curves.size());
    for (const GpsXCurve& curve : curves) {
        supports.push_back(support_state(curve));
        support_ids.push_back(
            supports.back().support_id);
    }
    const std::vector<VertexState> vertices =
        build_vertices(curves, support_ids);
    std::vector<BoundaryFeatureRecord2> records;
    records.reserve(curves.size());
    for (std::size_t index = 0; index < curves.size(); ++index) {
        const GpsXCurve& curve = curves[index];
        const VertexState& source =
            vertices[find_vertex(vertices, curve.source())];
        const VertexState& target =
            vertices[find_vertex(vertices, curve.target())];
        const std::string trim = tagged_record(
            "trimmed-domain-v1",
            {
                source.exact_record,
                target.exact_record,
                curve.is_directed_right()
                    ? "directed-right"
                    : "directed-left",
            });
        std::size_t overlap_multiplicity = 0;
        for (const GpsXCurve& candidate : curves) {
            if (curve.equals(candidate)) {
                ++overlap_multiplicity;
            }
        }
        const std::string feature_id = ccan_tagged(
            "boundary-feature-id-v1",
            ccan_map(
                {
                    {
                        "material-side",
                        ccan_tagged(
                            "material-side-inside-v1",
                            {}),
                    },
                    {
                        "overlap-multiplicity",
                        ccan_integer(
                            overlap_multiplicity),
                    },
                    {"source-vertex", source.vertex_id},
                    {"support-id", support_ids[index]},
                    {"target-vertex", target.vertex_id},
                }));
        records.push_back(
            {
                curve,
                supports[index].kind,
                supports[index].exact_coefficients,
                supports[index].primitive_coefficients,
                support_ids[index],
                source.exact_record,
                target.exact_record,
                source.vertex_id,
                target.vertex_id,
                incident_record(
                    support_ids[index],
                    "trim-leaving-v1"),
                incident_record(
                    support_ids[index],
                    "trim-entering-v1"),
                "left",
                trim,
                feature_id,
                overlap_multiplicity,
            });
    }
    std::sort(
        records.begin(),
        records.end(),
        [](const BoundaryFeatureRecord2& left,
           const BoundaryFeatureRecord2& right) {
            return left.feature_id < right.feature_id;
        });
    return records;
}

std::string point_vertex_id(
    const GpsPoint& point,
    const BoundaryFeatureRecord2& first,
    const BoundaryFeatureRecord2& second,
    const std::vector<IntersectionResult>& intersections)
{
    for (const BoundaryFeatureRecord2* feature :
         {&first, &second}) {
        if (point_equal(point, feature->curve.source())) {
            return feature->source_vertex_id;
        }
        if (point_equal(point, feature->curve.target())) {
            return feature->target_vertex_id;
        }
    }
    std::vector<GpsPoint> points;
    for (const IntersectionResult& intersection : intersections) {
        if (const auto* hit =
                std::get_if<IntersectionPoint>(&intersection)) {
            points.push_back(hit->first);
        }
    }
    std::sort(
        points.begin(),
        points.end(),
        [](const GpsPoint& left, const GpsPoint& right) {
            return GpsTraits().compare_xy_2_object()(left, right)
                == CGAL::SMALLER;
        });
    const auto found = std::find_if(
        points.begin(),
        points.end(),
        [&point](const GpsPoint& candidate) {
            return point_equal(candidate, point);
        });
    if (found == points.end()) {
        throw MissingBoundaryIntersectionError(
            "intersection point is absent from exact pair results");
    }
    std::vector<std::string> incidents{
        first.support_id,
        second.support_id,
    };
    std::sort(incidents.begin(), incidents.end());
    return ccan_tagged(
        "boundary-event-intersection-v1",
        ccan_map(
            {
                {
                    "intersection-ordinal",
                    ccan_integer(
                        static_cast<std::size_t>(
                            std::distance(points.begin(), found))),
                },
                {
                    "support-ids",
                    ccan_sequence(incidents),
                },
            }));
}

std::vector<BoundaryEvent2> classify_pair(
    const BoundaryFeatureRecord2& first,
    const BoundaryFeatureRecord2& second)
{
    std::vector<IntersectionResult> intersections;
    GpsTraits traits;
    traits.intersect_2_object()(
        first.curve,
        second.curve,
        std::back_inserter(intersections));
    std::vector<BoundaryEvent2> events;
    const unsigned int overlap_multiplicity =
        static_cast<unsigned int>(
            std::count_if(
                intersections.begin(),
                intersections.end(),
                [](const IntersectionResult& result) {
                    return std::holds_alternative<GpsXCurve>(
                        result);
                }));
    for (const IntersectionResult& intersection :
         intersections) {
        if (const auto* overlap =
                std::get_if<GpsXCurve>(&intersection)) {
            const std::string overlap_record = ccan_tagged(
                "exact-boundary-overlap-v1",
                ccan_map(
                    {
                        {
                            "first-support",
                            first.support_id,
                        },
                        {
                            "second-support",
                            second.support_id,
                        },
                        {
                            "source",
                            exact_point_record(
                                overlap->source()),
                        },
                        {
                            "target",
                            exact_point_record(
                                overlap->target()),
                        },
                    }));
            std::vector<std::string> overlap_endpoints{
                point_vertex_id(
                    overlap->source(),
                    first,
                    second,
                    intersections),
                point_vertex_id(
                    overlap->target(),
                    first,
                    second,
                    intersections),
            };
            std::sort(
                overlap_endpoints.begin(),
                overlap_endpoints.end());
            const std::string overlap_id = ccan_tagged(
                "boundary-overlap-id-v1",
                ccan_map(
                    {
                        {
                            "endpoint-vertices",
                            ccan_sequence(
                                overlap_endpoints),
                        },
                        {
                            "feature-ids",
                            ccan_sequence(
                                {
                                    std::min(
                                        first.feature_id,
                                        second.feature_id),
                                    std::max(
                                        first.feature_id,
                                        second.feature_id),
                                }),
                        },
                    }));
            events.push_back(
                {
                    "overlap",
                    std::min(
                        first.feature_id,
                        second.feature_id),
                    std::max(
                        first.feature_id,
                        second.feature_id),
                    {},
                    overlap_id,
                    overlap_record,
                    overlap_multiplicity,
                });
            continue;
        }
        const IntersectionPoint& hit =
            std::get<IntersectionPoint>(intersection);
        const bool endpoint =
            point_equal(hit.first, first.curve.source())
            || point_equal(hit.first, first.curve.target())
            || point_equal(hit.first, second.curve.source())
            || point_equal(hit.first, second.curve.target());
        const bool support_seam =
            endpoint
            && first.support_id == second.support_id;
        events.push_back(
            {
                support_seam
                    ? "seam"
                    : hit.second > 1
                    ? "tangent"
                    : endpoint
                    ? "vertex"
                    : "transverse",
                std::min(
                    first.feature_id,
                    second.feature_id),
                std::max(
                    first.feature_id,
                    second.feature_id),
                point_vertex_id(
                    hit.first,
                    first,
                    second,
                    intersections),
                {},
                {},
                hit.second,
            });
    }
    std::sort(
        events.begin(),
        events.end(),
        [](const BoundaryEvent2& left,
           const BoundaryEvent2& right) {
            return std::tie(
                       left.kind,
                       left.vertex_id,
                       left.overlap_id,
                       left.first_feature_id,
                       left.second_feature_id,
                       left.multiplicity)
                < std::tie(
                       right.kind,
                       right.vertex_id,
                       right.overlap_id,
                       right.first_feature_id,
                       right.second_feature_id,
                       right.multiplicity);
        });
    return events;
}

} // namespace

std::vector<BoundaryFeatureRecord2> extract_boundary_records(
    const Stock2& stock)
{
    return records_from_curves(stock_curves(stock));
}

std::string canonical_stock_boundary_identity(
    const Stock2& stock)
{
    const std::vector<BoundaryFeatureRecord2> records =
        extract_boundary_records(stock);
    std::vector<std::string> feature_ids;
    feature_ids.reserve(records.size());
    for (const BoundaryFeatureRecord2& record : records) {
        feature_ids.push_back(record.feature_id);
    }
    return ccan_tagged(
        "exact-stock-boundary-v1",
        ccan_sequence(feature_ids));
}

std::vector<BoundaryEvent2> classify_boundary_pair(
    const Stock2& stock,
    std::size_t first_index,
    std::size_t second_index)
{
    const std::vector<BoundaryFeatureRecord2> records =
        extract_boundary_records(stock);
    if (first_index >= records.size()
        || second_index >= records.size()) {
        throw BoundaryFeatureIndexError(
            "boundary feature index is out of range");
    }
    return classify_pair(
        records[first_index],
        records[second_index]);
}

std::vector<BoundaryEvent2> extract_boundary_events(
    const Stock2& stock)
{
    const std::vector<BoundaryFeatureRecord2> records =
        extract_boundary_records(stock);
    std::vector<BoundaryEvent2> events;
    for (std::size_t first = 0;
         first < records.size();
         ++first) {
        for (std::size_t second = first + 1;
             second < records.size();
             ++second) {
            std::vector<BoundaryEvent2> pair =
                classify_pair(
                    records[first],
                    records[second]);
            events.insert(
                events.end(),
                std::make_move_iterator(pair.begin()),
                std::make_move_iterator(pair.end()));
        }
    }
    std::sort(
        events.begin(),
        events.end(),
        [](const BoundaryEvent2& left,
           const BoundaryEvent2& right) {
            return std::tie(
                       left.kind,
                       left.vertex_id,
                       left.overlap_id,
                       left.first_feature_id,
                       left.second_feature_id,
                       left.multiplicity)
                < std::tie(
                       right.kind,
                       right.vertex_id,
                       right.overlap_id,
                       right.first_feature_id,
                       right.second_feature_id,
                       right.multiplicity);
        });
    events.erase(
        std::unique(
            events.begin(),
            events.end(),
            [](const BoundaryEvent2& left,
               const BoundaryEvent2& right) {
                return left.kind == right.kind
                    && left.vertex_id == right.vertex_id
                    && left.overlap_id
                        == right.overlap_id
                    && left.first_feature_id
                        == right.first_feature_id
                    && left.second_feature_id
                        == right.second_feature_id
                    && left.multiplicity
                        == right.multiplicity;
            }),
        events.end());
    return events;
}
