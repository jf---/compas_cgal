#include "boundary_events.h"

#include <algorithm>
#include <cstdint>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <tuple>
#include <utility>
#include <variant>

#include <CGAL/number_utils.h>

namespace {

using IntersectionPoint = std::pair<GpsPoint, unsigned int>;
using IntersectionResult = std::variant<IntersectionPoint, GpsXCurve>;

struct VertexState {
    GpsPoint point;
    std::string exact_record;
    std::vector<std::string> incidents;
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

std::vector<std::string> normalized_support_coefficients(
    const GpsXCurve& curve)
{
    std::vector<Epeck::FT> coefficients;
    if (curve.is_linear()) {
        const Epeck::Line_2 line = curve.supporting_line();
        coefficients = {line.a(), line.b(), line.c()};
        const auto pivot = std::find_if(
            coefficients.begin(),
            coefficients.end(),
            [](const Epeck::FT& value) {
                return CGAL::sign(value) != CGAL::ZERO;
            });
        if (pivot == coefficients.end()) {
            throw std::runtime_error(
                "boundary line support has zero coefficients");
        }
        const Epeck::FT scale = *pivot;
        for (Epeck::FT& coefficient : coefficients) {
            coefficient /= scale;
        }
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
    std::vector<std::string> records;
    records.reserve(coefficients.size());
    for (const Epeck::FT& coefficient : coefficients) {
        records.push_back(exact_rational_record(coefficient));
    }
    return records;
}

std::string support_record(
    const std::string& kind,
    const std::vector<std::string>& coefficients)
{
    std::vector<std::string> fields{kind};
    fields.insert(
        fields.end(),
        coefficients.begin(),
        coefficients.end());
    return tagged_record("incident-support-id-v1", fields);
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
        throw std::runtime_error(
            "boundary endpoint is absent from vertex state");
    }
    return static_cast<std::size_t>(
        std::distance(vertices.begin(), found));
}

std::vector<VertexState> build_vertices(
    const std::vector<GpsXCurve>& curves,
    const std::vector<std::string>& support_ids)
{
    std::vector<VertexState> vertices;
    for (std::size_t index = 0; index < curves.size(); ++index) {
        for (const GpsPoint& point :
             {curves[index].source(), curves[index].target()}) {
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
            found->incidents.push_back(support_ids[index]);
        }
    }
    for (VertexState& vertex : vertices) {
        std::sort(
            vertex.incidents.begin(),
            vertex.incidents.end());
        vertex.incidents.erase(
            std::unique(
                vertex.incidents.begin(),
                vertex.incidents.end()),
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
        for (std::size_t previous = 0;
             previous < index;
             ++previous) {
            if (vertices[previous].incidents
                == vertices[index].incidents) {
                ++ordinal;
            }
        }
        std::vector<std::string> fields =
            vertices[index].incidents;
        fields.push_back(u64_record(ordinal));
        vertices[index].vertex_id = tagged_record(
            "boundary-vertex-id-v1",
            {
                tagged_record(
                    "support-intersection-v1",
                    fields),
            });
    }
    return vertices;
}

std::vector<BoundaryFeatureRecord2> records_from_curves(
    const std::vector<GpsXCurve>& curves)
{
    std::vector<std::vector<std::string>> coefficients;
    std::vector<std::string> support_ids;
    coefficients.reserve(curves.size());
    support_ids.reserve(curves.size());
    for (const GpsXCurve& curve : curves) {
        coefficients.push_back(
            normalized_support_coefficients(curve));
        support_ids.push_back(
            support_record(
                curve.is_linear() ? "line" : "circle",
                coefficients.back()));
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
        const std::string feature_id = tagged_record(
            "boundary-feature-id-v1",
            {
                support_ids[index],
                source.vertex_id,
                target.vertex_id,
                "material-left",
                u64_record(1),
            });
        records.push_back(
            {
                curve,
                curve.is_linear() ? "line" : "circle",
                coefficients[index],
                support_ids[index],
                source.exact_record,
                target.exact_record,
                source.vertex_id,
                target.vertex_id,
                "left",
                trim,
                feature_id,
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
        throw std::runtime_error(
            "intersection point is absent from exact pair results");
    }
    std::vector<std::string> incidents{
        first.support_id,
        second.support_id,
    };
    std::sort(incidents.begin(), incidents.end());
    return tagged_record(
        "boundary-vertex-id-v1",
        {
            tagged_record(
                "support-intersection-v1",
                {
                    incidents[0],
                    incidents[1],
                    u64_record(
                        static_cast<std::size_t>(
                            std::distance(points.begin(), found))),
                }),
        });
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
    for (const IntersectionResult& intersection :
         intersections) {
        if (std::holds_alternative<GpsXCurve>(intersection)) {
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
                    0,
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
        events.push_back(
            {
                endpoint
                    ? "vertex"
                    : hit.second > 1
                    ? "tangent"
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
                       left.first_feature_id,
                       left.second_feature_id)
                < std::tie(
                       right.kind,
                       right.vertex_id,
                       right.first_feature_id,
                       right.second_feature_id);
        });
    return events;
}

} // namespace

std::vector<BoundaryFeatureRecord2> extract_boundary_records(
    const Stock2& stock)
{
    return records_from_curves(stock_curves(stock));
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
        throw std::out_of_range(
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
                       left.first_feature_id,
                       left.second_feature_id)
                < std::tie(
                       right.kind,
                       right.vertex_id,
                       right.first_feature_id,
                       right.second_feature_id);
        });
    events.erase(
        std::unique(
            events.begin(),
            events.end(),
            [](const BoundaryEvent2& left,
               const BoundaryEvent2& right) {
                return left.kind == right.kind
                    && left.vertex_id == right.vertex_id;
            }),
        events.end());
    return events;
}
