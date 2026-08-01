#include "segment_site_catalog_neck.h"

#include "segment_site_catalog.h"
#include "segment_site_catalog_graph.h"
#include "segment_site_clipping.h"
#include "segment_site_rational_sources.h"

#include <algorithm>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace {

const CanonicalMatRationalPointSource2 *
point_source(const CanonicalMatRationalSources2 &sources,
             const std::string &stable_id) {
  const auto found = std::lower_bound(
      sources.points().begin(), sources.points().end(), stable_id,
      [](const auto &point, const std::string &identity) {
        return point.stable_site_id < identity;
      });
  if (found == sources.points().end() || found->stable_site_id != stable_id) {
    return nullptr;
  }
  return &*found;
}

const CanonicalMatRationalOpenSegmentSource2 *
segment_source(const CanonicalMatRationalSources2 &sources,
               const std::string &stable_id) {
  const auto found = std::lower_bound(
      sources.segments().begin(), sources.segments().end(), stable_id,
      [](const auto &segment, const std::string &identity) {
        return segment.stable_site_id < identity;
      });
  if (found == sources.segments().end() || found->stable_site_id != stable_id) {
    return nullptr;
  }
  return &*found;
}

MatExactPointSiteSource2
exact_point_source(const CanonicalMatRationalPointSource2 &point) {
  return {
      point.stable_site_id,
      {point.x, 0},
      {point.y, 0},
      1,
  };
}

CORE::BigRat
squared_distance_to_support(const RationalPrimitiveParameterization2 &primitive,
                            const MatExactOpenSegmentSource2 &segment) {
  const CORE::BigRat line_constant =
      segment.line_a * primitive.x_coefficients[0] +
      segment.line_b * primitive.y_coefficients[0] + segment.line_c;
  const CORE::BigRat line_direction =
      segment.line_a * primitive.x_coefficients[1] +
      segment.line_b * primitive.y_coefficients[1];
  const CORE::BigRat line_norm =
      segment.line_a * segment.line_a + segment.line_b * segment.line_b;
  if (line_direction != 0 || line_norm == 0) {
    throw UnsupportedCanonicalMatClearanceProfileError(
        "parallel MAT profile is not constant on its source supports");
  }
  return line_constant * line_constant / line_norm;
}

MatSquaredClearanceFunction2
parallel_squared_clearance(const MatExactOpenSegmentSource2 &first,
                           const MatExactOpenSegmentSource2 &second) {
  const RationalPrimitiveParameterization2 primitive =
      parallel_segment_bisector_parameterization(first, second);
  const CORE::BigRat first_distance =
      squared_distance_to_support(primitive, first);
  const CORE::BigRat second_distance =
      squared_distance_to_support(primitive, second);
  if (first_distance != second_distance) {
    throw UnsupportedCanonicalMatClearanceProfileError(
        "parallel MAT profile is not equidistant from its sources");
  }
  const ClearanceRootBoundary2 zero = parallel_segment_clearance_boundary(
      primitive, first, second, first_distance);
  const ClearanceRootBoundary2 negative = parallel_segment_clearance_boundary(
      primitive, first, second, first_distance + 1);
  if (zero.constant_sign != CGAL::ZERO ||
      negative.constant_sign != CGAL::NEGATIVE) {
    throw UnsupportedCanonicalMatClearanceProfileError(
        "parallel MAT profile disagrees with clearance clipping");
  }
  return MatSquaredClearanceFunction2::build({first_distance});
}

MatSquaredClearanceFunction2
nonparallel_squared_clearance(const MatExactOpenSegmentSource2 &first,
                              const MatExactOpenSegmentSource2 &second) {
  const CORE::BigRat determinant =
      first.line_a * second.line_b - second.line_a * first.line_b;
  const CORE::BigRat second_norm =
      second.line_a * second.line_a + second.line_b * second.line_b;
  const CORE::BigRat quadratic = second_norm * determinant * determinant;
  if (determinant == 0 || quadratic <= 0) {
    throw UnsupportedCanonicalMatClearanceProfileError(
        "nonparallel MAT profile has no positive quadratic clearance");
  }
  return MatSquaredClearanceFunction2::build({0, 0, quadratic});
}

MatSquaredClearanceFunction2 point_segment_squared_clearance(
    const CanonicalMatRationalPointSource2 &point,
    const CanonicalMatRationalOpenSegmentSource2 &segment) {
  const MatExactPointSiteSource2 exact_point = exact_point_source(point);
  const SourceParabolaParameterization2 primitive =
      source_parameterization(exact_point, segment.support);
  const auto has_radical = [](const std::vector<CORE::BigRat> &coefficients) {
    return std::any_of(
        coefficients.begin(), coefficients.end(),
        [](const CORE::BigRat &coefficient) { return coefficient != 0; });
  };
  if (has_radical(primitive.x_radical) || has_radical(primitive.y_radical) ||
      primitive.x_rational.empty() || primitive.y_rational.empty()) {
    throw UnsupportedCanonicalMatClearanceProfileError(
        "canonical P-S profile is not rational");
  }
  RationalPolynomial dx = primitive.x_rational;
  RationalPolynomial dy = primitive.y_rational;
  dx.front() -= point.x;
  dy.front() -= point.y;
  trim(dx);
  trim(dy);
  RationalPolynomial clearance = multiply(dx, dx);
  add_in_place(clearance, multiply(dy, dy));
  const ClearanceRootBoundary2 boundary =
      source_parabola_clearance_boundary(exact_point, segment.support, 0);
  if (boundary.constant_sign.has_value() ||
      boundary.primitive_coefficients !=
          primitive_integer_coefficients(clearance)) {
    throw UnsupportedCanonicalMatClearanceProfileError(
        "canonical P-S profile disagrees with clearance clipping");
  }
  return MatSquaredClearanceFunction2::build(std::move(clearance));
}

MatSquaredClearanceFunction2
squared_clearance(const MatExactGraphEdge2 &edge,
                  const CanonicalMatRationalSources2 &sources) {
  const auto *first_point = point_source(sources, edge.generator_site_ids[0]);
  const auto *second_point = point_source(sources, edge.generator_site_ids[1]);
  const auto *first_segment =
      segment_source(sources, edge.generator_site_ids[0]);
  const auto *second_segment =
      segment_source(sources, edge.generator_site_ids[1]);
  if (first_segment != nullptr && second_segment != nullptr &&
      first_point == nullptr && second_point == nullptr &&
      edge.primitive_kind == "LINE") {
    const CORE::BigRat determinant =
        first_segment->support.line_a * second_segment->support.line_b -
        second_segment->support.line_a * first_segment->support.line_b;
    return determinant == 0
               ? parallel_squared_clearance(first_segment->support,
                                            second_segment->support)
               : nonparallel_squared_clearance(first_segment->support,
                                               second_segment->support);
  }
  if (edge.primitive_kind == "PARABOLA" &&
      ((first_point != nullptr && second_segment != nullptr &&
        first_segment == nullptr && second_point == nullptr) ||
       (second_point != nullptr && first_segment != nullptr &&
        second_segment == nullptr && first_point == nullptr))) {
    return first_point != nullptr
               ? point_segment_squared_clearance(*first_point, *second_segment)
               : point_segment_squared_clearance(*second_point, *first_segment);
  }
  throw UnsupportedCanonicalMatClearanceProfileError(
      "canonical MAT edge has an unsupported clearance-source pair");
}

bool endpoints_equal(const MatParameterEndpoint2 &lhs,
                     const MatParameterEndpoint2 &rhs) {
  return lhs.parameter.has_value() == rhs.parameter.has_value() &&
         (!lhs.parameter.has_value() ||
          mat_endpoint_root_identity_v1(lhs) ==
              mat_endpoint_root_identity_v1(rhs)) &&
         lhs.provenance_ids == rhs.provenance_ids &&
         lhs.exact_evidence == rhs.exact_evidence;
}

} // namespace

MatClearanceProfileGraph2::MatClearanceProfileGraph2(
    MatExactGraph2 graph, std::vector<MatClearanceEdgeProfile2> profiles)
    : graph_(std::move(graph)), profiles_(std::move(profiles)) {}

MatClearanceProfileGraph2 MatClearanceProfileGraph2::build(
    MatExactGraph2 graph, std::vector<MatClearanceEdgeProfile2> profiles) {
  if (graph.edges.size() != profiles.size()) {
    throw IncompleteMatClearanceProfileGraphError(
        "MAT graph does not have one clearance profile per edge");
  }
  for (std::size_t index = 0; index < graph.edges.size(); ++index) {
    const MatExactGraphEdge2 &edge = graph.edges[index];
    const MatClearanceEdgeProfile2 &profile = profiles[index];
    if ((index > 0 && edge.edge_id <= graph.edges[index - 1].edge_id) ||
        profile.edge_id() != edge.edge_id ||
        profile.defining_site_ids() != edge.generator_site_ids ||
        !endpoints_equal(profile.lower(), edge.source_endpoint) ||
        !endpoints_equal(profile.upper(), edge.target_endpoint)) {
      throw IncompleteMatClearanceProfileGraphError(
          "MAT clearance profile is not bound to its canonical edge");
    }
  }
  return MatClearanceProfileGraph2(std::move(graph), std::move(profiles));
}

const MatExactGraph2 &MatClearanceProfileGraph2::graph() const noexcept {
  return graph_;
}

const std::vector<MatClearanceEdgeProfile2> &
MatClearanceProfileGraph2::profiles() const noexcept {
  return profiles_;
}

MatClearanceProfileGraph2
canonical_l_shape_mat_clearance_graph(const CanonicalReachInput2 &input,
                                      const CORE::BigRat &radius_squared) {
  MatExactGraph2 graph = canonical_l_shape_mat_graph(input, radius_squared);
  const CanonicalMatRationalSources2 sources =
      CanonicalMatRationalSources2::build(canonical_mat_site_catalog(input));
  std::vector<MatClearanceEdgeProfile2> profiles;
  profiles.reserve(graph.edges.size());
  for (const MatExactGraphEdge2 &edge : graph.edges) {
    profiles.push_back(MatClearanceEdgeProfile2::build(
        edge, squared_clearance(edge, sources)));
  }
  return MatClearanceProfileGraph2::build(std::move(graph),
                                          std::move(profiles));
}
