#include "segment_site_catalog_sampling.h"

#include "segment_site_catalog.h"
#include "segment_site_parameterization.h"
#include "segment_site_provenance.h"
#include "segment_site_rational_sources.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
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
  return found != sources.points().end() && found->stable_site_id == stable_id
             ? &*found
             : nullptr;
}

const CanonicalMatRationalOpenSegmentSource2 *
segment_source(const CanonicalMatRationalSources2 &sources,
               const std::string &stable_id) {
  const auto found = std::lower_bound(
      sources.segments().begin(), sources.segments().end(), stable_id,
      [](const auto &segment, const std::string &identity) {
        return segment.stable_site_id < identity;
      });
  return found != sources.segments().end() && found->stable_site_id == stable_id
             ? &*found
             : nullptr;
}

bool is_zero(const RationalPolynomial &coefficients) {
  return std::all_of(
      coefficients.begin(), coefficients.end(),
      [](const CORE::BigRat &coefficient) { return coefficient == 0; });
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

int nonparallel_branch_sign(const MatExactGraphEdge2 &edge) {
  const std::string negative = stable_dual_identity_v1(
      "segment-segment/branch-negative", edge.generator_site_ids);
  const std::string positive = stable_dual_identity_v1(
      "segment-segment/branch-positive", edge.generator_site_ids);
  if (edge.original_dual_id == negative) {
    return -1;
  }
  if (edge.original_dual_id == positive) {
    return 1;
  }
  throw UnsupportedCanonicalMatSamplingCurveError(
      "canonical nonparallel sampling edge has no exact branch identity");
}

MatTraits::Segment_2 branch_segment(const MatExactOpenSegmentSource2 &first,
                                    const MatExactOpenSegmentSource2 &second,
                                    const int branch_sign) {
  const MatExactOpenSegmentSource2 *ordered_first = &first;
  const MatExactOpenSegmentSource2 *ordered_second = &second;
  if (ordered_second->stable_site_id < ordered_first->stable_site_id) {
    std::swap(ordered_first, ordered_second);
  }
  const CORE::BigRat first_norm =
      ordered_first->line_a * ordered_first->line_a +
      ordered_first->line_b * ordered_first->line_b;
  const CORE::BigRat second_norm =
      ordered_second->line_a * ordered_second->line_a +
      ordered_second->line_b * ordered_second->line_b;
  const CORE::Expr radical = CORE::sqrt(CORE::Expr(first_norm * second_norm));
  const MatTraits::Line_2 branch(
      CORE::Expr(second_norm * ordered_first->line_a) +
          CORE::Expr(branch_sign * ordered_second->line_a) * radical,
      CORE::Expr(second_norm * ordered_first->line_b) +
          CORE::Expr(branch_sign * ordered_second->line_b) * radical,
      CORE::Expr(second_norm * ordered_first->line_c) +
          CORE::Expr(branch_sign * ordered_second->line_c) * radical);
  return MatTraits::Segment_2(branch.point(0), branch.point(1));
}

MatRationalSamplingCurve2
segment_segment_curve(const MatExactGraphEdge2 &edge,
                      const MatClearanceEdgeProfile2 &profile,
                      const CanonicalMatRationalOpenSegmentSource2 &first,
                      const CanonicalMatRationalOpenSegmentSource2 &second) {
  const CORE::BigRat determinant =
      first.support.line_a * second.support.line_b -
      first.support.line_b * second.support.line_a;
  if (determinant == 0) {
    RationalPrimitiveParameterization2 primitive =
        parallel_segment_bisector_parameterization(first.support,
                                                   second.support);
    return MatRationalSamplingCurve2::build(
        edge, profile, std::move(primitive.x_coefficients),
        std::move(primitive.y_coefficients));
  }

  const NonparallelSegmentBisectorParameterization2 primitive =
      nonparallel_segment_bisector_parameterization(
          first.support, second.support,
          branch_segment(first.support, second.support,
                         nonparallel_branch_sign(edge)));
  if (!is_zero(primitive.x_radical) || !is_zero(primitive.y_radical)) {
    throw UnsupportedCanonicalMatSamplingCurveError(
        "canonical nonparallel sampling chart is not rational");
  }
  return MatRationalSamplingCurve2::build(edge, profile, primitive.x_rational,
                                          primitive.y_rational);
}

MatRationalSamplingCurve2
point_segment_curve(const MatExactGraphEdge2 &edge,
                    const MatClearanceEdgeProfile2 &profile,
                    const CanonicalMatRationalPointSource2 &point,
                    const CanonicalMatRationalOpenSegmentSource2 &segment) {
  const SourceParabolaParameterization2 primitive =
      source_parameterization(exact_point_source(point), segment.support);
  if (!is_zero(primitive.x_radical) || !is_zero(primitive.y_radical)) {
    throw UnsupportedCanonicalMatSamplingCurveError(
        "canonical point-segment sampling chart is not rational");
  }
  return MatRationalSamplingCurve2::build(edge, profile, primitive.x_rational,
                                          primitive.y_rational);
}

MatRationalSamplingCurve2
sampling_curve(const MatExactGraphEdge2 &edge,
               const MatClearanceEdgeProfile2 &profile,
               const CanonicalMatRationalSources2 &sources) {
  if (edge.generator_site_ids.size() != 2) {
    throw UnsupportedCanonicalMatSamplingCurveError(
        "canonical sampling edge does not have two generators");
  }
  const auto *first_point = point_source(sources, edge.generator_site_ids[0]);
  const auto *second_point = point_source(sources, edge.generator_site_ids[1]);
  const auto *first_segment =
      segment_source(sources, edge.generator_site_ids[0]);
  const auto *second_segment =
      segment_source(sources, edge.generator_site_ids[1]);
  if (edge.primitive_kind == "LINE" && first_segment != nullptr &&
      second_segment != nullptr && first_point == nullptr &&
      second_point == nullptr) {
    return segment_segment_curve(edge, profile, *first_segment,
                                 *second_segment);
  }
  if (edge.primitive_kind == "PARABOLA" && first_point != nullptr &&
      second_segment != nullptr && first_segment == nullptr &&
      second_point == nullptr) {
    return point_segment_curve(edge, profile, *first_point, *second_segment);
  }
  if (edge.primitive_kind == "PARABOLA" && second_point != nullptr &&
      first_segment != nullptr && second_segment == nullptr &&
      first_point == nullptr) {
    return point_segment_curve(edge, profile, *second_point, *first_segment);
  }
  throw UnsupportedCanonicalMatSamplingCurveError(
      "canonical sampling edge has an unsupported source pair");
}

} // namespace

MatProposalSamplingGraph2::MatProposalSamplingGraph2(
    MatClearanceProfileGraph2 profile_graph,
    std::vector<std::int64_t> sample_offsets,
    std::vector<MatWorldXYProposalSample2> samples)
    : profile_graph_(std::move(profile_graph)),
      sample_offsets_(std::move(sample_offsets)), samples_(std::move(samples)) {
}

MatProposalSamplingGraph2
MatProposalSamplingGraph2::build(MatClearanceProfileGraph2 profile_graph,
                                 std::vector<MatProposalSamplingRun2> runs) {
  const auto &edges = profile_graph.graph().edges;
  if (runs.size() != edges.size()) {
    throw IncompleteMatProposalSamplingGraphError(
        "MAT graph does not have one proposal sampling run per edge");
  }

  std::vector<std::int64_t> offsets{
      0,
  };
  std::vector<MatWorldXYProposalSample2> samples;
  for (std::size_t index = 0; index < runs.size(); ++index) {
    if (runs[index].edge_id() != edges[index].edge_id) {
      throw IncompleteMatProposalSamplingGraphError(
          "MAT proposal sampling run is not in canonical edge order");
    }
    if (runs[index].samples().size() >
        static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max() -
                                 offsets.back())) {
      throw MatSamplingOffsetOverflowError(
          "MAT proposal sample offset exceeds int64");
    }
    offsets.push_back(offsets.back() +
                      static_cast<std::int64_t>(runs[index].samples().size()));
    samples.insert(samples.end(), runs[index].samples().begin(),
                   runs[index].samples().end());
  }
  return MatProposalSamplingGraph2(std::move(profile_graph), std::move(offsets),
                                   std::move(samples));
}

const MatClearanceProfileGraph2 &
MatProposalSamplingGraph2::profile_graph() const noexcept {
  return profile_graph_;
}

const std::vector<std::int64_t> &
MatProposalSamplingGraph2::sample_offsets() const noexcept {
  return sample_offsets_;
}

const std::vector<MatWorldXYProposalSample2> &
MatProposalSamplingGraph2::samples() const noexcept {
  return samples_;
}

MatProposalSamplingGraph2 canonical_l_shape_mat_proposal_graph(
    const CanonicalReachInput2 &input, const CORE::BigRat &radius_squared,
    const MatStationSpacingMm2 &station_spacing,
    const MatSagittaBoundMm2 &max_sagitta,
    const std::size_t max_refinement_depth) {
  MatClearanceProfileGraph2 profile_graph =
      canonical_l_shape_mat_clearance_graph(input, radius_squared);
  const CanonicalMatRationalSources2 sources =
      CanonicalMatRationalSources2::build(canonical_mat_site_catalog(input));
  std::vector<MatProposalSamplingRun2> runs;
  runs.reserve(profile_graph.graph().edges.size());
  for (std::size_t index = 0; index < profile_graph.graph().edges.size();
       ++index) {
    const MatRationalSamplingCurve2 curve =
        sampling_curve(profile_graph.graph().edges[index],
                       profile_graph.profiles()[index], sources);
    runs.push_back(proposal_sampling_run(curve, station_spacing, max_sagitta,
                                         max_refinement_depth));
  }
  return MatProposalSamplingGraph2::build(std::move(profile_graph),
                                          std::move(runs));
}
