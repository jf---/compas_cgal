#include "segment_site_mat_bundle.h"

#include "reachable_certificate_encoding_2.h"
#include "reachable_domain_2.h"
#include "segment_site_mat_certificate.h"
#include "segment_site_neck_classification.h"

#include <utility>

SegmentSiteMatBundle2::SegmentSiteMatBundle2(
    MatProposalSamplingGraph2 sampled, CanonicalMatSiteCatalog2 catalog,
    std::string center_domain_digest, MatNumericMatTable2 numeric_table,
    std::vector<std::string> sample_parameter_ids,
    std::vector<std::string> neck_owner_ids,
    std::vector<std::vector<std::string>> neck_defining_site_ids)
    : sampled_(std::move(sampled)), catalog_(std::move(catalog)),
      center_domain_digest_(std::move(center_domain_digest)),
      numeric_table_(std::move(numeric_table)),
      sample_parameter_ids_(std::move(sample_parameter_ids)),
      neck_owner_ids_(std::move(neck_owner_ids)),
      neck_defining_site_ids_(std::move(neck_defining_site_ids)) {}

SegmentSiteMatBundle2
SegmentSiteMatBundle2::build(const CanonicalReachInput2 &input,
                             const MatStationSpacingMm2 &station_spacing,
                             const MatSagittaBoundMm2 &max_sagitta,
                             const std::size_t max_refinement_depth) {
  const ReachableDomain2 reachable = ReachableDomain2::build(input);
  const ReachableDomainCertificateIdentity2 domain_identity =
      reachable_domain_certificate_identity(input, reachable);
  const MatToolRadiusMm2 tool_radius =
      MatToolRadiusMm2::build(input.binary64_radius);
  MatProposalSamplingGraph2 sampled = canonical_l_shape_mat_proposal_graph(
      input, tool_radius.squared_exact(), station_spacing, max_sagitta,
      max_refinement_depth);
  CanonicalMatSiteCatalog2 catalog = canonical_mat_site_catalog(input);
  MatCertifiedExactProjection2 exact = certified_mat_exact_projection_v1(
      sampled.profile_graph(), catalog, domain_identity.center_domain_digest());

  std::vector<std::string> sample_parameter_ids;
  sample_parameter_ids.reserve(sampled.samples().size());
  for (const MatWorldXYProposalSample2 &sample : sampled.samples()) {
    sample_parameter_ids.push_back(sample.exact_parameter_id());
  }

  std::vector<std::string> neck_owner_ids;
  neck_owner_ids.reserve(exact.neck_evidence.size());
  std::vector<std::vector<std::string>> neck_defining_site_ids;
  neck_defining_site_ids.reserve(exact.neck_evidence.size());
  for (const MatNeckEvidenceV1 &evidence : exact.neck_evidence) {
    neck_owner_ids.push_back(evidence.evidence().owner_id());
    neck_defining_site_ids.push_back(evidence.evidence().defining_site_ids());
  }

  const std::string center_domain_digest =
      domain_identity.center_domain_digest();
  MatNumericProposalTable2 proposal{
      std::move(exact.graph),
      numeric_sample_table(sampled, tool_radius),
  };
  MatNumericMatTable2 numeric_table{
      numeric_node_reporting_coordinates(sampled),
      std::move(proposal),
      std::move(exact.necks.neck_evidence),
      std::move(exact.necks.neck_cut_offsets),
      std::move(exact.necks.neck_cut_edge_ids),
      center_domain_digest,
      std::move(exact.certificate).release_canonical_bytes(),
  };
  return SegmentSiteMatBundle2(
      std::move(sampled), std::move(catalog), center_domain_digest,
      std::move(numeric_table), std::move(sample_parameter_ids),
      std::move(neck_owner_ids), std::move(neck_defining_site_ids));
}

const MatProposalSamplingGraph2 &
SegmentSiteMatBundle2::sampled() const noexcept {
  return sampled_;
}

const CanonicalMatSiteCatalog2 &
SegmentSiteMatBundle2::catalog() const noexcept {
  return catalog_;
}

const std::string &
SegmentSiteMatBundle2::center_domain_digest() const noexcept {
  return center_domain_digest_;
}

const MatNumericMatTable2 &
SegmentSiteMatBundle2::numeric_table() const noexcept {
  return numeric_table_;
}

std::vector<std::string> SegmentSiteMatBundle2::edge_ids() const {
  std::vector<std::string> result;
  result.reserve(sampled_.profile_graph().graph().edges.size());
  for (const MatExactGraphEdge2 &edge :
       sampled_.profile_graph().graph().edges) {
    result.push_back(edge.edge_id);
  }
  return result;
}

const std::vector<std::string> &
SegmentSiteMatBundle2::sample_parameter_ids() const noexcept {
  return sample_parameter_ids_;
}

const std::vector<std::string> &
SegmentSiteMatBundle2::neck_owner_ids() const noexcept {
  return neck_owner_ids_;
}

const std::vector<std::vector<std::string>> &
SegmentSiteMatBundle2::neck_defining_site_ids() const noexcept {
  return neck_defining_site_ids_;
}

std::pair<std::vector<std::int64_t>, std::vector<std::string>>
SegmentSiteMatBundle2::validate_and_classify_necks(
    const std::string &mat_certificate,
    const std::vector<std::string> &neck_evidence,
    const std::vector<std::string> &squared_width_boundaries) const {
  static_cast<void>(replay_mat_certificate_v1(sampled_.profile_graph(),
                                              catalog_, center_domain_digest_,
                                              mat_certificate));
  return validate_and_classify_necks_v1(sampled_.profile_graph(), neck_evidence,
                                        squared_width_boundaries);
}

MatNumericMatTable2 SegmentSiteMatBundle2::release_numeric_table() && noexcept {
  return std::move(numeric_table_);
}

MatNumericMatTable2
canonical_l_shape_mat_numeric_table(const CanonicalReachInput2 &input,
                                    const MatStationSpacingMm2 &station_spacing,
                                    const MatSagittaBoundMm2 &max_sagitta,
                                    const std::size_t max_refinement_depth) {
  return SegmentSiteMatBundle2::build(input, station_spacing, max_sagitta,
                                      max_refinement_depth)
      .release_numeric_table();
}
