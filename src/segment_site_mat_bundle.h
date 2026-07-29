#pragma once

#include "segment_site_catalog.h"
#include "segment_site_mat_numeric_table.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

class SegmentSiteMatBundle2 {
public:
  static SegmentSiteMatBundle2
  build(const CanonicalReachInput2 &input,
        const MatStationSpacingMm2 &station_spacing,
        const MatSagittaBoundMm2 &max_sagitta,
        std::size_t max_refinement_depth);

  const MatProposalSamplingGraph2 &sampled() const noexcept;
  const CanonicalMatSiteCatalog2 &catalog() const noexcept;
  const std::string &center_domain_digest() const noexcept;
  const MatNumericMatTable2 &numeric_table() const noexcept;
  std::vector<std::string> edge_ids() const;
  const std::vector<MatNeckEvidenceV1> &neck_evidence() const noexcept;
  const std::vector<std::string> &sample_parameter_ids() const noexcept;
  const std::vector<std::string> &neck_owner_ids() const noexcept;
  const std::vector<std::vector<std::string>> &
  neck_defining_site_ids() const noexcept;

  std::pair<std::vector<std::int64_t>, std::vector<std::string>>
  validate_and_classify_necks(
      const std::string &mat_certificate,
      const std::vector<std::string> &neck_evidence,
      const std::vector<std::string> &squared_width_boundaries) const;

  MatNumericMatTable2 release_numeric_table() && noexcept;

private:
  SegmentSiteMatBundle2(
      MatProposalSamplingGraph2 sampled, CanonicalMatSiteCatalog2 catalog,
      std::string center_domain_digest, MatNumericMatTable2 numeric_table,
      std::vector<MatNeckEvidenceV1> neck_evidence,
      std::vector<std::string> sample_parameter_ids,
      std::vector<std::string> neck_owner_ids,
      std::vector<std::vector<std::string>> neck_defining_site_ids);

  MatProposalSamplingGraph2 sampled_;
  CanonicalMatSiteCatalog2 catalog_;
  std::string center_domain_digest_;
  MatNumericMatTable2 numeric_table_;
  std::vector<MatNeckEvidenceV1> neck_evidence_;
  std::vector<std::string> sample_parameter_ids_;
  std::vector<std::string> neck_owner_ids_;
  std::vector<std::vector<std::string>> neck_defining_site_ids_;
};
