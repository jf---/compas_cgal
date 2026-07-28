#include "segment_site_mat_certificate.h"
#include "segment_site_mat_numeric_table.h"

#include "segment_site_catalog.h"
#include "segment_site_catalog_neck.h"

#include <array>
#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

namespace {

std::string bytes_from_hex(const std::string_view hex) {
  const auto nibble = [](const char digit) {
    if (digit >= '0' && digit <= '9') {
      return digit - '0';
    }
    return digit - 'a' + 10;
  };
  std::string bytes;
  bytes.reserve(hex.size() / 2);
  for (std::size_t index = 0; index < hex.size(); index += 2) {
    bytes.push_back(
        static_cast<char>((nibble(hex[index]) << 4) | nibble(hex[index + 1])));
  }
  return bytes;
}

compas::RowMatrixXd matrix(const std::vector<std::array<double, 2>> &points) {
  compas::RowMatrixXd result(static_cast<Eigen::Index>(points.size()), 2);
  for (std::size_t index = 0; index < points.size(); ++index) {
    result(static_cast<Eigen::Index>(index), 0) = points[index][0];
    result(static_cast<Eigen::Index>(index), 1) = points[index][1];
  }
  return result;
}

CanonicalReachInput2 l_shape_input(const bool reversed = false) {
  return canonical_reach_input(reversed ? matrix({
                                              {0.0, 6.0},
                                              {2.0, 6.0},
                                              {2.0, 2.0},
                                              {6.0, 2.0},
                                              {6.0, 0.0},
                                              {0.0, 0.0},
                                          })
                                        : matrix({
                                              {0.0, 0.0},
                                              {6.0, 0.0},
                                              {6.0, 2.0},
                                              {2.0, 2.0},
                                              {2.0, 6.0},
                                              {0.0, 6.0},
                                          }),
                               {}, 0.5);
}

MatNumericMatTable2 table(const CanonicalReachInput2 &input,
                          const double station_spacing) {
  return canonical_l_shape_mat_numeric_table(
      input, MatStationSpacingMm2::build(station_spacing),
      MatSagittaBoundMm2::build(0.02), 32);
}

MatCertificateV1 replay(const CanonicalReachInput2 &input,
                        const MatNumericMatTable2 &numeric,
                        const std::string &record) {
  const MatClearanceProfileGraph2 exact =
      canonical_l_shape_mat_clearance_graph(input, CORE::BigRat(1, 4));
  return replay_mat_certificate_v1(exact, canonical_mat_site_catalog(input),
                                   numeric.center_domain_digest, record);
}

bool field_twenty_replays_exactly() {
  const CanonicalReachInput2 input = l_shape_input();
  const MatNumericMatTable2 numeric = table(input, 0.75);
  const MatCertificateV1 verified =
      replay(input, numeric, numeric.mat_certificate);
  return numeric.mat_certificate.size() == 124796 &&
         verified.canonical_bytes() == numeric.mat_certificate &&
         verified.canonical_digest() ==
             bytes_from_hex("59ccf145c161b819ee91b047edba687e"
                            "b9b4f0a697a4ed49574a8692f3e986dd");
}

bool certificate_is_input_reversal_invariant() {
  return table(l_shape_input(false), 0.75).mat_certificate ==
         table(l_shape_input(true), 0.75).mat_certificate;
}

bool certificate_is_refinement_invariant() {
  const CanonicalReachInput2 input = l_shape_input();
  const MatNumericMatTable2 coarse = table(input, 0.75);
  const MatNumericMatTable2 fine = table(input, 0.25);
  return coarse.proposal.samples.sample_centers.size() !=
             fine.proposal.samples.sample_centers.size() &&
         coarse.mat_certificate == fine.mat_certificate;
}

bool mutated_certificate_fails_loudly() {
  const CanonicalReachInput2 input = l_shape_input();
  const MatNumericMatTable2 numeric = table(input, 0.75);
  std::string mutated = numeric.mat_certificate;
  mutated[mutated.size() / 2] ^= 1;
  try {
    static_cast<void>(replay(input, numeric, mutated));
  } catch (const InvalidMatCertificateReplayError &) {
    return true;
  }
  return false;
}

bool truncated_certificate_fails_loudly() {
  const CanonicalReachInput2 input = l_shape_input();
  const MatNumericMatTable2 numeric = table(input, 0.75);
  std::string truncated = numeric.mat_certificate;
  truncated.pop_back();
  try {
    static_cast<void>(replay(input, numeric, truncated));
  } catch (const InvalidMatCertificateReplayError &) {
    return true;
  }
  return false;
}

bool mismatched_center_domain_fails_loudly() {
  const CanonicalReachInput2 input = l_shape_input();
  MatNumericMatTable2 numeric = table(input, 0.75);
  numeric.center_domain_digest.front() ^= 1;
  try {
    static_cast<void>(replay(input, numeric, numeric.mat_certificate));
  } catch (const InvalidMatCertificateReplayError &) {
    return true;
  }
  return false;
}

bool invalid_center_domain_digest_fails_loudly() {
  const CanonicalReachInput2 input = l_shape_input();
  const MatClearanceProfileGraph2 exact =
      canonical_l_shape_mat_clearance_graph(input, CORE::BigRat(1, 4));
  try {
    static_cast<void>(certified_mat_exact_projection_v1(
        exact, canonical_mat_site_catalog(input), std::string(31, '\0')));
  } catch (const InvalidMatCenterDomainDigestError &) {
    return true;
  }
  return false;
}

} // namespace

bool mat_certificate_gate() {
  return field_twenty_replays_exactly() &&
         certificate_is_input_reversal_invariant() &&
         certificate_is_refinement_invariant() &&
         mutated_certificate_fails_loudly() &&
         truncated_certificate_fails_loudly() &&
         mismatched_center_domain_fails_loudly() &&
         invalid_center_domain_digest_fails_loudly();
}
