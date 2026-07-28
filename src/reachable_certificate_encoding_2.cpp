#include "reachable_certificate_encoding_2.h"

#include "canonical_encoding.h"
#include "continuous_tea_2/sha256.h"
#include "reachable_domain_2.h"

#include <algorithm>
#include <array>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

constexpr std::string_view REACHABLE_STRATEGY_VERSION =
    "exact-reachable-arrangement-v2";

void require_sorted_unique_nonempty(const std::vector<std::string> &records,
                                    const std::string_view name) {
  if (records.empty() || !std::is_sorted(records.begin(), records.end()) ||
      std::adjacent_find(records.begin(), records.end()) != records.end()) {
    throw InvalidReachableDomainCertificateIdentityError(
        std::string(name) + " must be nonempty, sorted, and unique");
  }
}

std::string canonical_point_bytes(const std::array<double, 2> &point) {
  return canonical_encode_tagged_union("point2-world-xy-v1",
                                       canonical_encode_sequence({
                                           canonical_encode_binary64(point[0]),
                                           canonical_encode_binary64(point[1]),
                                       }));
}

std::string canonical_ring_bytes(const CanonicalReachRing2 &ring) {
  std::vector<std::string> points;
  points.reserve(ring.binary64_points.size());
  for (const std::array<double, 2> &point : ring.binary64_points) {
    points.push_back(canonical_point_bytes(point));
  }
  return canonical_encode_tagged_union(ring.outer ? "outer-ring-ccw-v1"
                                                  : "hole-ring-cw-v1",
                                       canonical_encode_sequence(points));
}

std::string canonical_polygon_bytes(const CanonicalReachInput2 &input) {
  std::vector<std::string> holes;
  holes.reserve(input.holes.size());
  for (const CanonicalReachRing2 &hole : input.holes) {
    holes.push_back(canonical_ring_bytes(hole));
  }
  std::sort(holes.begin(), holes.end());
  if (std::adjacent_find(holes.begin(), holes.end()) != holes.end()) {
    throw InvalidReachableDomainCertificateIdentityError(
        "canonical polygon holes must be unique");
  }
  return canonical_encode_tagged_union(
      "polygon-with-holes-world-xy-v1",
      canonical_encode_component_map({
          {
              "holes",
              canonical_encode_sequence(holes),
          },
          {
              "outer-ring",
              canonical_ring_bytes(input.outer),
          },
      }));
}

void validate_certificate_identity_inputs(
    const CanonicalReachInput2 &input,
    const ReachableDomainCertificate2 &certificate) {
  validate_canonical_reach_input(input);
  if (certificate.input_recipe_record != input.recipe_record ||
      certificate.strategy_version != REACHABLE_STRATEGY_VERSION ||
      !certificate.exact_cell_selection ||
      !certificate.complete_source_provenance ||
      !certificate.reachable_subset_of_design) {
    throw InvalidReachableDomainCertificateIdentityError(
        "reachable-domain certificate does not bind its exact input");
  }
  require_sorted_unique_nonempty(certificate.source_curve_records,
                                 "source curve records");
  require_sorted_unique_nonempty(certificate.arrangement_vertex_records,
                                 "arrangement vertex records");
  require_sorted_unique_nonempty(certificate.selected_cell_records,
                                 "selected cell records");
  if (certificate.component_records.size() != 1) {
    throw InvalidReachableDomainCertificateIdentityError(
        "reachable-domain certificate must bind one component");
  }
  require_sorted_unique_nonempty(certificate.component_records,
                                 "component records");
}

std::string center_domain_digest(const std::string &certificate_digest) {
  return sha256_bytes(canonical_encode_tagged_union(
      "exact-region-center-domain-v1",
      canonical_encode_component_map({
          {
              "reachable-domain-certificate-digest",
              canonical_encode_bytes(certificate_digest),
          },
      })));
}

} // namespace

ReachableDomainCertificateIdentity2::ReachableDomainCertificateIdentity2(
    std::string canonical_bytes, std::string certificate_digest,
    std::string center_domain_digest)
    : canonical_bytes_(std::move(canonical_bytes)),
      certificate_digest_(std::move(certificate_digest)),
      center_domain_digest_(std::move(center_domain_digest)) {}

const std::string &
ReachableDomainCertificateIdentity2::canonical_bytes() const noexcept {
  return canonical_bytes_;
}

const std::string &
ReachableDomainCertificateIdentity2::certificate_digest() const noexcept {
  return certificate_digest_;
}

const std::string &
ReachableDomainCertificateIdentity2::center_domain_digest() const noexcept {
  return center_domain_digest_;
}

ReachableDomainCertificateIdentity2
reachable_domain_certificate_identity(const CanonicalReachInput2 &input,
                                      const ReachableDomain2 &domain) {
  const ReachableDomainCertificate2 certificate = domain.certificate();
  validate_certificate_identity_inputs(input, certificate);
  std::string canonical_bytes = canonical_encode_tagged_union(
      "reachable-domain-certificate-v2",
      canonical_encode_component_map({
          {
              "complete-source-provenance",
              canonical_encode_boolean(certificate.complete_source_provenance),
          },
          {
              "component-records",
              canonical_encode_sequence(certificate.component_records),
          },
          {
              "design",
              canonical_polygon_bytes(input),
          },
          {
              "exact-cell-selection",
              canonical_encode_boolean(certificate.exact_cell_selection),
          },
          {
              "arrangement-vertex-records",
              canonical_encode_sequence(certificate.arrangement_vertex_records),
          },
          {
              "reachable-subset-of-design",
              canonical_encode_boolean(certificate.reachable_subset_of_design),
          },
          {
              "selected-cell-records",
              canonical_encode_sequence(certificate.selected_cell_records),
          },
          {
              "source-curve-records",
              canonical_encode_sequence(certificate.source_curve_records),
          },
          {
              "strategy-version",
              canonical_encode_bytes(certificate.strategy_version),
          },
          {
              "tool-radius",
              canonical_encode_binary64(input.binary64_radius),
          },
      }));
  std::string certificate_digest = sha256_bytes(canonical_bytes);
  return ReachableDomainCertificateIdentity2(
      std::move(canonical_bytes), certificate_digest,
      center_domain_digest(certificate_digest));
}
