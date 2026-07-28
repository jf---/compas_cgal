#include "reachable_certificate_encoding_2.h"
#include "reachable_domain_2.h"

#include <array>
#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

namespace {

compas::RowMatrixXd matrix(const std::vector<std::array<double, 2>> &points) {
  compas::RowMatrixXd result(static_cast<Eigen::Index>(points.size()), 2);
  for (std::size_t index = 0; index < points.size(); ++index) {
    result(static_cast<Eigen::Index>(index), 0) = points[index][0];
    result(static_cast<Eigen::Index>(index), 1) = points[index][1];
  }
  return result;
}

CanonicalReachInput2 l_shape_input(const bool reversed,
                                   const double tool_radius = 0.5) {
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
                               {}, tool_radius);
}

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

ReachableDomainCertificateIdentity2
l_shape_identity(const CanonicalReachInput2 &input) {
  const ReachableDomain2 domain = ReachableDomain2::build(input);
  return reachable_domain_certificate_identity(input, domain);
}

bool identity_matches_python_golden() {
  const ReachableDomainCertificateIdentity2 identity =
      l_shape_identity(l_shape_input(false));
  return identity.canonical_bytes().size() == 48947 &&
         identity.certificate_digest() ==
             bytes_from_hex("77b7b15003f2d5de9fbee7a495d4d490"
                            "04d958a58d10a9c73a9265e3a98d5874") &&
         identity.center_domain_digest() ==
             bytes_from_hex("98aaa87a6fdc0589ef364e7ae3edad561"
                            "208ff8d4245b340d616c63878be12de");
}

bool identity_is_input_reversal_invariant() {
  return l_shape_identity(l_shape_input(false)) ==
         l_shape_identity(l_shape_input(true));
}

bool mismatched_input_fails_loudly() {
  const CanonicalReachInput2 input = l_shape_input(false);
  const ReachableDomain2 domain = ReachableDomain2::build(input);
  try {
    static_cast<void>(reachable_domain_certificate_identity(
        l_shape_input(false, 0.75), domain));
  } catch (const InvalidReachableDomainCertificateIdentityError &) {
    return true;
  }
  return false;
}

} // namespace

bool domain_identity_gate() {
  return identity_matches_python_golden() &&
         identity_is_input_reversal_invariant() &&
         mismatched_input_fails_loudly();
}
