#pragma once

#include "reachable_input_2.h"

#include <stdexcept>
#include <string>

class InvalidReachableDomainCertificateIdentityError
    : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

class ReachableDomain2;

class ReachableDomainCertificateIdentity2 {
public:
  const std::string &canonical_bytes() const noexcept;
  const std::string &certificate_digest() const noexcept;
  const std::string &center_domain_digest() const noexcept;

  bool operator==(const ReachableDomainCertificateIdentity2 &) const = default;

private:
  ReachableDomainCertificateIdentity2(std::string canonical_bytes,
                                      std::string certificate_digest,
                                      std::string center_domain_digest);

  std::string canonical_bytes_;
  std::string certificate_digest_;
  std::string center_domain_digest_;

  friend ReachableDomainCertificateIdentity2
  reachable_domain_certificate_identity(const CanonicalReachInput2 &input,
                                        const ReachableDomain2 &domain);
};

ReachableDomainCertificateIdentity2
reachable_domain_certificate_identity(const CanonicalReachInput2 &input,
                                      const ReachableDomain2 &domain);
