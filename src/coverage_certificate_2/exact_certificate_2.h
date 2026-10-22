#pragma once

// The certificate against the exact target T = Omega* (+) D_r (exact_target_2.h).
//
// Theorem (docs/held_coverage_certificate.md, "Exact target Omega*"): let the
// swept pieces (cut circles, capsule caps, connector flanks) be cut at every
// swept boundary and at every target cut support (declared boundary lines and
// circles, vertex circles |x - w| = r), which contain dT. Then
//
//     T subset U   <=>   no swept edge is exposed   and   p in U,
//
// for any one point p in int T. An edge e on dS_i is exposed when no other sweep
// covers its outer side and that side lies in T. int T is connected (Omega* is
// one component, a neck splitting it is a hard error upstream, and
// Omega* (+) open D_r is connected and dense in int T), so an uncovered face
// either borders a swept boundary (an exposed edge) or the uncovered set is all
// of int T (then p is uncovered). No target piece is swept.
//
// The target side of an exposed edge is decided at its exact witness w by closed
// membership w in T: w lies on no cut support (the edge is cut at all of them),
// so w is not on dT, and w in T means the whole neighbourhood is. A swept piece
// may lie ON a target cut support (a coincident support: exact connectors along
// dOmega* do by design); then its target side is decided at the exact side probe
// beside w (side_probe_2.h).
//
// Exact dOmega* connectors (exact_connectors_2.h) enter as composite regions:
// their supports cut every piece, their membership is decided exactly at the
// edge witness, and their flank segments and support-circle arcs are checked
// pieces.

#include "coverage_certificate_2/certificate_2.h"
#include "coverage_certificate_2/exact_connectors_2.h"
#include "coverage_certificate_2/exact_target_2.h"
#include "coverage_certificate_2/lines_2.h"
#include "coverage_certificate_2/regions_2.h"
#include "coverage_certificate_2/target_2.h"

#include <vector>

namespace compas_cgal::coverage_certificate {

/// Run the certificate against the exact target.
///
/// Args:
///     cuts: full-circle cuts (capsule caps and target cut circles are appended).
///     path_domain: the polygonal centre domain the connectors follow (their
///         sites); it does not define the target.
///     connectors: connectors along path_domain.
///     exact_connectors: exact dOmega* connectors (sub-pieces of the target's pieces).
///     target: the exact target T = Omega* (+) D_r.
///     audit: also confirm every edge the counters call covered (slow).
///
/// Raises:
///     ConnectorBoundaryCoverageError: the exact connectors do not cover dOmega*.
///     CoincidentTargetSupportError: an irrational (polygon connector) flank lies
///         on an exact connector support, where no rational side probe exists.
///     CoverageCertificateNativeError: an internal contract breach.
[[nodiscard]] CertificateResult run_exact_target_certificate(CircleRegionSet cuts, const CentreDomain& path_domain,
                                                             const std::vector<ConnectorSpec>& connectors,
                                                             const std::vector<ExactConnectorSpec>& exact_connectors,
                                                             const ExactTarget& target, FilterMode mode, bool audit);

}  // namespace compas_cgal::coverage_certificate
