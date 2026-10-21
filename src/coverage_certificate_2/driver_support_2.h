#pragma once

// Shared plumbing of the certificate drivers (certificate_2.cpp,
// exact_certificate_2.cpp): box overlaps, flank-piece boxes, piece kinds.

#include "coverage_certificate_2/certificate_2.h"
#include "coverage_certificate_2/lines_2.h"
#include "coverage_certificate_2/regions_2.h"

#include <CGAL/Bbox_2.h>

#include <cstdint>
#include <vector>

namespace compas_cgal::coverage_certificate {

/// Every (a, b) whose certified boxes overlap (a superset filter), sorted so
/// that the candidate lists are a function of the input alone.
[[nodiscard]] std::vector<std::vector<std::uint32_t>> overlaps(const std::vector<CGAL::Bbox_2>& as,
                                                               const std::vector<CGAL::Bbox_2>& bs);

/// Certified box of the flank segment x([t_lo, t_hi]), widened by the margin.
[[nodiscard]] CGAL::Bbox_2 flank_box(const LineSupport& line, const FT& t_lo, const FT& t_hi);

/// The certificate piece kind of a boundary circle.
[[nodiscard]] CertificatePieceKind circle_kind(PieceKind kind);

}  // namespace compas_cgal::coverage_certificate
