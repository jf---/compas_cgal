#include "coverage_certificate_2/driver_support_2.h"

#include "coverage_certificate_2/errors_2.h"
#include "coverage_certificate_2/line_sweep_2.h"

#include <CGAL/box_intersection_d.h>

#include <algorithm>

namespace compas_cgal::coverage_certificate {

namespace {
using Box = CGAL::Box_intersection_d::Box_with_info_d<double, 2, std::uint32_t>;
}  // namespace

std::vector<std::vector<std::uint32_t>> overlaps(const std::vector<CGAL::Bbox_2>& as, const std::vector<CGAL::Bbox_2>& bs)
{
    std::vector<std::vector<std::uint32_t>> out(as.size());
    if (as.empty() || bs.empty()) return out;
    std::vector<Box> ba, bb;
    ba.reserve(as.size());
    bb.reserve(bs.size());
    for (std::uint32_t i = 0; i < as.size(); ++i) ba.emplace_back(as[i], i);
    for (std::uint32_t i = 0; i < bs.size(); ++i) bb.emplace_back(bs[i], i);
    CGAL::box_intersection_d(ba.begin(), ba.end(), bb.begin(), bb.end(),
                             [&](const Box& a, const Box& b) { out[a.info()].push_back(b.info()); });
    for (auto& list : out) {
        std::sort(list.begin(), list.end());
        list.erase(std::unique(list.begin(), list.end()), list.end());
    }
    return out;
}

/// Certified box of the flank segment x([t_lo, t_hi]), widened by the margin.
CGAL::Bbox_2 flank_box(const LineSupport& line, const FT& t_lo, const FT& t_hi)
{
    const Interval dx = line.ny_i, dy = -line.nx_i;
    const Interval bx = line.b0x_i + line.b1x_i * line.sqrt_l_i;
    const Interval by = line.b0y_i + line.b1y_i * line.sqrt_l_i;
    const Interval x0 = bx + interval_of(t_lo) * dx, x1 = bx + interval_of(t_hi) * dx;
    const Interval y0 = by + interval_of(t_lo) * dy, y1 = by + interval_of(t_hi) * dy;
    const Interval margin(LINE_PIECE_BOX_MARGIN);
    return CGAL::Bbox_2((Interval(std::min(x0.inf(), x1.inf())) - margin).inf(),
                        (Interval(std::min(y0.inf(), y1.inf())) - margin).inf(),
                        (Interval(std::max(x0.sup(), x1.sup())) + margin).sup(),
                        (Interval(std::max(y0.sup(), y1.sup())) + margin).sup());
}

CertificatePieceKind circle_kind(const PieceKind kind)
{
    switch (kind) {
        case PieceKind::AnnulusOuter: return CertificatePieceKind::AnnulusOuter;
        case PieceKind::AnnulusInner: return CertificatePieceKind::AnnulusInner;
        case PieceKind::Disk: return CertificatePieceKind::Disk;
        case PieceKind::CapsuleCap: return CertificatePieceKind::CapsuleCap;
        case PieceKind::TargetCircle: return CertificatePieceKind::TargetBoundary;
    }
    throw CoverageCertificateNativeError("unknown piece kind");
}

}  // namespace compas_cgal::coverage_certificate
