#include "coverage_certificate_2/regions_2.h"

#include "coverage_certificate_2/errors_2.h"

#include <string>

namespace compas_cgal::coverage_certificate {

namespace {
/// Compare an interval-shadowed squared distance with an interval-shadowed
/// squared radius; exact on overlap (or always, in EXACT_ONLY mode).
CGAL::Comparison_result compare_squared(const Interval& d2, const Interval& r2, const Point_2& q, const Point_2& c,
                                        const FT& exact_r2, const FilterMode mode, DecisionCounters& counters)
{
    if (mode == FilterMode::Lazy) {
        if (d2.sup() < r2.inf()) return CGAL::SMALLER;
        if (d2.inf() > r2.sup()) return CGAL::LARGER;
    }
    ++counters.exact_membership;
    return CGAL::compare(CGAL::squared_distance(q, c), exact_r2);
}
}  // namespace

CircleRegionSet CircleRegionSet::build(const std::vector<std::array<double, 2>>& centres,
                                       const std::vector<double>& path_radii, const double tool_radius)
{
    if (centres.size() != path_radii.size()) {
        throw InvalidCertificateInputError("full-circle centres and path radii differ in length");
    }
    const FT r = inject_binary64(tool_radius);
    if (CGAL::sign(r) != CGAL::POSITIVE) throw InvalidCertificateInputError("tool radius must be positive");
    CircleRegionSet set;
    set.regions.reserve(centres.size());
    set.circles.reserve(2 * centres.size());
    set.first_circle.reserve(centres.size() + 1);
    for (std::size_t i = 0; i < centres.size(); ++i) {
        const FT x = inject_binary64(centres[i][0]);
        const FT y = inject_binary64(centres[i][1]);
        const FT big_r = inject_binary64(path_radii[i]);
        if (CGAL::sign(big_r) != CGAL::POSITIVE) {
            throw InvalidCertificateInputError("full-circle path radius " + std::to_string(i) + " must be positive");
        }
        const Point_2 c(x, y);
        const bool annulus = CGAL::compare(big_r, r) == CGAL::LARGER;
        const FT outer = big_r + r;
        const FT inner = big_r - r;
        const FT outer_sq = outer * outer;
        const FT inner_sq = annulus ? inner * inner : FT(0);
        set.regions.push_back({c, outer_sq, inner_sq, annulus, interval_of(x), interval_of(y), interval_of(outer_sq),
                               interval_of(inner_sq), static_cast<std::uint32_t>(i),
                               annulus ? PieceKind::AnnulusOuter : PieceKind::Disk});
        set.first_circle.push_back(static_cast<std::uint32_t>(set.circles.size()));
        const auto add = [&](const PieceKind kind, const FT& radius) {
            const FT radius_sq = radius * radius;
            const Circle_2 circle(c, radius_sq);
            set.circles.push_back({static_cast<std::uint32_t>(i), kind, circle, radius, interval_of(x), interval_of(y),
                                   interval_of(radius_sq), circle.bbox()});
        };
        add(annulus ? PieceKind::AnnulusOuter : PieceKind::Disk, outer);
        if (annulus) add(PieceKind::AnnulusInner, inner);
    }
    set.first_circle.push_back(static_cast<std::uint32_t>(set.circles.size()));
    return set;
}

void CircleRegionSet::add_capsule_caps(const std::vector<Point_2>& centres, const std::vector<std::uint32_t>& report_index,
                                       const FT& tool_radius)
{
    if (centres.size() != report_index.size()) throw InvalidCertificateInputError("capsule caps and report indices differ");
    const FT r_sq = tool_radius * tool_radius;
    first_circle.pop_back();
    for (std::size_t i = 0; i < centres.size(); ++i) {
        const Point_2& c = centres[i];
        const auto region = static_cast<std::uint32_t>(regions.size());
        regions.push_back({c, r_sq, FT(0), false, interval_of(c.x()), interval_of(c.y()), interval_of(r_sq), interval_of(FT(0)),
                           report_index[i], PieceKind::CapsuleCap});
        first_circle.push_back(static_cast<std::uint32_t>(circles.size()));
        const Circle_2 circle(c, r_sq);
        circles.push_back({region, PieceKind::CapsuleCap, circle, tool_radius, interval_of(c.x()), interval_of(c.y()),
                           interval_of(r_sq), circle.bbox()});
    }
    first_circle.push_back(static_cast<std::uint32_t>(circles.size()));
}

void CircleRegionSet::add_target_circles(const std::vector<Point_2>& centres, const std::vector<std::uint32_t>& member,
                                         const FT& tool_radius)
{
    if (centres.size() != member.size()) throw InvalidCertificateInputError("target circles and members differ");
    const FT r_sq = tool_radius * tool_radius;
    for (std::size_t i = 0; i < centres.size(); ++i) {
        const Point_2& c = centres[i];
        const Circle_2 circle(c, r_sq);
        circles.push_back({member[i], PieceKind::TargetCircle, circle, tool_radius, interval_of(c.x()), interval_of(c.y()),
                           interval_of(r_sq), circle.bbox()});
    }
}

std::uint32_t CircleRegionSet::add_target_circle(const Point_2& centre, const FT& radius, const std::uint32_t member)
{
    if (CGAL::sign(radius) != CGAL::POSITIVE) throw InvalidCertificateInputError("target circle radius must be positive");
    const FT r_sq = radius * radius;
    const Circle_2 circle(centre, r_sq);
    const auto id = static_cast<std::uint32_t>(circles.size());
    circles.push_back({member, PieceKind::TargetCircle, circle, radius, interval_of(centre.x()), interval_of(centre.y()),
                       interval_of(r_sq), circle.bbox()});
    return id;
}

std::uint32_t CircleRegionSet::add_support_region(const Point_2& centre, const FT& outer, const FT& inner,
                                                const std::uint32_t report_index)
{
    if (CGAL::sign(inner) == CGAL::NEGATIVE || CGAL::compare(outer, inner) != CGAL::LARGER) {
        throw InvalidCertificateInputError("a support region needs 0 <= inner < outer");
    }
    const bool annulus = CGAL::sign(inner) == CGAL::POSITIVE;
    const FT outer_sq = outer * outer;
    const FT inner_sq = annulus ? inner * inner : FT(0);
    const auto region = static_cast<std::uint32_t>(regions.size());
    CircleRegion added{centre, outer_sq, inner_sq, annulus, interval_of(centre.x()), interval_of(centre.y()),
                       interval_of(outer_sq), interval_of(inner_sq), report_index,
                       annulus ? PieceKind::AnnulusOuter : PieceKind::Disk};
    added.covers = false;
    regions.push_back(std::move(added));
    first_circle.pop_back();
    first_circle.push_back(static_cast<std::uint32_t>(circles.size()));
    const auto add = [&](const PieceKind kind, const FT& radius, const FT& radius_sq) {
        const Circle_2 circle(centre, radius_sq);
        circles.push_back({region, kind, circle, radius, interval_of(centre.x()), interval_of(centre.y()),
                           interval_of(radius_sq), circle.bbox()});
    };
    add(annulus ? PieceKind::AnnulusOuter : PieceKind::Disk, outer, outer_sq);
    if (annulus) add(PieceKind::AnnulusInner, inner, inner_sq);
    first_circle.push_back(static_cast<std::uint32_t>(circles.size()));
    return region;
}

Side membership(const CircleRegion& region, const Point_2& q, const Interval& qx, const Interval& qy,
                const FilterMode mode, DecisionCounters& counters)
{
    const Interval dx = qx - region.centre_x;
    const Interval dy = qy - region.centre_y;
    const Interval d2 = dx * dx + dy * dy;
    const auto outer = compare_squared(d2, region.outer_squared_i, q, region.centre, region.outer_squared, mode, counters);
    if (outer == CGAL::LARGER) return Side::Exterior;
    if (outer == CGAL::EQUAL) return Side::Boundary;
    if (!region.annulus) return Side::Interior;
    const auto inner = compare_squared(d2, region.inner_squared_i, q, region.centre, region.inner_squared, mode, counters);
    if (inner == CGAL::SMALLER) return Side::Exterior;
    if (inner == CGAL::EQUAL) return Side::Boundary;
    return Side::Interior;
}

Point_2 stereographic_point(const BoundaryCircle& circle, const FT& s)
{
    const FT one(1);
    const FT den = one + s * s;
    const auto& c = circle.circle.center();
    return Point_2(c.x() + circle.radius * (one - s * s) / den, c.y() + circle.radius * (FT(2) * s) / den);
}

}  // namespace compas_cgal::coverage_certificate
