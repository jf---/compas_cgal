#include "cutter_centre_domain_2.h"

#include "reachable_arrangement_2.h"
#include "reachable_errors_2.h"

#include <cmath>
#include <utility>
#include <vector>

#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/enum.h>
#include <CGAL/number_utils.h>
#include <CGAL/squared_distance_2.h>

namespace {

CGAL::Bounded_side ring_side(
    const CanonicalReachRing2& ring,
    const ReachKernelPoint& query)
{
    return CGAL::bounded_side_2(
        ring.points.begin(),
        ring.points.end(),
        query,
        ReachKernel());
}

bool clears_ring_boundary(
    const CanonicalReachRing2& ring,
    const ReachKernelPoint& query,
    const ReachFT& squared_radius)
{
    for (std::size_t index = 0; index < ring.points.size(); ++index) {
        const ReachKernel::Segment_2 edge(
            ring.points[index],
            ring.points[(index + 1) % ring.points.size()]);
        if (CGAL::compare(
                CGAL::squared_distance(query, edge),
                squared_radius)
            == CGAL::SMALLER) {
            return false;
        }
    }
    return true;
}

} // namespace

CutterCentreDomain2::CutterCentreDomain2(
    Eigen::Ref<const compas::RowMatrixXd> design_boundary,
    const std::vector<compas::RowMatrixXd>& holes,
    double tool_radius)
    : CutterCentreDomain2(
        canonical_reach_input(
            design_boundary,
            holes,
            tool_radius))
{
}

CutterCentreDomain2::CutterCentreDomain2(CanonicalReachInput2 input)
    : input_(std::move(input))
    , squared_radius_(input_.radius * input_.radius)
{
    validate_canonical_reach_input(input_);
    static_cast<void>(reachable_design_polygon(input_));
}

bool CutterCentreDomain2::contains(double x, double y) const
{
    if (!std::isfinite(x) || !std::isfinite(y)) {
        throw InvalidReachableDomainInputError(
            "cutter-centre query coordinates must be finite binary64 values");
    }
    const ReachKernelPoint query{ReachFT(x), ReachFT(y)};
    if (ring_side(input_.outer, query) == CGAL::ON_UNBOUNDED_SIDE) {
        return false;
    }
    for (const CanonicalReachRing2& hole : input_.holes) {
        if (ring_side(hole, query) == CGAL::ON_BOUNDED_SIDE) {
            return false;
        }
    }
    if (!clears_ring_boundary(input_.outer, query, squared_radius_)) {
        return false;
    }
    for (const CanonicalReachRing2& hole : input_.holes) {
        if (!clears_ring_boundary(hole, query, squared_radius_)) {
            return false;
        }
    }
    return true;
}
