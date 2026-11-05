#pragma once

// CCW order of one-root points on a tool rim, decided exactly: the upper half
// (y above the centre, plus the rightmost point) comes first in decreasing x,
// then the lower half (plus the leftmost point) in increasing x. compare_xy alone
// does not give angular order. Points are compared only, never combined, so two
// points may carry different roots.

#include "exact_motion_2.h"
#include "stock_2.h"

#include <CGAL/enum.h>
#include <CGAL/number_utils.h>

inline bool rim_upper(const GpsPoint& point, const EPoint& center)
{
    using CoordNT = GpsPoint::CoordNT;
    const CoordNT center_x(center.x());
    const CoordNT center_y(center.y());
    const CGAL::Comparison_result y = CGAL::compare(point.y(), center_y);
    if (y == CGAL::LARGER) return true;
    if (y == CGAL::SMALLER) return false;
    return CGAL::compare(point.x(), center_x) == CGAL::LARGER;
}

inline bool rim_ccw_less(const GpsPoint& lhs, const GpsPoint& rhs, const EPoint& center)
{
    const bool lhs_upper = rim_upper(lhs, center);
    const bool rhs_upper = rim_upper(rhs, center);
    if (lhs_upper != rhs_upper) return lhs_upper;
    const CGAL::Comparison_result x = CGAL::compare(lhs.x(), rhs.x());
    return lhs_upper ? x == CGAL::LARGER : x == CGAL::SMALLER;
}

inline bool strictly_inside_ccw_arc(
    const GpsPoint& point,
    const GpsPoint& start,
    const GpsPoint& end,
    const EPoint& center)
{
    if (point == start || point == end) return false;
    if (!rim_ccw_less(end, start, center)) {
        return rim_ccw_less(start, point, center)
            && rim_ccw_less(point, end, center);
    }
    return rim_ccw_less(start, point, center)
        || rim_ccw_less(point, end, center);
}
