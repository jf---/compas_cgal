#pragma once

// Number and traits types of the exact Held centre domain (CORE-free).

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/Gps_circle_segment_traits_2.h>

namespace compas_cgal::held_centre_domain {

using Kernel = CGAL::Epeck;
using FT = Kernel::FT;
using RationalPoint = Kernel::Point_2;
using Traits = CGAL::Gps_circle_segment_traits_2<Kernel>;
using Gps = CGAL::General_polygon_set_2<Traits>;
using Polygon = Traits::Polygon_2;
using OneRootPoint = Traits::Point_2;
using XCurve = Traits::X_monotone_curve_2;
using Curve = Traits::Curve_2;

}  // namespace compas_cgal::held_centre_domain
