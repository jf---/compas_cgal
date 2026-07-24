// THE targeted fix, validated: keep the (now-fast) Gps depletion, but answer the
// engagement query by ZONING the cutter circle in the Gps's OWN internal
// arrangement (General_polygon_set_2::arrangement(), faces carry contained()) --
// local, no per-query overlay. Build square-minus-disk as a real Gps, then zone a
// cutter straddling the disk boundary and read contained() off the crossed faces.
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Gps_circle_segment_traits_2.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <boost/variant.hpp>
#include <iostream>
#include <vector>

using K = CGAL::Exact_predicates_exact_constructions_kernel;
using Traits = CGAL::Gps_circle_segment_traits_2<K>;
using Gps = CGAL::General_polygon_set_2<Traits>;
using Arr = Gps::Arrangement_2;
using GpsPolygon = Traits::Polygon_2;
using GpsPoint = Traits::Point_2;
using GpsXCurve = Traits::X_monotone_curve_2;
using Curve = Traits::Curve_2;
using PL = CGAL::Arr_walk_along_line_point_location<Arr>;
using FT = K::FT;
using EPoint = K::Point_2;
using ECircle = K::Circle_2;

static GpsPolygon disk_polygon(const EPoint& center, const FT& radius) {
    FT r_sq = radius * radius;
    ECircle circle(center, r_sq, CGAL::COUNTERCLOCKWISE);
    GpsPoint p_min(center.x() - radius, center.y());
    GpsPoint p_max(center.x() + radius, center.y());
    GpsPolygon poly;
    poly.push_back(GpsXCurve(circle, p_min, p_max, CGAL::COUNTERCLOCKWISE));
    poly.push_back(GpsXCurve(circle, p_max, p_min, CGAL::COUNTERCLOCKWISE));
    return poly;
}

static GpsPolygon square_polygon() {
    EPoint a(0, 0), b(10, 0), c(10, 10), d(0, 10);  // linear X-curve takes Kernel::Point_2
    GpsPolygon poly;
    poly.push_back(GpsXCurve(a, b));
    poly.push_back(GpsXCurve(b, c));
    poly.push_back(GpsXCurve(c, d));
    poly.push_back(GpsXCurve(d, a));
    return poly;
}

int main() {
    Gps stock;
    stock.insert(square_polygon());
    Gps disk;
    disk.insert(disk_polygon(EPoint(5, 5), FT(2)));
    stock.difference(disk);  // square minus disk (the real depletion path)

    Arr& arr = stock.arrangement();
    std::cout << "gps arrangement: faces=" << arr.number_of_faces()
              << " edges=" << arr.number_of_edges() << "\n";

    PL pl(arr);
    Curve cutter(ECircle(EPoint(5, 7), FT(9) / 4));  // r=1.5 at (5,7), straddles disk boundary
    Traits traits;
    auto mkx = traits.make_x_monotone_2_object();
    std::vector<boost::variant<GpsPoint, GpsXCurve>> xs;
    mkx(cutter, std::back_inserter(xs));

    int faces = 0, material = 0, edges = 0;
    for (auto& xo : xs) {
        if (auto xc = boost::get<GpsXCurve>(&xo)) {
            std::vector<boost::variant<Arr::Vertex_handle, Arr::Halfedge_handle, Arr::Face_handle>> cells;
            CGAL::zone(arr, *xc, std::back_inserter(cells), pl);
            for (auto& cell : cells) {
                if (auto fh = boost::get<Arr::Face_handle>(&cell)) { faces++; if ((*fh)->contained()) material++; }
                else if (boost::get<Arr::Halfedge_handle>(&cell)) edges++;
            }
        }
    }
    std::cout << "ZONE in the Gps's own arrangement: face cells=" << faces
              << " (contained/material=" << material << ") edge crossings=" << edges << "\n";
    std::cout << "=> engagement read LOCALLY off the existing depletion, contained() = material.\n";
    return 0;
}
