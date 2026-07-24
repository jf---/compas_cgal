// THE harvest, validated: drive Arrangement_zone_2 over the cutter's x-monotone
// arcs with a visitor that collects sub-arcs whose containing face is material
// (contained()); sum their spans = engaged TEA -- read LOCALLY off the Gps's own
// arrangement, no overlay. Validate against known answers (full immersion = 2*pi,
// edge-straddle = pi).
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Gps_circle_segment_traits_2.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/Arrangement_zone_2.h>
#include <boost/variant.hpp>
#include <cmath>
#include <numbers>
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

template <class ArrT>
struct EngagementVisitor {
    using X_monotone_curve_2 = typename ArrT::X_monotone_curve_2;
    using Vertex_handle = typename ArrT::Vertex_handle;
    using Halfedge_handle = typename ArrT::Halfedge_handle;
    using Face_handle = typename ArrT::Face_handle;
    using Result = std::pair<Halfedge_handle, bool>;
    std::vector<X_monotone_curve_2> engaged;
    void init(ArrT*) {}
    Result found_subcurve(const X_monotone_curve_2& cv, Face_handle face,
                          Vertex_handle, Halfedge_handle, Vertex_handle, Halfedge_handle) {
        if (face->contained()) engaged.push_back(cv);   // cutter rim in material
        return Result(Halfedge_handle(), false);
    }
    Result found_overlap(const X_monotone_curve_2&, Halfedge_handle, Vertex_handle, Vertex_handle) {
        return Result(Halfedge_handle(), false);         // grazing boundary: measure-zero
    }
};

static GpsPolygon disk_polygon(const EPoint& c, const FT& r) {
    ECircle circle(c, r * r, CGAL::COUNTERCLOCKWISE);
    GpsPoint pmin(c.x() - r, c.y()), pmax(c.x() + r, c.y());
    GpsPolygon p;
    p.push_back(GpsXCurve(circle, pmin, pmax, CGAL::COUNTERCLOCKWISE));
    p.push_back(GpsXCurve(circle, pmax, pmin, CGAL::COUNTERCLOCKWISE));
    return p;
}
static GpsPolygon square_polygon() {
    EPoint a(0, 0), b(10, 0), c(10, 10), d(0, 10);
    GpsPolygon p;
    p.push_back(GpsXCurve(a, b)); p.push_back(GpsXCurve(b, c));
    p.push_back(GpsXCurve(c, d)); p.push_back(GpsXCurve(d, a));
    return p;
}

static double engaged_tea(Gps& stock, double cx, double cy, double r) {
    Arr& arr = stock.arrangement();
    PL pl(arr);
    Curve cutter(ECircle(EPoint(cx, cy), FT(r * r)));
    Traits traits;
    auto mkx = traits.make_x_monotone_2_object();
    std::vector<boost::variant<GpsPoint, GpsXCurve>> xs;
    mkx(cutter, std::back_inserter(xs));
    EngagementVisitor<Arr> vis;
    for (auto& xo : xs)
        if (auto xc = boost::get<GpsXCurve>(&xo)) {
            CGAL::Arrangement_zone_2<Arr, EngagementVisitor<Arr>> zd(arr, &vis);
            zd.init(*xc, pl);
            zd.compute_zone();
        }
    double tea = 0.0;
    for (auto& cv : vis.engaged) {
        double sx = CGAL::to_double(cv.source().x()), sy = CGAL::to_double(cv.source().y());
        double tx = CGAL::to_double(cv.target().x()), ty = CGAL::to_double(cv.target().y());
        if (cv.orientation() == CGAL::CLOCKWISE) { std::swap(sx, tx); std::swap(sy, ty); }
        double span = std::atan2(ty - cy, tx - cx) - std::atan2(sy - cy, sx - cx);
        if (span <= 0.0) span += 2.0 * std::numbers::pi;
        tea += span;
    }
    return tea;
}

int main() {
    Gps stock;
    stock.insert(square_polygon());
    Gps disk;
    disk.insert(disk_polygon(EPoint(5, 5), FT(2)));
    stock.difference(disk);

    double pi = std::numbers::pi;
    std::cout << "full immersion (2,2 r0.5):  TEA=" << engaged_tea(stock, 2, 2, 0.5) << "  (want " << 2 * pi << ")\n";
    std::cout << "straddle left edge (0,5 r1): TEA=" << engaged_tea(stock, 0, 5, 1.0) << "  (want " << pi << ")\n";
    std::cout << "over the disk centre (5,5 r0.5): TEA=" << engaged_tea(stock, 5, 5, 0.5) << "  (want 0 -- fully cleared)\n";
    return 0;
}
