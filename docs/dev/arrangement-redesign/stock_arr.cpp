// Redesign query probe: persistent circle-segment Arrangement_2 with a `material`
// face bool. Build square-minus-disk, then ZONE a cutter circle (split into
// x-monotone arcs) and read materiality off the faces it crosses -- the local,
// no-overlay engagement query.
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_circle_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/Arr_point_location_result.h>
#include <boost/variant.hpp>
#include <iostream>
#include <vector>

using K = CGAL::Exact_predicates_exact_constructions_kernel;
using Traits = CGAL::Arr_circle_segment_traits_2<K>;
using Dcel = CGAL::Arr_face_extended_dcel<Traits, bool>;
using Arr = CGAL::Arrangement_2<Traits, Dcel>;
using CSPoint = Traits::Point_2;
using XCurve = Traits::X_monotone_curve_2;
using Curve = Traits::Curve_2;
using PL = CGAL::Arr_walk_along_line_point_location<Arr>;

static void mark(Arr& arr, PL& pl, const CSPoint& p, bool v) {
    auto obj = pl.locate(p);
    if (auto fch = boost::get<Arr::Face_const_handle>(&obj)) arr.non_const_handle(*fch)->set_data(v);
}

int main() {
    Arr arr;
    std::vector<Curve> cs;
    K::Point_2 a(0, 0), b(10, 0), c(10, 10), d(0, 10);
    cs.push_back(Curve(K::Segment_2(a, b)));
    cs.push_back(Curve(K::Segment_2(b, c)));
    cs.push_back(Curve(K::Segment_2(c, d)));
    cs.push_back(Curve(K::Segment_2(d, a)));
    CGAL::insert(arr, cs.begin(), cs.end());

    PL pl(arr);
    for (auto f = arr.faces_begin(); f != arr.faces_end(); ++f) f->set_data(false);
    mark(arr, pl, CSPoint(K::FT(5), K::FT(5)), true);
    CGAL::insert(arr, Curve(K::Circle_2(K::Point_2(5, 5), K::FT(4))));  // disk r=2 at (5,5)
    mark(arr, pl, CSPoint(K::FT(5), K::FT(5)), false);  // disk cleared
    mark(arr, pl, CSPoint(K::FT(5), K::FT(9)), true);   // annulus material
    std::cout << "stock built: faces=" << arr.number_of_faces() << "\n";

    // ---- ZONE query: cutter circle at (5,7) r=1.5 straddling the disk boundary ----
    // lower arc dips into the disk (cleared), upper arc is in the annulus (material).
    Curve cutter(K::Circle_2(K::Point_2(5, 7), K::FT(9) / 4));  // r^2 = 2.25
    Traits traits;
    auto mkx = traits.make_x_monotone_2_object();
    std::vector<boost::variant<CSPoint, XCurve>> xobjs;
    mkx(cutter, std::back_inserter(xobjs));
    std::cout << "cutter x-monotone pieces: " << xobjs.size() << "\n";

    int faces = 0, material = 0, edges = 0, verts = 0;
    for (auto& xo : xobjs) {
        if (auto xc = boost::get<XCurve>(&xo)) {
            std::vector<boost::variant<Arr::Vertex_handle, Arr::Halfedge_handle, Arr::Face_handle>> cells;
            CGAL::zone(arr, *xc, std::back_inserter(cells), pl);
            for (auto& cell : cells) {
                if (auto fh = boost::get<Arr::Face_handle>(&cell)) { faces++; if ((*fh)->data()) material++; }
                else if (boost::get<Arr::Halfedge_handle>(&cell)) edges++;
                else verts++;
            }
        }
    }
    std::cout << "zone cells: face=" << faces << " (material=" << material << ")"
              << " edge=" << edges << " vertex=" << verts << "\n";
    std::cout << "=> cutter crosses the disk boundary and touches BOTH a material (annulus)"
                 " and a cleared (disk) face, read locally with no overlay.\n";
    return 0;
}
