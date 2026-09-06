#include "held_disk_contour_2.h"

#include <algorithm>
#include <cmath>

#include <nanobind/nanobind.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/pair.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace {
double report_coordinate(const GpsPoint::CoordNT& value)
{
    if (!value.is_extended() || CGAL::is_zero(value.a0())
        || CGAL::sign(value.a0()) == CGAL::sign(value.a1())) {
        return CGAL::to_double(value);
    }
    // Opposite-sign terms cancel in a0 + a1*sqrt(root). Rationalize before
    // reporting: (a0² - a1²*root) / (a0 - a1*sqrt(root)). The numerator is
    // evaluated exactly in CGAL; the denominator's terms have the same sign.
    const Epeck::FT rational_squared = CGAL::square(value.a0());
    const Epeck::FT radical_squared = CGAL::square(value.a1()) * value.root();
    const Epeck::FT numerator = rational_squared - radical_squared;
    const GpsPoint::CoordNT conjugate(value.a0(), -value.a1(), value.root());
    // Normalize BEFORE conversion: a squared numerator can overflow/underflow
    // while the final coordinate is representable. Twice the larger term makes
    // |numerator/scale| <= |value| and |conjugate/scale| lie in [1/2, 1].
    // All operations remain in this one root extension; neither term cancels.
    using Coord = GpsPoint::CoordNT;
    const Coord dominant = CGAL::compare(rational_squared, radical_squared) != CGAL::SMALLER
        ? Coord(value.a0()) : Coord(Epeck::FT(0), value.a1(), value.root());
    const Coord scale = Coord(Epeck::FT(2)) * dominant;
    return CGAL::to_double(Coord(numerator) / scale)
        / CGAL::to_double(conjugate / scale);
}
} // namespace

void register_held_disk_contour_2(nb::module_& module)
{
    nb::exception<InvalidHeldContourInputError>(module, "InvalidHeldContourInputError", PyExc_ValueError);
    nb::exception<UndefinedPredecessorDirectionError>(module, "UndefinedPredecessorDirectionError", PyExc_ValueError);
    nb::exception<NoExposedPredecessorArcError>(module, "NoExposedPredecessorArcError", PyExc_RuntimeError);
    nb::exception<BrokenHeldContourBoundaryError>(module, "BrokenHeldContourBoundaryError", PyExc_RuntimeError);
    nb::class_<HeldDiskContour2>(module, "HeldDiskContour2")
        .def(nb::init<const HeldDiskContour2::XY&, double, double>(),
             "center"_a, "guide_radius"_a, "tool_radius"_a)
        .def("append", &HeldDiskContour2::append, "center"_a, "guide_radius"_a)
        .def("engagement_bound", [](const HeldDiskContour2& contour,
                                     const HeldDiskContour2::XY& center,
                                     double guide_radius, double cap_chord_ratio) {
            const auto [cosine, exceeded] = contour.engagement_bound(center, guide_radius, cap_chord_ratio);
            // The native cosine lies exactly in [-1,1]. Clamp only its display
            // conversion for acos; the exact cap decision above never uses it.
            const double angle = std::acos(std::clamp(report_coordinate(cosine), -1.0, 1.0));
            return std::pair{angle, exceeded};
        }, "candidate_center"_a, "guide_radius"_a, "cap_chord_ratio"_a)
        .def("contact_toward", [](const HeldDiskContour2& contour,
                                   const HeldDiskContour2::XY& center) {
            const auto [point, moved] = contour.contact_toward(center);
            // Display/reporting boundary only: no rounded point re-enters CGAL.
            return std::pair{std::pair{report_coordinate(point.x()),
                                       report_coordinate(point.y())}, moved};
        }, "candidate_center"_a);
}
