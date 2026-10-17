#include "held_boundary_clearance_2.h"

#include <ios>
#include <string>
#include <tuple>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace {
using held_boundary_clearance::AdmissibleSite;
using held_boundary_clearance::Integer;

// Exact integers cross as signed hex text (int(text, 16) in Python); no double
// ever leaves this boundary for these values.
using SiteRecord = std::tuple<std::size_t, std::string, std::string, std::string, std::string, std::string, int>;

std::string integer_text(const Integer& value)
{
    const bool negative = value < 0;
    const Integer magnitude = negative ? Integer(-value) : value;
    std::string text = magnitude.str(0, std::ios_base::hex);
    return negative ? "-" + text : text;
}

SiteRecord site_record(const AdmissibleSite& site)
{
    return {site.segment, integer_text(site.along), integer_text(site.length_squared), integer_text(site.foot_x),
            integer_text(site.foot_y), integer_text(site.distance_numerator), site.exponent};
}
} // namespace

void register_held_boundary_clearance_2(nb::module_& module)
{
    nb::exception<InvalidHeldBoundaryClearanceInputError>(module, "InvalidHeldBoundaryClearanceInputError", PyExc_ValueError);
    nb::exception<DegenerateHeldBoundarySegmentError>(module, "DegenerateHeldBoundarySegmentError", PyExc_ValueError);
    nb::exception<BrokenHeldBoundaryClearancePruningError>(module, "BrokenHeldBoundaryClearancePruningError", PyExc_RuntimeError);
    nb::class_<HeldBoundaryClearance2>(module, "HeldBoundaryClearance2")
        .def(nb::init<const std::vector<HeldBoundaryClearance2::XY>&>(), "ring"_a)
        .def_prop_ro("segment_count", &HeldBoundaryClearance2::segment_count)
        .def("distance_admissible_sites", [](const HeldBoundaryClearance2& index,
                                              const std::vector<HeldBoundaryClearance2::XY>& stations,
                                              double gap) {
            std::vector<std::vector<SiteRecord>> records;
            {
                nb::gil_scoped_release release;
                const auto sites = index.distance_admissible_sites(stations, gap);
                records.reserve(sites.size());
                for (const auto& station : sites) {
                    auto& row = records.emplace_back();
                    row.reserve(station.size());
                    for (const auto& site : station) row.push_back(site_record(site));
                }
            }
            return records;
        }, "stations"_a, "gap"_a);
}
