#include "annulus_view_2.h"
#include "engagement_2.h"
#include "stock_2.h"

#include <cstddef>
#include <tuple>
#include <vector>

#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace {

using SampleTuple = std::tuple<double, double, bool>;

SampleTuple sample_tuple(const EngagementSample& sample)
{
    return {sample.total_tea, sample.max_run_tea, sample.cap_exceeded};
}

std::vector<SampleTuple> sample_tuples(const std::vector<EngagementSample>& samples)
{
    std::vector<SampleTuple> tuples;
    tuples.reserve(samples.size());
    for (const EngagementSample& sample : samples) {
        tuples.push_back(sample_tuple(sample));
    }
    return tuples;
}

} // namespace

void register_annulus_view_2(nb::module_& m)
{
    nb::exception<InvalidAnnulusViewPolygonError>(
        m, "InvalidAnnulusViewPolygonError", PyExc_ValueError);
    nb::exception<NonFiniteAnnulusViewInputError>(
        m, "NonFiniteAnnulusViewInputError", PyExc_ValueError);
    nb::exception<InvalidAnnulusViewRadiiError>(
        m, "InvalidAnnulusViewRadiiError", PyExc_ValueError);
    nb::exception<InvalidAnnulusViewRangeError>(
        m, "InvalidAnnulusViewRangeError", PyExc_ValueError);
    nb::exception<AnnulusViewWitnessError>(
        m, "AnnulusViewWitnessError", PyExc_RuntimeError);

    nb::class_<AnnulusView2>(m, "AnnulusView2")
        .def_static("build", &AnnulusView2::build, "polygon"_a, "annuli"_a, "tool_radius"_a)
        .def_prop_ro("size", &AnnulusView2::size)
        .def("station",
             [](const AnnulusView2& view, std::size_t first, std::size_t end,
                double cx, double cy, double cap_chord_ratio, double gap_close_ratio) {
                 return sample_tuple(view.station(
                     {first, end}, cx, cy, cap_chord_ratio, gap_close_ratio));
             },
             "first"_a, "end"_a, "cx"_a, "cy"_a, "cap_chord_ratio"_a,
             "gap_close_ratio"_a = 0.0)
        .def("circle_prefixes",
             [](const AnnulusView2& view, std::size_t first, std::size_t end,
                double guide_cx, double guide_cy, double guide_radius,
                double entry_x, double entry_y, bool clockwise,
                const std::vector<double>& parameters, double cap_chord_ratio,
                double gap_close_ratio) {
                 return sample_tuples(view.circle_prefixes(
                     {first, end}, guide_cx, guide_cy, guide_radius, entry_x, entry_y,
                     clockwise, parameters, cap_chord_ratio, gap_close_ratio));
             },
             "first"_a, "end"_a, "guide_cx"_a, "guide_cy"_a, "guide_radius"_a,
             "entry_x"_a, "entry_y"_a, "clockwise"_a, "parameters"_a,
             "cap_chord_ratio"_a, "gap_close_ratio"_a = 0.0)
        .def("circle_exceeded",
             [](const AnnulusView2& view, std::size_t first, std::size_t end,
                double guide_cx, double guide_cy, double guide_radius,
                double entry_x, double entry_y, bool clockwise,
                const std::vector<double>& parameters, double cap_chord_ratio,
                double gap_close_ratio) {
                 return view.circle_exceeded(
                     {first, end}, guide_cx, guide_cy, guide_radius, entry_x, entry_y,
                     clockwise, parameters, cap_chord_ratio, gap_close_ratio);
             },
             "first"_a, "end"_a, "guide_cx"_a, "guide_cy"_a, "guide_radius"_a,
             "entry_x"_a, "entry_y"_a, "clockwise"_a, "parameters"_a,
             "cap_chord_ratio"_a, "gap_close_ratio"_a = 0.0)
        .def("station_matches_stock",
             [](const AnnulusView2& view, const Stock2& stock, std::size_t first,
                std::size_t end, double cx, double cy, double cap_chord_ratio,
                double gap_close_ratio) {
                 const StationAgreement2 agreement = view.station_agrees_with_stock(
                     stock, {first, end}, cx, cy, cap_chord_ratio, gap_close_ratio);
                 return agreement.disposition && agreement.runs;
             },
             "stock"_a, "first"_a, "end"_a, "cx"_a, "cy"_a, "cap_chord_ratio"_a,
             "gap_close_ratio"_a = 0.0)
        .def("circle_matches_stock",
             [](const AnnulusView2& view, const Stock2& stock, std::size_t first,
                std::size_t end, double guide_cx, double guide_cy,
                double guide_radius, double entry_x, double entry_y, bool clockwise,
                const std::vector<double>& parameters, double cap_chord_ratio,
                double gap_close_ratio) {
                 std::vector<std::tuple<bool, bool>> rows;
                 for (const StationAgreement2& agreement : view.circle_agrees_with_stock(
                          stock, {first, end}, guide_cx, guide_cy, guide_radius,
                          entry_x, entry_y, clockwise, parameters, cap_chord_ratio,
                          gap_close_ratio)) {
                     rows.emplace_back(agreement.disposition, agreement.runs);
                 }
                 return rows;
             },
             "stock"_a, "first"_a, "end"_a, "guide_cx"_a, "guide_cy"_a,
             "guide_radius"_a, "entry_x"_a, "entry_y"_a, "clockwise"_a, "parameters"_a,
             "cap_chord_ratio"_a, "gap_close_ratio"_a = 0.0);
}
