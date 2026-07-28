#include "containment_2.h"

#include <nanobind/nanobind.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace {

nb::bytes bytes_value(const std::string &value) {
  return nb::bytes(value.data(), value.size());
}

std::string from_bytes(const nb::bytes &value) {
  char *data = nullptr;
  Py_ssize_t size = 0;
  if (PyBytes_AsStringAndSize(value.ptr(), &data, &size) != 0) {
    throw nb::python_error();
  }
  return std::string(data, static_cast<std::size_t>(size));
}

} // namespace

NB_MODULE(_containment_2, m) {
  nb::exception<ContainmentConstructionError>(m, "ContainmentConstructionError",
                                              PyExc_RuntimeError);

  nb::class_<ContainmentRecord2>(m, "ContainmentRecord2")
      .def_ro("contained", &ContainmentRecord2::contained)
      .def_ro("guide_anchor_in_center_domain",
              &ContainmentRecord2::guide_anchor_in_center_domain)
      .def_ro("outer_disk_contained", &ContainmentRecord2::outer_disk_contained)
      .def_ro("disk_sweep", &ContainmentRecord2::disk_sweep)
      .def_prop_ro("strategy_version",
                   [](const ContainmentRecord2 &record) {
                     return bytes_value(record.strategy_version);
                   })
      .def_prop_ro("structural_record",
                   [](const ContainmentRecord2 &record) {
                     return bytes_value(record.structural_record);
                   })
      .def(
          "matches_exact_segment",
          [](const ContainmentRecord2 &record,
             const nb::bytes &authority_digest, const double x0,
             const double y0, const double x1, const double y1,
             const double tool_radius) {
            return record.matches_exact_segment(from_bytes(authority_digest),
                                                x0, y0, x1, y1, tool_radius);
          },
          "authority_digest"_a, "x0"_a, "y0"_a, "x1"_a, "y1"_a, "tool_radius"_a)
      .def(
          "matches_exact_full_circle",
          [](const ContainmentRecord2 &record,
             const nb::bytes &authority_digest, const double center_x,
             const double center_y, const double phase_x, const double phase_y,
             const double tool_radius) {
            return record.matches_exact_full_circle(
                from_bytes(authority_digest), center_x, center_y, phase_x,
                phase_y, tool_radius);
          },
          "authority_digest"_a, "center_x"_a, "center_y"_a, "phase_x"_a,
          "phase_y"_a, "tool_radius"_a)
      .def(
          "matches_exact_entry_disk",
          [](const ContainmentRecord2 &record,
             const nb::bytes &authority_digest, const double center_x,
             const double center_y, const double entry_radius,
             const double tool_radius) {
            return record.matches_exact_entry_disk(from_bytes(authority_digest),
                                                   center_x, center_y,
                                                   entry_radius, tool_radius);
          },
          "authority_digest"_a, "center_x"_a, "center_y"_a, "entry_radius"_a,
          "tool_radius"_a)
      .def(
          "matches_exact_full_circle_in_disk",
          [](const ContainmentRecord2 &record,
             const nb::bytes &authority_digest, const double entry_center_x,
             const double entry_center_y, const double entry_radius,
             const double circle_center_x, const double circle_center_y,
             const double phase_x, const double phase_y,
             const double tool_radius) {
            return record.matches_exact_full_circle_in_disk(
                from_bytes(authority_digest), entry_center_x, entry_center_y,
                entry_radius, circle_center_x, circle_center_y, phase_x,
                phase_y, tool_radius);
          },
          "authority_digest"_a, "entry_center_x"_a, "entry_center_y"_a,
          "entry_radius"_a, "circle_center_x"_a, "circle_center_y"_a,
          "phase_x"_a, "phase_y"_a, "tool_radius"_a);

  m.def(
      "evaluate_exact_segment_containment",
      [](const ExactRegion2 &design, const ExactRegion2 &center_domain,
         const nb::bytes &authority_digest, const double x0, const double y0,
         const double x1, const double y1, const double tool_radius) {
        return evaluate_exact_segment_containment(design, center_domain,
                                                  from_bytes(authority_digest),
                                                  x0, y0, x1, y1, tool_radius);
      },
      "design"_a, "center_domain"_a, "authority_digest"_a, "x0"_a, "y0"_a,
      "x1"_a, "y1"_a, "tool_radius"_a);
  m.def(
      "evaluate_exact_full_circle_containment",
      [](const ExactRegion2 &design, const ExactRegion2 &center_domain,
         const nb::bytes &authority_digest, const double center_x,
         const double center_y, const double phase_x, const double phase_y,
         const double tool_radius) {
        return evaluate_exact_full_circle_containment(
            design, center_domain, from_bytes(authority_digest), center_x,
            center_y, phase_x, phase_y, tool_radius);
      },
      "design"_a, "center_domain"_a, "authority_digest"_a, "center_x"_a,
      "center_y"_a, "phase_x"_a, "phase_y"_a, "tool_radius"_a);
  m.def(
      "evaluate_exact_entry_disk_containment",
      [](const ExactRegion2 &design, const ExactRegion2 &center_domain,
         const nb::bytes &authority_digest, const double center_x,
         const double center_y, const double entry_radius,
         const double tool_radius) {
        return evaluate_exact_entry_disk_containment(
            design, center_domain, from_bytes(authority_digest), center_x,
            center_y, entry_radius, tool_radius);
      },
      "design"_a, "center_domain"_a, "authority_digest"_a, "center_x"_a,
      "center_y"_a, "entry_radius"_a, "tool_radius"_a);
  m.def(
      "evaluate_exact_full_circle_in_disk",
      [](const nb::bytes &authority_digest, const double entry_center_x,
         const double entry_center_y, const double entry_radius,
         const double circle_center_x, const double circle_center_y,
         const double phase_x, const double phase_y, const double tool_radius) {
        return evaluate_exact_full_circle_in_disk(
            from_bytes(authority_digest), entry_center_x, entry_center_y,
            entry_radius, circle_center_x, circle_center_y, phase_x, phase_y,
            tool_radius);
      },
      "authority_digest"_a, "entry_center_x"_a, "entry_center_y"_a,
      "entry_radius"_a, "circle_center_x"_a, "circle_center_y"_a, "phase_x"_a,
      "phase_y"_a, "tool_radius"_a);
  m.def("containment_strategy_version",
        []() { return bytes_value(::containment_strategy_version()); });
}
