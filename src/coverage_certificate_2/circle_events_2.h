#pragma once

// The events on one boundary circle, both kinds behind one interface: crossings
// with other circles (Circular_kernel_2 points, lazy_crossings_2.h) and crossings
// with line supports (tower points, line_crossings_2.h).
//
// The circle sweep orders events CCW from the leftmost point: half (lower,
// upper), then x. Every comparison runs the certified interval front first and
// falls back to an exact predicate:
//   * circle-circle vs circle-circle: CK compare_x (cross-root, exact);
//   * any pair with a rational-base side (CK points are rational-base towers):
//     the tower sign rule, sign(t1 - t2) (tower_2.h);
//   * two line events over two different irrational bases (two flanks of
//     non-parallel edges crossing the circle near one point): equality is
//     decided by provenance -- X1 lies on line 2 and on line 2's branch --
//     and distinct points are separated by exact rational bisection on x (they
//     are distinct, and x is strictly monotone on one half, so it terminates).

#include "coverage_certificate_2/counters_2.h"
#include "coverage_certificate_2/lazy_crossings_2.h"
#include "coverage_certificate_2/line_crossings_2.h"
#include "coverage_certificate_2/lines_2.h"
#include "coverage_certificate_2/regions_2.h"
#include "coverage_certificate_2/vertex_events_2.h"

#include <cstdint>
#include <utility>
#include <vector>

namespace compas_cgal::coverage_certificate {

class CircleEvents {
public:
    CircleEvents(const CircleRegionSet& set, const LineTable& lines, CrossingTable& table,
                 const std::vector<VertexCircle>* vertices = nullptr)
        : set_(&set), lines_(&lines), table_(&table), vertices_(vertices)
    {
    }

    /// Add the crossings of circle k with line `line` to table.events[k] (source 1).
    void add_line(std::uint32_t k, std::uint32_t line, DecisionCounters& counters);

    /// Add the crossings of circle k with vertex circle `vertex` (source 2; cut
    /// supports, they flip nothing).
    void add_vertex(std::uint32_t k, std::uint32_t vertex, DecisionCounters& counters);

    /// 0 lower half (y < c_y, or y == c_y and x < c_x), 1 upper half.
    [[nodiscard]] int half(const CrossingEvent& event, const BoundaryCircle& circle, DecisionCounters& counters);

    /// Exact comparison of the x coordinates of two events on one circle.
    [[nodiscard]] CGAL::Comparison_result compare_x(const CrossingEvent& a, const CrossingEvent& b,
                                                    DecisionCounters& counters);

    /// Exact comparison of an event's x with the x of a rational point.
    [[nodiscard]] CGAL::Comparison_result compare_x(const CrossingEvent& a, const Point_2& q, const Interval& qx,
                                                    DecisionCounters& counters);

    /// The event point's coordinates as doubles from its EXACT value (mode-
    /// independent; witness proposals only, never a decision).
    [[nodiscard]] std::pair<double, double> proposal_xy(const CrossingEvent& event, DecisionCounters& counters);

    /// The event point's exact coordinates as towers (materialised on demand).
    [[nodiscard]] std::pair<Tower, Tower> exact_xy(const CrossingEvent& event, DecisionCounters& counters);

    [[nodiscard]] FilterMode mode() const { return table_->mode(); }
    [[nodiscard]] CrossingTable& table() { return *table_; }
    [[nodiscard]] const VertexCircle& vertex(std::uint32_t id) const { return (*vertices_)[id]; }

private:
    [[nodiscard]] const CircleLineExact& exact_line(const CrossingEvent& event, DecisionCounters& counters);
    [[nodiscard]] const CircleVertexExact& exact_vertex(const CrossingEvent& event, DecisionCounters& counters);
    /// Whether event a's point is event b's point (b of source 1 or 2), decided
    /// exactly by provenance: a lies on b's support and on b's branch.
    [[nodiscard]] bool same_point_by_provenance(const CrossingEvent& a, const CrossingEvent& b, DecisionCounters& counters);
    [[nodiscard]] Tower exact_x(const CrossingEvent& event, DecisionCounters& counters);
    [[nodiscard]] CGAL::Comparison_result separate_cross_base(const CrossingEvent& a, const CrossingEvent& b,
                                                              DecisionCounters& counters);

    const CircleRegionSet* set_;
    const LineTable* lines_;
    CrossingTable* table_;
    std::vector<CircleLineRecord> records_;
    const std::vector<VertexCircle>* vertices_;
    std::vector<CircleVertexRecord> vertex_records_;
};


}  // namespace compas_cgal::coverage_certificate
