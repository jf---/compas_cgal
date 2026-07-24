#pragma once

#include "compas.h"
#include "exact_depletion_2.h"
#include "exact_motion_2.h"

// compas.h supplies nanobind/stl/bind_vector.h (opaque nb::bind_vector) but NOT
// the automatic std::vector<T> <-> Python-list type caster. The Stock2
// constructor takes `holes` as a plain Python list, so this module needs the
// caster in its own (NB_STATIC) translation unit.
#include <nanobind/stl/vector.h>

#include <memory>
#include <utility>

#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/Gps_circle_segment_traits_2.h>
#include <CGAL/Polygon_2.h>

// Exact constructions kernel: boolean set operations on the 2D stock model
// build new geometry (sweep boundaries, intersection vertices) whose exact
// coordinates must be representable, so Epeck (not the repo-default Epick) is
// mandatory here. See CLAUDE.md: exact predicates, no epsilon decisions.
typedef CGAL::Gps_circle_segment_traits_2<Epeck> GpsTraits;
typedef CGAL::General_polygon_set_2<GpsTraits> Gps;
typedef GpsTraits::Polygon_2 GpsPolygon;             // circle-segment general polygon
typedef GpsTraits::Polygon_with_holes_2 GpsPolygonWithHoles;
typedef GpsTraits::Point_2 GpsPoint;                 // one-root coordinates
typedef GpsTraits::X_monotone_curve_2 GpsXCurve;

// Full disk of the given radius centred at `center`, as a two-arc CCW general
// polygon (split at its x-extreme vertical-tangency points). Shared by the
// stock-subtraction paths and the engagement query (which intersects the stock
// with this exact cutter disk), so it carries external linkage here rather than
// hiding in stock_2.cpp's anonymous namespace.
GpsPolygon disk_polygon(const EPoint& center, const Epeck::FT& radius);

// Exact 2D stock model: the remaining material as a general polygon set of
// linear (this task) and circular (later tasks) boundary arcs. Later tasks
// subtract tool sweeps and query engagement; this task owns init + point-in.
class Stock2 {
public:
    Stock2(Eigen::Ref<const compas::RowMatrixXd> boundary,
           const std::vector<compas::RowMatrixXd>& holes);
    Stock2(Stock2&&) noexcept = default;
    Stock2& operator=(Stock2&&) noexcept = default;
    Stock2(const Stock2&) = delete;
    Stock2& operator=(const Stock2&) = delete;

    bool contains(double x, double y) const;
    bool is_empty() const;
    Stock2 clone() const;
    bool is_subset_of(const Stock2& other) const;
    bool exactly_equals(const Stock2& other) const;

    // Remove the exact disk of the given radius centred at (cx, cy).
    void subtract_disk(double cx, double cy, double radius);

    // Remove the tool sweep along segment (x0,y0)->(x1,y1) as a certified
    // under-covering disk chain (the exact oriented capsule has irrational
    // side lines, so it is not representable in the circle-segment traits).
    void subtract_capsule(double x0, double y0, double x1, double y1, double radius);

    // Remove the tool sweep along the circular guide arc from (sx,sy) to
    // (ex,ey) about (cx,cy) — cw selects the sweep direction, start == end
    // means the full circle — as a certified under-covering disk chain.
    void subtract_arc_sweep(double cx, double cy, double sx, double sy,
                            double ex, double ey, bool cw, double tool_radius);

    DepletionTrace subtract_exact_segment(
        const ExactSegmentMotion2& motion,
        const Epeck::FT& tool_radius,
        const Epeck::FT& max_chord,
        std::size_t center_count_limit);

    DepletionTrace subtract_exact_full_circle(
        const ExactCircleMotion2& motion,
        const Epeck::FT& tool_radius,
        const Epeck::FT& max_chord,
        std::size_t center_count_limit);

    const Gps& set() const { return *set_; }          // engagement kernel reads this
    Gps& set() { return *set_; }

private:
    explicit Stock2(std::unique_ptr<Gps> set) noexcept;

    // Subtract the union of exact tool disks of the given radius centred at the
    // listed points — the one chain implementation shared by capsule and arc.
    void subtract_point_chain(const std::vector<std::pair<double, double>>& centers,
                              double radius);

    std::unique_ptr<Gps> set_;
};
