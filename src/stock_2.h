#pragma once

#include "compas.h"

// compas.h supplies nanobind/stl/bind_vector.h (opaque nb::bind_vector) but NOT
// the automatic std::vector<T> <-> Python-list type caster. The Stock2
// constructor takes `holes` as a plain Python list, so this module needs the
// caster in its own (NB_STATIC) translation unit.
#include <nanobind/stl/vector.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/Gps_circle_segment_traits_2.h>
#include <CGAL/Polygon_2.h>

// Exact constructions kernel: boolean set operations on the 2D stock model
// build new geometry (sweep boundaries, intersection vertices) whose exact
// coordinates must be representable, so Epeck (not the repo-default Epick) is
// mandatory here. See CLAUDE.md: exact predicates, no epsilon decisions.
typedef CGAL::Exact_predicates_exact_constructions_kernel Epeck;
typedef CGAL::Gps_circle_segment_traits_2<Epeck> GpsTraits;
typedef CGAL::General_polygon_set_2<GpsTraits> Gps;
typedef GpsTraits::Polygon_2 GpsPolygon;             // circle-segment general polygon
typedef GpsTraits::Polygon_with_holes_2 GpsPolygonWithHoles;
typedef GpsTraits::Point_2 GpsPoint;                 // one-root coordinates
typedef GpsTraits::X_monotone_curve_2 GpsXCurve;
typedef Epeck::Point_2 EPoint;
typedef Epeck::Circle_2 ECircle;
typedef Epeck::Segment_2 ESegment;

// Exact 2D stock model: the remaining material as a general polygon set of
// linear (this task) and circular (later tasks) boundary arcs. Later tasks
// subtract tool sweeps and query engagement; this task owns init + point-in.
class Stock2 {
public:
    Stock2(Eigen::Ref<const compas::RowMatrixXd> boundary,
           const std::vector<compas::RowMatrixXd>& holes);

    bool contains(double x, double y) const;
    bool is_empty() const;

    // Remove the exact disk of the given radius centred at (cx, cy).
    void subtract_disk(double cx, double cy, double radius);

    // Remove the tool sweep along segment (x0,y0)->(x1,y1) as a certified
    // under-covering disk chain (the exact oriented capsule has irrational
    // side lines, so it is not representable in the circle-segment traits).
    void subtract_capsule(double x0, double y0, double x1, double y1, double radius);

    const Gps& set() const { return set_; }          // engagement kernel reads this
    Gps& set() { return set_; }

private:
    Gps set_;
};
