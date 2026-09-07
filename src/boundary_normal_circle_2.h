#pragma once

#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Polygon_2.h>

#include <array>
#include <cstdint>
#include <stdexcept>
#include <vector>

class NativeBoundary2;

namespace boundary_normal {

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt;
using Point = Kernel::Point_2;
using FT = Kernel::FT;
using XY = std::array<double, 2>;

class InvalidBoundaryPolygonError : public std::runtime_error {
public: using std::runtime_error::runtime_error;
};
class InvalidBoundaryNormalInputError : public std::runtime_error {
public: using std::runtime_error::runtime_error;
};
class BoundaryVertexQueryUnsupportedError : public std::runtime_error {
public: using std::runtime_error::runtime_error;
};
class InvalidBoundaryVertexDirectionError : public std::runtime_error {
public: using std::runtime_error::runtime_error;
};
class NoPositiveBoundaryCircleError : public std::runtime_error {
public: using std::runtime_error::runtime_error;
};
class InvalidBoundaryCircleContactError : public std::runtime_error {
public: using std::runtime_error::runtime_error;
};
class BoundaryNormalConstructionError : public std::runtime_error {
public: using std::runtime_error::runtime_error;
};

class BoundaryNormalCircle2;

// Exact world-XY geometry in mm. Only reporting properties cross into doubles.
class BoundaryNormalCircleProposal2 {
public:
    // Native consumer seam: retain the exact cutter-centre boundary contact.
    const Point& exact_contact() const noexcept { return q_; }
    const Point& exact_center() const noexcept { return center_; }
    const FT& exact_guide_radius() const noexcept { return guide_radius_; }
    FT exact_tool_radius() const { return clearance_ - FT(2) * guide_radius_; }
    bool is_stationary() const { return CGAL::is_zero(guide_radius_); }
    XY p_mm() const;
    XY m_mm() const;
    XY q_mm() const;
    XY center_mm() const;
    double guide_radius_mm() const;
    double clearance_mm() const;
    const std::vector<std::size_t>& competing_vertex_indices() const;
    const std::vector<std::size_t>& competing_segment_indices() const;
    const std::vector<std::size_t>& competing_arc_indices() const { return arcs_; }

private:
    BoundaryNormalCircleProposal2(
        Point p, Point m, Point q, Point center, FT guide_radius, FT clearance,
        std::vector<std::size_t> vertices, std::vector<std::size_t> segments,
        std::vector<std::size_t> arcs = {});
    Point p_, m_, q_, center_;
    FT guide_radius_, clearance_;
    std::vector<std::size_t> vertices_, segments_, arcs_;
    friend class BoundaryNormalCircle2;
    friend class ::NativeBoundary2;
};

// One validated simple polygon. Segment indices retain the supplied order;
// either winding is accepted. Holes are excluded. Reflex vertices have an
// explicit point-site normal-cone query, separate from segment parameters.
class BoundaryNormalCircle2 {
public:
    explicit BoundaryNormalCircle2(const std::vector<XY>& boundary);
    BoundaryNormalCircleProposal2 query(
        std::int64_t source_segment, double parameter, double tool_radius) const;
    BoundaryNormalCircleProposal2 query_vertex(
        std::int64_t vertex, const XY& inward_direction, double tool_radius) const;

    BoundaryNormalCircleProposal2 at_contact(const Point& q, double tool_radius) const;

private:
    BoundaryNormalCircleProposal2 construct(
        const Point& p, const Kernel::Vector_2& normal, const FT& tool, const FT& normal_length, bool allow_stationary = false) const;
    CGAL::Polygon_2<Kernel> polygon_;
};

} // namespace boundary_normal
