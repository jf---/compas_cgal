#pragma once

#include <CGAL/AABB_segment_primitive_2.h>
#include <CGAL/AABB_traits_2.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>

#include <boost/multiprecision/cpp_int.hpp>

#include <array>
#include <cstddef>
#include <stdexcept>
#include <vector>

class InvalidHeldBoundaryClearanceInputError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};
class DegenerateHeldBoundarySegmentError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};
// Internal contract breach: the pruning box must contain the hint segment.
class BrokenHeldBoundaryClearancePruningError : public std::logic_error {
public:
    using std::logic_error::logic_error;
};

namespace held_boundary_clearance {

// Every coordinate is an injected double, so Epick points and segments carry the
// exact input; Epick predicates and the Filtered_predicate wrappers in the .cpp
// decide exactly (interval filter, exact-rational fallback).
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = Kernel::Point_2;
using Segment = Kernel::Segment_2;
using Segments = std::vector<Segment>;
using FilterPrimitive = CGAL::AABB_segment_primitive_2<Kernel, Segments::const_iterator>;
using FilterTree = CGAL::AABB_tree<CGAL::AABB_traits_2<Kernel, FilterPrimitive>>;
// Nearest-segment hints only. CGAL 6.0.1's AABB_traits_2 distance query does not
// compile for statically filtered kernels (it calls a three-argument
// do_intersect(Circle_2, Bbox_2, bool) that exists only in 3D), so hints use a
// plain double kernel. A hint's quality sizes the pruning box; its distance is
// re-bounded with certified intervals, so hint rounding never reaches a decision.
using HintKernel = CGAL::Simple_cartesian<double>;
using HintSegments = std::vector<HintKernel::Segment_2>;
using HintPrimitive = CGAL::AABB_segment_primitive_2<HintKernel, HintSegments::const_iterator>;
using HintTree = CGAL::AABB_tree<CGAL::AABB_traits_2<HintKernel, HintPrimitive>>;

// Arbitrary-precision integer with 512 inline bits: typical mm inputs never
// allocate, extreme exponents (subnormals) grow without loss.
using Integer = boost::multiprecision::number<boost::multiprecision::cpp_int_backend<
    512, 0, boost::multiprecision::signed_magnitude, boost::multiprecision::unchecked,
    std::allocator<boost::multiprecision::limb_type>>>;

// One distance-admissible site on the closed ring, in exact integer closed form.
// With scale 2^exponent applied to every input coordinate:
//   parameter        = along / length_squared                       (in [0, 1))
//   footpoint        = (foot_x, foot_y) / length_squared * 2^exponent
//   squared_distance = distance_numerator / length_squared * 2^(2*exponent)
// `segment` is canonical: a foot at a segment's target is reported as parameter
// zero on the following segment, so every ring vertex has one identity.
struct AdmissibleSite {
    std::size_t segment;
    Integer along;
    Integer length_squared;
    Integer foot_x;
    Integer foot_y;
    Integer distance_numerator;
    int exponent;
};

} // namespace held_boundary_clearance

// Exact nearest-boundary query over a closed polygonal ring in world XY mm.
//
// For a station p with exact squared distances d_i to the ring's sites and the
// exact minimum m, the admissible set is {i : sqrt(d_i) - sqrt(m) <= gap},
// decided without square roots. Doubles enter once, by exact injection.
//
// The AABB trees only prune. An approximate nearest segment yields a certified
// interval upper bound U >= m; every admissible segment meets the disk of
// radius sqrt(U) + gap about p, which lies inside an outward-rounded box. The
// box-segment test is an exact Epick predicate, so pruning never drops an
// admissible segment.
class HeldBoundaryClearance2 {
public:
    using XY = std::array<double, 2>;

    explicit HeldBoundaryClearance2(const std::vector<XY>& ring);
    // The trees hold iterators into the owned segment vectors.
    HeldBoundaryClearance2(const HeldBoundaryClearance2&) = delete;
    HeldBoundaryClearance2& operator=(const HeldBoundaryClearance2&) = delete;

    std::vector<std::vector<held_boundary_clearance::AdmissibleSite>> distance_admissible_sites(
        const std::vector<XY>& stations, double gap) const;

    std::size_t segment_count() const { return segments_.size(); }

private:
    std::vector<held_boundary_clearance::AdmissibleSite> station_sites(const XY& station, double gap) const;

    held_boundary_clearance::Segments segments_;
    held_boundary_clearance::FilterTree tree_;
    held_boundary_clearance::HintSegments hint_segments_;
    held_boundary_clearance::HintTree hint_tree_;
};
