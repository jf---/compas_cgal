#pragma once

// AABB tree over the declared primitives' outward-rounded boxes: pruning only.
//
// A primitive is returned unless its box is provably more than r from the
// element's box along some axis. Boxes are rounded outward from the exact
// data (an arc's box is its full circle's), the query box is the element box
// grown by r with outward rounding (Interval_nt), and box overlap is
// CGAL::do_overlap on doubles (comparisons only, exact), so a primitive is
// never dropped while it can be within r. Every returned pair is decided by
// the filter or the exact stage.
//
// The tree's geometric traits are Epick with two functors replaced: a datum
// is its own Bbox_2, and Bbox_2-Bbox_2 intersection is do_overlap (the
// kernel's Do_intersect_2 has no Bbox_2-Bbox_2 overload).

#include <CGAL/AABB_traits_2.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <cstddef>
#include <vector>

namespace compas_cgal::held_declared_clearance {

/// An axis-aligned box with binary64 ends (rounded outward from exact bounds).
struct Box {
    double xmin;
    double ymin;
    double xmax;
    double ymax;
};

/// Epick, with boxes as data and box overlap as the intersection test.
struct BoxTraits : CGAL::Exact_predicates_inexact_constructions_kernel {
    struct Do_intersect_2 {
        using result_type = bool;
        bool operator()(const CGAL::Bbox_2& first, const CGAL::Bbox_2& second) const
        {
            return CGAL::do_overlap(first, second);
        }
    };
    struct Construct_bbox_2 {
        using result_type = CGAL::Bbox_2;
        CGAL::Bbox_2 operator()(const CGAL::Bbox_2& box) const { return box; }
    };
    [[nodiscard]] Do_intersect_2 do_intersect_2_object() const { return {}; }
    [[nodiscard]] Construct_bbox_2 construct_bbox_2_object() const { return {}; }
};

class PrimitiveIndex {
public:
    using BoxKernel = CGAL::Exact_predicates_inexact_constructions_kernel;

    /// One indexed box.
    struct Entry {
        Box box;
        std::size_t index;
    };

    /// AABB primitive over an Entry (Datum: its Bbox_2, Id: primitive index).
    class BoxPrimitive {
    public:
        using Id = std::size_t;
        using Datum = CGAL::Bbox_2;
        using Point = BoxKernel::Point_2;

        BoxPrimitive() = default;
        explicit BoxPrimitive(std::vector<Entry>::const_iterator entry);

        [[nodiscard]] Datum datum() const;
        [[nodiscard]] Id id() const { return index_; }
        [[nodiscard]] Point reference_point() const;

    private:
        Box box_{0, 0, 0, 0};
        std::size_t index_ = 0;
    };

    using Traits = CGAL::AABB_traits_2<BoxTraits, BoxPrimitive>;
    using Tree = CGAL::AABB_tree<Traits>;

    explicit PrimitiveIndex(std::vector<Box> boxes);
    PrimitiveIndex(const PrimitiveIndex&) = delete;
    PrimitiveIndex& operator=(const PrimitiveIndex&) = delete;

    /// Indices of the primitives not provably farther than `reach` from `box`, ascending.
    void near(const Box& box, double reach, std::vector<std::size_t>& out) const;

private:
    std::vector<Entry> entries_;
    Tree tree_;
};

}  // namespace compas_cgal::held_declared_clearance
