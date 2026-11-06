#pragma once

// The append-only store of emitted swept annuli a refinement decision reads: per
// circle its injected (c, rho), outward-rounded enclosures of its exact radii, and
// its outward-rounded box in a uniform grid.
//
// A view is an emission range [first, end). Truncation drops the tail in O(1) and
// leaves its references in the grid: every query filters by range and box, so a
// stale reference costs one test and never selects anything. The grid is rebuilt
// once stale references outnumber live ones.
//
// Selection is never a decision: a box miss proves an exact miss, a hit only makes a
// candidate. Cell indices are always computed under upward rounding, so the same
// coordinate falls in the same cell at insertion and at query, and the index
// function is monotone: two boxes that meet share a cell.
//
// Why a grid, not Point_set_2: a Point_set_2 query radius has to reach the largest
// annulus, and the set merges coincident centres, which would lose a concentric annulus.

#include <CGAL/Interval_nt.h>

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

class NonFiniteAnnulusViewInputError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

class InvalidAnnulusViewRadiiError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

class InvalidAnnulusViewRangeError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

// The half-open emission range [first, end) of annuli a decision reads.
struct AnnulusViewRange2 {
    std::size_t first;
    std::size_t end;
};

// The nested views of one decision: the last 1 and the last 16 emissions of its range,
// then all of it. A tier's size is cost only: a pass on any of them is a pass on the
// range, and only a refusal on the range refuses.
inline constexpr std::size_t DECISION_TIER_LAST_ONE = 1;
inline constexpr std::size_t DECISION_TIER_LAST_SIXTEEN = 16;

inline AnnulusViewRange2 last_of_range(AnnulusViewRange2 range, std::size_t last)
{
    return {range.end - range.first > last ? range.end - last : range.first, range.end};
}

// An axis-aligned box in doubles. Every coordinate is an injected double or an
// outward-rounded bound.
struct AnnulusSearchBox2 {
    double min_x;
    double min_y;
    double max_x;
    double max_y;
};

inline bool search_boxes_meet(const AnnulusSearchBox2& lhs, const AnnulusSearchBox2& rhs)
{
    return !(lhs.max_x < rhs.min_x || lhs.min_x > rhs.max_x
             || lhs.max_y < rhs.min_y || lhs.min_y > rhs.max_y);
}

// One emitted circle's closed annulus {max(rho - r, 0) <= |x - c| <= rho + r}.
struct StoredAnnulus2 {
    double cx;
    double cy;
    double guide_radius;
    bool has_inner;  // rho > r: compared on the doubles, which are the exact values
    CGAL::Interval_nt<false> inner;
    CGAL::Interval_nt<false> outer;
    CGAL::Interval_nt<false> inner_squared;
    CGAL::Interval_nt<false> outer_squared;
    AnnulusSearchBox2 box;
};

class AnnulusStore2 {
public:
    // `domain` only sizes the grid: boxes outside it clamp to the border cells.
    AnnulusStore2(double tool_radius, const AnnulusSearchBox2& domain);

    std::size_t size() const noexcept { return annuli_.size(); }
    const StoredAnnulus2& operator[](std::size_t id) const { return annuli_[id]; }

    // Raises:
    //   NonFiniteAnnulusViewInputError: a non-finite centre or guide radius.
    //   InvalidAnnulusViewRadiiError: a negative guide radius.
    std::size_t append(double cx, double cy, double guide_radius);

    // Raises InvalidAnnulusViewRangeError: `size` exceeds the store.
    void truncate(std::size_t size);

    // Every annulus of `range` whose box meets `box`, each once, in no fixed order.
    void query(
        const AnnulusSearchBox2& box, AnnulusViewRange2 range,
        std::vector<std::uint32_t>& out) const;

    std::size_t stale_references() const noexcept { return stale_references_; }
    std::size_t rebuilds() const noexcept { return rebuilds_; }

private:
    struct CellRange {
        std::size_t x0;
        std::size_t x1;
        std::size_t y0;
        std::size_t y1;
    };
    CellRange cells_of(const AnnulusSearchBox2& box) const;
    void insert_references(std::size_t id);
    void rebuild();

    CGAL::Interval_nt<false> tool_radius_;
    double origin_x_;
    double origin_y_;
    double cell_;
    std::size_t columns_;
    std::size_t rows_;
    std::vector<StoredAnnulus2> annuli_;
    std::vector<std::vector<std::uint32_t>> cells_;
    std::size_t live_references_ = 0;
    std::size_t stale_references_ = 0;
    std::size_t rebuilds_ = 0;
    mutable std::vector<std::uint32_t> stamps_;
    mutable std::uint32_t generation_ = 0;
};
