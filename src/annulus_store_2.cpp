#include "annulus_store_2.h"

#include <CGAL/FPU.h>

#include <algorithm>
#include <cmath>

namespace {

using Interval = CGAL::Interval_nt<false>;

// Grid cell edge in tool radii: on Figure 5 and Monstera a rim's box then reads about
// four cells. Measured, not derived; it only moves cost.
constexpr double CELL_TOOL_RADII = 4.0;
// Upper bound on cells per axis, so a tiny tool on a large part cannot size the grid
// past 65,536 cell headers.
constexpr double MAX_CELLS_PER_AXIS = 256.0;

std::size_t clamped_cell(double coordinate, double origin, double cell, std::size_t count)
{
    const double index = std::floor((coordinate - origin) / cell);
    if (!(index > 0.0)) return 0;
    const double last = static_cast<double>(count - 1);
    return index >= last ? count - 1 : static_cast<std::size_t>(index);
}

} // namespace

AnnulusStore2::AnnulusStore2(double tool_radius, const AnnulusSearchBox2& domain)
{
    const CGAL::Protect_FPU_rounding<true> rounding;
    tool_radius_ = Interval(tool_radius);
    origin_x_ = domain.min_x;
    origin_y_ = domain.min_y;
    const double width = domain.max_x - domain.min_x;
    const double height = domain.max_y - domain.min_y;
    cell_ = std::max({
        CELL_TOOL_RADII * tool_radius,
        width / MAX_CELLS_PER_AXIS,
        height / MAX_CELLS_PER_AXIS,
    });
    columns_ = static_cast<std::size_t>(std::ceil(width / cell_)) + 1;
    rows_ = static_cast<std::size_t>(std::ceil(height / cell_)) + 1;
    cells_.resize(columns_ * rows_);
}

AnnulusStore2::CellRange AnnulusStore2::cells_of(const AnnulusSearchBox2& box) const
{
    const CGAL::Protect_FPU_rounding<true> rounding;
    return {
        clamped_cell(box.min_x, origin_x_, cell_, columns_),
        clamped_cell(box.max_x, origin_x_, cell_, columns_),
        clamped_cell(box.min_y, origin_y_, cell_, rows_),
        clamped_cell(box.max_y, origin_y_, cell_, rows_),
    };
}

void AnnulusStore2::insert_references(std::size_t id)
{
    const CellRange range = cells_of(annuli_[id].box);
    for (std::size_t y = range.y0; y <= range.y1; ++y) {
        for (std::size_t x = range.x0; x <= range.x1; ++x) {
            cells_[y * columns_ + x].push_back(static_cast<std::uint32_t>(id));
            ++live_references_;
        }
    }
}

std::size_t AnnulusStore2::append(double cx, double cy, double guide_radius)
{
    if (!std::isfinite(cx) || !std::isfinite(cy) || !std::isfinite(guide_radius)) {
        throw NonFiniteAnnulusViewInputError(
            "annulus view centres and guide radii must be finite");
    }
    if (guide_radius < 0.0) {
        throw InvalidAnnulusViewRadiiError(
            "annulus view guide radii must be nonnegative");
    }
    StoredAnnulus2 stored;
    {
        const CGAL::Protect_FPU_rounding<true> rounding;
        const Interval guide(guide_radius);
        stored.cx = cx;
        stored.cy = cy;
        stored.guide_radius = guide_radius;
        stored.has_inner = guide_radius > tool_radius_.inf();
        stored.outer = guide + tool_radius_;
        stored.inner = stored.has_inner ? Interval(guide - tool_radius_) : Interval(0.0);
        stored.outer_squared = CGAL::square(stored.outer);
        stored.inner_squared = CGAL::square(stored.inner);
        const Interval x(cx);
        const Interval y(cy);
        const Interval min_x = x - stored.outer;
        const Interval min_y = y - stored.outer;
        const Interval max_x = x + stored.outer;
        const Interval max_y = y + stored.outer;
        stored.box = {min_x.inf(), min_y.inf(), max_x.sup(), max_y.sup()};
    }
    annuli_.push_back(stored);
    stamps_.resize(annuli_.size(), 0);
    const std::size_t id = annuli_.size() - 1;
    insert_references(id);
    return id;
}

void AnnulusStore2::truncate(std::size_t size)
{
    if (size > annuli_.size()) {
        throw InvalidAnnulusViewRangeError("annulus store truncation beyond its size");
    }
    for (std::size_t id = size; id < annuli_.size(); ++id) {
        const CellRange range = cells_of(annuli_[id].box);
        const std::size_t references = (range.x1 - range.x0 + 1) * (range.y1 - range.y0 + 1);
        live_references_ -= references;
        stale_references_ += references;
    }
    annuli_.resize(size);
    if (stale_references_ > live_references_) rebuild();
}

void AnnulusStore2::rebuild()
{
    for (auto& cell : cells_) cell.clear();
    live_references_ = 0;
    stale_references_ = 0;
    for (std::size_t id = 0; id < annuli_.size(); ++id) insert_references(id);
    ++rebuilds_;
}

void AnnulusStore2::query(
    const AnnulusSearchBox2& box, AnnulusViewRange2 range,
    std::vector<std::uint32_t>& out) const
{
    ++generation_;
    if (generation_ == 0) {
        std::fill(stamps_.begin(), stamps_.end(), 0);
        generation_ = 1;
    }
    const CellRange cells = cells_of(box);
    for (std::size_t y = cells.y0; y <= cells.y1; ++y) {
        for (std::size_t x = cells.x0; x <= cells.x1; ++x) {
            for (const std::uint32_t id : cells_[y * columns_ + x]) {
                if (id < range.first || id >= range.end) continue;
                if (stamps_[id] == generation_) continue;
                stamps_[id] = generation_;
                if (search_boxes_meet(annuli_[id].box, box)) out.push_back(id);
            }
        }
    }
}
