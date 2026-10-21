#include "held_declared_clearance_2/primitive_index_2.h"

#include <CGAL/Interval_nt.h>

#include <algorithm>
#include <iterator>
#include <utility>

namespace compas_cgal::held_declared_clearance {

PrimitiveIndex::BoxPrimitive::BoxPrimitive(std::vector<Entry>::const_iterator entry)
    : box_(entry->box), index_(entry->index)
{
}

PrimitiveIndex::BoxPrimitive::Datum PrimitiveIndex::BoxPrimitive::datum() const
{
    return CGAL::Bbox_2(box_.xmin, box_.ymin, box_.xmax, box_.ymax);
}

PrimitiveIndex::BoxPrimitive::Point PrimitiveIndex::BoxPrimitive::reference_point() const
{
    return Point(box_.xmin, box_.ymin);
}

PrimitiveIndex::PrimitiveIndex(std::vector<Box> boxes)
{
    entries_.reserve(boxes.size());
    for (std::size_t index = 0; index < boxes.size(); ++index) entries_.push_back({boxes[index], index});
    tree_.insert(entries_.cbegin(), entries_.cend());
    tree_.build();
}

void PrimitiveIndex::near(const Box& box, double reach, std::vector<std::size_t>& out) const
{
    // Outward-rounded growth by `reach` (self-protected intervals).
    using Grow = CGAL::Interval_nt<true>;
    const CGAL::Bbox_2 query((Grow(box.xmin) - reach).inf(), (Grow(box.ymin) - reach).inf(),
                             (Grow(box.xmax) + reach).sup(), (Grow(box.ymax) + reach).sup());
    out.clear();
    tree_.all_intersected_primitives(query, std::back_inserter(out));
    std::sort(out.begin(), out.end());
}

}  // namespace compas_cgal::held_declared_clearance
