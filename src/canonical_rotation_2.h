#pragma once

#include <algorithm>
#include <cstddef>
#include <vector>

#include <CGAL/enum.h>

template <typename Value, typename Compare>
std::size_t minimal_rotation_index(
    const std::vector<Value>& values,
    Compare compare,
    std::size_t& comparison_count)
{
    const std::size_t size = values.size();
    std::size_t left = 0;
    std::size_t right = 1;
    std::size_t offset = 0;
    while (left < size && right < size && offset < size) {
        ++comparison_count;
        const CGAL::Comparison_result order = compare(
            values[(left + offset) % size],
            values[(right + offset) % size]);
        if (order == CGAL::EQUAL) {
            ++offset;
            continue;
        }
        if (order == CGAL::LARGER) {
            left += offset + 1;
            if (left == right) {
                ++left;
            }
        }
        else {
            right += offset + 1;
            if (left == right) {
                ++right;
            }
        }
        offset = 0;
    }
    return std::min(left, right);
}
