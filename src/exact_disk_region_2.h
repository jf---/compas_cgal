#pragma once

#include "stock_2.h"

#include <stdexcept>
#include <vector>

class ExactDiskRegionEmptyCentersError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

class ExactDiskRegionRadiusError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

// The `_into` form joins the disk union into a set the caller already owns, and
// the by-value form delegates to it. A Gps has no move constructor (its base
// declares a virtual destructor), so a set returned by value can only reach
// long-lived storage through CGAL's copy constructor -- which leaves the copy's
// arrangement reading the traits of the returned object, and that object then
// dies. A caller that must KEEP the union alive builds it with the `_into` form.
void build_exact_disk_union_region_2_into(
    Gps& target,
    const std::vector<EPoint>& centers,
    const Epeck::FT& radius);
Gps build_exact_disk_union_region_2(
    const std::vector<EPoint>& centers,
    const Epeck::FT& radius);
