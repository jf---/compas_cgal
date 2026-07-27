#pragma once

#include "segment_site_delaunay.h"

#include <stdexcept>

class ConsumedCanonicalMatDelaunaySourceError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class InvalidCanonicalMatVoronoiTransferError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class CanonicalMatVoronoiSource2 {
public:
    CanonicalMatVoronoiSource2(
        const CanonicalMatVoronoiSource2&) = delete;
    CanonicalMatVoronoiSource2&
    operator=(
        const CanonicalMatVoronoiSource2&) = delete;
    CanonicalMatVoronoiSource2(
        CanonicalMatVoronoiSource2&&) = delete;
    CanonicalMatVoronoiSource2&
    operator=(
        CanonicalMatVoronoiSource2&&) = delete;
    static CanonicalMatVoronoiSource2 build(
        CanonicalMatDelaunaySource2&& source);
    const SegmentSiteVoronoi2&
    voronoi() const noexcept;
    const CanonicalMatSiteGeometryIndex2&
    site_index() const noexcept;

private:
    explicit CanonicalMatVoronoiSource2(
        CanonicalMatDelaunaySource2&& source);

    SegmentSiteVoronoi2 voronoi_;
    CanonicalMatSiteGeometryIndex2 site_index_;
};
