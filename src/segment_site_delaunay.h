#pragma once

#include "segment_site_catalog.h"

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

class DuplicateCanonicalMatSiteGeometryError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class UnknownCanonicalMatSiteGeometryError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class CanonicalMatDelaunayBijectionError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class CanonicalMatSiteGeometryIndex2 {
public:
    std::size_t size() const noexcept;
    const std::string& stable_id(
        const MatTraits::Site_2& site) const;

private:
    struct PointRecord {
        MatTraits::Point_2 point;
        std::string stable_id;
    };

    struct SegmentRecord {
        MatTraits::Point_2 first;
        MatTraits::Point_2 second;
        std::string stable_id;
    };

    static bool segment_key_less(
        const SegmentRecord& left,
        const SegmentRecord& right);
    static bool segment_key_equal(
        const SegmentRecord& left,
        const SegmentRecord& right);
    explicit CanonicalMatSiteGeometryIndex2(
        const CanonicalMatSiteCatalog2& catalog);

    std::vector<PointRecord> points_;
    std::vector<SegmentRecord> segments_;

    friend class CanonicalMatDelaunaySource2;
};

class CanonicalMatDelaunaySource2 {
public:
    CanonicalMatDelaunaySource2(
        const CanonicalMatDelaunaySource2&) = delete;
    CanonicalMatDelaunaySource2&
    operator=(
        const CanonicalMatDelaunaySource2&) = delete;
    CanonicalMatDelaunaySource2(
        CanonicalMatDelaunaySource2&&) = delete;
    CanonicalMatDelaunaySource2&
    operator=(
        CanonicalMatDelaunaySource2&&) = delete;
    static CanonicalMatDelaunaySource2 build(
        const CanonicalMatSiteCatalog2& catalog);
    const SegmentSiteDelaunay2&
    delaunay() const noexcept;
    const CanonicalMatSiteGeometryIndex2&
    site_index() const noexcept;

private:
    explicit CanonicalMatDelaunaySource2(
        const CanonicalMatSiteCatalog2& catalog);

    CanonicalMatSiteGeometryIndex2 site_index_;
    SegmentSiteDelaunay2 delaunay_;
};
