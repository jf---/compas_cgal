#pragma once

#include "reachable_input_2.h"
#include "segment_site_provenance.h"

#include <array>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

enum class MatSiteKind2 : std::int64_t {
    Point = 0,
    OpenSegment = 1,
};

struct MatSiteProvenance2 {
    MatSiteKind2 kind;
    std::int64_t ring;
    std::int64_t feature;

    bool operator==(
        const MatSiteProvenance2&) const = default;
};

class MatSiteCatalogOverflowError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class InvalidCanonicalMatSiteCatalogError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class UnknownMatSiteIdentityError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class CanonicalMatSiteCatalog2;

class CanonicalMatSite2 {
public:
    const std::string& stable_id() const noexcept;
    const MatSiteProvenance2& provenance() const noexcept;
    const MatTraits::Site_2& site() const noexcept;
    const std::optional<std::array<std::string, 2>>&
    endpoint_site_ids() const noexcept;

private:
    static CanonicalMatSite2 build_point(
        std::int64_t ring,
        std::int64_t feature,
        const MatTraits::Point_2& point);
    static CanonicalMatSite2 build_open_segment(
        std::int64_t ring,
        std::int64_t feature,
        std::int64_t target_feature,
        const MatTraits::Point_2& source,
        const MatTraits::Point_2& target);
    CanonicalMatSite2(
        std::string stable_id,
        MatSiteProvenance2 provenance,
        MatTraits::Site_2 site,
        std::optional<std::array<std::string, 2>>
            endpoint_site_ids);

    std::string stable_id_;
    MatSiteProvenance2 provenance_;
    MatTraits::Site_2 site_;
    std::optional<std::array<std::string, 2>>
        endpoint_site_ids_;

    friend CanonicalMatSiteCatalog2
    canonical_mat_site_catalog(
        const CanonicalReachInput2& input);
};

class CanonicalMatSiteCatalog2 {
public:
    const std::vector<CanonicalMatSite2>&
    sites() const noexcept;
    std::vector<std::array<std::int64_t, 3>>
    site_provenance() const;
    std::int64_t index_of(
        const std::string& stable_id) const;

private:
    explicit CanonicalMatSiteCatalog2(
        std::vector<CanonicalMatSite2> sites);

    std::vector<CanonicalMatSite2> sites_;

    friend CanonicalMatSiteCatalog2
    canonical_mat_site_catalog(
        const CanonicalReachInput2& input);
};

CanonicalMatSiteCatalog2
canonical_mat_site_catalog(
    const CanonicalReachInput2& input);
