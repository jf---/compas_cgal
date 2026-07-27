#pragma once

#include "segment_site_catalog.h"

#include <stdexcept>
#include <string>
#include <vector>

class NonRationalCanonicalMatCoordinateError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class MissingCanonicalMatEndpointSourceError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class MismatchedCanonicalMatEndpointGeometryError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

CORE::BigRat exact_mat_input_rational(
    const MatTraits::FT& coordinate);

class CanonicalMatRationalPointSource2 {
public:
    const std::string stable_site_id;
    const CORE::BigRat x;
    const CORE::BigRat y;

private:
    CanonicalMatRationalPointSource2(
        std::string stable_site_id,
        CORE::BigRat x,
        CORE::BigRat y);

    friend class CanonicalMatRationalSources2;
};

class CanonicalMatRationalOpenSegmentSource2 {
public:
    const std::string stable_site_id;
    const std::string source_point_id;
    const std::string target_point_id;
    const MatExactOpenSegmentSource2 support;

private:
    CanonicalMatRationalOpenSegmentSource2(
        std::string stable_site_id,
        std::string source_point_id,
        std::string target_point_id,
        MatExactOpenSegmentSource2 support);

    friend class CanonicalMatRationalSources2;
};

class CanonicalMatRationalSources2 {
public:
    static CanonicalMatRationalSources2 build(
        const CanonicalMatSiteCatalog2& catalog);
    const std::vector<
        CanonicalMatRationalPointSource2>&
    points() const noexcept;
    const std::vector<
        CanonicalMatRationalOpenSegmentSource2>&
    segments() const noexcept;

private:
    explicit CanonicalMatRationalSources2(
        const CanonicalMatSiteCatalog2& catalog);

    std::vector<
        CanonicalMatRationalPointSource2>
        points_;
    std::vector<
        CanonicalMatRationalOpenSegmentSource2>
        segments_;
};
