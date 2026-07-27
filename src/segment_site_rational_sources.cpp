#include "segment_site_rational_sources.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <exception>
#include <optional>
#include <string>
#include <utility>

namespace {

const CanonicalMatRationalPointSource2&
point_source(
    const std::vector<
        CanonicalMatRationalPointSource2>& points,
    const std::string& stable_id)
{
    const auto found = std::lower_bound(
        points.begin(),
        points.end(),
        stable_id,
        [](const CanonicalMatRationalPointSource2&
               point,
           const std::string& identity) {
            return point.stable_site_id < identity;
        });
    if (found == points.end()
        || found->stable_site_id != stable_id) {
        throw MissingCanonicalMatEndpointSourceError(
            "canonical MAT segment endpoint has no rational point source");
    }
    return *found;
}

MatTraits::Point_2 exact_point(
    const CanonicalMatRationalPointSource2& point)
{
    return {
        MatTraits::FT(point.x),
        MatTraits::FT(point.y),
    };
}

} // namespace

CORE::BigRat exact_mat_input_rational(
    const MatTraits::FT& coordinate)
{
    CORE::BigRat candidate;
    try {
        candidate = coordinate.BigRatValue();
    } catch (const std::exception&) {
        throw NonRationalCanonicalMatCoordinateError(
            "canonical MAT input coordinate has no rational value");
    }
    if (MatTraits::FT(candidate) != coordinate) {
        throw NonRationalCanonicalMatCoordinateError(
            "canonical MAT input coordinate is not exactly rational");
    }
    return candidate;
}

CanonicalMatRationalPointSource2::
    CanonicalMatRationalPointSource2(
        std::string stable_site_id_,
        CORE::BigRat x_,
        CORE::BigRat y_)
    : stable_site_id(
          std::move(stable_site_id_)),
      x(std::move(x_)),
      y(std::move(y_))
{
}

CanonicalMatRationalOpenSegmentSource2::
    CanonicalMatRationalOpenSegmentSource2(
        std::string stable_site_id_,
        std::string source_point_id_,
        std::string target_point_id_,
        MatExactOpenSegmentSource2 support_)
    : stable_site_id(
          std::move(stable_site_id_)),
      source_point_id(
          std::move(source_point_id_)),
      target_point_id(
          std::move(target_point_id_)),
      support(std::move(support_))
{
}

CanonicalMatRationalSources2
CanonicalMatRationalSources2::build(
    const CanonicalMatSiteCatalog2& catalog)
{
    return CanonicalMatRationalSources2(
        catalog);
}

CanonicalMatRationalSources2::
    CanonicalMatRationalSources2(
        const CanonicalMatSiteCatalog2& catalog)
{
    points_.reserve(
        catalog.sites().size() / 2);
    segments_.reserve(
        catalog.sites().size() / 2);
    for (const CanonicalMatSite2& site :
         catalog.sites()) {
        if (!site.site().is_point()) {
            continue;
        }
        points_.push_back(
            CanonicalMatRationalPointSource2(
                site.stable_id(),
                exact_mat_input_rational(
                    site.site().point().x()),
                exact_mat_input_rational(
                    site.site().point().y())));
    }

    for (const CanonicalMatSite2& site :
         catalog.sites()) {
        if (!site.site().is_segment()) {
            continue;
        }
        const std::optional<
            std::array<std::string, 2>>&
            endpoint_ids =
                site.endpoint_site_ids();
        if (!endpoint_ids.has_value()) {
            throw MissingCanonicalMatEndpointSourceError(
                "canonical MAT segment has no endpoint ownership");
        }
        const CanonicalMatRationalPointSource2&
            source =
                point_source(
                    points_,
                    (*endpoint_ids)[0]);
        const CanonicalMatRationalPointSource2&
            target =
                point_source(
                    points_,
                    (*endpoint_ids)[1]);
        if (exact_point(source)
                != site.site().source()
            || exact_point(target)
                != site.site().target()) {
            throw MismatchedCanonicalMatEndpointGeometryError(
                "canonical MAT endpoint ownership contradicts exact segment geometry");
        }
        segments_.push_back(
            CanonicalMatRationalOpenSegmentSource2(
                site.stable_id(),
                source.stable_site_id,
                target.stable_site_id,
                canonical_open_segment_source(
                    site.stable_id(),
                    source.x,
                    source.y,
                    target.x,
                    target.y)));
    }
}

const std::vector<
    CanonicalMatRationalPointSource2>&
CanonicalMatRationalSources2::points() const noexcept
{
    return points_;
}

const std::vector<
    CanonicalMatRationalOpenSegmentSource2>&
CanonicalMatRationalSources2::segments() const noexcept
{
    return segments_;
}
