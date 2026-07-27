#include "segment_site_mat.h"

#include <algorithm>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace {

struct SegmentLimiterFixture2 {
    MatTraits::Point_2 segment_source;
    MatTraits::Point_2 segment_target;
    MatTraits::Point_2 focus;
    MatTraits::Point_2 limiter_source;
    MatTraits::Point_2 limiter_target;
    MatExactPointSiteSource2 focus_record;
    MatExactOpenSegmentSource2 segment_record;
    MatExactPointSiteSource2 segment_source_record;
    MatExactPointSiteSource2 segment_target_record;
    MatExactOpenSegmentSource2 limiter_record;
    MatExactPointSiteSource2 limiter_source_record;
    MatExactPointSiteSource2 limiter_target_record;
};

SegmentSiteDelaunay2 segment_limiter_delaunay(
    const SegmentLimiterFixture2& source)
{
    SegmentSiteDelaunay2 delaunay;
    delaunay.insert(
        source.segment_source,
        source.segment_target);
    delaunay.insert(source.focus);
    delaunay.insert(
        source.limiter_source,
        source.limiter_target);
    return delaunay;
}

struct LiveSegmentLimiterFixture2 {
    explicit LiveSegmentLimiterFixture2(
        SegmentLimiterFixture2 source)
        : source(std::move(source)),
          delaunay(segment_limiter_delaunay(this->source)),
          voronoi(delaunay)
    {
    }

    SegmentLimiterFixture2 source;
    SegmentSiteDelaunay2 delaunay;
    SegmentSiteVoronoi2 voronoi;
};

struct BoundEndpointObservation2 {
    MatParameterEndpoint2 endpoint;
    SegmentSiteVoronoi2::Halfedge_handle halfedge;
    MatTraits::Point_2 point;
};

bool exact_site_equal(
    const MatTraits::Site_2& lhs,
    const MatTraits::Site_2& rhs)
{
    if (lhs.is_point() != rhs.is_point()) {
        return false;
    }
    if (lhs.is_point()) {
        return lhs.point() == rhs.point();
    }
    return (lhs.source() == rhs.source()
            && lhs.target() == rhs.target())
        || (lhs.source() == rhs.target()
            && lhs.target() == rhs.source());
}

std::vector<BoundEndpointObservation2>
bind_segment_limiter_fixture(
    LiveSegmentLimiterFixture2& fixture)
{
    const MatTraits::Site_2 focus_site =
        MatTraits::Site_2::construct_site_2(
            fixture.source.focus);
    const MatTraits::Site_2 segment_site =
        MatTraits::Site_2::construct_site_2(
            fixture.source.segment_source,
            fixture.source.segment_target);
    const MatTraits::Site_2 limiter_site =
        MatTraits::Site_2::construct_site_2(
            fixture.source.limiter_source,
            fixture.source.limiter_target);
    std::vector<BoundEndpointObservation2> result;
    for (auto halfedge = fixture.voronoi.halfedges_begin();
         halfedge != fixture.voronoi.halfedges_end();
         ++halfedge)
    {
        const MatTraits::Site_2 up = halfedge->up()->site();
        const MatTraits::Site_2 down = halfedge->down()->site();
        if (!((exact_site_equal(up, focus_site)
               && exact_site_equal(down, segment_site))
              || (exact_site_equal(up, segment_site)
                  && exact_site_equal(down, focus_site))))
        {
            continue;
        }
        const bool limiter_is_left =
            halfedge->has_source()
            && exact_site_equal(
                halfedge->left()->site(),
                limiter_site);
        const bool limiter_is_right =
            halfedge->has_target()
            && exact_site_equal(
                halfedge->right()->site(),
                limiter_site);
        if (limiter_is_left == limiter_is_right) {
            continue;
        }
        const MatTraits::Point_2 point =
            limiter_is_left
            ? halfedge->source()->point()
            : halfedge->target()->point();
        SegmentSiteParabola2 parabola;
        if (!CGAL::assign(
                parabola,
                fixture.voronoi.dual().primal(
                    halfedge->dual()))
            || (point != parabola.p1
                && point != parabola.p2))
        {
            continue;
        }
        result.push_back(
            {
                bind_segment_limiter_parabola_endpoint(
                    fixture.source.focus_record,
                    fixture.source.segment_record,
                    fixture.source.segment_source_record,
                    fixture.source.segment_target_record,
                    fixture.source.limiter_record,
                    fixture.source.limiter_source_record,
                    fixture.source.limiter_target_record,
                    fixture.voronoi,
                    halfedge),
                halfedge,
                point,
            });
    }
    return result;
}

std::string canonical_endpoint_bytes(
    const MatParameterEndpoint2& endpoint)
{
    if (!endpoint.parameter.has_value()) {
        return "unbounded";
    }
    std::string result =
        algebraic_root_identity_v1(*endpoint.parameter);
    for (const std::string& provenance : endpoint.provenance_ids) {
        result += ":";
        result += std::to_string(provenance.size());
        result += ":";
        result += provenance;
    }
    return result;
}

SegmentLimiterFixture2 real_segment_limiter_fixture()
{
    return {
        {-100, 0},
        {100, 0},
        {1, 3},
        {0, 1},
        {0, 5},
        {
            "segment-focus",
            {1, 0},
            {3, 0},
            1,
        },
        canonical_open_segment_source(
            "segment-open-segment",
            -100,
            0,
            100,
            0),
        {
            "segment-source",
            {-100, 0},
            {0, 0},
            1,
        },
        {
            "segment-target",
            {100, 0},
            {0, 0},
            1,
        },
        canonical_open_segment_source(
            "segment-limiter",
            0,
            1,
            0,
            5),
        {
            "limiter-source",
            {0, 0},
            {1, 0},
            1,
        },
        {
            "limiter-target",
            {0, 0},
            {5, 0},
            1,
        },
    };
}

SegmentLimiterFixture2 endpoint_reduction_fixture()
{
    return {
        {-100, 0},
        {100, 0},
        {1, 3},
        {CORE::BigRat(-1, 2), CORE::BigRat(3, 2)},
        {CORE::BigRat(-1, 2), 5},
        {
            "endpoint-focus",
            {1, 0},
            {3, 0},
            1,
        },
        canonical_open_segment_source(
            "endpoint-open-segment",
            -100,
            0,
            100,
            0),
        {
            "endpoint-segment-source",
            {-100, 0},
            {0, 0},
            1,
        },
        {
            "endpoint-segment-target",
            {100, 0},
            {0, 0},
            1,
        },
        canonical_open_segment_source(
            "endpoint-limiter-segment",
            CORE::BigRat(-1, 2),
            CORE::BigRat(3, 2),
            CORE::BigRat(-1, 2),
            5),
        {
            "endpoint-limiter",
            {CORE::BigRat(-1, 2), 0},
            {CORE::BigRat(3, 2), 0},
            1,
        },
        {
            "endpoint-limiter-target",
            {CORE::BigRat(-1, 2), 0},
            {5, 0},
            1,
        },
    };
}

SegmentLimiterFixture2 identical_support_fixture()
{
    return {
        {-2, 0},
        {2, 0},
        {2, 3},
        {2, 0},
        {4, 0},
        {
            "identical-focus",
            {2, 0},
            {3, 0},
            1,
        },
        canonical_open_segment_source(
            "identical-open-segment",
            -2,
            0,
            2,
            0),
        {
            "identical-segment-source",
            {-2, 0},
            {0, 0},
            1,
        },
        {
            "identical-segment-target",
            {2, 0},
            {0, 0},
            1,
        },
        canonical_open_segment_source(
            "identical-limiter",
            2,
            0,
            4,
            0),
        {
            "identical-limiter-source",
            {2, 0},
            {0, 0},
            1,
        },
        {
            "identical-limiter-target",
            {4, 0},
            {0, 0},
            1,
        },
    };
}

bool segment_limiter_binds_square_free_endpoint()
{
    LiveSegmentLimiterFixture2 fixture(
        real_segment_limiter_fixture());
    const std::vector<BoundEndpointObservation2> observations =
        bind_segment_limiter_fixture(fixture);
    if (observations.size() != 2
        || !observations[0].endpoint.parameter.has_value()) {
        return false;
    }
    const std::vector<ExactAlgebraicInteger1>
        governing_factor{45, 72, -18, 0, 1};
    const std::string expected_root_id =
        algebraic_root_id_v1(governing_factor, 1);
    return observations[0].point == observations[1].point
        && observations[0].halfedge->opposite()
            == observations[1].halfedge
        && std::find(
               observations[0].endpoint.provenance_ids.begin(),
               observations[0].endpoint.provenance_ids.end(),
               expected_root_id)
            != observations[0].endpoint.provenance_ids.end()
        && std::find(
               observations[0].endpoint.provenance_ids.begin(),
               observations[0].endpoint.provenance_ids.end(),
               "segment-limiter/norm-factor-multiplicity/2")
            != observations[0].endpoint.provenance_ids.end()
        && canonical_endpoint_bytes(observations[0].endpoint)
            == canonical_endpoint_bytes(
                observations[1].endpoint);
}

bool segment_endpoint_reduces_to_point_limiter()
{
    LiveSegmentLimiterFixture2 fixture(
        endpoint_reduction_fixture());
    const MatTraits::Site_2 focus_site =
        MatTraits::Site_2::construct_site_2(
            fixture.source.focus);
    const MatTraits::Site_2 segment_site =
        MatTraits::Site_2::construct_site_2(
            fixture.source.segment_source,
            fixture.source.segment_target);
    const MatTraits::Site_2 endpoint_site =
        MatTraits::Site_2::construct_site_2(
            fixture.source.limiter_source);
    std::vector<MatParameterEndpoint2> endpoints;
    std::size_t segment_rejections = 0;
    for (auto halfedge = fixture.voronoi.halfedges_begin();
         halfedge != fixture.voronoi.halfedges_end();
         ++halfedge)
    {
        const MatTraits::Site_2 up = halfedge->up()->site();
        const MatTraits::Site_2 down = halfedge->down()->site();
        if (!((exact_site_equal(up, focus_site)
               && exact_site_equal(down, segment_site))
              || (exact_site_equal(up, segment_site)
                  && exact_site_equal(down, focus_site))))
        {
            continue;
        }
        const bool endpoint_is_left =
            halfedge->has_source()
            && exact_site_equal(
                halfedge->left()->site(),
                endpoint_site);
        const bool endpoint_is_right =
            halfedge->has_target()
            && exact_site_equal(
                halfedge->right()->site(),
                endpoint_site);
        if (!endpoint_is_left
            && !endpoint_is_right)
        {
            continue;
        }
        endpoints.push_back(
            bind_point_limiter_parabola_endpoint(
                fixture.source.focus_record,
                fixture.source.segment_record,
                fixture.source.segment_source_record,
                fixture.source.segment_target_record,
                fixture.source.limiter_source_record,
                fixture.voronoi,
                halfedge));
        try {
            bind_segment_limiter_parabola_endpoint(
                fixture.source.focus_record,
                fixture.source.segment_record,
                fixture.source.segment_source_record,
                fixture.source.segment_target_record,
                fixture.source.limiter_record,
                fixture.source.limiter_source_record,
                fixture.source.limiter_target_record,
                fixture.voronoi,
                halfedge);
        }
        catch (const MismatchedLiveParabolaBridgeError&) {
            ++segment_rejections;
        }
    }
    return endpoints.size() == 2
        && segment_rejections == 2
        && canonical_endpoint_bytes(endpoints[0])
            == canonical_endpoint_bytes(endpoints[1]);
}

bool unrelated_live_segment_limiter_is_rejected()
{
    LiveSegmentLimiterFixture2 fixture(
        real_segment_limiter_fixture());
    const std::vector<BoundEndpointObservation2> observations =
        bind_segment_limiter_fixture(fixture);
    if (observations.empty()) {
        return false;
    }
    MatExactPointSiteSource2 unrelated_source =
        fixture.source.limiter_source_record;
    MatExactPointSiteSource2 unrelated_target =
        fixture.source.limiter_target_record;
    unrelated_source.y = {9, 0};
    unrelated_target.y = {10, 0};
    try {
        bind_segment_limiter_parabola_endpoint(
            fixture.source.focus_record,
            fixture.source.segment_record,
            fixture.source.segment_source_record,
            fixture.source.segment_target_record,
            fixture.source.limiter_record,
            unrelated_source,
            unrelated_target,
            fixture.voronoi,
            observations[0].halfedge);
    }
    catch (const MismatchedLiveParabolaBridgeError&) {
        return true;
    }
    return false;
}

bool detached_identical_segment_limiter_fails_provenance()
{
    LiveSegmentLimiterFixture2 fixture(
        real_segment_limiter_fixture());
    const std::vector<BoundEndpointObservation2> observations =
        bind_segment_limiter_fixture(fixture);
    if (observations.empty()) {
        return false;
    }
    try {
        bind_segment_limiter_parabola_endpoint(
            fixture.source.focus_record,
            fixture.source.segment_record,
            fixture.source.segment_source_record,
            fixture.source.segment_target_record,
            fixture.source.segment_record,
            fixture.source.segment_source_record,
            fixture.source.segment_target_record,
            fixture.voronoi,
            observations[0].halfedge);
    }
    catch (const MismatchedLiveParabolaBridgeError&) {
        return true;
    }
    return false;
}

bool identical_support_has_no_segment_limiter_owner()
{
    LiveSegmentLimiterFixture2 fixture(
        identical_support_fixture());
    const MatTraits::Site_2 focus_site =
        MatTraits::Site_2::construct_site_2(
            fixture.source.focus);
    const MatTraits::Site_2 segment_site =
        MatTraits::Site_2::construct_site_2(
            fixture.source.segment_source,
            fixture.source.segment_target);
    const MatTraits::Site_2 limiter_site =
        MatTraits::Site_2::construct_site_2(
            fixture.source.limiter_source,
            fixture.source.limiter_target);
    std::size_t owned_halfedges = 0;
    for (auto halfedge = fixture.voronoi.halfedges_begin();
         halfedge != fixture.voronoi.halfedges_end();
         ++halfedge)
    {
        const MatTraits::Site_2 up =
            halfedge->up()->site();
        const MatTraits::Site_2 down =
            halfedge->down()->site();
        const bool generators_match =
            (exact_site_equal(up, focus_site)
             && exact_site_equal(down, segment_site))
            || (exact_site_equal(up, segment_site)
                && exact_site_equal(down, focus_site));
        const bool limiter_owns_endpoint =
            (halfedge->has_source()
             && exact_site_equal(
                 halfedge->left()->site(),
                 limiter_site))
            || (halfedge->has_target()
                && exact_site_equal(
                    halfedge->right()->site(),
                    limiter_site));
        if (generators_match
            && limiter_owns_endpoint)
        {
            ++owned_halfedges;
        }
    }
    return owned_halfedges == 0;
}

bool square_field_radical_segment_is_rejected()
{
    LiveSegmentLimiterFixture2 fixture(
        real_segment_limiter_fixture());
    const std::vector<BoundEndpointObservation2> observations =
        bind_segment_limiter_fixture(fixture);
    if (observations.empty()) {
        return false;
    }
    fixture.source.limiter_source_record.y.radical = 1;
    try {
        bind_segment_limiter_parabola_endpoint(
            fixture.source.focus_record,
            fixture.source.segment_record,
            fixture.source.segment_source_record,
            fixture.source.segment_target_record,
            fixture.source.limiter_record,
            fixture.source.limiter_source_record,
            fixture.source.limiter_target_record,
            fixture.voronoi,
            observations[0].halfedge);
    }
    catch (const NonCanonicalQuadraticFieldError&) {
        return true;
    }
    return false;
}

bool coincident_limiter_segment_is_rejected()
{
    LiveSegmentLimiterFixture2 fixture(
        real_segment_limiter_fixture());
    const std::vector<BoundEndpointObservation2> observations =
        bind_segment_limiter_fixture(fixture);
    if (observations.empty()) {
        return false;
    }
    try {
        bind_segment_limiter_parabola_endpoint(
            fixture.source.focus_record,
            fixture.source.segment_record,
            fixture.source.segment_source_record,
            fixture.source.segment_target_record,
            fixture.source.limiter_record,
            fixture.source.limiter_source_record,
            fixture.source.limiter_source_record,
            fixture.voronoi,
            observations[0].halfedge);
    }
    catch (const CoincidentSegmentEndpointsError&) {
        return true;
    }
    return false;
}

bool degenerate_limiter_is_rejected_at_construction()
{
    try {
        static_cast<void>(
            canonical_open_segment_source(
                "zero-limiter-line",
                CORE::BigRat(-2, 5),
                CORE::BigRat(7, 11),
                CORE::BigRat(-2, 5),
                CORE::BigRat(7, 11)));
    }
    catch (const DegenerateOpenSegmentSourceError&) {
        return true;
    }
    return false;
}

} // namespace

bool segment_limiter_gate()
{
    return segment_limiter_binds_square_free_endpoint()
        && segment_endpoint_reduces_to_point_limiter()
        && unrelated_live_segment_limiter_is_rejected()
        && detached_identical_segment_limiter_fails_provenance()
        && identical_support_has_no_segment_limiter_owner()
        && square_field_radical_segment_is_rejected()
        && coincident_limiter_segment_is_rejected()
        && degenerate_limiter_is_rejected_at_construction();
}
