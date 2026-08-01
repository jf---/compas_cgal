#pragma once

#include "segment_site_parameterization.h"

#include <compare>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <CGAL/Segment_Delaunay_graph_2.h>
#include <CGAL/Segment_Delaunay_graph_adaptation_policies_2.h>
#include <CGAL/Segment_Delaunay_graph_adaptation_traits_2.h>
#include <CGAL/Voronoi_diagram_2.h>

using SegmentSiteDelaunay2 =
    CGAL::Segment_Delaunay_graph_2<MatTraits>;
using SegmentSiteAdaptationTraits2 =
    CGAL::Segment_Delaunay_graph_adaptation_traits_2<
        SegmentSiteDelaunay2>;
using SegmentSiteAdaptationPolicy2 =
    CGAL::Segment_Delaunay_graph_degeneracy_removal_policy_2<
        SegmentSiteDelaunay2>;
using SegmentSiteVoronoi2 = CGAL::Voronoi_diagram_2<
    SegmentSiteDelaunay2,
    SegmentSiteAdaptationTraits2,
    SegmentSiteAdaptationPolicy2>;

struct GeneratorSite2 {
    std::string stable_id;
    MatTraits::Site_2 site;
};

enum class MatExactCurveKind2 : std::int64_t {
    Line = 0,
    Parabola = 1,
    Circle = 2,
};

enum class MatEndpointDomainKind2 : std::int64_t {
    Design = 0,
    Clearance = 1,
};

inline constexpr std::int64_t MAT_CURVE_LINE =
    static_cast<std::int64_t>(
        MatExactCurveKind2::Line);
inline constexpr std::int64_t MAT_CURVE_PARABOLA =
    static_cast<std::int64_t>(
        MatExactCurveKind2::Parabola);
inline constexpr std::int64_t MAT_CURVE_CIRCLE =
    static_cast<std::int64_t>(
        MatExactCurveKind2::Circle);
inline constexpr std::int64_t MAT_ENDPOINT_DOMAIN_DESIGN =
    static_cast<std::int64_t>(
        MatEndpointDomainKind2::Design);
inline constexpr std::int64_t MAT_ENDPOINT_DOMAIN_CLEARANCE =
    static_cast<std::int64_t>(
        MatEndpointDomainKind2::Clearance);

struct MatEndpointBoundaryFeature2 {
    MatEndpointDomainKind2 domain_kind;
    std::int64_t component;
    MatExactCurveKind2 curve_kind;
    std::int64_t source_site_or_ring;
    std::int64_t derived_feature_index;

    auto operator<=>(
        const MatEndpointBoundaryFeature2&) const = default;
};

struct MatEndpointExactEvidence2 {
    bool original_voronoi_vertex = false;
    bool domain_boundary = false;
    bool clearance_root = false;
    std::vector<MatEndpointBoundaryFeature2>
        boundary_features;

    bool operator==(
        const MatEndpointExactEvidence2&) const = default;
};

class MatEndpointRootIdentity2 {
public:
    static MatEndpointRootIdentity2 build(
        const ExactAlgebraicKernel1::Algebraic_real_1& root);

    const std::string& root_id() const noexcept
    {
        return root_id_;
    }

    void require_matches(
        const ExactAlgebraicKernel1::Algebraic_real_1& root) const;

private:
    MatEndpointRootIdentity2(
        std::vector<ExactAlgebraicInteger1> factor_coefficients,
        std::size_t root_ordinal,
        std::string root_id)
        : factor_coefficients_(
              std::move(factor_coefficients)),
          root_ordinal_(root_ordinal),
          root_id_(std::move(root_id))
    {
    }

    std::vector<ExactAlgebraicInteger1>
        factor_coefficients_;
    std::size_t root_ordinal_;
    std::string root_id_;
};

struct MatParameterEndpoint2 {
    std::optional<ExactAlgebraicKernel1::Algebraic_real_1>
        parameter;
    std::vector<std::string> provenance_ids;
    MatEndpointExactEvidence2 exact_evidence;
    std::optional<MatEndpointRootIdentity2>
        parameter_root_identity;
};

struct MatAdmissibleComponent2 {
    std::string component_id;
    MatParameterEndpoint2 lower;
    MatParameterEndpoint2 upper;
    std::int64_t component_index = -1;
};

struct RationalDomainRoot2 {
    CORE::BigRat parameter;
    std::vector<std::string> provenance_ids;
    std::vector<MatEndpointBoundaryFeature2>
        boundary_features;
};

struct AlgebraicDomainRoot2 {
    ExactAlgebraicKernel1::Algebraic_real_1 parameter;
    std::vector<std::string> provenance_ids;
    std::vector<MatEndpointBoundaryFeature2>
        boundary_features;
};

struct TaggedEndpoint2 {
    MatParameterEndpoint2 endpoint;
    bool clearance_root;
};

class GeneratorSiteBijectionError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class InvalidDualIdentityError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class MatEndpointFeatureOverflowError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class MatGraphComponentIndexOverflowError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class InvalidMatEndpointRootIdentityError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

std::size_t matched_generator_site_count(
    const SegmentSiteDelaunay2& delaunay,
    const std::vector<GeneratorSite2>& generators);
std::size_t require_generator_site_bijection(
    const SegmentSiteDelaunay2& delaunay,
    const std::vector<GeneratorSite2>& generators);
std::string stable_generator_site_id(
    const MatTraits::Site_2& site,
    const std::vector<GeneratorSite2>& generators);
std::vector<std::string> ordered_generator_site_ids(
    std::string first,
    std::string second);
void union_stable_ids(
    std::vector<std::string>& target,
    const std::vector<std::string>& source);
MatEndpointBoundaryFeature2 design_line_endpoint_feature(
    std::size_t ring,
    std::size_t feature);
void union_endpoint_boundary_features(
    std::vector<MatEndpointBoundaryFeature2>& target,
    const std::vector<MatEndpointBoundaryFeature2>& source);
void union_endpoint_evidence(
    MatParameterEndpoint2& target,
    const MatParameterEndpoint2& source);
std::string stable_dual_identity_v1(
    const std::string& dual_kind,
    const std::vector<std::string>& ordered_generator_ids);
std::string stable_voronoi_node_identity_v1(
    const std::vector<std::string>& ordered_generator_ids);
std::string stable_normalized_voronoi_node_identity_v1(
    const std::vector<std::string>& ordered_generator_ids);
std::string algebraic_root_identity_v1(
    const ExactAlgebraicKernel1::Algebraic_real_1& root);
const std::string& mat_endpoint_root_identity_v1(
    const MatParameterEndpoint2& endpoint);
MatParameterEndpoint2 exact_graph_endpoint_binding(
    const MatParameterEndpoint2& endpoint);
std::string stable_endpoint_node_identity_v1(
    const std::string& dual_id,
    const MatParameterEndpoint2& endpoint);
MatParameterEndpoint2 domain_endpoint(
    const std::string& dual_id,
    const char* side,
    const std::optional<CORE::BigRat>& value,
    const ExactAlgebraicKernel1& kernel);
void append_root_endpoint(
    std::vector<TaggedEndpoint2>& endpoints,
    TaggedEndpoint2& upper,
    const ExactAlgebraicKernel1::Algebraic_real_1& root,
    const std::string& root_id,
    const ExactAlgebraicKernel1& kernel);
std::vector<AlgebraicDomainRoot2> parabola_domain_roots(
    const std::string& dual_id,
    const RationalPrimitiveParameterization2& primitive,
    const MatDomainPolygonWithHoles2& domain);
std::vector<AlgebraicDomainRoot2>
nonparallel_segment_domain_roots(
    const std::string& dual_id,
    const NonparallelSegmentBisectorParameterization2& primitive,
    const NonparallelSegmentFeatureDomain2& feature_domain,
    const MatDomainPolygonWithHoles2& domain);
