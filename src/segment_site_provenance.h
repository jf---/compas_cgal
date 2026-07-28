#pragma once

#include "segment_site_parameterization.h"

#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string>
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

struct MatParameterEndpoint2 {
    std::optional<ExactAlgebraicKernel1::Algebraic_real_1>
        parameter;
    std::vector<std::string> provenance_ids;
};

struct MatAdmissibleComponent2 {
    std::string component_id;
    MatParameterEndpoint2 lower;
    MatParameterEndpoint2 upper;
};

struct RationalDomainRoot2 {
    CORE::BigRat parameter;
    std::vector<std::string> provenance_ids;
};

struct AlgebraicDomainRoot2 {
    ExactAlgebraicKernel1::Algebraic_real_1 parameter;
    std::vector<std::string> provenance_ids;
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
std::string stable_dual_identity_v1(
    const std::string& dual_kind,
    const std::vector<std::string>& ordered_generator_ids);
std::string stable_voronoi_node_identity_v1(
    const std::vector<std::string>& ordered_generator_ids);
std::string stable_normalized_voronoi_node_identity_v1(
    const std::vector<std::string>& ordered_generator_ids);
std::string algebraic_root_identity_v1(
    const ExactAlgebraicKernel1::Algebraic_real_1& root);
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
