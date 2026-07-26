#pragma once

#include "segment_site_parameterization.h"

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

#include <CGAL/Segment_Delaunay_graph_2.h>

using SegmentSiteDelaunay2 =
    CGAL::Segment_Delaunay_graph_2<MatTraits>;

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

std::size_t matched_generator_site_count(
    const SegmentSiteDelaunay2& delaunay,
    const std::vector<GeneratorSite2>& generators);
std::string algebraic_root_identity_v1(
    const ExactAlgebraicKernel1::Algebraic_real_1& root);
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
