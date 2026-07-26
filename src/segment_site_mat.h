#pragma once

#include "exact_algebraic_1.h"

#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/CORE/BigRat.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Parabola_segment_2.h>
#include <CGAL/Segment_Delaunay_graph_2.h>
#include <CGAL/Segment_Delaunay_graph_adaptation_policies_2.h>
#include <CGAL/Segment_Delaunay_graph_adaptation_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_traits_2.h>
#include <CGAL/Voronoi_diagram_2.h>

using MatKernel =
    CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt;
using MatTraits =
    CGAL::Segment_Delaunay_graph_traits_without_intersections_2<
        MatKernel,
        CGAL::Field_with_sqrt_tag>;
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
using SegmentSiteParabola2 =
    CGAL::Parabola_segment_2<MatTraits>;

struct MatParameterDomain2 {
    std::optional<MatTraits::FT> lower;
    std::optional<MatTraits::FT> upper;
};

MatParameterDomain2
exact_parameter_domain(const MatTraits::Line_2& line);

MatParameterDomain2
exact_parameter_domain(const MatTraits::Ray_2& ray);

MatParameterDomain2
exact_parameter_domain(const MatTraits::Segment_2& segment);

MatParameterDomain2
exact_parameter_domain(const SegmentSiteParabola2& parabola);

struct RationalPrimitiveParameterization2 {
    std::vector<CORE::BigRat> x_coefficients;
    std::vector<CORE::BigRat> y_coefficients;
    std::optional<CORE::BigRat> domain_lower;
    std::optional<CORE::BigRat> domain_upper;
};

class InvalidRationalPrimitiveError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

struct ClearanceRootBoundary2 {
    std::optional<CGAL::Sign> constant_sign;
    std::vector<ExactAlgebraicInteger1> primitive_coefficients;
    std::vector<ExactAlgebraicKernel1::Algebraic_real_1> roots;
};

ClearanceRootBoundary2 point_clearance_boundary(
    const RationalPrimitiveParameterization2& primitive,
    const CORE::BigRat& site_x,
    const CORE::BigRat& site_y,
    const CORE::BigRat& radius_squared);

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

std::vector<MatAdmissibleComponent2>
maximal_clearance_components(
    const std::string& original_dual_id,
    const RationalPrimitiveParameterization2& primitive,
    const ClearanceRootBoundary2& boundary);

struct SegmentSiteMatCompileEvidence2 {
    bool delaunay_valid;
    std::size_t assigned_dual_primitives;
    std::size_t exact_clearance_roots;
    std::size_t matched_generator_sites;
    std::size_t delaunay_vertices;
};

SegmentSiteMatCompileEvidence2
segment_site_mat_compile_spike();
