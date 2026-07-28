#include "segment_site_neck_clearance.h"

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

namespace {

MatParameterEndpoint2 endpoint(
    const CORE::BigRat& parameter)
{
    ExactAlgebraicKernel1 kernel;
    const auto algebraic =
        kernel.construct_algebraic_real_1_object()(
            parameter);
    return {
        algebraic,
        {
            algebraic_root_identity_v1(
                algebraic),
        },
        {
            true,
            false,
            false,
            {},
        },
    };
}

MatExactGraphEdge2 edge(
    MatParameterEndpoint2 lower,
    MatParameterEndpoint2 upper)
{
    return {
        "edge-a",
        "LINE",
        "dual-a",
        "node-a",
        "node-b",
        std::move(lower),
        std::move(upper),
        {"site-a", "site-b"},
        {"parent-a", "parent-b"},
        0,
        true,
    };
}

MatClearanceEdgeProfile2 profile(
    RationalPolynomial coefficients,
    const CORE::BigRat& lower,
    const CORE::BigRat& upper)
{
    return MatClearanceEdgeProfile2::build(
        edge(
            endpoint(lower),
            endpoint(upper)),
        MatSquaredClearanceFunction2::build(
            std::move(coefficients)));
}

bool quadratic_minimum_is_exact()
{
    const auto minima =
        strict_edge_clearance_minima(
            profile(
                {
                    CORE::BigRat(5, 2),
                    CORE::BigRat(-1),
                    CORE::BigRat(1),
                },
                -1,
                2));
    ExactAlgebraicKernel1 kernel;
    return minima.size() == 1
        && minima.front().edge_id()
            == "edge-a"
        && minima.front()
               .defining_site_ids()
            == std::vector<std::string>{
                "site-a",
                "site-b",
            }
        && kernel.compare_1_object()(
               minima.front().parameter(),
               CORE::BigRat(1, 2))
            == CGAL::EQUAL
        && minima.front().parameter_root_id()
            == algebraic_root_id_v1(
                {-1, 2},
                0)
        && minima.front()
               .left_derivative_sign()
            == CGAL::NEGATIVE
        && minima.front()
               .right_derivative_sign()
            == CGAL::POSITIVE
        && kernel.compare_1_object()(
               minima.front()
                   .squared_width()
                   .value(),
               CORE::BigRat(9))
            == CGAL::EQUAL
        && minima.front()
               .squared_width()
               .root_id()
            == algebraic_root_id_v1(
                {-9, 1},
                0);
}

bool quartic_has_two_algebraic_minima()
{
    const auto minima =
        strict_edge_clearance_minima(
            profile(
                {
                    5,
                    0,
                    -4,
                    0,
                    1,
                },
                -2,
                2));
    if (minima.size() != 2) {
        return false;
    }
    return minima[0].parameter_root_id()
            == algebraic_root_id_v1(
                {-2, 0, 1},
                0)
        && minima[1].parameter_root_id()
            == algebraic_root_id_v1(
                {-2, 0, 1},
                1)
        && minima[0]
               .squared_width()
               .root_id()
            == algebraic_root_id_v1(
                {-4, 1},
                0)
        && minima[1]
               .squared_width()
               .root_id()
            == minima[0]
                   .squared_width()
                   .root_id();
}

bool nonminima_are_not_emitted()
{
    return strict_edge_clearance_minima(
               profile({4}, -1, 1))
               .empty()
        && strict_edge_clearance_minima(
               profile({4, 4, 1}, 0, 1))
               .empty()
        && strict_edge_clearance_minima(
               profile({4, 0, -1}, -1, 1))
               .empty();
}

bool malformed_profiles_are_rejected()
{
    bool empty_function_rejected = false;
    try {
        static_cast<void>(
            MatSquaredClearanceFunction2::build(
                {}));
    } catch (
        const InvalidMatSquaredClearanceFunctionError&) {
        empty_function_rejected = true;
    }

    bool degree_rejected = false;
    try {
        static_cast<void>(
            MatSquaredClearanceFunction2::build(
                {1, 0, 0, 0, 0, 1}));
    } catch (
        const UnsupportedMatSquaredClearanceDegreeError&) {
        degree_rejected = true;
    }

    bool negative_rejected = false;
    try {
        static_cast<void>(
            profile({-1, 0, 1}, -1, 1));
    } catch (
        const NegativeMatSquaredClearanceError&) {
        negative_rejected = true;
    }
    bool negative_leading_rejected = false;
    try {
        static_cast<void>(
            profile({-1, 0, -1}, -1, 1));
    } catch (
        const NegativeMatSquaredClearanceError&) {
        negative_leading_rejected = true;
    }

    bool reversed_rejected = false;
    try {
        static_cast<void>(
            profile({1}, 1, -1));
    } catch (
        const InvalidMatClearanceEdgeProfileError&) {
        reversed_rejected = true;
    }

    bool unbounded_rejected = false;
    try {
        MatParameterEndpoint2 unbounded{
            std::nullopt,
            {},
            {},
        };
        static_cast<void>(
            MatClearanceEdgeProfile2::build(
                edge(
                    std::move(unbounded),
                    endpoint(1)),
                MatSquaredClearanceFunction2::build(
                    {1})));
    } catch (
        const InvalidMatClearanceEdgeProfileError&) {
        unbounded_rejected = true;
    }

    bool identity_rejected = false;
    try {
        MatExactGraphEdge2 malformed =
            edge(
                endpoint(-1),
                endpoint(1));
        malformed.edge_id.clear();
        std::swap(
            malformed.generator_site_ids[0],
            malformed.generator_site_ids[1]);
        static_cast<void>(
            MatClearanceEdgeProfile2::build(
                malformed,
                MatSquaredClearanceFunction2::build(
                    {1})));
    } catch (
        const InvalidMatClearanceEdgeProfileError&) {
        identity_rejected = true;
    }
    bool evidence_rejected = false;
    try {
        MatExactGraphEdge2 malformed =
            edge(
                endpoint(-1),
                endpoint(1));
        malformed.source_endpoint
            .exact_evidence = {};
        static_cast<void>(
            MatClearanceEdgeProfile2::build(
                malformed,
                MatSquaredClearanceFunction2::build(
                    {1})));
    } catch (
        const InvalidMatClearanceEdgeProfileError&) {
        evidence_rejected = true;
    }
    return empty_function_rejected
        && degree_rejected
        && negative_rejected
        && negative_leading_rejected
        && reversed_rejected
        && unbounded_rejected
        && identity_rejected
        && evidence_rejected;
}

} // namespace

bool neck_clearance_gate()
{
    return quadratic_minimum_is_exact()
        && quartic_has_two_algebraic_minima()
        && nonminima_are_not_emitted()
        && malformed_profiles_are_rejected();
}
