#pragma once

#include "segment_site_mat.h"

#include <stdexcept>
#include <string>
#include <vector>

class InvalidMatSquaredClearanceFunctionError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class UnsupportedMatSquaredClearanceDegreeError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class NegativeMatSquaredClearanceError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class InvalidMatClearanceEdgeProfileError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class IncompleteMatSquaredWidthError
    : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class MatStrictEdgeClearanceMinimum2;

class MatSquaredClearanceFunction2 {
public:
    static MatSquaredClearanceFunction2 build(
        RationalPolynomial coefficients);

    const RationalPolynomial& coefficients()
        const noexcept;
    bool is_constant() const noexcept;

private:
    explicit MatSquaredClearanceFunction2(
        RationalPolynomial coefficients);

    RationalPolynomial coefficients_;
};

class MatClearanceEdgeProfile2 {
public:
    static MatClearanceEdgeProfile2 build(
        const MatExactGraphEdge2& edge,
        MatSquaredClearanceFunction2
            squared_clearance);

    const std::string& edge_id() const noexcept;
    const std::vector<std::string>&
    defining_site_ids() const noexcept;
    const MatParameterEndpoint2& lower()
        const noexcept;
    const MatParameterEndpoint2& upper()
        const noexcept;
    const MatSquaredClearanceFunction2&
    squared_clearance() const noexcept;

private:
    MatClearanceEdgeProfile2(
        std::string edge_id,
        std::vector<std::string>
            defining_site_ids,
        MatParameterEndpoint2 lower,
        MatParameterEndpoint2 upper,
        MatSquaredClearanceFunction2
            squared_clearance);

    std::string edge_id_;
    std::vector<std::string>
        defining_site_ids_;
    MatParameterEndpoint2 lower_;
    MatParameterEndpoint2 upper_;
    MatSquaredClearanceFunction2
        squared_clearance_;
};

class MatExactSquaredWidth2 {
public:
    static MatExactSquaredWidth2 from_clearance(
        const MatSquaredClearanceFunction2&
            clearance,
        const ExactAlgebraicKernel1::
            Algebraic_real_1& parameter);

    const ExactAlgebraicKernel1::Algebraic_real_1&
    value() const noexcept;
    const std::string& root_id() const noexcept;

private:
    MatExactSquaredWidth2(
        ExactAlgebraicKernel1::Algebraic_real_1
            value,
        std::string root_id);

    ExactAlgebraicKernel1::Algebraic_real_1
        value_;
    std::string root_id_;
};

class MatClearanceEndpointBehavior2 {
public:
    CGAL::Sign inward_clearance_sign()
        const noexcept;
    const MatExactSquaredWidth2&
    squared_width() const noexcept;

private:
    MatClearanceEndpointBehavior2(
        CGAL::Sign inward_clearance_sign,
        MatExactSquaredWidth2 squared_width);

    CGAL::Sign inward_clearance_sign_;
    MatExactSquaredWidth2 squared_width_;

    friend MatClearanceEndpointBehavior2
    lower_endpoint_clearance_behavior(
        const MatClearanceEdgeProfile2&
            profile);
    friend MatClearanceEndpointBehavior2
    upper_endpoint_clearance_behavior(
        const MatClearanceEdgeProfile2&
            profile);
};

MatClearanceEndpointBehavior2
lower_endpoint_clearance_behavior(
    const MatClearanceEdgeProfile2& profile);

MatClearanceEndpointBehavior2
upper_endpoint_clearance_behavior(
    const MatClearanceEdgeProfile2& profile);

class MatStrictEdgeClearanceMinimum2 {
public:
    const std::string& edge_id() const noexcept;
    const std::vector<std::string>&
    defining_site_ids() const noexcept;
    const ExactAlgebraicKernel1::Algebraic_real_1&
    parameter() const noexcept;
    const std::string& parameter_root_id()
        const noexcept;
    CGAL::Sign left_derivative_sign()
        const noexcept;
    CGAL::Sign right_derivative_sign()
        const noexcept;
    const MatExactSquaredWidth2&
    squared_width() const noexcept;

private:
    MatStrictEdgeClearanceMinimum2(
        std::string edge_id,
        std::vector<std::string>
            defining_site_ids,
        ExactAlgebraicKernel1::Algebraic_real_1
            parameter,
        std::string parameter_root_id,
        MatExactSquaredWidth2 squared_width);

    std::string edge_id_;
    std::vector<std::string>
        defining_site_ids_;
    ExactAlgebraicKernel1::Algebraic_real_1
        parameter_;
    std::string parameter_root_id_;
    MatExactSquaredWidth2 squared_width_;

    friend std::vector<
        MatStrictEdgeClearanceMinimum2>
    strict_edge_clearance_minima(
        const MatClearanceEdgeProfile2&
            profile);
};

std::vector<MatStrictEdgeClearanceMinimum2>
strict_edge_clearance_minima(
    const MatClearanceEdgeProfile2& profile);
