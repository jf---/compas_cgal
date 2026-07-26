#include "segment_fibre.h"

#include "../exact_algebraic_1.h"
#include "event_certificate.h"
#include "exact_polynomial_sign_2.h"
#include "segment_partition.h"

#include <algorithm>
#include <iterator>
#include <map>
#include <optional>
#include <string_view>
#include <tuple>
#include <utility>

#include <CGAL/CORE/BigRat.h>
#include <CGAL/Fraction_traits.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>

namespace {

using Integer = CORE::BigInt;
using Rational = CORE::BigRat;
using Kernel1 = ExactAlgebraicKernel1;
using Kernel2 = ExactAlgebraicKernel2;
using Polynomial1 = Kernel1::Polynomial_1;
using Polynomial2 = Kernel2::Polynomial_2;
using AlgebraicReal1 = Kernel1::Algebraic_real_1;
using AlgebraicReal2 = Kernel2::Algebraic_real_2;

struct CoordinateParts {
    Rational base;
    Rational radical_coefficient;
    Rational radicand;
};

struct RadicalPredicate {
    Polynomial2 rational_part;
    Polynomial2 radical_part;
    Rational radicand;
};

struct PointNumerators {
    Polynomial2 x;
    Polynomial2 y;
    Polynomial2 denominator;
};

Rational parse_rational(
    const std::string& text,
    std::string_view role)
{
    const std::size_t separator = text.find('/');
    try {
        if (separator == std::string::npos) {
            return Rational(Integer(text));
        }
        if (text.find('/', separator + 1)
            != std::string::npos) {
            throw IncompleteSegmentPartitionError(
                std::string(role)
                + " is not an exact rational");
        }
        const Integer numerator(
            text.substr(0, separator));
        const Integer denominator(
            text.substr(separator + 1));
        if (denominator == 0) {
            throw IncompleteSegmentPartitionError(
                std::string(role)
                + " has zero denominator");
        }
        return Rational(numerator, denominator);
    } catch (const EventSubstrateError&) {
        throw;
    } catch (const std::exception&) {
        throw IncompleteSegmentPartitionError(
            std::string(role)
            + " is not an exact rational");
    }
}

Rational exact_rational(const Epeck::FT& value)
{
    return CGAL::exact(value);
}

CoordinateParts coordinate_parts(
    const GpsPoint::CoordNT& coordinate)
{
    return {
        exact_rational(coordinate.a0()),
        exact_rational(coordinate.a1()),
        exact_rational(coordinate.root()),
    };
}

Integer greatest_common_divisor(
    Integer left,
    Integer right)
{
    left = left < 0 ? -left : left;
    right = right < 0 ? -right : right;
    while (right != 0) {
        const Integer remainder = left % right;
        left = right;
        right = remainder;
    }
    return left;
}

Integer least_common_multiple(
    const Integer& left,
    const Integer& right)
{
    return left / greatest_common_divisor(left, right)
        * right;
}

Polynomial1 polynomial_1(
    const std::vector<std::string>& coefficients)
{
    std::vector<Integer> parsed;
    parsed.reserve(coefficients.size());
    for (const std::string& coefficient : coefficients) {
        parsed.emplace_back(coefficient);
    }
    return typename CGAL::Polynomial_traits_d<
        Polynomial1>::Construct_polynomial()(
        parsed.begin(),
        parsed.end());
}

Polynomial2 variable(std::size_t index)
{
    return CGAL::shift(
        Polynomial2(Integer(1)),
        1,
        index);
}

Polynomial2 polynomial_in_x(
    const std::vector<std::string>& coefficients)
{
    const Polynomial2 x = variable(0);
    Polynomial2 result(Integer(0));
    Polynomial2 power(Integer(1));
    for (const std::string& coefficient : coefficients) {
        result += Integer(coefficient) * power;
        power *= x;
    }
    return result;
}

Polynomial2 polynomial_from_rows(
    const std::vector<std::vector<std::string>>& rows)
{
    const Polynomial2 x = variable(0);
    const Polynomial2 y = variable(1);
    Polynomial2 result(Integer(0));
    Polynomial2 x_power(Integer(1));
    for (const auto& row : rows) {
        Polynomial2 y_power(Integer(1));
        for (const std::string& coefficient : row) {
            result += Integer(coefficient)
                * x_power * y_power;
            y_power *= y;
        }
        x_power *= x;
    }
    return result;
}

std::vector<std::string> primitive_coefficients(
    const Polynomial1& polynomial)
{
    const Polynomial1 canonical =
        typename CGAL::Polynomial_traits_d<
            Polynomial1>::Canonicalize()(
            polynomial);
    const int degree = CGAL::degree(canonical);
    std::vector<std::string> result;
    result.reserve(
        static_cast<std::size_t>(degree + 1));
    for (int index = 0; index <= degree; ++index) {
        result.push_back(
            canonical[index].convert_to<std::string>());
    }
    return result;
}

AlgebraicReal1 reconstruct_root(
    const AlgebraicRootRecord2& record)
{
    Kernel1 kernel;
    const Polynomial1 factor =
        polynomial_1(record.factor_coefficients);
    std::vector<AlgebraicReal1> roots;
    kernel.solve_1_object()(
        factor,
        true,
        std::back_inserter(roots));
    if (record.root_ordinal >= roots.size()) {
        throw IncompleteSegmentPartitionError(
            "segment fibre root ordinal is out of range");
    }
    if (algebraic_root_id_v1(
            [&record]() {
                std::vector<Integer> values;
                values.reserve(
                    record.factor_coefficients.size());
                for (const std::string& coefficient :
                     record.factor_coefficients) {
                    values.emplace_back(coefficient);
                }
                return values;
            }(),
            record.root_ordinal)
        != record.root_id) {
        throw IncompleteSegmentPartitionError(
            "segment fibre root identity is not canonical");
    }
    return roots[record.root_ordinal];
}

std::vector<std::string> pullback_fields(
    const ProjectionRecord2& pullback)
{
    const std::vector<std::string> fields =
        decode_string_sequence(
            pullback.projection_id);
    if (fields.size() != 3
        || fields.front()
            != "segment-pullback-v1") {
        throw IncompleteSegmentPartitionError(
            "segment pullback identity is malformed");
    }
    return fields;
}

bool same_event_identity(
    const PartitionEvent2& left,
    const PartitionEvent2& right)
{
    return std::tie(
               left.kind,
               left.feature_id,
               left.support_id,
               left.trim_id,
               left.vertex_id,
               left.branch_id,
               left.first_feature_id,
               left.second_feature_id,
               left.first_chart_id,
               left.second_chart_id)
        == std::tie(
               right.kind,
               right.feature_id,
               right.support_id,
               right.trim_id,
               right.vertex_id,
               right.branch_id,
               right.first_feature_id,
               right.second_feature_id,
               right.first_chart_id,
               right.second_chart_id);
}

const SegmentBranchState2* branch_for_id(
    const std::vector<SegmentBranchState2>& branches,
    const std::string& branch_id)
{
    const auto found = std::find_if(
        branches.begin(),
        branches.end(),
        [&branch_id](const SegmentBranchState2& branch) {
            return branch.branch.branch_id == branch_id;
        });
    return found == branches.end() ? nullptr : &*found;
}

const BranchPairDisposition2* event_pair(
    const PartitionEvent2& event,
    const std::vector<SegmentBranchState2>& branches,
    const std::vector<BranchPairDisposition2>& pairs)
{
    const auto found = std::find_if(
        pairs.begin(),
        pairs.end(),
        [&event, &branches](
            const BranchPairDisposition2& pair) {
            const SegmentBranchState2* first =
                branch_for_id(
                    branches,
                    pair.first_branch_id);
            const SegmentBranchState2* second =
                branch_for_id(
                    branches,
                    pair.second_branch_id);
            return first != nullptr
                && second != nullptr
                && first->branch.feature_id
                    == event.first_feature_id
                && second->branch.feature_id
                    == event.second_feature_id
                && first->branch.rim_chart_id
                    == event.first_chart_id
                && second->branch.rim_chart_id
                    == event.second_chart_id;
        });
    return found == pairs.end() ? nullptr : &*found;
}

void bind_pair(
    PartitionEvent2& event,
    const BranchPairDisposition2& pair)
{
    event.pair_sheet_id = pair.pair_sheet_id;
    event.first_branch_id = pair.first_branch_id;
    event.second_branch_id = pair.second_branch_id;
    event.branch_id = pair.first_branch_id;
}

const BoundaryFeatureRecord2& record_for_feature(
    const std::vector<BoundaryFeatureRecord2>& records,
    const std::string& feature_id)
{
    const auto found = std::find_if(
        records.begin(),
        records.end(),
        [&feature_id](
            const BoundaryFeatureRecord2& record) {
            return record.feature_id == feature_id;
        });
    if (found == records.end()) {
        throw IncompleteSegmentPartitionError(
            "segment event references an unknown feature");
    }
    return *found;
}

PointNumerators
point_numerators(
    const SegmentEventSource2& source,
    const std::string& chart)
{
    const Rational x0 =
        parse_rational(source.x0().text(), "segment x0");
    const Rational y0 =
        parse_rational(source.y0().text(), "segment y0");
    const Rational dx =
        parse_rational(source.x1().text(), "segment x1")
        - x0;
    const Rational dy =
        parse_rational(source.y1().text(), "segment y1")
        - y0;
    const Rational radius =
        parse_rational(
            source.tool_radius().text(),
            "tool radius");
    const Rational chart_sign =
        chart == "rim-half-0-v1"
        ? Rational(1)
        : Rational(-1);
    if (chart != "rim-half-0-v1"
        && chart != "rim-half-1-v1") {
        throw ChartCoverageError(
            "segment fibre uses an unknown rim chart");
    }
    const Polynomial2 x = variable(0);
    const Polynomial2 y = variable(1);
    const Polynomial2 unit_denominator =
        Polynomial2(Integer(1)) + y * y;
    Integer coefficient_denominator = 1;
    for (const Rational& value :
         {x0, y0, dx, dy, radius}) {
        coefficient_denominator =
            least_common_multiple(
                coefficient_denominator,
                CORE::denominator(value));
    }
    const auto scaled =
        [&coefficient_denominator](
            const Rational& value) {
            return CORE::numerator(
                value
                * Rational(coefficient_denominator));
        };
    const Polynomial2 center_x =
        scaled(x0) + scaled(dx) * x;
    const Polynomial2 center_y =
        scaled(y0) + scaled(dy) * x;
    return {
        center_x * unit_denominator
            + scaled(radius * chart_sign)
                * (
                    Polynomial2(Integer(1))
                    - y * y),
        center_y * unit_denominator
            + scaled(radius * chart_sign)
                * Integer(2) * y,
        coefficient_denominator
            * unit_denominator,
    };
}

Polynomial2 rational_polynomial(
    const std::vector<
        std::pair<Polynomial2, Rational>>& terms)
{
    Integer denominator_lcm = 1;
    for (const auto& [polynomial, coefficient] : terms) {
        static_cast<void>(polynomial);
        denominator_lcm = least_common_multiple(
            denominator_lcm,
            CORE::denominator(coefficient));
    }
    Polynomial2 result(Integer(0));
    for (const auto& [polynomial, coefficient] : terms) {
        result += CORE::numerator(
                      coefficient
                      * Rational(denominator_lcm))
            * polynomial;
    }
    return result;
}

RadicalPredicate trim_predicate(
    const Polynomial2& point_x,
    const Polynomial2& point_y,
    const Polynomial2& point_denominator,
    const GpsPoint& endpoint,
    const Rational& tangent_x,
    const Rational& tangent_y,
    bool source_bound)
{
    const CoordinateParts x =
        coordinate_parts(endpoint.x());
    const CoordinateParts y_coordinate =
        coordinate_parts(endpoint.y());
    Rational radicand = 0;
    if (x.radical_coefficient != 0) {
        radicand = x.radicand;
    }
    if (y_coordinate.radical_coefficient != 0) {
        if (radicand != 0
            && radicand
                != y_coordinate.radicand) {
            throw UnsupportedAlgebraicVertexProjectionError(
                "trim endpoint coordinates use distinct radicals");
        }
        radicand = y_coordinate.radicand;
    }
    const Rational direction =
        source_bound ? Rational(1) : Rational(-1);
    const Polynomial2 rational_part =
        rational_polynomial(
            {
                {
                    point_x,
                    direction * tangent_x,
                },
                {
                    point_y,
                    direction * tangent_y,
                },
                {
                    point_denominator,
                    -direction
                        * (
                            tangent_x * x.base
                            + tangent_y
                                * y_coordinate.base),
                },
            });
    const Polynomial2 radical_part =
        rational_polynomial(
            {
                {
                    point_denominator,
                    -direction
                        * (
                            tangent_x
                                * x.radical_coefficient
                            + tangent_y
                                * y_coordinate
                                      .radical_coefficient),
                },
            });
    return {
        rational_part,
        radical_part,
        radicand,
    };
}

CGAL::Sign radical_sign(
    const RadicalPredicate& predicate,
    const AlgebraicReal2& point)
{
    const CGAL::Sign rational_sign =
        exact_polynomial_sign_at_2(
            predicate.rational_part,
            point);
    const CGAL::Sign radical_sign =
        exact_polynomial_sign_at_2(
            predicate.radical_part,
            point);
    if (predicate.radicand == 0
        || radical_sign == CGAL::ZERO) {
        return rational_sign;
    }
    if (rational_sign == CGAL::ZERO) {
        return radical_sign;
    }
    if (rational_sign == radical_sign) {
        return rational_sign;
    }
    const Polynomial2 comparison =
        CORE::denominator(predicate.radicand)
            * predicate.rational_part
            * predicate.rational_part
        - CORE::numerator(predicate.radicand)
            * predicate.radical_part
            * predicate.radical_part;
    const CGAL::Sign magnitude =
        exact_polynomial_sign_at_2(
            comparison,
            point);
    if (magnitude == CGAL::ZERO) {
        return CGAL::ZERO;
    }
    return magnitude == CGAL::POSITIVE
        ? rational_sign
        : radical_sign;
}

bool line_trim_accepts(
    const BoundaryFeatureRecord2& record,
    const SegmentEventSource2& source,
    const std::string& chart,
    const AlgebraicReal2& point)
{
    const Rational a = parse_rational(
        record.primitive_coefficients[0],
        "line support");
    const Rational b = parse_rational(
        record.primitive_coefficients[1],
        "line support");
    Rational tangent_x = -b;
    Rational tangent_y = a;
    const GpsPoint::CoordNT source_lambda =
        GpsPoint::CoordNT(
            Epeck::FT(tangent_x))
            * record.curve.source().x()
        + GpsPoint::CoordNT(
            Epeck::FT(tangent_y))
            * record.curve.source().y();
    const GpsPoint::CoordNT target_lambda =
        GpsPoint::CoordNT(
            Epeck::FT(tangent_x))
            * record.curve.target().x()
        + GpsPoint::CoordNT(
            Epeck::FT(tangent_y))
            * record.curve.target().y();
    if (CGAL::compare(
            source_lambda,
            target_lambda)
        == CGAL::LARGER) {
        tangent_x = -tangent_x;
        tangent_y = -tangent_y;
    }
    const PointNumerators numerators =
        point_numerators(source, chart);
    const RadicalPredicate source_predicate =
        trim_predicate(
            numerators.x,
            numerators.y,
            numerators.denominator,
            record.curve.source(),
            tangent_x,
            tangent_y,
            true);
    const RadicalPredicate target_predicate =
        trim_predicate(
            numerators.x,
            numerators.y,
            numerators.denominator,
            record.curve.target(),
            tangent_x,
            tangent_y,
            false);
    return radical_sign(source_predicate, point)
            != CGAL::NEGATIVE
        && radical_sign(target_predicate, point)
            != CGAL::NEGATIVE;
}

bool circle_trim_accepts(
    const BoundaryFeatureRecord2& record,
    const SegmentEventSource2& source,
    const std::string& chart,
    const AlgebraicReal2& point)
{
    const PointNumerators numerators =
        point_numerators(source, chart);
    const RadicalPredicate left_predicate =
        trim_predicate(
            numerators.x,
            numerators.y,
            numerators.denominator,
            record.curve.left(),
            Rational(1),
            Rational(0),
            true);
    const RadicalPredicate right_predicate =
        trim_predicate(
            numerators.x,
            numerators.y,
            numerators.denominator,
            record.curve.right(),
            Rational(1),
            Rational(0),
            false);
    if (radical_sign(left_predicate, point)
            == CGAL::NEGATIVE
        || radical_sign(right_predicate, point)
            == CGAL::NEGATIVE) {
        return false;
    }
    const Rational quadratic = parse_rational(
        record.primitive_coefficients[0],
        "circle support");
    const Rational center_y =
        -parse_rational(
            record.primitive_coefficients[2],
            "circle support")
        / (Rational(2) * quadratic);
    const bool upper =
        (
            record.curve.orientation()
                == CGAL::COUNTERCLOCKWISE
            && !record.curve.is_directed_right())
        || (
            record.curve.orientation()
                == CGAL::CLOCKWISE
            && record.curve.is_directed_right());
    const Polynomial2 hemisphere =
        upper
        ? rational_polynomial(
              {
                  {
                      numerators.y,
                      Rational(1),
                  },
                  {
                      numerators.denominator,
                      -center_y,
                  },
              })
        : rational_polynomial(
              {
                  {
                      numerators.denominator,
                      center_y,
                  },
                  {
                      numerators.y,
                      Rational(-1),
                  },
              });
    return exact_polynomial_sign_at_2(
               hemisphere,
               point)
        != CGAL::NEGATIVE;
}

std::string exact_vertex_id(
    const BoundaryFeatureRecord2& record,
    const SegmentEventSource2& source,
    const std::string& chart,
    const AlgebraicReal2& point)
{
    const PointNumerators numerators =
        point_numerators(source, chart);
    const auto equals =
        [&numerators, &point](
            const GpsPoint& endpoint) {
            const RadicalPredicate x_difference =
                trim_predicate(
                    numerators.x,
                    numerators.y,
                    numerators.denominator,
                    endpoint,
                    Rational(1),
                    Rational(0),
                    true);
            const RadicalPredicate y_difference =
                trim_predicate(
                    numerators.x,
                    numerators.y,
                    numerators.denominator,
                    endpoint,
                    Rational(0),
                    Rational(1),
                    true);
            return radical_sign(
                       x_difference,
                       point)
                    == CGAL::ZERO
                && radical_sign(
                       y_difference,
                       point)
                    == CGAL::ZERO;
        };
    if (equals(record.curve.source())) {
        return record.source_vertex_id;
    }
    if (equals(record.curve.target())) {
        return record.target_vertex_id;
    }
    return {};
}

bool chart_accepts(
    const AlgebraicReal1& parameter)
{
    Kernel1 kernel;
    return kernel.compare_1_object()(
               parameter,
               Rational(-1))
            != CGAL::SMALLER
        && kernel.compare_1_object()(
               parameter,
               Rational(1))
            != CGAL::LARGER;
}

struct ExactBranches {
    std::vector<SegmentBranchState2> finite;
    std::vector<std::string> active_ids;
};

ExactBranches exact_fibre_branches(
    const std::vector<BoundaryFeatureRecord2>& records,
    const SegmentEventSource2& source,
    const std::vector<ProjectionRecord2>& pullbacks,
    const AlgebraicRootRecord2& root_record,
    const AlgebraicReal1& root,
    const std::vector<SegmentBranchState2>& side_states,
    const std::vector<PartitionEvent2>& events)
{
    const Polynomial2 root_polynomial =
        polynomial_in_x(
            root_record.factor_coefficients);
    Kernel1 kernel1;
    Kernel2 kernel2;
    ExactBranches result;
    for (const ProjectionRecord2& pullback :
         pullbacks) {
        const std::vector<std::string> fields =
            pullback_fields(pullback);
        const BoundaryFeatureRecord2& record =
            record_for_feature(records, fields[1]);
        const Polynomial2 support =
            polynomial_from_rows(
                pullback.coefficient_rows);
        if (!kernel2.is_coprime_2_object()(
                root_polynomial,
                support)) {
            const bool overlap_rechecked =
                std::any_of(
                    events.begin(),
                    events.end(),
                    [&record](
                        const PartitionEvent2& event) {
                        return event.kind
                                == "support-overlap"
                            && event.feature_id
                                == record.feature_id;
                    });
            if (!overlap_rechecked) {
                throw IncompleteSegmentPartitionError(
                    "unrecorded positive-dimensional fibre component");
            }
            for (const SegmentBranchState2& state :
                 side_states) {
                if (state.branch.feature_id
                        == record.feature_id
                    && std::find(
                           result.active_ids.begin(),
                           result.active_ids.end(),
                           state.branch.branch_id)
                        == result.active_ids.end()) {
                    result.active_ids.push_back(
                        state.branch.branch_id);
                }
            }
            continue;
        }
        std::vector<
            std::pair<
                AlgebraicReal2,
                Kernel2::Multiplicity_type>>
            solutions;
        kernel2.solve_2_object()(
            root_polynomial,
            support,
            std::back_inserter(solutions));
        struct AcceptedParameter {
            AlgebraicReal1 parameter;
            std::string vertex_id;
        };
        std::vector<AcceptedParameter> parameters;
        for (const auto& [point, multiplicity] :
             solutions) {
            static_cast<void>(multiplicity);
            if (kernel1.compare_1_object()(
                    kernel2.compute_x_2_object()(
                        point),
                    root)
                    != CGAL::EQUAL
                || exact_polynomial_sign_at_2(
                       support,
                       point)
                    != CGAL::ZERO) {
                continue;
            }
            const AlgebraicReal1 parameter =
                kernel2.compute_y_2_object()(point);
            if (!chart_accepts(parameter)
                || !(
                    record.support_kind == "line"
                    ? line_trim_accepts(
                          record,
                          source,
                          fields[2],
                          point)
                    : circle_trim_accepts(
                          record,
                          source,
                          fields[2],
                          point))) {
                continue;
            }
            parameters.push_back(
                {
                    parameter,
                    exact_vertex_id(
                        record,
                        source,
                        fields[2],
                        point),
                });
        }
        std::sort(
            parameters.begin(),
            parameters.end(),
            [&kernel1](
                const AcceptedParameter& left,
                const AcceptedParameter& right) {
                return kernel1.compare_1_object()(
                           left.parameter,
                           right.parameter)
                    == CGAL::SMALLER;
            });
        std::vector<SegmentBoundaryBranch2>
            identities;
        for (const SegmentBranchState2& state :
             side_states) {
            if (state.branch.feature_id
                    == record.feature_id
                && state.branch.rim_chart_id
                    == fields[2]) {
                identities.push_back(state.branch);
            }
        }
        std::sort(
            identities.begin(),
            identities.end(),
            [](const SegmentBoundaryBranch2& left,
               const SegmentBoundaryBranch2& right) {
                return left.rim_sheet_ordinal
                    < right.rim_sheet_ordinal;
            });
        if (parameters.empty() || identities.empty()) {
            continue;
        }
        if (parameters.size() != identities.size()
            && parameters.size() != 1) {
            throw IncompleteSegmentPartitionError(
                "exact fibre branch multiplicity does not match adjacent sheets");
        }
        for (std::size_t index = 0;
             index < identities.size();
             ++index) {
            const AcceptedParameter& accepted =
                parameters.size() == 1
                ? parameters.front()
                : parameters[index];
            const AlgebraicReal1& parameter =
                accepted.parameter;
            const Polynomial1 factor =
                kernel1.compute_polynomial_1_object()(
                    parameter);
            SegmentBoundaryBranch2 identity =
                identities[index];
            identity.vertex_id =
                accepted.vertex_id;
            result.finite.push_back(
                {
                    std::move(identity),
                    parameter,
                    factor,
                    index,
                });
        }
    }
    std::sort(
        result.finite.begin(),
        result.finite.end(),
        [&kernel1](
            const SegmentBranchState2& left,
            const SegmentBranchState2& right) {
            if (left.branch.rim_chart_id
                != right.branch.rim_chart_id) {
                return left.branch.rim_chart_id
                    < right.branch.rim_chart_id;
            }
            const CGAL::Comparison_result comparison =
                kernel1.compare_1_object()(
                    left.rim_parameter,
                    right.rim_parameter);
            if (comparison != CGAL::EQUAL) {
                return comparison == CGAL::SMALLER;
            }
            return left.branch.rim_sheet_ordinal
                < right.branch.rim_sheet_ordinal;
        });
    for (const SegmentBranchState2& state :
         result.finite) {
        if (std::find(
                result.active_ids.begin(),
                result.active_ids.end(),
                state.branch.branch_id)
            == result.active_ids.end()) {
            result.active_ids.push_back(
                state.branch.branch_id);
        }
    }
    return result;
}

Polynomial1 derivative(
    const std::vector<std::string>& coefficients)
{
    std::vector<Integer> values;
    if (coefficients.size() < 2) {
        return Polynomial1(Integer(0));
    }
    values.reserve(coefficients.size() - 1);
    for (std::size_t index = 1;
         index < coefficients.size();
         ++index) {
        values.push_back(
            Integer(coefficients[index])
            * Integer(index));
    }
    return typename CGAL::Polynomial_traits_d<
        Polynomial1>::Construct_polynomial()(
        values.begin(),
        values.end());
}

bool vertex_original_holds(
    const BoundaryFeatureRecord2& record,
    const SegmentEventSource2& source,
    const PartitionEvent2& event,
    const AlgebraicReal1& root)
{
    const GpsPoint* vertex = nullptr;
    if (event.vertex_id == record.source_vertex_id) {
        vertex = &record.curve.source();
    } else if (
        event.vertex_id == record.target_vertex_id) {
        vertex = &record.curve.target();
    }
    if (vertex == nullptr) {
        return false;
    }
    const CoordinateParts x =
        coordinate_parts(vertex->x());
    const CoordinateParts y =
        coordinate_parts(vertex->y());
    Rational radicand = 0;
    if (x.radical_coefficient != 0) {
        radicand = x.radicand;
    }
    if (y.radical_coefficient != 0) {
        if (radicand != 0
            && radicand != y.radicand) {
            throw UnsupportedAlgebraicVertexProjectionError(
                "vertex replay uses distinct radicals");
        }
        radicand = y.radicand;
    }
    const Rational x0 =
        parse_rational(source.x0().text(), "segment x0");
    const Rational y0 =
        parse_rational(source.y0().text(), "segment y0");
    const Rational dx =
        parse_rational(source.x1().text(), "segment x1")
        - x0;
    const Rational dy =
        parse_rational(source.y1().text(), "segment y1")
        - y0;
    const Rational radius =
        parse_rational(
            source.tool_radius().text(),
            "tool radius");
    using RationalPolynomial =
        CGAL::Polynomial<Rational>;
    const RationalPolynomial center_x(
        {x0 - x.base, dx});
    const RationalPolynomial center_y(
        {y0 - y.base, dy});
    const RationalPolynomial rational_part =
        center_x * center_x
        + center_y * center_y
        + RationalPolynomial(
            radicand
                    * (
                        x.radical_coefficient
                            * x.radical_coefficient
                        + y.radical_coefficient
                            * y.radical_coefficient)
                - radius * radius);
    const RationalPolynomial radical_part =
        Rational(-2)
        * (
            x.radical_coefficient * center_x
            + y.radical_coefficient * center_y);
    const auto integer_polynomial =
        [](const RationalPolynomial& polynomial) {
            using FractionTraits =
                CGAL::Fraction_traits<
                    RationalPolynomial>;
            typename FractionTraits::Numerator_type
                numerator;
            typename FractionTraits::Denominator_type
                denominator;
            typename FractionTraits::Decompose()(
                polynomial,
                numerator,
                denominator);
            if (denominator <= 0) {
                throw IncompleteSegmentPartitionError(
                    "vertex replay denominator is not positive");
            }
            return numerator;
        };
    const Polynomial1 a =
        integer_polynomial(rational_part);
    const Polynomial1 b =
        integer_polynomial(radical_part);
    Kernel1 kernel;
    const CGAL::Sign a_sign =
        kernel.sign_at_1_object()(a, root);
    const CGAL::Sign b_sign =
        CGAL::is_zero(b)
        ? CGAL::ZERO
        : kernel.sign_at_1_object()(b, root);
    if (radicand == 0 || b_sign == CGAL::ZERO) {
        return a_sign == CGAL::ZERO;
    }
    if (a_sign == CGAL::ZERO
        || a_sign == b_sign) {
        return false;
    }
    const Polynomial1 resultant =
        integer_polynomial(
            rational_part * rational_part
            - radicand
                * radical_part * radical_part);
    return kernel.sign_at_1_object()(
               resultant,
               root)
        == CGAL::ZERO;
}

bool support_overlap_holds(
    const BoundaryFeatureRecord2& record,
    const SegmentEventSource2& source,
    const AlgebraicReal1& root)
{
    if (record.support_kind != "circle") {
        return false;
    }
    const Rational quadratic = parse_rational(
        record.primitive_coefficients[0],
        "circle support");
    const Rational center_x =
        -parse_rational(
            record.primitive_coefficients[1],
            "circle support")
        / (Rational(2) * quadratic);
    const Rational center_y =
        -parse_rational(
            record.primitive_coefficients[2],
            "circle support")
        / (Rational(2) * quadratic);
    const Rational constant = parse_rational(
        record.primitive_coefficients[3],
        "circle support");
    const Rational radius_squared =
        center_x * center_x
        + center_y * center_y
        - constant / quadratic;
    const Rational tool_radius = parse_rational(
        source.tool_radius().text(),
        "tool radius");
    if (radius_squared
        != tool_radius * tool_radius) {
        return false;
    }
    const std::vector<std::string> x_coefficients{
        source.x0().text(),
        (
            parse_rational(source.x1().text(), "segment x1")
            - parse_rational(source.x0().text(), "segment x0"))
            .convert_to<std::string>(),
    };
    const std::vector<std::string> y_coefficients{
        source.y0().text(),
        (
            parse_rational(source.y1().text(), "segment y1")
            - parse_rational(source.y0().text(), "segment y0"))
            .convert_to<std::string>(),
    };
    const auto shifted =
        [](const std::vector<std::string>& values,
           const Rational& center) {
            std::vector<Rational> coefficients{
                parse_rational(values[0], "segment coordinate")
                    - center,
                parse_rational(values[1], "segment coordinate"),
            };
            using RationalPolynomial =
                CGAL::Polynomial<Rational>;
            using FractionTraits =
                CGAL::Fraction_traits<
                    RationalPolynomial>;
            const RationalPolynomial polynomial(
                coefficients.begin(),
                coefficients.end());
            typename FractionTraits::Numerator_type numerator;
            typename FractionTraits::Denominator_type denominator;
            typename FractionTraits::Decompose()(
                polynomial,
                numerator,
                denominator);
            return numerator;
        };
    Kernel1 kernel;
    return kernel.sign_at_1_object()(
               shifted(x_coefficients, center_x),
               root)
            == CGAL::ZERO
        && kernel.sign_at_1_object()(
               shifted(y_coefficients, center_y),
               root)
            == CGAL::ZERO;
}

std::vector<PartitionEvent2> replay_events(
    const std::vector<BoundaryFeatureRecord2>& records,
    const SegmentEventSource2& source,
    const std::vector<ProjectionInput2>& projections,
    const AlgebraicReal1& root,
    const std::vector<SegmentBranchState2>& branches,
    const std::vector<BranchPairDisposition2>& pairs,
    const std::vector<PartitionEvent2>& events,
    std::size_t left_count,
    std::size_t right_count,
    bool& candidates_rechecked)
{
    Kernel1 kernel;
    candidates_rechecked = true;
    std::vector<PartitionEvent2> result;
    result.reserve(events.size());
    for (PartitionEvent2 event : events) {
        const BoundaryFeatureRecord2& record =
            record_for_feature(records, event.feature_id);
        bool projection_rechecked = false;
        CGAL::Sign crossing_sign = CGAL::ZERO;
        for (const ProjectionInput2& projection :
             projections) {
            const bool owns_event = std::any_of(
                projection.events.begin(),
                projection.events.end(),
                [&event](const PartitionEvent2& candidate) {
                    return same_event_identity(
                        event,
                        candidate);
                });
            if (!owns_event) {
                continue;
            }
            const Polynomial1 polynomial =
                polynomial_1(projection.coefficients);
            if (kernel.sign_at_1_object()(
                    polynomial,
                    root)
                != CGAL::ZERO) {
                continue;
            }
            projection_rechecked = true;
            const std::vector<std::string> fields =
                decode_string_sequence(
                    projection.projection_id);
            if (!fields.empty()
                && (
                    fields.front()
                        == "segment-endpoint-order-v1"
                    || fields.front()
                        == "segment-vertex-passage-v1")) {
                const Polynomial1 slope =
                    derivative(
                        projection.coefficients);
                if (!CGAL::is_zero(slope)) {
                    crossing_sign =
                        kernel.sign_at_1_object()(
                            slope,
                            root);
                }
            }
        }
        bool original_rechecked =
            projection_rechecked;
        bool orientation_rechecked =
            original_rechecked;
        bool trim_rechecked = std::any_of(
            branches.begin(),
            branches.end(),
            [&event](const SegmentBranchState2& branch) {
                return branch.branch.feature_id
                    == event.feature_id;
            });
        if (event.kind == "endpoint-order"
            && !event.first_feature_id.empty()) {
            const BranchPairDisposition2* pair =
                event_pair(event, branches, pairs);
            const SegmentBranchState2* first =
                pair == nullptr
                ? nullptr
                : branch_for_id(
                      branches,
                      pair->first_branch_id);
            const SegmentBranchState2* second =
                pair == nullptr
                ? nullptr
                : branch_for_id(
                      branches,
                      pair->second_branch_id);
            const bool common_vertex =
                first != nullptr
                && second != nullptr
                && !first->branch.vertex_id.empty()
                && first->branch.vertex_id
                    == second->branch.vertex_id
                && pair->orientation_disposition
                    == "ccw-zero";
            if (common_vertex) {
                bind_pair(event, *pair);
                event.vertex_id =
                    first->branch.vertex_id;
            }
            original_rechecked =
                projection_rechecked
                && common_vertex;
            trim_rechecked = original_rechecked;
            event.incidence_permutation_rechecked =
                original_rechecked
                && crossing_sign != CGAL::ZERO;
            orientation_rechecked =
                event.incidence_permutation_rechecked;
            if (crossing_sign == CGAL::NEGATIVE) {
                event.disposition = "split";
            } else if (
                crossing_sign == CGAL::POSITIVE) {
                event.disposition = "merge";
            }
        } else if (event.kind == "endpoint-order") {
            const auto exact_vertex_branch =
                std::find_if(
                    branches.begin(),
                    branches.end(),
                    [&event](
                        const SegmentBranchState2& branch) {
                        return branch.branch.feature_id
                                == event.feature_id
                            && branch.branch.vertex_id
                                == event.vertex_id;
                    });
            const bool has_exact_vertex_branch =
                exact_vertex_branch != branches.end();
            if (has_exact_vertex_branch) {
                event.branch_id =
                    exact_vertex_branch
                        ->branch.branch_id;
            }
            const bool original_vertex =
                vertex_original_holds(
                    record,
                    source,
                    event,
                    root);
            if (!has_exact_vertex_branch
                && original_vertex) {
                const auto incident_branch =
                    std::find_if(
                        branches.begin(),
                        branches.end(),
                        [&event](
                            const SegmentBranchState2& branch) {
                            return branch.branch.feature_id
                                    == event.feature_id
                                && branch.branch.support_id
                                    == event.support_id;
                        });
                if (incident_branch != branches.end()) {
                    event.branch_id =
                        incident_branch
                            ->branch.branch_id;
                }
            }
            original_rechecked =
                projection_rechecked
                && (
                    has_exact_vertex_branch
                    || original_vertex);
            trim_rechecked = original_rechecked;
            event.incidence_permutation_rechecked =
                crossing_sign != CGAL::ZERO;
            orientation_rechecked =
                original_rechecked
                && event
                       .incidence_permutation_rechecked;
            if (crossing_sign == CGAL::NEGATIVE) {
                event.disposition = "split";
            } else if (
                crossing_sign == CGAL::POSITIVE) {
                event.disposition = "merge";
            }
        } else if (event.kind == "tangent") {
            event.disposition =
                right_count > left_count
                ? "birth"
                : "death";
            const std::size_t feature_branch_count =
                static_cast<std::size_t>(
                    std::count_if(
                        branches.begin(),
                        branches.end(),
                        [&event](
                            const SegmentBranchState2& branch) {
                            return branch.branch.feature_id
                                == event.feature_id;
                        }));
            orientation_rechecked =
                original_rechecked
                && (
                    feature_branch_count == 1
                    || std::any_of(
                        pairs.begin(),
                        pairs.end(),
                        [](const BranchPairDisposition2& pair) {
                            return pair
                                       .orientation_disposition
                                == "ccw-zero";
                        }));
        } else if (event.kind == "cap-crossing") {
            if (!event.first_feature_id.empty()) {
                const BranchPairDisposition2* pair =
                    event_pair(event, branches, pairs);
                const bool own_pair_equal =
                    pair != nullptr
                    && pair->cap_disposition
                        == "equal-cap";
                if (own_pair_equal) {
                    bind_pair(event, *pair);
                }
                original_rechecked =
                    projection_rechecked
                    && own_pair_equal;
                trim_rechecked = original_rechecked;
                orientation_rechecked =
                    original_rechecked;
            } else {
                orientation_rechecked =
                    original_rechecked
                    && std::any_of(
                        pairs.begin(),
                        pairs.end(),
                        [](const BranchPairDisposition2& pair) {
                            return pair.cap_disposition
                                == "equal-cap";
                        });
            }
        } else if (
            event.kind == "support-overlap") {
            original_rechecked =
                projection_rechecked
                && support_overlap_holds(
                    record,
                    source,
                    root);
            trim_rechecked = original_rechecked;
            orientation_rechecked =
                original_rechecked;
        }
        event.left_active_count = left_count;
        event.right_active_count = right_count;
        event.original_equations_rechecked =
            original_rechecked;
        event.orientation_rechecked =
            orientation_rechecked;
        event.trim_disposition =
            trim_rechecked ? "accepted" : "rejected";
        candidates_rechecked =
            candidates_rechecked
            && projection_rechecked;
        if (original_rechecked
            && orientation_rechecked
            && trim_rechecked) {
            result.push_back(std::move(event));
        }
    }
    return result;
}

} // namespace

SegmentFibreEvaluation2 evaluate_segment_fibre(
    const std::vector<BoundaryFeatureRecord2>& records,
    const SegmentEventSource2& source,
    const std::vector<ProjectionRecord2>& pullbacks,
    const std::vector<ProjectionInput2>& event_projections,
    const AlgebraicRootRecord2& root_record,
    const std::vector<PartitionEvent2>& events,
    const std::vector<SegmentBranchState2>& left_states,
    const std::vector<SegmentBranchState2>& right_states)
{
    const AlgebraicReal1 root =
        reconstruct_root(root_record);
    std::vector<SegmentBranchState2> side_states =
        left_states.empty()
        ? right_states
        : left_states;
    for (const SegmentBranchState2& state :
         right_states) {
        const auto found = std::find_if(
            side_states.begin(),
            side_states.end(),
            [&state](const SegmentBranchState2& candidate) {
                return candidate.branch.branch_id
                    == state.branch.branch_id;
            });
        if (found == side_states.end()) {
            side_states.push_back(state);
        }
    }
    ExactBranches branches =
        exact_fibre_branches(
            records,
            source,
            pullbacks,
            root_record,
            root,
            side_states,
            events);
    std::vector<BranchPairDisposition2> pairs;
    if (!branches.finite.empty()) {
        pairs = segment_branch_pair_dispositions(
            branches.finite,
            source);
    }
    bool candidates_rechecked = false;
    std::vector<PartitionEvent2> replayed =
        replay_events(
            records,
            source,
            event_projections,
            root,
            branches.finite,
            pairs,
            events,
            left_states.size(),
            right_states.size(),
            candidates_rechecked);
    const bool original_rechecked =
        (
            replayed.empty()
            && candidates_rechecked)
        || (
            !replayed.empty()
            && std::all_of(
            replayed.begin(),
            replayed.end(),
            [](const PartitionEvent2& event) {
                return event.original_equations_rechecked;
            }));
    const bool trim_rechecked =
        replayed.empty()
        || std::all_of(
            replayed.begin(),
            replayed.end(),
            [](const PartitionEvent2& event) {
                return event.trim_disposition
                    == "accepted";
            });
    const bool orientation_rechecked =
        replayed.empty()
        || std::all_of(
            replayed.begin(),
            replayed.end(),
            [](const PartitionEvent2& event) {
                return event.orientation_rechecked;
            });
    return {
        std::move(branches.finite),
        std::move(branches.active_ids),
        std::move(replayed),
        std::move(pairs),
        true,
        original_rechecked,
        orientation_rechecked,
        trim_rechecked,
    };
}
