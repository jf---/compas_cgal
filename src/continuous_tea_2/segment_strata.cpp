#include "segment_strata.h"

#include "exact_polynomial_sign_2.h"
#include "event_certificate.h"

#include <algorithm>
#include <iterator>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <variant>

#include <CGAL/CORE/BigRat.h>
#include <CGAL/Fraction_traits.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/number_utils.h>

namespace {

using Integer = CORE::BigInt;
using Rational = CORE::BigRat;
using Kernel1 = ExactAlgebraicKernel1;
using Kernel2 = ExactAlgebraicKernel2;
using Polynomial1 = Kernel1::Polynomial_1;
using Polynomial2 = Kernel2::Polynomial_2;
using AlgebraicReal1 = Kernel1::Algebraic_real_1;
using AlgebraicReal2 = Kernel2::Algebraic_real_2;
using IntersectionPoint =
    std::pair<GpsPoint, unsigned int>;
using IntersectionResult =
    std::variant<IntersectionPoint, GpsXCurve>;

struct RimRoot {
    AlgebraicReal1 value;
    Polynomial1 factor;
    std::vector<Integer> coefficients;
    std::size_t ordinal;
};

Rational parse_rational(const std::string& text)
{
    const std::size_t separator = text.find('/');
    try {
        if (separator == std::string::npos) {
            return Rational(Integer(text));
        }
        if (text.find('/', separator + 1)
            != std::string::npos) {
            throw IncompleteSegmentPartitionError(
                "exact cell witness is not rational");
        }
        const Integer numerator(
            text.substr(0, separator));
        const Integer denominator(
            text.substr(separator + 1));
        if (denominator == 0) {
            throw IncompleteSegmentPartitionError(
                "exact cell witness has zero denominator");
        }
        return Rational(numerator, denominator);
    } catch (const EventSubstrateError&) {
        throw;
    } catch (const std::exception&) {
        throw IncompleteSegmentPartitionError(
            "exact cell witness is not rational");
    }
}

Rational exact_rational(const Epeck::FT& value)
{
    return CGAL::exact(value);
}

Epeck::FT exact_ft(const Rational& value)
{
    static_assert(
        std::is_same_v<Rational, CGAL::Epeck_ft>);
    return Epeck::FT(value);
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

std::vector<Integer> primitive_coefficients(
    const std::vector<Rational>& values)
{
    Integer denominator_lcm = 1;
    for (const Rational& value : values) {
        denominator_lcm = least_common_multiple(
            denominator_lcm,
            CORE::denominator(value));
    }
    std::vector<Integer> result;
    result.reserve(values.size());
    Integer divisor = 0;
    for (const Rational& value : values) {
        const Integer integer = CORE::numerator(
            value * Rational(denominator_lcm));
        result.push_back(integer);
        divisor = greatest_common_divisor(
            divisor,
            integer);
    }
    for (Integer& value : result) {
        value /= divisor;
    }
    if (result.back() < 0) {
        for (Integer& value : result) {
            value = -value;
        }
    }
    return result;
}

std::optional<Rational> rational_square_root(
    const Rational& value)
{
    if (value < 0) {
        return std::nullopt;
    }
    using FractionTraits =
        CGAL::Fraction_traits<Rational>;
    Integer numerator;
    Integer denominator;
    typename FractionTraits::Decompose()(
        value,
        numerator,
        denominator);
    Integer numerator_root;
    Integer denominator_root;
    if (!CGAL::is_square(
            numerator,
            numerator_root)
        || !CGAL::is_square(
            denominator,
            denominator_root)) {
        return std::nullopt;
    }
    return typename FractionTraits::Compose()(
        numerator_root,
        denominator_root);
}

Polynomial1 polynomial_from_integers(
    const std::vector<Integer>& coefficients)
{
    return typename CGAL::Polynomial_traits_d<
        Polynomial1>::Construct_polynomial()(
        coefficients.begin(),
        coefficients.end());
}

std::vector<std::string> coefficient_text(
    const std::vector<Integer>& coefficients)
{
    std::vector<std::string> result;
    result.reserve(coefficients.size());
    for (const Integer& coefficient : coefficients) {
        result.push_back(
            coefficient.convert_to<std::string>());
    }
    return result;
}

bool same_point(
    const GpsPoint& left,
    const GpsPoint& right)
{
    return GpsTraits().compare_xy_2_object()(
               left,
               right)
        == CGAL::EQUAL;
}

GpsPolygon cutter_polygon(
    const EPoint& center,
    const Epeck::FT& radius)
{
    const ECircle circle(
        center,
        radius * radius,
        CGAL::COUNTERCLOCKWISE);
    const GpsPoint left(
        center.x() - radius,
        center.y());
    const GpsPoint right(
        center.x() + radius,
        center.y());
    GpsPolygon polygon;
    polygon.push_back(
        GpsXCurve(
            circle,
            left,
            right,
            CGAL::COUNTERCLOCKWISE));
    polygon.push_back(
        GpsXCurve(
            circle,
            right,
            left,
            CGAL::COUNTERCLOCKWISE));
    return polygon;
}

std::vector<GpsPoint> trimmed_hits(
    const BoundaryFeatureRecord2& record,
    const GpsPolygon& cutter)
{
    std::vector<GpsPoint> hits;
    const auto intersect =
        GpsTraits().intersect_2_object();
    for (auto cutter_curve = cutter.curves_begin();
         cutter_curve != cutter.curves_end();
         ++cutter_curve) {
        std::vector<IntersectionResult> intersections;
        intersect(
            record.curve,
            *cutter_curve,
            std::back_inserter(intersections));
        for (const IntersectionResult& intersection :
             intersections) {
            const auto* point =
                std::get_if<IntersectionPoint>(
                    &intersection);
            if (point == nullptr) {
                throw IncompleteSegmentPartitionError(
                    "positive-width cutter/support overlap requires an overlap stratum");
            }
            const auto found = std::find_if(
                hits.begin(),
                hits.end(),
                [&point](const GpsPoint& candidate) {
                    return same_point(
                        candidate,
                        point->first);
                });
            if (found == hits.end()) {
                hits.push_back(point->first);
            }
        }
    }
    return hits;
}

std::string chart_for_point(
    const GpsPoint& point,
    const Epeck::FT& center_x,
    const Epeck::FT& center_y)
{
    const GpsPoint::CoordNT exact_center_x(center_x);
    const CGAL::Comparison_result x_side =
        CGAL::compare(point.x(), exact_center_x);
    if (x_side == CGAL::LARGER) {
        return "rim-half-0-v1";
    }
    if (x_side == CGAL::SMALLER) {
        return "rim-half-1-v1";
    }
    const CGAL::Comparison_result y_side =
        CGAL::compare(
            point.y(),
            GpsPoint::CoordNT(center_y));
    return y_side == CGAL::SMALLER
        ? "rim-half-0-v1"
        : "rim-half-1-v1";
}

GpsPoint::CoordNT rim_parameter(
    const GpsPoint& point,
    const Epeck::FT& center_x,
    const Epeck::FT& center_y,
    const Epeck::FT& radius,
    const std::string& chart)
{
    const GpsPoint::CoordNT relative_x =
        point.x() - GpsPoint::CoordNT(center_x);
    const GpsPoint::CoordNT relative_y =
        point.y() - GpsPoint::CoordNT(center_y);
    if (chart == "rim-half-0-v1") {
        return relative_y
            / (
                GpsPoint::CoordNT(radius)
                + relative_x);
    }
    return -relative_y
        / (
            GpsPoint::CoordNT(radius)
            - relative_x);
}

RimRoot identify_rim_root(
    const GpsPoint::CoordNT& parameter)
{
    const Rational base =
        exact_rational(parameter.a0());
    const Rational coefficient =
        exact_rational(parameter.a1());
    const Rational radicand =
        exact_rational(parameter.root());
    const std::optional<Rational>
        rational_radical =
            rational_square_root(radicand);
    if (coefficient == 0
        || rational_radical.has_value()) {
        const Rational value =
            base
            + coefficient
                * rational_radical.value_or(
                    Rational(0));
        const std::vector<Integer> primitive =
            primitive_coefficients(
                {
                    -value,
                    Rational(1),
                });
        const Polynomial1 factor =
            polynomial_from_integers(primitive);
        Kernel1 kernel;
        return {
            kernel.construct_algebraic_real_1_object()(
                value),
            factor,
            primitive,
            0,
        };
    }
    const std::vector<Rational> rational_factor =
        {
            base * base
                - coefficient * coefficient
                    * radicand,
            Rational(-2) * base,
            Rational(1),
        };
    const std::vector<Integer> primitive =
        primitive_coefficients(rational_factor);
    const Polynomial1 factor =
        polynomial_from_integers(primitive);
    Kernel1 kernel;
    std::vector<AlgebraicReal1> roots;
    kernel.solve_1_object()(
        factor,
        true,
        std::back_inserter(roots));
    for (std::size_t ordinal = 0;
         ordinal < roots.size();
         ++ordinal) {
        const auto interval =
            kernel.isolate_1_object()(
                roots[ordinal],
                factor);
        if (CGAL::compare(
                parameter,
                GpsPoint::CoordNT(
                    exact_ft(interval.first)))
                != CGAL::SMALLER
            && CGAL::compare(
                   parameter,
                   GpsPoint::CoordNT(
                       exact_ft(interval.second)))
                != CGAL::LARGER) {
            return {
                roots[ordinal],
                factor,
                primitive,
                ordinal,
            };
        }
    }
    throw IncompleteSegmentPartitionError(
        "trimmed rim parameter did not match its exact root");
}

Polynomial2 polynomial_in_x(
    const std::vector<Integer>& coefficients)
{
    const Polynomial2 x =
        CGAL::shift(
            Polynomial2(Integer(1)),
            1,
            0);
    Polynomial2 result(Integer(0));
    Polynomial2 power(Integer(1));
    for (const Integer& coefficient : coefficients) {
        result += coefficient * power;
        power *= x;
    }
    return result;
}

Polynomial2 polynomial_in_y(
    const std::vector<Integer>& coefficients)
{
    const Polynomial2 y =
        CGAL::shift(
            Polynomial2(Integer(1)),
            1,
            1);
    Polynomial2 result(Integer(0));
    Polynomial2 power(Integer(1));
    for (const Integer& coefficient : coefficients) {
        result += coefficient * power;
        power *= y;
    }
    return result;
}

AlgebraicReal2 algebraic_pair(
    const SegmentBranchState2& first,
    const SegmentBranchState2& second)
{
    Kernel2 kernel;
    return kernel.construct_algebraic_real_2_object()(
        first.rim_parameter,
        second.rim_parameter);
}

int chart_sign(const std::string& chart)
{
    if (chart == "rim-half-0-v1") {
        return 1;
    }
    if (chart == "rim-half-1-v1") {
        return -1;
    }
    throw ChartCoverageError(
        "branch does not use a frozen rim chart");
}

struct PairPolynomials {
    Polynomial2 determinant;
    Polynomial2 dot;
    Polynomial2 cap_difference;
};

PairPolynomials pair_polynomials(
    const std::string& first_chart,
    const std::string& second_chart,
    const Rational& cap)
{
    const Polynomial2 x =
        CGAL::shift(
            Polynomial2(Integer(1)),
            1,
            0);
    const Polynomial2 y =
        CGAL::shift(
            Polynomial2(Integer(1)),
            1,
            1);
    const Polynomial2 one(Integer(1));
    const Polynomial2 first_denominator =
        one + x * x;
    const Polynomial2 second_denominator =
        one + y * y;
    const Integer first_sign =
        chart_sign(first_chart);
    const Integer second_sign =
        chart_sign(second_chart);
    const Polynomial2 first_x =
        first_sign * (one - x * x);
    const Polynomial2 first_y =
        Integer(2) * first_sign * x;
    const Polynomial2 second_x =
        second_sign * (one - y * y);
    const Polynomial2 second_y =
        Integer(2) * second_sign * y;
    const Polynomial2 determinant =
        first_x * second_y
        - first_y * second_x;
    const Polynomial2 dot =
        first_x * second_x
        + first_y * second_y;
    const Polynomial2 chord_x =
        first_x * second_denominator
        - second_x * first_denominator;
    const Polynomial2 chord_y =
        first_y * second_denominator
        - second_y * first_denominator;
    const Polynomial2 cap_difference =
        CORE::denominator(cap)
            * (
                chord_x * chord_x
                + chord_y * chord_y)
        - CORE::numerator(cap)
            * first_denominator * first_denominator
            * second_denominator * second_denominator;
    return {
        determinant,
        dot,
        cap_difference,
    };
}

std::string sign_text(CGAL::Sign sign)
{
    if (sign == CGAL::NEGATIVE) {
        return "negative";
    }
    if (sign == CGAL::ZERO) {
        return "zero";
    }
    return "positive";
}

} // namespace

std::vector<SegmentBranchState2> segment_branches_at(
    const std::vector<BoundaryFeatureRecord2>& records,
    const SegmentEventSource2& source,
    const std::string& witness_numerator,
    const std::string& witness_denominator)
{
    const Rational parameter{
        Integer(witness_numerator),
        Integer(witness_denominator)};
    const Rational x0 =
        parse_rational(source.x0().text());
    const Rational y0 =
        parse_rational(source.y0().text());
    const Rational x1 =
        parse_rational(source.x1().text());
    const Rational y1 =
        parse_rational(source.y1().text());
    const Epeck::FT center_x = exact_ft(
        x0 + parameter * (x1 - x0));
    const Epeck::FT center_y = exact_ft(
        y0 + parameter * (y1 - y0));
    const Epeck::FT radius = exact_ft(
        parse_rational(source.tool_radius().text()));
    const GpsPolygon cutter = cutter_polygon(
        EPoint(center_x, center_y),
        radius);

    std::vector<SegmentBranchState2> branches;
    for (const BoundaryFeatureRecord2& record :
         records) {
        std::vector<SegmentBranchState2>
            feature_branches;
        for (const GpsPoint& point :
             trimmed_hits(record, cutter)) {
            const std::string chart =
                chart_for_point(
                    point,
                    center_x,
                    center_y);
            const RimRoot root = identify_rim_root(
                rim_parameter(
                    point,
                    center_x,
                    center_y,
                    radius,
                    chart));
            const std::string root_id =
                algebraic_root_id_v1(
                    root.coefficients,
                    root.ordinal);
            std::string vertex_id;
            if (same_point(
                    point,
                    record.curve.source())) {
                vertex_id = record.source_vertex_id;
            } else if (
                same_point(
                    point,
                    record.curve.target())) {
                vertex_id = record.target_vertex_id;
            }
            feature_branches.push_back(
                {
                    {
                        {},
                        record.feature_id,
                        record.support_id,
                        record.support_kind,
                        record.trim_predicate,
                        std::move(vertex_id),
                        record.material_side,
                        "accepted",
                        chart,
                        0,
                        root_id,
                        coefficient_text(
                            root.coefficients),
                        root.ordinal,
                    },
                    root.value,
                    root.factor,
                    0,
                });
        }
        Kernel1 feature_kernel;
        std::sort(
            feature_branches.begin(),
            feature_branches.end(),
            [&feature_kernel](
                const SegmentBranchState2& first,
                const SegmentBranchState2& second) {
                if (first.branch.rim_chart_id
                    != second.branch.rim_chart_id) {
                    return first.branch.rim_chart_id
                        < second.branch.rim_chart_id;
                }
                return feature_kernel.compare_1_object()(
                           first.rim_parameter,
                           second.rim_parameter)
                    == CGAL::SMALLER;
            });
        std::string current_chart;
        std::size_t sheet_ordinal = 0;
        for (SegmentBranchState2& state :
             feature_branches) {
            if (state.branch.rim_chart_id
                != current_chart) {
                current_chart =
                    state.branch.rim_chart_id;
                sheet_ordinal = 0;
            } else {
                ++sheet_ordinal;
            }
            state.branch.rim_sheet_ordinal =
                sheet_ordinal;
            state.branch.branch_id =
                encode_string_sequence(
                    {
                        "segment-boundary-branch-v1",
                        record.feature_id,
                        current_chart,
                        std::to_string(
                            sheet_ordinal),
                    });
            branches.push_back(
                std::move(state));
        }
    }
    Kernel1 kernel;
    std::sort(
        branches.begin(),
        branches.end(),
        [&kernel](
            const SegmentBranchState2& first,
            const SegmentBranchState2& second) {
            if (first.branch.rim_chart_id
                != second.branch.rim_chart_id) {
                return first.branch.rim_chart_id
                    < second.branch.rim_chart_id;
            }
            const CGAL::Comparison_result comparison =
                kernel.compare_1_object()(
                    first.rim_parameter,
                    second.rim_parameter);
            if (comparison != CGAL::EQUAL) {
                return comparison == CGAL::SMALLER;
            }
            return first.branch.branch_id
                < second.branch.branch_id;
        });
    for (std::size_t index = 0;
         index < branches.size();
         ++index) {
        branches[index].rim_order = index;
    }
    return branches;
}

std::vector<BranchPairDisposition2>
segment_branch_pair_dispositions(
    const std::vector<SegmentBranchState2>& branches,
    const SegmentEventSource2& source)
{
    Kernel2 kernel;
    const Rational cap_chord_ratio =
        parse_rational(
            source.cap_chord_ratio().text());
    std::vector<BranchPairDisposition2> result;
    for (std::size_t first_index = 0;
         first_index < branches.size();
         ++first_index) {
        for (std::size_t second_index = 0;
             second_index < branches.size();
             ++second_index) {
            if (first_index == second_index) {
                continue;
            }
            const SegmentBranchState2& first =
                branches[first_index];
            const SegmentBranchState2& second =
                branches[second_index];
            const AlgebraicReal2 pair =
                algebraic_pair(first, second);
            const PairPolynomials polynomials =
                pair_polynomials(
                    first.branch.rim_chart_id,
                    second.branch.rim_chart_id,
                    cap_chord_ratio);
            const CGAL::Sign determinant =
                exact_polynomial_sign_at_2(
                    polynomials.determinant,
                    pair);
            const CGAL::Sign dot =
                exact_polynomial_sign_at_2(
                    polynomials.dot,
                    pair);
            const CGAL::Sign cap =
                exact_polynomial_sign_at_2(
                    polynomials.cap_difference,
                    pair);
            std::string orientation;
            std::string cap_disposition;
            if (determinant == CGAL::POSITIVE) {
                orientation = "ccw-minor";
                cap_disposition =
                    cap == CGAL::NEGATIVE
                    ? "below-cap"
                    : (
                          cap == CGAL::ZERO
                              ? "equal-cap"
                              : "above-cap");
            } else if (
                determinant == CGAL::NEGATIVE) {
                orientation = "ccw-major";
                cap_disposition = "above-cap";
            } else if (dot == CGAL::POSITIVE) {
                orientation = "ccw-zero";
                cap_disposition = "below-cap";
            } else {
                orientation = "ccw-pi";
                cap_disposition =
                    cap_chord_ratio == Rational(4)
                    ? "equal-cap"
                    : "above-cap";
            }
            result.push_back(
                {
                    encode_string_sequence(
                        {
                            "segment-ordered-branch-pair-v1",
                            first.branch.branch_id,
                            second.branch.branch_id,
                        }),
                    first.branch.branch_id,
                    second.branch.branch_id,
                    orientation,
                    cap_disposition,
                });
        }
    }
    return result;
}

SegmentCellStratum2 construct_segment_cell_stratum(
    const Stock2& stock,
    const SegmentEventSource2& source,
    const std::string& witness_numerator,
    const std::string& witness_denominator)
{
    const std::vector<BoundaryFeatureRecord2> records =
        extract_boundary_records(stock);
    if (records.empty()) {
        throw IncompleteSegmentPartitionError(
            "segment cell requires a nonempty stock boundary");
    }
    const std::vector<SegmentBranchState2> states =
        segment_branches_at(
            records,
            source,
            witness_numerator,
            witness_denominator);
    std::vector<SegmentBoundaryBranch2> branches;
    std::vector<std::string> active_branch_ids;
    branches.reserve(states.size());
    active_branch_ids.reserve(states.size());
    for (const SegmentBranchState2& state : states) {
        branches.push_back(state.branch);
        active_branch_ids.push_back(
            state.branch.branch_id);
    }
    return {
        std::move(branches),
        {
            "cell",
            {},
            {},
            {},
            "segment-linear-v1",
            witness_numerator,
            witness_denominator,
            {},
            0,
            std::move(active_branch_ids),
            {},
            {},
            segment_branch_pair_dispositions(
                states,
                source),
            {},
            false,
            false,
            false,
            false,
        },
    };
}

std::vector<std::string>
segment_pair_literal_signs()
{
    const std::vector<Integer> first_factor{
        Integer(-1),
        Integer(3),
    };
    const std::vector<Integer> second_factor{
        Integer(1),
        Integer(3),
    };
    Kernel2 kernel;
    std::vector<
        std::pair<
            AlgebraicReal2,
            Kernel2::Multiplicity_type>>
        solutions;
    kernel.solve_2_object()(
        polynomial_in_x(first_factor),
        polynomial_in_y(second_factor),
        std::back_inserter(solutions));
    if (solutions.size() != 1) {
        throw IncompleteSegmentPartitionError(
            "literal common-domain diagnostic did not isolate one pair");
    }
    const PairPolynomials polynomials =
        pair_polynomials(
            "rim-half-0-v1",
            "rim-half-1-v1",
            Rational(4));
    const Polynomial2 x =
        CGAL::shift(
            Polynomial2(Integer(1)),
            1,
            0);
    const Polynomial2 y =
        CGAL::shift(
            Polynomial2(Integer(1)),
            1,
            1);
    const Polynomial2 negative_factor =
        x * y - Polynomial2(Integer(1));
    const Polynomial2 zero_factor =
        Integer(3) * x
        - Polynomial2(Integer(1));
    const AlgebraicReal2& point =
        solutions.front().first;
    return {
        sign_text(
            exact_polynomial_sign_at_2(
                polynomials.determinant,
                point)),
        sign_text(
            exact_polynomial_sign_at_2(
                -polynomials.determinant,
                point)),
        sign_text(
            exact_polynomial_sign_at_2(
                polynomials.dot,
                point)),
        sign_text(
            exact_polynomial_sign_at_2(
                -polynomials.dot,
                point)),
        sign_text(
            exact_polynomial_sign_at_2(
                negative_factor
                    * negative_factor,
                point)),
        sign_text(
            exact_polynomial_sign_at_2(
                negative_factor
                    * negative_factor
                    * negative_factor,
                point)),
        sign_text(
            exact_polynomial_sign_at_2(
                polynomials.cap_difference,
                point)),
        sign_text(
            exact_polynomial_sign_at_2(
                zero_factor
                    * zero_factor,
                point)),
    };
}

std::vector<std::string>
segment_rational_square_root_cases()
{
    const std::optional<Rational> rational =
        rational_square_root(
            Rational(4, 9));
    const std::optional<Rational> irrational =
        rational_square_root(
            Rational(2));
    const std::optional<Rational> zero =
        rational_square_root(
            Rational(0));
    const std::optional<Rational> negative =
        rational_square_root(
            Rational(-1));
    std::string rational_text =
        "missing-rational";
    if (rational.has_value()) {
        Integer numerator;
        Integer denominator;
        typename CGAL::Fraction_traits<
            Rational>::Decompose()(
            rational.value(),
            numerator,
            denominator);
        rational_text =
            numerator.convert_to<std::string>()
            + "/"
            + denominator.convert_to<
                std::string>();
    }
    return {
        rational_text,
        irrational.has_value()
            ? "accepted-non-square"
            : "non-square",
        zero.has_value()
            && zero.value() == 0
            ? "zero"
            : "missing-zero",
        negative.has_value()
            ? "accepted-negative"
            : "negative",
    };
}
