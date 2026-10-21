#include "held_declared_clearance_2/wire_2.h"

#include "held_declared_clearance_2/errors_2.h"

#include <CGAL/Fraction_traits.h>

#include <boost/multiprecision/cpp_int.hpp>

#include <cmath>

namespace compas_cgal::held_declared_clearance {

namespace {

using RationalTraits = CGAL::Fraction_traits<FT::ET>;
using Integer = RationalTraits::Numerator_type;

constexpr std::size_t LINE_RATIONALS = 7;
constexpr std::size_t ARC_RATIONALS = 7;
constexpr std::size_t SEGMENT_ELEMENT_RATIONALS = 5;
constexpr std::size_t ARC_ELEMENT_RATIONALS = 3;

Integer integer_from_text(const std::string& text)
{
    if (text.empty()) throw InvalidClearanceWireError("empty integer text");
    const std::size_t first = text[0] == '-' ? 1 : 0;
    if (first == text.size()) throw InvalidClearanceWireError("integer text has no digits");
    for (std::size_t i = first; i < text.size(); ++i) {
        if (text[i] < '0' || text[i] > '9') throw InvalidClearanceWireError("integer text is not decimal: " + text);
    }
    return Integer(text.c_str());
}

}  // namespace

FT rational_from_wire(const RationalWire& wire)
{
    const Integer numerator = integer_from_text(wire.first);
    const Integer denominator = integer_from_text(wire.second);
    if (denominator <= 0) throw InvalidClearanceWireError("a rational's denominator is not positive");
    return FT(RationalTraits::Compose()(numerator, denominator));
}

Root root_from_wire(const RootWire& wire)
{
    return make_root(rational_from_wire(std::get<0>(wire)), rational_from_wire(std::get<1>(wire)),
                     rational_from_wire(std::get<2>(wire)));
}

OneRootPoint point_from_wire(const PointWire& wire)
{
    return make_point(root_from_wire(wire.first), root_from_wire(wire.second));
}

FT binary64_to_exact(double value)
{
    if (!std::isfinite(value)) throw InvalidClearanceWireError("binary64 inputs must be finite");
    return FT(value);  // exact: every finite binary64 is a dyadic rational
}

Primitive primitive_from_wire(std::size_t index, const PrimitiveWire& wire)
{
    const auto& [kind, sense, values] = wire;
    if (kind == LINE_KIND) {
        if (values.size() != LINE_RATIONALS) throw InvalidClearanceWireError("a line needs 7 rationals");
        return Primitive::line(index, RationalPoint(rational_from_wire(values[0]), rational_from_wire(values[1])),
                               RationalPoint(rational_from_wire(values[2]), rational_from_wire(values[3])),
                               Vector(rational_from_wire(values[4]), rational_from_wire(values[5])),
                               rational_from_wire(values[6]));
    }
    if (kind == ARC_KIND) {
        if (values.size() != ARC_RATIONALS) throw InvalidClearanceWireError("an arc needs 7 rationals");
        return Primitive::arc(index, RationalPoint(rational_from_wire(values[0]), rational_from_wire(values[1])),
                              RationalPoint(rational_from_wire(values[2]), rational_from_wire(values[3])),
                              RationalPoint(rational_from_wire(values[4]), rational_from_wire(values[5])),
                              rational_from_wire(values[6]), sense);
    }
    throw InvalidClearanceWireError("unknown primitive kind");
}

DecodedElement element_from_wire(const ElementWire& wire, std::size_t primitive_count)
{
    const auto& [kind, sense, owners, start_wire, end_wire, values] = wire;
    DecodedElement decoded;
    for (const std::size_t owner : owners) {
        if (owner >= primitive_count) throw InvalidClearanceWireError("an owner index is out of range");
    }
    decoded.owners = owners;
    const OneRootPoint start = point_from_wire(start_wire);
    const OneRootPoint end = point_from_wire(end_wire);
    if (kind == LINE_KIND) {
        if (values.size() != SEGMENT_ELEMENT_RATIONALS) throw InvalidClearanceWireError("a segment element needs 5 rationals");
        decoded.segment = SegmentElement::build(start, end, Vector(rational_from_wire(values[0]), rational_from_wire(values[1])),
                                                Vector(rational_from_wire(values[2]), rational_from_wire(values[3])),
                                                rational_from_wire(values[4]));
        return decoded;
    }
    if (kind == ARC_KIND) {
        if (values.size() != ARC_ELEMENT_RATIONALS) throw InvalidClearanceWireError("an arc element needs 3 rationals");
        decoded.is_arc = true;
        decoded.arc = ArcElement::build(RationalPoint(rational_from_wire(values[0]), rational_from_wire(values[1])),
                                        rational_from_wire(values[2]), sense, start, end);
        return decoded;
    }
    throw InvalidClearanceWireError("unknown element kind");
}

}  // namespace compas_cgal::held_declared_clearance
