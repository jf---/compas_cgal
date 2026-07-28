#include "canonical_encoding.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

#include <boost/multiprecision/cpp_int/import_export.hpp>

namespace {

inline constexpr std::size_t CANONICAL_PREFIX_SIZE = 6;
inline constexpr std::size_t CANONICAL_SIZE_SIZE = 8;
inline constexpr std::size_t CANONICAL_HEADER_SIZE =
    CANONICAL_PREFIX_SIZE + 1 + CANONICAL_SIZE_SIZE;

struct ParsedCanonicalNode {
  char kind;
  std::string_view payload;
  std::size_t end;
};

std::string unsigned_64_record(const std::size_t value) {
  static_assert(sizeof(std::size_t) <= sizeof(std::uint64_t));
  std::string result;
  result.reserve(sizeof(std::uint64_t));
  for (int shift = 56; shift >= 0; shift -= 8) {
    result.push_back(static_cast<char>(
        (static_cast<std::uint64_t>(value) >> shift) & 0xffU));
  }
  return result;
}

std::string canonical_node(const char kind, const std::string &payload) {
  std::string result("CCAN\0\1", 6);
  result.push_back(kind);
  result += unsigned_64_record(payload.size());
  result += payload;
  return result;
}

} // namespace

std::string canonical_encode_bytes(const std::string &value) {
  return canonical_node('B', value);
}

std::string canonical_encode_integer(const ExactAlgebraicInteger1 &value) {
  const ExactAlgebraicInteger1 magnitude = value < 0 ? -value : value;
  std::vector<unsigned char> magnitude_bytes;
  if (magnitude != 0) {
    export_bits(magnitude, std::back_inserter(magnitude_bytes), 8, true);
  }

  std::string payload(1,
                      value < 0 ? static_cast<char>(1) : static_cast<char>(0));
  for (const unsigned char byte : magnitude_bytes) {
    payload.push_back(static_cast<char>(byte));
  }
  return canonical_node('I', payload);
}

std::string canonical_encode_binary64(double value) {
  if (!std::isfinite(value)) {
    throw InvalidCanonicalEncodingError(
        "canonical binary64 value must be finite");
  }
  if (value == 0.0) {
    value = 0.0;
  }
  const std::uint64_t bits = std::bit_cast<std::uint64_t>(value);
  std::string payload;
  payload.reserve(sizeof(bits));
  for (int shift = 56; shift >= 0; shift -= 8) {
    payload.push_back(static_cast<char>((bits >> shift) & 0xffU));
  }
  return canonical_node('D', payload);
}

namespace {

std::uint64_t decode_unsigned_64(const std::string_view value,
                                 const std::size_t offset) {
  std::uint64_t result = 0;
  for (std::size_t index = 0; index < CANONICAL_SIZE_SIZE; ++index) {
    result = (result << 8U) | static_cast<unsigned char>(value[offset + index]);
  }
  return result;
}

ParsedCanonicalNode parse_canonical_node(const std::string_view value,
                                         const std::size_t offset) {
  if (offset > value.size() || value.size() - offset < CANONICAL_HEADER_SIZE ||
      value.substr(offset, CANONICAL_PREFIX_SIZE) !=
          std::string_view("CCAN\0\1", CANONICAL_PREFIX_SIZE)) {
    throw InvalidCanonicalEncodingError(
        "canonical bytes must contain one complete CCAN record");
  }
  const std::size_t kind_offset = offset + CANONICAL_PREFIX_SIZE;
  const std::size_t size_offset = kind_offset + 1;
  const std::uint64_t payload_size = decode_unsigned_64(value, size_offset);
  if (payload_size > std::numeric_limits<std::size_t>::max()) {
    throw InvalidCanonicalEncodingError(
        "canonical node payload exceeds native address space");
  }
  const std::size_t payload_offset = size_offset + CANONICAL_SIZE_SIZE;
  const std::size_t native_payload_size =
      static_cast<std::size_t>(payload_size);
  if (native_payload_size > value.size() - payload_offset) {
    throw InvalidCanonicalEncodingError("canonical node payload is truncated");
  }
  return {
      value[kind_offset],
      value.substr(payload_offset, native_payload_size),
      payload_offset + native_payload_size,
  };
}

ExactAlgebraicInteger1
decode_canonical_integer(const ParsedCanonicalNode &node) {
  if (node.kind != 'I' || node.payload.empty()) {
    throw InvalidCanonicalEncodingError(
        "canonical rational requires two integer child nodes");
  }
  const unsigned char sign = static_cast<unsigned char>(node.payload.front());
  if (sign > 1U) {
    throw InvalidCanonicalEncodingError(
        "canonical integer has an invalid sign byte");
  }
  const std::string_view magnitude_bytes = node.payload.substr(1);
  if (!magnitude_bytes.empty() && magnitude_bytes.front() == '\0') {
    throw InvalidCanonicalEncodingError(
        "canonical integer magnitude is not minimal");
  }
  ExactAlgebraicInteger1 magnitude = 0;
  for (const char byte : magnitude_bytes) {
    magnitude <<= 8;
    magnitude += static_cast<unsigned char>(byte);
  }
  if (sign == 1U && magnitude == 0) {
    throw InvalidCanonicalEncodingError(
        "canonical integer cannot encode negative zero");
  }
  return sign == 1U ? -magnitude : magnitude;
}

} // namespace

std::string canonical_encode_rational(const CORE::BigRat &value) {
  return canonical_node('R',
                        canonical_encode_integer(CORE::numerator(value)) +
                            canonical_encode_integer(CORE::denominator(value)));
}

CORE::BigRat canonical_decode_rational(const std::string &value) {
  const ParsedCanonicalNode rational = parse_canonical_node(value, 0);
  if (rational.kind != 'R' || rational.end != value.size()) {
    throw InvalidCanonicalEncodingError(
        "canonical rational must be one complete rational record");
  }
  const ParsedCanonicalNode numerator_node =
      parse_canonical_node(rational.payload, 0);
  const ParsedCanonicalNode denominator_node =
      parse_canonical_node(rational.payload, numerator_node.end);
  if (denominator_node.end != rational.payload.size()) {
    throw InvalidCanonicalEncodingError(
        "canonical rational contains trailing child bytes");
  }
  const ExactAlgebraicInteger1 numerator =
      decode_canonical_integer(numerator_node);
  const ExactAlgebraicInteger1 denominator =
      decode_canonical_integer(denominator_node);
  const ExactAlgebraicInteger1 magnitude =
      numerator < 0 ? -numerator : numerator;
  if (denominator <= 0 || CORE::gcd(magnitude, denominator) != 1) {
    throw InvalidCanonicalEncodingError(
        "canonical rational must be reduced with positive denominator");
  }
  return CORE::BigRat(numerator, denominator);
}

std::string canonical_encode_boolean(const bool value) {
  return canonical_node(
      '?', std::string(1, value ? static_cast<char>(1) : static_cast<char>(0)));
}

std::string canonical_encode_sequence(const std::vector<std::string> &values) {
  std::string payload = unsigned_64_record(values.size());
  for (const std::string &value : values) {
    payload += canonical_encode_bytes(value);
  }
  return canonical_node('S', payload);
}

std::string canonical_encode_component_map(
    const std::vector<std::pair<std::string, std::string>> &fields) {
  std::vector<std::pair<std::string, std::string>> encoded;
  encoded.reserve(fields.size());
  for (const auto &[key, value] : fields) {
    if (key.empty()) {
      throw InvalidCanonicalEncodingError(
          "canonical component name must be nonempty");
    }
    encoded.emplace_back(canonical_encode_bytes(key),
                         canonical_encode_bytes(value));
  }
  std::sort(encoded.begin(), encoded.end());
  if (std::adjacent_find(encoded.begin(), encoded.end(),
                         [](const auto &left, const auto &right) {
                           return left.first == right.first;
                         }) != encoded.end()) {
    throw InvalidCanonicalEncodingError(
        "canonical component names must be unique");
  }

  std::string payload = unsigned_64_record(encoded.size());
  for (const auto &[key, value] : encoded) {
    payload += key;
    payload += value;
  }
  return canonical_node('M', payload);
}

std::string canonical_encode_tagged_union(const std::string &tag,
                                          const std::string &payload) {
  if (tag.empty()) {
    throw InvalidCanonicalEncodingError("canonical union tag must be nonempty");
  }
  return canonical_node('T', canonical_encode_bytes(tag) +
                                 canonical_encode_bytes(payload));
}
