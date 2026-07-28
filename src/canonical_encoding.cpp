#include "canonical_encoding.h"

#include <algorithm>
#include <cstdint>
#include <iterator>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_int/import_export.hpp>

namespace {

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
