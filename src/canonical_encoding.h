#pragma once

#include "exact_algebraic_1.h"

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <CGAL/CORE/BigRat.h>

class InvalidCanonicalEncodingError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

std::string canonical_encode_bytes(const std::string &value);
std::string canonical_encode_integer(const ExactAlgebraicInteger1 &value);
std::string canonical_encode_rational(const CORE::BigRat &value);
CORE::BigRat canonical_decode_rational(const std::string &value);
std::string canonical_encode_boolean(bool value);
std::string canonical_encode_sequence(const std::vector<std::string> &values);
std::string canonical_encode_component_map(
    const std::vector<std::pair<std::string, std::string>> &fields);
std::string canonical_encode_tagged_union(const std::string &tag,
                                          const std::string &payload);
