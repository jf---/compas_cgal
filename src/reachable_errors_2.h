#pragma once

#include <stdexcept>

class InvalidReachableDomainInputError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class ReachableArrangementTopologyError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class PocketNotMachinableError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class ReachableMaterialContainmentError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};
