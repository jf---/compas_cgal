#pragma once

#include "audit_digest_2.h"
#include "compas_matrix.h"

#include <stdexcept>
#include <string>
#include <vector>

class AuditNativeStockIdentityError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

class AuditNativeStockShapeError : public AuditNativeStockIdentityError {
public:
    using AuditNativeStockIdentityError::AuditNativeStockIdentityError;
};

class AuditNativeStockNonFiniteInputError : public AuditNativeStockIdentityError {
public:
    using AuditNativeStockIdentityError::AuditNativeStockIdentityError;
};

class AuditNativeStockRingError : public AuditNativeStockIdentityError {
public:
    using AuditNativeStockIdentityError::AuditNativeStockIdentityError;
};

class AuditNativeStockDuplicateHoleError : public AuditNativeStockIdentityError {
public:
    using AuditNativeStockIdentityError::AuditNativeStockIdentityError;
};

class AuditNativeStockIdentity2 {
public:
    static AuditNativeStockIdentity2 build(
        Eigen::Ref<const compas::RowMatrixXd> boundary,
        const std::vector<compas::RowMatrixXd>& holes);

    const std::string& canonical_bytes() const noexcept;
    const AuditNativeStockDigest2& digest() const noexcept;

private:
    AuditNativeStockIdentity2(
        std::string canonical_bytes,
        AuditNativeStockDigest2 digest);

    std::string canonical_bytes_;
    AuditNativeStockDigest2 digest_;
};
