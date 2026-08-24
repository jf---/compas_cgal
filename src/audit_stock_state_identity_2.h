#pragma once

#include "audit_digest_2.h"
#include "stock_2.h"

#include <string>

class AuditStockStateIdentity2 {
public:
    static AuditStockStateIdentity2 build(const Stock2& stock);
    const std::string& canonical_bytes() const noexcept;
    const AuditStockStateDigest2& digest() const noexcept;

private:
    AuditStockStateIdentity2(
        std::string canonical_bytes,
        AuditStockStateDigest2 digest);
    std::string canonical_bytes_;
    AuditStockStateDigest2 digest_;
};
