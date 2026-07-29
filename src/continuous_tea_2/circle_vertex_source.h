#pragma once

#include "../stock_2.h"

#include <string>

#include <CGAL/CORE/BigRat.h>

class FullCircleCoordinate2 {
public:
    static FullCircleCoordinate2 from_exact(
        const GpsPoint::CoordNT& coordinate);

    static FullCircleCoordinate2 decode(
        const std::string& source);

    std::string canonical_source() const;

    const CORE::BigRat& base() const
    {
        return base_;
    }

    const CORE::BigRat& radical_coefficient() const
    {
        return radical_coefficient_;
    }

    const CORE::BigRat& radicand() const
    {
        return radicand_;
    }

private:
    FullCircleCoordinate2(
        CORE::BigRat base,
        CORE::BigRat radical_coefficient,
        CORE::BigRat radicand);

    CORE::BigRat base_;
    CORE::BigRat radical_coefficient_;
    CORE::BigRat radicand_;
};
