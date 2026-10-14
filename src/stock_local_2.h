#pragma once

#include "stock_2.h"

#include <stdexcept>
#include <vector>

// The flood that marks the removed faces escaped the removal region: it reached
// the unbounded face, which is never inside a bounded region. Either a seed was
// not strictly interior or a boundary circle was not fully inserted. Raised
// rather than leaving a stock whose face containment is silently wrong.
class LocalDepletionEscapedError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// A bounded removal region described the way a LOCAL subtraction needs it,
// rather than as the point set a global difference consumes.
//
// The two descriptions are not interchangeable: a global difference only needs
// the region's boundary, while the local update additionally needs to be able to
// (a) RECOGNISE an arrangement edge that lies on that boundary, so the flood
// that re-marks containment never crosses out of the region, and (b) START that
// flood, once per connected component of the region's interior.
struct LocalRemovalRegion {
    // The x-monotone curves of the region boundary, each traversed with the
    // material the removal KEEPS on its left, in the exact form inserted into
    // the stock's own arrangement. They are also the boundary DESCRIPTION: an
    // arrangement edge lies on the boundary iff it is a SUB-CURVE of one of
    // them -- same supporting curve (exact rational coefficients, which survive
    // every split: CGAL 6.0.1 Circle_segment_2.h m_first/m_second/m_third), the
    // same half of a supporting circle, and both endpoints on the closed curve.
    // Inserting a curve makes its endpoints arrangement vertices, so an edge on
    // its support either lies inside it or misses its relative interior: the
    // test is exact and complete, for whole circles and for capsule pieces alike.
    std::vector<GpsXCurve> boundary_curves;

    // At least one exact point STRICTLY inside every connected component of the
    // region's interior. Strictness matters: a seed on the boundary would sit on
    // the separating curve and could seed the wrong side.
    std::vector<GpsPoint> seeds;
};

// The exact closed disk of the given radius about `center`.
LocalRemovalRegion local_disk_region(const EPoint& center, const Epeck::FT& radius);

// The exact closed annulus between `inner_radius` and `outer_radius` about
// `center`. `inner_radius == 0` degenerates to the disk, decided by an exact
// sign test, not a tolerance.
LocalRemovalRegion local_annulus_region(const EPoint& center,
                                        const Epeck::FT& inner_radius,
                                        const Epeck::FT& outer_radius);

// Remove `region` from `set` by LOCAL arrangement surgery, leaving the exact
// same point set -- and the same canonical Gps representation -- that
// `Gps::difference(region)` would produce, without overlaying the whole stock.
//
// See docs/local_depletion.md for the soundness argument; the invariants it
// relies on are restated at each step in stock_local_2.cpp.
void subtract_region_local(Gps& set, const LocalRemovalRegion& region);

void local_depletion_phase_reset();
void local_depletion_phase_report();
