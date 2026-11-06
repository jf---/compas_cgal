#pragma once

// The interval filter of the rim engine: certifies one circle's stations WITHIN the
// cap on an emission range of swept annuli, without one exact number.
//
// It is one-sided, on purpose. For each station it builds M+, a superset of the exact
// material on the rim: whatever it removes from the rim is certainly not material,
// and whatever it cannot prove stays. A station is certified when every maximal arc
// of M+ is certainly within the cap. Clip monotonicity (swept_prefix_engagement_2.h)
// and the monotone cap test make that a proof about the exact engine's verdict.
// Nothing is ever certified exceeded: a station the filter cannot certify goes to the
// next view, then to the exact path, whose verdict stands. A refused decision costs
// one exact station.
//
// The covers (what is removed from the rim), each an inner bound of the true set:
//   own prefix   the rear half-rim and the start disk. On a circle's station k > 0 a
//                rim point is covered by the prefix iff it is on the rear half or in
//                the start disk, at every extent: a rear rim point is in the sector or
//                nearer the start than the station (the squared distance to a guide
//                point grows with the angular difference), and the same holds for a
//                forward rim point the sector reaches past the start ray.
//   annuli       on a rim at distance d from c_j, the disk of radius R covers the arc
//                of half-angle acos((d^2 + r^2 - R^2) / (2 d r)) about the direction
//                of c_j: outer arc minus inner arc.
//   walls        the arcs between a certain leaving crossing and the next certain
//                entering one. Any crossing that is not certain (an enclosure around a
//                tangency or a vertex) either leaves the arcs it cannot touch alone, or
//                drops every wall cover of the station.
//
// Positions on the rim are compared by a key: r - x on the upper half, 3 r + x on the
// lower half, so half a turn is exactly 2 r of key. Each cover keeps the far side of
// its enclosures (start's upper bound, end's lower bound), so the comparison of two
// keys, a comparison of doubles, can only leave more material.
//
// The cap, for an arc of M+ from end point p to start point q of its two covers:
//   orientation certainly positive and the chord squared certainly below the cap (the
//   two legs of run_exceeds_cap, exact_cap_certificate_2.h); or, when the two ends are
//   too close for an orientation, the key length l: within one half the chord squared
//   is at most 2 r l (with x = r cos t, l >= 2 r sin^2(dt/2)), over two halves at most
//   4 r l, and l < 2 r is less than half a turn.

#include "annulus_store_2.h"

#include <CGAL/Interval_nt.h>

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

// The certifying view of one station; UNCERTIFIED goes to the exact path.
enum class AnnulusFilterTier2 : std::uint8_t {
    LAST_ONE = 0,
    LAST_SIXTEEN = 1,
    RANGE = 2,
    UNCERTIFIED = 3,
};

// One station's result: its tier and an upper bound of the chord squared of its largest
// M+ arc (0 when M+ is empty), which an exact run on the same view never exceeds.
struct AnnulusFilterStation2 {
    AnnulusFilterTier2 tier;
    double chord_squared_bound;
};

// The exact rotation of a motion's phase at one parameter, enclosed. The entry station
// (parameter 0) has none.
struct AnnulusFilterRotation2 {
    CGAL::Interval_nt<false> cosine;
    CGAL::Interval_nt<false> sine;
    bool entry;
};

// A part edge with the part on its left.
struct AnnulusFilterEdge2 {
    double ax;
    double ay;
    double bx;
    double by;
};

// One circle: guide centre and radius, entry. Injected doubles.
struct AnnulusFilterCircle2 {
    double cx;
    double cy;
    double guide_radius;
    double entry_x;
    double entry_y;
};

class AnnulusFilter2 {
public:
    AnnulusFilter2(double tool_radius, const std::vector<AnnulusFilterEdge2>& edges);

    // `unit_t` and `unit_reflected` are rational_unit_steering's construction of the start
    // phase's unit. `rotations` holds one entry per station, and `out` is resized to it.
    void certify_circle(
        const AnnulusStore2& store, AnnulusViewRange2 range,
        const AnnulusFilterCircle2& circle, double unit_t, bool unit_reflected,
        const std::vector<AnnulusFilterRotation2>& rotations, double cap_chord_ratio,
        std::vector<AnnulusFilterStation2>& out) const;

private:
    using Interval = CGAL::Interval_nt<false>;
    struct Point {
        Interval x;
        Interval y;
    };
    struct Edge {
        double ax;
        double ay;
        double min_x;
        double min_y;
        double max_x;
        double max_y;
        Interval dx;
        Interval dy;
        Interval squared_length;
        Interval inverse_squared_length;
    };
    struct Piece {
        double low;
        double high;
        int low_point;
        int high_point;
    };
    struct WallEvent {
        double low;
        double high;
        bool entering;
        int point;
    };
    struct RimBox {
        double min_x;
        double min_y;
        double max_x;
        double max_y;
    };

    RimBox rim_box(const Point& station) const;
    bool key_of(const Point& point, double& low, double& high) const;
    int add_point(const Point& point) const;
    void add_cover(const Point& start, const Point& end, bool at_most_half_turn) const;
    void add_own_prefix(const Point& phase, const Point& offset, const Interval& scale) const;
    void add_walls(const Point& station) const;
    bool add_annulus(const Point& station, const StoredAnnulus2& annulus) const;
    bool add_tier(
        int tier, const AnnulusStore2& store, AnnulusViewRange2 range, const Point& station,
        const AnnulusSearchBox2& read_box, bool& gathered) const;
    bool key_length_within(double length) const;
    bool arc_within(int from, int to, double length) const;
    bool certify() const;

    Interval tool_radius_;
    Interval rim_squared_;
    double half_turn_key_;
    double full_turn_key_;
    double key_width_limit_;
    std::vector<Edge> edges_;

    // Scratch, reused between decisions; a filter belongs to one view and one thread.
    mutable double cap_squared_ = 0.0;
    mutable double chord_bound_ = 0.0;
    mutable std::vector<const Edge*> wall_candidates_;
    mutable std::vector<std::uint32_t> range_candidates_;
    mutable std::vector<Point> points_;
    mutable std::vector<Piece> pieces_;
    mutable std::vector<Piece> sorted_;
    mutable std::vector<WallEvent> events_;
    mutable std::vector<std::pair<double, double>> tangencies_;
};
