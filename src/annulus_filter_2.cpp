#include "annulus_filter_2.h"

#include <CGAL/FPU.h>

#include <algorithm>
#include <cmath>

namespace {

using Interval = CGAL::Interval_nt<false>;

// An enclosure wider than this many tool radii of key is useless for ordering, and its
// cover is dropped. Dropping is always sound: it only leaves more material.
constexpr double KEY_WIDTH_LIMIT_RADII = 0.5;

inline Interval square(const Interval& value)
{
    return CGAL::square(value);
}

} // namespace

AnnulusFilter2::AnnulusFilter2(double tool_radius, const std::vector<AnnulusFilterEdge2>& edges)
{
    const CGAL::Protect_FPU_rounding<true> rounding;
    tool_radius_ = Interval(tool_radius);
    rim_squared_ = square(tool_radius_);
    half_turn_key_ = (Interval(2.0) * tool_radius_).sup();
    full_turn_key_ = (Interval(4.0) * tool_radius_).sup();
    key_width_limit_ = KEY_WIDTH_LIMIT_RADII * tool_radius;
    edges_.reserve(edges.size());
    for (const AnnulusFilterEdge2& edge : edges) {
        const Interval dx = Interval(edge.bx) - Interval(edge.ax);
        const Interval dy = Interval(edge.by) - Interval(edge.ay);
        const Interval squared_length = square(dx) + square(dy);
        edges_.push_back({
            edge.ax, edge.ay,
            std::min(edge.ax, edge.bx), std::min(edge.ay, edge.by),
            std::max(edge.ax, edge.bx), std::max(edge.ay, edge.by),
            dx, dy, squared_length, Interval(1.0) / squared_length,
        });
    }
}

AnnulusFilter2::RimBox AnnulusFilter2::rim_box(const Point& station) const
{
    const Interval min_x = station.x - tool_radius_;
    const Interval min_y = station.y - tool_radius_;
    const Interval max_x = station.x + tool_radius_;
    const Interval max_y = station.y + tool_radius_;
    return {min_x.inf(), min_y.inf(), max_x.sup(), max_y.sup()};
}

// The enclosure of a rim point's key. Over a half line the enclosure of y straddles 0
// and the point may be on either half; at the right extreme the key is taken below 0 on
// the lower half's side, so the enclosure stays connected over the seam. false: the
// enclosure is too wide to order anything.
bool AnnulusFilter2::key_of(const Point& point, double& low, double& high) const
{
    Interval key;
    if (point.y.inf() > 0.0) {
        key = tool_radius_ - point.x;
    } else if (point.y.sup() < 0.0) {
        key = Interval(3.0) * tool_radius_ + point.x;
    } else if (point.x.inf() > 0.0) {
        const Interval upper = tool_radius_ - point.x;
        const Interval lower = point.x - tool_radius_;
        key = Interval(std::min(upper.inf(), lower.inf()), std::max(upper.sup(), lower.sup()));
    } else if (point.x.sup() < 0.0) {
        const Interval upper = tool_radius_ - point.x;
        const Interval lower = Interval(3.0) * tool_radius_ + point.x;
        key = Interval(std::min(upper.inf(), lower.inf()), std::max(upper.sup(), lower.sup()));
    } else {
        return false;
    }
    low = key.inf();
    high = key.sup();
    return high - low < key_width_limit_;
}

int AnnulusFilter2::add_point(const Point& point) const
{
    points_.push_back(point);
    return static_cast<int>(points_.size()) - 1;
}

// The cover of the CCW rim arc from `start` to `end`, shrunk by their enclosures. A true
// span of at most half a turn whose shrunk ends cross reads as nearly a full turn of key:
// it is near empty, and dropped.
void AnnulusFilter2::add_cover(const Point& start, const Point& end, bool at_most_half_turn) const
{
    double start_low;
    double start_high;
    double end_low;
    double end_high;
    if (!key_of(start, start_low, start_high) || !key_of(end, end_low, end_high)) return;
    double from = start_high;
    double to = end_low;
    if (at_most_half_turn) {
        double forward = to - from;
        if (forward < 0.0) forward += full_turn_key_;
        if (forward > half_turn_key_) return;
    }
    const int start_point = add_point(start);
    const int end_point = add_point(end);
    if (to < 0.0) to = (Interval(to) + Interval(full_turn_key_)).inf();
    if (from > full_turn_key_) from = 0.0;
    if (from <= to) {
        pieces_.push_back({from, to, start_point, end_point});
    } else {
        pieces_.push_back({from, full_turn_key_, start_point, -1});
        if (to > 0.0) pieces_.push_back({0.0, to, -1, end_point});
    }
}

// The rear half-rim, CCW from the inward radial to the outward one (a CCW motion), and the
// start disk. `offset` is station - guide, `phase` is start - guide, so neither pays the
// centre's magnitude. Both circles have radius r, so the start disk's a is exactly 1/2.
void AnnulusFilter2::add_own_prefix(const Point& phase, const Point& offset, const Interval& scale) const
{
    const Point outward{scale * offset.x, scale * offset.y};
    add_cover({-outward.x, -outward.y}, outward, false);
    const Point to_start{phase.x - offset.x, phase.y - offset.y};
    const Interval distance_squared = square(to_start.x) + square(to_start.y);
    if (!(distance_squared.inf() > 0.0)) return;
    const Interval radicand = rim_squared_ / distance_squared - Interval(0.25);
    if (!(radicand.inf() > 0.0)) return;
    const Interval root = CGAL::sqrt(radicand);
    const Interval half(0.5);
    const Point middle{half * to_start.x, half * to_start.y};
    const Point side{-root * to_start.y, root * to_start.x};
    add_cover({middle.x - side.x, middle.y - side.y}, {middle.x + side.x, middle.y + side.y}, true);
}

// The arcs outside P. With the part on each edge's left, the crossing on the positive root
// enters it: there the rim's CCW tangent points left of the edge. See the header.
void AnnulusFilter2::add_walls(const Point& station) const
{
    const RimBox rim = rim_box(station);
    events_.clear();
    tangencies_.clear();
    for (const Edge* candidate : wall_candidates_) {
        const Edge& edge = *candidate;
        if (edge.max_x < rim.min_x || edge.min_x > rim.max_x
            || edge.max_y < rim.min_y || edge.min_y > rim.max_y) {
            continue;
        }
        const Point offset{Interval(edge.ax) - station.x, Interval(edge.ay) - station.y};
        const Interval side = edge.dx * offset.y - edge.dy * offset.x;
        const Interval discriminant = edge.squared_length * rim_squared_ - square(side);
        if (discriminant.sup() < 0.0) continue;
        const Interval along = edge.dx * offset.x + edge.dy * offset.y;
        const bool crossing = discriminant.inf() > 0.0;
        const Interval root =
            CGAL::sqrt(crossing ? discriminant : Interval(0.0, discriminant.sup()));
        const auto wall_point = [&](const Interval& parameter) {
            const Interval scale = parameter * edge.inverse_squared_length;
            return Point{offset.x + scale * edge.dx, offset.y + scale * edge.dy};
        };
        if (crossing) {
            for (const bool entering : {false, true}) {
                const Interval parameter = entering ? root - along : -root - along;
                if (parameter.sup() < 0.0 || parameter.inf() > edge.squared_length.sup()) continue;
                if (parameter.inf() <= 0.0 || parameter.sup() >= edge.squared_length.inf()) return;
                const Point point = wall_point(parameter);
                double low;
                double high;
                if (!key_of(point, low, high)) return;
                events_.push_back({low, high, entering, add_point(point)});
            }
        } else {
            const Interval parameter = Interval(-root.sup(), root.sup()) - along;
            if (parameter.sup() < 0.0 || parameter.inf() > edge.squared_length.sup()) continue;
            if (parameter.inf() <= 0.0 || parameter.sup() >= edge.squared_length.inf()) return;
            double low;
            double high;
            if (!key_of(wall_point(parameter), low, high) || low < 0.0) return;
            tangencies_.emplace_back(low, high);
        }
    }
    if (events_.empty()) return;
    std::sort(events_.begin(), events_.end(), [](const WallEvent& lhs, const WallEvent& rhs) {
        return lhs.low < rhs.low;
    });
    const std::size_t count = events_.size();
    for (std::size_t index = 0; index < count; ++index) {
        const WallEvent& event = events_[index];
        const WallEvent& next = events_[(index + 1) % count];
        if (event.entering == next.entering || event.low < 0.0 || event.high > full_turn_key_) return;
        if (index + 1 < count && event.high >= next.low) return;
    }
    const std::size_t first_wall_piece = pieces_.size();
    for (std::size_t index = 0; index < count; ++index) {
        const WallEvent& leaving = events_[index];
        if (leaving.entering) continue;
        const WallEvent& entering = events_[(index + 1) % count];
        if (leaving.high <= entering.low) {
            pieces_.push_back({leaving.high, entering.low, leaving.point, entering.point});
        } else {
            pieces_.push_back({leaving.high, full_turn_key_, leaving.point, -1});
            pieces_.push_back({0.0, entering.low, -1, entering.point});
        }
    }
    for (const auto& [low, high] : tangencies_) {
        for (std::size_t index = first_wall_piece; index < pieces_.size(); ++index) {
            if (low <= pieces_[index].high && high >= pieces_[index].low) {
                pieces_.resize(first_wall_piece);
                return;
            }
        }
    }
}

// The cover of one annulus. Returns true when the whole rim is certainly inside it.
bool AnnulusFilter2::add_annulus(const Point& station, const StoredAnnulus2& annulus) const
{
    const Point to_centre{Interval(annulus.cx) - station.x, Interval(annulus.cy) - station.y};
    const Interval distance_squared = square(to_centre.x) + square(to_centre.y);
    if (!(distance_squared.inf() > 0.0)) return false;
    const Interval inverse = Interval(1.0) / distance_squared;
    const Interval ratio = rim_squared_ * inverse;
    const Interval shifted = distance_squared + rim_squared_;
    const Interval half_inverse = Interval(0.5) * inverse;
    enum class Disk { NONE, WHOLE, ARC };

    // The outer disk: a cover only shrinks, so a near-tangency that may be no cover is none.
    const Interval outer_a = (shifted - annulus.outer_squared) * half_inverse;
    const Interval outer_q = ratio - square(outer_a);
    Disk outer;
    Interval outer_root(0.0);
    if (outer_q.sup() < 0.0) {
        if (!(outer_a.sup() < 0.0)) return false;
        outer = Disk::WHOLE;
    } else if (outer_q.inf() > 0.0) {
        outer = Disk::ARC;
        outer_root = CGAL::sqrt(outer_q);
    } else {
        if (!(outer_a.sup() < 0.0)) return false;
        outer = Disk::ARC;
        outer_root = CGAL::sqrt(Interval(0.0, outer_q.sup()));
    }

    // The inner disk grows out of the cover: a near-tangency keeps the largest arc it can
    // be, and a rim that may be inside the hole has no cover at all.
    Disk inner = Disk::NONE;
    Interval inner_a(0.0);
    Interval inner_root(0.0);
    if (annulus.has_inner) {
        inner_a = (shifted - annulus.inner_squared) * half_inverse;
        const Interval inner_q = ratio - square(inner_a);
        if (inner_q.sup() < 0.0) {
            if (!(inner_a.inf() > 0.0)) return false;
        } else if (inner_q.inf() > 0.0) {
            inner = Disk::ARC;
            inner_root = CGAL::sqrt(inner_q);
        } else {
            if (!(inner_a.inf() > 0.0)) return false;
            inner = Disk::ARC;
            inner_root = CGAL::sqrt(Interval(0.0, inner_q.sup()));
        }
    }

    const auto rim_point = [&to_centre](const Interval& a, const Interval& root) {
        return Point{a * to_centre.x - root * to_centre.y, a * to_centre.y + root * to_centre.x};
    };
    if (outer == Disk::WHOLE) {
        if (inner == Disk::NONE) return true;
        add_cover(rim_point(inner_a, inner_root), rim_point(inner_a, -inner_root), inner_a.sup() < 0.0);
        return false;
    }
    if (inner == Disk::NONE) {
        add_cover(rim_point(outer_a, -outer_root), rim_point(outer_a, outer_root), outer_a.inf() > 0.0);
        return false;
    }
    add_cover(rim_point(outer_a, -outer_root), rim_point(inner_a, -inner_root), true);
    add_cover(rim_point(inner_a, inner_root), rim_point(outer_a, outer_root), true);
    return false;
}

// The covers `tier` adds to the tiers before it. Returns true when the whole rim is covered.
bool AnnulusFilter2::add_tier(
    int tier, const AnnulusStore2& store, AnnulusViewRange2 range, const Point& station,
    const AnnulusSearchBox2& read_box, bool& gathered) const
{
    const AnnulusViewRange2 one = last_of_range(range, DECISION_TIER_LAST_ONE);
    const AnnulusViewRange2 sixteen = last_of_range(range, DECISION_TIER_LAST_SIXTEEN);
    if (tier == 0) {
        for (std::size_t id = one.first; id < one.end; ++id) {
            if (add_annulus(station, store[id])) return true;
        }
        return false;
    }
    if (tier == 1) {
        for (std::size_t id = sixteen.first; id < one.first; ++id) {
            if (add_annulus(station, store[id])) return true;
        }
        return false;
    }
    if (!gathered) {
        range_candidates_.clear();
        store.query(read_box, {range.first, sixteen.first}, range_candidates_);
        gathered = true;
    }
    const RimBox rim = rim_box(station);
    const AnnulusSearchBox2 rim_search{rim.min_x, rim.min_y, rim.max_x, rim.max_y};
    for (const std::uint32_t id : range_candidates_) {
        const StoredAnnulus2& annulus = store[id];
        if (search_boxes_meet(annulus.box, rim_search) && add_annulus(station, annulus)) return true;
    }
    return false;
}

bool AnnulusFilter2::key_length_within(double length) const
{
    const Interval bound = Interval(4.0) * tool_radius_ * Interval(length);
    return length < half_turn_key_ && bound.sup() < cap_squared_;
}

// The two exact legs, on enclosures: an arc is within when its orientation is certainly
// positive and its chord squared certainly below the cap; and the key-length bound for
// two ends too close to orient.
bool AnnulusFilter2::arc_within(int from, int to, double length) const
{
    const Point& start = points_[static_cast<std::size_t>(from)];
    const Point& end = points_[static_cast<std::size_t>(to)];
    const Interval turn = start.x * end.y - start.y * end.x;
    if (turn.inf() > 0.0) {
        const Interval chord = square(end.x - start.x) + square(end.y - start.y);
        chord_bound_ = std::max(chord_bound_, chord.sup());
        return chord.sup() < cap_squared_;
    }
    if (key_length_within(length)) {
        chord_bound_ = std::max(chord_bound_, (Interval(4.0) * tool_radius_ * Interval(length)).sup());
        return true;
    }
    return false;
}

// Whether every arc of the rim left by the covers is within the cap. No cover at all leaves
// the whole rim, a closed loop.
bool AnnulusFilter2::certify() const
{
    chord_bound_ = 0.0;
    if (pieces_.empty()) return false;
    sorted_ = pieces_;
    std::sort(sorted_.begin(), sorted_.end(), [](const Piece& lhs, const Piece& rhs) {
        return lhs.low < rhs.low;
    });
    double reach = sorted_.front().high;
    int reach_point = sorted_.front().high_point;
    for (std::size_t index = 1; index < sorted_.size(); ++index) {
        const Piece& piece = sorted_[index];
        if (piece.low > reach) {
            if (reach_point < 0 || piece.low_point < 0) return false;
            if (!arc_within(reach_point, piece.low_point, piece.low - reach)) return false;
        }
        if (piece.high > reach) {
            reach = piece.high;
            reach_point = piece.high_point;
        }
    }
    const double first_low = sorted_.front().low;
    if (reach < full_turn_key_ || first_low > 0.0) {
        const int first_point = sorted_.front().low_point;
        if (reach_point < 0 || first_point < 0) return false;
        if (!arc_within(reach_point, first_point, (full_turn_key_ - reach) + first_low)) return false;
    }
    return true;
}

void AnnulusFilter2::certify_circle(
    const AnnulusStore2& store, AnnulusViewRange2 range,
    const AnnulusFilterCircle2& circle, double unit_t, bool unit_reflected,
    const std::vector<AnnulusFilterRotation2>& rotations, double cap_chord_ratio,
    std::vector<AnnulusFilterStation2>& out) const
{
    const CGAL::Protect_FPU_rounding<true> rounding;
    cap_squared_ = (Interval(cap_chord_ratio) * rim_squared_).inf();
    const Interval centre_x(circle.cx);
    const Interval centre_y(circle.cy);
    const Interval guide(circle.guide_radius);

    // phase = rho * ((1 - t^2), 2 t) / (1 + t^2), negated when reflected: rational_unit_along.
    const Interval t(unit_t);
    const Interval t_squared = square(t);
    const Interval inverse = Interval(1.0) / (Interval(1.0) + t_squared);
    Point unit{(Interval(1.0) - t_squared) * inverse, Interval(2.0) * t * inverse};
    if (unit_reflected) unit = {-unit.x, -unit.y};
    const Point phase{guide * unit.x, guide * unit.y};
    const Interval rear_scale = tool_radius_ / guide;

    // Every rim of the circle, the entry's included, lies in the read disk.
    const Interval entry_x = Interval(circle.entry_x) - centre_x;
    const Interval entry_y = Interval(circle.entry_y) - centre_y;
    const Interval reach =
        CGAL::sqrt(CGAL::max(square(guide), square(entry_x) + square(entry_y))) + tool_radius_;
    const Interval read_min_x = centre_x - reach;
    const Interval read_min_y = centre_y - reach;
    const Interval read_max_x = centre_x + reach;
    const Interval read_max_y = centre_y + reach;
    const AnnulusSearchBox2 read_box{
        read_min_x.inf(), read_min_y.inf(), read_max_x.sup(), read_max_y.sup()};
    wall_candidates_.clear();
    for (const Edge& edge : edges_) {
        if (edge.max_x < read_box.min_x || edge.min_x > read_box.max_x
            || edge.max_y < read_box.min_y || edge.min_y > read_box.max_y) {
            continue;
        }
        wall_candidates_.push_back(&edge);
    }

    out.resize(rotations.size());
    bool gathered = false;
    for (std::size_t k = 0; k < rotations.size(); ++k) {
        const AnnulusFilterRotation2& rotation = rotations[k];
        points_.clear();
        pieces_.clear();
        Point station{Interval(circle.entry_x), Interval(circle.entry_y)};
        if (!rotation.entry) {
            const Point offset{
                rotation.cosine * phase.x - rotation.sine * phase.y,
                rotation.cosine * phase.y + rotation.sine * phase.x};
            station = {centre_x + offset.x, centre_y + offset.y};
            add_own_prefix(phase, offset, rear_scale);
        }
        add_walls(station);
        AnnulusFilterStation2 result{AnnulusFilterTier2::UNCERTIFIED, 0.0};
        for (int tier = 0; tier <= static_cast<int>(AnnulusFilterTier2::RANGE); ++tier) {
            if (add_tier(tier, store, range, station, read_box, gathered)) {
                result = {static_cast<AnnulusFilterTier2>(tier), 0.0};
                break;
            }
            if (certify()) {
                result = {static_cast<AnnulusFilterTier2>(tier), chord_bound_};
                break;
            }
        }
        out[k] = result;
    }
}
