#pragma once

// The cover-counter sweep along one flank piece (a CAPSULE_FLANK of a merged
// connector rectangle, or a target flank checked by the target-edge rule).
//
// A flank piece is the open chart interval (t_lo, t_hi) of a line support
// (lines_2.h). Its events are the crossings with
//   * the boundary circles of candidate circle regions (two tower roots, or a
//     tangency, line_crossings_2.h),
//   * the lines of candidate rectangles (one tower root each), and
//   * for a swept piece, the target cut supports (vertex circles, target flank
//     lines), which flip nothing;
// restricted exactly to (t_lo, t_hi). All event parameters lie in towers over
// the piece's own field Q(sqrt(L)), so every comparison is the tower sign rule
// (no cross-base case on a line).
//
// Coincident rectangle lines (P5/P6): a candidate rectangle with a line equal to
// the piece's support has that line's side identically zero along the piece. If
// the rectangle lies on the piece's checked side there, that half-plane counts
// as satisfied for the whole piece (the rectangle then covers wherever its other
// three sides hold: a coincident support); otherwise the rectangle can never
// cover the checked side and is excluded. This is what makes a connector flank
// on edge e cover the target flank of edge e exactly.
//
// The start state is decided at a rational chart parameter strictly inside the
// first sub-edge (exact P4 at a one-root point); flips follow in order.

#include "coverage_certificate_2/counters_2.h"
#include "coverage_certificate_2/line_crossings_2.h"
#include "coverage_certificate_2/lines_2.h"
#include "coverage_certificate_2/regions_2.h"
#include "coverage_certificate_2/tower_2.h"
#include "coverage_certificate_2/vertex_events_2.h"

#include <CGAL/Bbox_2.h>

#include <cstdint>
#include <vector>

namespace compas_cgal::coverage_certificate {

/// One flank piece.
struct LinePiece {
    std::uint32_t line;     ///< support id in the LineTable
    FT t_lo, t_hi;          ///< open chart interval
    bool target;            ///< target flank (else a connector flank)
    std::uint32_t owner;    ///< rectangle index (swept or target list)
    std::uint32_t report_index;  ///< connector index, or target member index
    CGAL::Sign check_side;  ///< side of the line's normal whose cover is checked
    CGAL::Bbox_2 box;       ///< certified box of the segment, widened by LINE_PIECE_BOX_MARGIN
};

/// Superset filter margin of a flank piece's candidate box, mm (2^-10). A flank
/// along an axis has a degenerate box; widening it lets P7 place its uncovered
/// point beside the piece while the candidate set stays complete for it. Policy
/// of the box filter only: every decision is exact.
inline constexpr double LINE_PIECE_BOX_MARGIN = 1.0 / 1024.0;

/// Candidates of one flank piece.
struct LineCandidates {
    std::vector<std::uint32_t> regions;         ///< swept circle regions, sorted
    std::vector<std::uint32_t> rects;           ///< swept rectangles (not the owner), sorted
    std::vector<std::uint32_t> target_circles;  ///< target circle ids (cuts), sorted; empty for a target piece
    std::vector<std::uint32_t> target_lines;    ///< target flank lines (cuts), sorted; empty for a target piece
    std::vector<std::uint32_t> target_vertices; ///< one-root vertex circles (cuts), sorted
};

/// One event on a flank piece.
struct LineEvent {
    std::uint8_t source;        ///< 0 circle crossing, 1 line crossing, 2 vertex-circle root (a cut, flips nothing)
    std::uint32_t record;       ///< index into the sweep's circle or line records
    std::uint8_t branch;
    std::uint8_t multiplicity;  ///< 1 transversal, 2 tangent
    std::uint32_t other;        ///< circle id or line id
    Interval t;
};

/// The arrangement of one flank piece after the sweep.
class LineSweep {
public:
    std::vector<LineEvent> events;            ///< strictly inside (t_lo, t_hi), ascending t
    std::vector<std::uint32_t> group_starts;  ///< vertex v spans events [group_starts[v], group_starts[v + 1])
    /// Sub-edge e runs from vertex e - 1 (t_lo for e = 0) to vertex e (t_hi for the last).
    std::vector<std::uint8_t> edge_exposed;
    std::vector<std::uint32_t> regions;   ///< counted circle regions
    std::vector<std::uint32_t> rects;     ///< counted rectangles
    std::vector<std::uint8_t> rect_coincident_slot;  ///< per counted rect: slot of its coincident covering line, or 4
    std::vector<std::uint32_t> excluded_rects;       ///< coincident on the unchecked side: never cover

    [[nodiscard]] std::size_t vertex_count() const { return group_starts.empty() ? 0 : group_starts.size() - 1; }

    /// Exact chart parameter of a source-0/1 event (materialised on demand).
    [[nodiscard]] Tower exact_t(const LineEvent& event, DecisionCounters& counters);

    /// The exact roots of a source-2 event (materialised on demand).
    [[nodiscard]] const LineVertexExact& exact_vertex(const LineEvent& event, DecisionCounters& counters);

    /// Exact sign(t_event - q) for rational q, any source.
    [[nodiscard]] CGAL::Sign minus_rational(const LineEvent& event, const FT& q, DecisionCounters& counters);

    /// Exact order of two events, any sources (towers by the tower sign rule;
    /// a vertex root by provenance equality and rational separation).
    [[nodiscard]] CGAL::Comparison_result compare_events(const LineEvent& a, const LineEvent& b, DecisionCounters& counters);

    std::vector<CircleLineRecord> circle_records;
    std::vector<LineLineRecord> line_records;
    std::vector<LineVertexRecord> vertex_records;
    const std::vector<VertexCircle>* vertices = nullptr;
    const CircleRegionSet* set = nullptr;
    const LineTable* lines = nullptr;
    const LinePiece* piece = nullptr;
};

/// Run the sweep on one flank piece.
///
/// Raises:
///     CoverageCertificateNativeError: the start witness lies on a support, or
///         a filter enclosure misses its exact value.
[[nodiscard]] LineSweep sweep_line(const LinePiece& piece, const CircleRegionSet& set, const LineTable& lines,
                                   const std::vector<RectRegion>& rects, const LineCandidates& candidates, FilterMode mode,
                                   DecisionCounters& counters, const std::vector<VertexCircle>* vertices = nullptr);

/// One end of a sub-edge: a rational chart bound, an event's exact tower, or a
/// vertex-circle root (sign oracle, vertex_events_2.h).
struct ChartBound {
    bool rational;
    FT q;
    Tower t;
    const LineVertexExact* vertex = nullptr;  ///< set for a vertex-circle root
    std::uint8_t branch = 0;
    Interval vertex_t = Interval(0);

    [[nodiscard]] Interval interval() const
    {
        if (vertex != nullptr) return vertex_t;
        return rational ? interval_of(q) : interval_of(t);
    }
    [[nodiscard]] Tower tower() const { return rational ? Tower::rational(q) : t; }
    /// Exact sign(bound - q).
    [[nodiscard]] CGAL::Sign minus(const FT& value) const;
};

/// The two ends of sub-edge e of a swept line.
[[nodiscard]] std::pair<ChartBound, ChartBound> sub_edge_bounds(LineSweep& sweep, std::size_t edge,
                                                                DecisionCounters& counters);

/// P3 on a line chart: a short rational strictly between lo < hi (exact).
[[nodiscard]] FT rational_between(const ChartBound& lo, const ChartBound& hi, DecisionCounters& counters);

/// Exact closed side of a one-root point relative to a circle region (P4).
[[nodiscard]] Side membership(const CircleRegion& region, const ExactPoint& p);

}  // namespace compas_cgal::coverage_certificate
