#pragma once

// Native declared-boundary clearance oracle: the exact sign of
// dist(element, dP) - r for binary64 machining circles and exact segment and
// arc elements, against the declared (Pythagorean) boundary
// (benchmarks/held_declared_clearance.py is the Python oracle it reproduces).
//
// DECIDING: every verdict is exact. Primitives are pruned only by the AABB
// index (primitive_index_2.h); every remaining pair is decided by exact,
// lazily filtered Epeck / Sqrt_extension predicates: a crossing test, then the
// closed-form distance candidates (candidates_2.h), or, for a declared owner
// whose lemma hypotheses hold, structurally (owner_lemma_2.h).
//
// REPORTING: per element, the worst sign, then among the primitives decided
// with that sign the one with the smallest reported distance, then the
// smallest index (the Python oracle's selection, on the same binary64 views,
// reporting_2.h). The depth is r - distance for GOUGE, 0 otherwise. Reported
// values never feed back into a verdict.

#include "held_declared_clearance_2/elements_2.h"
#include "held_declared_clearance_2/primitive_index_2.h"
#include "held_declared_clearance_2/primitives_2.h"

#include <cstddef>
#include <memory>
#include <vector>

namespace compas_cgal::held_declared_clearance {

enum class Clearance : int { Gouge = -1, Touch = 0, Clear = 1 };

/// One element's verdict; `primitive` is -1 when every primitive is CLEAR.
struct ClearanceVerdict {
    Clearance clearance;
    std::ptrdiff_t primitive;
    double depth;
};

/// Reported counters of the work done (never a decision input).
struct ClearanceCounters {
    std::size_t elements = 0;
    std::size_t pairs = 0;          // primitive pairs that survived the AABB pruning
    std::size_t owner_touches = 0;  // pairs decided TOUCH by an owner lemma
    std::size_t crossings = 0;      // pairs decided GOUGE by an exact crossing
};

/// Scratch storage reused across elements of one batch.
struct ClearanceWorkspace {
    std::vector<std::size_t> near;
    ClearanceCounters counters;
};

class DeclaredClearance2 {
public:
    /// Index the declared primitives (primitive i has index i).
    ///
    /// Raises:
    ///     NonPositiveClearanceRadiusError: the tool radius is not a positive binary64.
    ///     InvalidDeclaredPrimitiveError: no primitives, or a primitive's index is not its position.
    [[nodiscard]] static std::unique_ptr<DeclaredClearance2> build(std::vector<Primitive> primitives, const FT& tool_radius);

    /// Exact sign(dist(circle, dP) - r).
    ///
    /// Raises:
    ///     NonPositiveClearanceRadiusError: the circle radius is not positive.
    [[nodiscard]] ClearanceVerdict circle(const CircleElement& element, ClearanceWorkspace& workspace) const;
    /// Exact sign(dist(segment, dP) - r); `owners` are declared primitives the
    /// element may offset (owner_lemma_2.h; each is still checked).
    [[nodiscard]] ClearanceVerdict segment(const SegmentElement& element, const std::vector<std::size_t>& owners,
                                           ClearanceWorkspace& workspace) const;
    /// Exact sign(dist(arc, dP) - r), as segment().
    [[nodiscard]] ClearanceVerdict arc(const ArcElement& element, const std::vector<std::size_t>& owners,
                                       ClearanceWorkspace& workspace) const;

    [[nodiscard]] std::size_t primitive_count() const { return primitives_.size(); }
    [[nodiscard]] const Primitive& primitive(std::size_t index) const { return primitives_.at(index); }
    [[nodiscard]] const FT& tool_radius() const { return tool_radius_; }

private:
    DeclaredClearance2(std::vector<Primitive> primitives, const FT& tool_radius, double reach);

    std::vector<Primitive> primitives_;
    FT tool_radius_;
    double reach_;  // the tool radius, an exact binary64
    PrimitiveIndex index_;
};

}  // namespace compas_cgal::held_declared_clearance
