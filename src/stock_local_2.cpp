#include "stock_local_2.h"

#include <algorithm>
#include <functional>
#include <string>
#include <variant>

#include <CGAL/Arr_observer.h>
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/Arrangement_2.h>

// ----------------------------------------------------------------------------
// LOCAL DEPLETION -- removing a bounded region from the stock by editing the
// Gps's OWN arrangement, instead of overlaying the whole stock with the region.
//
// `General_polygon_set_2::difference` rebuilds the entire arrangement: it
// sweeps the stock against the region and allocates a fresh vertex, halfedge
// and face record for every feature, whether or not that feature is anywhere
// near the removed material. This is the depletion twin of the engagement
// query's local zone walk (docs/engagement_zone_query.md,
// engagement_2.cpp::engaged_arcs_zone) -- with the difference that a query only
// READS the arrangement, whereas a removal must edit it and then restore every
// invariant the boolean engine maintains.
//
// SOUNDNESS. Write S for the stock point set, R for the removal region, and A
// for the arrangement the Gps holds (its faces carry contained() = "in S").
//
//   I1  After inserting the x-monotone arcs of dR into A, every face of the
//       refined arrangement A' still satisfies f subset S <=> f.contained().
//       Insertion never moves a point across a boundary; it only SPLITS faces,
//       and both pieces of a split are subsets of one original face. CGAL
//       creates the split face with a default-constructed DCEL record
//       (Arrangement_on_surface_2_impl.h:2927 -- a bare `_dcel().new_face()`
//       with no `assign`), so contained() would silently read FALSE on every
//       piece dR cuts off; ContainmentPropagator restores the invariant through
//       the `after_split_face` notification. On a bounded planar topology that
//       notification is the ONLY face-creating event insertion can raise -- the
//       remaining `new_face()` sites in vendored CGAL 6.0.1 are DCEL
//       initialisation and the unbounded-planar topology traits.
//
//   I2  dR is now part of A', so no face of A' straddles dR: each face lies
//       entirely inside R or entirely outside it.
//
//   I3  The flood marks exactly the faces inside R. int(R) is open, and
//       dR ^ int(R) is empty, so a path between two faces inside R crosses only
//       NON-boundary edges: flooding from a seed in each component of int(R),
//       never crossing an edge that lies on dR, reaches every face inside R and
//       no face outside it. Setting contained(false) on exactly those faces
//       turns the characteristic function of S into that of S \ R (I1 + I2).
//
//   I4  An edge is redundant (its two incident faces agree on contained()) only
//       if the marking changed one of their flags, i.e. only if it bounds a
//       flooded face -- every other face kept its flag, and the Gps invariant
//       held before. So the candidate set is exactly the boundary of the
//       flooded faces, and removing the redundant ones restores the canonical
//       Gps representation locally. Containment values do not change during
//       removal (a merge joins two faces that already agreed), so the predicate
//       is stable across the loop.
//
// EXACTNESS. Every decision above is an exact predicate on exact quantities:
// point location and curve insertion are CGAL's own exact machinery, the
// boundary test is exact rational equality of a circle's centre and squared
// radius, and the redundancy test is a comparison of two booleans. No epsilon,
// no tolerance, no ulp-nudge anywhere in the path.
// ----------------------------------------------------------------------------

namespace {

using Arr = Gps::Arrangement_2;
using PL = CGAL::Arr_walk_along_line_point_location<Arr>;

// Restores I1 across the face splits incremental insertion causes. Both pieces
// of a split are subsets of one original face and therefore carry its
// containment, so copying is the complete rule and it does not matter which of
// the two CGAL hands back as `new_face`.
class ContainmentPropagator : public CGAL::Arr_observer<Arr> {
public:
    explicit ContainmentPropagator(Arr& arrangement)
        : CGAL::Arr_observer<Arr>(arrangement)
    {
    }

    void after_split_face(Face_handle face, Face_handle new_face, bool) override
    {
        new_face->set_contained(face->contained());
    }
};

// Does `curve` lie on the region boundary? Exact: a circular arc stores its
// supporting circle's centre and SQUARED radius as rationals, and both survive
// every split, so this is rational equality -- not a distance test.
bool lies_on_region_boundary(const GpsXCurve& curve,
                             const std::vector<std::pair<EPoint, Epeck::FT>>& circles)
{
    if (!curve.is_circular()) {
        return false;
    }
    const ECircle support = curve.supporting_circle();
    for (const auto& [center, squared_radius] : circles) {
        if (support.center() == center && support.squared_radius() == squared_radius) {
            return true;
        }
    }
    return false;
}

// The face a strictly-interior seed belongs to. By I2 no face straddles dR, so
// a seed that lands on an edge or a vertex instead of in a face interior is
// still surrounded by faces inside R: any incident face will do.
Arr::Face_handle seed_face(Arr& arrangement, const PL& point_location, const GpsPoint& seed)
{
    const auto located = point_location.locate(seed);
    if (const auto* face = std::get_if<Arr::Face_const_handle>(&located)) {
        return arrangement.non_const_handle(*face);
    }
    if (const auto* halfedge = std::get_if<Arr::Halfedge_const_handle>(&located)) {
        return arrangement.non_const_handle(*halfedge)->face();
    }
    const Arr::Vertex_handle vertex =
        arrangement.non_const_handle(std::get<Arr::Vertex_const_handle>(located));
    return vertex->is_isolated() ? vertex->face() : vertex->incident_halfedges()->face();
}

// Flood (I3). Fills `inside` with every face of the refined arrangement that
// lies inside the region and `incident` with every halfedge on their boundaries
// -- the redundancy candidate set of I4, collected here because the walk that
// finds the faces already visits exactly those halfedges.
//
// The Gps face `visited` bit is the flood mark. That is the flag's purpose --
// General_polygon_set_2's own traversals use it the same way and clear it with
// `_reset_faces` afterwards; this operation upholds the same contract, clearing
// every bit it sets, including on the throwing path.
void flood_region_faces(Arr& arrangement,
                        const LocalRemovalRegion& region,
                        const PL& point_location,
                        std::vector<Arr::Face_handle>& inside,
                        std::vector<Arr::Halfedge_handle>& incident)
{
    std::vector<Arr::Face_handle> pending;
    for (const GpsPoint& seed : region.seeds) {
        const Arr::Face_handle start = seed_face(arrangement, point_location, seed);
        if (start->visited()) {
            continue;
        }
        start->set_visited(true);
        inside.push_back(start);
        pending.push_back(start);
    }

    const auto walk_ccb = [&](Arr::Ccb_halfedge_circulator start) {
        Arr::Ccb_halfedge_circulator curr = start;
        do {
            incident.push_back(curr);
            if (lies_on_region_boundary(curr->curve(), region.boundary_circles)) {
                continue;   // dR separates inside from outside: never cross it
            }
            const Arr::Face_handle neighbour = curr->twin()->face();
            if (neighbour->visited()) {
                continue;
            }
            neighbour->set_visited(true);
            inside.push_back(neighbour);
            pending.push_back(neighbour);
        } while (++curr != start);
    };

    while (!pending.empty()) {
        const Arr::Face_handle face = pending.back();
        pending.pop_back();
        // A bounded region encloses no unbounded face, so reaching one means the
        // flood escaped: fail loudly instead of depleting the wrong material.
        if (!face->has_outer_ccb()) {
            throw LocalDepletionEscapedError(
                "local depletion flood reached the unbounded face; the removal "
                "region's boundary does not separate its interior.");
        }
        walk_ccb(face->outer_ccb());
        for (auto hole = face->inner_ccbs_begin(); hole != face->inner_ccbs_end(); ++hole) {
            walk_ccb(*hole);
        }
    }
}

// A stable per-EDGE identity for the two halfedges of a twin pair, so the
// candidate set can be de-duplicated: the flood sees an interior edge from both
// of its faces.
const void* edge_address(Arr::Halfedge_handle halfedge)
{
    const Arr::Halfedge* self = &*halfedge;
    const Arr::Halfedge* twin = &*halfedge->twin();
    return std::less<const void*>{}(self, twin) ? static_cast<const void*>(self)
                                                : static_cast<const void*>(twin);
}

// Restore the Gps invariant on the candidate set (I4). Removing one edge never
// invalidates another candidate's handle: `Arrangement_2::remove_edge`
// deallocates only that twin pair, its end vertices if they become isolated,
// and one of the two merged faces -- it does not merge the surviving curves.
void remove_redundant_local_edges(Arr& arrangement,
                                  std::vector<Arr::Halfedge_handle>& candidates)
{
    const auto by_address = [](Arr::Halfedge_handle lhs, Arr::Halfedge_handle rhs) {
        return std::less<const void*>{}(edge_address(lhs), edge_address(rhs));
    };
    const auto same_edge = [](Arr::Halfedge_handle lhs, Arr::Halfedge_handle rhs) {
        return edge_address(lhs) == edge_address(rhs);
    };
    std::sort(candidates.begin(), candidates.end(), by_address);
    candidates.erase(std::unique(candidates.begin(), candidates.end(), same_edge),
                     candidates.end());

    for (const Arr::Halfedge_handle edge : candidates) {
        if (edge->face()->contained() == edge->twin()->face()->contained()) {
            arrangement.remove_edge(edge, true, true);
        }
    }
}

// Traversal orientations named by which side of the circle keeps its material.
// Clockwise puts a circle's EXTERIOR on the left, counterclockwise its interior.
constexpr CGAL::Orientation MATERIAL_OUTSIDE = CGAL::CLOCKWISE;
constexpr CGAL::Orientation MATERIAL_INSIDE = CGAL::COUNTERCLOCKWISE;

// The two x-monotone arcs of one boundary circle, split at its x-extreme
// vertical-tangency points -- the same subdivision disk_polygon uses -- but
// traversed so the side that SURVIVES the removal lies to the left of each arc.
//
// The orientation is not cosmetic. General_polygon_set_2's representation
// invariant (Gps_on_surface_base_2::_is_valid) requires every stored curve to
// carry the contained side on its left, and CGAL::insert stores the curve
// exactly as handed over -- Split_2 preserves that direction through every
// intersection. A removal region's interior is by construction NOT contained, so
// its OUTER circle is traversed clockwise: the mirror of disk_polygon, which
// builds the counterclockwise circle of a region being ADDED. An annulus's INNER
// circle bounds material the removal keeps, so it is traversed counterclockwise.
// Arcs whose surviving side turns out to be void anyway are redundant and are
// dropped in the last phase, so the rule only has to be right where it shows.
void append_boundary_circle(const EPoint& center,
                            const Epeck::FT& radius,
                            CGAL::Orientation traversal,
                            LocalRemovalRegion& region)
{
    const Epeck::FT squared_radius = radius * radius;
    const ECircle circle(center, squared_radius, traversal);
    const GpsPoint leftmost(center.x() - radius, center.y());
    const GpsPoint rightmost(center.x() + radius, center.y());
    // Clockwise is decreasing angle, so it runs rightmost -> bottom -> leftmost
    // and then leftmost -> top -> rightmost; counterclockwise is the reverse.
    const GpsPoint& first = (traversal == CGAL::CLOCKWISE) ? rightmost : leftmost;
    const GpsPoint& second = (traversal == CGAL::CLOCKWISE) ? leftmost : rightmost;
    region.boundary_curves.emplace_back(circle, first, second, traversal);
    region.boundary_curves.emplace_back(circle, second, first, traversal);
    region.boundary_circles.emplace_back(center, squared_radius);
}

} // namespace

LocalRemovalRegion local_disk_region(const EPoint& center, const Epeck::FT& radius)
{
    LocalRemovalRegion region;
    append_boundary_circle(center, radius, MATERIAL_OUTSIDE, region);
    // The centre is strictly inside a disk of positive radius, and a disk's
    // interior is connected, so one seed is enough.
    region.seeds.emplace_back(center.x(), center.y());
    return region;
}

LocalRemovalRegion local_annulus_region(const EPoint& center,
                                        const Epeck::FT& inner_radius,
                                        const Epeck::FT& outer_radius)
{
    // A guide no wider than the tool sweeps a filled disk, with no hole to punch:
    // a structural branch on an exact sign, not an epsilon. A NEGATIVE radius is a
    // contract violation, rejected upstream at the one boundary that owns the
    // error model (Stock2::subtract_annulus_exact_local); folding it in here keeps
    // this builder total rather than half-defined.
    if (CGAL::sign(inner_radius) != CGAL::POSITIVE) {
        return local_disk_region(center, outer_radius);
    }
    LocalRemovalRegion region;
    append_boundary_circle(center, outer_radius, MATERIAL_OUTSIDE, region);
    append_boundary_circle(center, inner_radius, MATERIAL_INSIDE, region);
    // Mid-radius on the positive x axis: an exact rational strictly between the
    // two radii, and an open annulus is connected, so one seed is enough.
    const Epeck::FT mid_radius = (inner_radius + outer_radius) / Epeck::FT(2);
    region.seeds.emplace_back(center.x() + mid_radius, center.y());
    return region;
}

void subtract_region_local(Gps& set, const LocalRemovalRegion& region)
{
    Arr& arrangement = set.arrangement();

    {
        // Scoped so the propagator is detached before the marking phase: it must
        // observe insertion (I1) and nothing else.
        ContainmentPropagator propagator(arrangement);
        PL point_location(arrangement);
        for (const GpsXCurve& curve : region.boundary_curves) {
            CGAL::insert(arrangement, curve, point_location);
        }
    }

    const PL point_location(arrangement);
    std::vector<Arr::Face_handle> inside;
    std::vector<Arr::Halfedge_handle> incident;
    try {
        flood_region_faces(arrangement, region, point_location, inside, incident);
    } catch (...) {
        for (const Arr::Face_handle face : inside) {
            face->set_visited(false);
        }
        throw;
    }

    // Clear the flood marks BEFORE any removal: edge removal merges faces, so
    // these handles stop being dereferenceable once the next phase starts.
    for (const Arr::Face_handle face : inside) {
        face->set_visited(false);
        face->set_contained(false);
    }

    remove_redundant_local_edges(arrangement, incident);
}
