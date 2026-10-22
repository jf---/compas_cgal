#pragma once

#include <CGAL/Arr_accessor.h>
#include <CGAL/enum.h>

#include <array>
#include <cstddef>
#include <stdexcept>
#include <vector>

// Restore the planar one-outer-CCB invariant on a Gps arrangement.
//
// CGAL 6.0.1 `Gps_on_surface_base_2::_remove_redundant_edges` (and CGAL main as
// of 2026-10) tags the CCBs of a merged face from twin evidence only: a cycle
// whose twin lies on an inner CCB becomes outer, and its last loop files every
// still-undecided cycle as OUTER on the stated assumption that it is "the only
// ccb of the face" (Gps_on_surface_base_2.h:1110-1119; it never consults
// face_outer_ccb_set there). That assumption fails when one outer CCB of the
// pre-merge faces splits into two cycles neither of which has inner-CCB
// evidence -- e.g. an island of material left inside a removed region whose
// boundary touched it only through edges that turned redundant. The result is
// a bounded face carrying TWO outer CCBs, one of which is really a hole.
//
// In a bounded planar arrangement every bounded face has exactly one outer CCB
// (a boundary pinched at a vertex is still ONE CCB whose next-chain visits the
// vertex twice), so this state exists only through the misfile.
// `Arrangement_2::is_valid` and `General_polygon_set_2`'s own operations
// tolerate it, but every consumer that relies on the planar invariant does not:
// `Face::outer_ccb()` (precondition compiled out in release), the face-split
// hole relocation of incremental insertion, and `Arrangement_2::remove_edge`
// (which then calls `erase_inner_ccb(nullptr)` on a two-halfedge "outer" cycle
// -- SIGSEGV). Call this right after the Gps boolean whose arrangement such a
// consumer reads.
//
// The repair is exact. Every curve of the arrangement is x-monotone, so each
// cycle's xy-lexicographically smallest point is one of its vertices. A bounded
// face's true outer boundary encloses every hole, holes are disjoint from it
// (a touching hole would be the same CCB), and the closed region it bounds has a
// unique xy-minimum lying on it: so the outer cycle is the one whose smallest
// vertex is strictly smallest, decided by the traits' exact `Compare_xy_2`
// (cross-root comparison of one-root points is exact and supported).
//
// Cost: one pass over the faces reading a container size, plus an exact
// comparison per vertex of a misfiled face only.

// The three error types below are templated on a module tag (`Domain`). Each
// nanobind module that runs the repair registers its own instantiation as its
// own named Python exception. One shared C++ type registered in two modules
// would let whichever module registered last translate the other's throws:
// translators live in nanobind's shared internals, and hidden-visibility type
// identity across shared objects falls back to name comparison on this target.

// Two distinct CCBs of one face share their xy-smallest vertex. Distinct CCBs are
// distinct connected components of the edge graph and can never share a vertex,
// so this is a corrupted arrangement, reported rather than guessed around.
template <typename Domain>
class AmbiguousOuterCcbError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// A bounded face with no outer CCB: impossible in a bounded planar arrangement
// (and never produced by the misfile above), so the arrangement is corrupt.
template <typename Domain>
class BoundedFaceWithoutOuterCcbError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// The unbounded face carries an outer CCB. Bounded planar topology gives it none,
// and the misfile above cannot produce one (CGAL forces face_outer_ccb_set for
// the unbounded set, and the twin of an inner-CCB halfedge faces a bounded
// face), so this is an unmodelled corruption, reported rather than repaired.
template <typename Domain>
class UnboundedFaceOuterCcbError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

namespace gps_ccb_repair_detail {

// The xy-smallest vertex of the cycle through `start`, by exact comparison.
template <typename Arr, typename CompareXy>
typename Arr::Vertex_handle lowest_vertex(typename Arr::Ccb_halfedge_circulator start, const CompareXy& compare_xy)
{
    typename Arr::Vertex_handle best = start->source();
    typename Arr::Ccb_halfedge_circulator curr = start;
    do {
        if (compare_xy(curr->source()->point(), best->point()) == CGAL::SMALLER) {
            best = curr->source();
        }
    } while (++curr != start);
    return best;
}

// Move the cycle through `start` from the outer-CCB list of `face` to its
// inner-CCB list, mirroring the DCEL bookkeeping CGAL's own Gps repair performs
// (Gps_on_surface_base_2.h:1226-1241: new Inner_ccb record, every halfedge
// re-pointed, record attached to the face, the emptied Outer_ccb record freed).
template <typename Arr>
void refile_as_inner(Arr& arrangement, typename Arr::Face_handle face, typename Arr::Halfedge_handle start)
{
    using Dcel = typename Arr::Dcel;
    CGAL::Arr_accessor<Arr> accessor(arrangement);
    typename Dcel::Face* dcel_face = &(*face);
    typename Dcel::Halfedge* first = &(*start);
    typename Dcel::Outer_ccb* outer = first->outer_ccb();
    typename Dcel::Inner_ccb* inner = accessor.new_inner_ccb();

    dcel_face->erase_outer_ccb(outer);
    typename Dcel::Halfedge* curr = first;
    do {
        curr->set_inner_ccb(inner);
        curr = curr->next();
    } while (curr != first);
    dcel_face->add_inner_ccb(inner, first);
    inner->set_halfedge(first);
    inner->set_face(dcel_face);
    accessor.delete_outer_ccbs(std::array<typename Dcel::Outer_ccb*, 1>{outer});
}

}  // namespace gps_ccb_repair_detail

// Re-file every misfiled outer CCB of `set`'s arrangement as an inner CCB.
// Returns the number of CCBs re-filed (0 on an already planar arrangement).
//
// `Domain` is the calling module's tag; it selects which instantiation of the
// three errors is thrown (see above).
//
// Raises:
//     AmbiguousOuterCcbError<Domain>, BoundedFaceWithoutOuterCcbError<Domain>,
//     UnboundedFaceOuterCcbError<Domain>: the arrangement is corrupt beyond the
//     misfile.
template <typename Domain, typename GpsSet>
std::size_t refile_misfiled_outer_ccbs(GpsSet& set)
{
    using Arr = typename GpsSet::Arrangement_2;
    Arr& arrangement = set.arrangement();
    const auto compare_xy = set.traits().compare_xy_2_object();

    std::size_t refiled = 0;
    std::vector<typename Arr::Halfedge_handle> holes;
    for (auto face = arrangement.faces_begin(); face != arrangement.faces_end(); ++face) {
        const std::size_t outer_count = face->number_of_outer_ccbs();
        if (face->is_unbounded()) {
            if (outer_count != 0) {
                throw UnboundedFaceOuterCcbError<Domain>("the unbounded face of a Gps arrangement carries an outer CCB.");
            }
            continue;
        }
        if (outer_count == 0) {
            throw BoundedFaceWithoutOuterCcbError<Domain>("a bounded face of a Gps arrangement has no outer CCB.");
        }
        if (outer_count == 1) continue;

        holes.clear();
        typename Arr::Halfedge_handle keep = *face->outer_ccbs_begin();
        typename Arr::Vertex_handle keep_lowest = gps_ccb_repair_detail::lowest_vertex<Arr>(keep->ccb(), compare_xy);
        for (auto outer = face->outer_ccbs_begin(); outer != face->outer_ccbs_end(); ++outer) {
            const typename Arr::Halfedge_handle candidate = *outer;
            if (candidate == keep) continue;
            const typename Arr::Vertex_handle lowest =
                gps_ccb_repair_detail::lowest_vertex<Arr>(candidate->ccb(), compare_xy);
            const CGAL::Comparison_result order = compare_xy(lowest->point(), keep_lowest->point());
            if (order == CGAL::EQUAL) {
                throw AmbiguousOuterCcbError<Domain>(
                    "two outer CCBs of one face share their lowest vertex; the arrangement is corrupt.");
            }
            if (order == CGAL::SMALLER) {
                holes.push_back(keep);
                keep = candidate;
                keep_lowest = lowest;
            } else {
                holes.push_back(candidate);
            }
        }
        for (const typename Arr::Halfedge_handle hole : holes) {
            gps_ccb_repair_detail::refile_as_inner(arrangement, face, hole);
            ++refiled;
        }
    }
    return refiled;
}
