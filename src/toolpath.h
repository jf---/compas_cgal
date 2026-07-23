#pragma once

#include "compas.h"

// CGAL straight skeleton / offset
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/create_straight_skeleton_2.h>
#include <CGAL/create_straight_skeleton_from_polygon_with_holes_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/intersections.h>

/**
 * @brief Interior straight-skeleton vertices with exact boundary clearance.
 *
 * The straight skeleton coincides with the medial axis for convex polygons
 * only; for non-convex polygons the loci differ near reflex vertices, but the
 * returned radii are exact clearances (largest gouge-free disk) at each locus.
 *
 * @param vertices Outer boundary as Nx3 row-major matrix (float64)
 * @param holes Hole polygons, each Mx3 (float64); may be empty
 * @return std::tuple<compas::RowMatrixXd, compas::RowMatrixXd> containing:
 *         - skeleton points as Mx3 matrix (float64)
 *         - exact clearance radii as Mx1 matrix (float64)
 */
std::tuple<compas::RowMatrixXd, compas::RowMatrixXd>
pmp_polygon_skeleton_clearance(
    Eigen::Ref<const compas::RowMatrixXd> vertices,
    const std::vector<compas::RowMatrixXd>& holes
);

/**
 * @brief Certified gouge-free trochoidal pocket toolpaths.
 *
 * Trochoid stations are sampled along straight-skeleton edges at
 * engagement-bounded spacing; every station radius derives from an EXACT
 * clearance query at its own center (no interpolation), and every bridge
 * tangent is certified against the boundary at tool radius.
 *
 * @param vertices Outer boundary as Nx3 row-major matrix (float64)
 * @param holes Hole polygons, each Mx3 (float64); may be empty
 * @param tool_diameter Cutter diameter
 * @param stepover Maximum radial engagement per trochoid cycle
 * @param pitch Desired trochoid advance per cycle (capped by stepover)
 * @param min_trochoid_radius Minimum trochoid radius
 * @param max_trochoid_radius Maximum trochoid radius
 * @param mat_scale Scale factor for clearance-derived radius
 * @param radial_clearance Safety clearance subtracted from available radius
 * @param samples_per_cycle Polyline samples per arc primitive (>= 4)
 * @param max_passes Maximum number of emitted toolpaths
 * @param climb true: climb milling (CW arcs); false: conventional (CCW)
 * @return (list of Kx3 polylines, number of skeleton edges skipped by max_passes)
 */
std::tuple<std::vector<compas::RowMatrixXd>, int>
pmp_trochoidal_mat_toolpath(
    Eigen::Ref<const compas::RowMatrixXd> vertices,
    const std::vector<compas::RowMatrixXd>& holes,
    double tool_diameter,
    double stepover,
    double pitch,
    double min_trochoid_radius,
    double max_trochoid_radius,
    double mat_scale,
    double radial_clearance,
    int samples_per_cycle,
    int max_passes,
    bool climb
);

/**
 * @brief Create linked circular arc toolpath with path ordering, leads, links, retract/plunge.
 *
 * Meta matrix columns: [path_index, type(0=line,1=arc,2=circle), clockwise, operation]
 * Operation codes: 0=cut, 1=lead_in, 2=lead_out, 3=link, 4=retract, 5=plunge
 *
 * Leads are certified against the boundary and shrink (or drop) to stay
 * gouge-free.  Flat links (no clearance_z) that would gouge raise instead of
 * being emitted.
 *
 * @return std::tuple of (meta Nx4, starts Nx3, ends Nx3, centers Nx3, radii Nx1,
 *         polyline Mx3, start_tangents Nx3, end_tangents Nx3,
 *         skipped edge count)
 */
std::tuple<compas::RowMatrixXd, compas::RowMatrixXd, compas::RowMatrixXd, compas::RowMatrixXd, compas::RowMatrixXd, compas::RowMatrixXd, compas::RowMatrixXd, compas::RowMatrixXd, int>
pmp_trochoidal_mat_toolpath_circular(
    Eigen::Ref<const compas::RowMatrixXd> vertices,
    const std::vector<compas::RowMatrixXd>& holes,
    double tool_diameter,
    double stepover,
    double pitch,
    double min_trochoid_radius,
    double max_trochoid_radius,
    double mat_scale,
    double radial_clearance,
    int samples_per_cycle,
    int max_passes,
    bool climb,
    double lead_in,
    double lead_out,
    bool link_paths,
    bool optimize_order,
    double cut_z,
    double clearance_z,
    bool has_clearance_z,
    bool retract_at_end,
    double samples_per_radian
);
