#include "containment_2.h"

#include "exact_sweep_2.h"

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include <CGAL/enum.h>
#include <CGAL/number_utils.h>

namespace {

const std::string STRATEGY_VERSION = "exact-sweep-containment-v2";
constexpr std::size_t SHA256_DIGEST_BYTES = 32;

double finite_value(const double value, const std::string &name) {
  if (!std::isfinite(value)) {
    throw ContainmentConstructionError(name + " must be finite.");
  }
  return value == 0.0 ? 0.0 : value;
}

ReachFT positive_length(const double value, const std::string &name) {
  const double canonical = finite_value(value, name);
  if (!(canonical > 0.0)) {
    throw ContainmentConstructionError(name + " must be positive.");
  }
  return ReachFT(canonical);
}

ReachKernelPoint exact_point(const double x, const double y,
                             const std::string &name) {
  return ReachKernelPoint(finite_value(x, name + " x"),
                          finite_value(y, name + " y"));
}

ReachKernelVector exact_vector(const double x, const double y,
                               const std::string &name) {
  const ReachKernelVector result(finite_value(x, name + " x"),
                                 finite_value(y, name + " y"));
  if (result == ReachKernelVector(CGAL::NULL_VECTOR)) {
    throw ContainmentConstructionError(name + " must be nonzero.");
  }
  return result;
}

const std::string &
validated_authority_digest(const std::string &authority_digest) {
  if (authority_digest.size() != SHA256_DIGEST_BYTES) {
    throw ContainmentConstructionError(
        "containment authority must be one SHA-256 digest.");
  }
  return authority_digest;
}

void require_domain_roles(const ExactRegion2 &design,
                          const ExactRegion2 &center_domain) {
  if (design.role() != ExactRegionRole2::Design) {
    throw ContainmentConstructionError(
        "containment design must have the exact design-region role.");
  }
  if (center_domain.role() != ExactRegionRole2::CenterDomain) {
    throw ContainmentConstructionError(
        "containment center domain must have the exact center-domain role.");
  }
}

bool contains_exact_point(const ExactRegion2 &region,
                          const ReachKernelPoint &point) {
  return region.set().oriented_side(ReachPoint(point.x(), point.y())) !=
         CGAL::ON_NEGATIVE_SIDE;
}

ReachSet disk_set(const ReachKernelPoint &center, const ReachFT &radius) {
  ReachSet result;
  result.insert(reach_disk_polygon(center, radius));
  return result;
}

ReachSet capsule_set(const ReachKernelPoint &start, const ReachKernelPoint &end,
                     const ReachFT &radius) {
  return reach_join_parts(reach_capsule_parts(start, end, radius), {});
}

std::string bool_record(const bool value) {
  return reach_u64_record(value ? 1U : 0U);
}

std::string segment_record(const std::string &authority_digest, const double x0,
                           const double y0, const double x1, const double y1,
                           const double tool_radius, const bool contained,
                           const bool anchors) {
  return reach_tagged_record("exact-segment-containment-v1",
                             {
                                 STRATEGY_VERSION,
                                 authority_digest,
                                 reach_binary64_record(x0),
                                 reach_binary64_record(y0),
                                 reach_binary64_record(x1),
                                 reach_binary64_record(y1),
                                 reach_binary64_record(tool_radius),
                                 bool_record(anchors),
                                 bool_record(contained),
                             });
}

std::string circle_record(const std::string &authority_digest,
                          const double center_x, const double center_y,
                          const double phase_x, const double phase_y,
                          const double tool_radius, const bool contained,
                          const bool anchor, const bool outer_disk_contained,
                          const bool disk_sweep) {
  return reach_tagged_record("exact-full-circle-containment-v1",
                             {
                                 STRATEGY_VERSION,
                                 authority_digest,
                                 reach_binary64_record(center_x),
                                 reach_binary64_record(center_y),
                                 reach_binary64_record(phase_x),
                                 reach_binary64_record(phase_y),
                                 reach_binary64_record(tool_radius),
                                 bool_record(anchor),
                                 bool_record(outer_disk_contained),
                                 bool_record(disk_sweep),
                                 bool_record(contained),
                             });
}

std::string entry_record(const std::string &authority_digest,
                         const double center_x, const double center_y,
                         const double entry_radius, const double tool_radius,
                         const bool contained, const bool center_in_domain) {
  return reach_tagged_record("exact-entry-disk-containment-v1",
                             {
                                 STRATEGY_VERSION,
                                 authority_digest,
                                 reach_binary64_record(center_x),
                                 reach_binary64_record(center_y),
                                 reach_binary64_record(entry_radius),
                                 reach_binary64_record(tool_radius),
                                 bool_record(center_in_domain),
                                 bool_record(contained),
                             });
}

std::string
circle_in_disk_record(const std::string &authority_digest,
                      const double entry_center_x, const double entry_center_y,
                      const double entry_radius, const double circle_center_x,
                      const double circle_center_y, const double phase_x,
                      const double phase_y, const double tool_radius,
                      const bool contained, const bool disk_sweep) {
  return reach_tagged_record("exact-full-circle-in-entry-disk-v1",
                             {
                                 STRATEGY_VERSION,
                                 authority_digest,
                                 reach_binary64_record(entry_center_x),
                                 reach_binary64_record(entry_center_y),
                                 reach_binary64_record(entry_radius),
                                 reach_binary64_record(circle_center_x),
                                 reach_binary64_record(circle_center_y),
                                 reach_binary64_record(phase_x),
                                 reach_binary64_record(phase_y),
                                 reach_binary64_record(tool_radius),
                                 bool_record(disk_sweep),
                                 bool_record(contained),
                             });
}

} // namespace

const std::string &containment_strategy_version() { return STRATEGY_VERSION; }

bool ContainmentRecord2::matches_exact_segment(
    const std::string &authority_digest, const double x0, const double y0,
    const double x1, const double y1, const double tool_radius) const {
  return kind == ContainmentKind2::Segment &&
         strategy_version == STRATEGY_VERSION &&
         structural_record == segment_record(authority_digest, x0, y0, x1, y1,
                                             tool_radius, contained,
                                             guide_anchor_in_center_domain);
}

bool ContainmentRecord2::matches_exact_full_circle(
    const std::string &authority_digest, const double center_x,
    const double center_y, const double phase_x, const double phase_y,
    const double tool_radius) const {
  return kind == ContainmentKind2::FullCircle &&
         strategy_version == STRATEGY_VERSION &&
         structural_record == circle_record(authority_digest, center_x,
                                            center_y, phase_x, phase_y,
                                            tool_radius, contained,
                                            guide_anchor_in_center_domain,
                                            outer_disk_contained, disk_sweep);
}

bool ContainmentRecord2::matches_exact_entry_disk(
    const std::string &authority_digest, const double center_x,
    const double center_y, const double entry_radius,
    const double tool_radius) const {
  return kind == ContainmentKind2::EntryDisk &&
         strategy_version == STRATEGY_VERSION &&
         structural_record == entry_record(authority_digest, center_x, center_y,
                                           entry_radius, tool_radius, contained,
                                           guide_anchor_in_center_domain);
}

bool ContainmentRecord2::matches_exact_full_circle_in_disk(
    const std::string &authority_digest, const double entry_center_x,
    const double entry_center_y, const double entry_radius,
    const double circle_center_x, const double circle_center_y,
    const double phase_x, const double phase_y,
    const double tool_radius) const {
  return kind == ContainmentKind2::FullCircleInDisk &&
         strategy_version == STRATEGY_VERSION &&
         structural_record ==
             circle_in_disk_record(authority_digest, entry_center_x,
                                   entry_center_y, entry_radius,
                                   circle_center_x, circle_center_y, phase_x,
                                   phase_y, tool_radius, contained, disk_sweep);
}

ContainmentRecord2 evaluate_exact_segment_containment(
    const ExactRegion2 &design, const ExactRegion2 &center_domain,
    const std::string &authority_digest, const double x0, const double y0,
    const double x1, const double y1, const double tool_radius) {
  require_domain_roles(design, center_domain);
  validated_authority_digest(authority_digest);
  const ReachKernelPoint start = exact_point(x0, y0, "segment start");
  const ReachKernelPoint end = exact_point(x1, y1, "segment end");
  if (start == end) {
    throw ContainmentConstructionError(
        "exact segment containment requires nonzero progress.");
  }
  const ReachFT radius = positive_length(tool_radius, "segment tool radius");
  const bool anchors = contains_exact_point(center_domain, start) &&
                       contains_exact_point(center_domain, end);
  // Regularized area erosion omits valid lower-dimensional centre loci, such
  // as the axis of a corridor exactly one tool diameter wide. The exact sweep
  // subset is the authoritative gouge predicate; anchor membership remains
  // independent diagnostic evidence for full-dimensional centre domains.
  const bool contained =
      reach_exact_subset(capsule_set(start, end, radius), design.set());
  return {
      ContainmentKind2::Segment,
      contained,
      anchors,
      false,
      false,
      STRATEGY_VERSION,
      segment_record(authority_digest, x0, y0, x1, y1, tool_radius, contained,
                     anchors),
  };
}

ContainmentRecord2 evaluate_exact_full_circle_containment(
    const ExactRegion2 &design, const ExactRegion2 &center_domain,
    const std::string &authority_digest, const double center_x,
    const double center_y, const double phase_x, const double phase_y,
    const double tool_radius) {
  require_domain_roles(design, center_domain);
  validated_authority_digest(authority_digest);
  const ReachKernelPoint center =
      exact_point(center_x, center_y, "full-circle center");
  const ReachKernelVector phase =
      exact_vector(phase_x, phase_y, "full-circle phase vector");
  const ReachFT radius =
      positive_length(tool_radius, "full-circle tool radius");
  const ReachKernelPoint anchor = center + phase;
  const bool anchor_in_domain = contains_exact_point(center_domain, anchor);
  const ReachFT guide_radius = CGAL::sqrt(phase.squared_length());
  const bool disk_sweep = CGAL::compare(guide_radius, radius) != CGAL::LARGER;
  const ReachSet outer_disk = disk_set(center, guide_radius + radius);
  const bool outer_disk_contained =
      reach_exact_subset(outer_disk, design.set());
  const ReachSet sweep = reach_full_circle_sweep(center, phase, radius);
  const bool contained =
      anchor_in_domain && reach_exact_subset(sweep, design.set());
  return {
      ContainmentKind2::FullCircle,
      contained,
      anchor_in_domain,
      outer_disk_contained,
      disk_sweep,
      STRATEGY_VERSION,
      circle_record(authority_digest, center_x, center_y, phase_x, phase_y,
                    tool_radius, contained, anchor_in_domain,
                    outer_disk_contained, disk_sweep),
  };
}

ContainmentRecord2 evaluate_exact_entry_disk_containment(
    const ExactRegion2 &design, const ExactRegion2 &center_domain,
    const std::string &authority_digest, const double center_x,
    const double center_y, const double entry_radius,
    const double tool_radius) {
  require_domain_roles(design, center_domain);
  validated_authority_digest(authority_digest);
  const ReachKernelPoint center =
      exact_point(center_x, center_y, "entry center");
  const ReachFT exact_entry_radius =
      positive_length(entry_radius, "entry radius");
  static_cast<void>(positive_length(tool_radius, "entry tool radius"));
  const bool center_in_domain = contains_exact_point(center_domain, center);
  const bool contained =
      center_in_domain &&
      reach_exact_subset(disk_set(center, exact_entry_radius), design.set());
  return {
      ContainmentKind2::EntryDisk,
      contained,
      center_in_domain,
      contained,
      true,
      STRATEGY_VERSION,
      entry_record(authority_digest, center_x, center_y, entry_radius,
                   tool_radius, contained, center_in_domain),
  };
}

ContainmentRecord2 evaluate_exact_full_circle_in_disk(
    const std::string &authority_digest, const double entry_center_x,
    const double entry_center_y, const double entry_radius,
    const double circle_center_x, const double circle_center_y,
    const double phase_x, const double phase_y, const double tool_radius) {
  validated_authority_digest(authority_digest);
  const ReachKernelPoint entry_center =
      exact_point(entry_center_x, entry_center_y, "entry disk center");
  const ReachFT exact_entry_radius =
      positive_length(entry_radius, "entry disk radius");
  const ReachKernelPoint circle_center =
      exact_point(circle_center_x, circle_center_y, "full-circle center");
  const ReachKernelVector phase =
      exact_vector(phase_x, phase_y, "full-circle phase vector");
  const ReachFT exact_tool_radius =
      positive_length(tool_radius, "full-circle tool radius");
  const ReachFT guide_radius = CGAL::sqrt(phase.squared_length());
  const bool disk_sweep =
      CGAL::compare(guide_radius, exact_tool_radius) != CGAL::LARGER;
  const bool contained = reach_exact_subset(
      reach_full_circle_sweep(circle_center, phase, exact_tool_radius),
      disk_set(entry_center, exact_entry_radius));
  return {
      ContainmentKind2::FullCircleInDisk,
      contained,
      true,
      contained,
      disk_sweep,
      STRATEGY_VERSION,
      circle_in_disk_record(authority_digest, entry_center_x, entry_center_y,
                            entry_radius, circle_center_x, circle_center_y,
                            phase_x, phase_y, tool_radius, contained,
                            disk_sweep),
  };
}
