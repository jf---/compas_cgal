#pragma once

#include "exact_region_2.h"

#include <stdexcept>
#include <string>

class ContainmentConstructionError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

enum class ContainmentKind2 {
  Segment,
  FullCircle,
  EntryDisk,
  FullCircleInDisk,
};

struct ContainmentRecord2 {
  ContainmentKind2 kind;
  bool contained;
  bool guide_anchor_in_center_domain;
  bool outer_disk_contained;
  bool disk_sweep;
  std::string strategy_version;
  std::string structural_record;

  bool matches_exact_segment(const std::string &authority_digest, double x0,
                             double y0, double x1, double y1,
                             double tool_radius) const;
  bool matches_exact_full_circle(const std::string &authority_digest,
                                 double center_x, double center_y,
                                 double phase_x, double phase_y,
                                 double tool_radius) const;
  bool matches_exact_entry_disk(const std::string &authority_digest,
                                double center_x, double center_y,
                                double entry_radius, double tool_radius) const;
  bool matches_exact_full_circle_in_disk(
      const std::string &authority_digest, double entry_center_x,
      double entry_center_y, double entry_radius, double circle_center_x,
      double circle_center_y, double phase_x, double phase_y,
      double tool_radius) const;
};

ContainmentRecord2 evaluate_exact_segment_containment(
    const ExactRegion2 &design, const ExactRegion2 &center_domain,
    const std::string &authority_digest, double x0, double y0, double x1,
    double y1, double tool_radius);

ContainmentRecord2 evaluate_exact_full_circle_containment(
    const ExactRegion2 &design, const ExactRegion2 &center_domain,
    const std::string &authority_digest, double center_x, double center_y,
    double phase_x, double phase_y, double tool_radius);

ContainmentRecord2 evaluate_exact_entry_disk_containment(
    const ExactRegion2 &design, const ExactRegion2 &center_domain,
    const std::string &authority_digest, double center_x, double center_y,
    double entry_radius, double tool_radius);

ContainmentRecord2 evaluate_exact_full_circle_in_disk(
    const std::string &authority_digest, double entry_center_x,
    double entry_center_y, double entry_radius, double circle_center_x,
    double circle_center_y, double phase_x, double phase_y, double tool_radius);

const std::string &containment_strategy_version();
