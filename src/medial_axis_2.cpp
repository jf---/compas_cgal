#include "reachable_errors_2.h"
#include "segment_site_catalog_graph.h"
#include "segment_site_mat_bundle.h"
#include "segment_site_mat_certificate.h"
#include "segment_site_mat_numeric_table.h"
#include "segment_site_mat_sampling.h"
#include "segment_site_neck_classification.h"
#include "segment_site_neck_evidence_bytes.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <Eigen/Core>
#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace {

class MedialAxisConstructionError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

class BindingInvalidReachableDomainInputError
    : public MedialAxisConstructionError {
public:
  using MedialAxisConstructionError::MedialAxisConstructionError;
};

class BindingReachableArrangementTopologyError
    : public MedialAxisConstructionError {
public:
  using MedialAxisConstructionError::MedialAxisConstructionError;
};

class BindingPocketNotMachinableError : public MedialAxisConstructionError {
public:
  using MedialAxisConstructionError::MedialAxisConstructionError;
};

class BindingReachableMaterialContainmentError
    : public MedialAxisConstructionError {
public:
  using MedialAxisConstructionError::MedialAxisConstructionError;
};

class BindingInvalidMatSamplingPolicyError
    : public MedialAxisConstructionError {
public:
  using MedialAxisConstructionError::MedialAxisConstructionError;
};

class BindingConicSamplingLimitError : public MedialAxisConstructionError {
public:
  using MedialAxisConstructionError::MedialAxisConstructionError;
};

class BindingUnsupportedCanonicalMatLShapeGraphError
    : public MedialAxisConstructionError {
public:
  using MedialAxisConstructionError::MedialAxisConstructionError;
};

class BindingInvalidMatCertificateReplayError
    : public MedialAxisConstructionError {
public:
  using MedialAxisConstructionError::MedialAxisConstructionError;
};

class BindingInvalidMatNeckEvidenceError : public MedialAxisConstructionError {
public:
  using MedialAxisConstructionError::MedialAxisConstructionError;
};

class BindingInvalidMatNeckWidthBoundariesError
    : public MedialAxisConstructionError {
public:
  using MedialAxisConstructionError::MedialAxisConstructionError;
};

template <typename Scalar, std::size_t ColumnCount>
using RowMatrix = Eigen::Matrix<Scalar, Eigen::Dynamic,
                                static_cast<int>(ColumnCount), Eigen::RowMajor>;

template <typename Scalar, std::size_t ColumnCount>
RowMatrix<Scalar, ColumnCount>
row_matrix(const std::vector<std::array<Scalar, ColumnCount>> &records) {
  RowMatrix<Scalar, ColumnCount> result(
      static_cast<Eigen::Index>(records.size()),
      static_cast<Eigen::Index>(ColumnCount));
  for (std::size_t row = 0; row < records.size(); ++row) {
    for (std::size_t column = 0; column < ColumnCount; ++column) {
      result(static_cast<Eigen::Index>(row),
             static_cast<Eigen::Index>(column)) = records[row][column];
    }
  }
  return result;
}

template <typename Scalar>
Eigen::Matrix<Scalar, Eigen::Dynamic, 1>
column_vector(const std::vector<Scalar> &records) {
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> result(
      static_cast<Eigen::Index>(records.size()));
  for (std::size_t index = 0; index < records.size(); ++index) {
    result(static_cast<Eigen::Index>(index)) = records[index];
  }
  return result;
}

nb::bytes bytes_value(const std::string &value) {
  return nb::bytes(value.data(), value.size());
}

nb::tuple bytes_tuple(const std::vector<std::string> &records) {
  nb::list result;
  for (const std::string &record : records) {
    result.append(bytes_value(record));
  }
  return nb::tuple(result);
}

nb::tuple catalog_site_ids(const CanonicalMatSiteCatalog2 &catalog) {
  nb::list result;
  for (const CanonicalMatSite2 &site : catalog.sites()) {
    result.append(bytes_value(site.stable_id()));
  }
  return nb::tuple(result);
}

nb::tuple
nested_bytes_tuple(const std::vector<std::vector<std::string>> &records) {
  nb::list result;
  for (const std::vector<std::string> &record : records) {
    result.append(bytes_tuple(record));
  }
  return nb::tuple(result);
}

std::string from_bytes(const nb::bytes &value) {
  char *data = nullptr;
  Py_ssize_t size = 0;
  if (PyBytes_AsStringAndSize(value.ptr(), &data, &size) != 0) {
    throw nb::python_error();
  }
  return std::string(data, static_cast<std::size_t>(size));
}

std::vector<std::string> from_bytes_tuple(const nb::tuple &values,
                                          const char *name) {
  std::vector<std::string> result;
  result.reserve(values.size());
  for (nb::handle value : values) {
    if (!PyBytes_Check(value.ptr())) {
      throw nb::type_error(
          (std::string(name) + " must contain only exact bytes").c_str());
    }
    result.push_back(from_bytes(nb::borrow<nb::bytes>(value)));
  }
  return result;
}

auto validate_and_classify_necks(const SegmentSiteMatBundle2 &owner,
                                 const nb::bytes &mat_certificate,
                                 const nb::tuple &neck_evidence,
                                 const nb::tuple &squared_width_boundaries) {
  try {
    auto [classes, certificates] = owner.validate_and_classify_necks(
        from_bytes(mat_certificate),
        from_bytes_tuple(neck_evidence, "neck evidence"),
        from_bytes_tuple(squared_width_boundaries, "squared-width boundaries"));
    return std::make_tuple(column_vector(classes), bytes_tuple(certificates));
  } catch (const InvalidMatCertificateReplayError &error) {
    throw BindingInvalidMatCertificateReplayError(error.what());
  } catch (const InvalidMatNeckEvidenceBytesError &error) {
    throw BindingInvalidMatNeckEvidenceError(error.what());
  } catch (const InvalidMatNeckWidthBoundariesError &error) {
    throw BindingInvalidMatNeckWidthBoundariesError(error.what());
  } catch (const std::exception &error) {
    throw MedialAxisConstructionError(error.what());
  }
}

SegmentSiteMatBundle2
build_bundle(Eigen::Ref<const compas::RowMatrixXd> vertices,
             const std::vector<compas::RowMatrixXd> &holes,
             const double tool_radius, const double station_spacing,
             const double max_sagitta, const std::size_t max_refinement_depth) {
  try {
    const CanonicalReachInput2 input =
        canonical_reach_input(vertices, holes, tool_radius);
    return SegmentSiteMatBundle2::build(
        input, MatStationSpacingMm2::build(station_spacing),
        MatSagittaBoundMm2::build(max_sagitta), max_refinement_depth);
  } catch (const InvalidReachableDomainInputError &error) {
    throw BindingInvalidReachableDomainInputError(error.what());
  } catch (const ReachableArrangementTopologyError &error) {
    throw BindingReachableArrangementTopologyError(error.what());
  } catch (const PocketNotMachinableError &error) {
    throw BindingPocketNotMachinableError(error.what());
  } catch (const ReachableMaterialContainmentError &error) {
    throw BindingReachableMaterialContainmentError(error.what());
  } catch (const InvalidMatSamplingPolicyError &error) {
    throw BindingInvalidMatSamplingPolicyError(error.what());
  } catch (const ConicSamplingLimitError &error) {
    throw BindingConicSamplingLimitError(error.what());
  } catch (const UnsupportedCanonicalMatLShapeGraphError &error) {
    throw BindingUnsupportedCanonicalMatLShapeGraphError(error.what());
  } catch (const std::exception &error) {
    throw MedialAxisConstructionError(error.what());
  }
}

auto projection_tuple(const MatNumericMatTable2 &table) {
  const MatNumericGraphTable2 &graph = table.proposal.graph;
  const MatNumericSampleTable2 &samples = table.proposal.samples;
  return std::make_tuple(
      row_matrix(table.nodes), row_matrix(graph.edges),
      column_vector(graph.node_site_offsets),
      column_vector(graph.node_site_ids), row_matrix(graph.site_provenance),
      row_matrix(graph.edge_endpoint_provenance_flags),
      column_vector(graph.endpoint_feature_offsets),
      row_matrix(graph.endpoint_features), row_matrix(graph.edge_exact_flags),
      row_matrix(samples.sample_centers),
      column_vector(samples.sample_clearance),
      column_vector(samples.sample_guide_radius),
      row_matrix(samples.sample_flags),
      column_vector(samples.edge_sample_offsets),
      column_vector(samples.sample_parameter), bytes_tuple(table.neck_evidence),
      column_vector(table.neck_cut_offsets),
      column_vector(table.neck_cut_edge_ids),
      bytes_value(table.center_domain_digest),
      bytes_value(table.mat_certificate));
}

auto segment_site_medial_axis(Eigen::Ref<const compas::RowMatrixXd> vertices,
                              const std::vector<compas::RowMatrixXd> &holes,
                              const double tool_radius,
                              const double station_spacing,
                              const double max_sagitta,
                              const std::size_t max_refinement_depth) {
  const SegmentSiteMatBundle2 owner =
      build_bundle(vertices, holes, tool_radius, station_spacing, max_sagitta,
                   max_refinement_depth);
  return projection_tuple(owner.numeric_table());
}

} // namespace

NB_MODULE(_medial_axis_2, m) {
  nb::exception<MedialAxisConstructionError> construction_error(
      m, "MedialAxisConstructionError");
  nb::exception<BindingInvalidReachableDomainInputError>(
      m, "InvalidReachableDomainInputError", construction_error.ptr());
  nb::exception<BindingReachableArrangementTopologyError>(
      m, "ReachableArrangementTopologyError", construction_error.ptr());
  nb::exception<BindingPocketNotMachinableError>(m, "PocketNotMachinableError",
                                                 construction_error.ptr());
  nb::exception<BindingReachableMaterialContainmentError>(
      m, "ReachableMaterialContainmentError", construction_error.ptr());
  nb::exception<BindingInvalidMatSamplingPolicyError>(
      m, "InvalidMatSamplingPolicyError", construction_error.ptr());
  nb::exception<BindingConicSamplingLimitError>(m, "ConicSamplingLimitError",
                                                construction_error.ptr());
  nb::exception<BindingUnsupportedCanonicalMatLShapeGraphError>(
      m, "UnsupportedCanonicalMatLShapeGraphError", construction_error.ptr());
  nb::exception<BindingInvalidMatCertificateReplayError>(
      m, "InvalidMatCertificateReplayError", construction_error.ptr());
  nb::exception<BindingInvalidMatNeckEvidenceError>(
      m, "InvalidMatNeckEvidenceError", construction_error.ptr());
  nb::exception<BindingInvalidMatNeckWidthBoundariesError>(
      m, "InvalidMatNeckWidthBoundariesError", construction_error.ptr());

  nb::class_<SegmentSiteMatBundle2>(m, "SegmentSiteMedialAxis")
      .def_static("build", &build_bundle, "vertices"_a, "holes"_a,
                  "tool_radius"_a, "station_spacing"_a, "max_sagitta"_a,
                  "max_refinement_depth"_a)
      .def_prop_ro("projection",
                   [](const SegmentSiteMatBundle2 &owner) {
                     return projection_tuple(owner.numeric_table());
                   })
      .def_prop_ro("node_ids",
                   [](const SegmentSiteMatBundle2 &owner) {
                     return bytes_tuple(
                         owner.numeric_table().proposal.graph.node_ids);
                   })
      .def_prop_ro("edge_ids",
                   [](const SegmentSiteMatBundle2 &owner) {
                     return bytes_tuple(owner.edge_ids());
                   })
      .def_prop_ro("site_ids",
                   [](const SegmentSiteMatBundle2 &owner) {
                     return catalog_site_ids(owner.catalog());
                   })
      .def_prop_ro(
          "original_dual_ids",
          [](const SegmentSiteMatBundle2 &owner) {
            return bytes_tuple(
                owner.numeric_table().proposal.graph.original_dual_ids);
          })
      .def_prop_ro("sample_parameter_ids",
                   [](const SegmentSiteMatBundle2 &owner) {
                     return bytes_tuple(owner.sample_parameter_ids());
                   })
      .def_prop_ro("neck_owner_ids",
                   [](const SegmentSiteMatBundle2 &owner) {
                     return bytes_tuple(owner.neck_owner_ids());
                   })
      .def_prop_ro("neck_defining_site_ids",
                   [](const SegmentSiteMatBundle2 &owner) {
                     return nested_bytes_tuple(owner.neck_defining_site_ids());
                   })
      .def("validate_and_classify_necks", &validate_and_classify_necks,
           "mat_certificate"_a, "neck_evidence"_a,
           "squared_width_boundaries"_a);

  m.def("segment_site_medial_axis", &segment_site_medial_axis, "vertices"_a,
        "holes"_a, "tool_radius"_a, "station_spacing"_a, "max_sagitta"_a,
        "max_refinement_depth"_a);
}
