cmake_minimum_required(VERSION 3.15...3.26)

project(cgal_wrapper LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_DEBUG "")

# =====================================================================
# Set this flag to ON for developing to reduce build time.
# Set this flag to OFF for publishing for file size reduction.
# =====================================================================
option(ENABLE_PRECOMPILED_HEADERS "Enable precompiled headers for the build" OFF)


# =====================================================================
# Set maximum heap size for MSVC
# =====================================================================

if(MSVC)
set(CMAKE_GENERATOR_PLATFORM x64)
    add_compile_options(/Zm1200) 
endif()

# =====================================================================
# Build size reduction.
# =====================================================================

if (NOT ENABLE_PRECOMPILED_HEADERS)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
    if(MSVC)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /O1") # Optimize for size on MSVC
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Os") # Optimize for size on GCC/Clang
    endif()
endif()

# =====================================================================
# Dependencies
# =====================================================================
include(ExternalProject)
include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/NativeDependencies.cmake")

# Define source directories for external dependencies
set(EXTERNAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external")
set(CGAL_SOURCE_DIR "${EXTERNAL_DIR}/cgal")
set(BOOST_SOURCE_DIR "${EXTERNAL_DIR}/boost")
set(EIGEN_SOURCE_DIR "${EXTERNAL_DIR}/eigen")

# Create a custom target for all external dependencies
add_custom_target(external_downloads ALL)

native_dependencies_validate_lock(
    "${CMAKE_CURRENT_SOURCE_DIR}/docs/build/native-dependency-lock.json"
)

native_add_locked_dependency(
    NAME cgal
    SOURCE_DIR "${CGAL_SOURCE_DIR}"
    VERSION "${CGAL_VERSION}"
    URL "${CGAL_URL}"
    ARCHIVE_SHA256 "${CGAL_ARCHIVE_SHA256}"
    SOURCE_TREE_SHA256 "${CGAL_SOURCE_TREE_SHA256}"
)
native_add_locked_dependency(
    NAME boost
    SOURCE_DIR "${BOOST_SOURCE_DIR}"
    VERSION "${BOOST_VERSION}"
    URL "${BOOST_URL}"
    ARCHIVE_SHA256 "${BOOST_ARCHIVE_SHA256}"
    SOURCE_TREE_SHA256 "${BOOST_SOURCE_TREE_SHA256}"
)
native_add_locked_dependency(
    NAME eigen
    SOURCE_DIR "${EIGEN_SOURCE_DIR}"
    VERSION "${EIGEN_VERSION}"
    URL "${EIGEN_URL}"
    ARCHIVE_SHA256 "${EIGEN_ARCHIVE_SHA256}"
    SOURCE_TREE_SHA256 "${EIGEN_SOURCE_TREE_SHA256}"
)

# Add include directories
set(CGAL_INCLUDE_DIR "${CGAL_SOURCE_DIR}/include")
set(BOOST_INCLUDE_DIR "${BOOST_SOURCE_DIR}")
set(EIGEN_INCLUDE_DIR "${EIGEN_SOURCE_DIR}")

# Print include directories for verification
message(STATUS "====================================================")
message(STATUS "CGAL INCLUDE DIRECTORY: ${CGAL_INCLUDE_DIR}")
message(STATUS "BOOST INCLUDE DIRECTORY: ${BOOST_INCLUDE_DIR}")
message(STATUS "EIGEN INCLUDE DIRECTORY: ${EIGEN_INCLUDE_DIR}")
message(STATUS "====================================================")

include_directories(
    ${CGAL_INCLUDE_DIR}
    ${BOOST_INCLUDE_DIR}
    ${EIGEN_INCLUDE_DIR}
)

# Incase linking is not working check the includes
# add_compile_options(/showIncludes)

# Define libraries as header only
add_definitions(
    -DCGAL_HEADER_ONLY
    -DBOOST_ALL_NO_LIB
    -DBOOST_ALL_DYN_LINK=0
    -DCGAL_DISABLE_GMP
    -DCGAL_USE_BOOST_MP
)

# =====================================================================
# BOOST and CGAL flags to configure multiprecision
# =====================================================================
set(CGAL_DO_NOT_USE_MPZF ON CACHE BOOL "Do not use MPZF")
set(CGAL_USE_GMPXX OFF CACHE BOOL "Disable GMPXX")
set(CGAL_DISABLE_GMP ON CACHE BOOL "Disable GMP in CGAL")
set(CMAKE_DISABLE_FIND_PACKAGE_GMP ON CACHE BOOL "Disable CMake find package for GMP")
set(CMAKE_DISABLE_FIND_PACKAGE_MPFR ON CACHE BOOL "Disable CMake find package for MPFR")
set(CGAL_NT_USE_BOOST_MP ON CACHE BOOL "Use Boost Multiprecision")

# Eigen is header-only by default, no special defines needed

# =====================================================================
# Binding
# =====================================================================
if(NOT SKBUILD)
    message(WARNING "\
  This CMake file is meant to be executed using 'scikit-build'. Running
  it directly will almost certainly not produce the desired result. If
  you are a user trying to install this package, please use the command
  below, which will install all necessary build dependencies, compile
  the package in an isolated environment, and then install it.
  =====================================================================
   $ pip install .
  =====================================================================
  If you are a software developer, and this is your own package, then
  it is usually much more efficient to install the build dependencies
  in your environment once and use the following command that avoids
  a costly creation of a new virtual environment at every compilation:
  =====================================================================
   $ pip install nanobind scikit-build-core[pyproject]
   $ conda install cmake
   $ pip install --no-build-isolation -ve .
   $ pip install --no-build-isolation -ve . -Ceditable.rebuild=true
  =====================================================================
  You may optionally add -Ceditable.rebuild=true to auto-rebuild when
  the package is imported. Otherwise, you need to re-run the above
  after editing C++ files.")
endif()

# Try to import all Python components potentially needed by nanobind
find_package(Python 3.8
    REQUIRED COMPONENTS Interpreter Development.Module
    OPTIONAL_COMPONENTS Development.SABIModule)

# Import nanobind through CMake's find_package mechanism
find_package(nanobind CONFIG REQUIRED)

# We are now ready to compile the actual extension module

# Automatically include all C++ source files from the main src directory
# Note: Any new source files must follow nanobind type conversion guidelines
# for proper Python/C++ type handling, especially for:
# - Eigen::Ref for matrix passing
# - std::tuple return types
# - numpy array handling (float64, int32)
# - row-major order enforcement
# See: https://github.com/wjakob/nanobind/tree/master/tests
#
# Important: If you add a new .cpp file, you must rebuild the project with:
# pip install --no-build-isolation -ve . -Ceditable.rebuild=true

# Enhanced PCH configuration
if (ENABLE_PRECOMPILED_HEADERS)
    set(CMAKE_PCH_INSTANTIATE_TEMPLATES ON)
    set(CMAKE_PCH_WARN_INVALID ON)
endif() 

# Function to add a nanobind module with include directories
function(add_nanobind_module module_name source_file)
    nanobind_add_module(${module_name} STABLE_ABI NB_STATIC ${source_file})

    target_include_directories(${module_name} PRIVATE
        ${CGAL_INCLUDE_DIR}
        ${BOOST_INCLUDE_DIR}
        ${EIGEN_INCLUDE_DIR}
    )

    add_dependencies(${module_name} external_downloads)

    if (ENABLE_PRECOMPILED_HEADERS)
        target_precompile_headers(${module_name} PRIVATE src/compas.h)
    endif()

    install(TARGETS ${module_name} LIBRARY DESTINATION compas_cgal)
endfunction()

add_library(coverage_exact_core STATIC
    src/exact_region_2.cpp
    src/exact_sweep_2.cpp
    src/reachable_input_2.cpp
    src/reachable_arrangement_2.cpp
    src/reachable_certificate_2.cpp
    src/reachable_domain_2.cpp
    src/coverage_2.cpp
)
set_target_properties(coverage_exact_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_definitions(coverage_exact_core PRIVATE CGAL_USE_CORE=1)
target_include_directories(coverage_exact_core
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src
    PRIVATE
        ${CGAL_INCLUDE_DIR}
        ${BOOST_INCLUDE_DIR}
        ${EIGEN_INCLUDE_DIR}
)

add_executable(task3_algorithm_gate EXCLUDE_FROM_ALL
    tests/native/task3_algorithm_gate.cpp
    tests/native/task3_certificate_gate.cpp
)
target_link_libraries(task3_algorithm_gate PRIVATE coverage_exact_core)
target_compile_definitions(task3_algorithm_gate PRIVATE CGAL_USE_CORE=1)

add_executable(task9_mat_compile_gate EXCLUDE_FROM_ALL
    tests/native/task9_mat_compile_gate.cpp
    tests/native/task9_point_graph_gate.cpp
    tests/native/task9_segment_limiter_gate.cpp
    tests/native/task9_segment_segment_gate.cpp
    tests/native/task9_graph_csr_gate.cpp
    tests/native/task9_neck_gate.cpp
    tests/native/task9_neck_clearance_gate.cpp
    tests/native/task9_neck_evidence_gate.cpp
    tests/native/task9_neck_classification_gate.cpp
    tests/native/task9_mat_sampling_gate.cpp
    tests/native/task9_catalog_neck_gate.cpp
    tests/native/task9_catalog_sampling_gate.cpp
    tests/native/task9_mat_proposal_table_gate.cpp
    tests/native/task9_mat_numeric_table_gate.cpp
    tests/native/task9_domain_identity_gate.cpp
    tests/native/task9_mat_certificate_gate.cpp
    tests/native/task9_site_catalog_gate.cpp
    tests/native/task9_catalog_delaunay_gate.cpp
    tests/native/task9_catalog_voronoi_gate.cpp
    tests/native/task9_rational_sources_gate.cpp
    tests/native/task9_indexed_endpoint_binding_gate.cpp
    tests/native/task9_catalog_graph_gate.cpp
    src/segment_site_mat.cpp
    src/segment_site_catalog.cpp
    src/segment_site_delaunay.cpp
    src/segment_site_voronoi.cpp
    src/segment_site_rational_sources.cpp
    src/segment_site_catalog_graph.cpp
    src/segment_site_graph_emission.cpp
    src/segment_site_graph_csr.cpp
    src/segment_site_neck.cpp
    src/segment_site_neck_clearance.cpp
    src/segment_site_neck_evidence.cpp
    src/segment_site_neck_evidence_bytes.cpp
    src/segment_site_neck_classification.cpp
    src/segment_site_mat_sampling.cpp
    src/segment_site_catalog_neck.cpp
    src/segment_site_catalog_sampling.cpp
    src/segment_site_mat_proposal_table.cpp
    src/segment_site_mat_numeric_table.cpp
    src/segment_site_mat_certificate.cpp
    src/reachable_certificate_encoding_2.cpp
    src/segment_site_mat_spikes.cpp
    src/segment_site_parameterization.cpp
    src/segment_site_clipping.cpp
    src/segment_site_endpoint_binding.cpp
    src/segment_site_provenance.cpp
    src/canonical_encoding.cpp
    src/exact_algebraic_1.cpp
    src/exact_region_2.cpp
    src/exact_sweep_2.cpp
    src/reachable_input_2.cpp
    src/reachable_arrangement_2.cpp
    src/reachable_certificate_2.cpp
    src/reachable_domain_2.cpp
    src/continuous_tea_2/exact_polynomial_sign_2.cpp
    src/continuous_tea_2/sha256.cpp
)
add_dependencies(task9_mat_compile_gate external_downloads)
target_include_directories(task9_mat_compile_gate PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CGAL_INCLUDE_DIR}
    ${BOOST_INCLUDE_DIR}
)
target_compile_definitions(
    task9_mat_compile_gate
    PRIVATE
        CGAL_USE_CORE=1
        CGAL_CORE_USE_BOOST_BACKEND=1
)

# Add new modules here
add_nanobind_module(_types_std src/types_std.cpp)
add_nanobind_module(_booleans src/booleans.cpp)
add_nanobind_module(_meshing src/meshing.cpp)
add_nanobind_module(_intersections src/intersections.cpp)
add_nanobind_module(_measure src/measure.cpp)
add_nanobind_module(_reconstruction src/reconstruction.cpp)
add_nanobind_module(_skeletonization src/skeletonization.cpp)
add_nanobind_module(_slicer src/slicer.cpp)
add_nanobind_module(_straight_skeleton_2 src/straight_skeleton_2.cpp)
add_nanobind_module(_triangulation src/triangulation.cpp)
add_nanobind_module(_subdivision src/subdivision.cpp)
add_nanobind_module(_polylines src/polylines.cpp)
add_nanobind_module(_toolpath src/toolpath.cpp)
add_nanobind_module(_stock_2 "src/stock_2.cpp;src/exact_depletion_2.cpp;src/engagement_2.cpp")
add_nanobind_module(
    _continuous_tea_2
    "src/continuous_tea_2.cpp;src/continuous_tea_2/boundary_events.cpp;src/continuous_tea_2/parameter_charts.cpp;src/exact_algebraic_1.cpp;src/continuous_tea_2/event_partition.cpp;src/continuous_tea_2/trimmed_line_branches.cpp;src/continuous_tea_2/cap_partition.cpp;src/continuous_tea_2/event_certificate.cpp;src/continuous_tea_2/sha256.cpp;src/continuous_tea_2/segment_source.cpp;src/continuous_tea_2/exact_polynomial_sign_2.cpp;src/continuous_tea_2/segment_projection.cpp;src/continuous_tea_2/segment_pair_projection.cpp;src/continuous_tea_2/segment_strata.cpp;src/continuous_tea_2/segment_fibre.cpp;src/continuous_tea_2/segment_partition.cpp"
)
target_compile_definitions(
    _continuous_tea_2
    PRIVATE
        CGAL_USE_CORE=1
        CGAL_CORE_USE_BOOST_BACKEND=1
)
add_nanobind_module(_coverage_2 src/coverage_bindings_2.cpp)
target_link_libraries(_coverage_2 PRIVATE coverage_exact_core)
target_compile_definitions(_coverage_2 PRIVATE CGAL_USE_CORE=1)

add_nanobind_module(_geodesics src/geodesics.cpp)
add_nanobind_module(_isolines src/isolines.cpp)
