include_guard(GLOBAL)

set(NATIVE_DEPENDENCIES_MODULE_DIR "${CMAKE_CURRENT_LIST_DIR}")

set(CGAL_VERSION "6.0.1")
set(CGAL_URL "https://github.com/CGAL/cgal/releases/download/v6.0.1/CGAL-6.0.1.zip")
set(CGAL_ARCHIVE_SHA256 "884da87edd0701932302aa36d24f330e62cab27e595567507cfd2df0cf82fdd2")
set(CGAL_SOURCE_TREE_SHA256 "37353282eb92a62dd412cfdc144f9ab8706607d278ef0e3e0aff5a351043f421")
set(CGAL_LICENSE_IDENTIFIER "LicenseRef-CGAL-Package-Dependent")

set(BOOST_VERSION "1.82.0")
set(BOOST_URL "https://archives.boost.io/release/1.82.0/source/boost_1_82_0.tar.gz")
set(BOOST_ARCHIVE_SHA256 "66a469b6e608a51f8347236f4912e27dc5c60c60d7d53ae9bfe4683316c6f04c")
set(BOOST_SOURCE_TREE_SHA256 "75a4827f862f947bedf60daa19385c43289cd8666525f37cbedc79ec32a945fe")
set(BOOST_LICENSE_IDENTIFIER "BSL-1.0")

set(EIGEN_VERSION "3.4.0")
set(EIGEN_URL "https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz")
set(EIGEN_ARCHIVE_SHA256 "8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72")
set(EIGEN_SOURCE_TREE_SHA256 "f8737c321ce0287bfcb8ff8beea83486cb6e1cd32d83631cd24768b5692c42f0")
set(EIGEN_LICENSE_IDENTIFIER "MPL-2.0")

function(native_dependencies_render_lock output_variable)
    string(
        CONCAT native_lock
        "{\n"
        "  \"schema_version\": 1,\n"
        "  \"source_tree_digest\": \"native-source-tree-v1\",\n"
        "  \"packages\": [\n"
        "    {\n"
        "      \"name\": \"cgal\",\n"
        "      \"version\": \"${CGAL_VERSION}\",\n"
        "      \"url\": \"${CGAL_URL}\",\n"
        "      \"archive_sha256\": \"${CGAL_ARCHIVE_SHA256}\",\n"
        "      \"source_tree_sha256\": \"${CGAL_SOURCE_TREE_SHA256}\",\n"
        "      \"license_identifier\": \"${CGAL_LICENSE_IDENTIFIER}\"\n"
        "    },\n"
        "    {\n"
        "      \"name\": \"boost\",\n"
        "      \"version\": \"${BOOST_VERSION}\",\n"
        "      \"url\": \"${BOOST_URL}\",\n"
        "      \"archive_sha256\": \"${BOOST_ARCHIVE_SHA256}\",\n"
        "      \"source_tree_sha256\": \"${BOOST_SOURCE_TREE_SHA256}\",\n"
        "      \"license_identifier\": \"${BOOST_LICENSE_IDENTIFIER}\"\n"
        "    },\n"
        "    {\n"
        "      \"name\": \"eigen\",\n"
        "      \"version\": \"${EIGEN_VERSION}\",\n"
        "      \"url\": \"${EIGEN_URL}\",\n"
        "      \"archive_sha256\": \"${EIGEN_ARCHIVE_SHA256}\",\n"
        "      \"source_tree_sha256\": \"${EIGEN_SOURCE_TREE_SHA256}\",\n"
        "      \"license_identifier\": \"${EIGEN_LICENSE_IDENTIFIER}\"\n"
        "    }\n"
        "  ]\n"
        "}\n"
    )
    set("${output_variable}" "${native_lock}" PARENT_SCOPE)
endfunction()

function(native_dependencies_validate_lock lock_path)
    native_dependencies_render_lock(expected_lock)
    string(SHA256 expected_lock_digest "${expected_lock}")

    if(NOT EXISTS "${lock_path}")
        message(
            FATAL_ERROR
            "NATIVE_DEPENDENCY_INTEGRITY_ERROR package=manifest "
            "expected=${expected_lock_digest} observed=missing"
        )
    endif()

    file(READ "${lock_path}" observed_lock)
    if(NOT "${observed_lock}" STREQUAL "${expected_lock}")
        string(SHA256 observed_lock_digest "${observed_lock}")
        message(
            FATAL_ERROR
            "NATIVE_DEPENDENCY_INTEGRITY_ERROR package=manifest "
            "expected=${expected_lock_digest} observed=${observed_lock_digest}"
        )
    endif()
endfunction()

function(native_dependency_validate_sha256 package_name field_name digest)
    string(LENGTH "${digest}" digest_length)
    if(NOT digest_length EQUAL 64 OR NOT "${digest}" MATCHES "^[0-9a-f]+$")
        message(
            FATAL_ERROR
            "NATIVE_DEPENDENCY_INTEGRITY_ERROR package=${package_name} "
            "expected=64-char-${field_name} observed=${digest}"
        )
    endif()
endfunction()

function(native_add_locked_dependency)
    set(single_value_arguments
        NAME
        SOURCE_DIR
        VERSION
        URL
        ARCHIVE_SHA256
        SOURCE_TREE_SHA256
    )
    cmake_parse_arguments(
        NATIVE
        ""
        "${single_value_arguments}"
        ""
        ${ARGN}
    )

    foreach(required_argument IN LISTS single_value_arguments)
        if(NOT DEFINED "NATIVE_${required_argument}" OR "${NATIVE_${required_argument}}" STREQUAL "")
            message(
                FATAL_ERROR
                "NATIVE_DEPENDENCY_INTEGRITY_ERROR package=${NATIVE_NAME} "
                "expected=${required_argument} observed=missing"
            )
        endif()
    endforeach()

    native_dependency_validate_sha256(
        "${NATIVE_NAME}"
        "archive-sha256"
        "${NATIVE_ARCHIVE_SHA256}"
    )
    native_dependency_validate_sha256(
        "${NATIVE_NAME}"
        "source-tree-sha256"
        "${NATIVE_SOURCE_TREE_SHA256}"
    )

    set(observed_digest_file
        "${CMAKE_BINARY_DIR}/native-source-digests/${NATIVE_NAME}.sha256"
    )
    set(verify_command
        "${CMAKE_COMMAND}"
        "-DNATIVE_PACKAGE=${NATIVE_NAME}"
        "-DNATIVE_SOURCE_DIR=${NATIVE_SOURCE_DIR}"
        "-DNATIVE_EXPECTED_DIGEST=${NATIVE_SOURCE_TREE_SHA256}"
        "-DNATIVE_OBSERVED_DIGEST_FILE=${observed_digest_file}"
        "-P"
        "${NATIVE_DEPENDENCIES_MODULE_DIR}/VerifyNativeSource.cmake"
    )

    if(IS_SYMLINK "${NATIVE_SOURCE_DIR}")
        message(
            FATAL_ERROR
            "NATIVE_DEPENDENCY_INTEGRITY_ERROR package=${NATIVE_NAME} "
            "expected=source-directory observed=root-symlink"
        )
    endif()
    if(EXISTS "${NATIVE_SOURCE_DIR}" AND NOT IS_DIRECTORY "${NATIVE_SOURCE_DIR}")
        message(
            FATAL_ERROR
            "NATIVE_DEPENDENCY_INTEGRITY_ERROR package=${NATIVE_NAME} "
            "expected=source-directory observed=non-directory"
        )
    endif()

    set(native_source_is_populated false)
    if(IS_DIRECTORY "${NATIVE_SOURCE_DIR}")
        file(
            GLOB native_source_top_level_entries
            LIST_DIRECTORIES true
            "${NATIVE_SOURCE_DIR}/*"
        )
        if(native_source_top_level_entries)
            set(native_source_is_populated true)
        endif()
    endif()

    if(native_source_is_populated)
        execute_process(
            COMMAND ${verify_command}
            RESULT_VARIABLE verify_result
            OUTPUT_VARIABLE verify_stdout
            ERROR_VARIABLE verify_stderr
        )
        if(NOT verify_result EQUAL 0)
            message(FATAL_ERROR "${verify_stdout}${verify_stderr}")
        endif()
    else()
        ExternalProject_Add(
            "${NATIVE_NAME}_download"
            URL "${NATIVE_URL}"
            URL_HASH "SHA256=${NATIVE_ARCHIVE_SHA256}"
            SOURCE_DIR "${NATIVE_SOURCE_DIR}"
            CONFIGURE_COMMAND ""
            BUILD_COMMAND ""
            INSTALL_COMMAND ""
            UPDATE_COMMAND ""
            PATCH_COMMAND ""
            LOG_DOWNLOAD ON
            TLS_VERIFY ON
        )
        add_custom_target(
            "verify_${NATIVE_NAME}_source"
            COMMAND ${verify_command}
            VERBATIM
        )
        add_dependencies(
            "verify_${NATIVE_NAME}_source"
            "${NATIVE_NAME}_download"
        )
        add_dependencies(
            external_downloads
            "verify_${NATIVE_NAME}_source"
        )
    endif()
endfunction()
