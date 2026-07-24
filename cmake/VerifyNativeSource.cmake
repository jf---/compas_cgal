cmake_minimum_required(VERSION 3.15)

foreach(required_variable NATIVE_PACKAGE NATIVE_SOURCE_DIR)
    if(NOT DEFINED "${required_variable}" OR "${${required_variable}}" STREQUAL "")
        message(
            FATAL_ERROR
            "NATIVE_DEPENDENCY_INTEGRITY_ERROR package=unknown "
            "expected=${required_variable} observed=missing"
        )
    endif()
endforeach()

if(IS_SYMLINK "${NATIVE_SOURCE_DIR}")
    message(
        FATAL_ERROR
        "NATIVE_DEPENDENCY_INTEGRITY_ERROR package=${NATIVE_PACKAGE} "
        "expected=source-directory observed=root-symlink"
    )
endif()

if(NOT IS_DIRECTORY "${NATIVE_SOURCE_DIR}")
    message(
        FATAL_ERROR
        "NATIVE_DEPENDENCY_INTEGRITY_ERROR package=${NATIVE_PACKAGE} "
        "expected=source-tree observed=missing"
    )
endif()

file(
    GLOB_RECURSE native_source_entries
    LIST_DIRECTORIES true
    RELATIVE "${NATIVE_SOURCE_DIR}"
    "${NATIVE_SOURCE_DIR}/*"
)
list(SORT native_source_entries)

set(native_source_manifest "native-source-tree-v1\n")
foreach(relative_path IN LISTS native_source_entries)
    set(absolute_path "${NATIVE_SOURCE_DIR}/${relative_path}")

    if(IS_SYMLINK "${absolute_path}")
        message(
            FATAL_ERROR
            "NATIVE_DEPENDENCY_INTEGRITY_ERROR package=${NATIVE_PACKAGE} "
            "expected=regular-file-tree observed=symlink:${relative_path}"
        )
    endif()
    if(IS_DIRECTORY "${absolute_path}")
        continue()
    endif()
    if("${relative_path}" MATCHES "[;\r\n]")
        message(
            FATAL_ERROR
            "NATIVE_DEPENDENCY_INTEGRITY_ERROR package=${NATIVE_PACKAGE} "
            "expected=portable-relative-path observed=${relative_path}"
        )
    endif()

    string(LENGTH "${relative_path}" relative_path_length)
    file(SIZE "${absolute_path}" entry_size)
    file(SHA256 "${absolute_path}" entry_digest)
    string(
        APPEND native_source_manifest
        "F:${relative_path_length}:${relative_path}:${entry_size}:${entry_digest}\n"
    )
endforeach()

string(SHA256 native_source_digest "${native_source_manifest}")

if(NATIVE_PRINT_DIGEST)
    list(LENGTH native_source_entries native_source_entry_count)
    message(
        STATUS
        "${NATIVE_PACKAGE} source-tree SHA256: ${native_source_digest} "
        "(${native_source_entry_count} entries)"
    )
endif()

if(
    DEFINED NATIVE_EXPECTED_DIGEST
    AND NOT "${NATIVE_EXPECTED_DIGEST}" STREQUAL ""
    AND NOT "${native_source_digest}" STREQUAL "${NATIVE_EXPECTED_DIGEST}"
)
    message(
        FATAL_ERROR
        "NATIVE_DEPENDENCY_INTEGRITY_ERROR package=${NATIVE_PACKAGE} "
        "expected=${NATIVE_EXPECTED_DIGEST} observed=${native_source_digest}"
    )
endif()

if(
    DEFINED NATIVE_OBSERVED_DIGEST_FILE
    AND NOT "${NATIVE_OBSERVED_DIGEST_FILE}" STREQUAL ""
)
    get_filename_component(
        native_observed_digest_directory
        "${NATIVE_OBSERVED_DIGEST_FILE}"
        DIRECTORY
    )
    file(MAKE_DIRECTORY "${native_observed_digest_directory}")
    file(WRITE "${NATIVE_OBSERVED_DIGEST_FILE}" "${native_source_digest}\n")
endif()
