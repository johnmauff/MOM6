cmake_minimum_required(VERSION 3.24)
project(MOM6 VERSION 6 LANGUAGES Fortran)

# Absolute path to this repo's root; used by all subdirectories to locate sources.
set(MOM6_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(MOM6Options)

# ── Standalone dependency resolution ─────────────────────────────────────────
# When MOM6 is an add_subdirectory() of turbo-stack, the parent defines
# TURBO::infra_r8 and MARBL::marbl before this file runs — these blocks are
# skipped.  When building MOM6 directly (cmake -S MOM6 -B build), we find or
# build those dependencies ourselves.

if(NOT TARGET TURBO::infra_r8)
    if(TURBO_INFRA STREQUAL "FMS2")
        find_package(FMS REQUIRED COMPONENTS R8)
        list(APPEND CMAKE_MODULE_PATH "${FMS_INSTALL_PREFIX}/lib/cmake/fms")
        add_library(turbo_infra_r8 INTERFACE)
        target_link_libraries(turbo_infra_r8 INTERFACE FMS::fms_r8)
    elseif(TURBO_INFRA STREQUAL "TIM")
        find_package(TIM REQUIRED COMPONENTS R8)
        list(APPEND CMAKE_MODULE_PATH "${TIM_DIR}")
        add_library(turbo_infra_r8 INTERFACE)
        target_link_libraries(turbo_infra_r8 INTERFACE TIM::tim_r8)
    else()
        message(FATAL_ERROR "Unknown TURBO_INFRA='${TURBO_INFRA}'. Valid values: FMS2, TIM")
    endif()
    add_library(TURBO::infra_r8 ALIAS turbo_infra_r8)
    find_package(MPI REQUIRED COMPONENTS Fortran)
    find_package(NetCDF REQUIRED COMPONENTS C Fortran)
endif()

if(NOT TARGET MARBL::marbl)
    if(NOT MARBL_SOURCE_DIR)
        set(MARBL_SOURCE_DIR "$ENV{MARBL_ROOT}" CACHE PATH
            "Path to MARBL source checkout (required for standalone MOM6 builds)")
    endif()
    if(NOT MARBL_SOURCE_DIR OR NOT IS_DIRECTORY "${MARBL_SOURCE_DIR}")
        message(FATAL_ERROR
            "MARBL::marbl not defined and MARBL_SOURCE_DIR is not set or invalid.\n"
            "Set MARBL_ROOT in the environment or pass -DMARBL_SOURCE_DIR=/path/to/MARBL")
    endif()
    include(MARBLTargets)
endif()

add_subdirectory(pkg)
add_subdirectory(src)
add_subdirectory(config_src)
