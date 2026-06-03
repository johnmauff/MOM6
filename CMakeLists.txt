cmake_minimum_required(VERSION 3.24)
project(MOM6 VERSION 6 LANGUAGES Fortran)

# project(MOM6) above already defines MOM6_SOURCE_DIR as this repo's root, used
# by all subdirectories to locate sources. We rely on it rather than
# CMAKE_SOURCE_DIR, which points at the parent project (turbo-stack) when MOM6
# is built via add_subdirectory().
list(APPEND CMAKE_MODULE_PATH "${MOM6_SOURCE_DIR}/cmake")
include(MOM6Options)

# ── Standalone dependency resolution ─────────────────────────────────────────
# When MOM6 is an add_subdirectory() of turbo-stack, the parent defines
# TURBO::infra_r8 and MARBL::marbl before this file runs — these blocks are
# skipped.  When building MOM6 directly (cmake -S MOM6 -B build), we find or
# build those dependencies ourselves.

if(NOT TARGET TURBO::infra_r8)
    # FMSConfig.cmake / TIMConfig.cmake both call find_dependency(NetCDF) on
    # the way in. Most NetCDF installs (incl. NCAR Derecho's spack module)
    # ship no NetCDFConfig.cmake — only nc-config — so a FindNetCDF.cmake
    # must be on CMAKE_MODULE_PATH *before* find_package(FMS|TIM). FMS and
    # TIM each ship their own copy; discover the bundled module via the
    # backend's install prefix (FMS_ROOT / TIM_ROOT env or var, or via
    # CMAKE_PREFIX_PATH) and prepend it so the transitive resolution works.
    if(TURBO_INFRA STREQUAL "FMS2")
        if(NOT DEFINED FMS_ROOT AND DEFINED ENV{FMS_ROOT})
            set(FMS_ROOT "$ENV{FMS_ROOT}")
        endif()
        find_path(_FMS_module_dir NAMES FindNetCDF.cmake
            HINTS ${FMS_ROOT} ${CMAKE_PREFIX_PATH}
            PATH_SUFFIXES lib/cmake/fms cmake)
        if(_FMS_module_dir)
            list(APPEND CMAKE_MODULE_PATH "${_FMS_module_dir}")
        endif()
        find_package(FMS REQUIRED COMPONENTS R8)
        add_library(turbo_infra_r8 INTERFACE)
        target_link_libraries(turbo_infra_r8 INTERFACE FMS::fms_r8)
    elseif(TURBO_INFRA STREQUAL "TIM")
        if(NOT DEFINED TIM_ROOT AND DEFINED ENV{TIM_ROOT})
            set(TIM_ROOT "$ENV{TIM_ROOT}")
        endif()
        find_path(_TIM_module_dir NAMES FindNetCDF.cmake
            HINTS ${TIM_ROOT} ${CMAKE_PREFIX_PATH}
            PATH_SUFFIXES lib/cmake/tim cmake)
        if(_TIM_module_dir)
            list(APPEND CMAKE_MODULE_PATH "${_TIM_module_dir}")
        endif()
        find_package(TIM REQUIRED COMPONENTS R8)
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

# pkg/CVMix-src and pkg/GSW-Fortran are git submodules. Fail noisily here if
# they aren't initialized — much better than the cryptic "file not found"
# cmake produces when target_sources() references missing files later.
foreach(_pkg_required
        pkg/CVMix-src/src/shared/cvmix_kinds_and_types.F90
        pkg/GSW-Fortran/modules/gsw_mod_kinds.f90)
    if(NOT EXISTS "${MOM6_SOURCE_DIR}/${_pkg_required}")
        message(FATAL_ERROR
            "Missing ${_pkg_required}. pkg/CVMix-src and pkg/GSW-Fortran are "
            "git submodules; run `git submodule update --init --recursive` "
            "in the MOM6 repo.")
    endif()
endforeach()

add_subdirectory(pkg)
add_subdirectory(src)
add_subdirectory(config_src)
