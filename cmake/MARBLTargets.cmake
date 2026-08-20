# MARBLTargets.cmake
# Creates the MARBL::marbl target from source.
# Requires MARBL_SOURCE_DIR to point to a MARBL source checkout.
# This file is safe to include() from any cmake context.

if(NOT MARBL_SOURCE_DIR OR NOT IS_DIRECTORY "${MARBL_SOURCE_DIR}")
    message(FATAL_ERROR
        "MARBLTargets.cmake requires MARBL_SOURCE_DIR to be set to a valid MARBL source directory.\n"
        "Set MARBL_ROOT in the environment or pass -DMARBL_SOURCE_DIR=/path/to/MARBL")
endif()

set(_marbl_src "${MARBL_SOURCE_DIR}/src")

add_library(mom6_marbl STATIC
    ${_marbl_src}/marbl_kinds_mod.F90
    ${_marbl_src}/marbl_constants_mod.F90
    ${_marbl_src}/marbl_utils_mod.F90
    ${_marbl_src}/marbl_logging.F90
    ${_marbl_src}/marbl_debug_mod.F90
    ${_marbl_src}/marbl_timing_mod.F90
    ${_marbl_src}/marbl_interface_constants.F90
    ${_marbl_src}/marbl_interface_public_types.F90
    ${_marbl_src}/marbl_interface_private_types.F90
    ${_marbl_src}/marbl_diagnostics_share_mod.F90
    ${_marbl_src}/marbl_interior_tendency_share_mod.F90
    ${_marbl_src}/marbl_surface_flux_share_mod.F90
    ${_marbl_src}/marbl_saved_state_mod.F90
    ${_marbl_src}/marbl_settings_mod.F90
    ${_marbl_src}/marbl_pft_mod.F90
    ${_marbl_src}/marbl_glo_avg_mod.F90
    ${_marbl_src}/marbl_oxygen.F90
    ${_marbl_src}/marbl_temperature.F90
    ${_marbl_src}/marbl_schmidt_number_mod.F90
    ${_marbl_src}/marbl_co2calc_mod.F90
    ${_marbl_src}/marbl_nhx_surface_emis_mod.F90
    ${_marbl_src}/marbl_diagnostics_mod.F90
    ${_marbl_src}/marbl_init_tracer_metadata_mod.F90
    ${_marbl_src}/marbl_abio_dic_diagnostics_mod.F90
    ${_marbl_src}/marbl_ciso_diagnostics_mod.F90
    ${_marbl_src}/marbl_abio_dic_interior_tendency_mod.F90
    ${_marbl_src}/marbl_abio_dic_surface_flux_mod.F90
    ${_marbl_src}/marbl_ciso_interior_tendency_mod.F90
    ${_marbl_src}/marbl_ciso_surface_flux_mod.F90
    ${_marbl_src}/marbl_interior_tendency_mod.F90
    ${_marbl_src}/marbl_surface_flux_mod.F90
    ${_marbl_src}/marbl_restore_mod.F90
    ${_marbl_src}/marbl_init_mod.F90
    ${_marbl_src}/marbl_interface.F90
)
add_library(MARBL::marbl ALIAS mom6_marbl)

set_target_properties(mom6_marbl PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/marbl/mod)
target_include_directories(mom6_marbl
    INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/marbl/mod>)

unset(_marbl_src)
