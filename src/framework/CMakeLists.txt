# ── Why multiple framework libraries? ─────────────────────────────────────────
#
# MOM6's src/framework/ contains files at several dependency levels. Splitting
# into layers avoids circular dependencies between libraries:
#
#   mom6_framework needs mom6_grid needs mom6_framework   (CYCLE if one library)
#
# Layer 1 (this file):     mom6_framework_base -- leaf utilities, no project deps.
# Layer 2 (config_src/infra): mom6_infra -- FMS2/TIM infra wrapper.
# Layer 3 (this file):     mom6_framework -- framework core, needs Layers 1 and 2.
# Layer 4 (src/core):      mom6_grid -- MOM_grid and MOM_verticalGrid.
# Layer 5 (this file):     mom6_io -- I/O modules that need mom6_grid.
#
# ─────────────────────────────────────────────────────────────────────────────

# Layer 1: pure utility files, no dependencies on anything in this project.
add_library(mom6_framework_base STATIC
    MOM_string_functions.F90
    MOM_array_transform.F90
    MOM_intrinsic_functions.F90
    MOM_murmur_hash.F90
    MOM_safe_alloc.F90
    posix.F90
    numerical_testing_type.F90
)
add_library(MOM6::framework_base ALIAS mom6_framework_base)

set_target_properties(mom6_framework_base PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod/base
)
target_include_directories(mom6_framework_base
    INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod/base>
              $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

# Layer 3: higher-level framework modules that depend on mom6_infra.
# Files added in dependency order via target_sources(), one block per tier.

# Tier A -- depend only on infra + framework_base
add_library(mom6_framework STATIC
    MOM_error_handler.F90
    MOM_cpu_clock.F90
    MOM_document.F90
    MOM_unique_scales.F90
    MOM_coms.F90
    MOM_coupler_types.F90
    MOM_data_override.F90
    MOM_ensemble_manager.F90
    MOM_unit_testing.F90
)

# Tier B
target_sources(mom6_framework PRIVATE MOM_file_parser.F90)

# Tier C
target_sources(mom6_framework PRIVATE
    MOM_unit_scaling.F90
    MOM_domains.F90
)

# Tier D
target_sources(mom6_framework PRIVATE MOM_hor_index.F90)

# Tier E
target_sources(mom6_framework PRIVATE
    MOM_dyn_horgrid.F90
    MOM_random.F90
    MOM_checksums.F90
    MOM_netcdf.F90
    MOM_io_file.F90
)

# Three src/framework/ files are compiled as part of mom6_ocean (src/CMakeLists.txt)
# because they depend on modules outside src/framework/ and src/core/:
#   MOM_diag_mediator       needs MOM_EOS, MOM_interface_heights, MOM_variables
#   MOM_diag_remap          needs MOM_EOS, MOM_regridding, MOM_remapping, coord_* modules
#   MOM_horizontal_regridding needs MOM_debugging (src/diagnostics/)

add_library(MOM6::framework ALIAS mom6_framework)

set_target_properties(mom6_framework PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod/framework
)
target_include_directories(mom6_framework
    INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod/framework>
)
target_link_libraries(mom6_framework
    PUBLIC  MOM6::infra MOM6::framework_base
)

# Layer 5 (mom6_io): I/O modules that need mom6_grid (Layer 4).
# mom6_grid is defined in src/core/CMakeLists.txt, processed first because
# core/ is listed before framework/ in src/CMakeLists.txt.

# Tier A
add_library(mom6_io STATIC MOM_io.F90)

# Tier B
target_sources(mom6_io PRIVATE
    MOM_get_input.F90
    MOM_write_cputime.F90
)

# Tier C
target_sources(mom6_io PRIVATE
    MOM_interpolate.F90
    MOM_restart.F90
)

add_library(MOM6::io ALIAS mom6_io)

set_target_properties(mom6_io PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod/io
)
target_include_directories(mom6_io
    INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod/io>
)
target_link_libraries(mom6_io
    PUBLIC MOM6::grid
)
