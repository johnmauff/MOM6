# ── mom6_framework: model-agnostic framework ──────────────────────────────────
#
# One library for the src/framework/ code that the model (mom6_ocean) builds on
# but that does not itself depend on the model. There are no sub-libraries or
# tiers: every file compiles into a single target and CMake's Fortran dependency
# scanner orders them automatically from their module graph.
#
# Adding a src/framework/ file: if it `use`s only mom6_infra and other
# mom6_framework modules, list it here. If it `use`s grid/io or any model module
# (MOM_EOS, MOM_regridding, diagnostics, core, ...), compile it into mom6_ocean
# instead (src/CMakeLists.txt). When unsure, grep its `use` statements.
# ──────────────────────────────────────────────────────────────────────────────

add_library(mom6_framework STATIC
    # leaf utilities — no dependencies on anything else in this project
    MOM_string_functions.F90
    MOM_array_transform.F90
    MOM_intrinsic_functions.F90
    MOM_murmur_hash.F90
    MOM_safe_alloc.F90
    posix.F90
    numerical_testing_type.F90
    # framework core — depend on mom6_infra and/or the utilities above
    MOM_error_handler.F90
    MOM_cpu_clock.F90
    MOM_document.F90
    MOM_unique_scales.F90
    MOM_coms.F90
    MOM_coupler_types.F90
    MOM_data_override.F90
    MOM_ensemble_manager.F90
    MOM_unit_testing.F90
    MOM_file_parser.F90
    MOM_unit_scaling.F90
    MOM_domains.F90
    MOM_hor_index.F90
    MOM_dyn_horgrid.F90
    MOM_random.F90
    MOM_checksums.F90
    MOM_netcdf.F90
    MOM_io_file.F90
)
add_library(MOM6::framework ALIAS mom6_framework)

set_target_properties(mom6_framework PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod
)
target_include_directories(mom6_framework
    PUBLIC    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
)
target_link_libraries(mom6_framework
    PUBLIC MOM6::infra
)
