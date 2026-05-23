set(_EXT    "${MOM6_SOURCE_DIR}/config_src/external")
set(_MEMDIR "${MOM6_SOURCE_DIR}/config_src/memory/${TURBO_MEMORY_MODE}")

# mom6_ocean: all MOM6 src/ files not already in the framework/grid/io/infra targets,
# plus config_src/external/ stubs (MARBL excluded — provided by MARBL::marbl).
add_library(mom6_ocean STATIC)
add_library(MOM6::ocean ALIAS mom6_ocean)

set_target_properties(mom6_ocean PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod
)
target_include_directories(mom6_ocean
    INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
    PRIVATE   ${_MEMDIR}
)
target_link_libraries(mom6_ocean
    PUBLIC MOM6::io MOM6::grid MOM6::CVMix MOM6::GSW MARBL::marbl
)

# Framework files deferred from the framework layer; they depend on modules
# defined in mom6_ocean (EOS, core, diagnostics) so they are added here.
target_sources(mom6_ocean PRIVATE
    framework/MOM_ANN.F90
    framework/MOM_diag_buffers.F90
    framework/MOM_diag_mediator.F90
    framework/MOM_diag_remap.F90
    framework/MOM_horizontal_regridding.F90
)

# External stubs (MARBL excluded — provided by MARBL::marbl).
target_sources(mom6_ocean PRIVATE
    ${_EXT}/GFDL_ocean_BGC/MOM_generic_tracer.F90
    ${_EXT}/ODA_hooks/kdtree.f90
    ${_EXT}/ODA_hooks/ocean_da_core.F90
    ${_EXT}/ODA_hooks/ocean_da_types.F90
    ${_EXT}/ODA_hooks/write_ocean_obs.F90
    ${_EXT}/database_comms/MOM_database_comms.F90
    ${_EXT}/database_comms/database_client_interface.F90
    ${_EXT}/drifters/MOM_particles.F90
    ${_EXT}/drifters/MOM_particles_types.F90
    ${_EXT}/stochastic_physics/get_stochy_pattern.F90
    ${_EXT}/stochastic_physics/stochastic_physics.F90
)

# core must come before framework so mom6_grid is defined before mom6_io.
add_subdirectory(core)
add_subdirectory(framework)
add_subdirectory(ALE)
add_subdirectory(diagnostics)
add_subdirectory(equation_of_state)
add_subdirectory(ice_shelf)
add_subdirectory(initialization)
add_subdirectory(ocean_data_assim)
add_subdirectory(parameterizations)
add_subdirectory(tracer)
add_subdirectory(user)
