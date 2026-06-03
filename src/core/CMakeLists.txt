# Remaining core files (MOM_grid and MOM_verticalGrid are in mom6_grid below).
target_sources(mom6_ocean PRIVATE
    MOM.F90
    MOM_CoriolisAdv.F90
    MOM_PressureForce.F90
    MOM_PressureForce_FV.F90
    MOM_PressureForce_Montgomery.F90
    MOM_barotropic.F90
    MOM_boundary_update.F90
    MOM_check_scaling.F90
    MOM_checksum_packages.F90
    MOM_continuity.F90
    MOM_continuity_PPM.F90
    MOM_density_integrals.F90
    MOM_dynamics_split_RK2.F90
    MOM_dynamics_split_RK2b.F90
    MOM_dynamics_unsplit.F90
    MOM_dynamics_unsplit_RK2.F90
    MOM_forcing_type.F90
    MOM_interface_heights.F90
    MOM_isopycnal_slopes.F90
    MOM_open_boundary.F90
    MOM_porous_barriers.F90
    MOM_stoch_eos.F90
    MOM_transcribe_grid.F90
    MOM_unit_tests.F90
    MOM_variables.F90
)

# Layer 4 (mom6_grid): MOM_grid and MOM_verticalGrid.
# These are the minimal src/core/ files needed to complete the framework layer.
# MOM_io (Layer 5) needs them; they in turn only need mom6_framework (Layer 3).
add_library(mom6_grid STATIC
    MOM_grid.F90
    MOM_verticalGrid.F90
)
add_library(MOM6::grid ALIAS mom6_grid)

set_target_properties(mom6_grid PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod
)
target_include_directories(mom6_grid
    INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
    PRIVATE   ${TURBO_MEMORY_DIR}
)
target_link_libraries(mom6_grid
    PUBLIC MOM6::framework
)
