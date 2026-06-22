# All core files compile into mom6_ocean. MOM_grid and MOM_verticalGrid used to
# be a separate mom6_grid library that existed only to order compilation; that
# ordering is now handled automatically by CMake's Fortran scanner within
# mom6_ocean, so they are just ordinary mom6_ocean sources here.
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
    MOM_grid.F90
    MOM_verticalGrid.F90
)
