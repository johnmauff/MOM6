# mom6_infra: the FMS2/TIM infrastructure wrapper. This is the one hard boundary
# in MOM6's build — these sources compile against the selected backend only and
# must not `use` any MOM6 module. Keeping infra a separate target with a minimal
# link surface is what enforces that at compile time.
add_library(mom6_infra STATIC
    ${TURBO_INFRA}/MOM_error_infra.F90
    ${TURBO_INFRA}/MOM_constants.F90
    ${TURBO_INFRA}/MOM_time_manager.F90
    ${TURBO_INFRA}/MOM_cpu_clock_infra.F90
    ${TURBO_INFRA}/MOM_ensemble_manager_infra.F90
    ${TURBO_INFRA}/array_mod.F90
    ${TURBO_INFRA}/MOM_coms_infra.F90
    ${TURBO_INFRA}/MOM_domain_infra.F90
    ${TURBO_INFRA}/MOM_data_override_infra.F90
    ${TURBO_INFRA}/MOM_diag_manager_infra.F90
    ${TURBO_INFRA}/MOM_couplertype_infra.F90
    ${TURBO_INFRA}/MOM_string_infra.F90
    ${TURBO_INFRA}/MOM_io_infra.F90
    ${TURBO_INFRA}/MOM_interp_infra.F90
)
# MOM_coms_helpers.F90 exists only in the TIM infra directory
if(TURBO_INFRA STREQUAL "TIM")
    target_sources(mom6_infra PRIVATE ${TURBO_INFRA}/MOM_coms_helpers.F90)
endif()
add_library(MOM6::infra ALIAS mom6_infra)

set_target_properties(mom6_infra PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod
)
target_include_directories(mom6_infra
    INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
)
target_link_libraries(mom6_infra
    PUBLIC TURBO::infra_r8 MPI::MPI_Fortran NetCDF::NetCDF_Fortran
)
