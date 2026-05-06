add_executable(MOM6
    MOM_driver.F90
    MOM_surface_forcing.F90
    MESO_surface_forcing.F90
    atmos_ocean_fluxes.F90
    user_surface_forcing.F90
)

set_target_properties(MOM6 PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod
)
target_include_directories(MOM6
    PRIVATE ${MOM6_SOURCE_DIR}/config_src/memory/${TURBO_MEMORY_MODE}
)
target_link_libraries(MOM6
    PRIVATE MOM6::ocean
)
