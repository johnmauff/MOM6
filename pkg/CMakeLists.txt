set(_CVMIX  "${MOM6_SOURCE_DIR}/pkg/CVMix-src/src/shared")
set(_GSW_MOD "${MOM6_SOURCE_DIR}/pkg/GSW-Fortran/modules")
set(_GSW_TB  "${MOM6_SOURCE_DIR}/pkg/GSW-Fortran/toolbox")

# ── CVMix ──────────────────────────────────────────────────────────────────
add_library(mom6_cvmix STATIC
    ${_CVMIX}/cvmix_kinds_and_types.F90
    ${_CVMIX}/cvmix_background.F90
    ${_CVMIX}/cvmix_convection.F90
    ${_CVMIX}/cvmix_ddiff.F90
    ${_CVMIX}/cvmix_kpp.F90
    ${_CVMIX}/cvmix_math.F90
    ${_CVMIX}/cvmix_put_get.F90
    ${_CVMIX}/cvmix_shear.F90
    ${_CVMIX}/cvmix_tidal.F90
    ${_CVMIX}/cvmix_utils.F90
)
add_library(MOM6::CVMix ALIAS mom6_cvmix)
set_target_properties(mom6_cvmix PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod/cvmix)
target_include_directories(mom6_cvmix
    INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod/cvmix>)

# ── GSW ────────────────────────────────────────────────────────────────────
add_library(mom6_gsw STATIC
    ${_GSW_MOD}/gsw_mod_kinds.f90
    ${_GSW_MOD}/gsw_mod_baltic_data.f90
    ${_GSW_MOD}/gsw_mod_error_functions.f90
    ${_GSW_MOD}/gsw_mod_freezing_poly_coefficients.f90
    ${_GSW_MOD}/gsw_mod_gibbs_ice_coefficients.f90
    ${_GSW_MOD}/gsw_mod_saar_data.f90
    ${_GSW_MOD}/gsw_mod_sp_coefficients.f90
    ${_GSW_MOD}/gsw_mod_specvol_coefficients.f90
    ${_GSW_MOD}/gsw_mod_teos10_constants.f90
    ${_GSW_MOD}/gsw_mod_netcdf.f90
    ${_GSW_MOD}/gsw_mod_toolbox.f90
)
add_library(MOM6::GSW ALIAS mom6_gsw)
set_target_properties(mom6_gsw PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod/gsw)
target_include_directories(mom6_gsw
    INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod/gsw>)
# gsw_mod_netcdf.f90 uses NetCDF; NetCDF_C needed because nf-config --flibs
# omits the C library search path.
target_link_libraries(mom6_gsw
    PRIVATE NetCDF::NetCDF_Fortran NetCDF::NetCDF_C)
target_sources(mom6_gsw PRIVATE
    ${_GSW_TB}/gsw_add_barrier.f90
    ${_GSW_TB}/gsw_add_mean.f90
    ${_GSW_TB}/gsw_adiabatic_lapse_rate_from_ct.f90
    ${_GSW_TB}/gsw_adiabatic_lapse_rate_ice.f90
    ${_GSW_TB}/gsw_alpha.f90
    ${_GSW_TB}/gsw_alpha_on_beta.f90
    ${_GSW_TB}/gsw_alpha_wrt_t_exact.f90
    ${_GSW_TB}/gsw_alpha_wrt_t_ice.f90
    ${_GSW_TB}/gsw_beta.f90
    ${_GSW_TB}/gsw_beta_const_t_exact.f90
    ${_GSW_TB}/gsw_c_from_sp.f90
    ${_GSW_TB}/gsw_cabbeling.f90
    ${_GSW_TB}/gsw_chem_potential_water_ice.f90
    ${_GSW_TB}/gsw_chem_potential_water_t_exact.f90
    ${_GSW_TB}/gsw_cp_ice.f90
    ${_GSW_TB}/gsw_ct_first_derivatives.f90
    ${_GSW_TB}/gsw_ct_first_derivatives_wrt_t_exact.f90
    ${_GSW_TB}/gsw_ct_freezing.f90
    ${_GSW_TB}/gsw_ct_freezing_exact.f90
    ${_GSW_TB}/gsw_ct_freezing_first_derivatives.f90
    ${_GSW_TB}/gsw_ct_freezing_first_derivatives_poly.f90
    ${_GSW_TB}/gsw_ct_freezing_poly.f90
    ${_GSW_TB}/gsw_ct_from_enthalpy.f90
    ${_GSW_TB}/gsw_ct_from_enthalpy_exact.f90
    ${_GSW_TB}/gsw_ct_from_entropy.f90
    ${_GSW_TB}/gsw_ct_from_pt.f90
    ${_GSW_TB}/gsw_ct_from_rho.f90
    ${_GSW_TB}/gsw_ct_from_t.f90
    ${_GSW_TB}/gsw_ct_maxdensity.f90
    ${_GSW_TB}/gsw_ct_second_derivatives.f90
    ${_GSW_TB}/gsw_deltasa_atlas.f90
    ${_GSW_TB}/gsw_deltasa_from_sp.f90
    ${_GSW_TB}/gsw_dilution_coefficient_t_exact.f90
    ${_GSW_TB}/gsw_dynamic_enthalpy.f90
    ${_GSW_TB}/gsw_enthalpy.f90
    ${_GSW_TB}/gsw_enthalpy_ct_exact.f90
    ${_GSW_TB}/gsw_enthalpy_diff.f90
    ${_GSW_TB}/gsw_enthalpy_first_derivatives.f90
    ${_GSW_TB}/gsw_enthalpy_first_derivatives_ct_exact.f90
    ${_GSW_TB}/gsw_enthalpy_ice.f90
    ${_GSW_TB}/gsw_enthalpy_second_derivatives.f90
    ${_GSW_TB}/gsw_enthalpy_second_derivatives_ct_exact.f90
    ${_GSW_TB}/gsw_enthalpy_sso_0.f90
    ${_GSW_TB}/gsw_enthalpy_t_exact.f90
    ${_GSW_TB}/gsw_entropy_first_derivatives.f90
    ${_GSW_TB}/gsw_entropy_from_pt.f90
    ${_GSW_TB}/gsw_entropy_from_t.f90
    ${_GSW_TB}/gsw_entropy_ice.f90
    ${_GSW_TB}/gsw_entropy_part.f90
    ${_GSW_TB}/gsw_entropy_part_zerop.f90
    ${_GSW_TB}/gsw_entropy_second_derivatives.f90
    ${_GSW_TB}/gsw_fdelta.f90
    ${_GSW_TB}/gsw_frazil_properties.f90
    ${_GSW_TB}/gsw_frazil_properties_potential.f90
    ${_GSW_TB}/gsw_frazil_properties_potential_poly.f90
    ${_GSW_TB}/gsw_frazil_ratios_adiabatic.f90
    ${_GSW_TB}/gsw_frazil_ratios_adiabatic_poly.f90
    ${_GSW_TB}/gsw_geo_strf_dyn_height.f90
    ${_GSW_TB}/gsw_geo_strf_dyn_height_pc.f90
    ${_GSW_TB}/gsw_gibbs.f90
    ${_GSW_TB}/gsw_gibbs_ice.f90
    ${_GSW_TB}/gsw_gibbs_ice_part_t.f90
    ${_GSW_TB}/gsw_gibbs_ice_pt0.f90
    ${_GSW_TB}/gsw_gibbs_ice_pt0_pt0.f90
    ${_GSW_TB}/gsw_gibbs_pt0_pt0.f90
    ${_GSW_TB}/gsw_grav.f90
    ${_GSW_TB}/gsw_helmholtz_energy_ice.f90
    ${_GSW_TB}/gsw_hill_ratio_at_sp2.f90
    ${_GSW_TB}/gsw_ice_fraction_to_freeze_seawater.f90
    ${_GSW_TB}/gsw_internal_energy.f90
    ${_GSW_TB}/gsw_internal_energy_ice.f90
    ${_GSW_TB}/gsw_ipv_vs_fnsquared_ratio.f90
    ${_GSW_TB}/gsw_kappa.f90
    ${_GSW_TB}/gsw_kappa_const_t_ice.f90
    ${_GSW_TB}/gsw_kappa_ice.f90
    ${_GSW_TB}/gsw_kappa_t_exact.f90
    ${_GSW_TB}/gsw_latentheat_evap_ct.f90
    ${_GSW_TB}/gsw_latentheat_evap_t.f90
    ${_GSW_TB}/gsw_latentheat_melting.f90
    ${_GSW_TB}/gsw_linear_interp_sa_ct.f90
    ${_GSW_TB}/gsw_melting_ice_equilibrium_sa_ct_ratio.f90
    ${_GSW_TB}/gsw_melting_ice_equilibrium_sa_ct_ratio_poly.f90
    ${_GSW_TB}/gsw_melting_ice_into_seawater.f90
    ${_GSW_TB}/gsw_melting_ice_sa_ct_ratio.f90
    ${_GSW_TB}/gsw_melting_ice_sa_ct_ratio_poly.f90
    ${_GSW_TB}/gsw_melting_seaice_equilibrium_sa_ct_ratio.f90
    ${_GSW_TB}/gsw_melting_seaice_equilibrium_sa_ct_ratio_poly.f90
    ${_GSW_TB}/gsw_melting_seaice_into_seawater.f90
    ${_GSW_TB}/gsw_melting_seaice_sa_ct_ratio.f90
    ${_GSW_TB}/gsw_melting_seaice_sa_ct_ratio_poly.f90
    ${_GSW_TB}/gsw_nsquared.f90
    ${_GSW_TB}/gsw_pot_enthalpy_from_pt_ice.f90
    ${_GSW_TB}/gsw_pot_enthalpy_from_pt_ice_poly.f90
    ${_GSW_TB}/gsw_pot_enthalpy_ice_freezing.f90
    ${_GSW_TB}/gsw_pot_enthalpy_ice_freezing_first_derivatives.f90
    ${_GSW_TB}/gsw_pot_enthalpy_ice_freezing_first_derivatives_poly.f90
    ${_GSW_TB}/gsw_pot_enthalpy_ice_freezing_poly.f90
    ${_GSW_TB}/gsw_pot_rho_t_exact.f90
    ${_GSW_TB}/gsw_pressure_coefficient_ice.f90
    ${_GSW_TB}/gsw_pressure_freezing_ct.f90
    ${_GSW_TB}/gsw_pt0_cold_ice_poly.f90
    ${_GSW_TB}/gsw_pt0_from_t.f90
    ${_GSW_TB}/gsw_pt0_from_t_ice.f90
    ${_GSW_TB}/gsw_pt_first_derivatives.f90
    ${_GSW_TB}/gsw_pt_from_ct.f90
    ${_GSW_TB}/gsw_pt_from_entropy.f90
    ${_GSW_TB}/gsw_pt_from_pot_enthalpy_ice.f90
    ${_GSW_TB}/gsw_pt_from_pot_enthalpy_ice_poly.f90
    ${_GSW_TB}/gsw_pt_from_pot_enthalpy_ice_poly_dh.f90
    ${_GSW_TB}/gsw_pt_from_t.f90
    ${_GSW_TB}/gsw_pt_from_t_ice.f90
    ${_GSW_TB}/gsw_pt_second_derivatives.f90
    ${_GSW_TB}/gsw_rho.f90
    ${_GSW_TB}/gsw_rho_alpha_beta.f90
    ${_GSW_TB}/gsw_rho_first_derivatives.f90
    ${_GSW_TB}/gsw_rho_first_derivatives_wrt_enthalpy.f90
    ${_GSW_TB}/gsw_rho_ice.f90
    ${_GSW_TB}/gsw_rho_second_derivatives.f90
    ${_GSW_TB}/gsw_rho_second_derivatives_wrt_enthalpy.f90
    ${_GSW_TB}/gsw_rho_t_exact.f90
    ${_GSW_TB}/gsw_rr68_interp_sa_ct.f90
    ${_GSW_TB}/gsw_sa_freezing_estimate.f90
    ${_GSW_TB}/gsw_sa_freezing_from_ct.f90
    ${_GSW_TB}/gsw_sa_freezing_from_ct_poly.f90
    ${_GSW_TB}/gsw_sa_freezing_from_t.f90
    ${_GSW_TB}/gsw_sa_freezing_from_t_poly.f90
    ${_GSW_TB}/gsw_sa_from_rho.f90
    ${_GSW_TB}/gsw_sa_from_sp.f90
    ${_GSW_TB}/gsw_sa_from_sp_baltic.f90
    ${_GSW_TB}/gsw_sa_from_sstar.f90
    ${_GSW_TB}/gsw_sa_p_inrange.f90
    ${_GSW_TB}/gsw_saar.f90
    ${_GSW_TB}/gsw_saar_init.f90
    ${_GSW_TB}/gsw_seaice_fraction_to_freeze_seawater.f90
    ${_GSW_TB}/gsw_sigma0.f90
    ${_GSW_TB}/gsw_sigma1.f90
    ${_GSW_TB}/gsw_sigma2.f90
    ${_GSW_TB}/gsw_sigma3.f90
    ${_GSW_TB}/gsw_sigma4.f90
    ${_GSW_TB}/gsw_sound_speed.f90
    ${_GSW_TB}/gsw_sound_speed_ice.f90
    ${_GSW_TB}/gsw_sound_speed_t_exact.f90
    ${_GSW_TB}/gsw_sp_from_c.f90
    ${_GSW_TB}/gsw_sp_from_sa.f90
    ${_GSW_TB}/gsw_sp_from_sa_baltic.f90
    ${_GSW_TB}/gsw_sp_from_sk.f90
    ${_GSW_TB}/gsw_sp_from_sr.f90
    ${_GSW_TB}/gsw_sp_from_sstar.f90
    ${_GSW_TB}/gsw_specvol.f90
    ${_GSW_TB}/gsw_specvol_alpha_beta.f90
    ${_GSW_TB}/gsw_specvol_anom_standard.f90
    ${_GSW_TB}/gsw_specvol_first_derivatives.f90
    ${_GSW_TB}/gsw_specvol_first_derivatives_wrt_enthalpy.f90
    ${_GSW_TB}/gsw_specvol_ice.f90
    ${_GSW_TB}/gsw_specvol_second_derivatives.f90
    ${_GSW_TB}/gsw_specvol_second_derivatives_wrt_enthalpy.f90
    ${_GSW_TB}/gsw_specvol_sso_0.f90
    ${_GSW_TB}/gsw_specvol_t_exact.f90
    ${_GSW_TB}/gsw_spiciness0.f90
    ${_GSW_TB}/gsw_spiciness1.f90
    ${_GSW_TB}/gsw_spiciness2.f90
    ${_GSW_TB}/gsw_sr_from_sp.f90
    ${_GSW_TB}/gsw_sstar_from_sa.f90
    ${_GSW_TB}/gsw_sstar_from_sp.f90
    ${_GSW_TB}/gsw_t_deriv_chem_potential_water_t_exact.f90
    ${_GSW_TB}/gsw_t_freezing.f90
    ${_GSW_TB}/gsw_t_freezing_exact.f90
    ${_GSW_TB}/gsw_t_freezing_first_derivatives.f90
    ${_GSW_TB}/gsw_t_freezing_first_derivatives_poly.f90
    ${_GSW_TB}/gsw_t_freezing_poly.f90
    ${_GSW_TB}/gsw_t_from_ct.f90
    ${_GSW_TB}/gsw_t_from_pt0_ice.f90
    ${_GSW_TB}/gsw_thermobaric.f90
    ${_GSW_TB}/gsw_turner_rsubrho.f90
    ${_GSW_TB}/gsw_util_indx.f90
    ${_GSW_TB}/gsw_util_interp1q_int.f90
    ${_GSW_TB}/gsw_util_sort_real.f90
    ${_GSW_TB}/gsw_util_xinterp1.f90
    ${_GSW_TB}/gsw_z_from_p.f90
)

# ── MARBL ──────────────────────────────────────────────────────────────────
# MARBL::marbl is NOT built here. It is a TURBO-stack-specific biogeochemistry
# plugin that the parent project (turbo-stack) builds and provides as a
# pre-defined CMake target before calling add_subdirectory on this repo —
# the same pattern as FMS::fms_r8.  src/CMakeLists.txt links against it
# directly.  If MARBL::marbl is not defined at configure time the build will
# fail with a clear "target not found" error.
