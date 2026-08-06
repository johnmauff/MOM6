! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Solve the layer continuity equation.
module MOM_continuity

use MOM_continuity_PPM, only : continuity_PPM
use MOM_continuity_PPM, only : continuity_stencil=>continuity_PPM_stencil
use MOM_continuity_PPM, only : continuity_init=>continuity_PPM_init
use MOM_continuity_PPM, only : continuity_CS=>continuity_PPM_CS
use MOM_continuity_PPM, only : continuity_fluxes, continuity_adjust_vel
use MOM_continuity_PPM, only : zonal_mass_flux, meridional_mass_flux
use MOM_continuity_PPM, only : zonal_edge_thickness, meridional_edge_thickness
use MOM_continuity_PPM, only : continuity_zonal_convergence, continuity_meridional_convergence
use MOM_continuity_PPM, only : zonal_flux_thickness, meridional_flux_thickness
use MOM_continuity_PPM, only : zonal_BT_mass_flux, meridional_BT_mass_flux
use MOM_continuity_PPM, only : set_continuity_loop_bounds, cont_loop_bounds_type
use MOM_continuity_PPM, only : set_continuity_box

use MOM_grid, only : ocean_grid_type
use MOM_open_boundary, only : ocean_OBC_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : BT_cont_type, porous_barrier_type
use MOM_verticalGrid, only : verticalGrid_type

use array_mod, only : RealArray_t

#include <MOM_memory.h>

implicit none ; private

! These are direct pass-throughs of routines in continuity_PPM
public continuity, continuity_init, continuity_stencil, continuity_CS
public continuity_fluxes, continuity_adjust_vel
public zonal_mass_flux, meridional_mass_flux
public zonal_edge_thickness, meridional_edge_thickness
public continuity_zonal_convergence, continuity_meridional_convergence
public zonal_flux_thickness, meridional_flux_thickness
public zonal_BT_mass_flux, meridional_BT_mass_flux
public set_continuity_loop_bounds, cont_loop_bounds_type
public set_continuity_box

contains

!> Time steps the layer thicknesses, using a monotonically limit, directionally split PPM scheme,
!! based on Lin (1994). A pure forwarding wrapper around continuity_PPM.
subroutine continuity(u, v, hin, h, uh, vh, dt, G, GV, US, CS, OBC, pbv, uhbt_a, vhbt_a, &
                      visc_rem_u_a, visc_rem_v_a, u_cor_a, v_cor_a, BT_cont, du_cor_a, dv_cor_a)
  type(ocean_grid_type),   intent(in)    :: G   !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV  !< Vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: u   !< Zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in)    :: v   !< Meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: hin !< Initial layer thickness [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: h   !< Final layer thickness [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(out)   :: uh  !< Zonal volume flux, u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(out)   :: vh  !< Meridional volume flux, v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1].
  real,                    intent(in)    :: dt  !< Time increment [T ~> s].
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  type(continuity_CS),     intent(in)    :: CS  !< Module's control structure.
  type(ocean_OBC_type),    pointer       :: OBC !< Open boundaries control structure.
  type(porous_barrier_type), intent(in)  :: pbv !< pointers to porous barrier fractional cell metrics
  type(RealArray_t),       intent(in)    :: uhbt_a !< The summed volume flux through zonal faces
                                                       !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(RealArray_t),       intent(in)    :: vhbt_a !< The summed volume flux through meridional
                                                       !! faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(RealArray_t),       intent(in)    :: visc_rem_u_a
                             !< The fraction of zonal momentum originally
                             !! in a layer that remains after a time-step of viscosity, and the
                             !! fraction of a time-step's worth of a barotropic acceleration that
                             !! a layer experiences after viscosity is applied [nondim].
                             !! Visc_rem_u is between 0 (at the bottom) and 1 (far above the bottom).
  type(RealArray_t),       intent(in)    :: visc_rem_v_a
                             !< The fraction of meridional momentum originally
                             !! in a layer that remains after a time-step of viscosity, and the
                             !! fraction of a time-step's worth of a barotropic acceleration that
                             !! a layer experiences after viscosity is applied [nondim].
                             !! Visc_rem_v is between 0 (at the bottom) and 1 (far above the bottom).
  type(RealArray_t),       intent(inout) :: u_cor_a
                             !< The zonal velocities that give uhbt as the depth-integrated
                             !! transport [L T-1 ~> m s-1].
  type(RealArray_t),       intent(inout) :: v_cor_a
                             !< The meridional velocities that give vhbt as the depth-integrated
                             !! transport [L T-1 ~> m s-1].
  type(BT_cont_type), optional, pointer  :: BT_cont !< A structure with elements that describe
                             !!  the effective open face areas as a function of barotropic flow.
  type(RealArray_t),       intent(inout) :: du_cor_a !< The zonal velocity increments from u
                                                         !! that give uhbt as the depth-integrated
                                                         !! transports [L T-1 ~> m s-1].
  type(RealArray_t),       intent(inout) :: dv_cor_a !< The meridional velocity increments from
                                                         !! v that give vhbt as the depth-integrated
                                                         !! transports [L T-1 ~> m s-1].

  call continuity_PPM(u, v, hin, h, uh, vh, dt, G, GV, US, CS, OBC, pbv, uhbt_a, vhbt_a, &
                      visc_rem_u_a, visc_rem_v_a, u_cor_a, v_cor_a, BT_cont, du_cor_a, dv_cor_a)

end subroutine continuity

end module MOM_continuity
