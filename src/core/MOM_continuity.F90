! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0
!!SKILLS: 0.3.1

!> Solve the layer continuity equation.
module MOM_continuity

use MOM_continuity_PPM, only : continuity_PPM
use MOM_continuity_PPM, only : continuity_stencil=>continuity_PPM_stencil
use MOM_continuity_PPM, only : continuity_init=>continuity_PPM_init
use MOM_continuity_PPM, only : continuity_CS=>continuity_PPM_CS
use MOM_continuity_PPM, only : continuity_PPM_3d_fluxes
use MOM_continuity_PPM, only : continuity_PPM_2d_fluxes
use MOM_continuity_PPM, only : continuity_PPM_adjust_vel
use MOM_continuity_PPM, only : zonal_mass_flux, meridional_mass_flux
use MOM_continuity_PPM, only : zonal_edge_thickness, meridional_edge_thickness
use MOM_continuity_PPM, only : continuity_zonal_convergence, continuity_meridional_convergence
use MOM_continuity_PPM, only : zonal_flux_thickness, meridional_flux_thickness
use MOM_continuity_PPM, only : zonal_BT_mass_flux, meridional_BT_mass_flux
use MOM_continuity_PPM, only : set_continuity_loop_bounds, cont_loop_bounds_type
use box_mod, only : Box_t
use MOM_grid, only : ocean_grid_type
use MOM_open_boundary, only : ocean_OBC_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : BT_cont_type, porous_barrier_type
use MOM_verticalGrid, only : verticalGrid_type
use array_mod, only : RealArray_t

implicit none ; private

#include <MOM_memory.h>

! These are direct pass-throughs of routines in continuity_PPM
public continuity, continuity_init, continuity_stencil, continuity_CS
public continuity_fluxes, continuity_adjust_vel
public zonal_mass_flux, meridional_mass_flux
public zonal_edge_thickness, meridional_edge_thickness
public continuity_zonal_convergence, continuity_meridional_convergence
public zonal_flux_thickness, meridional_flux_thickness
public zonal_BT_mass_flux, meridional_BT_mass_flux
public set_continuity_loop_bounds, cont_loop_bounds_type

!> Finds the thickness fluxes from the continuity solver or their vertical sum without
!! actually updating the layer thicknesses.
interface continuity_fluxes
  module procedure continuity_3d_fluxes, continuity_2d_fluxes
end interface continuity_fluxes

contains

!> Time steps the layer thicknesses, using a monotonically limit, directionally split PPM scheme,
!! based on Lin (1994).
subroutine continuity(u, v, hin, h, uh, vh, dt, G, GV, US, CS, OBC, pbv, uhbt, vhbt, &
                       visc_rem_u, visc_rem_v, u_cor, v_cor, BT_cont, du_cor, dv_cor)
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
  real, dimension(SZIB_(G),SZJ_(G)), &
                 optional, intent(in)    :: uhbt !< The summed volume flux through zonal faces
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G)), &
                 optional, intent(in)    :: vhbt !< The summed volume flux through meridional faces
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(in)    :: visc_rem_u
                             !< The fraction of zonal momentum originally
                             !! in a layer that remains after a time-step of viscosity, and the
                             !! fraction of a time-step's worth of a barotropic acceleration that
                             !! a layer experiences after viscosity is applied [nondim].
                             !! Visc_rem_u is between 0 (at the bottom) and 1 (far above the bottom).
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                 optional, intent(in)    :: visc_rem_v
                             !< The fraction of meridional momentum originally
                             !! in a layer that remains after a time-step of viscosity, and the
                             !! fraction of a time-step's worth of a barotropic acceleration that
                             !! a layer experiences after viscosity is applied [nondim].
                             !! Visc_rem_v is between 0 (at the bottom) and 1 (far above the bottom).
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(out)   :: u_cor
                             !< The zonal velocities that give uhbt as the depth-integrated transport [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                 optional, intent(out)   :: v_cor
                             !< The meridional velocities that give vhbt as the depth-integrated
                             !! transport [L T-1 ~> m s-1].
  type(BT_cont_type), optional, pointer  :: BT_cont !< A structure with elements that describe
                             !!  the effective open face areas as a function of barotropic flow.
  real, dimension(SZIB_(G),SZJ_(G)), &
                 optional, intent(out)   :: du_cor !< The zonal velocity increments from u that give uhbt
                                                 !! as the depth-integrated transports [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G)), &
                 optional, intent(out)   :: dv_cor !< The meridional velocity increments from v that give vhbt
                                                 !! as the depth-integrated transports [L T-1 ~> m s-1].

  type(RealArray_t) :: u_a, v_a, hin_a, h_a, uh_a, vh_a
  type(RealArray_t) :: uhbt_a, vhbt_a, visc_rem_u_a, visc_rem_v_a
  type(RealArray_t) :: u_cor_a, v_cor_a, du_cor_a, dv_cor_a
  type(RealArray_t) :: mask2dT_a, dy_Cu_a, IareaT_a, IdxT_a, areaT_a, dxT_a
  type(RealArray_t) :: mask2dCu_a, dxCu_a, dx_Cv_a, IdyT_a, dyT_a, mask2dCv_a, dyCv_a
  type(Box_t) :: bx0
  integer :: stencil
  logical :: x_first
  type(BT_cont_type), pointer :: BT_cont_local ! continuity_PPM's BT_cont is mandatory (a
                                               ! pointer, possibly disassociated); this caller's
                                               ! own BT_cont stays optional (fixed external API).

  call u_a%alloc(lb=LBOUND(u), ub=UBOUND(u), source=u)
  call v_a%alloc(lb=LBOUND(v), ub=UBOUND(v), source=v)
  call hin_a%alloc(lb=LBOUND(hin), ub=UBOUND(hin), source=hin)
  call h_a%alloc(lb=LBOUND(h), ub=UBOUND(h), source=h)
  call uh_a%alloc(lb=LBOUND(uh), ub=UBOUND(uh), source=uh)
  call vh_a%alloc(lb=LBOUND(vh), ub=UBOUND(vh), source=vh)

  call bx0%safe_alloc(ndims=3)
  call bx0%set(idxS=[G%isc,G%jsc,1], idxE=[G%iec,G%jec,GV%ke])
  stencil = continuity_stencil(CS)
  x_first = (MOD(G%first_direction,2) == 0)

  call mask2dT_a%alloc(lb=LBOUND(G%mask2dT), ub=UBOUND(G%mask2dT), source=G%mask2dT)
  call dy_Cu_a%alloc(lb=LBOUND(G%dy_Cu), ub=UBOUND(G%dy_Cu), source=G%dy_Cu)
  call IareaT_a%alloc(lb=LBOUND(G%IareaT), ub=UBOUND(G%IareaT), source=G%IareaT)
  call IdxT_a%alloc(lb=LBOUND(G%IdxT), ub=UBOUND(G%IdxT), source=G%IdxT)
  call areaT_a%alloc(lb=LBOUND(G%areaT), ub=UBOUND(G%areaT), source=G%areaT)
  call dxT_a%alloc(lb=LBOUND(G%dxT), ub=UBOUND(G%dxT), source=G%dxT)
  call mask2dCu_a%alloc(lb=LBOUND(G%mask2dCu), ub=UBOUND(G%mask2dCu), source=G%mask2dCu)
  call dxCu_a%alloc(lb=LBOUND(G%dxCu), ub=UBOUND(G%dxCu), source=G%dxCu)
  call dx_Cv_a%alloc(lb=LBOUND(G%dx_Cv), ub=UBOUND(G%dx_Cv), source=G%dx_Cv)
  call IdyT_a%alloc(lb=LBOUND(G%IdyT), ub=UBOUND(G%IdyT), source=G%IdyT)
  call dyT_a%alloc(lb=LBOUND(G%dyT), ub=UBOUND(G%dyT), source=G%dyT)
  call mask2dCv_a%alloc(lb=LBOUND(G%mask2dCv), ub=UBOUND(G%mask2dCv), source=G%mask2dCv)
  call dyCv_a%alloc(lb=LBOUND(G%dyCv), ub=UBOUND(G%dyCv), source=G%dyCv)

  if (present(uhbt)) call uhbt_a%alloc(lb=LBOUND(uhbt), ub=UBOUND(uhbt), source=uhbt)
  if (present(vhbt)) call vhbt_a%alloc(lb=LBOUND(vhbt), ub=UBOUND(vhbt), source=vhbt)
  if (present(visc_rem_u)) &
    call visc_rem_u_a%alloc(lb=LBOUND(visc_rem_u), ub=UBOUND(visc_rem_u), source=visc_rem_u)
  if (present(visc_rem_v)) &
    call visc_rem_v_a%alloc(lb=LBOUND(visc_rem_v), ub=UBOUND(visc_rem_v), source=visc_rem_v)
  if (present(u_cor)) call u_cor_a%alloc(lb=LBOUND(u_cor), ub=UBOUND(u_cor), source=u_cor)
  if (present(v_cor)) call v_cor_a%alloc(lb=LBOUND(v_cor), ub=UBOUND(v_cor), source=v_cor)
  if (present(du_cor)) call du_cor_a%alloc(lb=LBOUND(du_cor), ub=UBOUND(du_cor), source=du_cor)
  if (present(dv_cor)) call dv_cor_a%alloc(lb=LBOUND(dv_cor), ub=UBOUND(dv_cor), source=dv_cor)

  nullify(BT_cont_local)
  if (present(BT_cont)) BT_cont_local => BT_cont

  call continuity_PPM(u_a, v_a, hin_a, h_a, uh_a, vh_a, dt, bx0, stencil, x_first, &
                       mask2dT_a, dy_Cu_a, IareaT_a, IdxT_a, areaT_a, dxT_a, &
                       mask2dCu_a, dxCu_a, dx_Cv_a, IdyT_a, dyT_a, mask2dCv_a, dyCv_a, &
                       G%isd, G%ied, GV%Angstrom_H, GV%H_subroundoff, CS, OBC, pbv, &
                       uhbt_a=uhbt_a, vhbt_a=vhbt_a, visc_rem_u_a=visc_rem_u_a, &
                       visc_rem_v_a=visc_rem_v_a, u_cor_a=u_cor_a, v_cor_a=v_cor_a, &
                       BT_cont=BT_cont_local, du_cor_a=du_cor_a, dv_cor_a=dv_cor_a)

  call h_a%copy2F(h)
  call uh_a%copy2F(uh)
  call vh_a%copy2F(vh)

  call u_a%free()
  call v_a%free()
  call hin_a%free()
  call h_a%free()
  call uh_a%free()
  call vh_a%free()

  call mask2dT_a%free() ; call dy_Cu_a%free() ; call IareaT_a%free() ; call IdxT_a%free()
  call areaT_a%free() ; call dxT_a%free() ; call mask2dCu_a%free() ; call dxCu_a%free()
  call dx_Cv_a%free() ; call IdyT_a%free() ; call dyT_a%free() ; call mask2dCv_a%free()
  call dyCv_a%free()
  call bx0%free()

  if (u_cor_a%associated()) call u_cor_a%copy2F(u_cor)
  if (v_cor_a%associated()) call v_cor_a%copy2F(v_cor)
  if (du_cor_a%associated()) call du_cor_a%copy2F(du_cor)
  if (dv_cor_a%associated()) call dv_cor_a%copy2F(dv_cor)
  call u_cor_a%free() ; call v_cor_a%free() ; call du_cor_a%free() ; call dv_cor_a%free()
  call uhbt_a%free() ; call vhbt_a%free() ; call visc_rem_u_a%free() ; call visc_rem_v_a%free()

end subroutine continuity

!> Finds the thickness fluxes from the continuity solver without actually updating the
!! layer thicknesses.  Because the fluxes in the two directions are calculated based on the
!! input thicknesses, which are not updated between the direcitons, the fluxes returned here
!! are not the same as those that would be returned by a call to continuity.
subroutine continuity_3d_fluxes(u, v, h, uh, vh, dt, G, GV, US, CS, OBC, pbv)
  type(ocean_grid_type),   intent(inout) :: G   !< Ocean grid structure.
  type(verticalGrid_type), intent(in)    :: GV  !< Vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: u   !< Zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in)    :: v   !< Meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  &
                           intent(in)    :: h   !< Layer thickness [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(out)   :: uh  !< Thickness fluxes through zonal faces,
                                                !! u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(out)   :: vh  !< Thickness fluxes through meridional faces,
                                                !! v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1].
  real,                    intent(in)    :: dt  !< Time increment [T ~> s].
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  type(continuity_CS),     intent(in)    :: CS  !< Control structure for mom_continuity.
  type(ocean_OBC_type),    pointer       :: OBC !< Open boundaries control structure.
  type(porous_barrier_type), intent(in)  :: pbv !< porous barrier fractional cell metrics

  type(RealArray_t) :: u_a, v_a, h_a, uh_a, vh_a
  type(RealArray_t) :: mask2dT_a, dy_Cu_a, IareaT_a, IdxT_a, areaT_a, dxT_a
  type(RealArray_t) :: mask2dCu_a, dxCu_a, dx_Cv_a, IdyT_a, dyT_a, mask2dCv_a, dyCv_a
  type(Box_t) :: bxC

  call u_a%alloc(lb=LBOUND(u), ub=UBOUND(u), source=u)
  call v_a%alloc(lb=LBOUND(v), ub=UBOUND(v), source=v)
  call h_a%alloc(lb=LBOUND(h), ub=UBOUND(h), source=h)
  call uh_a%alloc(lb=LBOUND(uh), ub=UBOUND(uh), source=uh)
  call vh_a%alloc(lb=LBOUND(vh), ub=UBOUND(vh), source=vh)

  call bxC%safe_alloc(ndims=3)
  call bxC%set(idxS=[G%isc,G%jsc,1], idxE=[G%iec,G%jec,GV%ke])

  call mask2dT_a%alloc(lb=LBOUND(G%mask2dT), ub=UBOUND(G%mask2dT), source=G%mask2dT)
  call dy_Cu_a%alloc(lb=LBOUND(G%dy_Cu), ub=UBOUND(G%dy_Cu), source=G%dy_Cu)
  call IareaT_a%alloc(lb=LBOUND(G%IareaT), ub=UBOUND(G%IareaT), source=G%IareaT)
  call IdxT_a%alloc(lb=LBOUND(G%IdxT), ub=UBOUND(G%IdxT), source=G%IdxT)
  call areaT_a%alloc(lb=LBOUND(G%areaT), ub=UBOUND(G%areaT), source=G%areaT)
  call dxT_a%alloc(lb=LBOUND(G%dxT), ub=UBOUND(G%dxT), source=G%dxT)
  call mask2dCu_a%alloc(lb=LBOUND(G%mask2dCu), ub=UBOUND(G%mask2dCu), source=G%mask2dCu)
  call dxCu_a%alloc(lb=LBOUND(G%dxCu), ub=UBOUND(G%dxCu), source=G%dxCu)
  call dx_Cv_a%alloc(lb=LBOUND(G%dx_Cv), ub=UBOUND(G%dx_Cv), source=G%dx_Cv)
  call IdyT_a%alloc(lb=LBOUND(G%IdyT), ub=UBOUND(G%IdyT), source=G%IdyT)
  call dyT_a%alloc(lb=LBOUND(G%dyT), ub=UBOUND(G%dyT), source=G%dyT)
  call mask2dCv_a%alloc(lb=LBOUND(G%mask2dCv), ub=UBOUND(G%mask2dCv), source=G%mask2dCv)
  call dyCv_a%alloc(lb=LBOUND(G%dyCv), ub=UBOUND(G%dyCv), source=G%dyCv)

  call continuity_PPM_3d_fluxes(u_a, v_a, h_a, uh_a, vh_a, dt, bxC, &
                                mask2dT_a, dy_Cu_a, IareaT_a, IdxT_a, areaT_a, dxT_a, &
                                mask2dCu_a, dxCu_a, dx_Cv_a, IdyT_a, dyT_a, mask2dCv_a, &
                                dyCv_a, G%isd, G%ied, GV%Angstrom_H, GV%H_subroundoff, &
                                CS, OBC, pbv)

  call uh_a%copy2F(uh)
  call vh_a%copy2F(vh)

  call u_a%free()
  call v_a%free()
  call h_a%free()
  call uh_a%free()
  call vh_a%free()

  call mask2dT_a%free() ; call dy_Cu_a%free() ; call IareaT_a%free() ; call IdxT_a%free()
  call areaT_a%free() ; call dxT_a%free() ; call mask2dCu_a%free() ; call dxCu_a%free()
  call dx_Cv_a%free() ; call IdyT_a%free() ; call dyT_a%free() ; call mask2dCv_a%free()
  call dyCv_a%free()
  call bxC%free()

end subroutine continuity_3d_fluxes

!> Find the vertical sum of the thickness fluxes from the continuity solver without actually
!! updating the layer thicknesses.  Because the fluxes in the two directions are calculated
!! based on the input thicknesses, which are not updated between the directions, the fluxes
!! returned here are not the same as those that would be returned by a call to continuity.
subroutine continuity_2d_fluxes(u, v, h, uhbt, vhbt, dt, G, GV, US, CS, OBC, pbv)
  type(ocean_grid_type),   intent(inout) :: G   !< Ocean grid structure.
  type(verticalGrid_type), intent(in)    :: GV  !< Vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: u   !< Zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in)    :: v   !< Meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  &
                           intent(in)    :: h   !< Layer thickness [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G)), &
                           intent(out)   :: uhbt !< Vertically summed thickness flux through
                                                !! zonal faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G)), &
                           intent(out)   :: vhbt !< Vertically summed thickness flux through
                                                !! meridional faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  real,                    intent(in)    :: dt  !< Time increment [T ~> s].
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  type(continuity_CS),     intent(in)    :: CS  !< Control structure for mom_continuity.
  type(ocean_OBC_type),    pointer       :: OBC !< Open boundaries control structure.
  type(porous_barrier_type), intent(in)  :: pbv !< porous barrier fractional cell metrics

  type(RealArray_t) :: u_a, v_a, h_a, uhbt_a, vhbt_a
  type(RealArray_t) :: mask2dT_a, dy_Cu_a, IareaT_a, IdxT_a, dx_Cv_a, IdyT_a
  type(Box_t) :: bxC

  call u_a%alloc(lb=LBOUND(u), ub=UBOUND(u), source=u)
  call v_a%alloc(lb=LBOUND(v), ub=UBOUND(v), source=v)
  call h_a%alloc(lb=LBOUND(h), ub=UBOUND(h), source=h)
  call uhbt_a%alloc(lb=LBOUND(uhbt), ub=UBOUND(uhbt), source=uhbt)
  call vhbt_a%alloc(lb=LBOUND(vhbt), ub=UBOUND(vhbt), source=vhbt)

  call bxC%safe_alloc(ndims=3)
  call bxC%set(idxS=[G%isc,G%jsc,1], idxE=[G%iec,G%jec,GV%ke])

  call mask2dT_a%alloc(lb=LBOUND(G%mask2dT), ub=UBOUND(G%mask2dT), source=G%mask2dT)
  call dy_Cu_a%alloc(lb=LBOUND(G%dy_Cu), ub=UBOUND(G%dy_Cu), source=G%dy_Cu)
  call IareaT_a%alloc(lb=LBOUND(G%IareaT), ub=UBOUND(G%IareaT), source=G%IareaT)
  call IdxT_a%alloc(lb=LBOUND(G%IdxT), ub=UBOUND(G%IdxT), source=G%IdxT)
  call dx_Cv_a%alloc(lb=LBOUND(G%dx_Cv), ub=UBOUND(G%dx_Cv), source=G%dx_Cv)
  call IdyT_a%alloc(lb=LBOUND(G%IdyT), ub=UBOUND(G%IdyT), source=G%IdyT)

  call continuity_PPM_2d_fluxes(u_a, v_a, h_a, uhbt_a, vhbt_a, dt, bxC, &
                                mask2dT_a, dy_Cu_a, IareaT_a, IdxT_a, dx_Cv_a, IdyT_a, &
                                GV%Angstrom_H, CS, OBC, pbv)

  call uhbt_a%copy2F(uhbt)
  call vhbt_a%copy2F(vhbt)

  call u_a%free()
  call v_a%free()
  call h_a%free()
  call uhbt_a%free()
  call vhbt_a%free()

  call mask2dT_a%free() ; call dy_Cu_a%free() ; call IareaT_a%free() ; call IdxT_a%free()
  call dx_Cv_a%free() ; call IdyT_a%free()
  call bxC%free()

end subroutine continuity_2d_fluxes

!> Correct the velocities to give the specified depth-integrated transports by applying a
!! barotropic acceleration (subject to viscous drag) to the velocities.
subroutine continuity_adjust_vel(u, v, h, dt, G, GV, US, CS, OBC, pbv, uhbt, vhbt, visc_rem_u, visc_rem_v)
  type(ocean_grid_type),   intent(inout) :: G   !< Ocean grid structure.
  type(verticalGrid_type), intent(in)    :: GV  !< Vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(inout) :: u   !< Zonal velocity, which will be adjusted to
                                                !! give uhbt as the depth-integrated
                                                !! transport [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(inout) :: v   !< Meridional velocity, which will be adjusted
                                                !! to give vhbt as the depth-integrated
                                                !! transport [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  &
                           intent(in)    :: h   !< Layer thickness [H ~> m or kg m-2].
  real,                    intent(in)    :: dt  !< Time increment [T ~> s].
  type(unit_scale_type),   intent(in)    :: US  !< A dimensional unit scaling type
  type(continuity_CS),     intent(in)    :: CS  !< Control structure for mom_continuity.
  type(ocean_OBC_type),    pointer       :: OBC !< Open boundaries control structure.
  type(porous_barrier_type), intent(in)  :: pbv !< porous barrier fractional cell metrics
  real, dimension(SZIB_(G),SZJ_(G)), &
                           intent(in)    :: uhbt !< The vertically summed thickness flux through
                                                !! zonal faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G)), &
                           intent(in)    :: vhbt !< The vertically summed thickness flux through
                                                !! meridional faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(in)    :: visc_rem_u !< Both the fraction of the zonal momentum
                                                !! that remains after a time-step of viscosity, and
                                                !! the fraction of a time-step's worth of a barotropic
                                                !! acceleration that a layer experiences after viscosity
                                                !! is applied [nondim].  This goes between 0 (at the
                                                !! bottom) and 1 (far above the bottom).  When this
                                                !! column is under an ice shelf, this also goes to 0
                                                !! at the top due to the no-slip boundary condition there.
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                 optional, intent(in)    :: visc_rem_v !< Both the fraction of the meridional momentum
                                                !! that remains after a time-step of viscosity, and
                                                !! the fraction of a time-step's worth of a barotropic
                                                !! acceleration that a layer experiences after viscosity
                                                !! is applied [nondim].  This goes between 0 (at the
                                                !! bottom) and 1 (far above the bottom).  When this
                                                !! column is under an ice shelf, this also goes to 0
                                                !! at the top due to the no-slip boundary condition there.

  type(RealArray_t) :: u_a, v_a, h_a, uhbt_a, vhbt_a
  type(RealArray_t) :: visc_rem_u_a, visc_rem_v_a
  type(RealArray_t) :: mask2dT_a, dy_Cu_a, IareaT_a, IdxT_a, areaT_a, dxT_a
  type(RealArray_t) :: mask2dCu_a, dxCu_a, dx_Cv_a, IdyT_a, dyT_a, mask2dCv_a, dyCv_a
  type(Box_t) :: bxC

  call u_a%alloc(lb=LBOUND(u), ub=UBOUND(u), source=u)
  call v_a%alloc(lb=LBOUND(v), ub=UBOUND(v), source=v)
  call h_a%alloc(lb=LBOUND(h), ub=UBOUND(h), source=h)
  call uhbt_a%alloc(lb=LBOUND(uhbt), ub=UBOUND(uhbt), source=uhbt)
  call vhbt_a%alloc(lb=LBOUND(vhbt), ub=UBOUND(vhbt), source=vhbt)

  call bxC%safe_alloc(ndims=3)
  call bxC%set(idxS=[G%isc,G%jsc,1], idxE=[G%iec,G%jec,GV%ke])

  call mask2dT_a%alloc(lb=LBOUND(G%mask2dT), ub=UBOUND(G%mask2dT), source=G%mask2dT)
  call dy_Cu_a%alloc(lb=LBOUND(G%dy_Cu), ub=UBOUND(G%dy_Cu), source=G%dy_Cu)
  call IareaT_a%alloc(lb=LBOUND(G%IareaT), ub=UBOUND(G%IareaT), source=G%IareaT)
  call IdxT_a%alloc(lb=LBOUND(G%IdxT), ub=UBOUND(G%IdxT), source=G%IdxT)
  call areaT_a%alloc(lb=LBOUND(G%areaT), ub=UBOUND(G%areaT), source=G%areaT)
  call dxT_a%alloc(lb=LBOUND(G%dxT), ub=UBOUND(G%dxT), source=G%dxT)
  call mask2dCu_a%alloc(lb=LBOUND(G%mask2dCu), ub=UBOUND(G%mask2dCu), source=G%mask2dCu)
  call dxCu_a%alloc(lb=LBOUND(G%dxCu), ub=UBOUND(G%dxCu), source=G%dxCu)
  call dx_Cv_a%alloc(lb=LBOUND(G%dx_Cv), ub=UBOUND(G%dx_Cv), source=G%dx_Cv)
  call IdyT_a%alloc(lb=LBOUND(G%IdyT), ub=UBOUND(G%IdyT), source=G%IdyT)
  call dyT_a%alloc(lb=LBOUND(G%dyT), ub=UBOUND(G%dyT), source=G%dyT)
  call mask2dCv_a%alloc(lb=LBOUND(G%mask2dCv), ub=UBOUND(G%mask2dCv), source=G%mask2dCv)
  call dyCv_a%alloc(lb=LBOUND(G%dyCv), ub=UBOUND(G%dyCv), source=G%dyCv)

  if (present(visc_rem_u)) &
    call visc_rem_u_a%alloc(lb=LBOUND(visc_rem_u), ub=UBOUND(visc_rem_u), source=visc_rem_u)
  if (present(visc_rem_v)) &
    call visc_rem_v_a%alloc(lb=LBOUND(visc_rem_v), ub=UBOUND(visc_rem_v), source=visc_rem_v)

  call continuity_PPM_adjust_vel(u_a, v_a, h_a, dt, bxC, &
                                 mask2dT_a, dy_Cu_a, IareaT_a, IdxT_a, areaT_a, dxT_a, &
                                 mask2dCu_a, dxCu_a, dx_Cv_a, IdyT_a, dyT_a, mask2dCv_a, &
                                 dyCv_a, G%isd, G%ied, GV%Angstrom_H, GV%H_subroundoff, &
                                 CS, OBC, pbv, uhbt_a, vhbt_a, &
                                 visc_rem_u_a=visc_rem_u_a, visc_rem_v_a=visc_rem_v_a)

  call u_a%copy2F(u)
  call v_a%copy2F(v)

  call u_a%free()
  call v_a%free()
  call h_a%free()
  call uhbt_a%free()
  call vhbt_a%free()
  call visc_rem_u_a%free() ; call visc_rem_v_a%free()

  call mask2dT_a%free() ; call dy_Cu_a%free() ; call IareaT_a%free() ; call IdxT_a%free()
  call areaT_a%free() ; call dxT_a%free() ; call mask2dCu_a%free() ; call dxCu_a%free()
  call dx_Cv_a%free() ; call IdyT_a%free() ; call dyT_a%free() ; call mask2dCv_a%free()
  call dyCv_a%free()
  call bxC%free()

end subroutine continuity_adjust_vel

end module MOM_continuity
