! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0


!> Solve the layer continuity equation.
module MOM_continuity

use MOM_continuity_PPM, only : continuity_PPM
use MOM_continuity_PPM, only : continuity_stencil=>continuity_PPM_stencil
use MOM_continuity_PPM, only : continuity_init=>continuity_PPM_init
use MOM_continuity_PPM, only : continuity_CS=>continuity_PPM_CS
use MOM_continuity_PPM, only : continuity_2d_fluxes_PPM, continuity_3d_fluxes_PPM
use MOM_continuity_PPM, only : continuity_adjust_vel_PPM
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
use box_mod, only : Box_t

#include <MOM_memory.h>

implicit none ; private

!> Calculates the mass or volume fluxes through the zonal and meridional faces.
interface continuity_fluxes
  module procedure continuity_3d_fluxes, continuity_2d_fluxes
end interface continuity_fluxes

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
!! based on Lin (1994). A wrapper around continuity_PPM that presents MOM6's original raw/optional
!! interface to the dynamics core, marshalling into and out of continuity_PPM's array containers.
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

  ! Mandatory containers forwarded to continuity_PPM, built by copy-in from the raw arguments
  ! above; h_a/uh_a/vh_a are copied back out after the call since continuity_PPM writes into them.
  type(RealArray_t) :: u_a, v_a, hin_a, h_a, uh_a, vh_a
  ! Grid-metadata and porous-barrier containers, promoted out of continuity_PPM's own signature;
  ! sourced from G/pbv, which this wrapper already receives wholesale. Never written; freed with
  ! no copy-back.
  type(RealArray_t) :: mask2dT_a, IareaT_a, dy_Cu_a, IdxT_a, dxCu_a, areaT_a, dxT_a
  type(RealArray_t) :: mask2dCu_a, dx_Cv_a, IdyT_a, dyCv_a, dyT_a, mask2dCv_a
  type(RealArray_t) :: por_face_areaU_a, por_face_areaV_a
  integer :: stencil ! The continuity solver stencil size with the current settings
  ! Containers forwarded to continuity_PPM. Left unassociated (the callee's absence signal)
  ! for whichever of the optional raw arguments above the caller did not supply.
  type(RealArray_t) :: uhbt_a, vhbt_a, visc_rem_u_a, visc_rem_v_a
  type(RealArray_t) :: u_cor_a, v_cor_a, du_cor_a, dv_cor_a
  type(box_t) :: bxC ! The continuity solver's base iteration box, built once here

  call bxC%safe_alloc(ndims=3)
  call bxC%set(idxS=[G%isc, G%jsc, 1], idxE=[G%iec, G%jec, GV%ke])

  stencil = continuity_stencil(CS)

  call u_a%alloc(lb=LBOUND(u), ub=UBOUND(u), source=u)
  call v_a%alloc(lb=LBOUND(v), ub=UBOUND(v), source=v)
  call hin_a%alloc(lb=LBOUND(hin), ub=UBOUND(hin), source=hin)
  call h_a%alloc(lb=LBOUND(h), ub=UBOUND(h), source=h)
  call uh_a%alloc(lb=LBOUND(uh), ub=UBOUND(uh), source=uh)
  call vh_a%alloc(lb=LBOUND(vh), ub=UBOUND(vh), source=vh)

  call mask2dT_a%alloc(lb=LBOUND(G%mask2dT), ub=UBOUND(G%mask2dT), source=G%mask2dT)
  call IareaT_a%alloc(lb=LBOUND(G%IareaT), ub=UBOUND(G%IareaT), source=G%IareaT)
  call dy_Cu_a%alloc(lb=LBOUND(G%dy_Cu), ub=UBOUND(G%dy_Cu), source=G%dy_Cu)
  call IdxT_a%alloc(lb=LBOUND(G%IdxT), ub=UBOUND(G%IdxT), source=G%IdxT)
  call dxCu_a%alloc(lb=LBOUND(G%dxCu), ub=UBOUND(G%dxCu), source=G%dxCu)
  call areaT_a%alloc(lb=LBOUND(G%areaT), ub=UBOUND(G%areaT), source=G%areaT)
  call dxT_a%alloc(lb=LBOUND(G%dxT), ub=UBOUND(G%dxT), source=G%dxT)
  call mask2dCu_a%alloc(lb=LBOUND(G%mask2dCu), ub=UBOUND(G%mask2dCu), source=G%mask2dCu)
  call dx_Cv_a%alloc(lb=LBOUND(G%dx_Cv), ub=UBOUND(G%dx_Cv), source=G%dx_Cv)
  call IdyT_a%alloc(lb=LBOUND(G%IdyT), ub=UBOUND(G%IdyT), source=G%IdyT)
  call dyCv_a%alloc(lb=LBOUND(G%dyCv), ub=UBOUND(G%dyCv), source=G%dyCv)
  call dyT_a%alloc(lb=LBOUND(G%dyT), ub=UBOUND(G%dyT), source=G%dyT)
  call mask2dCv_a%alloc(lb=LBOUND(G%mask2dCv), ub=UBOUND(G%mask2dCv), source=G%mask2dCv)
  call por_face_areaU_a%alloc(lb=LBOUND(pbv%por_face_areaU), ub=UBOUND(pbv%por_face_areaU), &
                              source=pbv%por_face_areaU)
  call por_face_areaV_a%alloc(lb=LBOUND(pbv%por_face_areaV), ub=UBOUND(pbv%por_face_areaV), &
                              source=pbv%por_face_areaV)

  if (present(uhbt))       call uhbt_a%alloc(lb=LBOUND(uhbt), ub=UBOUND(uhbt), source=uhbt)
  if (present(vhbt))       call vhbt_a%alloc(lb=LBOUND(vhbt), ub=UBOUND(vhbt), source=vhbt)
  if (present(visc_rem_u)) &
    call visc_rem_u_a%alloc(lb=LBOUND(visc_rem_u), ub=UBOUND(visc_rem_u), source=visc_rem_u)
  if (present(visc_rem_v)) &
    call visc_rem_v_a%alloc(lb=LBOUND(visc_rem_v), ub=UBOUND(visc_rem_v), source=visc_rem_v)
  if (present(u_cor))  call u_cor_a%alloc(lb=LBOUND(u_cor), ub=UBOUND(u_cor), source=u_cor)
  if (present(v_cor))  call v_cor_a%alloc(lb=LBOUND(v_cor), ub=UBOUND(v_cor), source=v_cor)
  if (present(du_cor)) call du_cor_a%alloc(lb=LBOUND(du_cor), ub=UBOUND(du_cor), source=du_cor)
  if (present(dv_cor)) call dv_cor_a%alloc(lb=LBOUND(dv_cor), ub=UBOUND(dv_cor), source=dv_cor)

  call continuity_PPM(bxC, u_a, v_a, hin_a, h_a, uh_a, vh_a, dt, &
                      mask2dT_a, IareaT_a, dy_Cu_a, IdxT_a, dxCu_a, areaT_a, dxT_a, &
                      mask2dCu_a, dx_Cv_a, IdyT_a, dyCv_a, dyT_a, mask2dCv_a, &
                      por_face_areaU_a, por_face_areaV_a, OBC, &
                      G%first_direction, GV%Angstrom_H, GV%H_subroundoff, stencil, &
                      CS%initialized, CS%upwind_1st, CS%monotonic, CS%simple_2nd, &
                      CS%CFL_limit_adjust, CS%aggress_adjust, CS%vol_CFL, CS%better_iter, &
                      CS%use_visc_rem_max, CS%marginal_faces, CS%tol_eta, CS%tol_vel, &
                      uhbt_a, vhbt_a, visc_rem_u_a, visc_rem_v_a, u_cor_a, v_cor_a, &
                      BT_cont, du_cor_a, dv_cor_a)

  call bxC%free()

  call h_a%copy2F(h)   ; call h_a%free()
  call uh_a%copy2F(uh) ; call uh_a%free()
  call vh_a%copy2F(vh) ; call vh_a%free()
  call u_a%free()
  call v_a%free()
  call hin_a%free()

  call mask2dT_a%free()
  call IareaT_a%free()
  call dy_Cu_a%free()
  call IdxT_a%free()
  call dxCu_a%free()
  call areaT_a%free()
  call dxT_a%free()
  call mask2dCu_a%free()
  call dx_Cv_a%free()
  call IdyT_a%free()
  call dyCv_a%free()
  call dyT_a%free()
  call mask2dCv_a%free()
  call por_face_areaU_a%free()
  call por_face_areaV_a%free()

  if (present(uhbt))       call uhbt_a%free()
  if (present(vhbt))       call vhbt_a%free()
  if (present(visc_rem_u)) call visc_rem_u_a%free()
  if (present(visc_rem_v)) call visc_rem_v_a%free()
  if (present(u_cor))  then
    call u_cor_a%copy2F(u_cor) ; call u_cor_a%free()
  endif
  if (present(v_cor))  then
    call v_cor_a%copy2F(v_cor) ; call v_cor_a%free()
  endif
  if (present(du_cor)) then
    call du_cor_a%copy2F(du_cor) ; call du_cor_a%free()
  endif
  if (present(dv_cor)) then
    call dv_cor_a%copy2F(dv_cor) ; call dv_cor_a%free()
  endif

end subroutine continuity

!> Returns the barotropic mass fluxes for use in the barotropic solver. A wrapper around
!! continuity_2d_fluxes_PPM that presents the original raw interface, marshalling into and out
!! of its array containers.
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
  type(RealArray_t) :: por_face_areaU_a, por_face_areaV_a
  type(box_t) :: bxC

  call bxC%safe_alloc(ndims=3)
  call bxC%set(idxS=[G%isc, G%jsc, 1], idxE=[G%iec, G%jec, GV%ke])

  call u_a%alloc(lb=LBOUND(u), ub=UBOUND(u), source=u)
  call v_a%alloc(lb=LBOUND(v), ub=UBOUND(v), source=v)
  call h_a%alloc(lb=LBOUND(h), ub=UBOUND(h), source=h)
  call uhbt_a%alloc(lb=LBOUND(uhbt), ub=UBOUND(uhbt), source=uhbt)
  call vhbt_a%alloc(lb=LBOUND(vhbt), ub=UBOUND(vhbt), source=vhbt)

  call mask2dT_a%alloc(lb=LBOUND(G%mask2dT), ub=UBOUND(G%mask2dT), source=G%mask2dT)
  call dy_Cu_a%alloc(lb=LBOUND(G%dy_Cu), ub=UBOUND(G%dy_Cu), source=G%dy_Cu)
  call IareaT_a%alloc(lb=LBOUND(G%IareaT), ub=UBOUND(G%IareaT), source=G%IareaT)
  call IdxT_a%alloc(lb=LBOUND(G%IdxT), ub=UBOUND(G%IdxT), source=G%IdxT)
  call dx_Cv_a%alloc(lb=LBOUND(G%dx_Cv), ub=UBOUND(G%dx_Cv), source=G%dx_Cv)
  call IdyT_a%alloc(lb=LBOUND(G%IdyT), ub=UBOUND(G%IdyT), source=G%IdyT)
  call por_face_areaU_a%alloc(lb=LBOUND(pbv%por_face_areaU), ub=UBOUND(pbv%por_face_areaU), &
                              source=pbv%por_face_areaU)
  call por_face_areaV_a%alloc(lb=LBOUND(pbv%por_face_areaV), ub=UBOUND(pbv%por_face_areaV), &
                              source=pbv%por_face_areaV)

  call continuity_2d_fluxes_PPM(bxC, u_a, v_a, h_a, uhbt_a, vhbt_a, dt, &
                                mask2dT_a, dy_Cu_a, IareaT_a, IdxT_a, dx_Cv_a, IdyT_a, &
                                por_face_areaU_a, por_face_areaV_a, OBC, &
                                GV%Angstrom_H, CS%upwind_1st, CS%monotonic, CS%simple_2nd, &
                                CS%vol_CFL)

  call bxC%free()

  call uhbt_a%copy2F(uhbt) ; call uhbt_a%free()
  call vhbt_a%copy2F(vhbt) ; call vhbt_a%free()
  call u_a%free()
  call v_a%free()
  call h_a%free()
  call mask2dT_a%free()
  call dy_Cu_a%free()
  call IareaT_a%free()
  call IdxT_a%free()
  call dx_Cv_a%free()
  call IdyT_a%free()
  call por_face_areaU_a%free()
  call por_face_areaV_a%free()

end subroutine continuity_2d_fluxes

!> Returns the 3-D mass fluxes uh and vh. A wrapper around continuity_3d_fluxes_PPM that presents
!! the original raw interface, marshalling into and out of its array containers.
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
  type(RealArray_t) :: mask2dT_a, dy_Cu_a, IareaT_a, IdxT_a, dxCu_a, areaT_a, dxT_a
  type(RealArray_t) :: mask2dCu_a, dx_Cv_a, IdyT_a, dyCv_a, dyT_a, mask2dCv_a
  type(RealArray_t) :: por_face_areaU_a, por_face_areaV_a
  type(box_t) :: bxC

  call bxC%safe_alloc(ndims=3)
  call bxC%set(idxS=[G%isc, G%jsc, 1], idxE=[G%iec, G%jec, GV%ke])

  call u_a%alloc(lb=LBOUND(u), ub=UBOUND(u), source=u)
  call v_a%alloc(lb=LBOUND(v), ub=UBOUND(v), source=v)
  call h_a%alloc(lb=LBOUND(h), ub=UBOUND(h), source=h)
  call uh_a%alloc(lb=LBOUND(uh), ub=UBOUND(uh), source=uh)
  call vh_a%alloc(lb=LBOUND(vh), ub=UBOUND(vh), source=vh)

  call mask2dT_a%alloc(lb=LBOUND(G%mask2dT), ub=UBOUND(G%mask2dT), source=G%mask2dT)
  call dy_Cu_a%alloc(lb=LBOUND(G%dy_Cu), ub=UBOUND(G%dy_Cu), source=G%dy_Cu)
  call IareaT_a%alloc(lb=LBOUND(G%IareaT), ub=UBOUND(G%IareaT), source=G%IareaT)
  call IdxT_a%alloc(lb=LBOUND(G%IdxT), ub=UBOUND(G%IdxT), source=G%IdxT)
  call dxCu_a%alloc(lb=LBOUND(G%dxCu), ub=UBOUND(G%dxCu), source=G%dxCu)
  call areaT_a%alloc(lb=LBOUND(G%areaT), ub=UBOUND(G%areaT), source=G%areaT)
  call dxT_a%alloc(lb=LBOUND(G%dxT), ub=UBOUND(G%dxT), source=G%dxT)
  call mask2dCu_a%alloc(lb=LBOUND(G%mask2dCu), ub=UBOUND(G%mask2dCu), source=G%mask2dCu)
  call dx_Cv_a%alloc(lb=LBOUND(G%dx_Cv), ub=UBOUND(G%dx_Cv), source=G%dx_Cv)
  call IdyT_a%alloc(lb=LBOUND(G%IdyT), ub=UBOUND(G%IdyT), source=G%IdyT)
  call dyCv_a%alloc(lb=LBOUND(G%dyCv), ub=UBOUND(G%dyCv), source=G%dyCv)
  call dyT_a%alloc(lb=LBOUND(G%dyT), ub=UBOUND(G%dyT), source=G%dyT)
  call mask2dCv_a%alloc(lb=LBOUND(G%mask2dCv), ub=UBOUND(G%mask2dCv), source=G%mask2dCv)
  call por_face_areaU_a%alloc(lb=LBOUND(pbv%por_face_areaU), ub=UBOUND(pbv%por_face_areaU), &
                              source=pbv%por_face_areaU)
  call por_face_areaV_a%alloc(lb=LBOUND(pbv%por_face_areaV), ub=UBOUND(pbv%por_face_areaV), &
                              source=pbv%por_face_areaV)

  call continuity_3d_fluxes_PPM(bxC, u_a, v_a, h_a, uh_a, vh_a, dt, &
                                mask2dT_a, dy_Cu_a, IareaT_a, IdxT_a, dxCu_a, areaT_a, dxT_a, &
                                mask2dCu_a, dx_Cv_a, IdyT_a, dyCv_a, dyT_a, mask2dCv_a, &
                                por_face_areaU_a, por_face_areaV_a, OBC, &
                                GV%Angstrom_H, GV%H_subroundoff, &
                                CS%upwind_1st, CS%monotonic, CS%simple_2nd, CS%CFL_limit_adjust, &
                                CS%aggress_adjust, CS%vol_CFL, CS%better_iter, &
                                CS%use_visc_rem_max, CS%marginal_faces, CS%tol_eta, CS%tol_vel)

  call bxC%free()

  call uh_a%copy2F(uh) ; call uh_a%free()
  call vh_a%copy2F(vh) ; call vh_a%free()
  call u_a%free()
  call v_a%free()
  call h_a%free()
  call mask2dT_a%free()
  call dy_Cu_a%free()
  call IareaT_a%free()
  call IdxT_a%free()
  call dxCu_a%free()
  call areaT_a%free()
  call dxT_a%free()
  call mask2dCu_a%free()
  call dx_Cv_a%free()
  call IdyT_a%free()
  call dyCv_a%free()
  call dyT_a%free()
  call mask2dCv_a%free()
  call por_face_areaU_a%free()
  call por_face_areaV_a%free()

end subroutine continuity_3d_fluxes

!> Adjusts the velocities so that their depth-integrated transports match uhbt and vhbt. A
!! wrapper around continuity_adjust_vel_PPM that presents the original raw/container-mixed
!! interface, marshalling into and out of its array containers.
subroutine continuity_adjust_vel(u, v, h, dt, G, GV, US, CS, OBC, pbv, uhbt, vhbt, &
                                 visc_rem_u_a, visc_rem_v_a)
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
  type(RealArray_t),       intent(in) :: visc_rem_u_a !< Both the fraction of the zonal momentum
                                                !! that remains after a time-step of viscosity, and
                                                !! the fraction of a time-step's worth of a barotropic
                                                !! acceleration that a layer experiences after viscosity
                                                !! is applied [nondim].  This goes between 0 (at the
                                                !! bottom) and 1 (far above the bottom).  When this
                                                !! column is under an ice shelf, this also goes to 0
                                                !! at the top due to the no-slip boundary condition there.
  type(RealArray_t),       intent(in) :: visc_rem_v_a !< Both the fraction of the meridional
                                                !! momentum that remains after a time-step of
                                                !! viscosity, and the fraction of a time-step's
                                                !! worth of a barotropic acceleration that a layer
                                                !! experiences after viscosity is applied [nondim].
                                                !! This goes between 0 (at the bottom) and 1 (far
                                                !! above the bottom).  When this column is under an
                                                !! ice shelf, this also goes to 0 at the top due to
                                                !! the no-slip boundary condition there.

  type(RealArray_t) :: u_a, v_a, h_a, uhbt_a, vhbt_a
  type(RealArray_t) :: u_cor_a, v_cor_a
  type(RealArray_t) :: mask2dT_a, dy_Cu_a, IareaT_a, IdxT_a, dxCu_a, areaT_a, dxT_a
  type(RealArray_t) :: mask2dCu_a, dx_Cv_a, IdyT_a, dyCv_a, dyT_a, mask2dCv_a
  type(RealArray_t) :: por_face_areaU_a, por_face_areaV_a
  type(box_t) :: bxC

  call bxC%safe_alloc(ndims=3)
  call bxC%set(idxS=[G%isc, G%jsc, 1], idxE=[G%iec, G%jec, GV%ke])

  call u_a%alloc(lb=LBOUND(u), ub=UBOUND(u), source=u)
  call v_a%alloc(lb=LBOUND(v), ub=UBOUND(v), source=v)
  call h_a%alloc(lb=LBOUND(h), ub=UBOUND(h), source=h)
  call uhbt_a%alloc(lb=LBOUND(uhbt), ub=UBOUND(uhbt), source=uhbt)
  call vhbt_a%alloc(lb=LBOUND(vhbt), ub=UBOUND(vhbt), source=vhbt)
  call u_cor_a%alloc(lb=LBOUND(u), ub=UBOUND(u), source=u)
  call v_cor_a%alloc(lb=LBOUND(v), ub=UBOUND(v), source=v)

  call mask2dT_a%alloc(lb=LBOUND(G%mask2dT), ub=UBOUND(G%mask2dT), source=G%mask2dT)
  call dy_Cu_a%alloc(lb=LBOUND(G%dy_Cu), ub=UBOUND(G%dy_Cu), source=G%dy_Cu)
  call IareaT_a%alloc(lb=LBOUND(G%IareaT), ub=UBOUND(G%IareaT), source=G%IareaT)
  call IdxT_a%alloc(lb=LBOUND(G%IdxT), ub=UBOUND(G%IdxT), source=G%IdxT)
  call dxCu_a%alloc(lb=LBOUND(G%dxCu), ub=UBOUND(G%dxCu), source=G%dxCu)
  call areaT_a%alloc(lb=LBOUND(G%areaT), ub=UBOUND(G%areaT), source=G%areaT)
  call dxT_a%alloc(lb=LBOUND(G%dxT), ub=UBOUND(G%dxT), source=G%dxT)
  call mask2dCu_a%alloc(lb=LBOUND(G%mask2dCu), ub=UBOUND(G%mask2dCu), source=G%mask2dCu)
  call dx_Cv_a%alloc(lb=LBOUND(G%dx_Cv), ub=UBOUND(G%dx_Cv), source=G%dx_Cv)
  call IdyT_a%alloc(lb=LBOUND(G%IdyT), ub=UBOUND(G%IdyT), source=G%IdyT)
  call dyCv_a%alloc(lb=LBOUND(G%dyCv), ub=UBOUND(G%dyCv), source=G%dyCv)
  call dyT_a%alloc(lb=LBOUND(G%dyT), ub=UBOUND(G%dyT), source=G%dyT)
  call mask2dCv_a%alloc(lb=LBOUND(G%mask2dCv), ub=UBOUND(G%mask2dCv), source=G%mask2dCv)
  call por_face_areaU_a%alloc(lb=LBOUND(pbv%por_face_areaU), ub=UBOUND(pbv%por_face_areaU), &
                              source=pbv%por_face_areaU)
  call por_face_areaV_a%alloc(lb=LBOUND(pbv%por_face_areaV), ub=UBOUND(pbv%por_face_areaV), &
                              source=pbv%por_face_areaV)

  call continuity_adjust_vel_PPM(bxC, u_a, v_a, h_a, uhbt_a, vhbt_a, dt, &
                                 mask2dT_a, dy_Cu_a, IareaT_a, IdxT_a, dxCu_a, areaT_a, &
                                 dxT_a, mask2dCu_a, dx_Cv_a, IdyT_a, dyCv_a, dyT_a, &
                                 mask2dCv_a, por_face_areaU_a, por_face_areaV_a, OBC, &
                                 GV%Angstrom_H, GV%H_subroundoff, &
                                 CS%upwind_1st, CS%monotonic, CS%simple_2nd, CS%CFL_limit_adjust, &
                                 CS%aggress_adjust, CS%vol_CFL, CS%better_iter, &
                                 CS%use_visc_rem_max, CS%marginal_faces, CS%tol_eta, CS%tol_vel, &
                                 u_cor_a, v_cor_a, visc_rem_u_a, visc_rem_v_a)

  call bxC%free()

  call u_cor_a%copy2F(u) ; call u_cor_a%free()
  call v_cor_a%copy2F(v) ; call v_cor_a%free()
  call u_a%free()
  call v_a%free()
  call h_a%free()
  call uhbt_a%free()
  call vhbt_a%free()
  call mask2dT_a%free()
  call dy_Cu_a%free()
  call IareaT_a%free()
  call IdxT_a%free()
  call dxCu_a%free()
  call areaT_a%free()
  call dxT_a%free()
  call mask2dCu_a%free()
  call dx_Cv_a%free()
  call IdyT_a%free()
  call dyCv_a%free()
  call dyT_a%free()
  call mask2dCv_a%free()
  call por_face_areaU_a%free()
  call por_face_areaV_a%free()

end subroutine continuity_adjust_vel

end module MOM_continuity
