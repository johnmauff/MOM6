! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0
!!SKILLS: 0.3

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
use MOM_continuity_PPM, only : set_continuity_box
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
public set_continuity_box

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
  type(RealArray_t), pointer :: uhbt_a, vhbt_a, visc_rem_u_a, visc_rem_v_a
  type(RealArray_t), pointer :: u_cor_a, v_cor_a, du_cor_a, dv_cor_a

  call u_a%alloc(lb=LBOUND(u), ub=UBOUND(u), source=u)
  call v_a%alloc(lb=LBOUND(v), ub=UBOUND(v), source=v)
  call hin_a%alloc(lb=LBOUND(hin), ub=UBOUND(hin), source=hin)
  call h_a%alloc(lb=LBOUND(h), ub=UBOUND(h), source=h)
  call uh_a%alloc(lb=LBOUND(uh), ub=UBOUND(uh), source=uh)
  call vh_a%alloc(lb=LBOUND(vh), ub=UBOUND(vh), source=vh)

  nullify(uhbt_a, vhbt_a, visc_rem_u_a, visc_rem_v_a, u_cor_a, v_cor_a, du_cor_a, dv_cor_a)
  if (present(uhbt)) then
    allocate(uhbt_a) ; call uhbt_a%alloc(lb=LBOUND(uhbt), ub=UBOUND(uhbt), source=uhbt)
  endif
  if (present(vhbt)) then
    allocate(vhbt_a) ; call vhbt_a%alloc(lb=LBOUND(vhbt), ub=UBOUND(vhbt), source=vhbt)
  endif
  if (present(visc_rem_u)) then
    allocate(visc_rem_u_a)
    call visc_rem_u_a%alloc(lb=LBOUND(visc_rem_u), ub=UBOUND(visc_rem_u), source=visc_rem_u)
  endif
  if (present(visc_rem_v)) then
    allocate(visc_rem_v_a)
    call visc_rem_v_a%alloc(lb=LBOUND(visc_rem_v), ub=UBOUND(visc_rem_v), source=visc_rem_v)
  endif
  if (present(u_cor)) then
    allocate(u_cor_a) ; call u_cor_a%alloc(lb=LBOUND(u_cor), ub=UBOUND(u_cor), source=u_cor)
  endif
  if (present(v_cor)) then
    allocate(v_cor_a) ; call v_cor_a%alloc(lb=LBOUND(v_cor), ub=UBOUND(v_cor), source=v_cor)
  endif
  if (present(du_cor)) then
    allocate(du_cor_a) ; call du_cor_a%alloc(lb=LBOUND(du_cor), ub=UBOUND(du_cor), source=du_cor)
  endif
  if (present(dv_cor)) then
    allocate(dv_cor_a) ; call dv_cor_a%alloc(lb=LBOUND(dv_cor), ub=UBOUND(dv_cor), source=dv_cor)
  endif

  call continuity_PPM(u_a, v_a, hin_a, h_a, uh_a, vh_a, dt, G, GV, US, CS, OBC, pbv, &
                       uhbt_a=uhbt_a, vhbt_a=vhbt_a, visc_rem_u_a=visc_rem_u_a, &
                       visc_rem_v_a=visc_rem_v_a, u_cor_a=u_cor_a, v_cor_a=v_cor_a, &
                       BT_cont=BT_cont, du_cor_a=du_cor_a, dv_cor_a=dv_cor_a)

  call h_a%copy2F(h)
  call uh_a%copy2F(uh)
  call vh_a%copy2F(vh)

  call u_a%free()
  call v_a%free()
  call hin_a%free()
  call h_a%free()
  call uh_a%free()
  call vh_a%free()

  if (associated(u_cor_a)) then
    call u_cor_a%copy2F(u_cor) ; call u_cor_a%free() ; deallocate(u_cor_a)
  endif
  if (associated(v_cor_a)) then
    call v_cor_a%copy2F(v_cor) ; call v_cor_a%free() ; deallocate(v_cor_a)
  endif
  if (associated(du_cor_a)) then
    call du_cor_a%copy2F(du_cor) ; call du_cor_a%free() ; deallocate(du_cor_a)
  endif
  if (associated(dv_cor_a)) then
    call dv_cor_a%copy2F(dv_cor) ; call dv_cor_a%free() ; deallocate(dv_cor_a)
  endif
  if (associated(uhbt_a)) then ; call uhbt_a%free() ; deallocate(uhbt_a) ; endif
  if (associated(vhbt_a)) then ; call vhbt_a%free() ; deallocate(vhbt_a) ; endif
  if (associated(visc_rem_u_a)) then ; call visc_rem_u_a%free() ; deallocate(visc_rem_u_a) ; endif
  if (associated(visc_rem_v_a)) then ; call visc_rem_v_a%free() ; deallocate(visc_rem_v_a) ; endif

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

  call u_a%alloc(lb=LBOUND(u), ub=UBOUND(u), source=u)
  call v_a%alloc(lb=LBOUND(v), ub=UBOUND(v), source=v)
  call h_a%alloc(lb=LBOUND(h), ub=UBOUND(h), source=h)
  call uh_a%alloc(lb=LBOUND(uh), ub=UBOUND(uh), source=uh)
  call vh_a%alloc(lb=LBOUND(vh), ub=UBOUND(vh), source=vh)

  call continuity_PPM_3d_fluxes(u_a, v_a, h_a, uh_a, vh_a, dt, G, GV, US, CS, OBC, pbv)

  call uh_a%copy2F(uh)
  call vh_a%copy2F(vh)

  call u_a%free()
  call v_a%free()
  call h_a%free()
  call uh_a%free()
  call vh_a%free()

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

  call u_a%alloc(lb=LBOUND(u), ub=UBOUND(u), source=u)
  call v_a%alloc(lb=LBOUND(v), ub=UBOUND(v), source=v)
  call h_a%alloc(lb=LBOUND(h), ub=UBOUND(h), source=h)
  call uhbt_a%alloc(lb=LBOUND(uhbt), ub=UBOUND(uhbt), source=uhbt)
  call vhbt_a%alloc(lb=LBOUND(vhbt), ub=UBOUND(vhbt), source=vhbt)

  call continuity_PPM_2d_fluxes(u_a, v_a, h_a, uhbt_a, vhbt_a, dt, G, GV, US, CS, OBC, pbv)

  call uhbt_a%copy2F(uhbt)
  call vhbt_a%copy2F(vhbt)

  call u_a%free()
  call v_a%free()
  call h_a%free()
  call uhbt_a%free()
  call vhbt_a%free()

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
  type(RealArray_t), pointer :: visc_rem_u_a, visc_rem_v_a

  call u_a%alloc(lb=LBOUND(u), ub=UBOUND(u), source=u)
  call v_a%alloc(lb=LBOUND(v), ub=UBOUND(v), source=v)
  call h_a%alloc(lb=LBOUND(h), ub=UBOUND(h), source=h)
  call uhbt_a%alloc(lb=LBOUND(uhbt), ub=UBOUND(uhbt), source=uhbt)
  call vhbt_a%alloc(lb=LBOUND(vhbt), ub=UBOUND(vhbt), source=vhbt)

  nullify(visc_rem_u_a, visc_rem_v_a)
  if (present(visc_rem_u)) then
    allocate(visc_rem_u_a)
    call visc_rem_u_a%alloc(lb=LBOUND(visc_rem_u), ub=UBOUND(visc_rem_u), source=visc_rem_u)
  endif
  if (present(visc_rem_v)) then
    allocate(visc_rem_v_a)
    call visc_rem_v_a%alloc(lb=LBOUND(visc_rem_v), ub=UBOUND(visc_rem_v), source=visc_rem_v)
  endif

  call continuity_PPM_adjust_vel(u_a, v_a, h_a, dt, G, GV, US, CS, OBC, pbv, uhbt_a, vhbt_a, &
                                 visc_rem_u_a=visc_rem_u_a, visc_rem_v_a=visc_rem_v_a)

  call u_a%copy2F(u)
  call v_a%copy2F(v)

  call u_a%free()
  call v_a%free()
  call h_a%free()
  call uhbt_a%free()
  call vhbt_a%free()
  if (associated(visc_rem_u_a)) then ; call visc_rem_u_a%free() ; deallocate(visc_rem_u_a) ; endif
  if (associated(visc_rem_v_a)) then ; call visc_rem_v_a%free() ; deallocate(visc_rem_v_a) ; endif

end subroutine continuity_adjust_vel

end module MOM_continuity
