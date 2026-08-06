! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

#include "do_concurrent_compat.h"

!> Solve the layer continuity equation using the PPM method for layer fluxes.
module MOM_continuity_PPM

use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_diag_mediator, only : time_type, diag_ctrl
use MOM_error_handler, only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_open_boundary, only : ocean_OBC_type, OBC_segment_type, OBC_NONE
use MOM_open_boundary, only : OBC_DIRECTION_E, OBC_DIRECTION_W, OBC_DIRECTION_N, OBC_DIRECTION_S
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : BT_cont_type, porous_barrier_type
use MOM_verticalGrid, only : verticalGrid_type

use array_mod, only : RealArray_t, RealArray_c, LogicalArray_t
use box_mod, only : Box_t, Box_c
use iso_c_binding, only : c_double, c_int, c_ptr, c_loc, c_bool, c_null_char, c_null_ptr
use posix, only : mkdir_posix

use turbotmp_helperF, only : getenv_mode, io_recorder, already_recorded, mark_recorded
use turbotmp_helperF, only : TIMH_runAMREX, TIMH_capture, TIMH_runFORTRAN

implicit none ; private

  !----------------------------------------
  ! C interface (bridge to C++)
  !----------------------------------------
  interface
    !> Bridge for the PPM_limit_pos subroutine
    subroutine turbotmp_ppm_limit_pos_bridge(bx, h_in, h_L, h_R, h_min) bind(C)
       use iso_c_binding
       use array_mod, only : RealArray_c
       use box_mod,   only : Box_c
       implicit none
       type(Box_C), intent(in)          :: bx    !< Index space over which to iterate
       type(RealArray_C), intent(in)    :: h_in  !< Layer thickness [H ~> m or kg m-2].
       type(RealArray_C), intent(inout) :: h_L   !< Left thickness in the reconstruction
                                                 !! [H ~> m or kg m-2].
       type(RealArray_C), intent(inout) :: h_R   !< Right thickness in the reconstruction
                                                 !! [H ~> m or kg m-2].
       real(c_double), intent(in),value :: h_min !< The minimum thickness that can be obtained
                                                 !! by a concave parabolic fit [H ~> m or kg m-2]
    end subroutine turbotmp_ppm_limit_pos_bridge
  end interface

  interface
    !> Bridge for the PPM_limit_cw84 subroutine
    subroutine turbotmp_ppm_limit_cw84_bridge(bx, h_in, h_L, h_R) bind(C)
      use iso_c_binding
      use array_mod, only : RealArray_c
      use box_mod,   only : Box_c
      implicit none

      type(Box_C), intent(in)           :: bx   !< Index space over which to iterate
      type(RealArray_C),  intent(in)    :: h_in !< Layer thickness [H ~> m or kg m-2].
      type(RealArray_C),  intent(inout) :: h_L  !< Left thickness in the reconstruction
                                                !! [H ~> m or kg m-2].
      type(RealArray_C), intent(inout)  :: h_R  !< Right thickness in the reconstruction
                                                !! [H ~> m or kg m-2].
    end subroutine turbotmp_ppm_limit_cw84_bridge
  end interface

  interface
    !> Bridge for the PPM_reconstruction_y subroutine
    subroutine turbotmp_ppm_reconstruction_y_bridge(bx, h_in, h_S, h_N, mask2dT, &
                                           h_min, monotonic, simple_2nd, obc) bind(C)
      use iso_c_binding, only : c_double, c_bool, c_ptr
      use array_mod, only : RealArray_c
      use box_mod,   only : Box_c
      implicit none

      type(Box_C), intent(in)            :: bx         !< Index space over which to iterate
      type(RealArray_C), intent(in)      :: h_in       !< Layer thickness
      type(RealArray_C), intent(inout)   :: h_S        !< South edge thickness
      type(RealArray_C), intent(inout)   :: h_N        !< North edge thickness
      type(RealArray_C), intent(in)      :: mask2dT    !< Mask (0 land, 1 ocean)
      real(c_double),  intent(in), value :: h_min      !< Minimum thickness
      logical(c_bool), intent(in), value :: monotonic  !< Use CW84 limiter
      logical(c_bool), intent(in), value :: simple_2nd !< Use 2nd order scheme
      type(c_ptr),     intent(in), value :: obc        !< Pointer to OBC structure
    end subroutine turbotmp_ppm_reconstruction_y_bridge
  end interface

  interface
    !> Bridge for the PPM_reconstruction_x subroutine
    subroutine turbotmp_ppm_reconstruction_x_bridge(bx, h_in, h_W, h_E, mask2dT, &
                                           h_min, monotonic, simple_2nd, obc) bind(C)
      use iso_c_binding, only : c_double, c_bool, c_ptr
      use array_mod, only : RealArray_c
      use box_mod,   only : Box_c
      implicit none

      type(Box_C), intent(in)            :: bx         !< Index space over which to iterate
      type(RealArray_C), intent(in)      :: h_in       !< Layer thickness
      type(RealArray_C), intent(inout)   :: h_W        !< West edge thickness
      type(RealArray_C), intent(inout)   :: h_E        !< East edge thickness
      type(RealArray_C), intent(in)      :: mask2dT    !< Mask (0 land, 1 ocean)
      real(c_double),  intent(in), value :: h_min      !< Minimum thickness
      logical(c_bool), intent(in), value :: monotonic  !< Use CW84 limiter
      logical(c_bool), intent(in), value :: simple_2nd !< Use 2nd order scheme
      type(c_ptr),     intent(in), value :: obc        !< Pointer to OBC structure
    end subroutine turbotmp_ppm_reconstruction_x_bridge
  end interface

  interface
    !> Bridge for the zonal_edge_thickness subroutine
    subroutine turbotmp_zonal_edge_thickness_bridge(bx, h_in, h_W, h_E, mask2dT, &
                                           h_min, upwind_1st, monotonic, simple_2nd, obc) bind(C)
      use iso_c_binding, only : c_double, c_bool, c_ptr
      use array_mod, only : RealArray_c
      use box_mod,   only : Box_c
      implicit none

      type(Box_C),       intent(in)          :: bx          !< Index space over which to iterate
      type(RealArray_C), intent(in)          :: h_in        !< Layer thickness [H ~> m or kg m-2]
      type(RealArray_C), intent(inout)       :: h_W         !< Western edge thickness [H ~> m or kg m-2]
      type(RealArray_C), intent(inout)       :: h_E         !< Eastern edge thickness [H ~> m or kg m-2]
      type(RealArray_C), intent(in)          :: mask2dT     !< Cell land/ocean mask [nondim]
      real(c_double),    intent(in), value   :: h_min       !< Minimum layer thickness [H ~> m or kg m-2]
      logical(c_bool),   intent(in), value   :: upwind_1st  !< Use 1st-order upwind reconstruction
      logical(c_bool),   intent(in), value   :: monotonic   !< Use CW84 monotonic limiter
      logical(c_bool),   intent(in), value   :: simple_2nd  !< Use simple 2nd-order scheme
      type(c_ptr),       intent(in), value   :: obc         !< Pointer to OBC structure
    end subroutine turbotmp_zonal_edge_thickness_bridge
  end interface

  interface
    !> Bridge for the meridional_edge_thickness subroutine
    subroutine turbotmp_meridional_edge_thickness_bridge(bx, h_in, h_S, h_N, mask2dT, &
                                           h_min, upwind_1st, monotonic, simple_2nd, obc) bind(C)
      use iso_c_binding, only : c_double, c_bool, c_ptr
      use array_mod, only : RealArray_c
      use box_mod,   only : Box_c
      implicit none

      type(Box_C),       intent(in)          :: bx          !< Index space over which to iterate
      type(RealArray_C), intent(in)          :: h_in        !< Layer thickness [H ~> m or kg m-2]
      type(RealArray_C), intent(inout)       :: h_S         !< Southern edge thickness [H ~> m or kg m-2]
      type(RealArray_C), intent(inout)       :: h_N         !< Northern edge thickness [H ~> m or kg m-2]
      type(RealArray_C), intent(in)          :: mask2dT     !< Cell land/ocean mask [nondim]
      real(c_double),    intent(in), value   :: h_min       !< Minimum layer thickness [H ~> m or kg m-2]
      logical(c_bool),   intent(in), value   :: upwind_1st  !< Use 1st-order upwind reconstruction
      logical(c_bool),   intent(in), value   :: monotonic   !< Use CW84 monotonic limiter
      logical(c_bool),   intent(in), value   :: simple_2nd  !< Use simple 2nd-order scheme
      type(c_ptr),       intent(in), value   :: obc         !< Pointer to OBC structure
    end subroutine turbotmp_meridional_edge_thickness_bridge
  end interface

#include <MOM_memory.h>

public continuity_PPM, continuity_PPM_init, continuity_PPM_stencil
public continuity_fluxes, continuity_adjust_vel
public zonal_mass_flux, meridional_mass_flux
public zonal_edge_thickness, meridional_edge_thickness
public continuity_zonal_convergence, continuity_meridional_convergence
public zonal_flux_thickness, meridional_flux_thickness
public zonal_BT_mass_flux, meridional_BT_mass_flux
public set_continuity_loop_bounds
public set_continuity_box

!>@{ CPU time clock IDs
integer :: id_clock_reconstruct, id_clock_update, id_clock_correct
!>@}

!> Control structure for mom_continuity_ppm
type, public :: continuity_PPM_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  type(diag_ctrl), pointer :: diag !< Diagnostics control structure.
  logical :: upwind_1st      !< If true, use a first-order upwind scheme.
  logical :: monotonic       !< If true, use the Colella & Woodward monotonic
                             !! limiter; otherwise use a simple positive
                             !! definite limiter.
  logical :: simple_2nd      !< If true, use a simple second order (arithmetic
                             !! mean) interpolation of the edge values instead
                             !! of the higher order interpolation.
  real :: tol_eta            !< The tolerance for free-surface height
                             !! discrepancies between the barotropic solution and
                             !! the sum of the layer thicknesses [H ~> m or kg m-2].
  real :: tol_vel            !< The tolerance for barotropic velocity
                             !! discrepancies between the barotropic solution and
                             !! the sum of the layer thicknesses [L T-1 ~> m s-1].
  real :: CFL_limit_adjust   !< The maximum CFL of the adjusted velocities [nondim]
  real :: h_marg_min         !< Negligible floor on h_marg, the marginal thickness
                             !! used to calculate the partial derivative of transports
                             !! with velocities [H ~> m or kg m-2]
  logical :: aggress_adjust  !< If true, allow the adjusted velocities to have a
                             !! relative CFL change up to 0.5.  False by default.
  logical :: vol_CFL         !< If true, use the ratio of the open face lengths
                             !! to the tracer cell areas when estimating CFL
                             !! numbers.  Without aggress_adjust, the default is
                             !! false; it is always true with.
  logical :: better_iter     !< If true, stop corrective iterations using a
                             !! velocity-based criterion and only stop if the
                             !! iteration is better than all predecessors.
  logical :: use_visc_rem_max !< If true, use more appropriate limiting bounds
                             !! for corrections in strongly viscous columns.
  logical :: marginal_faces  !< If true, use the marginal face areas from the
                             !! continuity solver for use as the weights in the
                             !! barotropic solver.  Otherwise use the transport
                             !! averaged areas.
end type continuity_PPM_CS

!> A container for loop bounds
type, public :: cont_loop_bounds_type ; private
  !>@{ Loop bounds
  integer :: ish, ieh, jsh, jeh
  !>@}
end type cont_loop_bounds_type

!> Finds the thickness fluxes from the continuity solver or their vertical sum without
!! actually updating the layer thicknesses.
interface continuity_fluxes
  module procedure continuity_3d_fluxes, continuity_2d_fluxes
end interface continuity_fluxes

contains

!> Time steps the layer thicknesses, using a monotonically limit, directionally split PPM scheme,
!! based on Lin (1994).
subroutine continuity_PPM(u, v, hin, h, uh, vh, dt, G, GV, US, CS, OBC, pbv, uhbt_a, vhbt_a, &
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
  type(continuity_PPM_CS), intent(in)    :: CS  !< Module's control structure.
  type(ocean_OBC_type),    pointer       :: OBC !< Open boundaries control structure.
  type(porous_barrier_type), intent(in)  :: pbv !< pointers to porous barrier fractional cell metrics
  type(RealArray_t), optional, intent(in)    :: uhbt_a !< The summed volume flux through zonal faces
                                                       !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(RealArray_t), optional, intent(in)    :: vhbt_a !< The summed volume flux through meridional
                                                       !! faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(RealArray_t), optional, intent(in)    :: visc_rem_u_a
                             !< The fraction of zonal momentum originally
                             !! in a layer that remains after a time-step of viscosity, and the
                             !! fraction of a time-step's worth of a barotropic acceleration that
                             !! a layer experiences after viscosity is applied [nondim].
                             !! Visc_rem_u is between 0 (at the bottom) and 1 (far above the bottom).
  type(RealArray_t), optional, intent(in)    :: visc_rem_v_a
                             !< The fraction of meridional momentum originally
                             !! in a layer that remains after a time-step of viscosity, and the
                             !! fraction of a time-step's worth of a barotropic acceleration that
                             !! a layer experiences after viscosity is applied [nondim].
                             !! Visc_rem_v is between 0 (at the bottom) and 1 (far above the bottom).
  type(RealArray_t), optional, intent(inout) :: u_cor_a
                             !< The zonal velocities that give uhbt as the depth-integrated
                             !! transport [L T-1 ~> m s-1].
  type(RealArray_t), optional, intent(inout) :: v_cor_a
                             !< The meridional velocities that give vhbt as the depth-integrated
                             !! transport [L T-1 ~> m s-1].
  type(BT_cont_type), optional, pointer  :: BT_cont !< A structure with elements that describe
                             !!  the effective open face areas as a function of barotropic flow.
  type(RealArray_t), optional, intent(inout) :: du_cor_a !< The zonal velocity increments from u
                                                         !! that give uhbt as the depth-integrated
                                                         !! transports [L T-1 ~> m s-1].
  type(RealArray_t), optional, intent(inout) :: dv_cor_a !< The meridional velocity increments from
                                                         !! v that give vhbt as the depth-integrated
                                                         !! transports [L T-1 ~> m s-1].

  ! Local variables
  real :: h_W(SZI_(G),SZJ_(G),SZK_(GV)) ! West edge thicknesses in the zonal PPM reconstruction [H ~> m or kg m-2]
  real :: h_E(SZI_(G),SZJ_(G),SZK_(GV)) ! East edge thicknesses in the zonal PPM reconstruction [H ~> m or kg m-2]
  real :: h_S(SZI_(G),SZJ_(G),SZK_(GV)) ! South edge thicknesses in the meridional PPM reconstruction [H ~> m or kg m-2]
  real :: h_N(SZI_(G),SZJ_(G),SZK_(GV)) ! North edge thicknesses in the meridional PPM reconstruction [H ~> m or kg m-2]
  real :: h_min  ! The minimum layer thickness [H ~> m or kg m-2].  h_min could be 0.
  type(box_t) :: bxC                ! An iteration box
  logical :: x_first
  type(RealArray_t) :: h_in_a, h_W_a, h_E_a, mask2dT_a
  type(RealArray_t) :: h_S_a, h_N_a
  type(RealArray_t) :: h_a, uh_a, IareaT_a, hin_a
  type(RealArray_t) :: vh_a
  type(RealArray_t) :: u_a, por_face_areaU_a
  type(RealArray_t) :: v_a, por_face_areaV_a
  type(RealArray_t) :: dy_Cu_a, IdxT_a, dxCu_a, areaT_a, dxT_a, mask2dCu_a
  type(RealArray_t) :: dx_Cv_a, IdyT_a, dyCv_a, dyT_a, mask2dCv_a
  type(RealArray_t) :: no_hin_a ! Never allocated; unassociated data signals hin_a absent.
  real :: edge_h_min ! Minimum layer thickness (2*Angstrom_H) [H ~> m or kg m-2]

  h_min = GV%Angstrom_H

  if (.not.CS%initialized) call MOM_error(FATAL, &
         "MOM_continuity_PPM: Module must be initialized before it is used.")

  x_first = (MOD(G%first_direction,2) == 0)

  if (present(visc_rem_u_a) .neqv. present(visc_rem_v_a)) call MOM_error(FATAL, &
      "MOM_continuity_PPM: Either both visc_rem_u_a and visc_rem_v_a or neither "// &
      "one must be present in call to continuity_PPM.")

  edge_h_min = 2.0 * GV%Angstrom_H

  !$omp target enter data map(alloc: h_W, h_E, h_S, h_N)

  call h_W_a%alloc(lb=LBOUND(h_W), ub=UBOUND(h_W), source=h_W)
  call h_E_a%alloc(lb=LBOUND(h_E), ub=UBOUND(h_E), source=h_E)
  call h_S_a%alloc(lb=LBOUND(h_S), ub=UBOUND(h_S), source=h_S)
  call h_N_a%alloc(lb=LBOUND(h_N), ub=UBOUND(h_N), source=h_N)
  call mask2dT_a%alloc(lb=LBOUND(G%mask2dT), ub=UBOUND(G%mask2dT), source=G%mask2dT)
  call IareaT_a%alloc(lb=LBOUND(G%IareaT), ub=UBOUND(G%IareaT), source=G%IareaT)
  call hin_a%alloc(lb=LBOUND(hin), ub=UBOUND(hin), source=hin)
  call h_in_a%alloc(lb=LBOUND(hin), ub=UBOUND(hin), source=hin)
  call h_a%alloc(lb=LBOUND(h), ub=UBOUND(h), source=h)
  call u_a%alloc(lb=LBOUND(u), ub=UBOUND(u), source=u)
  call uh_a%alloc(lb=LBOUND(uh), ub=UBOUND(uh), source=uh)
  call por_face_areaU_a%alloc(lb=LBOUND(pbv%por_face_areaU), ub=UBOUND(pbv%por_face_areaU), &
                              source=pbv%por_face_areaU)
  call v_a%alloc(lb=LBOUND(v), ub=UBOUND(v), source=v)
  call vh_a%alloc(lb=LBOUND(vh), ub=UBOUND(vh), source=vh)
  call por_face_areaV_a%alloc(lb=LBOUND(pbv%por_face_areaV), ub=UBOUND(pbv%por_face_areaV), &
                              source=pbv%por_face_areaV)
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

  if (x_first) then
    !  First advect zonally, with loop bounds that accomodate the subsequent meridional advection.
    !LB  = set_continuity_loop_bounds(G, CS, i_stencil=.false., j_stencil=.true.)
    bxC = set_continuity_box(G,GV, CS, i_stencil=.false., j_stencil=.true.)
    call zonal_edge_thickness(bxC, h_in_a, h_W_a, h_E_a, mask2dT_a, &
                              edge_h_min, CS%upwind_1st, CS%monotonic, CS%simple_2nd, OBC)
    call zonal_mass_flux(bxC, u_a, h_in_a, h_W_a, h_E_a, uh_a, dt, OBC, &
                         por_face_areaU_a, uhbt_a=uhbt_a, visc_rem_u_a=visc_rem_u_a, &
                         u_cor_a=u_cor_a, BT_cont=BT_cont, du_cor_a=du_cor_a, &
                         dy_Cu_a=dy_Cu_a, IareaT_a=IareaT_a, IdxT_a=IdxT_a, dxCu_a=dxCu_a, &
                         areaT_a=areaT_a, dxT_a=dxT_a, mask2dCu_a=mask2dCu_a, &
                         H_subroundoff=GV%H_subroundoff, &
                         CFL_limit_adjust=CS%CFL_limit_adjust, aggress_adjust=CS%aggress_adjust, &
                         use_visc_rem_max=CS%use_visc_rem_max, vol_CFL=CS%vol_CFL, &
                         tol_vel=CS%tol_vel, tol_eta=CS%tol_eta, better_iter=CS%better_iter, &
                         marginal_faces=CS%marginal_faces)
    call continuity_zonal_convergence(bxC, h_a, uh_a, dt, IareaT_a, hin_a=hin_a)

    ! update host h from continuity_zonal_convergence

    !  Now advect meridionally, using the updated thicknesses to determine the fluxes.
    !LB  = set_continuity_loop_bounds(G, CS, i_stencil=.false., j_stencil=.false.)
    bxC = set_continuity_box(G, GV, CS, i_stencil=.false., j_stencil=.false.)
    call meridional_edge_thickness(bxC, h_a, h_S_a, h_N_a, mask2dT_a, &
                                   edge_h_min, CS%upwind_1st, CS%monotonic, CS%simple_2nd, OBC)
    call meridional_mass_flux(bxC, v_a, h_a, h_S_a, h_N_a, vh_a, dt, OBC, &
                              por_face_areaV_a, vhbt_a=vhbt_a, visc_rem_v_a=visc_rem_v_a, &
                              v_cor_a=v_cor_a, BT_cont=BT_cont, dv_cor_a=dv_cor_a, &
                              dx_Cv_a=dx_Cv_a, IareaT_a=IareaT_a, IdyT_a=IdyT_a, dyCv_a=dyCv_a, &
                              areaT_a=areaT_a, dyT_a=dyT_a, mask2dCv_a=mask2dCv_a, &
                              H_subroundoff=GV%H_subroundoff, &
                              CFL_limit_adjust=CS%CFL_limit_adjust, &
                              aggress_adjust=CS%aggress_adjust, &
                              use_visc_rem_max=CS%use_visc_rem_max, vol_CFL=CS%vol_CFL, &
                              tol_vel=CS%tol_vel, tol_eta=CS%tol_eta, better_iter=CS%better_iter, &
                              marginal_faces=CS%marginal_faces)
    call continuity_meridional_convergence(bxC, h_a, vh_a, dt, IareaT_a, hin_a=no_hin_a, hmin=h_min)

  else  ! .not. x_first
    !  First advect meridionally, with loop bounds that accomodate the subsequent zonal advection.
    !LB  = set_continuity_loop_bounds(G, CS, i_stencil=.true., j_stencil=.false.)
    bxC = set_continuity_box(G, GV, CS, i_stencil=.true., j_stencil=.false.)
    call meridional_edge_thickness(bxC, h_in_a, h_S_a, h_N_a, mask2dT_a, &
                                   edge_h_min, CS%upwind_1st, CS%monotonic, CS%simple_2nd, OBC)
    call meridional_mass_flux(bxC, v_a, h_in_a, h_S_a, h_N_a, vh_a, dt, OBC, &
                              por_face_areaV_a, vhbt_a=vhbt_a, visc_rem_v_a=visc_rem_v_a, &
                              v_cor_a=v_cor_a, BT_cont=BT_cont, dv_cor_a=dv_cor_a, &
                              dx_Cv_a=dx_Cv_a, IareaT_a=IareaT_a, IdyT_a=IdyT_a, dyCv_a=dyCv_a, &
                              areaT_a=areaT_a, dyT_a=dyT_a, mask2dCv_a=mask2dCv_a, &
                              H_subroundoff=GV%H_subroundoff, &
                              CFL_limit_adjust=CS%CFL_limit_adjust, &
                              aggress_adjust=CS%aggress_adjust, &
                              use_visc_rem_max=CS%use_visc_rem_max, vol_CFL=CS%vol_CFL, &
                              tol_vel=CS%tol_vel, tol_eta=CS%tol_eta, better_iter=CS%better_iter, &
                              marginal_faces=CS%marginal_faces)
    call continuity_meridional_convergence(bxC, h_a, vh_a, dt, IareaT_a, hin_a=hin_a)

    !  Now advect zonally, using the updated thicknesses to determine the fluxes.
    !LB  = set_continuity_loop_bounds(G, CS, i_stencil=.false., j_stencil=.false.)
    bxC = set_continuity_box(G, GV, CS, i_stencil=.false., j_stencil=.false.)
    call zonal_edge_thickness(bxC, h_a, h_W_a, h_E_a, mask2dT_a, &
                              edge_h_min, CS%upwind_1st, CS%monotonic, CS%simple_2nd, OBC)
    call zonal_mass_flux(bxC, u_a, h_a, h_W_a, h_E_a, uh_a, dt, OBC, &
                         por_face_areaU_a, uhbt_a=uhbt_a, visc_rem_u_a=visc_rem_u_a, &
                         u_cor_a=u_cor_a, BT_cont=BT_cont, du_cor_a=du_cor_a, &
                         dy_Cu_a=dy_Cu_a, IareaT_a=IareaT_a, IdxT_a=IdxT_a, dxCu_a=dxCu_a, &
                         areaT_a=areaT_a, dxT_a=dxT_a, mask2dCu_a=mask2dCu_a, &
                         H_subroundoff=GV%H_subroundoff, &
                         CFL_limit_adjust=CS%CFL_limit_adjust, aggress_adjust=CS%aggress_adjust, &
                         use_visc_rem_max=CS%use_visc_rem_max, vol_CFL=CS%vol_CFL, &
                         tol_vel=CS%tol_vel, tol_eta=CS%tol_eta, better_iter=CS%better_iter, &
                         marginal_faces=CS%marginal_faces)
    call continuity_zonal_convergence(bxC, h_a, uh_a, dt, IareaT_a, hin_a=no_hin_a, hmin=h_min)
  endif

  call h_W_a%copy2F(h_W)
  call h_E_a%copy2F(h_E)
  call h_S_a%copy2F(h_S)
  call h_N_a%copy2F(h_N)
  call uh_a%copy2F(uh)
  call vh_a%copy2F(vh)
  call h_a%copy2F(h)

  call h_W_a%free()
  call h_E_a%free()
  call h_S_a%free()
  call h_N_a%free()
  call mask2dT_a%free()
  call IareaT_a%free()
  call hin_a%free()
  call h_in_a%free()
  call h_a%free()
  call u_a%free()
  call uh_a%free()
  call por_face_areaU_a%free()
  call v_a%free()
  call vh_a%free()
  call por_face_areaV_a%free()
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

  ! Free the continuity solver iteration box
  call bxC%free()
  !$omp target exit data map(delete: h_W, h_E, h_S, h_N)

end subroutine continuity_PPM

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
  type(continuity_PPM_CS), intent(in)    :: CS  !< Control structure for mom_continuity.
  type(ocean_OBC_type),    pointer       :: OBC !< Open boundaries control structure.
  type(porous_barrier_type), intent(in)  :: pbv !< porous barrier fractional cell metrics

  ! Local variables
  real :: h_W(SZI_(G),SZJ_(G),SZK_(GV)) ! West edge thicknesses in the zonal PPM reconstruction [H ~> m or kg m-2]
  real :: h_E(SZI_(G),SZJ_(G),SZK_(GV)) ! East edge thicknesses in the zonal PPM reconstruction [H ~> m or kg m-2]
  real :: h_S(SZI_(G),SZJ_(G),SZK_(GV)) ! South edge thicknesses in the meridional PPM reconstruction [H ~> m or kg m-2]
  real :: h_N(SZI_(G),SZJ_(G),SZK_(GV)) ! North edge thicknesses in the meridional PPM reconstruction [H ~> m or kg m-2]
  type (box_t) :: bxC                   ! Iteration box for the continuity solver
  type(RealArray_t) :: h_in_a, h_W_a, h_E_a, mask2dT_a
  type(RealArray_t) :: h_S_a, h_N_a
  type(RealArray_t) :: u_a, uh_a, por_face_areaU_a
  type(RealArray_t) :: v_a, vh_a, por_face_areaV_a
  type(RealArray_t) :: dy_Cu_a, IareaT_a, IdxT_a, dxCu_a, areaT_a, dxT_a, mask2dCu_a
  type(RealArray_t) :: dx_Cv_a, IdyT_a, dyCv_a, dyT_a, mask2dCv_a
  real :: h_min                         ! Minimum layer thickness (2*Angstrom_H) [H ~> m or kg m-2]

  ! Construct the iteration box
  bxC = set_continuity_box(G,GV, CS)

  h_min = 2.0 * GV%Angstrom_H

  call h_in_a%alloc(lb=LBOUND(h), ub=UBOUND(h), source=h)
  call h_W_a%alloc(lb=LBOUND(h_W), ub=UBOUND(h_W), source=h_W)
  call h_E_a%alloc(lb=LBOUND(h_E), ub=UBOUND(h_E), source=h_E)
  call h_S_a%alloc(lb=LBOUND(h_S), ub=UBOUND(h_S), source=h_S)
  call h_N_a%alloc(lb=LBOUND(h_N), ub=UBOUND(h_N), source=h_N)
  call mask2dT_a%alloc(lb=LBOUND(G%mask2dT), ub=UBOUND(G%mask2dT), source=G%mask2dT)
  call u_a%alloc(lb=LBOUND(u), ub=UBOUND(u), source=u)
  call uh_a%alloc(lb=LBOUND(uh), ub=UBOUND(uh), source=uh)
  call por_face_areaU_a%alloc(lb=LBOUND(pbv%por_face_areaU), ub=UBOUND(pbv%por_face_areaU), &
                              source=pbv%por_face_areaU)
  call v_a%alloc(lb=LBOUND(v), ub=UBOUND(v), source=v)
  call vh_a%alloc(lb=LBOUND(vh), ub=UBOUND(vh), source=vh)
  call por_face_areaV_a%alloc(lb=LBOUND(pbv%por_face_areaV), ub=UBOUND(pbv%por_face_areaV), &
                              source=pbv%por_face_areaV)
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

  call zonal_edge_thickness(bxC, h_in_a, h_W_a, h_E_a, mask2dT_a, &
                            h_min, CS%upwind_1st, CS%monotonic, CS%simple_2nd, OBC)
  call zonal_mass_flux(bxC, u_a, h_in_a, h_W_a, h_E_a, uh_a, dt, OBC, &
                       por_face_areaU_a, dy_Cu_a=dy_Cu_a, IareaT_a=IareaT_a, IdxT_a=IdxT_a, &
                       dxCu_a=dxCu_a, areaT_a=areaT_a, dxT_a=dxT_a, mask2dCu_a=mask2dCu_a, &
                       H_subroundoff=GV%H_subroundoff, &
                       CFL_limit_adjust=CS%CFL_limit_adjust, aggress_adjust=CS%aggress_adjust, &
                       use_visc_rem_max=CS%use_visc_rem_max, vol_CFL=CS%vol_CFL, &
                       tol_vel=CS%tol_vel, tol_eta=CS%tol_eta, better_iter=CS%better_iter, &
                       marginal_faces=CS%marginal_faces)
  call meridional_edge_thickness(bxC, h_in_a, h_S_a, h_N_a, mask2dT_a, &
                                 h_min, CS%upwind_1st, CS%monotonic, CS%simple_2nd, OBC)
  call meridional_mass_flux(bxC, v_a, h_in_a, h_S_a, h_N_a, vh_a, dt, OBC, &
                            por_face_areaV_a, dx_Cv_a=dx_Cv_a, IareaT_a=IareaT_a, &
                            IdyT_a=IdyT_a, dyCv_a=dyCv_a, areaT_a=areaT_a, dyT_a=dyT_a, &
                            mask2dCv_a=mask2dCv_a, H_subroundoff=GV%H_subroundoff, &
                            CFL_limit_adjust=CS%CFL_limit_adjust, &
                            aggress_adjust=CS%aggress_adjust, &
                            use_visc_rem_max=CS%use_visc_rem_max, vol_CFL=CS%vol_CFL, &
                            tol_vel=CS%tol_vel, tol_eta=CS%tol_eta, better_iter=CS%better_iter, &
                            marginal_faces=CS%marginal_faces)

  call h_W_a%copy2F(h_W)
  call h_E_a%copy2F(h_E)
  call h_S_a%copy2F(h_S)
  call h_N_a%copy2F(h_N)
  call uh_a%copy2F(uh)
  call vh_a%copy2F(vh)

  call h_in_a%free()
  call h_W_a%free()
  call h_E_a%free()
  call h_S_a%free()
  call h_N_a%free()
  call mask2dT_a%free()
  call u_a%free()
  call uh_a%free()
  call por_face_areaU_a%free()
  call v_a%free()
  call vh_a%free()
  call por_face_areaV_a%free()
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

  ! Free the continuity solver iteration box
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
  type(continuity_PPM_CS), intent(in)    :: CS  !< Control structure for mom_continuity.
  type(ocean_OBC_type),    pointer       :: OBC !< Open boundaries control structure.
  type(porous_barrier_type), intent(in)  :: pbv !< porous barrier fractional cell metrics

  ! Local variables
  real :: h_W(SZI_(G),SZJ_(G),SZK_(GV)) ! West edge thicknesses in the zonal PPM reconstruction [H ~> m or kg m-2]
  real :: h_E(SZI_(G),SZJ_(G),SZK_(GV)) ! East edge thicknesses in the zonal PPM reconstruction [H ~> m or kg m-2]
  real :: h_S(SZI_(G),SZJ_(G),SZK_(GV)) ! South edge thicknesses in the meridional PPM reconstruction [H ~> m or kg m-2]
  real :: h_N(SZI_(G),SZJ_(G),SZK_(GV)) ! North edge thicknesses in the meridional PPM reconstruction [H ~> m or kg m-2]
  type (box_t) :: bxC                   ! Iteration box for the continuity solver
  type(RealArray_t) :: h_in_a, h_W_a, h_E_a, mask2dT_a
  type(RealArray_t) :: h_S_a, h_N_a
  type(RealArray_t) :: u_a, uhbt_a, por_face_areaU_a
  type(RealArray_t) :: v_a, vhbt_a, por_face_areaV_a
  type(RealArray_t) :: dy_Cu_a, IareaT_a, IdxT_a
  type(RealArray_t) :: dx_Cv_a, IdyT_a
  real :: h_min                         ! Minimum layer thickness (2*Angstrom_H) [H ~> m or kg m-2]

  ! Construct the iteration box
  bxC = set_continuity_box(G,GV, CS)

  call h_in_a%alloc(lb=LBOUND(h), ub=UBOUND(h), source=h)
  call mask2dT_a%alloc(lb=LBOUND(G%mask2dT), ub=UBOUND(G%mask2dT), source=G%mask2dT)
  h_min = 2.0 * GV%Angstrom_H
  call h_W_a%alloc(lb=LBOUND(h_W), ub=UBOUND(h_W), source=h_W)
  call h_E_a%alloc(lb=LBOUND(h_E), ub=UBOUND(h_E), source=h_E)
  call u_a%alloc(lb=LBOUND(u), ub=UBOUND(u), source=u)
  call uhbt_a%alloc(lb=LBOUND(uhbt), ub=UBOUND(uhbt))
  call por_face_areaU_a%alloc(lb=LBOUND(pbv%por_face_areaU), ub=UBOUND(pbv%por_face_areaU), &
                              source=pbv%por_face_areaU)
  call h_S_a%alloc(lb=LBOUND(h_S), ub=UBOUND(h_S), source=h_S)
  call h_N_a%alloc(lb=LBOUND(h_N), ub=UBOUND(h_N), source=h_N)
  call v_a%alloc(lb=LBOUND(v), ub=UBOUND(v), source=v)
  call vhbt_a%alloc(lb=LBOUND(vhbt), ub=UBOUND(vhbt))
  call por_face_areaV_a%alloc(lb=LBOUND(pbv%por_face_areaV), ub=UBOUND(pbv%por_face_areaV), &
                              source=pbv%por_face_areaV)
  call dy_Cu_a%alloc(lb=LBOUND(G%dy_Cu), ub=UBOUND(G%dy_Cu), source=G%dy_Cu)
  call IareaT_a%alloc(lb=LBOUND(G%IareaT), ub=UBOUND(G%IareaT), source=G%IareaT)
  call IdxT_a%alloc(lb=LBOUND(G%IdxT), ub=UBOUND(G%IdxT), source=G%IdxT)
  call dx_Cv_a%alloc(lb=LBOUND(G%dx_Cv), ub=UBOUND(G%dx_Cv), source=G%dx_Cv)
  call IdyT_a%alloc(lb=LBOUND(G%IdyT), ub=UBOUND(G%IdyT), source=G%IdyT)

  call zonal_edge_thickness(bxC, h_in_a, h_W_a, h_E_a, mask2dT_a, &
                            h_min, CS%upwind_1st, CS%monotonic, CS%simple_2nd, OBC)
  call zonal_BT_mass_flux(bxC, u_a, h_in_a, h_W_a, h_E_a, uhbt_a, dt, CS%vol_CFL, &
                          OBC, por_face_areaU_a, dy_Cu_a=dy_Cu_a, IareaT_a=IareaT_a, IdxT_a=IdxT_a)
  call meridional_edge_thickness(bxC, h_in_a, h_S_a, h_N_a, mask2dT_a, &
                                 h_min, CS%upwind_1st, CS%monotonic, CS%simple_2nd, OBC)
  call meridional_BT_mass_flux(bxC, v_a, h_in_a, h_S_a, h_N_a, vhbt_a, dt, &
                               CS%vol_CFL, OBC, por_face_areaV_a, dx_Cv_a=dx_Cv_a, &
                               IareaT_a=IareaT_a, IdyT_a=IdyT_a)

  call uhbt_a%copy2F(uhbt)
  call vhbt_a%copy2F(vhbt)
  call h_in_a%free()
  call mask2dT_a%free()
  call h_W_a%free()
  call h_E_a%free()
  call u_a%free()
  call uhbt_a%free()
  call por_face_areaU_a%free()
  call h_S_a%free()
  call h_N_a%free()
  call v_a%free()
  call vhbt_a%free()
  call por_face_areaV_a%free()
  call dy_Cu_a%free()
  call IareaT_a%free()
  call IdxT_a%free()
  call dx_Cv_a%free()
  call IdyT_a%free()

  ! Free the continuity solver iteration box
  call bxC%free()

end subroutine continuity_2d_fluxes

!> Correct the velocities to give the specified depth-integrated transports by applying a
!! barotropic acceleration (subject to viscous drag) to the velocities.
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
  type(continuity_PPM_CS), intent(in)    :: CS  !< Control structure for mom_continuity.
  type(ocean_OBC_type),    pointer       :: OBC !< Open boundaries control structure.
  type(porous_barrier_type), intent(in)  :: pbv !< porous barrier fractional cell metrics
  real, dimension(SZIB_(G),SZJ_(G)), &
                           intent(in)    :: uhbt !< The vertically summed thickness flux through
                                                !! zonal faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G)), &
                           intent(in)    :: vhbt !< The vertically summed thickness flux through
                                                !! meridional faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(RealArray_t), optional, intent(in) :: visc_rem_u_a !< Both the fraction of the zonal momentum
                                                !! that remains after a time-step of viscosity, and
                                                !! the fraction of a time-step's worth of a barotropic
                                                !! acceleration that a layer experiences after viscosity
                                                !! is applied [nondim].  This goes between 0 (at the
                                                !! bottom) and 1 (far above the bottom).  When this
                                                !! column is under an ice shelf, this also goes to 0
                                                !! at the top due to the no-slip boundary condition there.
  type(RealArray_t), optional, intent(in) :: visc_rem_v_a !< Both the fraction of the meridional
                                                !! momentum that remains after a time-step of
                                                !! viscosity, and the fraction of a time-step's
                                                !! worth of a barotropic acceleration that a layer
                                                !! experiences after viscosity is applied [nondim].
                                                !! This goes between 0 (at the bottom) and 1 (far
                                                !! above the bottom).  When this column is under an
                                                !! ice shelf, this also goes to 0 at the top due to
                                                !! the no-slip boundary condition there.

  ! Local variables
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: u_in  !< Input zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: v_in  !< Input meridional velocity [L T-1 ~> m s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: uh  !< Volume flux through zonal faces =
                                                !! u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: vh  !< Volume flux through meridional faces =
                                                !! v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: h_W(SZI_(G),SZJ_(G),SZK_(GV)) ! West edge thicknesses in the zonal PPM reconstruction [H ~> m or kg m-2]
  real :: h_E(SZI_(G),SZJ_(G),SZK_(GV)) ! East edge thicknesses in the zonal PPM reconstruction [H ~> m or kg m-2]
  real :: h_S(SZI_(G),SZJ_(G),SZK_(GV)) ! South edge thicknesses in the meridional PPM reconstruction [H ~> m or kg m-2]
  real :: h_N(SZI_(G),SZJ_(G),SZK_(GV)) ! North edge thicknesses in the meridional PPM reconstruction [H ~> m or kg m-2]
  type (box_t) :: bxC                   ! Iteration box for continuity solver
  type(RealArray_t) :: h_in_a, h_W_a, h_E_a, mask2dT_a
  type(RealArray_t) :: h_S_a, h_N_a
  type(RealArray_t) :: u_a, uh_a, por_face_areaU_a
  type(RealArray_t) :: uhbt_a, u_cor_a
  type(RealArray_t) :: v_a, vh_a, por_face_areaV_a
  type(RealArray_t) :: vhbt_a, v_cor_a
  type(RealArray_t) :: dy_Cu_a, IareaT_a, IdxT_a, dxCu_a, areaT_a, dxT_a, mask2dCu_a
  type(RealArray_t) :: dx_Cv_a, IdyT_a, dyCv_a, dyT_a, mask2dCv_a
  real :: h_min                         ! Minimum layer thickness (2*Angstrom_H) [H ~> m or kg m-2]

  ! It might not be necessary to separate the input velocity array from the adjusted velocities,
  ! but it seems safer to do so, even if it might be less efficient.
  u_in(:,:,:) = u(:,:,:)
  v_in(:,:,:) = v(:,:,:)

  bxC = set_continuity_box(G,GV, CS)

  h_min = 2.0 * GV%Angstrom_H

  call h_in_a%alloc(lb=LBOUND(h), ub=UBOUND(h), source=h)
  call h_W_a%alloc(lb=LBOUND(h_W), ub=UBOUND(h_W), source=h_W)
  call h_E_a%alloc(lb=LBOUND(h_E), ub=UBOUND(h_E), source=h_E)
  call h_S_a%alloc(lb=LBOUND(h_S), ub=UBOUND(h_S), source=h_S)
  call h_N_a%alloc(lb=LBOUND(h_N), ub=UBOUND(h_N), source=h_N)
  call mask2dT_a%alloc(lb=LBOUND(G%mask2dT), ub=UBOUND(G%mask2dT), source=G%mask2dT)
  call u_a%alloc(lb=LBOUND(u_in), ub=UBOUND(u_in), source=u_in)
  call uh_a%alloc(lb=LBOUND(uh), ub=UBOUND(uh), source=uh)
  call por_face_areaU_a%alloc(lb=LBOUND(pbv%por_face_areaU), ub=UBOUND(pbv%por_face_areaU), &
                              source=pbv%por_face_areaU)
  call uhbt_a%alloc(lb=LBOUND(uhbt), ub=UBOUND(uhbt), source=uhbt)
  call u_cor_a%alloc(lb=LBOUND(u), ub=UBOUND(u), source=u)
  call v_a%alloc(lb=LBOUND(v_in), ub=UBOUND(v_in), source=v_in)
  call vh_a%alloc(lb=LBOUND(vh), ub=UBOUND(vh), source=vh)
  call por_face_areaV_a%alloc(lb=LBOUND(pbv%por_face_areaV), ub=UBOUND(pbv%por_face_areaV), &
                              source=pbv%por_face_areaV)
  call vhbt_a%alloc(lb=LBOUND(vhbt), ub=UBOUND(vhbt), source=vhbt)
  call v_cor_a%alloc(lb=LBOUND(v), ub=UBOUND(v), source=v)
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

  call zonal_edge_thickness(bxC, h_in_a, h_W_a, h_E_a, mask2dT_a, &
                            h_min, CS%upwind_1st, CS%monotonic, CS%simple_2nd, OBC)
  call zonal_mass_flux(bxC, u_a, h_in_a, h_W_a, h_E_a, uh_a, dt, OBC, &
                       por_face_areaU_a, uhbt_a=uhbt_a, visc_rem_u_a=visc_rem_u_a, &
                       u_cor_a=u_cor_a, dy_Cu_a=dy_Cu_a, IareaT_a=IareaT_a, IdxT_a=IdxT_a, &
                       dxCu_a=dxCu_a, areaT_a=areaT_a, dxT_a=dxT_a, mask2dCu_a=mask2dCu_a, &
                       H_subroundoff=GV%H_subroundoff, &
                       CFL_limit_adjust=CS%CFL_limit_adjust, aggress_adjust=CS%aggress_adjust, &
                       use_visc_rem_max=CS%use_visc_rem_max, vol_CFL=CS%vol_CFL, &
                       tol_vel=CS%tol_vel, tol_eta=CS%tol_eta, better_iter=CS%better_iter, &
                       marginal_faces=CS%marginal_faces)
  call meridional_edge_thickness(bxC, h_in_a, h_S_a, h_N_a, mask2dT_a, &
                                 h_min, CS%upwind_1st, CS%monotonic, CS%simple_2nd, OBC)
  call meridional_mass_flux(bxC, v_a, h_in_a, h_S_a, h_N_a, vh_a, dt, OBC, &
                            por_face_areaV_a, vhbt_a=vhbt_a, visc_rem_v_a=visc_rem_v_a, &
                            v_cor_a=v_cor_a, dx_Cv_a=dx_Cv_a, IareaT_a=IareaT_a, &
                            IdyT_a=IdyT_a, dyCv_a=dyCv_a, areaT_a=areaT_a, dyT_a=dyT_a, &
                            mask2dCv_a=mask2dCv_a, H_subroundoff=GV%H_subroundoff, &
                            CFL_limit_adjust=CS%CFL_limit_adjust, &
                            aggress_adjust=CS%aggress_adjust, &
                            use_visc_rem_max=CS%use_visc_rem_max, vol_CFL=CS%vol_CFL, &
                            tol_vel=CS%tol_vel, tol_eta=CS%tol_eta, better_iter=CS%better_iter, &
                            marginal_faces=CS%marginal_faces)

  call u_cor_a%copy2F(u)
  call v_cor_a%copy2F(v)

  call h_in_a%free()
  call h_W_a%free()
  call h_E_a%free()
  call h_S_a%free()
  call h_N_a%free()
  call mask2dT_a%free()
  call u_a%free()
  call uh_a%free()
  call por_face_areaU_a%free()
  call uhbt_a%free()
  call u_cor_a%free()
  call v_a%free()
  call vh_a%free()
  call por_face_areaV_a%free()
  call vhbt_a%free()
  call v_cor_a%free()
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

  ! Free the continuity solver iteration box
  call bxC%free()

end subroutine continuity_adjust_vel


!> Updates the thicknesses due to zonal thickness fluxes.
subroutine continuity_zonal_convergence(bxC, h_a, uh_a, dt, IareaT_a, hin_a, hmin)
  type(box_t),        intent(in)    :: bxC      !< Iteration box for continuity solver
  type(RealArray_t),  intent(inout) :: h_a      !< Final layer thickness [H ~> m or kg m-2]
  type(RealArray_t),  intent(in)    :: uh_a     !< Zonal thickness flux, u*h*dy
                                                !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real,               intent(in)    :: dt       !< Time increment [T ~> s]
  type(RealArray_t),  intent(in)    :: IareaT_a !< The grid cell's 1/areaT [L-2 ~> m-2].
  type(RealArray_t), intent(in) :: hin_a !< Initial layer thickness [H ~> m or kg m-2].
                                              !! If hin is absent, h is also the initial thickness.
  real,              optional, intent(in)    :: hmin !< The minimum layer thickness
                                                     !! [H ~> m or kg m-2]

  real :: h_min  ! The minimum layer thickness [H ~> m or kg m-2].  h_min could be 0.
  integer :: i, j, k, ish, ieh, jsh, jeh, nz
  real, dimension(:,:,:), contiguous, pointer :: h, uh, hin
  real, dimension(:,:), contiguous, pointer :: IareaT

  call cpu_clock_begin(id_clock_update)

  h_min = 0.0 ; if (present(hmin)) h_min = hmin

  call h_a%view(h)
  call uh_a%view(uh)
  call IareaT_a%view(IareaT)

  if (hin_a%associated()) then
    call hin_a%view(hin)
    do concurrent(k=bxC%idxS(3):bxC%idxE(3), &
                  j=bxC%idxS(2):bxC%idxE(2), &
                  i=bxC%idxS(1):bxC%idxE(1))
      h(i,j,k) = max( hin(i,j,k) - dt * IareaT(i,j) * (uh(I,j,k) - uh(I-1,j,k)), h_min )
    enddo
  else
    ! untested
    do concurrent(k=bxC%idxS(3):bxC%idxE(3), &
                  j=bxC%idxS(2):bxC%idxE(2), &
                  i=bxC%idxS(1):bxC%idxE(1))
      h(i,j,k) = max( h(i,j,k) - dt * IareaT(i,j) * (uh(I,j,k) - uh(I-1,j,k)), h_min )
    enddo
  endif

  call cpu_clock_end(id_clock_update)

end subroutine continuity_zonal_convergence

!> Updates the thicknesses due to meridional thickness fluxes.
subroutine continuity_meridional_convergence(bxC, h_a, vh_a, dt, IareaT_a, hin_a, hmin)
  type(box_t),        intent(in)    :: bxC      !< Iteration box for continuity solver
  type(RealArray_t),  intent(inout) :: h_a      !< Final layer thickness [H ~> m or kg m-2]
  type(RealArray_t),  intent(in)    :: vh_a     !< Meridional thickness flux, v*h*dx
                                                !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real,               intent(in)    :: dt       !< Time increment [T ~> s]
  type(RealArray_t),  intent(in)    :: IareaT_a !< The grid cell's 1/areaT [L-2 ~> m-2].
  type(RealArray_t), intent(in) :: hin_a !< Initial layer thickness [H ~> m or kg m-2].
                                              !! If hin is absent, h is also the initial thickness.
  real,              optional, intent(in)    :: hmin !< The minimum layer thickness
                                                     !! [H ~> m or kg m-2]

  real :: h_min  ! The minimum layer thickness [H ~> m or kg m-2].  h_min could be 0.
  integer :: i, j, k, ish, ieh, jsh, jeh, nz
  real, dimension(:,:,:), contiguous, pointer :: h, vh, hin
  real, dimension(:,:), contiguous, pointer :: IareaT

  call cpu_clock_begin(id_clock_update)

  h_min = 0.0 ; if (present(hmin)) h_min = hmin

  call h_a%view(h)
  call vh_a%view(vh)
  call IareaT_a%view(IareaT)

  if (hin_a%associated()) then
    call hin_a%view(hin)
    ! untested
    do concurrent(k=bxC%idxS(3):bxC%idxE(3), &
                  j=bxC%idxS(2):bxC%idxE(2), &
                  i=bxC%idxS(1):bxC%idxE(1))
      h(i,j,k) = max( hin(i,j,k) - dt * IareaT(i,j) * (vh(i,J,k) - vh(i,J-1,k)), h_min )
    enddo
  else
    do concurrent(k=bxC%idxS(3):bxC%idxE(3), &
                  j=bxC%idxS(2):bxC%idxE(2), &
                  i=bxC%idxS(1):bxC%idxE(1))
      h(i,j,k) = max( h(i,j,k) - dt * IareaT(i,j) * (vh(i,J,k) - vh(i,J-1,k)), h_min )
    enddo
  endif

  call cpu_clock_end(id_clock_update)

end subroutine continuity_meridional_convergence


!> Original Fortran implementation of zonal_edge_thickness (renamed). Takes containers.
subroutine zonal_edge_thickness_fortran(bxC, h_in_a, h_W_a, h_E_a, mask2dT_a, &
                                        h_min, upwind_1st, monotonic, simple_2nd, OBC)
  type(Box_t),          intent(in)    :: bxC        !< Iteration box for continuity solver
  type(RealArray_t),    intent(in)    :: h_in_a     !< Tracer cell layer thickness [H ~> m or kg m-2]
  type(RealArray_t),    intent(inout) :: h_W_a      !< Western edge layer thickness [H ~> m or kg m-2]
  type(RealArray_t),    intent(inout) :: h_E_a      !< Eastern edge layer thickness [H ~> m or kg m-2]
  type(RealArray_t),    intent(in)    :: mask2dT_a  !< Cell land/ocean mask [nondim]
  real,                 intent(in)    :: h_min      !< Minimum layer thickness (2*Angstrom_H) [H ~> m or kg m-2]
  logical,              intent(in)    :: upwind_1st !< If true, use 1st-order upwind reconstruction
  logical,              intent(in)    :: monotonic  !< If true, use the CW84 monotonic limiter
  logical,              intent(in)    :: simple_2nd !< If true, use a simple 2nd-order scheme
  type(ocean_OBC_type), pointer       :: OBC        !< Open boundaries control structure

  integer :: i, j, k
  type(Box_t) :: bx
  real, dimension(:,:,:), contiguous, pointer :: h_in, h_W, h_E

  if (upwind_1st) then
    bx = bxC%grow(dim=1, n=1)
    call h_in_a%view(h_in)
    call h_W_a%view(h_W)
    call h_E_a%view(h_E)
    do concurrent (k=bx%idxS(3):bx%idxE(3), j=bx%idxS(2):bx%idxE(2), i=bx%idxS(1):bx%idxE(1))
      h_W(i,j,k) = h_in(i,j,k) ; h_E(i,j,k) = h_in(i,j,k)
    enddo
    call bx%free()
  else
    call PPM_reconstruction_x(bxC, h_in_a, h_W_a, h_E_a, mask2dT_a, h_min, monotonic, simple_2nd, OBC)
  endif

end subroutine zonal_edge_thickness_fortran

!> Shim for zonal_edge_thickness — dispatches via ZONAL_EDGE_THICKNESS_MODE env var.
subroutine zonal_edge_thickness(bxC, h_in_a, h_W_a, h_E_a, mask2dT_a, &
                                h_min, upwind_1st, monotonic, simple_2nd, OBC)
  type(box_t),           intent(in)    :: bxC        !< Iteration box for continuity solver
  type(RealArray_t),     intent(in)    :: h_in_a     !< Tracer cell layer thickness
                                                     !! [H ~> m or kg m-2].
  type(RealArray_t),     intent(inout) :: h_W_a      !< Western edge layer thickness
                                                     !! [H ~> m or kg m-2].
  type(RealArray_t),     intent(inout) :: h_E_a      !< Eastern edge layer thickness
                                                     !! [H ~> m or kg m-2].
  type(RealArray_t),     intent(in)    :: mask2dT_a  !< Cell land/ocean mask [nondim]
  real,                  intent(in)    :: h_min      !< Minimum layer thickness
                                                     !! (2*Angstrom_H) [H ~> m or kg m-2]
  logical,               intent(in)    :: upwind_1st !< If true, use 1st-order upwind reconstruction
  logical,               intent(in)    :: monotonic  !< If true, use the CW84 monotonic limiter
  logical,               intent(in)    :: simple_2nd !< If true, use a simple 2nd-order scheme
  type(ocean_OBC_type),  pointer       :: OBC        !< Open boundaries control structure.

  integer  :: mode, rc
  type(RealArray_C) :: h_in_c, h_W_c, h_E_c, mask2dT_c
  type(Box_C)        :: bxC_c
  type(c_ptr)        :: OBC_c
  logical(c_bool)    :: upwind_1st_c, monotonic_c, simple_2nd_c
  type(io_recorder)  :: rec
  logical            :: capture
  character(len=80)  :: kernel
  character(len=100) :: dir
  character(len=256) :: binFile, metaFile

  kernel = "zonal_edge_thickness"

  call cpu_clock_begin(id_clock_reconstruct)

  mode = getenv_mode("ZONAL_EDGE_THICKNESS_MODE", default=TIMH_runFORTRAN)

  select case (mode)

    case (TIMH_capture)
      capture = (.not. already_recorded(trim(kernel))) .and. is_root_pe()
      if (capture) then
        dir = "capture"
        rc = mkdir_posix(trim(dir) // c_null_char, int(o'755', c_int))
        binFile  = trim(dir) // "/" // trim(kernel) // ".bin"
        metaFile = trim(dir) // "/" // trim(kernel) // ".meta"
        call rec%open_write(binFile, metaFile)
        call rec%add("_bxC",        bxC)
        call rec%add("_h_in",       h_in_a)
        call rec%add("_h_W_before", h_W_a)
        call rec%add("_h_E_before", h_E_a)
        call rec%add("_mask2dT",    mask2dT_a)
        call rec%add("_h_min",      h_min)
        call rec%add("_upwind_1st", upwind_1st)
        call rec%add("_monotonic",  monotonic)
        call rec%add("_simple_2nd", simple_2nd)
      endif

      call zonal_edge_thickness_fortran(bxC, h_in_a, h_W_a, h_E_a, mask2dT_a, &
                                        h_min, upwind_1st, monotonic, simple_2nd, OBC)

      if (capture) then
        call rec%add("_h_W_after", h_W_a)
        call rec%add("_h_E_after", h_E_a)
        call rec%close()
        call mark_recorded(trim(kernel))
      endif

#ifdef _TIM
    case (TIMH_runAMREX)
      bxC_c        = bxC%to_c()
      h_in_c       = h_in_a%to_c()
      h_W_c        = h_W_a%to_c()
      h_E_c        = h_E_a%to_c()
      mask2dT_c    = mask2dT_a%to_c()
      upwind_1st_c = upwind_1st
      monotonic_c  = monotonic
      simple_2nd_c = simple_2nd
      if (associated(OBC)) then
        OBC_c = c_loc(OBC)
      else
        OBC_c = c_null_ptr
      endif
      call turbotmp_zonal_edge_thickness_bridge(bxC_c, h_in_c, h_W_c, h_E_c, mask2dT_c, &
               h_min, upwind_1st_c, monotonic_c, simple_2nd_c, OBC_c)
#endif

    case default
      call zonal_edge_thickness_fortran(bxC, h_in_a, h_W_a, h_E_a, mask2dT_a, &
                                        h_min, upwind_1st, monotonic, simple_2nd, OBC)

  end select

  call cpu_clock_end(id_clock_reconstruct)

end subroutine zonal_edge_thickness


!> Original Fortran implementation of meridional_edge_thickness (renamed). Takes containers.
subroutine meridional_edge_thickness_fortran(bxC, h_in_a, h_S_a, h_N_a, mask2dT_a, &
                                             h_min, upwind_1st, monotonic, simple_2nd, OBC)
  type(Box_t),          intent(in)    :: bxC        !< Iteration box for continuity solver
  type(RealArray_t),    intent(in)    :: h_in_a     !< Tracer cell layer thickness [H ~> m or kg m-2]
  type(RealArray_t),    intent(inout) :: h_S_a      !< Southern edge layer thickness [H ~> m or kg m-2]
  type(RealArray_t),    intent(inout) :: h_N_a      !< Northern edge layer thickness [H ~> m or kg m-2]
  type(RealArray_t),    intent(in)    :: mask2dT_a  !< Cell land/ocean mask [nondim]
  real,                 intent(in)    :: h_min      !< Minimum layer thickness (2*Angstrom_H) [H ~> m or kg m-2]
  logical,              intent(in)    :: upwind_1st !< If true, use 1st-order upwind reconstruction
  logical,              intent(in)    :: monotonic  !< If true, use the CW84 monotonic limiter
  logical,              intent(in)    :: simple_2nd !< If true, use a simple 2nd-order scheme
  type(ocean_OBC_type), pointer       :: OBC        !< Open boundaries control structure

  integer :: i, j, k
  type(Box_t) :: bx
  real, dimension(:,:,:), contiguous, pointer :: h_in, h_S, h_N

  if (upwind_1st) then
    bx = bxC%grow(dim=2, n=1)
    call h_in_a%view(h_in)
    call h_S_a%view(h_S)
    call h_N_a%view(h_N)
    do concurrent (k=bx%idxS(3):bx%idxE(3), j=bx%idxS(2):bx%idxE(2), i=bx%idxS(1):bx%idxE(1))
      h_S(i,j,k) = h_in(i,j,k) ; h_N(i,j,k) = h_in(i,j,k)
    enddo
    call bx%free()
  else
    call PPM_reconstruction_y(bxC, h_in_a, h_S_a, h_N_a, mask2dT_a, h_min, monotonic, simple_2nd, OBC)
  endif

end subroutine meridional_edge_thickness_fortran

!> Shim for meridional_edge_thickness — dispatches via MERIDIONAL_EDGE_THICKNESS_MODE env var.
subroutine meridional_edge_thickness(bxC, h_in_a, h_S_a, h_N_a, mask2dT_a, &
                                     h_min, upwind_1st, monotonic, simple_2nd, OBC)
  type(box_t),           intent(in)    :: bxC        !< Iteration box for continuity solver
  type(RealArray_t),     intent(in)    :: h_in_a     !< Tracer cell layer thickness
                                                     !! [H ~> m or kg m-2].
  type(RealArray_t),     intent(inout) :: h_S_a      !< Southern edge layer thickness
                                                     !! [H ~> m or kg m-2].
  type(RealArray_t),     intent(inout) :: h_N_a      !< Northern edge layer thickness
                                                     !! [H ~> m or kg m-2].
  type(RealArray_t),     intent(in)    :: mask2dT_a  !< Cell land/ocean mask [nondim]
  real,                  intent(in)    :: h_min      !< Minimum layer thickness
                                                     !! (2*Angstrom_H) [H ~> m or kg m-2]
  logical,               intent(in)    :: upwind_1st !< If true, use 1st-order upwind reconstruction
  logical,               intent(in)    :: monotonic  !< If true, use the CW84 monotonic limiter
  logical,               intent(in)    :: simple_2nd !< If true, use a simple 2nd-order scheme
  type(ocean_OBC_type),  pointer       :: OBC        !< Open boundaries control structure.

  integer  :: mode, rc
  type(RealArray_C) :: h_in_c, h_S_c, h_N_c, mask2dT_c
  type(Box_C)        :: bxC_c
  type(c_ptr)        :: OBC_c
  logical(c_bool)    :: upwind_1st_c, monotonic_c, simple_2nd_c
  type(io_recorder)  :: rec
  logical            :: capture
  character(len=80)  :: kernel
  character(len=100) :: dir
  character(len=256) :: binFile, metaFile

  kernel = "meridional_edge_thickness"

  call cpu_clock_begin(id_clock_reconstruct)

  mode = getenv_mode("MERIDIONAL_EDGE_THICKNESS_MODE", default=TIMH_runFORTRAN)

  select case (mode)

    case (TIMH_capture)
      capture = (.not. already_recorded(trim(kernel))) .and. is_root_pe()
      if (capture) then
        dir = "capture"
        rc = mkdir_posix(trim(dir) // c_null_char, int(o'755', c_int))
        binFile  = trim(dir) // "/" // trim(kernel) // ".bin"
        metaFile = trim(dir) // "/" // trim(kernel) // ".meta"
        call rec%open_write(binFile, metaFile)
        call rec%add("_bxC",        bxC)
        call rec%add("_h_in",       h_in_a)
        call rec%add("_h_S_before", h_S_a)
        call rec%add("_h_N_before", h_N_a)
        call rec%add("_mask2dT",    mask2dT_a)
        call rec%add("_h_min",      h_min)
        call rec%add("_upwind_1st", upwind_1st)
        call rec%add("_monotonic",  monotonic)
        call rec%add("_simple_2nd", simple_2nd)
      endif

      call meridional_edge_thickness_fortran(bxC, h_in_a, h_S_a, h_N_a, mask2dT_a, &
                                             h_min, upwind_1st, monotonic, simple_2nd, OBC)

      if (capture) then
        call rec%add("_h_S_after", h_S_a)
        call rec%add("_h_N_after", h_N_a)
        call rec%close()
        call mark_recorded(trim(kernel))
      endif

#ifdef _TIM
    case (TIMH_runAMREX)
      bxC_c        = bxC%to_c()
      h_in_c       = h_in_a%to_c()
      h_S_c        = h_S_a%to_c()
      h_N_c        = h_N_a%to_c()
      mask2dT_c    = mask2dT_a%to_c()
      upwind_1st_c = upwind_1st
      monotonic_c  = monotonic
      simple_2nd_c = simple_2nd
      if (associated(OBC)) then
        OBC_c = c_loc(OBC)
      else
        OBC_c = c_null_ptr
      endif
      call turbotmp_meridional_edge_thickness_bridge(bxC_c, h_in_c, h_S_c, h_N_c, mask2dT_c, &
               h_min, upwind_1st_c, monotonic_c, simple_2nd_c, OBC_c)
#endif

    case default
      call meridional_edge_thickness_fortran(bxC, h_in_a, h_S_a, h_N_a, mask2dT_a, &
                                             h_min, upwind_1st, monotonic, simple_2nd, OBC)

  end select

  call cpu_clock_end(id_clock_reconstruct)

end subroutine meridional_edge_thickness


!> Calculates the mass or volume fluxes through the zonal faces, and other related quantities.
subroutine zonal_mass_flux(bxC, u_a, h_in_a, h_W_a, h_E_a, uh_a, dt, OBC, &
                           por_face_areaU_a, uhbt_a, visc_rem_u_a, u_cor_a, BT_cont, du_cor_a, &
                           dy_Cu_a, IareaT_a, IdxT_a, dxCu_a, areaT_a, dxT_a, mask2dCu_a, &
                           H_subroundoff, CFL_limit_adjust, aggress_adjust, use_visc_rem_max, &
                           vol_CFL, tol_vel, tol_eta, better_iter, marginal_faces)
  type(Box_t),             intent(in)    :: bxC  !< Iteration box for continuity solver
  type(RealArray_t),       intent(in)    :: u_a    !< Zonal velocity [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in)    :: h_in_a !< Layer thickness used to calculate fluxes
                                                   !! [H ~> m or kg m-2].
  type(RealArray_t),       intent(in)    :: h_W_a !< Western edge thicknesses [H ~> m or kg m-2].
  type(RealArray_t),       intent(in)    :: h_E_a !< Eastern edge thicknesses [H ~> m or kg m-2].
  type(RealArray_t),       intent(inout) :: uh_a   !< Volume flux through zonal faces = u*h*dy
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real,                    intent(in)    :: dt   !< Time increment [T ~> s].
  type(ocean_OBC_type),    pointer       :: OBC  !< Open boundaries control structure.
  type(RealArray_t),       intent(in)    :: por_face_areaU_a !< fractional open area of U-faces
                                                             !! [nondim]
  type(RealArray_t), optional, intent(in) :: uhbt_a !< The summed volume flux through zonal faces
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(RealArray_t), optional, intent(in) :: visc_rem_u_a
                     !< The fraction of zonal momentum originally in a layer that remains after a
                     !! time-step of viscosity, and the fraction of a time-step's worth of a
                     !! barotropic acceleration that a layer experiences after viscosity is applied
                     !! [nondim]. Visc_rem_u is between 0 (at the bottom) and 1 (far above the
                     !! bottom).
  type(RealArray_t), optional, intent(inout) :: u_cor_a
                     !< The zonal velocities (u with a barotropic correction)
                     !! that give uhbt as the depth-integrated transport [L T-1 ~> m s-1]
  type(BT_cont_type), optional, pointer  :: BT_cont !< A structure with elements that describe the
                     !! effective open face areas as a function of barotropic flow.
  type(RealArray_t), optional, intent(inout) :: du_cor_a
                                     !< The zonal velocity increments from u that give uhbt
                                     !! as the depth-integrated transports [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in)    :: dy_Cu_a !< The grid cell's unblocked lengths of the
                                                    !! u/v-faces of the h-cell [L ~> m].
  type(RealArray_t),       intent(in)    :: IareaT_a !< The grid cell's 1/areaT [L-2 ~> m-2].
  type(RealArray_t),       intent(in)    :: IdxT_a !< The grid cell's 1/dxT [L-1 ~> m-1].
  type(RealArray_t),       intent(in)    :: dxCu_a !< The dx spacing at u points [L ~> m].
  type(RealArray_t),       intent(in)    :: areaT_a !< The grid cell's area [L2 ~> m2].
  type(RealArray_t),       intent(in)    :: dxT_a !< The dx spacing at h points [L ~> m].
  type(RealArray_t),       intent(in)    :: mask2dCu_a !< 0 for boundary points and 1 for ocean
                                                       !! points on the u grid [nondim].
  real,                    intent(in)    :: H_subroundoff !< A thickness that is so small that it
                     !! can be added to a thickness of Angstrom or larger without changing it at the
                     !! bit level [H ~> m or kg m-2].
  real,                    intent(in)    :: CFL_limit_adjust !< The maximum CFL of the adjusted
                                                             !! velocities [nondim]
  logical,                 intent(in)    :: aggress_adjust !< If true, allow the adjusted velocities
                                                           !! to have a relative CFL change up to
                                                           !! 0.5. False by default.
  logical,                 intent(in)    :: use_visc_rem_max !< If true, use more appropriate
                                                             !! limiting bounds for corrections in
                                                             !! strongly viscous columns.
  logical,                 intent(in)    :: vol_CFL !< If true, use the ratio of the open face
                                                     !! lengths to the tracer cell areas when
                                                     !! estimating CFL numbers. Without
                                                     !! aggress_adjust, the default is false; it is
                                                     !! always true with.
  real,                    intent(in)    :: tol_vel !< The tolerance for barotropic velocity
                                                     !! discrepancies between the barotropic
                                                     !! solution and the sum of the layer
                                                     !! thicknesses [L T-1 ~> m s-1].
  real,                    intent(in)    :: tol_eta !< The tolerance for free-surface height
                                                     !! discrepancies between the barotropic
                                                     !! solution and the sum of the layer
                                                     !! thicknesses [H ~> m or kg m-2].
  logical,                 intent(in)    :: better_iter !< If true, stop corrective iterations
                                                         !! using a velocity-based criterion and
                                                         !! only stop if the iteration is better
                                                         !! than all predecessors.
  logical,                 intent(in)    :: marginal_faces !< If true, use the marginal face areas
                          !! from the continuity solver for use as the weights in the barotropic
                          !! solver. Otherwise use the transport averaged areas.

  ! Local variables
  real :: FAuI  ! A sum of zonal face areas [H L ~> m2 or kg m-1].
  real :: FA_u    ! A sum of zonal face areas [H L ~> m2 or kg m-1].
  real :: I_vrm   ! 1.0 / visc_rem_max [nondim]
  real :: CFL_dt  ! The maximum CFL ratio of the adjusted velocities divided by
                  ! the time step [T-1 ~> s-1].
  real :: I_dt    ! 1.0 / dt [T-1 ~> s-1].
  real :: du_lim  ! The velocity change that give a relative CFL of 1 [L T-1 ~> m s-1].
  real :: dx_E, dx_W ! Effective x-grid spacings to the east and west [L ~> m].
  integer :: i, j, k, ish, ieh, jsh, jeh, n, nz
  integer :: l_seg ! The OBC segment number
  logical :: use_visc_rem, set_BT_cont
  logical :: local_specified_BC, local_open_BC, any_simple_OBC  ! OBC-related logicals
  real, dimension(u_a%lb(1):u_a%ub(1), u_a%lb(2):u_a%ub(2), u_a%lb(3):u_a%ub(3)) :: &
    duhdu ! Partial derivative of uh with u [H L ~> m2 or kg m-1].
  type(RealArray_t) :: uh_tot_0_a, duhdu_tot_0_a, du_a, du_max_CFL_a, du_min_CFL_a
  type(RealArray_t) :: visc_rem_max_a, visc_rem_u_tmp_a
  real, dimension(:,:,:), contiguous, pointer :: u, h_in, h_W, h_E, uh, por_face_areaU
  real, dimension(:,:,:), contiguous, pointer :: visc_rem_u
  real, dimension(:,:),   contiguous, pointer :: du_cor
  real, dimension(:,:),   contiguous, pointer :: &
    du, &         ! Corrective barotropic change in the velocity to give uhbt [L T-1 ~> m s-1].
    du_min_CFL, & ! Lower limit on du correction to avoid CFL violations [L T-1 ~> m s-1]
    du_max_CFL, & ! Upper limit on du correction to avoid CFL violations [L T-1 ~> m s-1]
    duhdu_tot_0, & ! Summed partial derivative of uh with u [H L ~> m2 or kg m-1].
    uh_tot_0, &   ! Summed transport with no barotropic correction [H L2 T-1 ~> m3 s-1 or kg s-1].
    visc_rem_max  ! The column maximum of visc_rem [nondim].
  real, dimension(:,:,:), contiguous, pointer :: &
    visc_rem_u_tmp      ! A 2-D copy of visc_rem_u or an array of 1's [nondim].
  real, dimension(:,:),   contiguous, pointer :: dy_Cu, IareaT, IdxT, areaT, dxT, mask2dCu

  nullify(visc_rem_u, du_cor)

  call cpu_clock_begin(id_clock_correct)

  call u_a%view(u)
  call h_in_a%view(h_in)
  call h_W_a%view(h_W)
  call h_E_a%view(h_E)
  call uh_a%view(uh)
  call por_face_areaU_a%view(por_face_areaU)
  if (present(visc_rem_u_a)) call visc_rem_u_a%view(visc_rem_u)
  if (present(du_cor_a)) call du_cor_a%view(du_cor)
  call dy_Cu_a%view(dy_Cu)
  call IareaT_a%view(IareaT)
  call IdxT_a%view(IdxT)
  call areaT_a%view(areaT)
  call dxT_a%view(dxT)
  call mask2dCu_a%view(mask2dCu)

  use_visc_rem = present(visc_rem_u_a)

  set_BT_cont = .false. ; if (present(BT_cont)) set_BT_cont = (associated(BT_cont))

  local_specified_BC = .false. ; local_open_BC = .false.
  if (associated(OBC)) then ; if (OBC%OBC_pe) then
    local_specified_BC = OBC%specified_u_BCs_exist_globally
    local_open_BC = OBC%open_u_BCs_exist_globally
  endif ; endif

  ish = bxC%idxS(1) ; ieh = bxC%idxE(1) ; jsh = bxC%idxS(2) ; jeh = bxC%idxE(2) ; nz  = bxC%idxE(3)

  CFL_dt = CFL_limit_adjust / dt
  I_dt = 1.0 / dt
  if (aggress_adjust) CFL_dt = I_dt

  call du_a%alloc(lb=[u_a%lb(1), u_a%lb(2)], ub=[u_a%ub(1), u_a%ub(2)])
  call du_min_CFL_a%alloc(lb=[u_a%lb(1), u_a%lb(2)], ub=[u_a%ub(1), u_a%ub(2)])
  call du_max_CFL_a%alloc(lb=[u_a%lb(1), u_a%lb(2)], ub=[u_a%ub(1), u_a%ub(2)])
  call duhdu_tot_0_a%alloc(lb=[u_a%lb(1), u_a%lb(2)], ub=[u_a%ub(1), u_a%ub(2)])
  call uh_tot_0_a%alloc(lb=[u_a%lb(1), u_a%lb(2)], ub=[u_a%ub(1), u_a%ub(2)])
  call visc_rem_max_a%alloc(lb=[u_a%lb(1), u_a%lb(2)], ub=[u_a%ub(1), u_a%ub(2)])
  call visc_rem_u_tmp_a%alloc(lb=[u_a%lb(1), u_a%lb(2), 1], ub=[u_a%ub(1), u_a%ub(2), u_a%ub(3)])
  call du_a%view(du)
  call du_min_CFL_a%view(du_min_CFL)
  call du_max_CFL_a%view(du_max_CFL)
  call duhdu_tot_0_a%view(duhdu_tot_0)
  call uh_tot_0_a%view(uh_tot_0)
  call visc_rem_max_a%view(visc_rem_max)
  call visc_rem_u_tmp_a%view(visc_rem_u_tmp)

  !$omp target enter data &
  !$omp   map(alloc: duhdu)

  do concurrent (j=jsh:jeh)

    if (present(du_cor_a)) then
      do concurrent (i=ish-1:ieh)
        du_cor(i,j) = 0.0
      enddo
    endif

    if (.not.use_visc_rem) then
      do concurrent (k=1:nz, i=ish-1:ieh)
        visc_rem_u_tmp(i,j,k) = 1.0
      enddo
    else
      ! this is expensive
      do concurrent (k=1:nz, i=ish-1:ieh)
        visc_rem_u_tmp(i,j,k) = visc_rem_u(i,j,k)
      enddo
    end if

    ! Set uh and duhdu.
    do concurrent (k=1:nz , I=ish-1:ieh)
      call flux_elem(u(I,j,k), h_in(I,j,k), h_in(I+1,j,k), h_W(I,j,k), h_W(I+1,j,k), h_E(I,j,k), &
                     h_E(I+1,j,k), uh(I,j,k), duhdu(I,j,k), visc_rem_u_tmp(I,j,k), dy_Cu(I,j), &
                     IareaT(I,j), IareaT(I+1,j), IdxT(I,j), IdxT(i+1,j), dt, &
                     vol_CFL, por_face_areaU(I,j,k))
    enddo
    if (local_open_BC) then
      do concurrent (k=1:nz, I=ish-1:ieh)
        call flux_elem_OBC(u(I,j,k), h_in(I,j,k), h_in(I+1,j,k), uh(I,j,k), duhdu(I,j,k), &
                           visc_rem_u_tmp(I,j,k), por_face_areaU(I,j,k), dy_Cu(I,j), &
                           OBC, OBC%segnum_u(I,j))
      enddo
    endif

    ! untested!
    if (local_specified_BC) then
      do concurrent (k=1:nz, I=ish-1:ieh, OBC%segnum_u(I,j) /= 0)
        l_seg = abs(OBC%segnum_u(I,j))
        if (OBC%segment(l_seg)%specified) uh(I,j,k) = OBC%segment(l_seg)%normal_trans(I,j,k)
      enddo
    endif

    if (present(uhbt_a) .or. set_BT_cont) then
      if (use_visc_rem.and.use_visc_rem_max) then
        ! poor performance for nvfortran + do concurrent if k is inside loop
        do concurrent (I=ish-1:ieh)
          visc_rem_max(I,j) = visc_rem_u_tmp(I,j,1)
        enddo
        do k=2,nz ; do concurrent (I=ish-1:ieh)
          visc_rem_max(I,j) = max(visc_rem_max(I,j), visc_rem_u_tmp(I,j,k))
        enddo ; enddo
      else
        do concurrent (i=ish-1:ieh)
          visc_rem_max(i, j) = 1.0
        enddo
      endif
      !   Set limits on du that will keep the CFL number between -1 and 1.
      ! This should be adequate to keep the root bracketed in all cases.
      do concurrent (I=ish-1:ieh)
        I_vrm = 0.0
        if (visc_rem_max(I,j) > 0.0) I_vrm = 1.0 / visc_rem_max(I,j)
        if (vol_CFL) then
          dx_W = ratio_max(areaT(i,j), dy_Cu(I,j), 1000.0*dxT(i,j))
          dx_E = ratio_max(areaT(i+1,j), dy_Cu(I,j), 1000.0*dxT(i+1,j))
        else ; dx_W = dxT(i,j) ; dx_E = dxT(i+1,j) ; endif
        du_max_CFL(I,j) = 2.0* (CFL_dt * dx_W) * I_vrm
        du_min_CFL(I,j) = -2.0 * (CFL_dt * dx_E) * I_vrm
        uh_tot_0(I,j) = 0.0 ; duhdu_tot_0(I,j) = 0.0
      enddo
      ! poor performance for nvfortran + do concurrent if k is inside loop
      do k=1,nz ; do concurrent (I=ish-1:ieh)
        duhdu_tot_0(I,j) = duhdu_tot_0(I,j) + duhdu(I, j, k)
        uh_tot_0(I,j) = uh_tot_0(I,j) + uh(I,j,k)
      enddo ; enddo

      if (use_visc_rem) then
        if (aggress_adjust) then
          ! untested!
          do k=1,nz ; do concurrent (I=ish-1:ieh)
            if (vol_CFL) then
              dx_W = ratio_max(areaT(i,j), dy_Cu(I,j), 1000.0*dxT(i,j))
              dx_E = ratio_max(areaT(i+1,j), dy_Cu(I,j), 1000.0*dxT(i+1,j))
            else ; dx_W = dxT(i,j) ; dx_E = dxT(i+1,j) ; endif

            du_lim = 0.499*((dx_W*I_dt - u(I,j,k)) + MIN(0.0,u(I-1,j,k)))
            if (du_max_CFL(I,j) * visc_rem_u_tmp(I,j,k) > du_lim) &
              du_max_CFL(I,j) = du_lim / visc_rem_u_tmp(I,j,k)

            du_lim = 0.499*((-dx_E*I_dt - u(I,j,k)) + MAX(0.0,u(I+1,j,k)))
            if (du_min_CFL(I,j) * visc_rem_u_tmp(I,j,k) < du_lim) &
              du_min_CFL(I,j) = du_lim / visc_rem_u_tmp(I,j,k)
          enddo ; enddo
        else
          do k=1,nz ; do concurrent (I=ish-1:ieh)
            if (vol_CFL) then
              dx_W = ratio_max(areaT(i,j), dy_Cu(I,j), 1000.0*dxT(i,j))
              dx_E = ratio_max(areaT(i+1,j), dy_Cu(I,j), 1000.0*dxT(i+1,j))
            else ; dx_W = dxT(i,j) ; dx_E = dxT(i+1,j) ; endif

            if (du_max_CFL(I,j) * visc_rem_u_tmp(I,j,k) > dx_W*CFL_dt - u(I,j,k)*mask2dCu(I,j)) &
              du_max_CFL(I,j) = (dx_W*CFL_dt - u(I,j,k)) / visc_rem_u_tmp(I,j,k)
            if (du_min_CFL(I,j) * visc_rem_u_tmp(I,j,k) < -dx_E*CFL_dt - u(I,j,k)*mask2dCu(I,j)) &
              du_min_CFL(I,j) = -(dx_E*CFL_dt + u(I,j,k)) / visc_rem_u_tmp(I,j,k)
          enddo ; enddo
        endif
      else
        ! untested!
        if (aggress_adjust) then
          do k=1,nz ; do concurrent (I=ish-1:ieh)
            if (vol_CFL) then
              dx_W = ratio_max(areaT(i,j), dy_Cu(I,j), 1000.0*dxT(i,j))
              dx_E = ratio_max(areaT(i+1,j), dy_Cu(I,j), 1000.0*dxT(i+1,j))
            else ; dx_W = dxT(i,j) ; dx_E = dxT(i+1,j) ; endif

            du_max_CFL(I,j) = MIN(du_max_CFL(I,j), 0.499 * &
                        ((dx_W*I_dt - u(I,j,k)) + MIN(0.0,u(I-1,j,k))) )
            du_min_CFL(I,j) = MAX(du_min_CFL(I,j), 0.499 * &
                        ((-dx_E*I_dt - u(I,j,k)) + MAX(0.0,u(I+1,j,k))) )
          enddo ; enddo
        else
          do k=1,nz ; do concurrent (I=ish-1:ieh)
            if (vol_CFL) then
              dx_W = ratio_max(areaT(i,j), dy_Cu(I,j), 1000.0*dxT(i,j))
              dx_E = ratio_max(areaT(i+1,j), dy_Cu(I,j), 1000.0*dxT(i+1,j))
            else ; dx_W = dxT(i,j) ; dx_E = dxT(i+1,j) ; endif

            du_max_CFL(I,j) = MIN(du_max_CFL(I,j), dx_W*CFL_dt - u(I,j,k))
            du_min_CFL(I,j) = MAX(du_min_CFL(I,j), -(dx_E*CFL_dt + u(I,j,k)))
          enddo ; enddo
        endif ! CS%agress_adjust
      endif ! use_visc_rem
      do concurrent (I=ish-1:ieh)
        du_max_CFL(I,j) = max(du_max_CFL(I,j),0.0)
        du_min_CFL(I,j) = min(du_min_CFL(I,j),0.0)
      enddo
    endif ! present(uhbt_a) .or. set_BT_cont
  enddo

  call present_uhbt_or_set_BT_cont(bxC, u_a, h_in_a, h_W_a, h_E_a, uh_tot_0_a, duhdu_tot_0_a, &
                                   du_a, du_max_CFL_a, du_min_CFL_a, &
                                   visc_rem_u_tmp_a, visc_rem_max_a, por_face_areaU_a, &
                                   uhbt_a=uhbt_a, uh_a=uh_a, u_cor_a=u_cor_a, &
                                   du_cor_a=du_cor_a, BT_cont=BT_cont, dt=dt, OBC=OBC, &
                                   dy_Cu_a=dy_Cu_a, H_subroundoff=H_subroundoff, &
                                   IareaT_a=IareaT_a, IdxT_a=IdxT_a, dxCu_a=dxCu_a, &
                                   tol_vel=tol_vel, tol_eta=tol_eta, &
                                   better_iter=better_iter, vol_CFL=vol_CFL, &
                                   marginal_faces=marginal_faces)
  call uh_tot_0_a%free()
  call duhdu_tot_0_a%free()
  call du_a%free()
  call du_max_CFL_a%free()
  call du_min_CFL_a%free()
  call visc_rem_max_a%free()
  call visc_rem_u_tmp_a%free()

  !$omp target exit data &
  !$omp   map(release: duhdu)

  call cpu_clock_end(id_clock_correct)

end subroutine zonal_mass_flux

subroutine present_uhbt_or_set_BT_cont(bxC, u_a, h_in_a, h_W_a, h_E_a, uh_tot_0_a, &
                                       duhdu_tot_0_a, du_a, du_max_CFL_a, du_min_CFL_a, &
                                       visc_rem_u_a, visc_rem_max_a, por_face_areaU_a, uhbt_a, &
                                       uh_a, u_cor_a, du_cor_a, BT_cont, dt, OBC, &
                                       dy_Cu_a, H_subroundoff, IareaT_a, IdxT_a, dxCu_a, &
                                       tol_vel, tol_eta, better_iter, vol_CFL, marginal_faces)
  type(Box_t),             intent(in)    :: bxC  !< Iteration box for continuity solver
  type(RealArray_t),       intent(in)    :: u_a    !< Zonal velocity [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in)    :: h_in_a !< Layer thickness used to calculate fluxes
                                                   !! [H ~> m or kg m-2].
  type(RealArray_t),       intent(in)    :: h_W_a !< Western edge thicknesses [H ~> m or kg m-2].
  type(RealArray_t),       intent(in)    :: h_E_a !< Eastern edge thicknesses [H ~> m or kg m-2].
  type(RealArray_t),       intent(in)    :: uh_tot_0_a !< Summed transport with no barotropic
                                                       !! correction [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(RealArray_t),       intent(in)    :: duhdu_tot_0_a !< Summed partial derivative of uh with u
                                                          !! [H L ~> m2 or kg m-1].
  type(RealArray_t),       intent(inout) :: du_a !< Corrective barotropic change in the velocity to
                                                 !! give uhbt [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in)    :: du_max_CFL_a !< Upper limit on du correction to avoid
                                                         !! CFL violations [L T-1 ~> m s-1]
  type(RealArray_t),       intent(in)    :: du_min_CFL_a !< Lower limit on du correction to avoid
                                                         !! CFL violations [L T-1 ~> m s-1]
  type(RealArray_t),       intent(inout) :: uh_a   !< Volume flux through zonal faces = u*h*dy
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(RealArray_t),       intent(in)    :: visc_rem_u_a
                     !< The fraction of zonal momentum originally in a layer that remains after a
                     !! time-step of viscosity, and the fraction of a time-step's worth of a barotropic
                     !! acceleration that a layer experiences after viscosity is applied [nondim].
                     !! Visc_rem_u is between 0 (at the bottom) and 1 (far above the bottom).
  type(RealArray_t),       intent(in)    :: visc_rem_max_a !< The column maximum of visc_rem
                                                           !! [nondim].
  type(RealArray_t),       intent(in)    :: por_face_areaU_a !< fractional open area of U-faces
                                                             !! [nondim]
  type(RealArray_t), optional, intent(inout) :: u_cor_a !< The zonal velocities (u with a barotropic
                                                         !! correction) that give uhbt as the
                                                         !! depth-integrated transport
                                                         !! [L T-1 ~> m s-1]
  type(RealArray_t), optional, intent(inout) :: du_cor_a !< The zonal velocity increments from u
                                                          !! that give uhbt as the depth-integrated
                                                          !! transports [L T-1 ~> m s-1].
  type(BT_cont_type), optional, pointer  :: BT_cont !< A structure with elements that describe the
                     !! effective open face areas as a function of barotropic flow.
  type(RealArray_t), optional, intent(in) :: uhbt_a !< The summed volume flux through zonal faces
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real,                    intent(in)    :: dt  !< Time increment [T ~> s].
  type(ocean_OBC_type),      pointer     :: OBC !< Open boundaries control structure.
  type(RealArray_t),       intent(in)    :: dy_Cu_a !< The grid cell's unblocked lengths of the
                                                    !! u/v-faces of the h-cell [L ~> m].
  real,                    intent(in)    :: H_subroundoff !< A thickness that is so small that it
                     !! can be added to a thickness of Angstrom or larger without changing it at the
                     !! bit level [H ~> m or kg m-2].
  type(RealArray_t),       intent(in)    :: IareaT_a !< The grid cell's 1/areaT [L-2 ~> m-2].
  type(RealArray_t),       intent(in)    :: IdxT_a !< The grid cell's 1/dxT [L-1 ~> m-1].
  type(RealArray_t),       intent(in)    :: dxCu_a !< The dx spacing at u points [L ~> m].
  real,                    intent(in)    :: tol_vel !< The tolerance for barotropic velocity
                                                     !! discrepancies between the barotropic
                                                     !! solution and the sum of the layer
                                                     !! thicknesses [L T-1 ~> m s-1].
  real,                    intent(in)    :: tol_eta !< The tolerance for free-surface height
                                                     !! discrepancies between the barotropic
                                                     !! solution and the sum of the layer
                                                     !! thicknesses [H ~> m or kg m-2].
  logical,                 intent(in)    :: better_iter !< If true, stop corrective iterations
                                                         !! using a velocity-based criterion and
                                                         !! only stop if the iteration is better
                                                         !! than all predecessors.
  logical,                 intent(in)    :: vol_CFL !< If true, use the ratio of the open face
                                                     !! lengths to the tracer cell areas when
                                                     !! estimating CFL numbers. Without
                                                     !! aggress_adjust, the default is false; it is
                                                     !! always true with.
  logical,                 intent(in)    :: marginal_faces !< If true, report the marginal face
                          !! thicknesses; otherwise report transport-averaged thicknesses.
  ! Local variables
  logical, dimension(u_a%lb(1):u_a%ub(1), u_a%lb(2):u_a%ub(2)) :: &
    simple_OBC_pt  ! Indicates points in a row with specified transport OBCs
  logical:: set_BT_cont
  logical:: local_specified_BC, local_Flather_OBC, local_open_BC, any_simple_OBC  ! OBC-related logicals
  integer:: l_seg, i, j, k, n, ish, ieh, jsh, jeh, nz
  real :: FAuI, FA_u
  real, dimension(:,:,:), contiguous, pointer :: u, h_in, h_W, h_E, uh, visc_rem_u, por_face_areaU
  real, dimension(:,:,:), contiguous, pointer :: u_cor
  real, dimension(:,:),   contiguous, pointer :: uh_tot_0, duhdu_tot_0, du, du_max_CFL, du_min_CFL
  real, dimension(:,:),   contiguous, pointer :: visc_rem_max, uhbt, du_cor
  real, dimension(:,:),   contiguous, pointer :: dy_Cu
  logical, dimension(:,:), contiguous, pointer :: do_I
  type(RealArray_t) :: h_u_a
  type(LogicalArray_t) :: do_I_a
  type(RealArray_t) :: no_uh_3d_a ! Never allocated; unassociated data signals uh_3d_a absent.

  call u_a%view(u)
  call h_in_a%view(h_in)
  call h_W_a%view(h_W)
  call h_E_a%view(h_E)
  call uh_tot_0_a%view(uh_tot_0)
  call duhdu_tot_0_a%view(duhdu_tot_0)
  call du_a%view(du)
  call du_max_CFL_a%view(du_max_CFL)
  call du_min_CFL_a%view(du_min_CFL)
  call uh_a%view(uh)
  call visc_rem_u_a%view(visc_rem_u)
  call visc_rem_max_a%view(visc_rem_max)
  call por_face_areaU_a%view(por_face_areaU)
  if (present(uhbt_a)) call uhbt_a%view(uhbt)
  if (present(u_cor_a)) call u_cor_a%view(u_cor)
  if (present(du_cor_a)) call du_cor_a%view(du_cor)
  call dy_Cu_a%view(dy_Cu)

  ish = bxC%idxS(1) ; ieh = bxC%idxE(1) ; jsh = bxC%idxS(2) ; jeh = bxC%idxE(2) ; nz  = bxC%idxE(3)

  set_BT_cont = .false. ; if (present(BT_cont)) set_BT_cont = (associated(BT_cont))

  local_specified_BC = .false. ; local_Flather_OBC = .false. ; local_open_BC = .false.
  if (associated(OBC)) then ; if (OBC%OBC_pe) then
    local_specified_BC = OBC%specified_u_BCs_exist_globally
    local_Flather_OBC = OBC%Flather_u_BCs_exist_globally
    local_open_BC = OBC%open_u_BCs_exist_globally
  endif ; endif

  if (present(uhbt_a) .or. set_BT_cont) then
    call do_I_a%alloc(lb=[u_a%lb(1), u_a%lb(2)], ub=[u_a%ub(1), u_a%ub(2)])
    call do_I_a%view(do_I)
    !$omp target enter data map(alloc: do_I, simple_OBC_pt)
    any_simple_OBC = .false.
    if (local_specified_BC .or. local_Flather_OBC) then
      do concurrent (j=jsh:jeh, I=ish-1:ieh)
        l_seg = abs(OBC%segnum_u(I,j))
        ! Avoid reconciling barotropic/baroclinic transports if transport is specified
        simple_OBC_pt(I,j) = .false.
        if (l_seg /= OBC_NONE) simple_OBC_pt(I,j) = OBC%segment(l_seg)%specified
        do_I(I, j) = .not.simple_OBC_pt(I,j)
        any_simple_OBC = any_simple_OBC .or. simple_OBC_pt(I,j)
      enddo
    else
      do concurrent (j=jsh:jeh, I=ish-1:ieh)
        do_I(I, j) = .true.
      enddo
    endif

    if (present(uhbt_a)) then
      ! Find du and uh.
      call zonal_flux_adjust(bxC, u_a, h_in_a, h_W_a, h_E_a, uh_tot_0_a, duhdu_tot_0_a, du_a, &
                            du_max_CFL_a, du_min_CFL_a, dt, tol_vel, tol_eta, &
                            better_iter, vol_CFL, visc_rem_u_a, do_I_a, por_face_areaU_a, &
                            uhbt_a, uh_a, OBC=OBC, dy_Cu_a=dy_Cu_a, IareaT_a=IareaT_a, &
                            IdxT_a=IdxT_a)

      do concurrent (j=jsh:jeh)
        if (present(u_cor_a)) then
          do concurrent (k=1:nz, I=ish-1:ieh)
            u_cor(I,j,k) = u(I,j,k) + du(I,j) * visc_rem_u(I,j,k)
          enddo
          if (any_simple_OBC) then
            ! untested
            do concurrent (k=1:nz, I=ish-1:ieh, simple_OBC_pt(I,j))
              u_cor(I,j,k) = OBC%segment(abs(OBC%segnum_u(I,j)))%normal_vel(I,j,k)
            enddo
          endif
        endif ! u-corrected

        if (present(du_cor_a)) then
          do concurrent (I=ish-1:ieh)
            du_cor(I,j) = du(I,j)
          enddo
        endif ! du-corrected
      enddo
    endif
    if (set_BT_cont) then
      ! Diagnose the zero-transport correction, du0.
      call zonal_flux_adjust(bxC, u_a, h_in_a, h_W_a, h_E_a, uh_tot_0_a, duhdu_tot_0_a, du_a, &
                            du_max_CFL_a, du_min_CFL_a, dt, tol_vel, tol_eta, &
                            better_iter, vol_CFL, visc_rem_u_a, do_I_a, por_face_areaU_a, &
                            uh_3d_a=no_uh_3d_a, dy_Cu_a=dy_Cu_a, IareaT_a=IareaT_a, IdxT_a=IdxT_a)
      call set_zonal_BT_cont(bxC, u_a, h_in_a, h_W_a, h_E_a, BT_cont, du_a, uh_tot_0_a, &
                             duhdu_tot_0_a, du_max_CFL_a, du_min_CFL_a, dt, vol_CFL, &
                             visc_rem_u_a, visc_rem_max_a, do_I_a, por_face_areaU_a, &
                             dxCu_a, dy_Cu_a, IareaT_a, IdxT_a)
      if (any_simple_OBC) then
        ! untested
        do concurrent (j=jsh:jeh, I=ish-1:ieh)
          ! NOTE: simple_OBC_pt(I, j) should prevent access to segment OBC_NONE
          if (simple_OBC_pt(I,j)) then
            FAuI = H_subroundoff*dy_Cu(I,j)
            do k=1,nz
              l_seg = abs(OBC%segnum_u(I,j))
              if ((abs(OBC%segment(l_seg)%normal_vel(I,j,k)) > 0.0) .and. (OBC%segment(l_seg)%specified)) &
                FAuI = FAuI + OBC%segment(l_seg)%normal_trans(I,j,k) / OBC%segment(l_seg)%normal_vel(I,j,k)
            enddo
            BT_cont%FA_u_W0(I,j) = FAuI ; BT_cont%FA_u_E0(I,j) = FAuI
            BT_cont%FA_u_WW(I,j) = FAuI ; BT_cont%FA_u_EE(I,j) = FAuI
            BT_cont%uBT_WW(I,j) = 0.0 ; BT_cont%uBT_EE(I,j) = 0.0
          endif
        enddo
      endif
    endif
    !$omp target exit data map(release: do_I, simple_OBC_pt)
    call do_I_a%free()
  endif

  ! untested!
  if (local_open_BC .and. set_BT_cont) then
    do n = 1, OBC%number_of_segments
      if (OBC%segment(n)%open .and. OBC%segment(n)%is_E_or_W) then
        I = OBC%segment(n)%HI%IsdB
        if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
          do concurrent (j = OBC%segment(n)%HI%Jsd:OBC%segment(n)%HI%Jed)
            FA_u = 0.0
            do k=1,nz ; FA_u = FA_u + h_in(i,j,k)*(dy_Cu(I,j)*por_face_areaU(I,j,k)) ; enddo
            BT_cont%FA_u_W0(I,j) = FA_u ; BT_cont%FA_u_E0(I,j) = FA_u
            BT_cont%FA_u_WW(I,j) = FA_u ; BT_cont%FA_u_EE(I,j) = FA_u
            BT_cont%uBT_WW(I,j) = 0.0 ; BT_cont%uBT_EE(I,j) = 0.0
          enddo
        else
          do concurrent (j = OBC%segment(n)%HI%Jsd:OBC%segment(n)%HI%Jed)
            FA_u = 0.0
            do k=1,nz ; FA_u = FA_u + h_in(i+1,j,k)*(dy_Cu(I,j)*por_face_areaU(I,j,k)) ; enddo
            BT_cont%FA_u_W0(I,j) = FA_u ; BT_cont%FA_u_E0(I,j) = FA_u
            BT_cont%FA_u_WW(I,j) = FA_u ; BT_cont%FA_u_EE(I,j) = FA_u
            BT_cont%uBT_WW(I,j) = 0.0 ; BT_cont%uBT_EE(I,j) = 0.0
          enddo
        endif
      endif
    enddo
  endif

  if  (set_BT_cont) then ; if (allocated(BT_cont%h_u)) then
    if (present(u_cor_a)) then
      call h_u_a%alloc(lb=LBOUND(BT_cont%h_u), ub=UBOUND(BT_cont%h_u), source=BT_cont%h_u)
      call zonal_flux_thickness(bxC, u_cor_a, h_in_a, h_W_a, h_E_a, h_u_a, dt, dy_Cu_a, IareaT_a, &
                                IdxT_a, vol_CFL, marginal_faces, por_face_areaU_a, OBC, &
                                visc_rem_u_a)
      call h_u_a%copy2F(BT_cont%h_u)
      call h_u_a%free()
    else
      call h_u_a%alloc(lb=LBOUND(BT_cont%h_u), ub=UBOUND(BT_cont%h_u), source=BT_cont%h_u)
      call zonal_flux_thickness(bxC, u_a, h_in_a, h_W_a, h_E_a, h_u_a, dt, dy_Cu_a, IareaT_a, &
                                IdxT_a, vol_CFL, marginal_faces, por_face_areaU_a, OBC, &
                                visc_rem_u_a)
      call h_u_a%copy2F(BT_cont%h_u)
      call h_u_a%free()
    endif
  endif ; endif

end subroutine present_uhbt_or_set_BT_cont


!> Calculates the vertically integrated mass or volume fluxes through the zonal faces.
subroutine zonal_BT_mass_flux(bxC, u_a, h_in_a, h_W_a, h_E_a, uhbt_a, dt, vol_CFL, &
                              OBC, por_face_areaU_a, dy_Cu_a, IareaT_a, IdxT_a)
  type(Box_t),                                intent(in)  :: bxC  !< Iteration box for continuity solver
  type(RealArray_t),       intent(in)  :: u_a  !< Zonal velocity [L T-1 ~> m s-1]
  type(RealArray_t),       intent(in)  :: h_in_a !< Layer thickness used to
                                                 !! calculate fluxes [H ~> m or kg m-2]
  type(RealArray_t),       intent(in)  :: h_W_a !< Western edge thickness in the PPM
                                                !! reconstruction [H ~> m or kg m-2].
  type(RealArray_t),       intent(in)  :: h_E_a !< Eastern edge thickness in the PPM
                                                !! reconstruction [H ~> m or kg m-2].
  type(RealArray_t),       intent(inout) :: uhbt_a !< The summed volume flux through zonal
                                                   !! faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  real,                                       intent(in)  :: dt   !< Time increment [T ~> s].
  logical,                 intent(in)  :: vol_CFL !< If true, use the ratio of the open face
                                                   !! lengths to the tracer cell areas when
                                                   !! estimating CFL numbers. Without
                                                   !! aggress_adjust, the default is false; it is
                                                   !! always true with.
  type(ocean_OBC_type),                       pointer     :: OBC  !< Open boundary condition type
                                                                  !! specifies whether, where, and what
                                                                  !! open boundary conditions are used.
  type(RealArray_t),       intent(in)  :: por_face_areaU_a !< fractional open area of U-faces
                                                           !! [nondim]
  type(RealArray_t),       intent(in)  :: dy_Cu_a !< The grid cell's unblocked lengths of the
                                                   !! u/v-faces of the h-cell [L ~> m].
  type(RealArray_t),       intent(in)  :: IareaT_a !< The grid cell's 1/areaT [L-2 ~> m-2].
  type(RealArray_t),       intent(in)  :: IdxT_a !< The grid cell's 1/dxT [L-1 ~> m-1].

  ! Local variables
  real, dimension(u_a%lb(1):u_a%ub(1), u_a%lb(2):u_a%ub(2), u_a%lb(3):u_a%ub(3)) :: &
    uh      ! Volume flux through zonal faces = u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(u_a%lb(1):u_a%ub(1), u_a%lb(2):u_a%ub(2), u_a%lb(3):u_a%ub(3)) :: &
    duhdu   ! Partial derivative of uh with u [H L ~> m2 or kg m-1].
  integer :: i, j, k, ish, ieh, jsh, jeh, nz, l_seg
  logical :: local_specified_BC
  logical, dimension(u_a%lb(2):u_a%ub(2)) :: OBC_in_row
  real, dimension(:,:,:), contiguous, pointer :: u, h_in, h_W, h_E, por_face_areaU
  real, dimension(:,:),   contiguous, pointer :: uhbt
  real, dimension(:,:),   contiguous, pointer :: dy_Cu, IareaT, IdxT

  call u_a%view(u)
  call h_in_a%view(h_in)
  call h_W_a%view(h_W)
  call h_E_a%view(h_E)
  call uhbt_a%view(uhbt)
  call por_face_areaU_a%view(por_face_areaU)
  call dy_Cu_a%view(dy_Cu)
  call IareaT_a%view(IareaT)
  call IdxT_a%view(IdxT)

  call cpu_clock_begin(id_clock_correct)

  local_specified_BC = .false.
  if (associated(OBC)) then ; if (OBC%OBC_pe) then
    local_specified_BC = OBC%specified_v_BCs_exist_globally
  endif ; endif

  ish = bxC%idxS(1) ; ieh = bxC%idxE(1) ; jsh = bxC%idxS(2) ; jeh = bxC%idxE(2) ; nz  = bxC%idxE(3)

  OBC_in_row(:) = .false.

  uhbt(:,:) = 0.0

  ! Determining whether there are any OBC points outside of the k-loop should be more efficient.
  if (local_specified_BC) then
    do j=jsh,jeh ; do I=ish-1,ieh ; if (OBC%segnum_u(I,j) /= 0) then
      if (OBC%segment(abs(OBC%segnum_u(I,j)))%specified) OBC_in_row(j) = .true.
    endif ; enddo ; enddo
  endif

  ! This sets uh and duhdu.
  do concurrent (k=1:nz, j=jsh:jeh, I=ish-1:ieh)
    call flux_elem(u(I,j,k), h_in(I,j,k), h_in(I+1,j,k), h_W(I,j,k), h_W(I+1,j,k), h_E(I,j,k), &
                   h_E(I+1,j,k), uh(I,j,k), duhdu(I,j,k), 1.0, dy_Cu(I,j), IareaT(I,j), &
                   IareaT(I+1,j), IdxT(I,j), IdxT(I+1,j), dt, vol_CFL, &
                   por_face_areaU(I,j,k))
    if (local_specified_BC) &
      call flux_elem_OBC(u(I,j,k), h_in(I,j,k), h_in(I+1,j,k), uh(I,j,k), duhdu(I,j,k), 1.0, &
                         por_face_areaU(I,j,k), dy_Cu(I,j), OBC, OBC%segnum_u(I,j))
  enddo

  do k=1,nz ; do j=jsh,jeh ; do i=ish-1,ieh
    if (OBC_in_row(j) .and. OBC%segnum_u(I,j) /= 0) then
      l_seg = abs(OBC%segnum_u(I,j))
      if (OBC%segment(l_seg)%specified) uh(I,j,k) = OBC%segment(l_seg)%normal_trans(I,j,k)
    endif
  enddo ; enddo ; enddo

  ! Accumulate the barotropic transport.
  do k=1,nz ; do j=jsh,jeh ; do I=ish-1,ieh
        uhbt(I,j) = uhbt(I,j) + uh(I,j,k)
  enddo ; enddo ; enddo ! j-loop

  call cpu_clock_end(id_clock_correct)

end subroutine zonal_BT_mass_flux

!> Evaluates the zonal mass or volume fluxes in an element.
elemental subroutine flux_elem(u, h, h_p1, h_L, h_L_p1, h_R, h_R_p1, uh, duhdu, visc_rem, &
                               G_dy_Cu, G_IareaT, G_IareaT_p1, G_IdxT, G_IdxT_p1, dt, &
                               vol_CFL, por_face_area)
  real,                    intent(in)  :: u        !< Zonal or meridional velocity [L T-1 ~> m s-1].
  real,                    intent(in)  :: visc_rem !< Both the fraction of the
                        !! momentum originally in a layer that remains after a time-step
                        !! of viscosity, and the fraction of a time-step's worth of a barotropic
                        !! acceleration that a layer experiences after viscosity is applied [nondim].
                        !! Visc_rem is between 0 (at the bottom) and 1 (far above the bottom).
  real,                    intent(in)  :: h        !< Layer thickness [H ~> m or kg m-2].
  real,                    intent(in)  :: h_p1     !< Layer thickness - offset by 1 [ H ~> m or kg m-2].
  real,                    intent(in)  :: h_L      !< West/South edge thickness [H ~> m or kg m-2].
  real,                    intent(in)  :: h_L_p1   !< West/South edge thickness - offset by 1 [H ~> m or kg m-2].
  real,                    intent(in)  :: h_R      !< East/North edge thickness [H ~> m or kg m-2].
  real,                    intent(in)  :: h_R_p1   !< East/North edge thickness - offset by 1 [H ~> m or kg m-2].
  real,                    intent(out) :: uh       !< Zonal or meridional mass or volume transport
                                                   !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real,                    intent(out) :: duhdu    !< Partial derivative of uh
                                                   !! with u [H L ~> m2 or kg m-1].
  real,                    intent(in)  :: dt       !< Time increment [T ~> s]
  logical,                 intent(in)  :: vol_CFL  !< If true, rescale the ratio of face areas to the
                                                   !! cell areas when estimating the CFL number.
  real,                    intent(in)  :: por_face_area !< fractional open area of U/V-faces [nondim].
  real,                    intent(in)  :: G_dy_Cu  !< The grid cell's unblocked lengths of the u/v-faces
                                                   !! of the h-cell [L ~> m].
  real,                    intent(in)  :: G_IareaT !< The grid cell's 1/areaT [L-2 ~> m-2].
  real,                    intent(in)  :: G_IareaT_p1 !< The grid cell's 1/areaT - offset by 1 [L-2 ~> m-2].
  real,                    intent(in)  :: G_IdxT   !< The grid cell's 1/dxT [L-1 ~> m-1].
  real,                    intent(in)  :: G_IdxT_p1 !< The grid cell's 1/dxT - offset by 1 [L-1 ~> m-1].
  ! Local variables
  real :: CFL  ! The CFL number based on the local velocity and grid spacing [nondim]
  real :: curv_3 ! A measure of the thickness curvature over a grid length [H ~> m or kg m-2]
  real :: h_marg ! The marginal thickness of a flux [H ~> m or kg m-2].
  real :: tmp ! temporary variable to store precalculted values
  real :: dh ! h differential between E/W

  ! Set new values of uh and duhdu.
  tmp = G_dy_Cu * por_face_area ! precalculate things
  if (u > 0.0) then
    if (vol_CFL) then ; CFL = (u * dt) * (G_dy_Cu * G_IareaT)
    else ; CFL = u * dt * G_IdxT ; endif
    curv_3 = (h_L + h_R) - 2.0*h
    dh = h_L - h_R
    uh = tmp * u * &
        (h_R + CFL * (0.5*dh + curv_3*(CFL - 1.5)))
    h_marg = h_R + CFL * (dh + 3.0*curv_3*(CFL - 1.0))
  elseif (u < 0.0) then
    if (vol_CFL) then ; CFL = (-u * dt) * (G_dy_Cu * G_IareaT_p1)
    else ; CFL = -u * dt * G_IdxT_p1 ; endif
    curv_3 = (h_L_p1 + h_R_p1) - 2.0*h_p1
    dh = h_R_p1-h_L_p1
    uh = tmp * u * &
        (h_L_p1 + CFL * (0.5*dh + curv_3*(CFL - 1.5)))
    h_marg = h_L_p1 + CFL * (dh + 3.0*curv_3*(CFL - 1.0))
  else
    uh = 0.0
    h_marg = 0.5 * (h_L_p1 + h_R)
  endif
  duhdu = tmp * h_marg * visc_rem

end subroutine flux_elem

elemental subroutine flux_elem_OBC(u, h, h_p1, uh, duhdu, visc_rem, por_face_area, &
                                     G_dy_Cu, OBC, l_seg)
  real,                     intent(in)    :: u        !< Zonal/meridional velocity [L T-1 ~> m s-1].
  real,                     intent(in)    :: visc_rem !< Both the fraction of the
                        !! momentum originally in a layer that remains after a time-step
                        !! of viscosity, and the fraction of a time-step's worth of a barotropic
                        !! acceleration that a layer experiences after viscosity is applied [nondim].
                        !! Visc_rem is between 0 (at the bottom) and 1 (far above the bottom).
  real,                     intent(in)    :: h        !< Layer thickness [H ~> m or kg m-2].
  real,                     intent(in)    :: h_p1     !< Layer thickness offset by 1 [H ~> m or kg m-2].
  real,                     intent(inout) :: uh       !< Zonal/meridional mass or volume
                                                      !! transport [H L2 T-1 ~> m3 s-1 or kg s-1].
  real,                     intent(inout) :: duhdu    !< Partial derivative of uh
                                                      !! with u [H L ~> m2 or kg m-1].
  real,                     intent(in)    :: por_face_area !< fractional open area of U/V-faces
                                                            !! [nondim].
  real,                     intent(in)    :: G_dy_Cu  !< The grid cell's unblocked lengths of the
                                                      !! u/v-faces of the h-cell [L ~> m].
          !! ratio of face areas to the cell areas when estimating the CFL number.
  type(ocean_OBC_type),     intent(in)    :: OBC !< Open boundaries control structure.
  integer, intent(in) :: l_seg !< Segment index.

  ! untested
  if (l_seg /= 0) then
    if (OBC%segment(abs(l_seg))%open) then
      if (l_seg > 0) then !  OBC_DIRECTION_E or OBC_DIRECTION_N
        uh = (G_dy_Cu * por_face_area) * u * h
        duhdu = (G_dy_Cu * por_face_area) * h * visc_rem
      else !  OBC_DIRECTION_W or OBC_DIRECTION_S
        uh = (G_dy_Cu * por_face_area) * u * h_p1
        duhdu = (G_dy_Cu* por_face_area) * h_p1 * visc_rem
      endif
    endif
  endif

end subroutine flux_elem_OBC


!> Sets the effective interface thickness associated with the fluxes at each zonal velocity point,
!! optionally scaling back these thicknesses to account for viscosity and fractional open areas.
subroutine zonal_flux_thickness(bxC, u_a, h_a, h_W_a, h_E_a, h_u_a, dt, dy_Cu_a, IareaT_a, &
                                IdxT_a, vol_CFL, marginal, por_face_areaU_a, OBC, visc_rem_u_a)
  type(box_t),                               intent(in)    :: bxC  !< Iteration box for continuity solver
  type(RealArray_t),       intent(in)    :: u_a  !< Zonal velocity [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in)    :: h_a  !< Layer thickness used to calculate fluxes
                                                 !! [H ~> m or kg m-2].
  type(RealArray_t),       intent(in)    :: h_W_a !< West edge thickness in the reconstruction
                                                  !! [H ~> m or kg m-2].
  type(RealArray_t),       intent(in)    :: h_E_a !< East edge thickness in the reconstruction
                                                  !! [H ~> m or kg m-2].
  type(RealArray_t),       intent(inout) :: h_u_a !< Effective thickness at zonal faces,
                                                  !! scaled down to account for the effects of
                                                  !! viscosity and the fractional open area
                                                  !! [H ~> m or kg m-2].
  real,                                      intent(in)    :: dt   !< Time increment [T ~> s].
  type(RealArray_t),       intent(in)    :: dy_Cu_a !< The grid cell's unblocked lengths of the
                                                    !! u/v-faces of the h-cell [L ~> m].
  type(RealArray_t),       intent(in)    :: IareaT_a !< The grid cell's 1/areaT [L-2 ~> m-2].
  type(RealArray_t),       intent(in)    :: IdxT_a !< The grid cell's 1/dxT [L-1 ~> m-1].
  logical,                                   intent(in)    :: vol_CFL !< If true, rescale the ratio
                          !! of face areas to the cell areas when estimating the CFL number.
  logical,                                   intent(in)    :: marginal !< If true, report the
                          !! marginal face thicknesses; otherwise report transport-averaged thicknesses.
  type(RealArray_t),       intent(in)    :: por_face_areaU_a !< fractional open area of U-faces
                                                             !! [nondim]
  type(ocean_OBC_type),                      pointer       :: OBC !< Open boundaries control structure.
  type(RealArray_t), intent(in) :: visc_rem_u_a
                          !< Both the fraction of the momentum originally in a layer that remains after
                          !! a time-step of viscosity, and the fraction of a time-step's worth of a
                          !! barotropic acceleration that a layer experiences after viscosity is applied [nondim].
                          !! Visc_rem_u is between 0 (at the bottom) and 1 (far above the bottom).

  ! Local variables
  real :: CFL  ! The CFL number based on the local velocity and grid spacing [nondim]
  real :: curv_3 ! A measure of the thickness curvature over a grid length [H ~> m or kg m-2]
  logical :: local_open_BC
  integer :: i, j, k, ish, ieh, jsh, jeh, nz, n
  real :: dh
  type(box_t) :: bxU
  real, dimension(:,:,:), contiguous, pointer :: u, h, h_W, h_E, h_u, por_face_areaU, visc_rem_u
  real, dimension(:,:),   contiguous, pointer :: dy_Cu, IareaT, IdxT

  call u_a%view(u)
  call h_a%view(h)
  call h_W_a%view(h_W)
  call h_E_a%view(h_E)
  call h_u_a%view(h_u)
  call dy_Cu_a%view(dy_Cu)
  call IareaT_a%view(IareaT)
  call IdxT_a%view(IdxT)
  call por_face_areaU_a%view(por_face_areaU)
  if (visc_rem_u_a%associated()) call visc_rem_u_a%view(visc_rem_u)

  ish = bxC%idxS(1) ; ieh = bxC%idxE(1) ; jsh = bxC%idxS(2) ; jeh = bxC%idxE(2) ; nz  = bxC%idxE(3)

  bxU = bxC%growLo(dim=1, n=1) !< Increase the lower extent of the x-dimension (U-grid)

  ! do concurrent (k=1:nz, j=jsh:jeh, I=ish-1:ieh)
  do concurrent(k=bxU%idxS(3):bxU%idxE(3), &
                j=bxU%idxS(2):bxU%idxE(2), &
                i=bxU%idxS(1):bxU%idxE(1)) ! U-grid
    if (u(I,j,k) > 0.0) then
      if (vol_CFL) then ; CFL = (u(I,j,k) * dt) * (dy_Cu(I,j) * IareaT(i,j))
      else ; CFL = u(I,j,k) * dt * IdxT(i,j) ; endif
      curv_3 = (h_W(i,j,k) + h_E(i,j,k)) - 2.0*h(i,j,k)
      dh = h_W(i,j,k) - h_E(i,j,k)
      if (marginal) then
        h_u(I,j,k) = h_E(i,j,k) + CFL * (dh + 3.0*curv_3*(CFL - 1.0))
      else
        h_u(I,j,k) = h_E(i,j,k) + CFL * (0.5*dh + curv_3*(CFL - 1.5))
      endif
    elseif (u(I,j,k) < 0.0) then
      if (vol_CFL) then ; CFL = (-u(I,j,k)*dt) * (dy_Cu(I,j) * IareaT(i+1,j))
      else ; CFL = -u(I,j,k) * dt * IdxT(i+1,j) ; endif
      curv_3 = (h_W(i+1,j,k) + h_E(i+1,j,k)) - 2.0*h(i+1,j,k)
      dh = h_E(i+1,j,k)-h_W(i+1,j,k)
      if (marginal) then
        h_u(I,j,k) = h_W(i+1,j,k) + CFL * (dh + 3.0*curv_3*(CFL - 1.0))
      else
        h_u(I,j,k) = h_W(i+1,j,k) + CFL * (0.5*dh + curv_3*(CFL - 1.5))
      endif
    else
      !   The choice to use the arithmetic mean here is somewhat arbitrarily, but
      ! it should be noted that h_W(i+1,j,k) and h_E(i,j,k) are usually the same.
      h_u(I,j,k) = 0.5 * (h_W(i+1,j,k) + h_E(i,j,k))
 !    h_marg = (2.0 * h_W(i+1,j,k) * h_E(i,j,k)) / &
 !             (h_W(i+1,j,k) + h_E(i,j,k) + GV%H_subroundoff)
    endif

    if (visc_rem_u_a%associated()) then
      ! Scale back the thickness to account for the effects of viscosity and the fractional open
      ! thickness to give an appropriate non-normalized weight for each layer in determining the
      ! barotropic acceleration.
      h_u(I,j,k) = h_u(I,j,k) * (visc_rem_u(I,j,k) * por_face_areaU(I,j,k))
    else
      h_u(I,j,k) = h_u(I,j,k) * por_face_areaU(I,j,k)
    endif
  enddo

  local_open_BC = .false.
  if (associated(OBC)) local_open_BC = OBC%open_u_BCs_exist_globally
  if (local_open_BC) then
    ! untested
    do n = 1, OBC%number_of_segments
      if (OBC%segment(n)%open .and. OBC%segment(n)%is_E_or_W) then
        I = OBC%segment(n)%HI%IsdB
        if (OBC%segment(n)%direction == OBC_DIRECTION_E) then
          if (visc_rem_u_a%associated()) then
            do concurrent (k=1:nz, j = OBC%segment(n)%HI%jsd:OBC%segment(n)%HI%jed)
              h_u(I,j,k) = h(i,j,k) * (visc_rem_u(I,j,k) * por_face_areaU(I,j,k))
            enddo
          else
            do concurrent (k=1:nz, j = OBC%segment(n)%HI%jsd:OBC%segment(n)%HI%jed)
              h_u(I,j,k) = h(i,j,k) * por_face_areaU(I,j,k)
            enddo
          endif
        else
          if (visc_rem_u_a%associated()) then
            do concurrent (k=1:nz, j = OBC%segment(n)%HI%jsd:OBC%segment(n)%HI%jed)
              h_u(I,j,k) = h(i+1,j,k) * (visc_rem_u(I,j,k) * por_face_areaU(I,j,k))
            enddo
          else
            do concurrent (k=1:nz, j = OBC%segment(n)%HI%jsd:OBC%segment(n)%HI%jed)
              h_u(I,j,k) = h(i+1,j,k) * por_face_areaU(I,j,k)
            enddo
          endif
        endif
      endif
    enddo
  endif

  call bxU%free()

end subroutine zonal_flux_thickness

!> Returns the barotropic velocity adjustment that gives the
!! desired barotropic (layer-summed) transport.
subroutine zonal_flux_adjust(bxC, u_a, h_in_a, h_W_a, h_E_a, uh_tot_0_a, duhdu_tot_0_a, &
                             du_a, du_max_CFL_a, du_min_CFL_a, dt, tol_vel_in, &
                             tol_eta_in, better_iter, vol_CFL, visc_rem_a, do_I_in_a, &
                             por_face_areaU_a, uhbt_a, uh_3d_a, OBC, dy_Cu_a, IareaT_a, IdxT_a)
  type(box_t),                                intent(in)    :: bxC  !< Iteration box for continuity solver
  type(RealArray_t),       intent(in)    :: u_a  !< Zonal velocity [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in)    :: h_in_a !< Layer thickness used to calculate fluxes
                                                   !! [H ~> m or kg m-2].
  type(RealArray_t),       intent(in)    :: h_W_a !< West edge thickness in the reconstruction
                                                  !! [H ~> m or kg m-2].
  type(RealArray_t),       intent(in)    :: h_E_a !< East edge thickness in the reconstruction
                                                  !! [H ~> m or kg m-2].
  type(RealArray_t),       intent(in)    :: visc_rem_a !< Both the fraction of the momentum
                                                       !! originally in a layer that remains after a
                                                       !! time-step of viscosity, and the fraction
                                                       !! of a time-step's worth of a barotropic
                                                       !! acceleration that a layer experiences
                                                       !! after viscosity is applied [nondim].
                                                       !! Visc_rem is between 0 (at the bottom) and
                                                       !! 1 (far above the bottom).
  type(RealArray_t), optional, intent(in)    :: uhbt_a !< The summed volume flux
                       !! through zonal faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(RealArray_t),       intent(in)    :: du_max_CFL_a  !< Maximum acceptable
                       !! value of du [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in)    :: du_min_CFL_a  !< Minimum acceptable
                       !! value of du [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in)    :: uh_tot_0_a !< The summed transport with 0 adjustment
                                                       !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(RealArray_t),       intent(in)    :: duhdu_tot_0_a !< The partial derivative of du_err with
                                                          !! du at 0 adjustment
                                                          !! [H L ~> m2 or kg m-1].
  type(RealArray_t),       intent(inout) :: du_a !<
                       !! The barotropic velocity adjustment [L T-1 ~> m s-1].
  real,                                       intent(in)    :: dt  !< Time increment [T ~> s].
  real,                    intent(in)    :: tol_vel_in !< The tolerance for barotropic velocity
                                                       !! discrepancies between the barotropic
                                                       !! solution and the sum of the layer
                                                       !! thicknesses [L T-1 ~> m s-1].
  real,                    intent(in)    :: tol_eta_in !< The tolerance for free-surface height
                                                       !! discrepancies between the barotropic
                                                       !! solution and the sum of the layer
                                                       !! thicknesses [H ~> m or kg m-2].
  logical,                 intent(in)    :: better_iter !< If true, stop corrective iterations
                                                        !! using a velocity-based criterion and
                                                        !! only stop if the iteration is better
                                                        !! than all predecessors.
  logical,                 intent(in)    :: vol_CFL !< If true, use the ratio of the open face
                                                     !! lengths to the tracer cell areas when
                                                     !! estimating CFL numbers. Without
                                                     !! aggress_adjust, the default is false; it is
                                                     !! always true with.


  type(LogicalArray_t),    intent(in)    :: do_I_in_a !< A logical flag indicating which I
                                                      !! values to work on.
  type(RealArray_t),       intent(in)    :: por_face_areaU_a !< fractional open area
                                                             !! of U-faces [nondim].
  type(RealArray_t), intent(inout) :: uh_3d_a !< Volume flux through zonal
                                                 !! faces = u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(ocean_OBC_type),             optional, pointer       :: OBC !< Open boundaries control structure.
  type(RealArray_t),       intent(in)    :: dy_Cu_a !< The grid cell's unblocked lengths of the
                                                    !! u/v-faces of the h-cell [L ~> m].
  type(RealArray_t),       intent(in)    :: IareaT_a !< The grid cell's 1/areaT [L-2 ~> m-2].
  type(RealArray_t),       intent(in)    :: IdxT_a !< The grid cell's 1/dxT [L-1 ~> m-1].
  ! Local variables
  real, dimension(u_a%lb(1):u_a%ub(1), u_a%lb(3):u_a%ub(3)) :: &
    uh_aux         ! An auxiliary zonal volume flux [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: &
    duhdu, &       ! Partial derivative of uh with u [H L ~> m2 or kg m-1].
    u_new          ! The velocity with the correction added [L T-1 ~> m s-1].
  real, dimension(u_a%lb(1):u_a%ub(1)) :: &
    uh_err, &      ! Difference between uhbt and the summed uh [H L2 T-1 ~> m3 s-1 or kg s-1].
    uh_err_best, & ! The smallest value of uh_err found so far [H L2 T-1 ~> m3 s-1 or kg s-1].
    duhdu_tot,&    ! Summed partial derivative of uh with u [H L ~> m2 or kg m-1].
    du_min, &      ! Lower limit on du correction based on CFL limits and previous iterations [L T-1 ~> m s-1]
    du_max         ! Upper limit on du correction based on CFL limits and previous iterations [L T-1 ~> m s-1]
  real :: du_prev  ! The previous value of du [L T-1 ~> m s-1].
  real :: ddu      ! The change in du from the previous iteration [L T-1 ~> m s-1].
  real :: tol_eta  ! The tolerance for the current iteration [H ~> m or kg m-2].
  real :: tol_vel  ! The tolerance for velocity in the current iteration [L T-1 ~> m s-1].
  integer :: i, j, k, nz, itt
  integer :: ish !< Start of i index range.
  integer :: jsh !< Start of j index range.
  integer :: ieh !< End of i index range.
  integer :: jeh !< End of j index range.
  logical :: do_I(u_a%lb(1):u_a%ub(1)), local_OBC, use_uhbt
  integer, parameter:: max_itts = 20
  real, dimension(:,:,:), contiguous, pointer :: u, h_in, h_W, h_E, visc_rem, por_face_areaU, uh_3d
  real, dimension(:,:),   contiguous, pointer :: uh_tot_0, duhdu_tot_0, du, du_max_CFL, du_min_CFL
  real, dimension(:,:),   contiguous, pointer :: uhbt
  real, dimension(:,:),   contiguous, pointer :: dy_Cu, IareaT, IdxT
  logical, dimension(:,:), contiguous, pointer :: do_I_in

  call u_a%view(u)
  call h_in_a%view(h_in)
  call h_W_a%view(h_W)
  call h_E_a%view(h_E)
  call visc_rem_a%view(visc_rem)
  call uh_tot_0_a%view(uh_tot_0)
  call duhdu_tot_0_a%view(duhdu_tot_0)
  call du_a%view(du)
  call du_max_CFL_a%view(du_max_CFL)
  call du_min_CFL_a%view(du_min_CFL)
  call por_face_areaU_a%view(por_face_areaU)
  if (present(uhbt_a)) call uhbt_a%view(uhbt)
  if (uh_3d_a%associated()) call uh_3d_a%view(uh_3d)
  call dy_Cu_a%view(dy_Cu)
  call IareaT_a%view(IareaT)
  call IdxT_a%view(IdxT)
  call do_I_in_a%view(do_I_in)


  local_OBC = .false.
  if (present(OBC)) then
    if (associated(OBC)) then
      local_OBC = OBC%open_u_BCs_exist_globally
    endif
  endif

  use_uhbt = present(uhbt_a)

  ish = bxC%idxS(1) ; ieh = bxC%idxE(1) ; jsh = bxC%idxS(2) ; jeh = bxC%idxE(2) ; nz  = bxC%idxE(3)

  tol_vel = tol_vel_in

  ! NVIDIA needs private arrays to be alloc'ed to prevent data transfers.
  ! GCC doesn't understand map(alloc: ...) for variables also marked private
  !$omp target enter data map(alloc: do_I, du_max, du_min, duhdu_tot, uh_err, uh_err_best, uh_aux)

  ! NVIDIA do concurrent doesn't work with private arrays (private scalars OK)
  !$omp target teams loop &
  !$omp   private(uh_err, uh_err_best, duhdu_tot, du_min, du_max, do_I, uh_aux, itt, tol_eta)
  do j=jsh,jeh

    if (uh_3d_a%associated()) then
      do concurrent (k=1:nz, I=ish-1:ieh)
        uh_aux(I,k) = uh_3d(I,j,k)
      enddo
    endif

    do concurrent (I=ish-1:ieh)
      du(I,j) = 0.0 ; do_I(I) = do_I_in(I,j)
      du_max(I) = du_max_CFL(I,j) ; du_min(I) = du_min_CFL(I,j)
      uh_err(I) = uh_tot_0(I,j)
      if (use_uhbt) uh_err(I) = uh_err(I) - uhbt(I,j)
      duhdu_tot(I) = duhdu_tot_0(I,j)
      uh_err_best(I) = abs(uh_err(I))
    enddo

    do itt=1,max_itts
      select case (itt)
        case (:1) ; tol_eta = 1e-6 * tol_eta_in
        case (2)  ; tol_eta = 1e-4 * tol_eta_in
        case (3)  ; tol_eta = 1e-2 * tol_eta_in
        case default ; tol_eta = tol_eta_in
      end select

      do concurrent (I=ish-1:ieh, do_I(I)) &
          & DO_LOCALITY(local(ddu, du_prev))
        if (uh_err(I) > 0.0) then ; du_max(I) = du(I,j)
        elseif (uh_err(I) < 0.0) then ; du_min(I) = du(I,j)
        else ; do_I(I) = .false. ; endif
        if ((dt * min(IareaT(i,j),IareaT(i+1,j))*abs(uh_err(I)) > tol_eta) .or. &
            (better_iter .and. ((abs(uh_err(I)) > tol_vel * duhdu_tot(I)) .or. &
                                  (abs(uh_err(I)) > uh_err_best(I))) )) then
        !   Use Newton's method, provided it stays bounded.  Otherwise bisect
        ! the value with the appropriate bound.
          ddu = -uh_err(I) / duhdu_tot(I)
          du_prev = du(I,j)
          du(I,j) = du(I,j) + ddu
          if (abs(ddu) < 1.0e-15*abs(du(I,j))) then
            do_I(I) = .false. ! ddu is small enough to quit.
          elseif (ddu > 0.0) then
            if (du(I,j) >= du_max(I)) then
              du(I,j) = 0.5*(du_prev + du_max(I))
              if (du_max(I) - du_prev < 1.0e-15*abs(du(I,j))) do_I(I) = .false.
            endif
          else ! ddu < 0.0
            if (du(I,j) <= du_min(I)) then
              du(I,j) = 0.5*(du_prev + du_min(I))
              if (du_prev - du_min(I) < 1.0e-15*abs(du(I,j))) do_I(I) = .false.
            endif
          endif
        else
          do_I(I) = .false.
        endif
      enddo

      ! Below conditional compilation is to control whether early exit happens when compiled with
      ! OpenMP - compiling with OpenMP prevents early exit. Without OpenMP, enables early exit.
      ! Early exit saves time on CPU, but causes other loops to be serialized on GPU.
      !$ if (.false.) then
      if (.not. any(do_I(ish-1:ieh))) exit
      !$ endif

      if ((itt < max_itts) .or. uh_3d_a%associated()) then
        do concurrent (I=ish-1:ieh)
          uh_err(I) = 0.0 ; duhdu_tot(I) = 0.0
          if (use_uhbt) uh_err(I) = -uhbt(I,j)
        enddo
        do k=1,nz ; do concurrent (I=ish-1:ieh, do_I(I)) DO_LOCALITY(local(u_new, duhdu))
          u_new = u(I,j,k) + du(I,j) * visc_rem(I,j,k)
          call flux_elem(u_new, h_in(I,j,k), h_in(I+1,j,k), h_W(I,j,k), h_W(I+1,j,k), h_E(I,j,k), &
                         h_E(I+1,j,k), uh_aux(I,k), duhdu, visc_rem(I,j,k), dy_Cu(I,j), &
                         IareaT(I,j), IareaT(I+1,j), IdxT(I,j), IdxT(i+1,j), dt, &
                         vol_CFL, por_face_areaU(I,j,k))
          ! Below if statement looks expensive in profiling results, but I believe it's
          ! masking the expensive update of uh_err beneath
          if (local_OBC) &
            call flux_elem_OBC(u_new, h_in(I,j,k), h_in(I+1,j,k), uh_aux(I,k), duhdu, &
                               visc_rem(I,j,k), por_face_areaU(I,j,k), dy_Cu(I,j), OBC, &
                               OBC%segnum_u(I,j))
          uh_err(I) = uh_err(I) + uh_aux(I,k)
          duhdu_tot(I) = duhdu_tot(I) + duhdu
        enddo ; enddo
        do concurrent (I=ish-1:ieh)
          uh_err_best(I) = min(uh_err_best(I), abs(uh_err(I)))
        enddo
      endif

    enddo ! itt-loop
    if (uh_3d_a%associated()) then
      do concurrent (k=1:nz, I=ish-1:ieh)
        uh_3d(I,j,k) = uh_aux(I,k)
      enddo
    endif
  enddo ! j-loop
  ! If there are any faces which have not converged to within the tolerance,
  ! so-be-it, or else use a final upwind correction?
  ! This never seems to happen with 20 iterations as max_itt.

  !$omp target exit data map(release: do_I, du_max, du_min, duhdu_tot, uh_err, uh_err_best, uh_aux)

end subroutine zonal_flux_adjust


!> Sets a structure that describes the zonal barotropic volume or mass fluxes as a
!! function of barotropic flow to agree closely with the sum of the layer's transports.
subroutine set_zonal_BT_cont(bxC, u_a, h_in_a, h_W_a, h_E_a, BT_cont, du0_a, uh_tot_0_a, &
                             duhdu_tot_0_a, du_max_CFL_a, du_min_CFL_a, dt, vol_CFL, &
                             visc_rem_a, visc_rem_max_a, do_I_a, por_face_areaU_a, &
                             dxCu_a, dy_Cu_a, IareaT_a, IdxT_a)
  type(box_t),             intent(in) :: bxC  !< Iteration box for continuity solver
  type(RealArray_t),       intent(in)    :: u_a    !< Zonal velocity [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in)    :: h_in_a !< Layer thickness used to calculate fluxes
                                                   !! [H ~> m or kg m-2].
  type(RealArray_t),       intent(in)    :: h_W_a !< West edge thickness in the reconstruction
                                                  !! [H ~> m or kg m-2].
  type(RealArray_t),       intent(in)    :: h_E_a !< East edge thickness in the reconstruction
                                                  !! [H ~> m or kg m-2].
  type(BT_cont_type),   intent(inout) :: BT_cont !< A structure with elements
                       !! that describe the effective open face areas as a function of barotropic flow.
  type(RealArray_t),       intent(in)    :: du0_a !< The barotropic velocity increment that gives 0
                                                  !! transport [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in)    :: uh_tot_0_a !< The summed transport with 0 adjustment
                                                       !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(RealArray_t),       intent(in)    :: duhdu_tot_0_a !< The partial derivative of du_err with
                                                          !! du at 0 adjustment
                                                          !! [H L ~> m2 or kg m-1].
  type(RealArray_t),       intent(in)    :: du_max_CFL_a !< Maximum acceptable value of du
                                                         !! [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in)    :: du_min_CFL_a !< Minimum acceptable value of du
                                                         !! [L T-1 ~> m s-1].
  real,                    intent(in) :: dt   !< Time increment [T ~> s].
  logical,                 intent(in)    :: vol_CFL !< If true, use the ratio of the open face
                                                     !! lengths to the tracer cell areas when
                                                     !! estimating CFL numbers. Without
                                                     !! aggress_adjust, the default is false; it is
                                                     !! always true with.
  type(RealArray_t),       intent(in)    :: visc_rem_a !< Both the fraction of the momentum
                                                       !! originally in a layer that remains after a
                                                       !! time-step of viscosity, and the fraction
                                                       !! of a time-step's worth of a barotropic
                                                       !! acceleration that a layer experiences
                                                       !! after viscosity is applied [nondim].
                                                       !! Visc_rem is between 0 (at the bottom) and
                                                       !! 1 (far above the bottom).
  type(RealArray_t),       intent(in)    :: visc_rem_max_a !< Maximum allowable visc_rem [nondim].
  type(LogicalArray_t),    intent(in) :: do_I_a  !< A logical flag indicating which I
                                                 !! values to work on.
  type(RealArray_t),       intent(in)    :: por_face_areaU_a !< fractional open area of U-faces
                                                             !! [nondim]
  type(RealArray_t),       intent(in)    :: dxCu_a !< The dx spacing at u points [L ~> m].
  type(RealArray_t),       intent(in)    :: dy_Cu_a !< The grid cell's unblocked lengths of the
                                                    !! u/v-faces of the h-cell [L ~> m].
  type(RealArray_t),       intent(in)    :: IareaT_a !< The grid cell's 1/areaT [L-2 ~> m-2].
  type(RealArray_t),       intent(in)    :: IdxT_a !< The grid cell's 1/dxT [L-1 ~> m-1].
  ! Local variables
  real, dimension(u_a%lb(1):u_a%ub(1)) :: &
    duL, duR, &       ! The barotropic velocity increments that give the westerly
    du_CFL, &         ! The velocity increment that corresponds to CFL_min [L T-1 ~> m s-1].
                      ! (duL) and easterly (duR) test velocities [L T-1 ~> m s-1].
    FAmt_L, FAmt_R, & ! The summed effective marginal face areas for the 3
    FAmt_0, &         ! test velocities [H L ~> m2 or kg m-1].
    uhtot_L, &        ! The summed transport with the westerly (uhtot_L) and
    uhtot_R           ! and easterly (uhtot_R) test velocities [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: &
    u_L, u_R, &   ! The westerly (u_L), easterly (u_R), and zero-barotropic
    u_0, &        ! transport (u_0) layer test velocities [L T-1 ~> m s-1].
    duhdu_L, &    ! The effective layer marginal face areas with the westerly
    duhdu_R, &    ! (_L), easterly (_R), and zero-barotropic (_0) test
    duhdu_0, &    ! velocities [H L ~> m2 or kg m-1].
    uh_L, uh_R, & ! The layer transports with the westerly (_L), easterly (_R),
    uh_0       ! and zero-barotropic (_0) test velocities [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: FA_0    ! The effective face area with 0 barotropic transport [L H ~> m2 or kg m-1].
  real :: FA_avg  ! The average effective face area [L H ~> m2 or kg m-1], nominally given by
                  ! the realized transport divided by the barotropic velocity.
  real :: visc_rem_lim ! The larger of visc_rem and min_visc_rem [nondim]. This
                       ! limiting is necessary to keep the inverse of visc_rem
                       ! from leading to large CFL numbers.
  real :: min_visc_rem ! The smallest permitted value for visc_rem that is used
                       ! in finding the barotropic velocity that changes the
                       ! flow direction [nondim].  This is necessary to keep the inverse
                       ! of visc_rem from leading to large CFL numbers.
  real :: CFL_min ! A minimal increment in the CFL to try to ensure that the
                  ! flow is truly upwind [nondim]
  real :: Idt     ! The inverse of the time step [T-1 ~> s-1].
  integer :: i, j, k, nz
  integer :: ish      !< Start of i index range.
  integer :: ieh      !< End of i index range.
  integer :: jsh      !< Start of j index range.
  integer :: jeh      !< End of j index range.
  real, dimension(:,:,:), contiguous, pointer :: u, h_in, h_W, h_E, visc_rem, por_face_areaU
  real, dimension(:,:),   contiguous, pointer :: du0, uh_tot_0, duhdu_tot_0, du_max_CFL, &
                                                 du_min_CFL, visc_rem_max
  real, dimension(:,:),   contiguous, pointer :: dxCu, dy_Cu, IareaT, IdxT
  logical, dimension(:,:), contiguous, pointer :: do_I

  call u_a%view(u)
  call h_in_a%view(h_in)
  call h_W_a%view(h_W)
  call h_E_a%view(h_E)
  call du0_a%view(du0)
  call uh_tot_0_a%view(uh_tot_0)
  call duhdu_tot_0_a%view(duhdu_tot_0)
  call du_max_CFL_a%view(du_max_CFL)
  call du_min_CFL_a%view(du_min_CFL)
  call visc_rem_a%view(visc_rem)
  call visc_rem_max_a%view(visc_rem_max)
  call por_face_areaU_a%view(por_face_areaU)
  call dxCu_a%view(dxCu)
  call dy_Cu_a%view(dy_Cu)
  call IareaT_a%view(IareaT)
  call IdxT_a%view(IdxT)
  call do_I_a%view(do_I)

  ish = bxC%idxS(1) ; ieh = bxC%idxE(1) ; jsh = bxC%idxS(2) ; jeh = bxC%idxE(2) ; nz  = bxC%idxE(3)
  Idt = 1.0 / dt
  min_visc_rem = 0.1 ; CFL_min = 1e-6

  !$omp target enter data map(alloc: duL, duR, du_CFL, FAmt_L, FAmT_R, FAmt_0, uhtot_L, uhtot_R)

  !$omp target teams loop private(duL, duR, du_CFL, FAmt_L, FAmt_R, FAmt_0, uhtot_L, uhtot_R)
  do j=jsh,jeh
    ! Determine the westerly- and easterly- fluxes.  Choose a sufficiently
    ! negative velocity correction for the easterly-flux, and a sufficiently
    ! positive correction for the westerly-flux.
    do concurrent (I=ish-1:ieh)
      du_CFL(I) = (CFL_min * Idt) * dxCu(I,j)
      duR(I) = min(0.0,du0(I,j) - du_CFL(I))
      duL(I) = max(0.0,du0(I,j) + du_CFL(I))
      FAmt_L(I) = 0.0 ; FAmt_R(I) = 0.0 ; FAmt_0(I) = 0.0
      uhtot_L(I) = 0.0 ; uhtot_R(I) = 0.0
    enddo

    do k=1,nz ; do concurrent (I=ish-1:ieh, do_I(I,j)) DO_LOCALITY(local(visc_rem_lim))
      visc_rem_lim = max(visc_rem(I,j,k), min_visc_rem*visc_rem_max(I,j))
      if (visc_rem_lim > 0.0) then ! This is almost always true for ocean points.
        if (u(I,j,k) + duR(I)*visc_rem_lim > -du_CFL(I)*visc_rem(I,j,k)) &
          duR(I) = -(u(I,j,k) + du_CFL(I)*visc_rem(I,j,k)) / visc_rem_lim
        if (u(I,j,k) + duL(I)*visc_rem_lim < du_CFL(I)*visc_rem(I,j,k)) &
          duL(I) = -(u(I,j,k) - du_CFL(I)*visc_rem(I,j,k)) / visc_rem_lim
      endif
    enddo ; enddo

    do k=1,nz ; do concurrent (I=ish-1:ieh, do_I(I,j)) &
        & DO_LOCALITY(local(u_0, u_L, u_R, uh_0, uh_L, uh_R, duhdu_0, duhdu_L, duhdu_R))
      u_L = u(I,j,k) + duL(I) * visc_rem(I,j,k)
      u_R = u(I,j,k) + duR(I) * visc_rem(I,j,k)
      u_0 = u(I,j,k) + du0(I,j) * visc_rem(I,j,k)
      call flux_elem(u_0, h_in(I,j,k), h_in(I+1,j,k), h_W(I,j,k), h_W(I+1,j,k), h_E(I,j,k), &
                     h_E(I+1,j,k), uh_0, duhdu_0, visc_rem(I,j,k), dy_Cu(I,j), &
                     IareaT(I,j), IareaT(I+1,j), IdxT(I,j), IdxT(i+1,j), dt, &
                     vol_CFL, por_face_areaU(I,j,k))
      call flux_elem(u_L, h_in(I,j,k), h_in(I+1,j,k), h_W(I,j,k), h_W(I+1,j,k), h_E(I,j,k), &
                     h_E(I+1,j,k), uh_L, duhdu_L, visc_rem(I,j,k), dy_Cu(I,j), &
                     IareaT(I,j), IareaT(I+1,j), IdxT(I,j), IdxT(i+1,j), dt, &
                     vol_CFL, por_face_areaU(I,j,k))
      call flux_elem(u_R, h_in(I,j,k), h_in(I+1,j,k), h_W(I,j,k), h_W(I+1,j,k), h_E(I,j,k), &
                     h_E(I+1,j,k), uh_R, duhdu_R, visc_rem(I,j,k), dy_Cu(I,j), &
                     IareaT(I,j), IareaT(I+1,j), IdxT(I,j), IdxT(i+1,j), dt, &
                     vol_CFL, por_face_areaU(I,j,k))
      FAmt_0(I) = FAmt_0(I) + duhdu_0
      FAmt_L(I) = FAmt_L(I) + duhdu_L
      FAmt_R(I) = FAmt_R(I) + duhdu_R
      uhtot_L(I) = uhtot_L(I) + uh_L
      uhtot_R(I) = uhtot_R(I) + uh_R
    enddo ; enddo

    do concurrent (I=ish-1:ieh) DO_LOCALITY(local(FA_0, FA_avg))
      if (do_I(I,j)) then
        FA_0 = FAmt_0(I) ; FA_avg = FAmt_0(I)
        if ((duL(I) - du0(I,j)) /= 0.0) &
          FA_avg = uhtot_L(I) / (duL(I) - du0(I,j))
        if (FA_avg > max(FA_0, FAmt_L(I))) then ; FA_avg = max(FA_0, FAmt_L(I))
        elseif (FA_avg < min(FA_0, FAmt_L(I))) then ; FA_0 = FA_avg ; endif

        BT_cont%FA_u_W0(I,j) = FA_0 ; BT_cont%FA_u_WW(I,j) = FAmt_L(I)
        if (abs(FA_0-FAmt_L(I)) <= 1e-12*FA_0) then ; BT_cont%uBT_WW(I,j) = 0.0 ; else
          BT_cont%uBT_WW(I,j) = (1.5 * (duL(I) - du0(I,j))) * &
                                ((FAmt_L(I) - FA_avg) / (FAmt_L(I) - FA_0))
        endif

        FA_0 = FAmt_0(I) ; FA_avg = FAmt_0(I)
        if ((duR(I) - du0(I,j)) /= 0.0) &
          FA_avg = uhtot_R(I) / (duR(I) - du0(I,j))
        if (FA_avg > max(FA_0, FAmt_R(I))) then ; FA_avg = max(FA_0, FAmt_R(I))
        elseif (FA_avg < min(FA_0, FAmt_R(I))) then ; FA_0 = FA_avg ; endif

        BT_cont%FA_u_E0(I,j) = FA_0 ; BT_cont%FA_u_EE(I,j) = FAmt_R(I)
        if (abs(FAmt_R(I) - FA_0) <= 1e-12*FA_0) then ; BT_cont%uBT_EE(I,j) = 0.0 ; else
          BT_cont%uBT_EE(I,j) = (1.5 * (duR(I) - du0(I,j))) * &
                                ((FAmt_R(I) - FA_avg) / (FAmt_R(I) - FA_0))
        endif
      else
        BT_cont%FA_u_W0(I,j) = 0.0 ; BT_cont%FA_u_WW(I,j) = 0.0
        BT_cont%FA_u_E0(I,j) = 0.0 ; BT_cont%FA_u_EE(I,j) = 0.0
        BT_cont%uBT_WW(I,j) = 0.0 ; BT_cont%uBT_EE(I,j) = 0.0
      endif
    enddo
  enddo

  !$omp target exit data map(release: duL, duR, du_CFL, FAmt_L, FAmT_R, FAmt_0, uhtot_L, uhtot_R)

end subroutine set_zonal_BT_cont

!> Calculates the mass or volume fluxes through the meridional faces, and other related quantities.
subroutine meridional_mass_flux(bxC, v_a, h_in_a, h_S_a, h_N_a, vh_a, dt, OBC, &
                                 por_face_areaV_a, vhbt_a, visc_rem_v_a, v_cor_a, BT_cont, &
                                 dv_cor_a, dx_Cv_a, IareaT_a, IdyT_a, dyCv_a, areaT_a, dyT_a, &
                                 mask2dCv_a, H_subroundoff, CFL_limit_adjust, aggress_adjust, &
                                 use_visc_rem_max, vol_CFL, tol_vel, tol_eta, better_iter, &
                                 marginal_faces)
  type(Box_t),             intent(in)    :: bxC  !< Iteration box for continuity solver
  type(RealArray_t),       intent(in)    :: v_a    !< Meridional velocity [L T-1 ~> m s-1]
  type(RealArray_t),       intent(in)    :: h_in_a !< Layer thickness used to calculate fluxes
                                                   !! [H ~> m or kg m-2]
  type(RealArray_t),       intent(in)    :: h_S_a !< South edge thickness in the reconstruction
                                                  !! [H ~> m or kg m-2].
  type(RealArray_t),       intent(in)    :: h_N_a !< North edge thickness in the reconstruction
                                                  !! [H ~> m or kg m-2].
  type(RealArray_t),       intent(inout) :: vh_a   !< Volume flux through meridional faces = v*h*dx
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real,                    intent(in)    :: dt   !< Time increment [T ~> s].
  type(ocean_OBC_type),    pointer       :: OBC  !< Open boundary condition type
                                                 !! specifies whether, where, and what
                                                 !! open boundary conditions are used.
  type(RealArray_t),       intent(in)    :: por_face_areaV_a !< fractional open area of V-faces
                                                             !! [nondim]
  type(RealArray_t), optional, intent(in) :: vhbt_a !< The summed volume flux through meridional
                                                 !! faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(RealArray_t), optional, intent(in) :: visc_rem_v_a
                     !< Both the fraction of the momentum originally in a layer that remains after a
                     !! time-step of viscosity, and the fraction of a time-step's worth of a
                     !! barotropic acceleration that a layer experiences after viscosity is applied
                     !! [nondim]. Visc_rem_v is between 0 (at the bottom) and 1 (far above the
                     !! bottom).
  type(RealArray_t), optional, intent(inout) :: v_cor_a
                     !< The meridional velocities (v with a barotropic correction)
                     !! that give vhbt as the depth-integrated transport [L T-1 ~> m s-1].
  type(BT_cont_type), optional, pointer  :: BT_cont !< A structure with elements that describe the
                     !! effective open face areas as a function of barotropic flow.
  type(RealArray_t), optional, intent(inout) :: dv_cor_a
                                     !< The meridional velocity increments from v that give vhbt
                                     !! as the depth-integrated transports [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in)    :: dx_Cv_a !< The grid cell's unblocked lengths of the
                                                    !! u/v-faces of the h-cell [L ~> m].
  type(RealArray_t),       intent(in)    :: IareaT_a !< The grid cell's 1/areaT [L-2 ~> m-2].
  type(RealArray_t),       intent(in)    :: IdyT_a !< The grid cell's 1/dyT [L-1 ~> m-1].
  type(RealArray_t),       intent(in)    :: dyCv_a !< The dy spacing at v points [L ~> m].
  type(RealArray_t),       intent(in)    :: areaT_a !< The grid cell's area [L2 ~> m2].
  type(RealArray_t),       intent(in)    :: dyT_a !< The dy spacing at h points [L ~> m].
  type(RealArray_t),       intent(in)    :: mask2dCv_a !< 0 for boundary points and 1 for ocean
                                                       !! points on the v grid [nondim].
  real,                    intent(in)    :: H_subroundoff !< A thickness that is so small that it
                     !! can be added to a thickness of Angstrom or larger without changing it at the
                     !! bit level [H ~> m or kg m-2].
  real,                    intent(in)    :: CFL_limit_adjust !< The maximum CFL of the adjusted
                                                             !! velocities [nondim]
  logical,                 intent(in)    :: aggress_adjust !< If true, allow the adjusted velocities
                                                           !! to have a relative CFL change up to
                                                           !! 0.5. False by default.
  logical,                 intent(in)    :: use_visc_rem_max !< If true, use more appropriate
                                                             !! limiting bounds for corrections in
                                                             !! strongly viscous columns.
  logical,                 intent(in)    :: vol_CFL !< If true, use the ratio of the open face
                                                     !! lengths to the tracer cell areas when
                                                     !! estimating CFL numbers. Without
                                                     !! aggress_adjust, the default is false; it is
                                                     !! always true with.
  real,                    intent(in)    :: tol_vel !< The tolerance for barotropic velocity
                                                     !! discrepancies between the barotropic
                                                     !! solution and the sum of the layer
                                                     !! thicknesses [L T-1 ~> m s-1].
  real,                    intent(in)    :: tol_eta !< The tolerance for free-surface height
                                                     !! discrepancies between the barotropic
                                                     !! solution and the sum of the layer
                                                     !! thicknesses [H ~> m or kg m-2].
  logical,                 intent(in)    :: better_iter !< If true, stop corrective iterations
                                                         !! using a velocity-based criterion and
                                                         !! only stop if the iteration is better
                                                         !! than all predecessors.
  logical,                 intent(in)    :: marginal_faces !< If true, use the marginal face areas
                          !! from the continuity solver for use as the weights in the barotropic
                          !! solver. Otherwise use the transport averaged areas.

  ! Local variables
  real, dimension(v_a%lb(1):v_a%ub(1), v_a%lb(2):v_a%ub(2), v_a%lb(3):v_a%ub(3)) :: &
    dvhdv         ! Partial derivative of vh with v [H L ~> m2 or kg m-1].
  real :: I_vrm   ! 1.0 / visc_rem_max [nondim]
  real :: CFL_dt  ! The maximum CFL ratio of the adjusted velocities divided by
                  ! the time step [T-1 ~> s-1].
  real :: I_dt    ! 1.0 / dt [T-1 ~> s-1].
  real :: dv_lim  ! The velocity change that give a relative CFL of 1 [L T-1 ~> m s-1].
  real :: dy_N, dy_S ! Effective y-grid spacings to the north and south [L ~> m].
  integer :: i, j, k, ish, ieh, jsh, jeh, n, nz
  integer :: l_seg ! The OBC segment number
  logical :: use_visc_rem, set_BT_cont
  logical :: local_specified_BC, local_open_BC, any_simple_OBC  ! OBC-related logicals
  type(RealArray_t) :: vh_tot_0_a, dvhdv_tot_0_a, dv_a, dv_max_CFL_a, dv_min_CFL_a
  type(RealArray_t) :: visc_rem_max_a, visc_rem_v_tmp_a
  real, dimension(:,:,:), contiguous, pointer :: v, h_in, h_S, h_N, vh, por_face_areaV
  real, dimension(:,:,:), contiguous, pointer :: visc_rem_v
  real, dimension(:,:),   contiguous, pointer :: dv_cor
  real, dimension(:,:),   contiguous, pointer :: dx_Cv, IareaT, IdyT, areaT, dyT, mask2dCv
  real, dimension(:,:),   contiguous, pointer :: &
    dv, &         ! Corrective barotropic change in the velocity to give vhbt [L T-1 ~> m s-1].
    dv_min_CFL, & ! Lower limit on dv correction to avoid CFL violations [L T-1 ~> m s-1]
    dv_max_CFL, & ! Upper limit on dv correction to avoid CFL violations [L T-1 ~> m s-1]
    dvhdv_tot_0, & ! Summed partial derivative of vh with v [H L ~> m2 or kg m-1].
    vh_tot_0, &   ! Summed transport with no barotropic correction [H L2 T-1 ~> m3 s-1 or kg s-1].
    visc_rem_max  ! The column maximum of visc_rem [nondim]
  real, dimension(:,:,:), contiguous, pointer :: &
    visc_rem_v_tmp ! A copy of visc_rem_v or an array of 1's [nondim]

  nullify(visc_rem_v, dv_cor)

  call cpu_clock_begin(id_clock_correct)

  call v_a%view(v)
  call h_in_a%view(h_in)
  call h_S_a%view(h_S)
  call h_N_a%view(h_N)
  call vh_a%view(vh)
  call por_face_areaV_a%view(por_face_areaV)
  if (present(visc_rem_v_a)) call visc_rem_v_a%view(visc_rem_v)
  if (present(dv_cor_a)) call dv_cor_a%view(dv_cor)
  call dx_Cv_a%view(dx_Cv)
  call IareaT_a%view(IareaT)
  call IdyT_a%view(IdyT)
  call areaT_a%view(areaT)
  call dyT_a%view(dyT)
  call mask2dCv_a%view(mask2dCv)

  use_visc_rem = present(visc_rem_v_a)

  set_BT_cont = .false. ; if (present(BT_cont)) set_BT_cont = (associated(BT_cont))

  local_specified_BC = .false. ; local_open_BC = .false.
  if (associated(OBC)) then ; if (OBC%OBC_pe) then
    local_specified_BC = OBC%specified_v_BCs_exist_globally
    local_open_BC = OBC%open_v_BCs_exist_globally
  endif ; endif

  ish = bxC%idxS(1) ; ieh = bxC%idxE(1) ; jsh = bxC%idxS(2) ; jeh = bxC%idxE(2) ; nz  = bxC%idxE(3)

  CFL_dt = CFL_limit_adjust / dt
  I_dt = 1.0 / dt
  if (aggress_adjust) CFL_dt = I_dt

  call dv_a%alloc(lb=[v_a%lb(1), v_a%lb(2)], ub=[v_a%ub(1), v_a%ub(2)])
  call dv_min_CFL_a%alloc(lb=[v_a%lb(1), v_a%lb(2)], ub=[v_a%ub(1), v_a%ub(2)])
  call dv_max_CFL_a%alloc(lb=[v_a%lb(1), v_a%lb(2)], ub=[v_a%ub(1), v_a%ub(2)])
  call dvhdv_tot_0_a%alloc(lb=[v_a%lb(1), v_a%lb(2)], ub=[v_a%ub(1), v_a%ub(2)])
  call vh_tot_0_a%alloc(lb=[v_a%lb(1), v_a%lb(2)], ub=[v_a%ub(1), v_a%ub(2)])
  call visc_rem_max_a%alloc(lb=[v_a%lb(1), v_a%lb(2)], ub=[v_a%ub(1), v_a%ub(2)])
  call visc_rem_v_tmp_a%alloc(lb=[v_a%lb(1), v_a%lb(2), 1], ub=[v_a%ub(1), v_a%ub(2), v_a%ub(3)])
  call dv_a%view(dv)
  call dv_min_CFL_a%view(dv_min_CFL)
  call dv_max_CFL_a%view(dv_max_CFL)
  call dvhdv_tot_0_a%view(dvhdv_tot_0)
  call vh_tot_0_a%view(vh_tot_0)
  call visc_rem_max_a%view(visc_rem_max)
  call visc_rem_v_tmp_a%view(visc_rem_v_tmp)

  !$omp target enter data &
  !$omp   map(alloc: dvhdv)

  do concurrent (J=jsh-1:jeh)

    if (present(dv_cor_a)) then
      do concurrent (i=ish:ieh)
        dv_cor(i,J) = 0.0
      enddo
    endif

    ! visc_rem_v_tmp must be valid over the full local domain, not just i=ish:ieh --
    ! meridional_flux_thickness's open-boundary-segment branch reads it beyond this box
    ! when this PE owns a segment of an open boundary.
    if (.not.use_visc_rem) then
      do concurrent (k=1:nz, i=LBOUND(visc_rem_v_tmp,1):UBOUND(visc_rem_v_tmp,1))
        visc_rem_v_tmp(i,J,k) = 1.0
      enddo
    else
      do concurrent (k=1:nz, i=LBOUND(visc_rem_v_tmp,1):UBOUND(visc_rem_v_tmp,1))
        visc_rem_v_tmp(i,J,k) = visc_rem_v(i,J,k)
      enddo
    endif
    do concurrent (k=1:nz, i=ish:ieh)
      call flux_elem(v(i,J,k), h_in(i,J,k), h_in(i,J+1,k), h_S(i,J,k), h_S(i,J+1,k), h_N(i,J,k), &
                     h_N(i,J+1,k), vh(i,J,k), dvhdv(i,J,k), visc_rem_v_tmp(i,J,k), dx_Cv(i,J), &
                     IareaT(i,J), IareaT(i,J+1), IdyT(i,J), IdyT(i,J+1), dt, &
                     vol_CFL, por_face_areaV(i,J,k))
    enddo
    if (local_open_BC) then
      do concurrent (k=1:nz, i=ish:ieh)
        ! untested!
        call flux_elem_OBC(v(i,J,k), h_in(i,J,k), h_in(i,J+1,k), vh(i,J,k), dvhdv(i,J,k), &
                           visc_rem_v_tmp(i,J,k), por_face_areaV(i,J,k), &
                           dx_Cv(i,J), OBC, OBC%segnum_v(i,J))
      enddo
    endif

    ! untested!
    if (local_specified_BC) then
      do concurrent (k=1:nz, i=ish:ieh, OBC%segnum_v(i,J) /= 0)
        l_seg = abs(OBC%segnum_v(i,J))
        if (OBC%segment(l_seg)%specified) vh(i,J,k) = OBC%segment(l_seg)%normal_trans(i,J,k)
      enddo
    endif

    if (present(vhbt_a) .or. set_BT_cont) then
      if (use_visc_rem .and. use_visc_rem_max) then
        do concurrent (i=ish:ieh)
          visc_rem_max(i,J) = 0.0
        enddo
        do k=1,nz ; do concurrent (i=ish:ieh)
          visc_rem_max(i,J) = max(visc_rem_max(i,J), visc_rem_v_tmp(i,J,k))
        enddo ; enddo
      else
        do concurrent (i=ish:ieh)
          visc_rem_max(i,J) = 1.0
        enddo
      endif
      !   Set limits on dv that will keep the CFL number between -1 and 1.
      ! This should be adequate to keep the root bracketed in all cases.
      do concurrent (i=ish:ieh)
        I_vrm = 0.0
        if (visc_rem_max(i,j) > 0.0) I_vrm = 1.0 / visc_rem_max(i,j)
        if (vol_CFL) then
          dy_S = ratio_max(areaT(i,j), dx_Cv(i,J), 1000.0*dyT(i,j))
          dy_N = ratio_max(areaT(i,j+1), dx_Cv(i,J), 1000.0*dyT(i,j+1))
        else ; dy_S = dyT(i,j) ; dy_N = dyT(i,j+1) ; endif
        dv_max_CFL(i,j) = 2.0 * (CFL_dt * dy_S) * I_vrm
        dv_min_CFL(i,j) = -2.0 * (CFL_dt * dy_N) * I_vrm
        vh_tot_0(i,j) = 0.0 ; dvhdv_tot_0(i,j) = 0.0
      enddo
      do k=1,nz ; do concurrent (i=ish:ieh)
        dvhdv_tot_0(i,j) = dvhdv_tot_0(i,j) + dvhdv(i,j,k)
        vh_tot_0(i,j) = vh_tot_0(i,j) + vh(i,J,k)
      enddo ; enddo


      if (use_visc_rem) then
        if (aggress_adjust) then
          ! untested
          do k=1,nz ; do concurrent (i=ish:ieh)
            if (vol_CFL) then
              dy_S = ratio_max(areaT(i,J), dx_Cv(i,J), 1000.0*dyT(i,J))
              dy_N = ratio_max(areaT(i,J+1), dx_Cv(i,J), 1000.0*dyT(i,J+1))
            else ; dy_S = dyT(i,J) ; dy_N = dyT(i,J+1) ; endif
            dv_lim = 0.499*((dy_S*I_dt - v(i,J,k)) + MIN(0.0,v(i,J-1,k)))
            if (dv_max_CFL(i,J) * visc_rem_v_tmp(i,J,k) > dv_lim) &
              dv_max_CFL(i,J) = dv_lim / visc_rem_v_tmp(i,J,k)

            dv_lim = 0.499*((-dy_N*CFL_dt - v(i,J,k)) + MAX(0.0,v(i,J+1,k)))
            if (dv_min_CFL(i,J) * visc_rem_v_tmp(i,J,k) < dv_lim) &
              dv_min_CFL(i,J) = dv_lim / visc_rem_v_tmp(i,J,k)
          enddo ; enddo
        else
          do k=1,nz ; do concurrent (i=ish:ieh)
            if (vol_CFL) then
              dy_S = ratio_max(areaT(i,J), dx_Cv(i,J), 1000.0*dyT(i,J))
              dy_N = ratio_max(areaT(i,J+1), dx_Cv(i,J), 1000.0*dyT(i,J+1))
            else ; dy_S = dyT(i,J) ; dy_N = dyT(i,J+1) ; endif
            if (dv_max_CFL(i,J) * visc_rem_v_tmp(i,J,k) > dy_S*CFL_dt - v(i,J,k)*mask2dCv(i,J)) &
              dv_max_CFL(i,J) = (dy_S*CFL_dt - v(i,J,k)) / visc_rem_v_tmp(i,J,k)
            if (dv_min_CFL(i,J) * visc_rem_v_tmp(i,J,k) < -dy_N*CFL_dt - v(i,J,k)*mask2dCv(i,J)) &
              dv_min_CFL(i,J) = -(dy_N*CFL_dt + v(i,J,k)) / visc_rem_v_tmp(i,J,k)
          enddo ; enddo
        endif ! CS%agress_adjust
      else
        if (aggress_adjust) then
          ! untested
          do k=1,nz ; do concurrent (i=ish:ieh)
            if (vol_CFL) then
              dy_S = ratio_max(areaT(i,J), dx_Cv(i,J), 1000.0*dyT(i,J))
              dy_N = ratio_max(areaT(i,J+1), dx_Cv(i,J), 1000.0*dyT(i,J+1))
            else ; dy_S = dyT(i,J) ; dy_N = dyT(i,J+1) ; endif
            dv_max_CFL(i,J) = min(dv_max_CFL(i,J), 0.499 * &
                        ((dy_S*I_dt - v(i,J,k)) + MIN(0.0,v(i,J-1,k))) )
            dv_min_CFL(i,J) = max(dv_min_CFL(i,J), 0.499 * &
                        ((-dy_N*I_dt - v(i,J,k)) + MAX(0.0,v(i,J+1,k))) )
          enddo ; enddo
        else
          do k=1,nz ; do concurrent (i=ish:ieh)
            if (vol_CFL) then
              dy_S = ratio_max(areaT(i,J), dx_Cv(i,J), 1000.0*dyT(i,J))
              dy_N = ratio_max(areaT(i,J+1), dx_Cv(i,J), 1000.0*dyT(i,J+1))
            else ; dy_S = dyT(i,J) ; dy_N = dyT(i,J+1) ; endif
            dv_max_CFL(i,J) = min(dv_max_CFL(i,J), dy_S*CFL_dt - v(i,J,k))
            dv_min_CFL(i,J) = max(dv_min_CFL(i,J), -(dy_N*CFL_dt + v(i,J,k)))
          enddo ; enddo
        endif ! CS%agress_adjust
      endif ! use_visc_rem
      do concurrent (i=ish:ieh)
        dv_max_CFL(i,J) = max(dv_max_CFL(i,J),0.0)
        dv_min_CFL(i,J) = min(dv_min_CFL(i,J),0.0)
      enddo
    endif ! present(vhbt_a) .or. set_BT_cont

  enddo

  call present_vhbt_or_set_BT_cont(bxC, v_a, h_in_a, h_S_a, h_N_a, vh_tot_0_a, dvhdv_tot_0_a, &
                                   dv_a, dv_max_CFL_a, dv_min_CFL_a, visc_rem_v_tmp_a, &
                                   visc_rem_max_a, por_face_areaV_a, vhbt_a, vh_a, v_cor_a, &
                                   dv_cor_a, BT_cont, dt, OBC, dx_Cv_a, H_subroundoff, &
                                   IareaT_a, IdyT_a, dyCv_a, tol_vel, tol_eta, &
                                   better_iter, vol_CFL, marginal_faces)
  call vh_tot_0_a%free()
  call dvhdv_tot_0_a%free()
  call dv_a%free()
  call dv_max_CFL_a%free()
  call dv_min_CFL_a%free()
  call visc_rem_max_a%free()
  call visc_rem_v_tmp_a%free()

  !$omp target exit data &
  !$omp   map(release: dvhdv)

  call cpu_clock_end(id_clock_correct)

end subroutine meridional_mass_flux

subroutine present_vhbt_or_set_BT_cont(bxC, v_a, h_in_a, h_S_a, h_N_a, vh_tot_0_a, &
                                       dvhdv_tot_0_a, dv_a, dv_max_CFL_a, dv_min_CFL_a, &
                                       visc_rem_v_a, visc_rem_max_a, por_face_areaV_a, vhbt_a, &
                                       vh_a, v_cor_a, dv_cor_a, BT_cont, dt, OBC, &
                                       dx_Cv_a, H_subroundoff, IareaT_a, IdyT_a, dyCv_a, &
                                       tol_vel, tol_eta, better_iter, vol_CFL, marginal_faces)
  type(box_t), intent(in) :: bxC                 !< Iteration box for continuity solver
  type(RealArray_t),       intent(in)    :: v_a    !< Meridional velocity [L T-1 ~> m s-1]
  type(RealArray_t),       intent(in)    :: h_in_a !< Layer thickness used to calculate fluxes
                                                   !! [H ~> m or kg m-2]
  type(RealArray_t),       intent(in)    :: h_S_a !< South edge thickness in the reconstruction
                                                  !! [H ~> m or kg m-2].
  type(RealArray_t),       intent(in)    :: h_N_a !< North edge thickness in the reconstruction
                                                  !! [H ~> m or kg m-2].
  type(RealArray_t),       intent(in)    :: vh_tot_0_a !< Summed transport with no barotropic
                                                       !! correction [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(RealArray_t),       intent(in)    :: dvhdv_tot_0_a !< Summed partial derivative of vh with v
                                                          !! [H L ~> m2 or kg m-1].
  type(RealArray_t),       intent(inout) :: dv_a !< Corrective barotropic change in the velocity to
                                                 !! give vhbt [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in)    :: dv_max_CFL_a !< Upper limit on dv correction to avoid
                                                         !! CFL violations [L T-1 ~> m s-1]
  type(RealArray_t),       intent(in)    :: dv_min_CFL_a !< Lower limit on dv correction to avoid
                                                         !! CFL violations [L T-1 ~> m s-1]
  type(RealArray_t),       intent(inout) :: vh_a   !< Volume flux through meridional faces = v*h*dx
                                                   !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  type(RealArray_t),       intent(in)    :: visc_rem_v_a !< Both the fraction of the momentum
                                   !! originally in a layer that remains after a time-step of viscosity,
                                   !! and the fraction of a time-step's worth of a barotropic acceleration
                                   !! that a layer experiences after viscosity is applied [nondim].
                                   !! Visc_rem_v is between 0 (at the bottom) and 1 (far above the bottom).
  type(RealArray_t),       intent(in)    :: visc_rem_max_a !< The column maximum of visc_rem
                                                           !! [nondim]
  type(RealArray_t),       intent(in)    :: por_face_areaV_a !< fractional open area of V-faces
                                                             !! [nondim]
  type(RealArray_t), optional, intent(inout) :: v_cor_a !< The meridional velocities (v with a
                                                         !! barotropic correction) that give vhbt as
                                                         !! the depth-integrated transport
                                                         !! [L T-1 ~> m s-1].
  type(RealArray_t), optional, intent(inout) :: dv_cor_a !< The meridional velocity increments from
                                                         !! v that give vhbt as the depth-integrated
                                                         !! transports [L T-1 ~> m s-1].
  type(BT_cont_type), optional, pointer :: BT_cont !< A structure with elements that describe the
                     !! effective open face areas as a function of barotropic flow.
  type(RealArray_t), optional, intent(in) :: vhbt_a !< The summed volume flux through meridional
                                                    !! faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  real,                    intent(in)    :: dt  !< Time increment [T ~> s].
  type(ocean_OBC_type), pointer :: OBC !< Open boundaries control structure.
  type(RealArray_t),       intent(in)    :: dx_Cv_a !< The grid cell's unblocked lengths of the
                                                    !! u/v-faces of the h-cell [L ~> m].
  real,                    intent(in)    :: H_subroundoff !< A thickness that is so small that it
                     !! can be added to a thickness of Angstrom or larger without changing it at the
                     !! bit level [H ~> m or kg m-2].
  type(RealArray_t),       intent(in)    :: IareaT_a !< The grid cell's 1/areaT [L-2 ~> m-2].
  type(RealArray_t),       intent(in)    :: IdyT_a !< The grid cell's 1/dyT [L-1 ~> m-1].
  type(RealArray_t),       intent(in)    :: dyCv_a !< The dy spacing at v points [L ~> m].
  real,                    intent(in)    :: tol_vel !< The tolerance for barotropic velocity
                                                     !! discrepancies between the barotropic
                                                     !! solution and the sum of the layer
                                                     !! thicknesses [L T-1 ~> m s-1].
  real,                    intent(in)    :: tol_eta !< The tolerance for free-surface height
                                                     !! discrepancies between the barotropic
                                                     !! solution and the sum of the layer
                                                     !! thicknesses [H ~> m or kg m-2].
  logical,                 intent(in)    :: better_iter !< If true, stop corrective iterations
                                                         !! using a velocity-based criterion and
                                                         !! only stop if the iteration is better
                                                         !! than all predecessors.
  logical,                 intent(in)    :: vol_CFL !< If true, use the ratio of the open face
                                                     !! lengths to the tracer cell areas when
                                                     !! estimating CFL numbers. Without
                                                     !! aggress_adjust, the default is false; it is
                                                     !! always true with.
  logical,                 intent(in)    :: marginal_faces !< If true, report the marginal face
                          !! thicknesses; otherwise report transport-averaged thicknesses.
  ! Local variables
  logical, dimension(v_a%lb(1):v_a%ub(1), v_a%lb(2):v_a%ub(2)) :: &
    simple_OBC_pt ! Indicates points in a row with specified transport OBCs
  type(OBC_segment_type), pointer :: segment => NULL()
  real :: FAvi, FA_v    ! A sum of meridional face areas [H L ~> m2 or kg m-1].
  logical :: set_BT_cont
  logical :: any_simple_OBC, local_specified_BC, local_Flather_OBC, local_open_BC  ! OBC-related logicals
  integer :: l_seg, i, j, k, n, ish, ieh, jsh, jeh, nz
  real, dimension(:,:,:), contiguous, pointer :: v, h_in, h_S, h_N, vh, visc_rem_v, por_face_areaV
  real, dimension(:,:,:), contiguous, pointer :: v_cor
  real, dimension(:,:),   contiguous, pointer :: vh_tot_0, dvhdv_tot_0, dv, dv_max_CFL, dv_min_CFL
  real, dimension(:,:),   contiguous, pointer :: visc_rem_max, vhbt, dv_cor
  real, dimension(:,:),   contiguous, pointer :: dx_Cv
  logical, dimension(:,:), contiguous, pointer :: do_I
  type(RealArray_t) :: h_v_a
  type(LogicalArray_t) :: do_I_a
  type(RealArray_t) :: no_vh_3d_a ! Never allocated; unassociated data signals vh_3d_a absent.

  call v_a%view(v)
  call h_in_a%view(h_in)
  call h_S_a%view(h_S)
  call h_N_a%view(h_N)
  call vh_tot_0_a%view(vh_tot_0)
  call dvhdv_tot_0_a%view(dvhdv_tot_0)
  call dv_a%view(dv)
  call dv_max_CFL_a%view(dv_max_CFL)
  call dv_min_CFL_a%view(dv_min_CFL)
  call vh_a%view(vh)
  call visc_rem_v_a%view(visc_rem_v)
  call visc_rem_max_a%view(visc_rem_max)
  call por_face_areaV_a%view(por_face_areaV)
  if (present(vhbt_a)) call vhbt_a%view(vhbt)
  if (present(v_cor_a)) call v_cor_a%view(v_cor)
  if (present(dv_cor_a)) call dv_cor_a%view(dv_cor)
  call dx_Cv_a%view(dx_Cv)

  ish = bxC%idxS(1) ; ieh = bxC%idxE(1) ; jsh = bxC%idxS(2) ; jeh = bxC%idxE(2) ; nz  = bxC%idxE(3)

  set_BT_cont = .false. ; if (present(BT_cont)) set_BT_cont = (associated(BT_cont))

  local_specified_BC = .false. ; local_Flather_OBC = .false. ; local_open_BC = .false.
  if (associated(OBC)) then ; if (OBC%OBC_pe) then
    local_specified_BC = OBC%specified_u_BCs_exist_globally
    local_Flather_OBC = OBC%Flather_u_BCs_exist_globally
    local_open_BC = OBC%open_u_BCs_exist_globally
  endif ; endif

  if (present(vhbt_a) .or. set_BT_cont) then
    call do_I_a%alloc(lb=[v_a%lb(1), v_a%lb(2)], ub=[v_a%ub(1), v_a%ub(2)])
    call do_I_a%view(do_I)
    !$omp target enter data map(alloc: do_I, simple_OBC_pt)
    any_simple_OBC = .false.
    if (local_specified_BC .or. local_Flather_OBC) then
      do concurrent (j=jsh-1:jeh, i=ish:ieh)
        l_seg = abs(OBC%segnum_v(i,J))

        ! Avoid reconciling barotropic/baroclinic transports if transport is specified
        simple_OBC_pt(i,J) = .false.
        if (l_seg /= 0) simple_OBC_pt(i,J) = OBC%segment(l_seg)%specified
        do_I(i,J) = .not.simple_OBC_pt(i,J)
        any_simple_OBC = any_simple_OBC .or. simple_OBC_pt(i,J)
      enddo
    else
      do concurrent (J=jsh-1:jeh, i=ish:ieh)
        do_I(i,J) = .true.
      enddo
    endif ! local_specified_BC .or. local_Flather_OBC

    if (present(vhbt_a)) then
      ! Find dv and vh.
      call meridional_flux_adjust(bxC, v_a, h_in_a, h_S_a, h_N_a, vh_tot_0_a, dvhdv_tot_0_a, &
                             dv_a, dv_max_CFL_a, dv_min_CFL_a, dt, tol_vel, &
                             tol_eta, better_iter, vol_CFL, visc_rem_v_a, do_I_a, &
                             por_face_areaV_a, dx_Cv_a, IareaT_a, IdyT_a, vhbt_a, vh_a, OBC=OBC)

      do concurrent (J=jsh-1:jeh)
        if (present(v_cor_a)) then
          do concurrent (k=1:nz, i=ish:ieh)
            v_cor(i,J,k) = v(i,J,k) + dv(i,J) * visc_rem_v(i,J,k)
          enddo
          if (any_simple_OBC) then
            ! untested
            do concurrent (k=1:nz, i=ish:ieh, simple_OBC_pt(i,J))
              v_cor(i,J,k) = OBC%segment(abs(OBC%segnum_v(i,J)))%normal_vel(i,J,k)
            enddo
          endif
        endif ! v-corrected

        if (present(dv_cor_a)) then
          do concurrent (i=ish:ieh)
            dv_cor(i,J) = dv(i,J)
          enddo
        endif ! dv-corrected
      enddo
    endif

    if (set_BT_cont) then
    ! Diagnose the zero-transport correction, dv0.
      call meridional_flux_adjust(bxC, v_a, h_in_a, h_S_a, h_N_a, vh_tot_0_a, dvhdv_tot_0_a, &
                            dv_a, dv_max_CFL_a, dv_min_CFL_a, dt, tol_vel, &
                            tol_eta, better_iter, vol_CFL, visc_rem_v_a, do_I_a, &
                            por_face_areaV_a, dx_Cv_a, IareaT_a, IdyT_a, vh_3d_a=no_vh_3d_a)
      call set_merid_BT_cont(bxC, v_a, h_in_a, h_S_a, h_N_a, BT_cont, dv_a, vh_tot_0_a, &
                             dvhdv_tot_0_a, dv_max_CFL_a, dv_min_CFL_a, dt, vol_CFL, &
                             visc_rem_v_a, visc_rem_max_a, do_I_a, por_face_areaV_a, &
                             dyCv_a, dx_Cv_a, IareaT_a, IdyT_a)

      if (any_simple_OBC) then
        ! untested
        ! NOTE: simple_OBC_pt(i,j) should prevent access to segment OBC_NONE
        do concurrent (J=jsh-1:jeh, i=ish:jeh, simple_OBC_pt(i,J))
          segment => OBC%segment(abs(OBC%segnum_v(i,J)))
          FAvi = H_subroundoff*dx_Cv(i,J)
          do k=1,nz
            if ((abs(segment%normal_vel(i,J,k)) > 0.0) .and. (segment%specified)) &
              FAvi = FAvi + segment%normal_trans(i,J,k) / segment%normal_vel(i,J,k)
          enddo
          BT_cont%FA_v_S0(i,J) = FAvi ; BT_cont%FA_v_N0(i,J) = FAvi
          BT_cont%FA_v_SS(i,J) = FAvi ; BT_cont%FA_v_NN(i,J) = FAvi
          BT_cont%vBT_SS(i,J) = 0.0 ; BT_cont%vBT_NN(i,J) = 0.0
        enddo
      endif ! any_simple_OBC
    endif ! set_BT_cont
    !$omp target exit data map(release: do_I, simple_OBC_pt)
    call do_I_a%free()
  endif ! present(vhbt_a) or set_BT_cont

  ! untested - probably needs to be refactored to be performant on GPU
  if (local_open_BC .and. set_BT_cont) then
    do n = 1, OBC%number_of_segments
      if (OBC%segment(n)%open .and. OBC%segment(n)%is_N_or_S) then
        J = OBC%segment(n)%HI%JsdB
        if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
          do concurrent (i = OBC%segment(n)%HI%Isd:OBC%segment(n)%HI%Ied)
            FA_v = 0.0
            do k=1,nz ; FA_v = FA_v + h_in(i,j,k)*(dx_Cv(i,J)*por_face_areaV(i,J,k)) ; enddo
            BT_cont%FA_v_S0(i,J) = FA_v ; BT_cont%FA_v_N0(i,J) = FA_v
            BT_cont%FA_v_SS(i,J) = FA_v ; BT_cont%FA_v_NN(i,J) = FA_v
            BT_cont%vBT_SS(i,J) = 0.0 ; BT_cont%vBT_NN(i,J) = 0.0
          enddo
        else
          do concurrent (i = OBC%segment(n)%HI%Isd:OBC%segment(n)%HI%Ied)
            FA_v = 0.0
            do k=1,nz ; FA_v = FA_v + h_in(i,j+1,k)*(dx_Cv(i,J)*por_face_areaV(i,J,k)) ; enddo
            BT_cont%FA_v_S0(i,J) = FA_v ; BT_cont%FA_v_N0(i,J) = FA_v
            BT_cont%FA_v_SS(i,J) = FA_v ; BT_cont%FA_v_NN(i,J) = FA_v
            BT_cont%vBT_SS(i,J) = 0.0 ; BT_cont%vBT_NN(i,J) = 0.0
          enddo
        endif
      endif
    enddo
  endif

  if (set_BT_cont) then ; if (allocated(BT_cont%h_v)) then
    if (present(v_cor_a)) then
      call h_v_a%alloc(lb=LBOUND(BT_cont%h_v), ub=UBOUND(BT_cont%h_v), source=BT_cont%h_v)
      call meridional_flux_thickness(bxC, v_cor_a, h_in_a, h_S_a, h_N_a, h_v_a, dt, dx_Cv_a, &
                                     IareaT_a, IdyT_a, vol_CFL, marginal_faces, &
                                     por_face_areaV_a, OBC, visc_rem_v_a)
      call h_v_a%copy2F(BT_cont%h_v)
      call h_v_a%free()
    else
      call h_v_a%alloc(lb=LBOUND(BT_cont%h_v), ub=UBOUND(BT_cont%h_v), source=BT_cont%h_v)
      call meridional_flux_thickness(bxC, v_a, h_in_a, h_S_a, h_N_a, h_v_a, dt, dx_Cv_a, &
                                     IareaT_a, IdyT_a, vol_CFL, marginal_faces, &
                                     por_face_areaV_a, OBC, visc_rem_v_a)
      call h_v_a%copy2F(BT_cont%h_v)
      call h_v_a%free()
    endif
  endif ; endif

end subroutine present_vhbt_or_set_BT_cont


!> Calculates the vertically integrated mass or volume fluxes through the meridional faces.
subroutine meridional_BT_mass_flux(bxC, v_a, h_in_a, h_S_a, h_N_a, vhbt_a, dt, &
                                   vol_CFL, OBC, por_face_areaV_a, dx_Cv_a, IareaT_a, IdyT_a)

  type(box_t),                                intent(in)  :: bxC  !< Iteration box for continuity solver
  type(RealArray_t),       intent(in)  :: v_a  !< Meridional velocity [L T-1 ~> m s-1]
  type(RealArray_t),       intent(in)  :: h_in_a !< Layer thickness used to
                                                 !! calculate fluxes [H ~> m or kg m-2]
  type(RealArray_t),       intent(in)  :: h_S_a !< Southern edge thickness in the PPM
                                                !! reconstruction [H ~> m or kg m-2].
  type(RealArray_t),       intent(in)  :: h_N_a !< Northern edge thickness in the PPM
                                                !! reconstruction [H ~> m or kg m-2].
  type(RealArray_t),       intent(inout) :: vhbt_a !< The summed volume flux through meridional
                                                   !! faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  real,                                       intent(in)  :: dt   !< Time increment [T ~> s].
  logical,                 intent(in)  :: vol_CFL !< If true, use the ratio of the open face
                                                   !! lengths to the tracer cell areas when
                                                   !! estimating CFL numbers. Without
                                                   !! aggress_adjust, the default is false; it is
                                                   !! always true with.
  type(ocean_OBC_type),                       pointer     :: OBC  !< Open boundary condition type
                                                                  !! specifies whether, where, and what
                                                                  !! open boundary conditions are used.
  type(RealArray_t),       intent(in)  :: por_face_areaV_a !< fractional open area of V-faces
                                                           !! [nondim]
  type(RealArray_t),       intent(in)  :: dx_Cv_a !< The grid cell's unblocked lengths of the
                                                   !! u/v-faces of the h-cell [L ~> m].
  type(RealArray_t),       intent(in)  :: IareaT_a !< The grid cell's 1/areaT [L-2 ~> m-2].
  type(RealArray_t),       intent(in)  :: IdyT_a !< The grid cell's 1/dyT [L-1 ~> m-1].

  ! Local variables
  real, dimension(v_a%lb(1):v_a%ub(1), v_a%lb(2):v_a%ub(2), v_a%lb(3):v_a%ub(3)) :: &
    vh      ! Volume flux through meridional faces = v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(v_a%lb(1):v_a%ub(1), v_a%lb(2):v_a%ub(2), v_a%lb(3):v_a%ub(3)) :: &
    dvhdv   ! Partial derivative of vh with v [H L ~> m2 or kg m-1].
  integer :: i, j, k, ish, ieh, jsh, jeh, nz, l_seg
  logical :: local_specified_BC, OBC_in_row(v_a%lb(2):v_a%ub(2))
  real, dimension(:,:,:), contiguous, pointer :: v, h_in, h_S, h_N, por_face_areaV
  real, dimension(:,:),   contiguous, pointer :: vhbt
  real, dimension(:,:),   contiguous, pointer :: dx_Cv, IareaT, IdyT

  call v_a%view(v)
  call h_in_a%view(h_in)
  call h_S_a%view(h_S)
  call h_N_a%view(h_N)
  call vhbt_a%view(vhbt)
  call por_face_areaV_a%view(por_face_areaV)
  call dx_Cv_a%view(dx_Cv)
  call IareaT_a%view(IareaT)
  call IdyT_a%view(IdyT)

  call cpu_clock_begin(id_clock_correct)

  local_specified_BC = .false.
  if (associated(OBC)) then ; if (OBC%OBC_pe) then
    local_specified_BC = OBC%specified_v_BCs_exist_globally
  endif ; endif

  ish = bxC%idxS(1) ; ieh = bxC%idxE(1) ; jsh = bxC%idxS(2) ; jeh = bxC%idxE(2) ; nz  = bxC%idxE(3)

  vhbt(:,:) = 0.0

  ! Determining whether there are any OBC points outside of the k-loop should be more efficient.
  OBC_in_row(:) = .false.
  if (local_specified_BC) then
    do j=jsh-1,jeh ; do i=ish,ieh ; if (OBC%segnum_v(i,J) /= 0) then
      if (OBC%segment(abs(OBC%segnum_v(i,J)))%specified) OBC_in_row(j) = .true.
    endif ; enddo ; enddo
  endif

  ! This sets vh and dvhdv.
  do concurrent (k=1:nz, J=jsh-1:jeh, i=ish:ieh)
    call flux_elem(v(i,J,k), h_in(i,J,k), h_in(i,J+1,k), h_S(i,J,k), h_S(i,J+1,k), &
                   h_N(i,J,k), h_N(i,J+1,k), vh(i,J,k), dvhdv(i,J,k), 1.0, dx_Cv(I,j), &
                   IareaT(i,J), IareaT(i,J+1), IdyT(i,J), IdyT(i,J+1), dt, &
                   vol_CFL, por_face_areaV(i,J,k))
    if (local_specified_BC) &
      call flux_elem_OBC(v(i,J,k), h_in(i,J,k), h_in(i,J+1,k), vh(i,J,k), dvhdv(i,J,k), 1.0, &
                         por_face_areaV(i,J,k), dx_Cv(i,J), OBC, OBC%segnum_v(i,J))
  enddo

  do k=1,nz ; do j=jsh-1,jeh ; do i=ish,ieh
    if (OBC_in_row(j) .and. OBC%segnum_v(i,J) /= 0) then
      l_seg = abs(OBC%segnum_v(i,J))
      if (OBC%segment(l_seg)%specified) vh(i,j,k) = OBC%segment(l_seg)%normal_trans(i,J,k)
    endif
  enddo ; enddo ; enddo


  ! Accumulate the barotropic transport.
  do k=1,nz ; do J=jsh-1,jeh ; do i=ish,ieh
    vhbt(i,J) = vhbt(i,J) + vh(i,J,k)
  enddo ; enddo ; enddo

  call cpu_clock_end(id_clock_correct)

end subroutine meridional_BT_mass_flux


!> Sets the effective interface thickness associated with the fluxes at each meridional velocity point,
!! optionally scaling back these thicknesses to account for viscosity and fractional open areas.
subroutine meridional_flux_thickness(bxC, v_a, h_a, h_S_a, h_N_a, h_v_a, dt, dx_Cv_a, IareaT_a, &
                                     IdyT_a, vol_CFL, marginal, por_face_areaV_a, OBC, visc_rem_v_a)
  type(box_t),                               intent(in)    :: bxC  !< Iteration box for continuity solver
  type(RealArray_t),       intent(in)    :: v_a  !< Meridional velocity [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in)    :: h_a  !< Layer thickness used to calculate fluxes,
                                                 !! [H ~> m or kg m-2].
  type(RealArray_t),       intent(in)    :: h_S_a !< South edge thickness in the reconstruction,
                                                  !! [H ~> m or kg m-2].
  type(RealArray_t),       intent(in)    :: h_N_a !< North edge thickness in the reconstruction,
                                                  !! [H ~> m or kg m-2].
  type(RealArray_t),       intent(inout) :: h_v_a !< Effective thickness at meridional faces,
                                                  !! scaled down to account for the effects of
                                                  !! viscosity and the fractional open area
                                                  !! [H ~> m or kg m-2].
  real,                                      intent(in)    :: dt   !< Time increment [T ~> s].
  type(RealArray_t),       intent(in)    :: dx_Cv_a !< The grid cell's unblocked lengths of the
                                                    !! u/v-faces of the h-cell [L ~> m].
  type(RealArray_t),       intent(in)    :: IareaT_a !< The grid cell's 1/areaT [L-2 ~> m-2].
  type(RealArray_t),       intent(in)    :: IdyT_a !< The grid cell's 1/dyT [L-1 ~> m-1].
  logical,                                   intent(in)    :: vol_CFL !< If true, rescale the ratio
                          !! of face areas to the cell areas when estimating the CFL number.
  logical,                                   intent(in)    :: marginal !< If true, report the marginal
                          !! face thicknesses; otherwise report transport-averaged thicknesses.
  type(RealArray_t),       intent(in)    :: por_face_areaV_a  !< fractional open area of V-faces
                                                              !! [nondim]
  type(ocean_OBC_type),                      pointer       :: OBC !< Open boundaries control structure.
  type(RealArray_t), intent(in) :: visc_rem_v_a !< Both the fraction
                          !! of the momentum originally in a layer that remains after a time-step of
                          !! viscosity, and the fraction of a time-step's worth of a barotropic
                          !! acceleration that a layer experiences after viscosity is applied [nondim].
                          !! Visc_rem_v is between 0 (at the bottom) and 1 (far above the bottom).

  ! Local variables
  real :: CFL ! The CFL number based on the local velocity and grid spacing [nondim]
  real :: curv_3 ! A measure of the thickness curvature over a grid length,
                 ! with the same units as h [H ~> m or kg m-2] .
  real :: h_avg  ! The average thickness of a flux [H ~> m or kg m-2].
  real :: h_marg ! The marginal thickness of a flux [H ~> m or kg m-2].
  logical :: local_open_BC
  integer :: i, j, k, ish, ieh, jsh, jeh, n, nz
  real :: dh
  type(box_t) :: bxV
  real, dimension(:,:,:), contiguous, pointer :: v, h, h_S, h_N, h_v, por_face_areaV, visc_rem_v
  real, dimension(:,:),   contiguous, pointer :: dx_Cv, IareaT, IdyT

  call v_a%view(v)
  call h_a%view(h)
  call h_S_a%view(h_S)
  call h_N_a%view(h_N)
  call h_v_a%view(h_v)
  call dx_Cv_a%view(dx_Cv)
  call IareaT_a%view(IareaT)
  call IdyT_a%view(IdyT)
  call por_face_areaV_a%view(por_face_areaV)
  if (visc_rem_v_a%associated()) call visc_rem_v_a%view(visc_rem_v)

  ish = bxC%idxS(1) ; ieh = bxC%idxE(1) ; jsh = bxC%idxS(2) ; jeh = bxC%idxE(2) ; nz  = bxC%idxE(3)

  bxV = bxC%growLo(dim=2, n=1) !< Increase the lower extent of the y-dimension (V-grid)

  ! do concurrent (k=1:nz, J=jsh-1:jeh, i=ish:ieh)
  do concurrent(k=bxV%idxS(3):bxV%idxE(3), &
                j=bxV%idxS(2):bxV%idxE(2), &
                i=bxV%idxS(1):bxV%idxE(1)) ! V-grid
    if (v(i,J,k) > 0.0) then
      if (vol_CFL) then ; CFL = (v(i,J,k) * dt) * (dx_Cv(i,J) * IareaT(i,j))
      else ; CFL = v(i,J,k) * dt * IdyT(i,j) ; endif
      curv_3 = (h_S(i,j,k) + h_N(i,j,k)) - 2.0*h(i,j,k)
      dh = h_S(i,J,k) - h_N(i,J,k)
      if (marginal) then
        h_v(i,J,k) = h_N(i,j,k) + CFL * (dh + 3.0*curv_3*(CFL - 1.0))
      else
        h_v(i,J,k) = h_N(i,j,k) + CFL * (0.5*dh + curv_3*(CFL - 1.5))
      endif
    elseif (v(i,J,k) < 0.0) then
      if (vol_CFL) then ; CFL = (-v(i,J,k)*dt) * (dx_Cv(i,J) * IareaT(i,j+1))
      else ; CFL = -v(i,J,k) * dt * IdyT(i,j+1) ; endif
      curv_3 = (h_S(i,j+1,k) + h_N(i,j+1,k)) - 2.0*h(i,j+1,k)
      dh = h_N(i,j+1,k)-h_S(i,j+1,k)
      if (marginal) then
        h_v(i,J,k) = h_S(i,j+1,k) + CFL * (dh + 3.0*curv_3*(CFL - 1.0))
      else
        h_v(i,J,k) = h_S(i,j+1,k) + CFL * (0.5*dh + curv_3*(CFL - 1.5))
      endif
    else
      !   The choice to use the arithmetic mean here is somewhat arbitrarily, but
      ! it should be noted that h_S(i+1,j,k) and h_N(i,j,k) are usually the same.
      h_v(i,J,k) = 0.5 * (h_S(i,j+1,k) + h_N(i,j,k))
 !    h_marg = (2.0 * h_S(i,j+1,k) * h_N(i,j,k)) / &
 !             (h_S(i,j+1,k) + h_N(i,j,k) + GV%H_subroundoff)
    endif

    if (visc_rem_v_a%associated()) then
      ! Scale back the thickness to account for the effects of viscosity and the fractional open
      ! thickness to give an appropriate non-normalized weight for each layer in determining the
      ! barotropic acceleration.
      h_v(i,J,k) = h_v(i,J,k) * (visc_rem_v(i,J,k) * por_face_areaV(i,J,k))
    else
      h_v(i,J,k) = h_v(i,J,k) * por_face_areaV(i,J,k)
    endif
  enddo

  local_open_BC = .false.
  if (associated(OBC)) local_open_BC = OBC%open_v_BCs_exist_globally
  ! untested - will need to be refactored to be performant on GPUs
  if (local_open_BC) then
    do n = 1, OBC%number_of_segments
      if (OBC%segment(n)%open .and. OBC%segment(n)%is_N_or_S) then
        J = OBC%segment(n)%HI%JsdB
        if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
          if (visc_rem_v_a%associated()) then
            do concurrent (k=1:nz, i = OBC%segment(n)%HI%isd:OBC%segment(n)%HI%ied)
              h_v(i,J,k) = h(i,J,k) * (visc_rem_v(i,J,k) * por_face_areaV(i,J,k))
            enddo
          else
            do concurrent (k=1:nz, i = OBC%segment(n)%HI%isd:OBC%segment(n)%HI%ied)
              h_v(i,J,k) = h(i,J,k) * por_face_areaV(i,J,k)
            enddo
          endif
        else
          if (visc_rem_v_a%associated()) then
            do concurrent (k=1:nz, i = OBC%segment(n)%HI%isd:OBC%segment(n)%HI%ied)
              h_v(i,J,k) = h(i,J+1,k) * (visc_rem_v(i,J,k) * por_face_areaV(i,J,k))
            enddo
          else
            do concurrent (k=1:nz, i = OBC%segment(n)%HI%isd:OBC%segment(n)%HI%ied)
              h_v(i,J,k) = h(i,J+1,k) * por_face_areaV(i,J,k)
            enddo
          endif
        endif
      endif
    enddo
  endif

  call bxV%free()

end subroutine meridional_flux_thickness


!> Returns the barotropic velocity adjustment that gives the desired barotropic (layer-summed) transport.
subroutine meridional_flux_adjust(bxC, v_a, h_in_a, h_S_a, h_N_a, vh_tot_0_a, dvhdv_tot_0_a, &
                                  dv_a, dv_max_CFL_a, dv_min_CFL_a, dt, tol_vel_in, &
                                  tol_eta_in, better_iter, vol_CFL, visc_rem_a, do_I_in_a, &
                                  por_face_areaV_a, dx_Cv_a, IareaT_a, IdyT_a, vhbt_a, vh_3d_a, OBC)
  type(box_t),             intent(in)    :: bxC  !< Iteration box for continuity solver
  type(RealArray_t),       intent(in)    :: v_a  !< Meridional velocity [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in)    :: h_in_a !< Layer thickness used to calculate fluxes
                                                   !! [H ~> m or kg m-2].
  type(RealArray_t),       intent(in)    :: h_S_a !< South edge thickness in the reconstruction
                                                  !! [H ~> m or kg m-2].
  type(RealArray_t),       intent(in)    :: h_N_a !< North edge thickness in the reconstruction
                                                  !! [H ~> m or kg m-2].
  type(RealArray_t),       intent(in)    :: visc_rem_a !< Both the fraction of the momentum
                                                       !! originally in a layer that remains after a
                                                       !! time-step of viscosity, and the fraction
                                                       !! of a time-step's worth of a barotropic
                                                       !! acceleration that a layer experiences
                                                       !! after viscosity is applied [nondim].
                                                       !! Visc_rem is between 0 (at the bottom) and
                                                       !! 1 (far above the bottom).
  type(RealArray_t), optional, intent(in)    :: vhbt_a !< The summed volume flux through meridional
                                                       !! faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(RealArray_t),       intent(in)    :: dv_max_CFL_a !< Maximum acceptable value of dv
                                                         !! [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in)    :: dv_min_CFL_a !< Minimum acceptable value of dv
                                                         !! [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in)    :: vh_tot_0_a !< The summed transport with 0 adjustment
                                                       !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(RealArray_t),       intent(in)    :: dvhdv_tot_0_a !< The partial derivative of dv_err with
                                                          !! dv at 0 adjustment
                                                          !! [H L ~> m2 or kg m-1].
  type(RealArray_t),       intent(inout) :: dv_a !< The barotropic velocity adjustment
                                                 !! [L T-1 ~> m s-1].
  real,                    intent(in)  :: dt      !< Time increment [T ~> s].
  real,                    intent(in)    :: tol_vel_in !< The tolerance for barotropic velocity
                                                       !! discrepancies between the barotropic
                                                       !! solution and the sum of the layer
                                                       !! thicknesses [L T-1 ~> m s-1].
  real,                    intent(in)    :: tol_eta_in !< The tolerance for free-surface height
                                                       !! discrepancies between the barotropic
                                                       !! solution and the sum of the layer
                                                       !! thicknesses [H ~> m or kg m-2].
  logical,                 intent(in)    :: better_iter !< If true, stop corrective iterations
                                                        !! using a velocity-based criterion and
                                                        !! only stop if the iteration is better
                                                        !! than all predecessors.
  logical,                 intent(in)    :: vol_CFL !< If true, use the ratio of the open face
                                                     !! lengths to the tracer cell areas when
                                                     !! estimating CFL numbers. Without
                                                     !! aggress_adjust, the default is false; it is
                                                     !! always true with.
  type(LogicalArray_t),    intent(in)  :: do_I_in_a !< A flag indicating which I values to work on.
  type(RealArray_t),       intent(in)  :: por_face_areaV_a !< fractional open area of V-faces
                                                           !! [nondim]
  type(RealArray_t),       intent(in)    :: dx_Cv_a !< The grid cell's unblocked lengths of the
                                                    !! u/v-faces of the h-cell [L ~> m].
  type(RealArray_t),       intent(in)    :: IareaT_a !< The grid cell's 1/areaT [L-2 ~> m-2].
  type(RealArray_t),       intent(in)    :: IdyT_a !< The grid cell's 1/dyT [L-1 ~> m-1].
  type(RealArray_t), intent(inout) :: vh_3d_a !< Volume flux through meridional
                             !! faces = v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(ocean_OBC_type), optional, pointer :: OBC !< Open boundaries control structure.
  ! Local variables
  real, dimension(v_a%lb(1):v_a%ub(1), v_a%lb(3):v_a%ub(3)) :: &
    vh_aux     ! An auxiliary meridional volume flux [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: &
    dvhdv, &   ! Partial derivative of vh with v [H L ~> m2 or kg m-1].
    v_new      ! The velocity with the correction added [L T-1 ~> m s-1].
  real, dimension(v_a%lb(1):v_a%ub(1)) :: &
    vh_err, &  ! Difference between vhbt and the summed vh [H L2 T-1 ~> m3 s-1 or kg s-1].
    vh_err_best, & ! The smallest value of vh_err found so far [H L2 T-1 ~> m3 s-1 or kg s-1].
    dvhdv_tot,&! Summed partial derivative of vh with u [H L ~> m2 or kg m-1].
    dv_min, &  ! Lower limit on dv correction based on CFL limits and previous iterations [L T-1 ~> m s-1]
    dv_max     ! Upper limit on dv correction based on CFL limits and previous iterations [L T-1 ~> m s-1]
  real :: dv_prev ! The previous value of dv [L T-1 ~> m s-1].
  real :: ddv     ! The change in dv from the previous iteration [L T-1 ~> m s-1].
  real :: tol_eta ! The tolerance for the current iteration [H ~> m or kg m-2].
  real :: tol_vel ! The tolerance for velocity in the current iteration [L T-1 ~> m s-1].
  integer :: i, j, k, nz, itt
  logical :: do_I(v_a%lb(1):v_a%ub(1)), local_OBC, use_vhbt
  integer, parameter :: max_itts = 20
  integer :: ish     !< Start of i index range.
  integer :: ieh     !< End of i index range.
  integer :: jsh     !< Start of j index range.
  integer :: jeh     !< End of j index range.
  real, dimension(:,:,:), contiguous, pointer :: v, h_in, h_S, h_N, visc_rem, por_face_areaV, vh_3d
  real, dimension(:,:),   contiguous, pointer :: vh_tot_0, dvhdv_tot_0, dv, dv_max_CFL, dv_min_CFL
  real, dimension(:,:),   contiguous, pointer :: vhbt, dx_Cv, IareaT, IdyT
  logical, dimension(:,:), contiguous, pointer :: do_I_in

  call v_a%view(v)
  call h_in_a%view(h_in)
  call h_S_a%view(h_S)
  call h_N_a%view(h_N)
  call visc_rem_a%view(visc_rem)
  call vh_tot_0_a%view(vh_tot_0)
  call dvhdv_tot_0_a%view(dvhdv_tot_0)
  call dv_a%view(dv)
  call dv_max_CFL_a%view(dv_max_CFL)
  call dv_min_CFL_a%view(dv_min_CFL)
  call por_face_areaV_a%view(por_face_areaV)
  call dx_Cv_a%view(dx_Cv)
  call IareaT_a%view(IareaT)
  call IdyT_a%view(IdyT)
  if (present(vhbt_a)) call vhbt_a%view(vhbt)
  if (vh_3d_a%associated()) call vh_3d_a%view(vh_3d)
  call do_I_in_a%view(do_I_in)

  local_OBC = .false.
  if (present(OBC)) then
    if (associated(OBC)) then
      local_OBC = OBC%open_u_BCs_exist_globally
    endif
  endif

  use_vhbt = present(vhbt_a)

  ish = bxC%idxS(1) ; ieh = bxC%idxE(1) ; jsh = bxC%idxS(2) ; jeh = bxC%idxE(2) ; nz  = bxC%idxE(3)

  tol_vel = tol_vel_in

  ! NVIDIA needs private arrays to be alloc'ed to prevent data transfers.
  ! GCC doesn't understand map(alloc: ...) for variables also marked private
  !$omp target enter data map(alloc: do_I, dv_max, dv_min, dvhdv_tot, vh_err, vh_err_best, vh_aux)

  !$omp target teams loop &
  !$omp   private(vh_err, vh_err_best, dvhdv_tot, dv_min, dv_max, do_I, vh_aux, itt, tol_eta)
  do J=jsh-1,jeh

    if (vh_3d_a%associated()) then
      do concurrent (k=1:nz, i=ish:ieh)
        vh_aux(i,k) = vh_3d(i,J,k)
      enddo
    endif

    do concurrent (i=ish:ieh)
      dv(i,J) = 0.0 ; do_I(i) = do_I_in(i,J)
      dv_max(i) = dv_max_CFL(i,J) ; dv_min(i) = dv_min_CFL(i,J)
      vh_err(i) = vh_tot_0(i,J)
      if (use_vhbt) vh_err(i) = vh_err(i) - vhbt(i,J)
      dvhdv_tot(i) = dvhdv_tot_0(i,J)
      vh_err_best(i) = abs(vh_err(i))
    enddo

    do itt=1,max_itts
      select case (itt)
        case (:1) ; tol_eta = 1e-6 * tol_eta_in
        case (2)  ; tol_eta = 1e-4 * tol_eta_in
        case (3)  ; tol_eta = 1e-2 * tol_eta_in
        case default ; tol_eta = tol_eta_in
      end select

      do concurrent (i=ish:ieh)
        if (vh_err(i) > 0.0) then ; dv_max(i) = dv(i,j)
        elseif (vh_err(i) < 0.0) then ; dv_min(i) = dv(i,j)
        else ; do_I(i) = .false. ; endif
      enddo

      do concurrent (i=ish:ieh, do_I(i)) &
          & DO_LOCALITY(local(ddv, dv_prev))
        if ((dt * min(IareaT(i,j),IareaT(i,j+1))*abs(vh_err(i)) > tol_eta) .or. &
            (better_iter .and. ((abs(vh_err(i)) > tol_vel * dvhdv_tot(i)) .or. &
                                  (abs(vh_err(i)) > vh_err_best(i))) )) then
          !   Use Newton's method, provided it stays bounded.  Otherwise bisect
          ! the value with the appropriate bound.
          ddv = -vh_err(i) / dvhdv_tot(i)
          dv_prev = dv(i,j)
          dv(i,j) = dv(i,j) + ddv
          if (abs(ddv) < 1.0e-15*abs(dv(i,j))) then
            do_I(i) = .false. ! ddv is small enough to quit.
          elseif (ddv > 0.0) then
            if (dv(i,j) >= dv_max(i)) then
              dv(i,j) = 0.5*(dv_prev + dv_max(i))
              if (dv_max(i) - dv_prev < 1.0e-15*abs(dv(i,j))) do_I(i) = .false.
            endif
          else ! ddv(i) < 0.0
            if (dv(i,j) <= dv_min(i)) then
              dv(i,j) = 0.5*(dv_prev + dv_min(i))
              if (dv_prev - dv_min(i) < 1.0e-15*abs(dv(i,j))) do_I(i) = .false.
            endif
          endif
        else
          do_I(i) = .false.
        endif
      enddo

      ! Below conditional compilation is to control whether early exit happens when compiled with
      ! OpenMP - compiling with OpenMP prevents early exit. Without OpenMP, enables early exit.
      ! Early exit saves time on CPU, but causes other loops to be serialized on GPU.
      !$ if (.false.) then
      if (.not. any(do_I(ish:ieh))) exit
      !$ endif

      if ((itt < max_itts) .or. vh_3d_a%associated()) then
        do concurrent (i=ish:ieh)
          vh_err(i) = 0.0 ; dvhdv_tot(i) = 0.0
          if (use_vhbt) vh_err(i) = -vhbt(i,J)
        enddo
        do k=1,nz ; do concurrent (i=ish:ieh, do_I(i)) DO_LOCALITY(local(v_new, dvhdv))
          v_new = v(i,J,k) + dv(i,j) * visc_rem(i,j,k)
          call flux_elem(v_new, h_in(i,J,k), h_in(i,J+1,k), h_S(i,J,k), h_S(i,J+1,k), &
                         h_N(i,J,k), h_N(i,J+1,k), vh_aux(i,k), dvhdv, visc_rem(i,J,k), &
                         dx_Cv(i,J), IareaT(i,J), IareaT(i,J+1), IdyT(i,J), IdyT(i,J+1), &
                         dt, vol_CFL, por_face_areaV(i,J,k))
          if (local_OBC) &
            call flux_elem_OBC(v_new, h_in(i,J,k), h_in(i,J+1,k), vh_aux(i,k), &
                               dvhdv, visc_rem(i,J,k), por_face_areaV(i,J,k), &
                               dx_Cv(i,J), OBC, OBC%segnum_v(i,J))
          vh_err(i) = vh_err(i) + vh_aux(i,k)
          dvhdv_tot(i) = dvhdv_tot(i) + dvhdv
        enddo ; enddo
        do concurrent (i=ish:ieh, do_I(i))
          vh_err_best(i) = min(vh_err_best(i), abs(vh_err(i)))
        enddo
      endif
    enddo ! itt-loop

    ! If there are any faces which have not converged to within the tolerance,
    ! so-be-it, or else use a final upwind correction?
    ! This never seems to happen with 20 iterations as max_itt.

    if (vh_3d_a%associated()) then
      do concurrent (k=1:nz, i=ish:ieh)
        vh_3d(i,J,k) = vh_aux(i,k)
      enddo
    endif
  enddo ! j-loop

  !$omp target exit data map(release: do_I, dv_max, dv_min, dvhdv_tot, vh_err, vh_err_best, vh_aux)

end subroutine meridional_flux_adjust


!> Sets of a structure that describes the meridional barotropic volume or mass fluxes as a
!! function of barotropic flow to agree closely with the sum of the layer's transports.
subroutine set_merid_BT_cont(bxC, v_a, h_in_a, h_S_a, h_N_a, BT_cont, dv0_a, vh_tot_0_a, &
                             dvhdv_tot_0_a, dv_max_CFL_a, dv_min_CFL_a, dt, vol_CFL, &
                             visc_rem_a, visc_rem_max_a, do_I_a, por_face_areaV_a, &
                             dyCv_a, dx_Cv_a, IareaT_a, IdyT_a)
  type(box_t),                                intent(in)    :: bxC  !< Iteration box for continuity solver
  type(RealArray_t),       intent(in)    :: v_a    !< Meridional velocity [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in)    :: h_in_a !< Layer thickness used to calculate fluxes,
                                                   !! [H ~> m or kg m-2].
  type(RealArray_t),       intent(in)    :: h_S_a !< South edge thickness in the reconstruction,
                                                  !! [H ~> m or kg m-2].
  type(RealArray_t),       intent(in)    :: h_N_a !< North edge thickness in the reconstruction,
                                                  !! [H ~> m or kg m-2].
  type(BT_cont_type),                         intent(inout) :: BT_cont !< A structure with elements
                       !! that describe the effective open face areas as a function of barotropic flow.
  type(RealArray_t),       intent(in)    :: dv0_a !< The barotropic velocity increment that gives 0
                                                  !! transport [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in)    :: vh_tot_0_a !< The summed transport with 0 adjustment
                                                       !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(RealArray_t),       intent(in)    :: dvhdv_tot_0_a !< The partial derivative of du_err with
                                                          !! dv at 0 adjustment
                                                          !! [H L ~> m2 or kg m-1].
  type(RealArray_t),       intent(in)    :: dv_max_CFL_a !< Maximum acceptable value of dv
                                                         !! [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in)    :: dv_min_CFL_a !< Minimum acceptable value of dv
                                                         !! [L T-1 ~> m s-1].
  real,                                       intent(in)    :: dt   !< Time increment [T ~> s].
  logical,                 intent(in)    :: vol_CFL !< If true, use the ratio of the open face
                                                     !! lengths to the tracer cell areas when
                                                     !! estimating CFL numbers. Without
                                                     !! aggress_adjust, the default is false; it is
                                                     !! always true with.
  type(RealArray_t),       intent(in)    :: visc_rem_a !< Both the fraction of the momentum
                                                       !! originally in a layer that remains after a
                                                       !! time-step of viscosity, and the fraction
                                                       !! of a time-step's worth of a barotropic
                                                       !! acceleration that a layer experiences
                                                       !! after viscosity is applied [nondim].
                                                       !! Visc_rem is between 0 (at the bottom) and
                                                       !! 1 (far above the bottom).
  type(RealArray_t),       intent(in)    :: visc_rem_max_a !< Maximum allowable visc_rem [nondim]
  type(LogicalArray_t),    intent(in) :: do_I_a  !< A logical flag indicating which I
                                                 !! values to work on.
  type(RealArray_t),       intent(in)    :: por_face_areaV_a !< fractional open area of V-faces
                                                             !! [nondim]
  type(RealArray_t),       intent(in)    :: dyCv_a !< The dy spacing at v points [L ~> m].
  type(RealArray_t),       intent(in)    :: dx_Cv_a !< The grid cell's unblocked lengths of the
                                                    !! u/v-faces of the h-cell [L ~> m].
  type(RealArray_t),       intent(in)    :: IareaT_a !< The grid cell's 1/areaT [L-2 ~> m-2].
  type(RealArray_t),       intent(in)    :: IdyT_a !< The grid cell's 1/dyT [L-1 ~> m-1].
  ! Local variables
  real, dimension(v_a%lb(1):v_a%ub(1)) :: &
    dvL, dvR, &        ! The barotropic velocity increments that give the southerly
                       ! (dvL) and northerly (dvR) test velocities [L T-1 ~> m s-1].
    dv_CFL, &          ! The velocity increment that corresponds to CFL_min [L T-1 ~> m s-1].
    FAmt_L, FAmt_R, &  ! The summed effective marginal face areas for the 3
    FAmt_0, &          ! test velocities [H L ~> m2 or kg m-1].
    vhtot_L, &         ! The summed transport with the southerly (vhtot_L) and
    vhtot_R            ! and northerly (vhtot_R) test velocities [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: &
    v_L, v_R, &        ! The southerly (v_L), northerly (v_R), and zero-barotropic
    v_0, &             ! transport (v_0) layer test velocities [L T-1 ~> m s-1].
    dvhdv_L, &         ! The effective layer marginal face areas with the southerly
    dvhdv_R, &         ! (_L), northerly (_R), and zero-barotropic (_0) test
    dvhdv_0, &         ! velocities [H L ~> m2 or kg m-1].
    vh_L, vh_R, &      ! The layer transports with the southerly (_L), northerly (_R)
    vh_0               ! and zero-barotropic (_0) test velocities [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: FA_0         ! The effective face area with 0 barotropic transport [H L ~> m2 or kg m-1].
  real :: FA_avg       ! The average effective face area [H L ~> m2 or kg m-1], nominally given by
                       ! the realized transport divided by the barotropic velocity.
  real :: visc_rem_lim ! The larger of visc_rem and min_visc_rem [nondim]  This
                       ! limiting is necessary to keep the inverse of visc_rem
                       ! from leading to large CFL numbers.
  real :: min_visc_rem ! The smallest permitted value for visc_rem that is used
                       ! in finding the barotropic velocity that changes the
                       ! flow direction [nondim].  This is necessary to keep the inverse
                       ! of visc_rem from leading to large CFL numbers.
  real :: CFL_min      ! A minimal increment in the CFL to try to ensure that the
                       ! flow is truly upwind [nondim]
  real :: Idt          ! The inverse of the time step [T-1 ~> s-1].
  integer :: i, j, k, nz
  integer :: ish  !< Start of i index range.
  integer :: ieh  !< End of i index range.
  integer :: jsh  !< Start of j index range.
  integer :: jeh  !< End of j index range.
  real, dimension(:,:,:), contiguous, pointer :: v, h_in, h_S, h_N, visc_rem, por_face_areaV
  real, dimension(:,:),   contiguous, pointer :: dv0, vh_tot_0, dvhdv_tot_0, dv_max_CFL, &
                                                 dv_min_CFL, visc_rem_max
  real, dimension(:,:),   contiguous, pointer :: dyCv, dx_Cv, IareaT, IdyT
  logical, dimension(:,:), contiguous, pointer :: do_I

  call v_a%view(v)
  call h_in_a%view(h_in)
  call h_S_a%view(h_S)
  call h_N_a%view(h_N)
  call dv0_a%view(dv0)
  call vh_tot_0_a%view(vh_tot_0)
  call dvhdv_tot_0_a%view(dvhdv_tot_0)
  call dv_max_CFL_a%view(dv_max_CFL)
  call dv_min_CFL_a%view(dv_min_CFL)
  call visc_rem_a%view(visc_rem)
  call visc_rem_max_a%view(visc_rem_max)
  call por_face_areaV_a%view(por_face_areaV)
  call dyCv_a%view(dyCv)
  call dx_Cv_a%view(dx_Cv)
  call IareaT_a%view(IareaT)
  call IdyT_a%view(IdyT)
  call do_I_a%view(do_I)

  ish = bxC%idxS(1) ; ieh = bxC%idxE(1) ; jsh = bxC%idxS(2) ; jeh = bxC%idxE(2) ; nz  = bxC%idxE(3)
  Idt = 1.0 / dt
  min_visc_rem = 0.1 ; CFL_min = 1e-6

  !$omp target enter data map(alloc: dvL, dvR, dv_CFL, FAmt_L, FAmt_R, FAmt_0, vhtot_L, vhtot_R)

  !$omp target teams loop private(dvL, dvR, dv_CFL, FAmt_L, FAmt_R, FAmt_0, vhtot_L, vhtot_R)
  do J=jsh-1,jeh
    ! Determine the southerly- and northerly- fluxes. Choose a sufficiently
    ! negative velocity correction for the northerly-flux, and a sufficiently
    ! positive correction for the southerly-flux.
    do concurrent (i=ish:ieh)
      dv_CFL(i) = (CFL_min * Idt) * dyCv(i,J)
      dvR(i) = min(0.0,dv0(i,J) - dv_CFL(i))
      dvL(i) = max(0.0,dv0(i,J) + dv_CFL(i))
      FAmt_L(i) = 0.0 ; FAmt_R(i) = 0.0 ; FAmt_0(i) = 0.0
      vhtot_L(i) = 0.0 ; vhtot_R(i) = 0.0
    enddo

    ! not parallelized on k because of dvR/L are calculated per column
    ! nvfortran do concurrent poor performance when k is inside
    do k=1,nz ; do concurrent (i=ish:ieh, do_I(i,J)) DO_LOCALITY(local(visc_rem_lim))
      visc_rem_lim = max(visc_rem(i,J,k), min_visc_rem*visc_rem_max(i,J))
      if (visc_rem_lim > 0.0) then ! This is almost always true for ocean points.
        if (v(i,J,k) + dvR(i)*visc_rem_lim > -dv_CFL(i)*visc_rem(i,J,k)) &
          dvR(i) = -(v(i,J,k) + dv_CFL(i)*visc_rem(i,J,k)) / visc_rem_lim
        if (v(i,J,k) + dvL(i)*visc_rem_lim < dv_CFL(i)*visc_rem(i,J,k)) &
          dvL(i) = -(v(i,J,k) - dv_CFL(i)*visc_rem(i,J,k)) / visc_rem_lim
      endif
    enddo ; enddo

    do k=1,nz ; do concurrent (i=ish:ieh, do_I(i,j)) &
        & DO_LOCALITY(local(v_0, v_L, v_R, dvhdv_0, dvhdv_L, dvhdv_R, vh_0, vh_L, vh_R))
      v_L = v(I,J,k) + dvL(i) * visc_rem(i,J,k)
      v_R = v(I,J,k) + dvR(i) * visc_rem(i,J,k)
      v_0 = v(I,J,k) + dv0(i,J) * visc_rem(i,J,k)
      call flux_elem(v_0, h_in(i,J,k), h_in(i,J+1,k), h_S(i,J,k), h_S(i,J+1,k), &
                     h_N(i,J,k), h_N(i,J+1,k), vh_0, dvhdv_0, visc_rem(i,J,k), &
                     dx_Cv(i,J), IareaT(i,J), IareaT(i,J+1), IdyT(i,J), &
                     IdyT(i,J+1), dt, vol_CFL, por_face_areaV(i,J,k))
      call flux_elem(v_L, h_in(i,J,k), h_in(i,J+1,k), h_S(i,J,k), h_S(i,J+1,k), &
                     h_N(i,J,k), h_N(i,J+1,k), vh_L, dvhdv_L, visc_rem(i,J,k), &
                     dx_Cv(i,J), IareaT(i,J), IareaT(i,J+1), IdyT(i,J), &
                     IdyT(i,J+1), dt, vol_CFL, por_face_areaV(i,J,k))
      call flux_elem(v_R, h_in(i,J,k), h_in(i,J+1,k), h_S(i,J,k), h_S(i,J+1,k), &
                     h_N(i,J,k), h_N(i,J+1,k), vh_R, dvhdv_R, visc_rem(i,J,k), &
                     dx_Cv(i,J), IareaT(i,J), IareaT(i,J+1), IdyT(i,J), &
                     IdyT(i,J+1), dt, vol_CFL, por_face_areaV(i,J,k))
      FAmt_0(i) = FAmt_0(i) + dvhdv_0
      FAmt_L(i) = FAmt_L(i) + dvhdv_L
      FAmt_R(i) = FAmt_R(i) + dvhdv_R
      vhtot_L(i) = vhtot_L(i) + vh_L
      vhtot_R(i) = vhtot_R(i) + vh_R
    enddo ; enddo

    do concurrent (i=ish:ieh) DO_LOCALITY(local(FA_0, FA_Avg))
      if (do_I(i,j)) then
        FA_0 = FAmt_0(i) ; FA_avg = FAmt_0(i)
        if ((dvL(i) - dv0(i,J)) /= 0.0) &
          FA_avg = vhtot_L(i) / (dvL(i) - dv0(i,J))
        if (FA_avg > max(FA_0, FAmt_L(i))) then ; FA_avg = max(FA_0, FAmt_L(i))
        elseif (FA_avg < min(FA_0, FAmt_L(i))) then ; FA_0 = FA_avg ; endif
        BT_cont%FA_v_S0(i,J) = FA_0 ; BT_cont%FA_v_SS(i,J) = FAmt_L(i)
        if (abs(FA_0-FAmt_L(i)) <= 1e-12*FA_0) then ; BT_cont%vBT_SS(i,J) = 0.0 ; else
          BT_cont%vBT_SS(i,J) = (1.5 * (dvL(i) - dv0(i,J))) * &
                      ((FAmt_L(i) - FA_avg) / (FAmt_L(i) - FA_0))
        endif

        FA_0 = FAmt_0(i) ; FA_avg = FAmt_0(i)
        if ((dvR(i) - dv0(i,j)) /= 0.0) &
          FA_avg = vhtot_R(i) / (dvR(i) - dv0(i,j))
        if (FA_avg > max(FA_0, FAmt_R(i))) then ; FA_avg = max(FA_0, FAmt_R(i))
        elseif (FA_avg < min(FA_0, FAmt_R(i))) then ; FA_0 = FA_avg ; endif
        BT_cont%FA_v_N0(i,J) = FA_0 ; BT_cont%FA_v_NN(i,J) = FAmt_R(i)
        if (abs(FAmt_R(i) - FA_0) <= 1e-12*FA_0) then ; BT_cont%vBT_NN(i,J) = 0.0 ; else
          BT_cont%vBT_NN(i,J) = (1.5 * (dvR(i) - dv0(i,j))) * &
                      ((FAmt_R(i) - FA_avg) / (FAmt_R(i) - FA_0))
        endif
      else
        BT_cont%FA_v_S0(i,J) = 0.0 ; BT_cont%FA_v_SS(i,J) = 0.0
        BT_cont%FA_v_N0(i,J) = 0.0 ; BT_cont%FA_v_NN(i,J) = 0.0
        BT_cont%vBT_SS(i,J) = 0.0 ; BT_cont%vBT_NN(i,J) = 0.0
      endif
    enddo
  enddo

  !$omp target exit data map(release: dvL, dvR, dv_CFL, FAmt_L, FAmt_R, FAmt_0, vhtot_L, vhtot_R)

end subroutine set_merid_BT_cont

!> Calculates left/right edge values for PPM reconstruction.
subroutine PPM_reconstruction_x_fortran(bxH, h_in_a, h_W_a, h_E_a, mask2dT_a, h_min, monotonic, simple_2nd, OBC)
  type(Box_t),                       intent(in)  :: bxH  !< H-grid iteration Box
  type(RealArray_t),  intent(in)    :: h_in_a !< Layer thickness [H ~> m or kg m-2].
  type(RealArray_t),  intent(inout) :: h_W_a  !< West edge thickness in the reconstruction
                                              !! [H ~> m or kg m-2].
  type(RealArray_t),  intent(inout) :: h_E_a  !< East edge thickness in the reconstruction
                                              !! [H ~> m or kg m-2].
  type(RealArray_t), intent(in)  :: mask2dT_a !< 0 for land points and 1 for ocean points
                                              !! on the h-grid [nondim]
  real,                              intent(in)  :: h_min !< The minimum thickness
                    !! that can be obtained by a concave parabolic fit [H ~> m or kg m-2]
  logical,                           intent(in)  :: monotonic !< If true, use the
                    !! Colella & Woodward monotonic limiter.
                    !! Otherwise use a simple positive-definite limiter.
  logical,                           intent(in)  :: simple_2nd !< If true, use the
                    !! arithmetic mean thicknesses as the default edge values
                    !! for a simple 2nd order scheme.
  type(ocean_OBC_type),              pointer     :: OBC !< Open boundaries control structure.
  integer :: k      !< vertical grid index

  ! Local variables with useful mnemonic names.
  real, dimension(:,:,:), allocatable  :: slp ! The slopes per grid point [H ~> m or kg m-2]
  real, parameter :: oneSixth = 1./6.  ! [nondim]
  real :: h_ip1, h_im1 ! Neighboring thicknesses or sensibly extrapolated values [H ~> m or kg m-2]
  real :: dMx, dMn     ! The difference between the local thickness and the maximum (dMx) or
                       ! minimum (dMn) of the surrounding values [H ~> m or kg m-2]
  integer :: i, j
  integer :: n
  logical :: local_open_BC
  type(OBC_segment_type), pointer :: segment => NULL()
  type(Box_t) :: bx, bxE

  real, dimension(:,:),   contiguous, pointer :: mask2dT
  real, dimension(:,:,:), contiguous, pointer :: h_in
  real, dimension(:,:,:), contiguous, pointer :: h_W
  real, dimension(:,:,:), contiguous, pointer :: h_E

  ! Get the views for containers (subroutine arguments)
  call h_W_a%view(h_W)
  call h_E_a%view(h_E)
  call h_in_a%view(h_in)
  call mask2dT_a%view(mask2dT)

  ! Allocate local slope array using bounds from h_in_a
  allocate(slp(h_in_a%lb(1):h_in_a%ub(1),h_in_a%lb(2):h_in_a%ub(2),h_in_a%lb(3):h_in_a%ub(3)))

  local_open_BC = .false.
  if (associated(OBC)) then
    local_open_BC = OBC%open_u_BCs_exist_globally
  endif

  ! The local iteration box is expanded by one element in the j-dimension
  bx = bxH%grow(dim=1,n=1)

  ! Create an second box that extent the iteration space by two in the i-dimension
  bxE = bxH%grow(dim=1,n=2)

  !$omp target enter data map(alloc: slp)

  if (simple_2nd) then
    ! untested
    do concurrent(k=bx%idxS(3):bx%idxE(3), &
                  j=bx%idxS(2):bx%idxE(2), &
                  i=bx%idxS(1):bx%idxE(1)) ! Local box (bx)
      h_im1 = mask2dT(i-1,j) * h_in(i-1,j,k) + (1.0-mask2dT(i-1,j)) * h_in(i,j,k)
      h_ip1 = mask2dT(i+1,j) * h_in(i+1,j,k) + (1.0-mask2dT(i+1,j)) * h_in(i,j,k)
      h_W(i,j,k) = 0.5*( h_im1 + h_in(i,j,k) )
      h_E(i,j,k) = 0.5*( h_ip1 + h_in(i,j,k) )
    enddo
  else
    do concurrent(k=bxE%idxS(3):bxE%idxE(3), &
                  j=bxE%idxS(2):bxE%idxE(2), &
                  i=bxE%idxS(1):bxE%idxE(1)) ! Expanded box (bxE)
      if ((mask2dT(i-1,j) * mask2dT(i,j) * mask2dT(i+1,j)) == 0.0) then
        slp(i,j,k) = 0.0
      else
        ! This uses a simple 2nd order slope.
        slp(i,j,k) = 0.5 * (h_in(i+1,j,k) - h_in(i-1,j,k))
        ! Monotonic constraint, see Eq. B2 in Lin 1994, MWR (132)
        dMx = max(h_in(i+1,j,k), h_in(i-1,j,k), h_in(i,j,k)) - h_in(i,j,k)
        dMn = h_in(i,j,k) - min(h_in(i+1,j,k), h_in(i-1,j,k), h_in(i,j,k))
        slp(i,j,k) = sign(1.,slp(i,j,k)) * min(abs(slp(i,j,k)), 2. * min(dMx, dMn))
                ! * (mask2dT(i-1,j) * mask2dT(i,j) * mask2dT(i+1,j))
      endif
    enddo

    if (local_open_BC) then
      ! untested
      do n=1, OBC%number_of_segments
        segment => OBC%segment(n)
        if (.not. segment%on_pe) cycle
        if (segment%is_E_or_W) then
          I=segment%HI%IsdB
          do concurrent(k=bx%idxS(3):bx%idxE(3),j=segment%HI%jsd:segment%HI%jed)
            slp(i+1,j,k) = 0.0
            slp(i,j,k) = 0.0
          enddo
        endif
      enddo
    endif

    do concurrent(k=bx%idxS(3):bx%idxE(3), &
                  j=bx%idxS(2):bx%idxE(2), &
                  i=bx%idxS(1):bx%idxE(1))
      ! Neighboring values should take into account any boundaries.  The 3
      ! following sets of expressions are equivalent.
    ! h_im1 = h_in(i-1,j,k) ; if (mask2dT(i-1,j) < 0.5) h_im1 = h_in(i,j)
    ! h_ip1 = h_in(i+1,j,k) ; if (mask2dT(i+1,j) < 0.5) h_ip1 = h_in(i,j)
      h_im1 = mask2dT(i-1,j) * h_in(i-1,j,k) + (1.0-mask2dT(i-1,j)) * h_in(i,j,k)
      h_ip1 = mask2dT(i+1,j) * h_in(i+1,j,k) + (1.0-mask2dT(i+1,j)) * h_in(i,j,k)
      ! Left/right values following Eq. B2 in Lin 1994, MWR (132)
      h_W(i,j,k) = 0.5*( h_im1 + h_in(i,j,k) ) + oneSixth*( slp(i-1,j,k) - slp(i,j,k) )
      h_E(i,j,k) = 0.5*( h_ip1 + h_in(i,j,k) ) + oneSixth*( slp(i,j,k) - slp(i+1,j,k) )
    enddo
  endif

  if (local_open_BC) then
    ! untested
    do n=1, OBC%number_of_segments
      segment => OBC%segment(n)
      if (.not. segment%on_pe) cycle
      if (segment%direction == OBC_DIRECTION_E) then
        I=segment%HI%IsdB
        do concurrent(k=bx%idxS(3):bx%idxE(3),j=segment%HI%jsd:segment%HI%jed)
          h_W(i+1,j,k) = h_in(i,j,k)
          h_E(i+1,j,k) = h_in(i,j,k)
          h_W(i,j,k) = h_in(i,j,k)
          h_E(i,j,k) = h_in(i,j,k)
        enddo
      elseif (segment%direction == OBC_DIRECTION_W) then
        I=segment%HI%IsdB
        do concurrent(k=bx%idxS(3):bx%idxE(3),j=segment%HI%jsd:segment%HI%jed)
          h_W(i,j,k) = h_in(i+1,j,k)
          h_E(i,j,k) = h_in(i+1,j,k)
          h_W(i+1,j,k) = h_in(i+1,j,k)
          h_E(i+1,j,k) = h_in(i+1,j,k)
        enddo
      endif
    enddo
  endif

  if (monotonic) then
    ! untested
    call PPM_limit_cw84(bx, h_in_a, h_W_a, h_E_a)
  else
    call PPM_limit_pos(bx, h_in_a, h_W_a, h_E_a, h_min)
  endif

  !$omp target exit data map(release: slp)

  ! Deallocate local temporary array
  if(allocated(slp)) deallocate(slp)

  ! Deallocate local iteration boxes
  call bx%free()
  call bxE%free()

  return

end subroutine PPM_reconstruction_x_fortran

!> Calculates left/right edge values for PPM reconstruction.
subroutine PPM_reconstruction_y_fortran(bxH, h_in_a, h_S_a, h_N_a, mask2dT_a, h_min, monotonic, simple_2nd, OBC)
  type(Box_t),                       intent(in)  :: bxH  !< H-grid iteration Box
  type(RealArray_t),  intent(in)    :: h_in_a !< Layer thickness [H ~> m or kg m-2].
  type(RealArray_t),  intent(inout) :: h_S_a  !< South edge thickness in the reconstruction
                                            !! [H ~> m or kg m-2].
  type(RealArray_t),  intent(inout) :: h_N_a  !< North edge thickness in the reconstruction,
                                            !! [H ~> m or kg m-2].
  type(RealArray_t), intent(in)  :: mask2dT_a !< 0 for land points and 1 for ocean points
                                              !! on the h-grid [nondim]
  real,                             intent(in)  :: h_min   !< The minimum thickness
                    !! that can be obtained by a concave parabolic fit [H ~> m or kg m-2]
  logical,                           intent(in)  :: monotonic !< If true, use the
                    !! Colella & Woodward monotonic limiter.
                    !! Otherwise use a simple positive-definite limiter.
  logical,                           intent(in)  :: simple_2nd !< If true, use the
                    !! arithmetic mean thicknesses as the default edge values
                    !! for a simple 2nd order scheme.
  type(ocean_OBC_type),              pointer     :: OBC !< Open boundaries control structure.
  integer :: k      !< vertical grid index

  ! Local variables with useful mnemonic names.
  !real, dimension(SZI_(G),SZJ_(G),SZK_(GV))  :: slp ! The slopes per grid point [H ~> m or kg m-2]
  real, dimension(:,:,:), allocatable  :: slp ! The slopes per grid point [H ~> m or kg m-2]

  real, parameter :: oneSixth = 1./6.      ! [nondim]
  real :: h_jp1, h_jm1 ! Neighboring thicknesses or sensibly extrapolated values [H ~> m or kg m-2]
  real :: dMx, dMn     ! The difference between the local thickness and the maximum (dMx) or
                       ! minimum (dMn) of the surrounding values [H ~> m or kg m-2]
  integer :: i, j
  integer :: n, ndims
  logical :: local_open_BC
  type(OBC_segment_type), pointer :: segment => NULL()

  real, dimension(:,:),   contiguous, pointer :: mask2dT
  real, dimension(:,:,:), contiguous, pointer :: h_in
  real, dimension(:,:,:), contiguous, pointer :: h_S
  real, dimension(:,:,:), contiguous, pointer :: h_N

  type(Box_t) :: bx, bxE

  ! Get the views for containers (subroutine arguments)
  call h_S_a%view(h_S)
  call h_N_a%view(h_N)
  call h_in_a%view(h_in)
  call mask2dT_a%view(mask2dT)

  ! Get the lower and upper bounds of the h_in_a container
  ! and allocate a local array
  !JMD KLUDGE: fix with more elegant solution
  allocate(slp(h_in_a%lb(1):h_in_a%ub(1),h_in_a%lb(2):h_in_a%ub(2),h_in_a%lb(3):h_in_a%ub(3)))

  local_open_BC = .false.
  if (associated(OBC)) then
    local_open_BC = OBC%open_v_BCs_exist_globally
  endif

  ! Local iteration box extends the h-grid by one element in the j-dimension
  bx = bxH%grow(dim=2,n=1)

  ! Extended iteration box extends the h-grid by two elements in the j-dimension
  bxE = bxH%grow(dim=2,n=2)

  !$omp target enter data map(alloc: slp)

  if (simple_2nd) then
    ! untested
    do concurrent(k=bx%idxS(3):bx%idxE(3), &
                  j=bx%idxS(2):bx%idxE(2), &
                  i=bx%idxS(1):bx%idxE(1))
      h_jm1 = mask2dT(i,j-1) * h_in(i,j-1,k) + (1.0-mask2dT(i,j-1)) * h_in(i,j,k)
      h_jp1 = mask2dT(i,j+1) * h_in(i,j+1,k) + (1.0-mask2dT(i,j+1)) * h_in(i,j,k)
      h_S(i,j,k) = 0.5*( h_jm1 + h_in(i,j,k) )
      h_N(i,j,k) = 0.5*( h_jp1 + h_in(i,j,k) )
    enddo
  else
    do concurrent(k=bxE%idxS(3):bxE%idxE(3), &
                  j=bxE%idxS(2):bxE%idxE(2), &
                  i=bxE%idxS(1):bxE%idxE(1)) ! Expanded box (bxE)
      if ((mask2dT(i,j-1) * mask2dT(i,j) * mask2dT(i,j+1)) == 0.0) then
        slp(i,j,k) = 0.0
      else
        ! This uses a simple 2nd order slope.
        slp(i,j,k) = 0.5 * (h_in(i,j+1,k) - h_in(i,j-1,k))
        ! Monotonic constraint, see Eq. B2 in Lin 1994, MWR (132)
        dMx = max(h_in(i,j+1,k), h_in(i,j-1,k), h_in(i,j,k)) - h_in(i,j,k)
        dMn = h_in(i,j,k) - min(h_in(i,j+1,k), h_in(i,j-1,k), h_in(i,j,k))
        slp(i,j,k) = sign(1.,slp(i,j,k)) * min(abs(slp(i,j,k)), 2. * min(dMx, dMn))
                ! * (mask2dT(i,j-1) * mask2dT(i,j) * mask2dT(i,j+1))
      endif
    enddo

    if (local_open_BC) then
      ! untested
      do n=1, OBC%number_of_segments
        segment => OBC%segment(n)
        if (.not. segment%on_pe) cycle
        if (segment%is_N_or_S) then
          J=segment%HI%JsdB
          do concurrent(k=bx%idxS(3):bx%idxE(3),i=segment%HI%isd:segment%HI%ied)
            slp(i,j+1,k) = 0.0
            slp(i,j,k) = 0.0
          enddo
        endif
      enddo
    endif

    do concurrent(k=bx%idxS(3):bx%idxE(3), &
                  j=bx%idxS(2):bx%idxE(2), &
                  i=bx%idxS(1):bx%idxE(1))
      ! Neighboring values should take into account any boundaries.  The 3
      ! following sets of expressions are equivalent.
      h_jm1 = mask2dT(i,j-1) * h_in(i,j-1,k) + (1.0-mask2dT(i,j-1)) * h_in(i,j,k)
      h_jp1 = mask2dT(i,j+1) * h_in(i,j+1,k) + (1.0-mask2dT(i,j+1)) * h_in(i,j,k)
      ! Left/right values following Eq. B2 in Lin 1994, MWR (132)
      h_S(i,j,k) = 0.5*( h_jm1 + h_in(i,j,k) ) + oneSixth*( slp(i,j-1,k) - slp(i,j,k) )
      h_N(i,j,k) = 0.5*( h_jp1 + h_in(i,j,k) ) + oneSixth*( slp(i,j,k) - slp(i,j+1,k) )
    enddo
  endif

  if (local_open_BC) then
    ! untested
    do n=1, OBC%number_of_segments
      segment => OBC%segment(n)
      if (.not. segment%on_pe) cycle
      if (segment%direction == OBC_DIRECTION_N) then
        J=segment%HI%JsdB
        do concurrent(k=bx%idxS(3):bx%idxE(3),i=segment%HI%isd:segment%HI%ied)
          h_S(i,j+1,k) = h_in(i,j,k)
          h_N(i,j+1,k) = h_in(i,j,k)
          h_S(i,j,k) = h_in(i,j,k)
          h_N(i,j,k) = h_in(i,j,k)
        enddo
      elseif (segment%direction == OBC_DIRECTION_S) then
        J=segment%HI%JsdB
        do concurrent(k=bx%idxS(3):bx%idxE(3),i=segment%HI%isd:segment%HI%ied)
          h_S(i,j,k) = h_in(i,j+1,k)
          h_N(i,j,k) = h_in(i,j+1,k)
          h_S(i,j+1,k) = h_in(i,j+1,k)
          h_N(i,j+1,k) = h_in(i,j+1,k)
        enddo
      endif
    enddo
  endif

  if (monotonic) then
    ! untested
    call PPM_limit_cw84(bx, h_in_a, h_S_a, h_N_a)
  else
    call PPM_limit_pos(bx, h_in_a, h_S_a, h_N_a, h_min)
  endif

  !$omp target exit data map(release: slp)

  ! Deallocate local temporary array
  if(allocated(slp)) deallocate(slp)

  ! Deallocate local iteration boxes
  call bx%free()
  call bxE%free()

end subroutine PPM_reconstruction_y_fortran

!> This subroutine limits the left/right edge values of the PPM reconstruction
!! to give a reconstruction that is positive-definite.  Here this is
!! reinterpreted as giving a constant thickness if the mean thickness is less
!! than h_min, with a minimum of h_min otherwise.
subroutine PPM_limit_pos_fortran(bx, h_in_a, h_L_a, h_R_a, h_min)
  type(Box_t),        intent(in)    :: bx     !< Box over which to iterate
  type(RealArray_t),  intent(in)    :: h_in_a !< Layer thickness [H ~> m or kg m-2].
  type(RealArray_t),  intent(inout) :: h_L_a  !< Left thickness in the reconstruction [H ~> m or kg m-2].
  type(RealArray_t),  intent(inout) :: h_R_a  !< Right thickness in the reconstruction [H ~> m or kg m-2].
  real,               intent(in)    :: h_min  !< The minimum thickness

! Local variables
  real    :: curv  ! The grid-normalized curvature of the three thicknesses  [H ~> m or kg m-2]
  real    :: dh    ! The difference between the edge thicknesses             [H ~> m or kg m-2]
  real    :: scale ! A scaling factor to reduce the curvature of the fit               [nondim]
  integer :: i,j,k
  real, dimension(:,:,:), contiguous, pointer :: h_in, h_L, h_R  ! pointers to Fortran arrays

  ! Get the views
  call h_in_a%view(h_in)
  call h_L_a%view(h_L)
  call h_R_a%view(h_R)

  do concurrent (k=bx%idxS(3):bx%idxE(3), &
                 j=bx%idxS(2):bx%idxE(2), &
                 i=bx%idxS(1):bx%idxE(1))
    ! This limiter prevents undershooting minima within the domain with
    ! values less than h_min.
    curv = 3.0*((h_L(i,j,k) + h_R(i,j,k)) - 2.0*h_in(i,j,k))
    if (curv > 0.0) then ! Only minima are limited.
      dh = h_R(i,j,k) - h_L(i,j,k)
      if (abs(dh) < curv) then ! The parabola's minimum is within the cell.
        if (h_in(i,j,k) <= h_min) then
          h_L(i,j,k) = h_in(i,j,k) ; h_R(i,j,k) = h_in(i,j,k)
        elseif (12.0*curv*(h_in(i,j,k) - h_min) < (curv**2 + 3.0*dh**2)) then
          ! The minimum value is h_in - (curv^2 + 3*dh^2)/(12*curv), and must
          ! be limited in this case.  0 < scale < 1.
          scale = 12.0*curv*(h_in(i,j,k) - h_min) / (curv**2 + 3.0*dh**2)
          h_L(i,j,k) = h_in(i,j,k) + scale*(h_L(i,j,k) - h_in(i,j,k))
          h_R(i,j,k) = h_in(i,j,k) + scale*(h_R(i,j,k) - h_in(i,j,k))
        endif
      endif
    endif
  enddo

end subroutine PPM_limit_pos_fortran

!> This subroutine limits the left/right edge values of the PPM reconstruction
!! according to the monotonic prescription of Colella and Woodward, 1984.
subroutine PPM_limit_CW84_fortran(bx, h_in_a, h_L_a, h_R_a)
  type(box_t),         intent(in)   :: bx     !< Iteration box
  type(RealArray_t),  intent(in)    :: h_in_a !< Layer thickness [H ~> m or kg m-2].
  type(RealArray_t),  intent(inout) :: h_L_a  !< Left thickness in the reconstruction,
                                              !! [H ~> m or kg m-2].
  type(RealArray_t), intent(inout)   :: h_R_a !< Right thickness in the reconstruction,

  ! Local variables
  real    :: h_i      ! A copy of the cell-average layer thickness                [H ~> m or kg m-2]
  real    :: RLdiff   ! The difference between the input edge values              [H ~> m or kg m-2]
  real    :: RLdiff2  ! The squared difference between the input edge values   [H2 ~> m2 or kg2 m-4]
  real    :: RLmean   ! The average of the input edge thicknesses                 [H ~> m or kg m-2]
  real    :: FunFac   ! A curious product of the thickness slope and curvature [H2 ~> m2 or kg2 m-4]
  integer :: i, j, k
  real, dimension(:,:,:), contiguous, pointer :: h_in, h_L, h_R  ! pointers to Fortran arrays

  ! Get the views
  call h_in_a%view(h_in)
  call h_L_a%view(h_L)
  call h_R_a%view(h_R)

  ! untested
  do concurrent(k=bx%idxS(3):bx%idxE(3), &
                j=bx%idxS(2):bx%idxE(2), &
                i=bx%idxS(1):bx%idxE(1))
    ! This limiter monotonizes the parabola following
    ! Colella and Woodward, 1984, Eq. 1.10
    h_i = h_in(i,j,k)
    if ( ( h_R(i,j,k) - h_i ) * ( h_i - h_L(i,j,k) ) <= 0. ) then
      h_L(i,j,k) = h_i ; h_R(i,j,k) = h_i
    else
      RLdiff = h_R(i,j,k) - h_L(i,j,k)            ! Difference of edge values
      RLmean = 0.5 * ( h_R(i,j,k) + h_L(i,j,k) )  ! Mean of edge values
      FunFac = 6. * RLdiff * ( h_i - RLmean ) ! Some funny factor
      RLdiff2 = RLdiff * RLdiff               ! Square of difference
      if ( FunFac >  RLdiff2 ) h_L(i,j,k) = 3. * h_i - 2. * h_R(i,j,k)
      if ( FunFac < -RLdiff2 ) h_R(i,j,k) = 3. * h_i - 2. * h_L(i,j,k)
    endif
  enddo

  return
end subroutine PPM_limit_CW84_fortran

!> Return the maximum ratio of a/b or maxrat.
pure function ratio_max(a, b, maxrat) result(ratio)
  real, intent(in) :: a       !< Numerator, in arbitrary units [A]
  real, intent(in) :: b       !< Denominator, in arbitrary units [B]
  real, intent(in) :: maxrat  !< Maximum value of ratio [A B-1]
  real :: ratio               !< Return value [A B-1]

  if (abs(a) > abs(maxrat*b)) then
    ratio = maxrat
  else
    ratio = a / b
  endif
end function ratio_max

!> Initializes continuity_ppm_cs
subroutine continuity_PPM_init(Time, G, GV, US, param_file, diag, CS, OBC)
  type(time_type), target, intent(in)    :: Time !< The current model time.
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< Vertical grid structure.
  type(unit_scale_type),   intent(in)    :: US   !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file !< A structure indicating
                  !! the open file to parse for model parameter values.
  type(diag_ctrl), target, intent(inout) :: diag !< A structure that is used to
                  !! regulate diagnostic output.
  type(continuity_PPM_CS), intent(inout) :: CS   !< Module's control structure.
  type(ocean_OBC_type),    pointer       :: OBC  !< Open boundaries control structure.
  logical :: local_open_BC, use_h_marg_min
  type(OBC_segment_type), pointer :: segment => NULL()
  integer :: n

  !> This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_continuity_PPM" ! This module's name.
  character(len=256) :: mesg

  CS%initialized = .true.

  local_open_BC = .false.
  if (associated(OBC)) then
    local_open_BC = OBC%open_u_BCs_exist_globally
  endif

! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "MONOTONIC_CONTINUITY", CS%monotonic, &
                 "If true, CONTINUITY_PPM uses the Colella and Woodward "//&
                 "monotonic limiter.  The default (false) is to use a "//&
                 "simple positive definite limiter.", default=.false.)
  call get_param(param_file, mdl, "SIMPLE_2ND_PPM_CONTINUITY", CS%simple_2nd, &
                 "If true, CONTINUITY_PPM uses a simple 2nd order "//&
                 "(arithmetic mean) interpolation of the edge values. "//&
                 "This may give better PV conservation properties. While "//&
                 "it formally reduces the accuracy of the continuity "//&
                 "solver itself in the strongly advective limit, it does "//&
                 "not reduce the overall order of accuracy of the dynamic "//&
                 "core.", default=.false.)
  call get_param(param_file, mdl, "UPWIND_1ST_CONTINUITY", CS%upwind_1st, &
                 "If true, CONTINUITY_PPM becomes a 1st-order upwind "//&
                 "continuity solver.  This scheme is highly diffusive "//&
                 "but may be useful for debugging or in single-column "//&
                 "mode where its minimal stencil is useful.", default=.false.)
  call get_param(param_file, mdl, "ETA_TOLERANCE", CS%tol_eta, &
                 "The tolerance for the differences between the "//&
                 "barotropic and baroclinic estimates of the sea surface "//&
                 "height due to the fluxes through each face.  The total "//&
                 "tolerance for SSH is 4 times this value.  The default "//&
                 "is 0.5*NK*ANGSTROM, and this should not be set less "//&
                 "than about 10^-15*MAXIMUM_DEPTH.", units="m", scale=GV%m_to_H, &
                 default=0.5*GV%ke*GV%Angstrom_m)

  call get_param(param_file, mdl, "VELOCITY_TOLERANCE", CS%tol_vel, &
                 "The tolerance for barotropic velocity discrepancies "//&
                 "between the barotropic solution and  the sum of the "//&
                 "layer thicknesses.", units="m s-1", default=3.0e8, scale=US%m_s_to_L_T)
                 ! The speed of light is the default.

  call get_param(param_file, mdl, "CONT_PPM_AGGRESS_ADJUST", CS%aggress_adjust,&
                 "If true, allow the adjusted velocities to have a "//&
                 "relative CFL change up to 0.5.", default=.false.)
  CS%vol_CFL = CS%aggress_adjust
  call get_param(param_file, mdl, "CONT_PPM_VOLUME_BASED_CFL", CS%vol_CFL, &
                 "If true, use the ratio of the open face lengths to the "//&
                 "tracer cell areas when estimating CFL numbers.  The "//&
                 "default is set by CONT_PPM_AGGRESS_ADJUST.", &
                 default=CS%aggress_adjust, do_not_read=CS%aggress_adjust)
  call get_param(param_file, mdl, "CONTINUITY_CFL_LIMIT", CS%CFL_limit_adjust, &
                 "The maximum CFL of the adjusted velocities.", units="nondim", &
                 default=0.5)
  call get_param(param_file, mdl, "CONT_PPM_BETTER_ITER", CS%better_iter, &
                 "If true, stop corrective iterations using a velocity "//&
                 "based criterion and only stop if the iteration is "//&
                 "better than all predecessors.", default=.true.)
  call get_param(param_file, mdl, "CONT_PPM_USE_VISC_REM_MAX", CS%use_visc_rem_max, &
                 "If true, use more appropriate limiting bounds for "//&
                 "corrections in strongly viscous columns.", default=.true.)
  call get_param(param_file, mdl, "CONT_PPM_MARGINAL_FACE_AREAS", CS%marginal_faces, &
                 "If true, use the marginal face areas from the continuity "//&
                 "solver for use as the weights in the barotropic solver. "//&
                 "Otherwise use the transport averaged areas.", default=.true.)
  call get_param(param_file, mdl, "CONT_USE_H_MARG_MIN", use_h_marg_min, &
                 "If true, the marginal thickness used and returned from continuity "//&
                 "is bounded from below by a sub-roundoff value. Otherwise the "//&
                 "minimum is 0.", default=.false.)
  CS%diag => diag
  !$omp target update to(CS)

  id_clock_reconstruct = cpu_clock_id('(Ocean continuity reconstruction)', grain=CLOCK_ROUTINE)
  id_clock_update = cpu_clock_id('(Ocean continuity update)', grain=CLOCK_ROUTINE)
  id_clock_correct = cpu_clock_id('(Ocean continuity correction)', grain=CLOCK_ROUTINE)

  if (use_h_marg_min) then
    CS%h_marg_min = GV%H_subroundoff
  else
    CS%h_marg_min = 0.
  endif

  if (local_open_BC) then
    do n=1, OBC%number_of_segments
      segment => OBC%segment(n)
      if (associated(segment%h_Reg)) then
        if (.not. allocated(segment%h_Reg%h_res)) then
          write(mesg,'("In MOM_continuity_PPM, continuity_PPM_init called with ", &
                & "badly configured h_res.")')
          call MOM_error(FATAL, mesg)
        endif
      endif
    enddo
  endif

end subroutine continuity_PPM_init

!> continuity_PPM_stencil returns the continuity solver stencil size
function continuity_PPM_stencil(CS) result(stencil)
  type(continuity_PPM_CS), intent(in) :: CS   !< Module's control structure.
  integer ::  stencil !< The continuity solver stencil size with the current settings.

  stencil = 3 ; if (CS%simple_2nd) stencil = 2 ; if (CS%upwind_1st) stencil = 1

end function continuity_PPM_stencil

!> Set up a structure that stores the sizes of the i- and j-loops to work on in the continuity solver.
function set_continuity_loop_bounds(G, CS, i_stencil, j_stencil) result(LB)
  type(ocean_grid_type),   intent(in) :: G   !< The ocean's grid structure.
  type(continuity_PPM_CS), intent(in) :: CS  !< Module's control structure.
  logical,       optional, intent(in) :: i_stencil !< If present and true, extend the i-loop bounds
                                             !! by the stencil width of the continuity scheme.
  logical,       optional, intent(in) :: j_stencil !< If present and true, extend the j-loop bounds
                                             !! by the stencil width of the continuity scheme.
  type(cont_loop_bounds_type)         :: LB  !< A type storing the array sizes to work on in the continuity routines.

  ! Local variables
  logical :: add_i_stencil, add_j_stencil ! Local variables set based on i_stencil and j_stensil
  integer :: stencil    ! The continuity solver stencil size with the current continuity scheme.

  add_i_stencil = .false. ; if (present(i_stencil)) add_i_stencil = i_stencil
  add_j_stencil = .false. ; if (present(j_stencil)) add_j_stencil = j_stencil

  stencil = continuity_PPM_stencil(CS)

  if (add_i_stencil) then
    LB%ish = G%isc-stencil ; LB%ieh = G%iec+stencil
  else
    LB%ish = G%isc ; LB%ieh = G%iec
  endif

  if (add_j_stencil) then
    LB%jsh = G%jsc-stencil ; LB%jeh = G%jec+stencil
  else
    LB%jsh = G%jsc ; LB%jeh = G%jec
  endif

end function set_continuity_loop_bounds

!> Set up a structure that stores the sizes of the i- and j-loops to work on in the continuity solver.
function set_continuity_box(G, GV, CS, i_stencil, j_stencil) result(box)
  type(ocean_grid_type),   intent(in) :: G   !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV  !< Vertical grid structure.

  type(continuity_PPM_CS), intent(in) :: CS  !< Module's control structure.
  logical,       optional, intent(in) :: i_stencil !< If present and true, extend the i-loop bounds
                                             !! by the stencil width of the continuity scheme.
  logical,       optional, intent(in) :: j_stencil !< If present and true, extend the j-loop bounds
                                             !! by the stencil width of the continuity scheme.
  type(box_t) :: box                         !< The iteration box

  ! Local variables
  logical :: add_i_stencil, add_j_stencil ! Local variables set based on i_stencil and j_stensil
  integer :: stencil    ! The continuity solver stencil size with the current continuity scheme.
  integer :: is, ie, js, je

  add_i_stencil = .false. ; if (present(i_stencil)) add_i_stencil = i_stencil
  add_j_stencil = .false. ; if (present(j_stencil)) add_j_stencil = j_stencil

  stencil = continuity_PPM_stencil(CS)

  ! Allocate a 3 dimension iteration box
  call box%safe_alloc(ndims=3)

  if (add_i_stencil) then
    is = G%isc-stencil ; ie = G%iec+stencil
  else
    is = G%isc ; ie = G%iec
  endif

  if (add_j_stencil) then
    js = G%jsc-stencil ; je = G%jec+stencil
  else
    js = G%jsc ; je = G%jec
  endif

  ! Set the extents of the iteration space
  call box%set(idxS=[is,js,1],idxE=[ie,je,GV%ke])

end function set_continuity_box

!< shim for PPM_limit_pos
subroutine PPM_limit_pos(bx, h_in, h_L, h_R, h_min)
    implicit none

    type(box_t), intent(in)          :: bx   !< Box over which to iterate
    type(RealArray_t), intent(in)    :: h_in !< Layer thickness [H ~> m or kg m-2].
    type(RealArray_t), intent(inout) :: h_L  !< Left thickness in the reconstruction
                                             !! [H ~> m or kg m-2].
    type(RealArray_t), intent(inout) :: h_R  !< Right thickness in the reconstruction
                                             !! [H ~> m or kg m-2].
    real, intent(in)    :: h_min             !< The minimum thickness that can be obtain by a
                                             !! concave parabolic fit [H ~> m or kg m-2]

    ! local variables
    integer :: mode
    type(RealArray_C) :: h_in_c, h_L_c, h_R_c
    type(Box_c) :: bx_c
    integer :: rc
    type (io_recorder) :: rec
    logical :: capture
    character(len=80)  :: kernel
    character(len=100) :: dir
    character(len=256) :: binFile, metaFile

    kernel = "ppm_limit_pos"

    mode = getenv_mode("PPM_LIMIT_POS_MODE",default=TIMH_runFORTRAN)

    select case (mode)
       case (TIMH_capture)
           capture = .FALSE.
           if((.not. already_recorded(TRIM(kernel))) .and. is_root_pe()) capture = .TRUE.

           if(capture) then
             ! -----------WRITE DATA---------------------
             ! open a dump file subroutine arguments
             dir = "capture"
             rc = mkdir_posix(TRIM(dir) // c_null_char, int(o'755', c_int))

             binFile  = TRIM(dir) // "/" // TRIM(kernel) // ".bin"
             metaFile = TRIM(dir) // "/" // TRIM(kernel) // ".meta"

             call rec%open_write(binFile, metaFile)

             ! write out the input arguments
             call rec%add("_bx", bx)
             call rec%add("_h_in", h_in )
             call rec%add("_h_L_before", h_L )
             call rec%add("_h_R_before", h_R )
             call rec%add("_h_min", h_min)
           endif

           ! Run Fortran truth
           call ppm_limit_pos_fortran(bx,h_in, h_L, h_R, h_min)

           if(capture) then
             ! Write out the output arguments
             call rec%add("_h_L_after", h_L)
             call rec%add("_h_R_after", h_R)
             ! Close the file
             call rec%close()
             call mark_recorded(TRIM(kernel))
           endif
#ifdef _TIM
       case (TIMH_runAMREX)
           ! create C-compatible descriptors
           bx_c = bx%to_c(); h_in_c = h_in%to_c(); h_L_c  = h_L%to_c(); h_R_c  = h_R%to_c()
           ! Call C++ bridge to execute AMReX code
           call turbotmp_ppm_limit_pos_bridge(bx_c, h_in_c, h_L_c, h_R_c, h_min)
#endif
       case default
          ! Run Fortran code
          call ppm_limit_pos_fortran(bx,h_in, h_L, h_R, h_min)

    end select

end subroutine PPM_limit_pos

!< shim for PPM_limit_cw84
subroutine PPM_limit_cw84(bx, h_in, h_L, h_R)
    implicit none

    type(Box_t), intent(in)          :: bx   !< Box over which to iterate
    type(RealArray_t), intent(in)    :: h_in !< Layer thickness [H ~> m or kg m-2].
    type(RealArray_t), intent(inout) :: h_L  !< Left thickness in the
                                             !! reconstruction [H ~> m or kg m-2].
    type(RealArray_t), intent(inout) :: h_R  !< Right thickness in the
                                             !! reconstruction [H ~> m or kg m-2].
    ! local variables
    integer :: mode
    type(RealArray_C) :: h_in_c, h_L_c, h_R_c
    type(Box_c) :: bx_c
    integer :: rc
    type (io_recorder) :: rec
    logical :: capture
    character(len=80)  :: kernel
    character(len=100) :: dir
    character(len=256) :: binFile, metaFile

    kernel = "ppm_limit_cw84"

    mode = getenv_mode("PPM_LIMIT_CW84_MODE", default=TIMH_runFORTRAN)
    ! Call C++ bridge
    select case (mode)
       case (TIMH_capture)
          capture = .FALSE.
          if((.not. already_recorded(TRIM(kernel))) .and. is_root_pe()) capture = .TRUE.

          if(capture) then
            ! -----------WRITE DATA---------------------
            ! open a dump file to capture arguments
            dir = "capture"
            rc = mkdir_posix(TRIM(dir) // c_null_char, int(o'755', c_int))

            binFile  = TRIM(dir) // "/" // TRIM(kernel) // ".bin"
            metaFile = TRIM(dir) // "/" // TRIM(kernel) // ".meta"

            call rec%open_write(binFile, metaFile)

            ! write out the input arguments
            call rec%add("_bx", bx)
            call rec%add("_h_in", h_in )
            call rec%add("_h_L_before", h_L )
            call rec%add("_h_R_before", h_R )
          endif

          ! Run Fortran truth
          call ppm_limit_cw84_fortran(bx, h_in, h_L, h_R)

          if(capture) then
            ! Write out the output arguments
            call rec%add("_h_L_after", h_L)
            call rec%add("_h_R_after", h_R)
            ! Close the file
            call rec%close()
            call mark_recorded(TRIM(kernel))
          endif
#ifdef _TIM
       case (TIMH_runAMREX)
          ! Create C compatable descriptors
          bx_c = bx%to_c(); h_in_c = h_in%to_c(); h_L_c  = h_L%to_c(); h_R_c  = h_R%to_c()
          !  Call C+ bridge to execute AMReX code
          call turbotmp_ppm_limit_cw84_bridge(bx_c, h_in_c, h_L_c, h_R_c)
#endif
       case default
          ! Run Fortran code
          call ppm_limit_cw84_fortran(bx, h_in, h_L, h_R)
     end select

end subroutine PPM_limit_cw84

!< shim for PPM_reconstruction_y
subroutine PPM_reconstruction_y(bxH, h_in_a, h_S_a, h_N_a, mask2dT_a, h_min, monotonic, simple_2nd, OBC)
    implicit none

    type(Box_t), intent(in)          :: bxH        !< H-grid iteration Box
    type(RealArray_t), intent(in)    :: h_in_a     !< Layer thickness
    type(RealArray_t), intent(inout) :: h_S_a      !< South edge thickness
    type(RealArray_t), intent(inout) :: h_N_a      !< North edge thickness
    type(RealArray_t), intent(in)    :: mask2dT_a  !< Mask (0 land, 1 ocean)
    real,            intent(in)      :: h_min      !< Minimum thickness
    logical,         intent(in)      :: monotonic  !< Use CW84 limiter
    logical,         intent(in)      :: simple_2nd !< Use simple 2nd order
    type(ocean_OBC_type), pointer    :: OBC        !< Open boundary control

    ! local variables
    integer :: mode
    type(Box_C) :: bx_c
    type(RealArray_C) :: h_in_c, h_S_c, h_N_c, mask2dT_c
    type(c_ptr) :: OBC_c
    logical(c_bool) :: monotonic_c, simple_2nd_c
    integer :: rc
    type (io_recorder) :: rec
    logical :: capture
    character(len=80)  :: kernel
    character(len=100) :: dir
    character(len=256) :: binFile, metaFile

    kernel="ppm_reconstruction_y"

    mode = getenv_mode("PPM_RECONSTRUCTION_Y_MODE", default=TIMH_runFORTRAN)

    ! Call C++ bridge
    select case (mode)

       case (TIMH_capture)
          capture = .FALSE.
          if((.not. already_recorded(TRIM(kernel))) .and. is_root_pe()) capture = .TRUE.

          if(capture) then
            ! -----------WRITE DATA---------------------
            ! open a file to capture arguments
            dir = "capture"
            rc = mkdir_posix(TRIM(dir) // c_null_char, int(o'755', c_int))

            binFile  = TRIM(dir) // "/" // TRIM(kernel) // ".bin"
            metaFile = TRIM(dir) // "/" // TRIM(kernel) // ".meta"

            call rec%open_write(binFile, metaFile)

            ! Write out the input arguments
            call rec%add("_bxH", bxH)
            call rec%add("_h_in", h_in_a )
            call rec%add("_h_S_before", h_S_a )
            call rec%add("_h_N_before", h_N_a )
            call rec%add("_mask2d_t", mask2DT_a )
            call rec%add("_h_min", h_min)
            call rec%add("_monotonic", monotonic)
            call rec%add("_simple_2nd", simple_2nd)
            !call rec%add("OBC", OBC)
          endif

          ! Capture the Fortran results
          call ppm_reconstruction_y_fortran(bxH, h_in_a, h_S_a, h_N_a, &
                         mask2dT_a, h_min, monotonic, simple_2nd, OBC)
          if(capture) then
            ! Write out the output arguments
            call rec%add("_h_S_after", h_S_a)
            call rec%add("_h_N_after", h_N_a)
            ! Close the file
            call rec%close()
            call mark_recorded(TRIM(kernel))
          endif
#ifdef _TIM
       case (TIMH_runAMREX)

          ! create C-compatible descriptors
          bx_c = bxH%to_c(); h_in_c = h_in_a%to_c(); h_S_c = h_S_a%to_c();
          h_N_c = h_N_a%to_c(); mask2dT_c = mask2dT_a%to_c()
          monotonic_c = monotonic; simple_2nd_c = simple_2nd
          if(associated(OBC)) then; OBC_c = c_loc(OBC); else; OBC_c = c_null_ptr; endif
          ! Call C++ bridge to execute AMReX code
          call turbotmp_ppm_reconstruction_y_bridge(bx_c, h_in_c, h_S_c, h_N_c, mask2dT_c, &
                  h_min, monotonic_c, simple_2nd_c,OBC_c)
#endif
       case default
          ! Run Fortran code
          call ppm_reconstruction_y_fortran(bxH, h_in_a, h_S_a, h_N_a, &
                           mask2dT_a, h_min, monotonic, simple_2nd, OBC)
    end select

end subroutine PPM_reconstruction_y

!< shim for PPM_reconstruction_x
subroutine PPM_reconstruction_x(bxH, h_in_a, h_W_a, h_E_a, mask2dT_a, h_min, monotonic, simple_2nd, OBC)
    implicit none

    type(Box_t), intent(in)          :: bxH        !< H-grid iteration Box
    type(RealArray_t), intent(in)    :: h_in_a     !< Layer thickness
    type(RealArray_t), intent(inout) :: h_W_a      !< West edge thickness
    type(RealArray_t), intent(inout) :: h_E_a      !< East edge thickness
    type(RealArray_t), intent(in)    :: mask2dT_a  !< Mask (0 land, 1 ocean)
    real,            intent(in)      :: h_min      !< Minimum thickness
    logical,         intent(in)      :: monotonic  !< Use CW84 limiter
    logical,         intent(in)      :: simple_2nd !< Use simple 2nd order
    type(ocean_OBC_type), pointer    :: OBC        !< Open boundary control

    ! local variables
    integer :: mode
    type(Box_C) :: bx_c
    type(RealArray_C) :: h_in_c, h_W_c, h_E_c, mask2dT_c
    type(c_ptr) :: OBC_c
    logical(c_bool) :: monotonic_c, simple_2nd_c
    integer :: rc
    type (io_recorder) :: rec
    logical :: capture
    character(len=80)  :: kernel
    character(len=100) :: dir
    character(len=256) :: binFile, metaFile

    kernel="ppm_reconstruction_x"

    mode = getenv_mode("PPM_RECONSTRUCTION_X_MODE", default=TIMH_runFORTRAN)

    select case (mode)

       case (TIMH_capture)
          capture = .FALSE.
          if((.not. already_recorded(TRIM(kernel))) .and. is_root_pe()) capture = .TRUE.

          if(capture) then
            dir = "capture"
            rc = mkdir_posix(TRIM(dir) // c_null_char, int(o'755', c_int))

            binFile  = TRIM(dir) // "/" // TRIM(kernel) // ".bin"
            metaFile = TRIM(dir) // "/" // TRIM(kernel) // ".meta"

            call rec%open_write(binFile, metaFile)

            call rec%add("_bxH", bxH)
            call rec%add("_h_in", h_in_a )
            call rec%add("_h_W_before", h_W_a )
            call rec%add("_h_E_before", h_E_a )
            call rec%add("_mask2d_t", mask2dT_a )
            call rec%add("_h_min", h_min)
            call rec%add("_monotonic", monotonic)
            call rec%add("_simple_2nd", simple_2nd)
            !call rec%add("OBC", OBC)
          endif

          call ppm_reconstruction_x_fortran(bxH, h_in_a, h_W_a, h_E_a, &
                         mask2dT_a, h_min, monotonic, simple_2nd, OBC)
          if(capture) then
            call rec%add("_h_W_after", h_W_a)
            call rec%add("_h_E_after", h_E_a)
            call rec%close()
            call mark_recorded(TRIM(kernel))
          endif
#ifdef _TIM
       case (TIMH_runAMREX)

          ! create C-compatible descriptors
          bx_c = bxH%to_c(); h_in_c = h_in_a%to_c(); h_W_c = h_W_a%to_c();
          h_E_c = h_E_a%to_c(); mask2dT_c = mask2dT_a%to_c()
          monotonic_c = monotonic; simple_2nd_c = simple_2nd
          if(associated(OBC)) then; OBC_c = c_loc(OBC); else; OBC_c = c_null_ptr; endif
          ! Call C++ bridge to execute AMReX code
          call turbotmp_ppm_reconstruction_x_bridge(bx_c, h_in_c, h_W_c, h_E_c, mask2dT_c, &
                  h_min, monotonic_c, simple_2nd_c, OBC_c)
#endif
       case default
          call ppm_reconstruction_x_fortran(bxH, h_in_a, h_W_a, h_E_a, &
                           mask2dT_a, h_min, monotonic, simple_2nd, OBC)
    end select

end subroutine PPM_reconstruction_x

!> \namespace mom_continuity_ppm
!!
!! This module contains the subroutines that advect layer
!! thickness.  The scheme here uses a Piecewise-Parabolic method with
!! a positive definite limiter.

end module MOM_continuity_PPM
