! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0
!!SKILLS: 0.3

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

use array_mod, only : RealArray_t, RealArray_c, LogicalArray_t, LogicalArray_c
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

  interface
    !> Bridge for the continuity_zonal_convergence subroutine
    subroutine turbotmp_continuity_zonal_convergence_bridge(bxC, h, uh, dt, IareaT, hin, h_min) bind(C)
      use iso_c_binding, only : c_double
      use array_mod, only : RealArray_c
      use box_mod,   only : Box_c
      implicit none

      type(Box_C),       intent(in)        :: bxC    !< Iteration box for continuity solver
      type(RealArray_C), intent(inout)     :: h      !< Final layer thickness [H ~> m or kg m-2]
      type(RealArray_C), intent(in)        :: uh     !< Zonal thickness flux, u*h*dy
                                                      !! [H L2 T-1 ~> m3 s-1 or kg s-1]
      real(c_double),    intent(in), value :: dt     !< Time increment [T ~> s]
      type(RealArray_C), intent(in)        :: IareaT !< The grid cell's 1/areaT [L-2 ~> m-2]
      type(RealArray_C), intent(in)        :: hin    !< Initial layer thickness [H ~> m or kg m-2]
      real(c_double),    intent(in), value :: h_min  !< The minimum layer thickness [H ~> m or kg m-2]
    end subroutine turbotmp_continuity_zonal_convergence_bridge
  end interface

  interface
    !> Bridge for the continuity_meridional_convergence subroutine
    subroutine turbotmp_continuity_meridional_convergence_bridge(bxC, h, vh, dt, IareaT, hin, &
                                                                  h_min) bind(C)
      use iso_c_binding, only : c_double
      use array_mod, only : RealArray_c
      use box_mod,   only : Box_c
      implicit none

      type(Box_C),       intent(in)        :: bxC    !< Iteration box for continuity solver
      type(RealArray_C), intent(inout)     :: h      !< Final layer thickness [H ~> m or kg m-2]
      type(RealArray_C), intent(in)        :: vh     !< Meridional thickness flux, v*h*dx
                                                      !! [H L2 T-1 ~> m3 s-1 or kg s-1]
      real(c_double),    intent(in), value :: dt     !< Time increment [T ~> s]
      type(RealArray_C), intent(in)        :: IareaT !< The grid cell's 1/areaT [L-2 ~> m-2]
      type(RealArray_C), intent(in)        :: hin    !< Initial layer thickness [H ~> m or kg m-2]
      real(c_double),    intent(in), value :: h_min  !< The minimum layer thickness [H ~> m or kg m-2]
    end subroutine turbotmp_continuity_meridional_convergence_bridge
  end interface


#include <MOM_memory.h>

public continuity_PPM, continuity_PPM_init, continuity_PPM_stencil
public continuity_PPM_3d_fluxes, continuity_PPM_2d_fluxes, continuity_PPM_adjust_vel
public zonal_mass_flux, meridional_mass_flux
public zonal_edge_thickness, meridional_edge_thickness
public continuity_zonal_convergence, continuity_meridional_convergence
public zonal_flux_thickness, meridional_flux_thickness
public zonal_BT_mass_flux, meridional_BT_mass_flux
public set_continuity_loop_bounds

!>@{ CPU time clock IDs
integer :: id_clock_reconstruct, id_clock_update, id_clock_correct
!>@}
!> Options controlling the edge-value reconstruction scheme used by the continuity solver.
type, public :: reconstruction_CS
  logical :: upwind_1st      !< If true, use a first-order upwind scheme.
  logical :: monotonic       !< If true, use the Colella & Woodward monotonic
                             !! limiter; otherwise use a simple positive
                             !! definite limiter.
  logical :: simple_2nd      !< If true, use a simple second order (arithmetic
                             !! mean) interpolation of the edge values instead
                             !! of the higher order interpolation.
end type reconstruction_CS

!> bind(C) mirror of reconstruction_CS, field-for-field, same order.
type, bind(C) :: reconstruction_CS_C
  logical(c_bool) :: upwind_1st  !< If true, use a first-order upwind scheme.
  logical(c_bool) :: monotonic   !< If true, use the Colella & Woodward monotonic limiter.
  logical(c_bool) :: simple_2nd  !< If true, use a simple second order interpolation.
end type reconstruction_CS_C

!> Options controlling the transport adjustment and barotropic-consistency
!! iteration used by the continuity solver.
type, public :: transport_adjust_CS
  real :: tol_eta            !< The tolerance for free-surface height
                             !! discrepancies between the barotropic solution and
                             !! the sum of the layer thicknesses [H ~> m or kg m-2].
  real :: tol_vel            !< The tolerance for barotropic velocity
                             !! discrepancies between the barotropic solution and
                             !! the sum of the layer thicknesses [L T-1 ~> m s-1].
  real :: CFL_limit_adjust   !< The maximum CFL of the adjusted velocities [nondim]
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
end type transport_adjust_CS

!> bind(C) mirror of transport_adjust_CS, field-for-field, same order.
type, bind(C) :: transport_adjust_CS_C
  real(c_double)  :: tol_eta            !< The tolerance for free-surface height discrepancies.
  real(c_double)  :: tol_vel            !< The tolerance for barotropic velocity discrepancies.
  real(c_double)  :: CFL_limit_adjust   !< The maximum CFL of the adjusted velocities.
  logical(c_bool) :: aggress_adjust     !< If true, allow a larger relative CFL change.
  logical(c_bool) :: vol_CFL            !< If true, use the ratio of open face lengths to
                                        !! tracer cell areas when estimating CFL numbers.
  logical(c_bool) :: better_iter        !< If true, use a velocity-based iteration criterion.
  logical(c_bool) :: use_visc_rem_max   !< If true, use limiting bounds for viscous columns.
  logical(c_bool) :: marginal_faces     !< If true, use marginal face areas as barotropic weights.
end type transport_adjust_CS_C

interface
  !> Bridge for the zonal_BT_mass_flux subroutine
  subroutine turbotmp_zonal_bt_mass_flux_bridge(bxC, u, h_in, h_W, h_E, uhbt, dt, dy_Cu, &
                                                IareaT, IdxT, CS, obc, por_face_areaU) bind(C)
    use iso_c_binding, only : c_double, c_ptr
    use array_mod, only : RealArray_c
    use box_mod,   only : Box_c
    import :: transport_adjust_CS_C
    implicit none

    type(Box_C), intent(in) :: bxC !< Iteration box for continuity solver
    type(RealArray_C), intent(in) :: u !< Zonal velocity [L T-1 ~> m s-1]
    type(RealArray_C), intent(in) :: h_in !< Layer thickness used to
                                         !! calculate fluxes [H ~> m or kg m-2]
    type(RealArray_C), intent(in) :: h_W !< Western edge thickness in the PPM
                                        !! reconstruction [H ~> m or kg m-2].
    type(RealArray_C), intent(in) :: h_E !< Eastern edge thickness in the PPM
                                        !! reconstruction [H ~> m or kg m-2].
    type(RealArray_C), intent(inout) :: uhbt !< The summed volume flux through
                                            !! zonal faces [H L2 T-1 ~> m3 s-1 or kg s-1].
    real(c_double), intent(in), value :: dt !< Time increment [T ~> s].
    type(RealArray_C), intent(in) :: dy_Cu !< The grid cell's unblocked lengths of
                                          !! the u-faces of the h-cell [L ~> m].
    type(RealArray_C), intent(in) :: IareaT !< The grid cell's 1/areaT [L-2 ~> m-2].
    type(RealArray_C), intent(in) :: IdxT !< The grid cell's 1/dxT [L-1 ~> m-1].
    type(transport_adjust_CS_C), intent(in) :: CS !< Options
                       !! controlling the transport adjustment and barotropic-consistency iteration.
    type(c_ptr), intent(in), value :: obc !< Open boundary condition type
                                         !! specifies whether, where, and what
                                         !! open boundary conditions are used.
    type(RealArray_C), intent(in) :: por_face_areaU !< fractional open area of U-faces [nondim]
  end subroutine turbotmp_zonal_bt_mass_flux_bridge
end interface

interface
  !> Bridge for the meridional_BT_mass_flux subroutine
  subroutine turbotmp_meridional_bt_mass_flux_bridge(bxC, v, h_in, h_S, h_N, vhbt, dt, dx_Cv, &
                                                     IareaT, IdyT, CS, obc, por_face_areaV) bind(C)
    use iso_c_binding, only : c_double, c_ptr
    use array_mod, only : RealArray_c
    use box_mod,   only : Box_c
    import :: transport_adjust_CS_C
    implicit none

    type(Box_C), intent(in) :: bxC !< Iteration box for continuity solver
    type(RealArray_C), intent(in) :: v !< Meridional velocity [L T-1 ~> m s-1]
    type(RealArray_C), intent(in) :: h_in !< Layer thickness used to
                                         !! calculate fluxes [H ~> m or kg m-2]
    type(RealArray_C), intent(in) :: h_S !< Southern edge thickness in the PPM
                                        !! reconstruction [H ~> m or kg m-2].
    type(RealArray_C), intent(in) :: h_N !< Northern edge thickness in the PPM
                                        !! reconstruction [H ~> m or kg m-2].
    type(RealArray_C), intent(inout) :: vhbt !< The summed volume flux through
                                            !! meridional faces [H L2 T-1 ~> m3 s-1 or kg s-1].
    real(c_double), intent(in), value :: dt !< Time increment [T ~> s].
    type(RealArray_C), intent(in) :: dx_Cv !< The grid cell's unblocked lengths of
                                          !! the v-faces of the h-cell [L ~> m].
    type(RealArray_C), intent(in) :: IareaT !< The grid cell's 1/areaT [L-2 ~> m-2].
    type(RealArray_C), intent(in) :: IdyT !< The grid cell's 1/dyT [L-1 ~> m-1].
    type(transport_adjust_CS_C), intent(in) :: CS !< Options
                       !! controlling the transport adjustment and barotropic-consistency iteration.
    type(c_ptr), intent(in), value :: obc !< Open boundary condition type
                                         !! specifies whether, where, and what
                                         !! open boundary conditions are used.
    type(RealArray_C), intent(in) :: por_face_areaV !< fractional open area of V-faces [nondim]
  end subroutine turbotmp_meridional_bt_mass_flux_bridge
end interface

interface
  !> Bridge for the zonal_flux_thickness subroutine
  subroutine turbotmp_zonal_flux_thickness_bridge(bxC, u, h, h_W, h_E, h_u, dt, dy_Cu, IareaT, &
                                                  IdxT, vol_CFL, marginal, obc, por_face_areaU, &
                                                  visc_rem_u) bind(C)
    use iso_c_binding, only : c_double, c_bool, c_ptr
    use array_mod, only : RealArray_c
    use box_mod,   only : Box_c
    implicit none

    type(Box_C), intent(in) :: bxC !< Iteration box for continuity solver
    type(RealArray_C), intent(in) :: u !< Zonal velocity [L T-1 ~> m s-1].
    type(RealArray_C), intent(in) :: h !< Layer thickness used to
                                      !! calculate fluxes [H ~> m or kg m-2].
    type(RealArray_C), intent(in) :: h_W !< West edge thickness in the
                                        !! reconstruction [H ~> m or kg m-2].
    type(RealArray_C), intent(in) :: h_E !< East edge thickness in the
                                        !! reconstruction [H ~> m or kg m-2].
    type(RealArray_C), intent(inout) :: h_u !< Effective thickness at zonal faces,
                                           !! scaled down to account for the effects of
                                           !! viscosity and the fractional open area
                                           !! [H ~> m or kg m-2].
    real(c_double), intent(in), value :: dt !< Time increment [T ~> s].
    type(RealArray_C), intent(in) :: dy_Cu !< The grid cell's unblocked lengths of
                                          !! the u-faces of the h-cell [L ~> m].
    type(RealArray_C), intent(in) :: IareaT !< The grid cell's 1/areaT [L-2 ~> m-2].
    type(RealArray_C), intent(in) :: IdxT !< The grid cell's 1/dxT [L-1 ~> m-1].
    logical(c_bool), intent(in), value :: vol_CFL !< If true, rescale the ratio
                          !! of face areas to the cell areas when estimating the CFL number.
    logical(c_bool), intent(in), value :: marginal !< If true, report the
                          !! marginal face thicknesses; otherwise report transport-averaged thicknesses.
    type(c_ptr), intent(in), value :: obc !< Open boundaries control structure.
    type(RealArray_C), intent(in) :: por_face_areaU !< fractional open area of
                                                    !! U-faces [nondim]
    type(RealArray_C), intent(in) :: visc_rem_u
                          !< Both the fraction of the momentum originally in a layer that remains after
                          !! a time-step of viscosity, and the fraction of a time-step's worth of a
                          !! barotropic acceleration that a layer experiences after viscosity is applied [nondim].
                          !! Visc_rem_u is between 0 (at the bottom) and 1 (far above the bottom).
  end subroutine turbotmp_zonal_flux_thickness_bridge
end interface

interface
  !> Bridge for the meridional_flux_thickness subroutine
  subroutine turbotmp_meridional_flux_thickness_bridge(bxC, v, h, h_S, h_N, h_v, dt, &
                                                        dx_Cv, IareaT, IdyT, vol_CFL, marginal, &
                                                        obc, por_face_areaV, visc_rem_v) bind(C)
    use iso_c_binding, only : c_double, c_bool, c_ptr
    use array_mod, only : RealArray_c
    use box_mod,   only : Box_c
    implicit none

    type(Box_C), intent(in) :: bxC !< Iteration box for continuity solver
    type(RealArray_C), intent(in) :: v !< Meridional velocity [L T-1 ~> m s-1].
    type(RealArray_C), intent(in) :: h !< Layer thickness used to
                                      !! calculate fluxes, [H ~> m or kg m-2].
    type(RealArray_C), intent(in) :: h_S !< South edge thickness in the
                                        !! reconstruction, [H ~> m or kg m-2].
    type(RealArray_C), intent(in) :: h_N !< North edge thickness in the
                                        !! reconstruction, [H ~> m or kg m-2].
    type(RealArray_C), intent(inout) :: h_v !< Effective thickness at meridional faces,
                                           !! scaled down to account for the effects of
                                           !! viscosity and the fractional open area
                                           !! [H ~> m or kg m-2].
    real(c_double), intent(in), value :: dt !< Time increment [T ~> s].
    type(RealArray_C), intent(in) :: dx_Cv !< The grid cell's unblocked lengths of
                                          !! the v-faces of the h-cell [L ~> m].
    type(RealArray_C), intent(in) :: IareaT !< The grid cell's 1/areaT [L-2 ~> m-2].
    type(RealArray_C), intent(in) :: IdyT !< The grid cell's 1/dyT [L-1 ~> m-1].
    logical(c_bool), intent(in), value :: vol_CFL !< If true, rescale the ratio
                          !! of face areas to the cell areas when estimating the CFL number.
    logical(c_bool), intent(in), value :: marginal !< If true, report the marginal
                          !! face thicknesses; otherwise report transport-averaged thicknesses.
    type(c_ptr), intent(in), value :: obc !< Open boundaries control structure.
    type(RealArray_C), intent(in) :: por_face_areaV !< fractional open area of
                                                    !! V-faces [nondim]
    type(RealArray_C), intent(in) :: visc_rem_v !< Both the fraction
                          !! of the momentum originally in a layer that remains after a time-step of
                          !! viscosity, and the fraction of a time-step's worth of a barotropic
                          !! acceleration that a layer experiences after viscosity is applied [nondim].
                          !! Visc_rem_v is between 0 (at the bottom) and 1 (far above the bottom).
  end subroutine turbotmp_meridional_flux_thickness_bridge
end interface

interface
  !> Bridge for the zonal_flux_adjust subroutine
  subroutine turbotmp_zonal_flux_adjust_bridge(bxC, u, h_in, h_W, h_E, uh_tot_0, duhdu_tot_0, du, &
                                               du_max_CFL, du_min_CFL, dt, dy_Cu, IareaT, IdxT, &
                                               CS, visc_rem, do_I_in, por_face_areaU, uhbt, &
                                               uh_3d, obc) bind(C)
    use iso_c_binding, only : c_double, c_ptr
    use array_mod, only : RealArray_c, LogicalArray_c
    use box_mod,   only : Box_c
    import :: transport_adjust_CS_C
    implicit none

    type(Box_C), intent(in) :: bxC !< Iteration box for continuity solver
    type(RealArray_C), intent(in) :: u !< Zonal velocity [L T-1 ~> m s-1].
    type(RealArray_C), intent(in) :: h_in !< Layer thickness used to
                                         !! calculate fluxes [H ~> m or kg m-2].
    type(RealArray_C), intent(in) :: h_W !< West edge thickness in the
                                        !! reconstruction [H ~> m or kg m-2].
    type(RealArray_C), intent(in) :: h_E !< East edge thickness in the
                                        !! reconstruction [H ~> m or kg m-2].
    type(RealArray_C), intent(in) :: uh_tot_0 !< The summed transport
                       !! with 0 adjustment [H L2 T-1 ~> m3 s-1 or kg s-1].
    type(RealArray_C), intent(in) :: duhdu_tot_0 !< The partial derivative
                       !! of du_err with du at 0 adjustment [H L ~> m2 or kg m-1].
    type(RealArray_C), intent(inout) :: du !< The barotropic velocity adjustment [L T-1 ~> m s-1].
    type(RealArray_C), intent(in) :: du_max_CFL !< Maximum acceptable
                       !! value of du [L T-1 ~> m s-1].
    type(RealArray_C), intent(in) :: du_min_CFL !< Minimum acceptable
                       !! value of du [L T-1 ~> m s-1].
    real(c_double), intent(in), value :: dt !< Time increment [T ~> s].
    type(RealArray_C), intent(in) :: dy_Cu !< The grid cell's unblocked
                       !! lengths of the u-faces of the h-cell [L ~> m].
    type(RealArray_C), intent(in) :: IareaT !< The grid cell's 1/areaT
                       !! [L-2 ~> m-2].
    type(RealArray_C), intent(in) :: IdxT !< The grid cell's 1/dxT
                       !! [L-1 ~> m-1].
    type(transport_adjust_CS_C), intent(in) :: CS !< Options
                       !! controlling the transport adjustment and barotropic-consistency iteration.
    type(RealArray_C), intent(in) :: visc_rem !< Both the fraction of the
                       !! momentum originally in a layer that remains after a time-step of viscosity, and
                       !! the fraction of a time-step's worth of a barotropic acceleration that a layer
                       !! experiences after viscosity is applied [nondim].
                       !! Visc_rem is between 0 (at the bottom) and 1 (far above the bottom).
    type(LogicalArray_C), intent(in) :: do_I_in !< A logical flag indicating
                                                                       !! which I values to work on.
    type(RealArray_C), intent(in) :: por_face_areaU !< fractional open area
                                                                              !! of U-faces [nondim].
    type(RealArray_C), intent(in) :: uhbt !< The summed volume flux
                       !! through zonal faces [H L2 T-1 ~> m3 s-1 or kg s-1].
    type(RealArray_C), intent(inout) :: uh_3d !< Volume flux through zonal
                                                 !! faces = u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1].
    type(c_ptr), intent(in), value :: obc !< Open boundaries control structure.
  end subroutine turbotmp_zonal_flux_adjust_bridge
end interface

interface
  !> Bridge for the meridional_flux_adjust subroutine
  subroutine turbotmp_meridional_flux_adjust_bridge(bxC, v, h_in, h_S, h_N, vh_tot_0, &
                                                    dvhdv_tot_0, dv, dv_max_CFL, dv_min_CFL, dt, &
                                                    dx_Cv, IareaT, IdyT, CS, visc_rem, do_I_in, &
                                                    por_face_areaV, vhbt, vh_3d, obc) bind(C)
    use iso_c_binding, only : c_double, c_ptr
    use array_mod, only : RealArray_c, LogicalArray_c
    use box_mod,   only : Box_c
    import :: transport_adjust_CS_C
    implicit none

    type(Box_C), intent(in) :: bxC !< Iteration box for continuity solver
    type(RealArray_C), intent(in) :: v !< Meridional velocity [L T-1 ~> m s-1].
    type(RealArray_C), intent(in) :: h_in !< Layer thickness used to calculate
                                          !! fluxes [H ~> m or kg m-2].
    type(RealArray_C), intent(in) :: h_S !< South edge thickness in the
                                        !! reconstruction [H ~> m or kg m-2].
    type(RealArray_C), intent(in) :: h_N !< North edge thickness in the
                                        !! reconstruction [H ~> m or kg m-2].
    type(RealArray_C), intent(in) :: vh_tot_0 !< The summed transport with 0 adjustment
                                              !! [H L2 T-1 ~> m3 s-1 or kg s-1].
    type(RealArray_C), intent(in) :: dvhdv_tot_0 !< The partial derivative of dv_err with
                                                 !! dv at 0 adjustment [H L ~> m2 or kg m-1].
    type(RealArray_C), intent(inout) :: dv !< The barotropic velocity adjustment [L T-1 ~> m s-1].
    type(RealArray_C), intent(in) :: dv_max_CFL !< Maximum acceptable value of dv [L T-1 ~> m s-1].
    type(RealArray_C), intent(in) :: dv_min_CFL !< Minimum acceptable value of dv [L T-1 ~> m s-1].
    real(c_double), intent(in), value :: dt !< Time increment [T ~> s].
    type(RealArray_C), intent(in) :: dx_Cv !< The grid cell's unblocked lengths of the
                                           !! v-faces of the h-cell [L ~> m].
    type(RealArray_C), intent(in) :: IareaT !< The grid cell's 1/areaT [L-2 ~> m-2].
    type(RealArray_C), intent(in) :: IdyT !< The grid cell's 1/dyT [L-1 ~> m-1].
    type(transport_adjust_CS_C), intent(in) :: CS !< Options controlling the
                       !! transport adjustment and barotropic-consistency iteration.
    type(RealArray_C), intent(in) :: visc_rem
                             !< Both the fraction of the momentum originally
                             !! in a layer that remains after a time-step of viscosity, and the
                             !! fraction of a time-step's worth of a barotropic acceleration that
                             !! a layer experiences after viscosity is applied [nondim].
                             !! Visc_rem is between 0 (at the bottom) and 1 (far above the bottom).
    type(LogicalArray_C), intent(in) :: do_I_in !< A flag indicating which I values to work on.
    type(RealArray_C), intent(in) :: por_face_areaV !< fractional open area of V-faces [nondim]
    type(RealArray_C), intent(in) :: vhbt !< The summed volume flux through
                                         !! meridional faces [H L2 T-1 ~> m3 s-1 or kg s-1].
    type(RealArray_C), intent(inout) :: vh_3d !< Volume flux through meridional
                       !! faces = v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1].
    type(c_ptr), intent(in), value :: obc !< Pointer to OBC structure.
  end subroutine turbotmp_meridional_flux_adjust_bridge
end interface

interface
  !> Bridge for the set_zonal_BT_cont subroutine
  subroutine turbotmp_set_zonal_bt_cont_bridge(bxC, u, h_in, h_W, h_E, FA_u_W0, FA_u_E0, &
                                               FA_u_WW, FA_u_EE, uBT_WW, uBT_EE, du0, &
                                               uh_tot_0, duhdu_tot_0, du_max_CFL, du_min_CFL, &
                                               dt, dxCu, dy_Cu, IareaT, IdxT, CS, visc_rem, &
                                               visc_rem_max, do_I, por_face_areaU) bind(C)
    use iso_c_binding, only : c_double
    use array_mod, only : RealArray_c, LogicalArray_c
    use box_mod,   only : Box_c
    import :: transport_adjust_CS_C
    implicit none

    type(Box_C), intent(in) :: bxC !< Iteration box for continuity solver
    type(RealArray_C), intent(in) :: u !< Zonal velocity [L T-1 ~> m s-1].
    type(RealArray_C), intent(in) :: h_in !< Layer thickness used to calculate
                                         !! fluxes [H ~> m or kg m-2].
    type(RealArray_C), intent(in) :: h_W !< West edge thickness in the reconstruction
                                        !! [H ~> m or kg m-2].
    type(RealArray_C), intent(in) :: h_E !< East edge thickness in the reconstruction
                                        !! [H ~> m or kg m-2].
    type(RealArray_C), intent(inout) :: FA_u_W0 !< Effective open face area, west, 0 transport
    type(RealArray_C), intent(inout) :: FA_u_E0 !< Effective open face area, east, 0 transport
    type(RealArray_C), intent(inout) :: FA_u_WW !< Effective open face area, westerly test velocity
    type(RealArray_C), intent(inout) :: FA_u_EE !< Effective open face area, easterly test velocity
    type(RealArray_C), intent(inout) :: uBT_WW  !< Westerly correction to the barotropic velocity
    type(RealArray_C), intent(inout) :: uBT_EE  !< Easterly correction to the barotropic velocity
    type(RealArray_C), intent(in) :: du0 !< The barotropic velocity increment that gives 0
                                        !! transport [L T-1 ~> m s-1].
    type(RealArray_C), intent(in) :: uh_tot_0 !< The summed transport with 0 adjustment
                                             !! [H L2 T-1 ~> m3 s-1 or kg s-1].
    type(RealArray_C), intent(in) :: duhdu_tot_0 !< The partial derivative
                       !! of du_err with du at 0 adjustment [H L ~> m2 or kg m-1].
    type(RealArray_C), intent(in) :: du_max_CFL !< Maximum acceptable value of
                                                !! du [L T-1 ~> m s-1].
    type(RealArray_C), intent(in) :: du_min_CFL !< Minimum acceptable value of
                                                !! du [L T-1 ~> m s-1].
    real(c_double), intent(in), value :: dt !< Time increment [T ~> s].
    type(RealArray_C), intent(in) :: dxCu !< The grid cell's u-point x-extent [L ~> m].
    type(RealArray_C), intent(in) :: dy_Cu !< The grid cell's unblocked lengths of the
                                          !! u-faces of the h-cell [L ~> m].
    type(RealArray_C), intent(in) :: IareaT !< The grid cell's 1/areaT [L-2 ~> m-2].
    type(RealArray_C), intent(in) :: IdxT !< The grid cell's 1/dxT [L-1 ~> m-1].
    type(transport_adjust_CS_C), intent(in) :: CS !< Options controlling the
                       !! transport adjustment and barotropic-consistency iteration.
    type(RealArray_C), intent(in) :: visc_rem !< Both the fraction of the
                       !! momentum originally in a layer that remains after a time-step of viscosity, and
                       !! the fraction of a time-step's worth of a barotropic acceleration that a layer
                       !! experiences after viscosity is applied [nondim].
                       !! Visc_rem is between 0 (at the bottom) and 1 (far above the bottom).
    type(RealArray_C), intent(in) :: visc_rem_max !< Maximum allowable visc_rem [nondim].
    type(LogicalArray_C), intent(in) :: do_I !< A logical flag indicating
                                                        !! which I values to work on.
    type(RealArray_C), intent(in) :: por_face_areaU !< fractional open area
                                                        !! of U-faces [nondim]
  end subroutine turbotmp_set_zonal_bt_cont_bridge
end interface

interface
  !> Bridge for the set_merid_BT_cont subroutine
  subroutine turbotmp_set_merid_bt_cont_bridge(bxC, v, h_in, h_S, h_N, FA_v_S0, FA_v_N0, &
                                               FA_v_SS, FA_v_NN, vBT_SS, vBT_NN, dv0, &
                                               vh_tot_0, dvhdv_tot_0, dv_max_CFL, dv_min_CFL, &
                                               dt, dyCv, dx_Cv, IareaT, IdyT, CS, visc_rem, &
                                               visc_rem_max, do_I, por_face_areaV) bind(C)
    use iso_c_binding, only : c_double
    use array_mod, only : RealArray_c, LogicalArray_c
    use box_mod,   only : Box_c
    import :: transport_adjust_CS_C
    implicit none

    type(Box_C), intent(in) :: bxC !< Iteration box for continuity solver
    type(RealArray_C), intent(in) :: v !< Meridional velocity [L T-1 ~> m s-1].
    type(RealArray_C), intent(in) :: h_in !< Layer thickness used to
                                         !! calculate fluxes, [H ~> m or kg m-2].
    type(RealArray_C), intent(in) :: h_S !< South edge thickness in the
                                        !! reconstruction, [H ~> m or kg m-2].
    type(RealArray_C), intent(in) :: h_N !< North edge thickness in the
                                        !! reconstruction, [H ~> m or kg m-2].
    type(RealArray_C), intent(inout) :: FA_v_S0 !< Effective open face area, south, 0 transport
    type(RealArray_C), intent(inout) :: FA_v_N0 !< Effective open face area, north, 0 transport
    type(RealArray_C), intent(inout) :: FA_v_SS !< Effective open face area, southerly test velocity
    type(RealArray_C), intent(inout) :: FA_v_NN !< Effective open face area, northerly test velocity
    type(RealArray_C), intent(inout) :: vBT_SS  !< Southerly correction to the barotropic velocity
    type(RealArray_C), intent(inout) :: vBT_NN  !< Northerly correction to the barotropic velocity
    type(RealArray_C), intent(in) :: dv0 !< The barotropic velocity increment
                                        !! that gives 0 transport [L T-1 ~> m s-1].
    type(RealArray_C), intent(in) :: vh_tot_0 !< The summed transport
                       !! with 0 adjustment [H L2 T-1 ~> m3 s-1 or kg s-1].
    type(RealArray_C), intent(in) :: dvhdv_tot_0 !< The partial derivative
                       !! of du_err with dv at 0 adjustment [H L ~> m2 or kg m-1].
    type(RealArray_C), intent(in) :: dv_max_CFL !< Maximum acceptable value
                                                                        !!  of dv [L T-1 ~> m s-1].
    type(RealArray_C), intent(in) :: dv_min_CFL !< Minimum acceptable value
                                                                        !!  of dv [L T-1 ~> m s-1].
    real(c_double), intent(in), value :: dt !< Time increment [T ~> s].
    type(RealArray_C), intent(in) :: dyCv !< The grid cell's v-point y-extent [L ~> m].
    type(RealArray_C), intent(in) :: dx_Cv !< The grid cell's unblocked lengths of the
                                          !! v-faces of the h-cell [L ~> m].
    type(RealArray_C), intent(in) :: IareaT !< The grid cell's 1/areaT [L-2 ~> m-2].
    type(RealArray_C), intent(in) :: IdyT !< The grid cell's 1/dyT [L-1 ~> m-1].
    type(transport_adjust_CS_C), intent(in) :: CS !< Options
                       !! controlling the transport adjustment and barotropic-consistency iteration.
    type(RealArray_C), intent(in) :: visc_rem !< Both the fraction of the
                       !! momentum originally in a layer that remains after a time-step
                       !! of viscosity, and the fraction of a time-step's worth of a barotropic
                       !! acceleration that a layer experiences after viscosity is applied [nondim].
                       !! Visc_rem is between 0 (at the bottom) and 1 (far above the bottom).
    type(RealArray_C), intent(in) :: visc_rem_max !< Maximum allowable visc_rem [nondim]
    type(LogicalArray_C), intent(in) :: do_I !< A logical flag indicating
                                             !! which I values to work on.
    type(RealArray_C), intent(in) :: por_face_areaV !< fractional open
                                                    !! area of V-faces [nondim]
  end subroutine turbotmp_set_merid_bt_cont_bridge
end interface

!> Control structure for mom_continuity_ppm
type, public :: continuity_PPM_CS ; private
  logical :: initialized = .false. !< True if this control structure has been initialized.
  type(diag_ctrl), pointer :: diag !< Diagnostics control structure.
  type(reconstruction_CS) :: reconstruction_CS !< Options controlling the
                             !! edge-value reconstruction scheme.
  type(transport_adjust_CS) :: transport_adjust_CS !< Options controlling the
                             !! transport adjustment and barotropic-consistency iteration.
  real :: h_marg_min         !< Negligible floor on h_marg, the marginal thickness
                             !! used to calculate the partial derivative of transports
                             !! with velocities [H ~> m or kg m-2]
end type continuity_PPM_CS

!> A container for loop bounds
type, public :: cont_loop_bounds_type ; private
  !>@{ Loop bounds
  integer :: ish, ieh, jsh, jeh
  !>@}
end type cont_loop_bounds_type

contains

!> Converts a reconstruction_CS to its bind(C) mirror.
function reconstruction_CS_to_c(opts) result(cdesc)
  type(reconstruction_CS), intent(in) :: opts !< Options to convert
  type(reconstruction_CS_C) :: cdesc                !< bind(C) mirror of opts
  cdesc%upwind_1st = opts%upwind_1st
  cdesc%monotonic  = opts%monotonic
  cdesc%simple_2nd = opts%simple_2nd
end function reconstruction_CS_to_c

!> Converts a transport_adjust_CS to its bind(C) mirror.
function transport_adjust_CS_to_c(opts) result(cdesc)
  type(transport_adjust_CS), intent(in) :: opts !< Options to convert
  type(transport_adjust_CS_C) :: cdesc                !< bind(C) mirror of opts
  cdesc%tol_eta          = opts%tol_eta
  cdesc%tol_vel          = opts%tol_vel
  cdesc%CFL_limit_adjust = opts%CFL_limit_adjust
  cdesc%aggress_adjust   = opts%aggress_adjust
  cdesc%vol_CFL          = opts%vol_CFL
  cdesc%better_iter      = opts%better_iter
  cdesc%use_visc_rem_max = opts%use_visc_rem_max
  cdesc%marginal_faces   = opts%marginal_faces
end function transport_adjust_CS_to_c

!> Time steps the layer thicknesses, using a monotonically limit, directionally split PPM scheme,
!! based on Lin (1994).
subroutine continuity_PPM(u_a, v_a, hin_a, h_a, uh_a, vh_a, dt, bx0, stencil, x_first, &
                          mask2dT_a, dy_Cu_a, IareaT_a, IdxT_a, areaT_a, dxT_a, &
                          mask2dCu_a, dxCu_a, dx_Cv_a, IdyT_a, dyT_a, mask2dCv_a, dyCv_a, &
                          isd, ied, Angstrom_H, H_subroundoff, CS, OBC, pbv, &
                          uhbt_a, vhbt_a, visc_rem_u_a, visc_rem_v_a, u_cor_a, v_cor_a, BT_cont, &
                          du_cor_a, dv_cor_a)
  type(Box_t),              intent(in)    :: bx0  !< The core (unstencilled) iteration box
  integer,                  intent(in)    :: stencil !< The continuity solver stencil width
  logical,                  intent(in)    :: x_first !< If true, advect zonally before
                                                 !! meridionally.
  type(RealArray_t),       intent(in)    :: u_a   !< Zonal velocity [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in)    :: v_a   !< Meridional velocity [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in)    :: hin_a !< Initial layer thickness [H ~> m or kg m-2].
  type(RealArray_t),       intent(inout) :: h_a   !< Final layer thickness [H ~> m or kg m-2].
  type(RealArray_t),       intent(inout) :: uh_a  !< Zonal volume flux, u*h*dy
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(RealArray_t),       intent(inout) :: vh_a  !< Meridional volume flux, v*h*dx
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real,                    intent(in)    :: dt  !< Time increment [T ~> s].
  type(RealArray_t),       intent(in)    :: mask2dT_a !< Cell land/ocean mask [nondim].
  type(RealArray_t),       intent(in)    :: dy_Cu_a  !< The grid cell's unblocked lengths of
                                                 !! the u-faces of the h-cell [L ~> m].
  type(RealArray_t),       intent(in)    :: IareaT_a !< The grid cell's 1/areaT [L-2 ~> m-2].
  type(RealArray_t),       intent(in)    :: IdxT_a   !< The grid cell's 1/dxT [L-1 ~> m-1].
  type(RealArray_t),       intent(in)    :: areaT_a  !< The area of the h-cell [L2 ~> m2].
  type(RealArray_t),       intent(in)    :: dxT_a    !< The x-extent of the h-cell [L ~> m].
  type(RealArray_t),       intent(in)    :: mask2dCu_a !< 0 for land points, 1 for ocean points
                                                 !! at u-locations [nondim].
  type(RealArray_t),       intent(in)    :: dxCu_a   !< The grid cell's u-point x-extent [L ~> m].
  type(RealArray_t),       intent(in)    :: dx_Cv_a  !< The grid cell's unblocked lengths of
                                                 !! the v-faces of the h-cell [L ~> m].
  type(RealArray_t),       intent(in)    :: IdyT_a   !< The grid cell's 1/dyT [L-1 ~> m-1].
  type(RealArray_t),       intent(in)    :: dyT_a    !< The y-extent of the h-cell [L ~> m].
  type(RealArray_t),       intent(in)    :: mask2dCv_a !< 0 for land points, 1 for ocean points
                                                 !! at v-locations [nondim].
  type(RealArray_t),       intent(in)    :: dyCv_a   !< The grid cell's v-point y-extent [L ~> m].
  integer,                 intent(in)    :: isd  !< The start i-index of the data domain.
  integer,                 intent(in)    :: ied  !< The end i-index of the data domain.
  real,                    intent(in)    :: Angstrom_H !< A one-Angstrom thickness
                                                 !! [H ~> m or kg m-2].
  real,                    intent(in)    :: H_subroundoff !< A negligibly small thickness used
                                                 !! to avoid division by zero [H ~> m or kg m-2].
  type(continuity_PPM_CS), intent(in)    :: CS  !< Module's control structure.
  type(ocean_OBC_type),    pointer       :: OBC !< Open boundaries control structure.
  type(porous_barrier_type), intent(in)  :: pbv !< pointers to porous barrier fractional cell metrics
  type(RealArray_t), &
                           intent(in)    :: uhbt_a !< The summed volume flux through zonal faces
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(RealArray_t), &
                           intent(in)    :: vhbt_a !< The summed volume flux through meridional
                                                 !! faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(RealArray_t), &
                           intent(in)    :: visc_rem_u_a
                             !< The fraction of zonal momentum originally
                             !! in a layer that remains after a time-step of viscosity, and the
                             !! fraction of a time-step's worth of a barotropic acceleration that
                             !! a layer experiences after viscosity is applied [nondim].
                             !! Visc_rem_u is between 0 (at the bottom) and 1 (far above the bottom).
  type(RealArray_t), &
                           intent(in)    :: visc_rem_v_a
                             !< The fraction of meridional momentum originally
                             !! in a layer that remains after a time-step of viscosity, and the
                             !! fraction of a time-step's worth of a barotropic acceleration that
                             !! a layer experiences after viscosity is applied [nondim].
                             !! Visc_rem_v is between 0 (at the bottom) and 1 (far above the bottom).
  type(RealArray_t), &
                           intent(inout) :: u_cor_a
                             !< The zonal velocities that give uhbt as the depth-integrated transport [L T-1 ~> m s-1].
  type(RealArray_t), &
                           intent(inout) :: v_cor_a
                             !< The meridional velocities that give vhbt as the depth-integrated
                             !! transport [L T-1 ~> m s-1].
  type(BT_cont_type), optional, pointer  :: BT_cont !< A structure with elements that describe
                             !!  the effective open face areas as a function of barotropic flow.
  type(RealArray_t), &
                           intent(inout) :: du_cor_a !< The zonal velocity increments from u that
                                                 !! give uhbt as the depth-integrated
                                                 !! transports [L T-1 ~> m s-1].
  type(RealArray_t), &
                           intent(inout) :: dv_cor_a !< The meridional velocity increments from v
                                                 !! that give vhbt as the depth-integrated
                                                 !! transports [L T-1 ~> m s-1].

  ! Local variables
  real :: h_min  ! The minimum layer thickness [H ~> m or kg m-2].  h_min could be 0.
  type(Box_t) :: bxC                ! An iteration box
  real, dimension(:,:,:), contiguous, pointer :: u, v, hin, h, uh, vh
  type(RealArray_t) :: h_W_a, h_E_a, h_S_a, h_N_a
  type(RealArray_t) :: por_face_areaU_a, por_face_areaV_a
  type(RealArray_t) :: hin_a_none ! Never allocated -- the h_min branch uses h_min, not hin.

  call u_a%view(u)
  call v_a%view(v)
  call hin_a%view(hin)
  call h_a%view(h)
  call uh_a%view(uh)
  call vh_a%view(vh)

  h_min = Angstrom_H

  if (.not.CS%initialized) call MOM_error(FATAL, &
         "MOM_continuity_PPM: Module must be initialized before it is used.")

  if (visc_rem_u_a%associated() .neqv. visc_rem_v_a%associated()) call MOM_error(FATAL, &
      "MOM_continuity_PPM: Either both visc_rem_u and visc_rem_v or neither "// &
      "one must be present in call to continuity_PPM.")

  if (x_first) then
    !  First advect zonally, with loop bounds that accomodate the subsequent meridional advection.
    bxC = bx0%grow(dim=2, n=stencil)
    call h_W_a%alloc(lb=[h_a%lb(1),h_a%lb(2),h_a%lb(3)], ub=[h_a%ub(1),h_a%ub(2),h_a%ub(3)])
    call h_E_a%alloc(lb=[h_a%lb(1),h_a%lb(2),h_a%lb(3)], ub=[h_a%ub(1),h_a%ub(2),h_a%ub(3)])
    call zonal_edge_thickness(bxC, hin_a, h_W_a, h_E_a, mask2dT_a, Angstrom_H, &
                              CS%reconstruction_CS, OBC)
    call por_face_areaU_a%alloc(lb=LBOUND(pbv%por_face_areaU), ub=UBOUND(pbv%por_face_areaU), &
                                source=pbv%por_face_areaU)
    call zonal_mass_flux(bxC, u_a, hin_a, h_W_a, h_E_a, uh_a, dt, &
                         dy_Cu_a, IareaT_a, IdxT_a, areaT_a, dxT_a, mask2dCu_a, dxCu_a, &
                         H_subroundoff, &
                         CS%transport_adjust_CS, OBC, por_face_areaU_a, uhbt_a=uhbt_a, &
                         visc_rem_u_a=visc_rem_u_a, u_cor_a=u_cor_a, BT_cont=BT_cont, &
                         du_cor_a=du_cor_a)
    call h_W_a%free() ; call h_E_a%free() ; call por_face_areaU_a%free()
    call continuity_zonal_convergence(bxC, h_a, uh_a, dt, IareaT_a, hin_a=hin_a)

    ! update host h from continuity_zonal_convergence

    !  Now advect meridionally, using the updated thicknesses to determine the fluxes.
    bxC = bx0
    call h_S_a%alloc(lb=[h_a%lb(1),h_a%lb(2),h_a%lb(3)], ub=[h_a%ub(1),h_a%ub(2),h_a%ub(3)])
    call h_N_a%alloc(lb=[h_a%lb(1),h_a%lb(2),h_a%lb(3)], ub=[h_a%ub(1),h_a%ub(2),h_a%ub(3)])
    call meridional_edge_thickness(bxC, h_a, h_S_a, h_N_a, mask2dT_a, Angstrom_H, &
                                   CS%reconstruction_CS, OBC)
    call por_face_areaV_a%alloc(lb=LBOUND(pbv%por_face_areaV), ub=UBOUND(pbv%por_face_areaV), &
                                source=pbv%por_face_areaV)
    call meridional_mass_flux(bxC, v_a, h_a, h_S_a, h_N_a, vh_a, dt, &
                              dx_Cv_a, IareaT_a, IdyT_a, areaT_a, dyT_a, mask2dCv_a, dyCv_a, &
                              isd, ied, H_subroundoff, &
                              CS%transport_adjust_CS, OBC, por_face_areaV_a, vhbt_a=vhbt_a, &
                              visc_rem_v_a=visc_rem_v_a, v_cor_a=v_cor_a, BT_cont=BT_cont, &
                              dv_cor_a=dv_cor_a)
    call h_S_a%free() ; call h_N_a%free() ; call por_face_areaV_a%free()
    call continuity_meridional_convergence(bxC, h_a, vh_a, dt, IareaT_a, hin_a=hin_a_none, &
                                           hmin=h_min)

  else  ! .not. x_first
    !  First advect meridionally, with loop bounds that accomodate the subsequent zonal advection.
    bxC = bx0%grow(dim=1, n=stencil)
    call h_S_a%alloc(lb=[h_a%lb(1),h_a%lb(2),h_a%lb(3)], ub=[h_a%ub(1),h_a%ub(2),h_a%ub(3)])
    call h_N_a%alloc(lb=[h_a%lb(1),h_a%lb(2),h_a%lb(3)], ub=[h_a%ub(1),h_a%ub(2),h_a%ub(3)])
    call meridional_edge_thickness(bxC, hin_a, h_S_a, h_N_a, mask2dT_a, Angstrom_H, &
                                   CS%reconstruction_CS, OBC)
    call por_face_areaV_a%alloc(lb=LBOUND(pbv%por_face_areaV), ub=UBOUND(pbv%por_face_areaV), &
                                source=pbv%por_face_areaV)
    call meridional_mass_flux(bxC, v_a, hin_a, h_S_a, h_N_a, vh_a, dt, &
                              dx_Cv_a, IareaT_a, IdyT_a, areaT_a, dyT_a, mask2dCv_a, dyCv_a, &
                              isd, ied, H_subroundoff, &
                              CS%transport_adjust_CS, OBC, por_face_areaV_a, vhbt_a=vhbt_a, &
                              visc_rem_v_a=visc_rem_v_a, v_cor_a=v_cor_a, BT_cont=BT_cont, &
                              dv_cor_a=dv_cor_a)
    call h_S_a%free() ; call h_N_a%free() ; call por_face_areaV_a%free()
    call continuity_meridional_convergence(bxC, h_a, vh_a, dt, IareaT_a, hin_a=hin_a)

    !  Now advect zonally, using the updated thicknesses to determine the fluxes.
    bxC = bx0
    call h_W_a%alloc(lb=[h_a%lb(1),h_a%lb(2),h_a%lb(3)], ub=[h_a%ub(1),h_a%ub(2),h_a%ub(3)])
    call h_E_a%alloc(lb=[h_a%lb(1),h_a%lb(2),h_a%lb(3)], ub=[h_a%ub(1),h_a%ub(2),h_a%ub(3)])
    call zonal_edge_thickness(bxC, h_a, h_W_a, h_E_a, mask2dT_a, Angstrom_H, &
                              CS%reconstruction_CS, OBC)
    call por_face_areaU_a%alloc(lb=LBOUND(pbv%por_face_areaU), ub=UBOUND(pbv%por_face_areaU), &
                                source=pbv%por_face_areaU)
    call zonal_mass_flux(bxC, u_a, h_a, h_W_a, h_E_a, uh_a, dt, &
                         dy_Cu_a, IareaT_a, IdxT_a, areaT_a, dxT_a, mask2dCu_a, dxCu_a, &
                         H_subroundoff, CS%transport_adjust_CS, &
                         OBC, por_face_areaU_a, uhbt_a=uhbt_a, visc_rem_u_a=visc_rem_u_a, &
                         u_cor_a=u_cor_a, BT_cont=BT_cont, du_cor_a=du_cor_a)
    call h_W_a%free() ; call h_E_a%free() ; call por_face_areaU_a%free()
    call continuity_zonal_convergence(bxC, h_a, uh_a, dt, IareaT_a, hin_a=hin_a_none, hmin=h_min)
  endif

  ! Free the continuity solver iteration box
  call bxC%free()

end subroutine continuity_PPM

!> Finds the thickness fluxes from the continuity solver without actually updating the
!! layer thicknesses.  Because the fluxes in the two directions are calculated based on the
!! input thicknesses, which are not updated between the direcitons, the fluxes returned here
!! are not the same as those that would be returned by a call to continuity.
subroutine continuity_PPM_3d_fluxes(u_a, v_a, h_a, uh_a, vh_a, dt, bxC, &
                                    mask2dT_a, dy_Cu_a, IareaT_a, IdxT_a, areaT_a, dxT_a, &
                                    mask2dCu_a, dxCu_a, dx_Cv_a, IdyT_a, dyT_a, mask2dCv_a, &
                                    dyCv_a, isd, ied, Angstrom_H, H_subroundoff, CS, OBC, pbv)
  type(Box_t),              intent(in)    :: bxC  !< Iteration box for continuity solver
  type(RealArray_t),       intent(in)    :: u_a  !< Zonal velocity [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in)    :: v_a  !< Meridional velocity [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in)    :: h_a  !< Layer thickness [H ~> m or kg m-2].
  type(RealArray_t),       intent(inout) :: uh_a !< Thickness fluxes through zonal faces,
                                                !! u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(RealArray_t),       intent(inout) :: vh_a !< Thickness fluxes through meridional faces,
                                                !! v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1].
  real,                    intent(in)    :: dt  !< Time increment [T ~> s].
  type(RealArray_t),       intent(in)    :: mask2dT_a !< Cell land/ocean mask [nondim].
  type(RealArray_t),       intent(in)    :: dy_Cu_a  !< The grid cell's unblocked lengths of
                                                 !! the u-faces of the h-cell [L ~> m].
  type(RealArray_t),       intent(in)    :: IareaT_a !< The grid cell's 1/areaT [L-2 ~> m-2].
  type(RealArray_t),       intent(in)    :: IdxT_a   !< The grid cell's 1/dxT [L-1 ~> m-1].
  type(RealArray_t),       intent(in)    :: areaT_a  !< The area of the h-cell [L2 ~> m2].
  type(RealArray_t),       intent(in)    :: dxT_a    !< The x-extent of the h-cell [L ~> m].
  type(RealArray_t),       intent(in)    :: mask2dCu_a !< 0 for land points, 1 for ocean points
                                                 !! at u-locations [nondim].
  type(RealArray_t),       intent(in)    :: dxCu_a   !< The grid cell's u-point x-extent [L ~> m].
  type(RealArray_t),       intent(in)    :: dx_Cv_a  !< The grid cell's unblocked lengths of
                                                 !! the v-faces of the h-cell [L ~> m].
  type(RealArray_t),       intent(in)    :: IdyT_a   !< The grid cell's 1/dyT [L-1 ~> m-1].
  type(RealArray_t),       intent(in)    :: dyT_a    !< The y-extent of the h-cell [L ~> m].
  type(RealArray_t),       intent(in)    :: mask2dCv_a !< 0 for land points, 1 for ocean points
                                                 !! at v-locations [nondim].
  type(RealArray_t),       intent(in)    :: dyCv_a   !< The grid cell's v-point y-extent [L ~> m].
  integer,                 intent(in)    :: isd  !< The start i-index of the data domain.
  integer,                 intent(in)    :: ied  !< The end i-index of the data domain.
  real,                    intent(in)    :: Angstrom_H !< A one-Angstrom thickness
                                                 !! [H ~> m or kg m-2].
  real,                    intent(in)    :: H_subroundoff !< A negligibly small thickness used
                                                 !! to avoid division by zero [H ~> m or kg m-2].
  type(continuity_PPM_CS), intent(in)    :: CS  !< Control structure for mom_continuity.
  type(ocean_OBC_type),    pointer       :: OBC !< Open boundaries control structure.
  type(porous_barrier_type), intent(in)  :: pbv !< porous barrier fractional cell metrics

  ! Local variables
  real, dimension(:,:,:), contiguous, pointer :: u, v, h, uh, vh
  type(RealArray_t) :: h_W_a, h_E_a, h_S_a, h_N_a
  type(RealArray_t) :: por_face_areaU_a, por_face_areaV_a
  ! Never allocated -- this caller does not report the barotropic-consistency outputs.
  type(RealArray_t) :: uhbt_a, visc_rem_u_a, u_cor_a, du_cor_a
  type(RealArray_t) :: vhbt_a, visc_rem_v_a, v_cor_a, dv_cor_a

  call u_a%view(u)
  call v_a%view(v)
  call h_a%view(h)
  call uh_a%view(uh)
  call vh_a%view(vh)

  call h_W_a%alloc(lb=[h_a%lb(1),h_a%lb(2),h_a%lb(3)], ub=[h_a%ub(1),h_a%ub(2),h_a%ub(3)])
  call h_E_a%alloc(lb=[h_a%lb(1),h_a%lb(2),h_a%lb(3)], ub=[h_a%ub(1),h_a%ub(2),h_a%ub(3)])
  call zonal_edge_thickness(bxC, h_a, h_W_a, h_E_a, mask2dT_a, Angstrom_H, &
                            CS%reconstruction_CS, OBC)
  call por_face_areaU_a%alloc(lb=LBOUND(pbv%por_face_areaU), ub=UBOUND(pbv%por_face_areaU), &
                              source=pbv%por_face_areaU)
  call zonal_mass_flux(bxC, u_a, h_a, h_W_a, h_E_a, uh_a, dt, &
                       dy_Cu_a, IareaT_a, IdxT_a, areaT_a, dxT_a, mask2dCu_a, dxCu_a, &
                       H_subroundoff, CS%transport_adjust_CS, &
                       OBC, por_face_areaU_a, uhbt_a=uhbt_a, visc_rem_u_a=visc_rem_u_a, &
                       u_cor_a=u_cor_a, du_cor_a=du_cor_a)
  call h_W_a%free() ; call h_E_a%free() ; call por_face_areaU_a%free()

  call h_S_a%alloc(lb=[h_a%lb(1),h_a%lb(2),h_a%lb(3)], ub=[h_a%ub(1),h_a%ub(2),h_a%ub(3)])
  call h_N_a%alloc(lb=[h_a%lb(1),h_a%lb(2),h_a%lb(3)], ub=[h_a%ub(1),h_a%ub(2),h_a%ub(3)])
  call meridional_edge_thickness(bxC, h_a, h_S_a, h_N_a, mask2dT_a, Angstrom_H, &
                                 CS%reconstruction_CS, OBC)
  call por_face_areaV_a%alloc(lb=LBOUND(pbv%por_face_areaV), ub=UBOUND(pbv%por_face_areaV), &
                              source=pbv%por_face_areaV)
  call meridional_mass_flux(bxC, v_a, h_a, h_S_a, h_N_a, vh_a, dt, &
                            dx_Cv_a, IareaT_a, IdyT_a, areaT_a, dyT_a, mask2dCv_a, dyCv_a, &
                            isd, ied, H_subroundoff, &
                            CS%transport_adjust_CS, OBC, por_face_areaV_a, vhbt_a=vhbt_a, &
                            visc_rem_v_a=visc_rem_v_a, v_cor_a=v_cor_a, dv_cor_a=dv_cor_a)
  call h_S_a%free() ; call h_N_a%free() ; call por_face_areaV_a%free()

end subroutine continuity_PPM_3d_fluxes

!> Find the vertical sum of the thickness fluxes from the continuity solver without actually
!! updating the layer thicknesses.  Because the fluxes in the two directions are calculated
!! based on the input thicknesses, which are not updated between the directions, the fluxes
!! returned here are not the same as those that would be returned by a call to continuity.
subroutine continuity_PPM_2d_fluxes(u_a, v_a, h_a, uhbt_a, vhbt_a, dt, bxC, &
                                    mask2dT_a, dy_Cu_a, IareaT_a, IdxT_a, dx_Cv_a, IdyT_a, &
                                    Angstrom_H, CS, OBC, pbv)
  type(Box_t),              intent(in)    :: bxC  !< Iteration box for continuity solver
  type(RealArray_t),       intent(in)    :: u_a  !< Zonal velocity [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in)    :: v_a  !< Meridional velocity [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in)    :: h_a  !< Layer thickness [H ~> m or kg m-2].
  type(RealArray_t),       intent(inout) :: uhbt_a !< Vertically summed thickness flux through
                                                !! zonal faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(RealArray_t),       intent(inout) :: vhbt_a !< Vertically summed thickness flux through
                                                !! meridional faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  real,                    intent(in)    :: dt  !< Time increment [T ~> s].
  type(RealArray_t),       intent(in)    :: mask2dT_a !< Cell land/ocean mask [nondim].
  type(RealArray_t),       intent(in)    :: dy_Cu_a  !< The grid cell's unblocked lengths of
                                                !! the u-faces of the h-cell [L ~> m].
  type(RealArray_t),       intent(in)    :: IareaT_a !< The grid cell's 1/areaT [L-2 ~> m-2].
  type(RealArray_t),       intent(in)    :: IdxT_a   !< The grid cell's 1/dxT [L-1 ~> m-1].
  type(RealArray_t),       intent(in)    :: dx_Cv_a  !< The grid cell's unblocked lengths of
                                                !! the v-faces of the h-cell [L ~> m].
  type(RealArray_t),       intent(in)    :: IdyT_a   !< The grid cell's 1/dyT [L-1 ~> m-1].
  real,                    intent(in)    :: Angstrom_H !< A one-Angstrom thickness
                                                !! [H ~> m or kg m-2].
  type(continuity_PPM_CS), intent(in)    :: CS  !< Control structure for mom_continuity.
  type(ocean_OBC_type),    pointer       :: OBC !< Open boundaries control structure.
  type(porous_barrier_type), intent(in)  :: pbv !< porous barrier fractional cell metrics

  ! Local variables
  real, dimension(:,:,:), contiguous, pointer :: u, v, h
  type(RealArray_t) :: h_W_a, h_E_a, h_S_a, h_N_a
  type(RealArray_t) :: por_face_areaU_a, por_face_areaV_a

  call u_a%view(u)
  call v_a%view(v)
  call h_a%view(h)

  call h_W_a%alloc(lb=[h_a%lb(1),h_a%lb(2),h_a%lb(3)], ub=[h_a%ub(1),h_a%ub(2),h_a%ub(3)])
  call h_E_a%alloc(lb=[h_a%lb(1),h_a%lb(2),h_a%lb(3)], ub=[h_a%ub(1),h_a%ub(2),h_a%ub(3)])
  call zonal_edge_thickness(bxC, h_a, h_W_a, h_E_a, mask2dT_a, Angstrom_H, &
                            CS%reconstruction_CS, OBC)
  call por_face_areaU_a%alloc(lb=LBOUND(pbv%por_face_areaU), ub=UBOUND(pbv%por_face_areaU), &
                              source=pbv%por_face_areaU)
  call zonal_BT_mass_flux(bxC, u_a, h_a, h_W_a, h_E_a, uhbt_a, dt, &
                          dy_Cu_a, IareaT_a, IdxT_a, CS%transport_adjust_CS, &
                          OBC, por_face_areaU_a)
  call h_W_a%free() ; call h_E_a%free() ; call por_face_areaU_a%free()

  call h_S_a%alloc(lb=[h_a%lb(1),h_a%lb(2),h_a%lb(3)], ub=[h_a%ub(1),h_a%ub(2),h_a%ub(3)])
  call h_N_a%alloc(lb=[h_a%lb(1),h_a%lb(2),h_a%lb(3)], ub=[h_a%ub(1),h_a%ub(2),h_a%ub(3)])
  call meridional_edge_thickness(bxC, h_a, h_S_a, h_N_a, mask2dT_a, Angstrom_H, &
                                 CS%reconstruction_CS, OBC)
  call por_face_areaV_a%alloc(lb=LBOUND(pbv%por_face_areaV), ub=UBOUND(pbv%por_face_areaV), &
                              source=pbv%por_face_areaV)
  call meridional_BT_mass_flux(bxC, v_a, h_a, h_S_a, h_N_a, vhbt_a, dt, &
                               dx_Cv_a, IareaT_a, IdyT_a, CS%transport_adjust_CS, &
                               OBC, por_face_areaV_a)
  call h_S_a%free() ; call h_N_a%free() ; call por_face_areaV_a%free()

end subroutine continuity_PPM_2d_fluxes

!> Correct the velocities to give the specified depth-integrated transports by applying a
!! barotropic acceleration (subject to viscous drag) to the velocities.
subroutine continuity_PPM_adjust_vel(u_a, v_a, h_a, dt, bxC, &
                                     mask2dT_a, dy_Cu_a, IareaT_a, IdxT_a, areaT_a, dxT_a, &
                                     mask2dCu_a, dxCu_a, dx_Cv_a, IdyT_a, dyT_a, mask2dCv_a, &
                                     dyCv_a, isd, ied, Angstrom_H, H_subroundoff, &
                                     CS, OBC, pbv, uhbt_a, vhbt_a, &
                                     visc_rem_u_a, visc_rem_v_a)
  type(Box_t),              intent(in)    :: bxC  !< Iteration box for continuity solver
  type(RealArray_t),       intent(inout) :: u_a  !< Zonal velocity, which will be adjusted to
                                                !! give uhbt as the depth-integrated
                                                !! transport [L T-1 ~> m s-1].
  type(RealArray_t),       intent(inout) :: v_a  !< Meridional velocity, which will be adjusted
                                                !! to give vhbt as the depth-integrated
                                                !! transport [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in)    :: h_a  !< Layer thickness [H ~> m or kg m-2].
  real,                    intent(in)    :: dt  !< Time increment [T ~> s].
  type(RealArray_t),       intent(in)    :: mask2dT_a !< Cell land/ocean mask [nondim].
  type(RealArray_t),       intent(in)    :: dy_Cu_a  !< The grid cell's unblocked lengths of
                                                !! the u-faces of the h-cell [L ~> m].
  type(RealArray_t),       intent(in)    :: IareaT_a !< The grid cell's 1/areaT [L-2 ~> m-2].
  type(RealArray_t),       intent(in)    :: IdxT_a   !< The grid cell's 1/dxT [L-1 ~> m-1].
  type(RealArray_t),       intent(in)    :: areaT_a  !< The area of the h-cell [L2 ~> m2].
  type(RealArray_t),       intent(in)    :: dxT_a    !< The x-extent of the h-cell [L ~> m].
  type(RealArray_t),       intent(in)    :: mask2dCu_a !< 0 for land points, 1 for ocean points
                                                !! at u-locations [nondim].
  type(RealArray_t),       intent(in)    :: dxCu_a   !< The grid cell's u-point x-extent [L ~> m].
  type(RealArray_t),       intent(in)    :: dx_Cv_a  !< The grid cell's unblocked lengths of
                                                !! the v-faces of the h-cell [L ~> m].
  type(RealArray_t),       intent(in)    :: IdyT_a   !< The grid cell's 1/dyT [L-1 ~> m-1].
  type(RealArray_t),       intent(in)    :: dyT_a    !< The y-extent of the h-cell [L ~> m].
  type(RealArray_t),       intent(in)    :: mask2dCv_a !< 0 for land points, 1 for ocean points
                                                !! at v-locations [nondim].
  type(RealArray_t),       intent(in)    :: dyCv_a   !< The grid cell's v-point y-extent [L ~> m].
  integer,                 intent(in)    :: isd  !< The start i-index of the data domain.
  integer,                 intent(in)    :: ied  !< The end i-index of the data domain.
  real,                    intent(in)    :: Angstrom_H !< A one-Angstrom thickness
                                                !! [H ~> m or kg m-2].
  real,                    intent(in)    :: H_subroundoff !< A negligibly small thickness used
                                                !! to avoid division by zero [H ~> m or kg m-2].
  type(continuity_PPM_CS), intent(in)    :: CS  !< Control structure for mom_continuity.
  type(ocean_OBC_type),    pointer       :: OBC !< Open boundaries control structure.
  type(porous_barrier_type), intent(in)  :: pbv !< porous barrier fractional cell metrics
  type(RealArray_t),       intent(in)    :: uhbt_a !< The vertically summed thickness flux through
                                                !! zonal faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(RealArray_t),       intent(in)    :: vhbt_a !< The vertically summed thickness flux through
                                                !! meridional faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(RealArray_t), &
                           intent(in)    :: visc_rem_u_a !< Both the fraction of the zonal momentum
                                                !! that remains after a time-step of viscosity, and
                                                !! the fraction of a time-step's worth of a barotropic
                                                !! acceleration that a layer experiences after viscosity
                                                !! is applied [nondim].  This goes between 0 (at the
                                                !! bottom) and 1 (far above the bottom).  When this
                                                !! column is under an ice shelf, this also goes to 0
                                                !! at the top due to the no-slip boundary condition there.
  type(RealArray_t), &
                           intent(in)    :: visc_rem_v_a !< Both the fraction of the meridional
                                                !! momentum
                                                !! that remains after a time-step of viscosity, and
                                                !! the fraction of a time-step's worth of a barotropic
                                                !! acceleration that a layer experiences after viscosity
                                                !! is applied [nondim].  This goes between 0 (at the
                                                !! bottom) and 1 (far above the bottom).  When this
                                                !! column is under an ice shelf, this also goes to 0
                                                !! at the top due to the no-slip boundary condition there.

  ! Local variables
  type(RealArray_t) :: u_cor_a, v_cor_a
  type(RealArray_t) :: du_cor_a, dv_cor_a ! Never allocated -- this caller does not report these.
  real, dimension(:,:,:), contiguous, pointer :: u, v, h
  type(RealArray_t) :: h_W_a, h_E_a, h_S_a, h_N_a
  type(RealArray_t) :: u_in_a, v_in_a, uh_a, vh_a, por_face_areaU_a, por_face_areaV_a

  call u_a%view(u)
  call v_a%view(v)
  call h_a%view(h)

  ! It might not be necessary to separate the input velocity array from the adjusted velocities,
  ! but it seems safer to do so, even if it might be less efficient.

  call h_W_a%alloc(lb=[h_a%lb(1),h_a%lb(2),h_a%lb(3)], ub=[h_a%ub(1),h_a%ub(2),h_a%ub(3)])
  call h_E_a%alloc(lb=[h_a%lb(1),h_a%lb(2),h_a%lb(3)], ub=[h_a%ub(1),h_a%ub(2),h_a%ub(3)])
  call zonal_edge_thickness(bxC, h_a, h_W_a, h_E_a, mask2dT_a, Angstrom_H, &
                            CS%reconstruction_CS, OBC)
  call u_cor_a%alloc(lb=LBOUND(u), ub=UBOUND(u), source=u)
  call u_in_a%alloc(lb=LBOUND(u), ub=UBOUND(u), source=u)
  call uh_a%alloc(lb=[u_a%lb(1),u_a%lb(2),u_a%lb(3)], ub=[u_a%ub(1),u_a%ub(2),u_a%ub(3)])
  call por_face_areaU_a%alloc(lb=LBOUND(pbv%por_face_areaU), ub=UBOUND(pbv%por_face_areaU), &
                              source=pbv%por_face_areaU)
  call zonal_mass_flux(bxC, u_in_a, h_a, h_W_a, h_E_a, uh_a, dt, &
                       dy_Cu_a, IareaT_a, IdxT_a, areaT_a, dxT_a, mask2dCu_a, dxCu_a, &
                       H_subroundoff, CS%transport_adjust_CS, &
                       OBC, por_face_areaU_a, uhbt_a=uhbt_a, &
                       visc_rem_u_a=visc_rem_u_a, u_cor_a=u_cor_a, du_cor_a=du_cor_a)
  call u_cor_a%copy2F(u) ; call u_cor_a%free()
  call h_W_a%free() ; call h_E_a%free() ; call u_in_a%free() ; call uh_a%free()
  call por_face_areaU_a%free()

  call h_S_a%alloc(lb=[h_a%lb(1),h_a%lb(2),h_a%lb(3)], ub=[h_a%ub(1),h_a%ub(2),h_a%ub(3)])
  call h_N_a%alloc(lb=[h_a%lb(1),h_a%lb(2),h_a%lb(3)], ub=[h_a%ub(1),h_a%ub(2),h_a%ub(3)])
  call meridional_edge_thickness(bxC, h_a, h_S_a, h_N_a, mask2dT_a, Angstrom_H, &
                                 CS%reconstruction_CS, OBC)
  call v_cor_a%alloc(lb=LBOUND(v), ub=UBOUND(v), source=v)
  call v_in_a%alloc(lb=LBOUND(v), ub=UBOUND(v), source=v)
  call vh_a%alloc(lb=[v_a%lb(1),v_a%lb(2),v_a%lb(3)], ub=[v_a%ub(1),v_a%ub(2),v_a%ub(3)])
  call por_face_areaV_a%alloc(lb=LBOUND(pbv%por_face_areaV), ub=UBOUND(pbv%por_face_areaV), &
                              source=pbv%por_face_areaV)
  call meridional_mass_flux(bxC, v_in_a, h_a, h_S_a, h_N_a, vh_a, dt, &
                            dx_Cv_a, IareaT_a, IdyT_a, areaT_a, dyT_a, mask2dCv_a, dyCv_a, &
                            isd, ied, H_subroundoff, &
                            CS%transport_adjust_CS, OBC, por_face_areaV_a, vhbt_a=vhbt_a, &
                            visc_rem_v_a=visc_rem_v_a, v_cor_a=v_cor_a, dv_cor_a=dv_cor_a)
  call v_cor_a%copy2F(v) ; call v_cor_a%free()
  call h_S_a%free() ; call h_N_a%free() ; call v_in_a%free() ; call vh_a%free()
  call por_face_areaV_a%free()

end subroutine continuity_PPM_adjust_vel


!> Updates the thicknesses due to zonal thickness fluxes.
!> Original Fortran implementation of continuity_zonal_convergence (renamed). Takes containers.
subroutine continuity_zonal_convergence_fortran(bxC, h_a, uh_a, dt, IareaT_a, hin_a, hmin)
  type(box_t), intent(in) :: bxC                 !< Iteration box for continuity solver
  type(RealArray_t),          intent(inout) :: h_a  !< Final layer thickness [H ~> m or kg m-2]
  type(RealArray_t),          intent(in)    :: uh_a !< Zonal thickness flux, u*h*dy
                                                     !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real,                        intent(in)    :: dt   !< Time increment [T ~> s]
  type(RealArray_t),          intent(in)    :: IareaT_a !< The grid cell's 1/areaT [L-2 ~> m-2]
  type(RealArray_t), &
                               intent(in)    :: hin_a !< Initial layer thickness [H ~> m or kg m-2].
                                                     !! If hin is absent, h is also the initial thickness.
  real,              optional, intent(in)    :: hmin !< The minimum layer thickness [H ~> m or kg m-2]

  real :: h_min  ! The minimum layer thickness [H ~> m or kg m-2].  h_min could be 0.
  integer :: i, j, k, ish, ieh, jsh, jeh, nz
  real, dimension(:,:,:), contiguous, pointer :: hin, h, uh
  real, dimension(:,:), contiguous, pointer :: IareaT

  call h_a%view(h)
  call uh_a%view(uh)
  call IareaT_a%view(IareaT)
  nullify(hin)
  if (hin_a%associated()) call hin_a%view(hin)

  h_min = 0.0 ; if (present(hmin)) h_min = hmin

  if (hin_a%associated()) then
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

end subroutine continuity_zonal_convergence_fortran

!> Shim for continuity_zonal_convergence -- dispatches via CONTINUITY_ZONAL_CONVERGENCE_MODE env var.
subroutine continuity_zonal_convergence(bxC, h_a, uh_a, dt, IareaT_a, hin_a, hmin)
  type(box_t), intent(in) :: bxC                 !< Iteration box for continuity solver
  type(RealArray_t),          intent(inout) :: h_a  !< Final layer thickness [H ~> m or kg m-2]
  type(RealArray_t),          intent(in)    :: uh_a !< Zonal thickness flux, u*h*dy
                                                     !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real,                        intent(in)    :: dt   !< Time increment [T ~> s]
  type(RealArray_t),          intent(in)    :: IareaT_a !< The grid cell's 1/areaT [L-2 ~> m-2]
  type(RealArray_t), &
                               intent(in)    :: hin_a !< Initial layer thickness [H ~> m or kg m-2].
                                                     !! If hin is absent, h is also the initial thickness.
  real,              optional, intent(in)    :: hmin !< The minimum layer thickness [H ~> m or kg m-2]

  integer :: mode, rc
  real :: h_min
  type(RealArray_C) :: h_c, uh_c, IareaT_c, hin_c
  type(Box_C)        :: bxC_c
  type(io_recorder)  :: rec
  logical            :: capture
  character(len=80)  :: kernel
  character(len=100) :: dir
  character(len=256) :: binFile, metaFile

  kernel = "continuity_zonal_convergence"
  h_min  = 0.0 ; if (present(hmin)) h_min = hmin

  call cpu_clock_begin(id_clock_update)

  mode = getenv_mode("CONTINUITY_ZONAL_CONVERGENCE_MODE", default=TIMH_runFORTRAN)

  select case (mode)

    case (TIMH_capture)
      capture = (.not. already_recorded(trim(kernel))) .and. is_root_pe()
      if (capture) then
        dir = "capture"
        rc = mkdir_posix(trim(dir) // c_null_char, int(o'755', c_int))
        binFile  = trim(dir) // "/" // trim(kernel) // ".bin"
        metaFile = trim(dir) // "/" // trim(kernel) // ".meta"
        call rec%open_write(binFile, metaFile)
        call rec%add("_bxC",      bxC)
        call rec%add("_h_before", h_a)
        call rec%add("_uh",       uh_a)
        call rec%add("_dt",       dt)
        call rec%add("_IareaT",   IareaT_a)
        call rec%add("_hin",      hin_a)
        call rec%add("_h_min",    h_min)
      endif

      call continuity_zonal_convergence_fortran(bxC, h_a, uh_a, dt, IareaT_a, hin_a, hmin)

      if (capture) then
        call rec%add("_h_after", h_a)
        call rec%close()
        call mark_recorded(trim(kernel))
      endif

#ifdef _TIM
    case (TIMH_runAMREX)
      bxC_c    = bxC%to_c()
      h_c      = h_a%to_c()
      uh_c     = uh_a%to_c()
      IareaT_c = IareaT_a%to_c()
      hin_c    = hin_a%to_c()
      call turbotmp_continuity_zonal_convergence_bridge(bxC_c, h_c, uh_c, dt, IareaT_c, hin_c, &
                                                        h_min)
#endif

    case default
      call continuity_zonal_convergence_fortran(bxC, h_a, uh_a, dt, IareaT_a, hin_a, hmin)

  end select

  call cpu_clock_end(id_clock_update)

end subroutine continuity_zonal_convergence

!> Updates the thicknesses due to meridional thickness fluxes.
!> Original Fortran implementation of continuity_meridional_convergence (renamed). Takes containers.
subroutine continuity_meridional_convergence_fortran(bxC, h_a, vh_a, dt, IareaT_a, hin_a, hmin)
  type(box_t), intent(in) :: bxC                 !< Iteration box for continuity solver
  type(RealArray_t),          intent(inout) :: h_a  !< Final layer thickness [H ~> m or kg m-2]
  type(RealArray_t),          intent(in)    :: vh_a !< Meridional thickness flux, v*h*dx
                                                     !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real,                        intent(in)    :: dt   !< Time increment [T ~> s]
  type(RealArray_t),          intent(in)    :: IareaT_a !< The grid cell's 1/areaT [L-2 ~> m-2]
  type(RealArray_t), &
                               intent(in)    :: hin_a !< Initial layer thickness [H ~> m or kg m-2].
                                                     !! If hin is absent, h is also the initial thickness.
  real,              optional, intent(in)    :: hmin !< The minimum layer thickness [H ~> m or kg m-2]

  real :: h_min  ! The minimum layer thickness [H ~> m or kg m-2].  h_min could be 0.
  integer :: i, j, k, ish, ieh, jsh, jeh, nz
  real, dimension(:,:,:), contiguous, pointer :: hin, h, vh
  real, dimension(:,:), contiguous, pointer :: IareaT

  call h_a%view(h)
  call vh_a%view(vh)
  call IareaT_a%view(IareaT)
  nullify(hin)
  if (hin_a%associated()) call hin_a%view(hin)

  h_min = 0.0 ; if (present(hmin)) h_min = hmin

  if (hin_a%associated()) then
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

end subroutine continuity_meridional_convergence_fortran

!> Shim for continuity_meridional_convergence -- dispatches via CONTINUITY_MERIDIONAL_CONVERGENCE_MODE
!! env var.
subroutine continuity_meridional_convergence(bxC, h_a, vh_a, dt, IareaT_a, hin_a, hmin)
  type(box_t), intent(in) :: bxC                 !< Iteration box for continuity solver
  type(RealArray_t),          intent(inout) :: h_a  !< Final layer thickness [H ~> m or kg m-2]
  type(RealArray_t),          intent(in)    :: vh_a !< Meridional thickness flux, v*h*dx
                                                     !! [H L2 T-1 ~> m3 s-1 or kg s-1]
  real,                        intent(in)    :: dt   !< Time increment [T ~> s]
  type(RealArray_t),          intent(in)    :: IareaT_a !< The grid cell's 1/areaT [L-2 ~> m-2]
  type(RealArray_t), &
                               intent(in)    :: hin_a !< Initial layer thickness [H ~> m or kg m-2].
                                                     !! If hin is absent, h is also the initial thickness.
  real,              optional, intent(in)    :: hmin !< The minimum layer thickness [H ~> m or kg m-2]

  integer :: mode, rc
  real :: h_min
  type(RealArray_C) :: h_c, vh_c, IareaT_c, hin_c
  type(Box_C)        :: bxC_c
  type(io_recorder)  :: rec
  logical            :: capture
  character(len=80)  :: kernel
  character(len=100) :: dir
  character(len=256) :: binFile, metaFile

  kernel = "continuity_meridional_convergence"
  h_min  = 0.0 ; if (present(hmin)) h_min = hmin

  call cpu_clock_begin(id_clock_update)

  mode = getenv_mode("CONTINUITY_MERIDIONAL_CONVERGENCE_MODE", default=TIMH_runFORTRAN)

  select case (mode)

    case (TIMH_capture)
      capture = (.not. already_recorded(trim(kernel))) .and. is_root_pe()
      if (capture) then
        dir = "capture"
        rc = mkdir_posix(trim(dir) // c_null_char, int(o'755', c_int))
        binFile  = trim(dir) // "/" // trim(kernel) // ".bin"
        metaFile = trim(dir) // "/" // trim(kernel) // ".meta"
        call rec%open_write(binFile, metaFile)
        call rec%add("_bxC",      bxC)
        call rec%add("_h_before", h_a)
        call rec%add("_vh",       vh_a)
        call rec%add("_dt",       dt)
        call rec%add("_IareaT",   IareaT_a)
        call rec%add("_hin",      hin_a)
        call rec%add("_h_min",    h_min)
      endif

      call continuity_meridional_convergence_fortran(bxC, h_a, vh_a, dt, IareaT_a, hin_a, hmin)

      if (capture) then
        call rec%add("_h_after", h_a)
        call rec%close()
        call mark_recorded(trim(kernel))
      endif

#ifdef _TIM
    case (TIMH_runAMREX)
      bxC_c    = bxC%to_c()
      h_c      = h_a%to_c()
      vh_c     = vh_a%to_c()
      IareaT_c = IareaT_a%to_c()
      hin_c    = hin_a%to_c()
      call turbotmp_continuity_meridional_convergence_bridge(bxC_c, h_c, vh_c, dt, &
                                                              IareaT_c, hin_c, h_min)
#endif

    case default
      call continuity_meridional_convergence_fortran(bxC, h_a, vh_a, dt, IareaT_a, hin_a, hmin)

  end select

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
subroutine zonal_edge_thickness(bxC, h_in_a, h_W_a, h_E_a, mask2dT_a, Angstrom_H, CS, OBC)
  type(box_t), intent(in) :: bxC                 !< Iteration box for continuity solver
  type(RealArray_t),       intent(in)    :: h_in_a !< Tracer cell layer thickness
                                                 !! [H ~> m or kg m-2].
  type(RealArray_t),       intent(inout) :: h_W_a  !< Western edge layer thickness
                                                 !! [H ~> m or kg m-2].
  type(RealArray_t),       intent(inout) :: h_E_a  !< Eastern edge layer thickness
                                                 !! [H ~> m or kg m-2].
  type(RealArray_t),       intent(in)    :: mask2dT_a !< Cell land/ocean mask [nondim].
  real,                    intent(in)    :: Angstrom_H !< A one-Angstrom thickness
                                                 !! [H ~> m or kg m-2].
  type(reconstruction_CS), intent(in) :: CS !< Options controlling the
                                                 !! edge-value reconstruction scheme.
  type(ocean_OBC_type),    pointer       :: OBC  !< Open boundaries control structure.

  integer  :: mode, rc
  real     :: h_min
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
  h_min  = 2.0 * Angstrom_H

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
        call rec%add("_upwind_1st", CS%upwind_1st)
        call rec%add("_monotonic",  CS%monotonic)
        call rec%add("_simple_2nd", CS%simple_2nd)
      endif

      call zonal_edge_thickness_fortran(bxC, h_in_a, h_W_a, h_E_a, mask2dT_a, &
                                        h_min, CS%upwind_1st, &
                                        CS%monotonic, &
                                        CS%simple_2nd, OBC)

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
      upwind_1st_c = CS%upwind_1st
      monotonic_c  = CS%monotonic
      simple_2nd_c = CS%simple_2nd
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
                                        h_min, CS%upwind_1st, &
                                        CS%monotonic, &
                                        CS%simple_2nd, OBC)

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
subroutine meridional_edge_thickness(bxC, h_in_a, h_S_a, h_N_a, mask2dT_a, Angstrom_H, CS, OBC)
  type(box_t), intent(in) :: bxC                 !< Iteration box for continuity solver
  type(RealArray_t),       intent(in)    :: h_in_a !< Tracer cell layer thickness
                                                 !! [H ~> m or kg m-2].
  type(RealArray_t),       intent(inout) :: h_S_a  !< Southern edge layer thickness
                                                 !! [H ~> m or kg m-2].
  type(RealArray_t),       intent(inout) :: h_N_a  !< Northern edge layer thickness
                                                 !! [H ~> m or kg m-2].
  type(RealArray_t),       intent(in)    :: mask2dT_a !< Cell land/ocean mask [nondim].
  real,                    intent(in)    :: Angstrom_H !< A one-Angstrom thickness
                                                 !! [H ~> m or kg m-2].
  type(reconstruction_CS), intent(in) :: CS !< Options controlling the
                                                 !! edge-value reconstruction scheme.
  type(ocean_OBC_type),    pointer       :: OBC  !< Open boundaries control structure.

  integer  :: mode, rc
  real     :: h_min
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
  h_min  = 2.0 * Angstrom_H

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
        call rec%add("_upwind_1st", CS%upwind_1st)
        call rec%add("_monotonic",  CS%monotonic)
        call rec%add("_simple_2nd", CS%simple_2nd)
      endif

      call meridional_edge_thickness_fortran(bxC, h_in_a, h_S_a, h_N_a, mask2dT_a, &
                                             h_min, CS%upwind_1st, &
                                             CS%monotonic, &
                                             CS%simple_2nd, OBC)

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
      upwind_1st_c = CS%upwind_1st
      monotonic_c  = CS%monotonic
      simple_2nd_c = CS%simple_2nd
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
                                             h_min, CS%upwind_1st, &
                                             CS%monotonic, &
                                             CS%simple_2nd, OBC)

  end select

  call cpu_clock_end(id_clock_reconstruct)

end subroutine meridional_edge_thickness


!> Calculates the mass or volume fluxes through the zonal faces, and other related quantities.
subroutine zonal_mass_flux(bxC, u_a, h_in_a, h_W_a, h_E_a, uh_a, dt, &
                           dy_Cu_a, IareaT_a, IdxT_a, areaT_a, dxT_a, mask2dCu_a, dxCu_a, &
                           H_subroundoff, CS, OBC, &
                           por_face_areaU_a, uhbt_a, visc_rem_u_a, u_cor_a, BT_cont, du_cor_a)
  type(Box_t),             intent(in)    :: bxC  !< Iteration box for continuity solver
  type(RealArray_t),       intent(in)    :: u_a    !< Zonal velocity [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in)    :: h_in_a !< Layer thickness used to calculate
                                                 !! fluxes [H ~> m or kg m-2].
  type(RealArray_t),       intent(in)    :: h_W_a !< Western edge thicknesses [H ~> m or kg m-2].
  type(RealArray_t),       intent(in)    :: h_E_a !< Eastern edge thicknesses [H ~> m or kg m-2].
  type(RealArray_t),       intent(inout) :: uh_a   !< Volume flux through zonal faces = u*h*dy
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real,                    intent(in)    :: dt   !< Time increment [T ~> s].
  type(RealArray_t),       intent(in)    :: dy_Cu_a  !< The grid cell's unblocked lengths of
                                                 !! the u-faces of the h-cell [L ~> m].
  type(RealArray_t),       intent(in)    :: IareaT_a !< The grid cell's 1/areaT [L-2 ~> m-2].
  type(RealArray_t),       intent(in)    :: IdxT_a   !< The grid cell's 1/dxT [L-1 ~> m-1].
  type(RealArray_t),       intent(in)    :: areaT_a  !< The area of the h-cell [L2 ~> m2].
  type(RealArray_t),       intent(in)    :: dxT_a    !< The x-extent of the h-cell [L ~> m].
  type(RealArray_t),       intent(in)    :: mask2dCu_a !< 0 for land points, 1 for ocean points at
                                                 !! u-locations [nondim].
  type(RealArray_t),       intent(in)    :: dxCu_a   !< The grid cell's u-point x-extent [L ~> m].
  real,                    intent(in)    :: H_subroundoff !< A negligibly small thickness used to
                                                 !! avoid division by zero [H ~> m or kg m-2].
  type(transport_adjust_CS), intent(in) :: CS !< Options controlling the
                                                 !! transport adjustment and barotropic-consistency
                                                 !! iteration.
  type(ocean_OBC_type),    pointer       :: OBC  !< Open boundaries control structure.
  type(RealArray_t),       intent(in)    :: por_face_areaU_a !< fractional open area of
                                                 !! U-faces [nondim]
  type(RealArray_t), &
                           intent(in)    :: uhbt_a !< The summed volume flux through zonal faces
                                                 !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(RealArray_t), &
                           intent(in)    :: visc_rem_u_a
                     !< The fraction of zonal momentum originally in a layer that remains after a
                     !! time-step of viscosity, and the fraction of a time-step's worth of a barotropic
                     !! acceleration that a layer experiences after viscosity is applied [nondim].
                     !! Visc_rem_u is between 0 (at the bottom) and 1 (far above the bottom).
  type(RealArray_t), &
                           intent(inout) :: u_cor_a
                     !< The zonal velocities (u with a barotropic correction)
                     !! that give uhbt as the depth-integrated transport [L T-1 ~> m s-1]
  type(BT_cont_type), optional, pointer  :: BT_cont !< A structure with elements that describe the
                     !! effective open face areas as a function of barotropic flow.
  type(RealArray_t), &
                           intent(inout) :: du_cor_a !< The zonal velocity increments from u that give uhbt
                                                 !! as the depth-integrated transports [L T-1 ~> m s-1].

  ! Local variables
  real, dimension(u_a%lb(1):u_a%ub(1), u_a%lb(2):u_a%ub(2), u_a%lb(3):u_a%ub(3)) :: &
    duhdu ! Partial derivative of uh with u [H L ~> m2 or kg m-1].
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
  logical :: local_specified_BC, local_Flather_OBC, local_open_BC, any_simple_OBC  ! OBC-related logicals
  logical, dimension(u_a%lb(1):u_a%ub(1), u_a%lb(2):u_a%ub(2)) :: &
    simple_OBC_pt  ! Indicates points in a row with specified transport OBCs
  real, dimension(:,:), contiguous, pointer :: FA_u_W0, FA_u_E0, FA_u_WW, FA_u_EE, uBT_WW, uBT_EE
  real, dimension(:,:,:), contiguous, pointer :: visc_rem_u, u_cor
  real, dimension(:,:), contiguous, pointer :: du_cor
  real, dimension(:,:,:), contiguous, pointer :: u, h_in, h_W, h_E, uh, por_face_areaU
  real, dimension(:,:), contiguous, pointer :: dy_Cu, IareaT, IdxT, areaT, dxT, mask2dCu, dxCu
  real, dimension(:,:,:), contiguous, pointer :: visc_rem_u_tmp
  real, dimension(:,:), contiguous, pointer :: du_max_CFL, du_min_CFL, duhdu_tot_0, uh_tot_0
  real, dimension(:,:), contiguous, pointer :: visc_rem_max, du
  logical, dimension(:,:), contiguous, pointer :: do_I
  type(RealArray_t) :: visc_rem_u_tmp_a
  type(RealArray_t) :: uh_tot_0_a, duhdu_tot_0_a, du_max_CFL_a, du_min_CFL_a, du_a
  type(RealArray_t) :: visc_rem_max_a
  type(LogicalArray_t) :: do_I_a
  ! Never allocated -- the zero-transport (du0) correction never uses uhbt or reports uh_3d.
  type(RealArray_t) :: uhbt_none, uh_3d_none

  call cpu_clock_begin(id_clock_correct)

  call u_a%view(u)
  call h_in_a%view(h_in)
  call h_W_a%view(h_W)
  call h_E_a%view(h_E)
  call uh_a%view(uh)
  call por_face_areaU_a%view(por_face_areaU)
  call dy_Cu_a%view(dy_Cu)
  call IareaT_a%view(IareaT)
  call IdxT_a%view(IdxT)
  call areaT_a%view(areaT)
  call dxT_a%view(dxT)
  call mask2dCu_a%view(mask2dCu)
  call dxCu_a%view(dxCu)

  nullify(FA_u_W0, FA_u_E0, FA_u_WW, FA_u_EE, uBT_WW, uBT_EE)
  nullify(visc_rem_u, u_cor, du_cor)
  if (visc_rem_u_a%associated()) call visc_rem_u_a%view(visc_rem_u)
  if (u_cor_a%associated()) call u_cor_a%view(u_cor)
  if (du_cor_a%associated()) call du_cor_a%view(du_cor)

  use_visc_rem = visc_rem_u_a%associated()

  set_BT_cont = .false. ; if (present(BT_cont)) set_BT_cont = (associated(BT_cont))
  if (set_BT_cont) then
    call BT_cont%FA_u_W0%view(FA_u_W0) ; call BT_cont%FA_u_E0%view(FA_u_E0)
    call BT_cont%FA_u_WW%view(FA_u_WW) ; call BT_cont%FA_u_EE%view(FA_u_EE)
    call BT_cont%uBT_WW%view(uBT_WW)   ; call BT_cont%uBT_EE%view(uBT_EE)
  endif

  local_specified_BC = .false. ; local_Flather_OBC = .false. ; local_open_BC = .false.
  if (associated(OBC)) then ; if (OBC%OBC_pe) then
    local_specified_BC = OBC%specified_u_BCs_exist_globally
    local_Flather_OBC = OBC%Flather_u_BCs_exist_globally
    local_open_BC = OBC%open_u_BCs_exist_globally
  endif ; endif

  ish = bxC%idxS(1) ; ieh = bxC%idxE(1) ; jsh = bxC%idxS(2) ; jeh = bxC%idxE(2) ; nz  = bxC%idxE(3)

  CFL_dt = CS%CFL_limit_adjust / dt
  I_dt = 1.0 / dt
  if (CS%aggress_adjust) CFL_dt = I_dt

  call visc_rem_u_tmp_a%allocView(visc_rem_u_tmp, lb=[u_a%lb(1),u_a%lb(2),u_a%lb(3)], &
                                  ub=[u_a%ub(1),u_a%ub(2),u_a%ub(3)])
  if (uhbt_a%associated() .or. set_BT_cont) then
    call du_max_CFL_a%allocView(du_max_CFL, lb=[u_a%lb(1),u_a%lb(2)], ub=[u_a%ub(1),u_a%ub(2)])
    call du_min_CFL_a%allocView(du_min_CFL, lb=[u_a%lb(1),u_a%lb(2)], ub=[u_a%ub(1),u_a%ub(2)])
    call duhdu_tot_0_a%allocView(duhdu_tot_0, lb=[u_a%lb(1),u_a%lb(2)], ub=[u_a%ub(1),u_a%ub(2)])
    call uh_tot_0_a%allocView(uh_tot_0, lb=[u_a%lb(1),u_a%lb(2)], ub=[u_a%ub(1),u_a%ub(2)])
    call visc_rem_max_a%allocView(visc_rem_max, lb=[u_a%lb(1),u_a%lb(2)], ub=[u_a%ub(1),u_a%ub(2)])
  endif

  do concurrent (j=jsh:jeh)

    if (du_cor_a%associated()) then
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
                     CS%vol_CFL, por_face_areaU(I,j,k))
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

    if (uhbt_a%associated() .or. set_BT_cont) then
      if (use_visc_rem.and.CS%use_visc_rem_max) then
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
        if (CS%vol_CFL) then
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
        if (CS%aggress_adjust) then
          ! untested!
          do k=1,nz ; do concurrent (I=ish-1:ieh)
            if (CS%vol_CFL) then
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
            if (CS%vol_CFL) then
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
        if (CS%aggress_adjust) then
          do k=1,nz ; do concurrent (I=ish-1:ieh)
            if (CS%vol_CFL) then
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
            if (CS%vol_CFL) then
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
    endif ! uhbt_a%associated() .or. set_BT_cont
  enddo

  if (uhbt_a%associated() .or. set_BT_cont) then
    call do_I_a%allocView(do_I, lb=[u_a%lb(1),u_a%lb(2)], ub=[u_a%ub(1),u_a%ub(2)])
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

    call du_a%allocView(du, lb=[u_a%lb(1),u_a%lb(2)], ub=[u_a%ub(1),u_a%ub(2)], source=0.0)

    if (uhbt_a%associated()) then
      ! Find du and uh.
      call zonal_flux_adjust(bxC, u_a, h_in_a, h_W_a, h_E_a, uh_tot_0_a, duhdu_tot_0_a, du_a, &
                            du_max_CFL_a, du_min_CFL_a, dt, dy_Cu_a, IareaT_a, IdxT_a, CS, &
                            visc_rem_u_tmp_a, &
                            do_I_a, por_face_areaU_a, uhbt_a=uhbt_a, uh_3d_a=uh_a, OBC=OBC)

      do concurrent (j=jsh:jeh)
        if (u_cor_a%associated()) then
          do concurrent (k=1:nz, I=ish-1:ieh)
            u_cor(I,j,k) = u(I,j,k) + du(I,j) * visc_rem_u_tmp(I,j,k)
          enddo
          if (any_simple_OBC) then
            ! untested
            do concurrent (k=1:nz, I=ish-1:ieh, simple_OBC_pt(I,j))
              u_cor(I,j,k) = OBC%segment(abs(OBC%segnum_u(I,j)))%normal_vel(I,j,k)
            enddo
          endif
        endif ! u-corrected

        if (du_cor_a%associated()) then
          do concurrent (I=ish-1:ieh)
            du_cor(I,j) = du(I,j)
          enddo
        endif ! du-corrected
      enddo
    endif
    if (set_BT_cont) then
      ! Diagnose the zero-transport correction, du0.
      call zonal_flux_adjust(bxC, u_a, h_in_a, h_W_a, h_E_a, uh_tot_0_a, duhdu_tot_0_a, du_a, &
                            du_max_CFL_a, du_min_CFL_a, dt, dy_Cu_a, IareaT_a, IdxT_a, CS, &
                            visc_rem_u_tmp_a, &
                            do_I_a, por_face_areaU_a, uhbt_a=uhbt_none, uh_3d_a=uh_3d_none)
      call set_zonal_BT_cont(bxC, u_a, h_in_a, h_W_a, h_E_a, BT_cont, du_a, uh_tot_0_a, &
                              duhdu_tot_0_a, du_max_CFL_a, du_min_CFL_a, dt, &
                              dxCu_a, dy_Cu_a, IareaT_a, IdxT_a, CS, &
                              visc_rem_u_tmp_a, &
                              visc_rem_max_a, do_I_a, por_face_areaU_a)
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
            FA_u_W0(I,j) = FAuI ; FA_u_E0(I,j) = FAuI
            FA_u_WW(I,j) = FAuI ; FA_u_EE(I,j) = FAuI
            uBT_WW(I,j) = 0.0 ; uBT_EE(I,j) = 0.0
          endif
        enddo
      endif
    endif
    call du_a%free()
    call uh_tot_0_a%free() ; call duhdu_tot_0_a%free() ; call du_max_CFL_a%free()
    call du_min_CFL_a%free() ; call visc_rem_max_a%free() ; call do_I_a%free()
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
            FA_u_W0(I,j) = FA_u ; FA_u_E0(I,j) = FA_u
            FA_u_WW(I,j) = FA_u ; FA_u_EE(I,j) = FA_u
            uBT_WW(I,j) = 0.0 ; uBT_EE(I,j) = 0.0
          enddo
        else
          do concurrent (j = OBC%segment(n)%HI%Jsd:OBC%segment(n)%HI%Jed)
            FA_u = 0.0
            do k=1,nz ; FA_u = FA_u + h_in(i+1,j,k)*(dy_Cu(I,j)*por_face_areaU(I,j,k)) ; enddo
            FA_u_W0(I,j) = FA_u ; FA_u_E0(I,j) = FA_u
            FA_u_WW(I,j) = FA_u ; FA_u_EE(I,j) = FA_u
            uBT_WW(I,j) = 0.0 ; uBT_EE(I,j) = 0.0
          enddo
        endif
      endif
    enddo
  endif

  if  (set_BT_cont) then ; if (BT_cont%h_u%associated()) then
    if (u_cor_a%associated()) then
      call zonal_flux_thickness(bxC, u_cor_a, h_in_a, h_W_a, h_E_a, BT_cont%h_u, dt, &
                                dy_Cu_a, IareaT_a, IdxT_a, &
                                CS%vol_CFL, &
                                CS%marginal_faces, OBC, por_face_areaU_a, &
                                visc_rem_u_tmp_a)
    else
      call zonal_flux_thickness(bxC, u_a, h_in_a, h_W_a, h_E_a, BT_cont%h_u, dt, &
                                dy_Cu_a, IareaT_a, IdxT_a, &
                                CS%vol_CFL, &
                                CS%marginal_faces, OBC, por_face_areaU_a, &
                                visc_rem_u_tmp_a)
    endif
  endif ; endif

  call visc_rem_u_tmp_a%free()

  call cpu_clock_end(id_clock_correct)

end subroutine zonal_mass_flux


!> Calculates the vertically integrated mass or volume fluxes through the zonal faces.
!> Original Fortran implementation of zonal_BT_mass_flux (renamed). Takes containers.
subroutine zonal_BT_mass_flux_fortran(bxC, u_a, h_in_a, h_W_a, h_E_a, uhbt_a, dt, &
                              dy_Cu_a, IareaT_a, IdxT_a, CS, &
                              OBC, por_face_areaU_a)
  type(Box_t),                                intent(in)  :: bxC  !< Iteration box for continuity solver
  type(RealArray_t),  intent(in)  :: u_a    !< Zonal velocity [L T-1 ~> m s-1]
  type(RealArray_t),  intent(in)  :: h_in_a !< Layer thickness used to
                                                                  !! calculate fluxes [H ~> m or kg m-2]
  type(RealArray_t),  intent(in)  :: h_W_a !< Western edge thickness in the PPM
                                                                  !! reconstruction [H ~> m or kg m-2].
  type(RealArray_t),  intent(in)  :: h_E_a !< Eastern edge thickness in the PPM
                                                                  !! reconstruction [H ~> m or kg m-2].
  type(RealArray_t),  intent(inout) :: uhbt_a !< The summed volume flux through
                                                 !! zonal faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  real,                                       intent(in)  :: dt   !< Time increment [T ~> s].
  type(RealArray_t),  intent(in)  :: dy_Cu_a  !< The grid cell's unblocked lengths of
                                              !! the u-faces of the h-cell [L ~> m].
  type(RealArray_t),  intent(in)  :: IareaT_a !< The grid cell's 1/areaT [L-2 ~> m-2].
  type(RealArray_t),  intent(in)  :: IdxT_a   !< The grid cell's 1/dxT [L-1 ~> m-1].
  type(transport_adjust_CS),                  intent(in)  :: CS !< Options
                       !! controlling the transport adjustment and barotropic-consistency iteration.
  type(ocean_OBC_type),                       pointer     :: OBC  !< Open boundary condition type
                                                                  !! specifies whether, where, and what
                                                                  !! open boundary conditions are used.
  type(RealArray_t), &
                 intent(in)  :: por_face_areaU_a !< fractional open area of U-faces [nondim]

  ! Local variables
  real :: uh(u_a%lb(1):u_a%ub(1), u_a%lb(2):u_a%ub(2), u_a%lb(3):u_a%ub(3))
     ! Volume flux through zonal faces = u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1]
  real :: duhdu(u_a%lb(1):u_a%ub(1), u_a%lb(2):u_a%ub(2), u_a%lb(3):u_a%ub(3))
     ! Partial derivative of uh with u [H L ~> m2 or kg m-1].
  integer :: i, j, k, ish, ieh, jsh, jeh, nz, l_seg
  logical :: local_specified_BC
  logical, dimension(u_a%lb(2):u_a%ub(2)) :: OBC_in_row
  real, dimension(:,:,:), contiguous, pointer :: u, h_in, h_W, h_E, por_face_areaU
  real, dimension(:,:), contiguous, pointer :: uhbt
  real, dimension(:,:), contiguous, pointer :: dy_Cu, IareaT, IdxT

  call u_a%view(u)
  call h_in_a%view(h_in)
  call h_W_a%view(h_W)
  call h_E_a%view(h_E)
  call por_face_areaU_a%view(por_face_areaU)
  call uhbt_a%view(uhbt)
  call dy_Cu_a%view(dy_Cu)
  call IareaT_a%view(IareaT)
  call IdxT_a%view(IdxT)

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
                   IareaT(I+1,j), IdxT(I,j), IdxT(I+1,j), dt, &
                   CS%vol_CFL, por_face_areaU(I,j,k))
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

end subroutine zonal_BT_mass_flux_fortran

!> Shim for zonal_BT_mass_flux -- dispatches via ZONAL_BT_MASS_FLUX_MODE env var.
subroutine zonal_BT_mass_flux(bxC, u_a, h_in_a, h_W_a, h_E_a, uhbt_a, dt, &
                              dy_Cu_a, IareaT_a, IdxT_a, CS, &
                              OBC, por_face_areaU_a)
  type(Box_t),                                intent(in)  :: bxC  !< Iteration box for continuity solver
  type(RealArray_t),  intent(in)  :: u_a    !< Zonal velocity [L T-1 ~> m s-1]
  type(RealArray_t),  intent(in)  :: h_in_a !< Layer thickness used to
                                                                  !! calculate fluxes [H ~> m or kg m-2]
  type(RealArray_t),  intent(in)  :: h_W_a !< Western edge thickness in the PPM
                                                                  !! reconstruction [H ~> m or kg m-2].
  type(RealArray_t),  intent(in)  :: h_E_a !< Eastern edge thickness in the PPM
                                                                  !! reconstruction [H ~> m or kg m-2].
  type(RealArray_t),  intent(inout) :: uhbt_a !< The summed volume flux through
                                                 !! zonal faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  real,                                       intent(in)  :: dt   !< Time increment [T ~> s].
  type(RealArray_t),  intent(in)  :: dy_Cu_a  !< The grid cell's unblocked lengths of
                                              !! the u-faces of the h-cell [L ~> m].
  type(RealArray_t),  intent(in)  :: IareaT_a !< The grid cell's 1/areaT [L-2 ~> m-2].
  type(RealArray_t),  intent(in)  :: IdxT_a   !< The grid cell's 1/dxT [L-1 ~> m-1].
  type(transport_adjust_CS),                  intent(in)  :: CS !< Options
                       !! controlling the transport adjustment and barotropic-consistency iteration.
  type(ocean_OBC_type),                       pointer     :: OBC  !< Open boundary condition type
                                                                  !! specifies whether, where, and what
                                                                  !! open boundary conditions are used.
  type(RealArray_t), &
                 intent(in)  :: por_face_areaU_a !< fractional open area of U-faces [nondim]

  integer :: mode, rc
  type(RealArray_C) :: u_c, h_in_c, h_W_c, h_E_c, uhbt_c, dy_Cu_c, IareaT_c, IdxT_c
  type(RealArray_C) :: por_face_areaU_c
  type(Box_C)                 :: bxC_c
  type(transport_adjust_CS_C) :: CS_c
  type(c_ptr)                 :: OBC_c
  type(io_recorder)  :: rec
  logical            :: capture
  character(len=80)  :: kernel
  character(len=100) :: dir
  character(len=256) :: binFile, metaFile

  kernel = "zonal_bt_mass_flux"

  call cpu_clock_begin(id_clock_correct)

  mode = getenv_mode("ZONAL_BT_MASS_FLUX_MODE", default=TIMH_runFORTRAN)

  select case (mode)

    case (TIMH_capture)
      capture = (.not. already_recorded(trim(kernel))) .and. is_root_pe()
      if (capture) then
        dir = "capture"
        rc = mkdir_posix(trim(dir) // c_null_char, int(o'755', c_int))
        binFile  = trim(dir) // "/" // trim(kernel) // ".bin"
        metaFile = trim(dir) // "/" // trim(kernel) // ".meta"
        call rec%open_write(binFile, metaFile)
        call rec%add("_bxC",            bxC)
        call rec%add("_u",              u_a)
        call rec%add("_h_in",           h_in_a)
        call rec%add("_h_W",            h_W_a)
        call rec%add("_h_E",            h_E_a)
        call rec%add("_uhbt_before",    uhbt_a)
        call rec%add("_dt",             dt)
        call rec%add("_dy_Cu",          dy_Cu_a)
        call rec%add("_IareaT",         IareaT_a)
        call rec%add("_IdxT",           IdxT_a)
        call rec%add("_vol_CFL",        CS%vol_CFL)
        call rec%add("_por_face_areaU", por_face_areaU_a)
      endif

      call zonal_BT_mass_flux_fortran(bxC, u_a, h_in_a, h_W_a, h_E_a, uhbt_a, dt, &
                                      dy_Cu_a, IareaT_a, IdxT_a, CS, &
                                      OBC, por_face_areaU_a)

      if (capture) then
        call rec%add("_uhbt_after", uhbt_a)
        call rec%close()
        call mark_recorded(trim(kernel))
      endif

#ifdef _TIM
    case (TIMH_runAMREX)
      bxC_c           = bxC%to_c()
      u_c             = u_a%to_c()
      h_in_c          = h_in_a%to_c()
      h_W_c           = h_W_a%to_c()
      h_E_c           = h_E_a%to_c()
      uhbt_c          = uhbt_a%to_c()
      dy_Cu_c         = dy_Cu_a%to_c()
      IareaT_c        = IareaT_a%to_c()
      IdxT_c          = IdxT_a%to_c()
      CS_c            = transport_adjust_CS_to_c(CS)
      por_face_areaU_c = por_face_areaU_a%to_c()
      if (associated(OBC)) then
        OBC_c = c_loc(OBC)
      else
        OBC_c = c_null_ptr
      endif
      call turbotmp_zonal_bt_mass_flux_bridge(bxC_c, u_c, h_in_c, h_W_c, h_E_c, uhbt_c, dt, &
                                              dy_Cu_c, IareaT_c, IdxT_c, CS_c, OBC_c, &
                                              por_face_areaU_c)
#endif

    case default
      call zonal_BT_mass_flux_fortran(bxC, u_a, h_in_a, h_W_a, h_E_a, uhbt_a, dt, &
                                      dy_Cu_a, IareaT_a, IdxT_a, CS, &
                                      OBC, por_face_areaU_a)

  end select

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
!> Original Fortran implementation of zonal_flux_thickness (renamed). Takes containers.
subroutine zonal_flux_thickness_fortran(bxC, u_a, h_a, h_W_a, h_E_a, h_u_a, dt, &
                                dy_Cu_a, IareaT_a, IdxT_a, vol_CFL, &
                                marginal, OBC, por_face_areaU_a, visc_rem_u_a)
  type(box_t),                               intent(in)    :: bxC  !< Iteration box for continuity solver
  type(RealArray_t),  intent(in)  :: u_a   !< Zonal velocity [L T-1 ~> m s-1].
  type(RealArray_t),  intent(in)  :: h_a  !< Layer thickness used to
                                                                   !! calculate fluxes [H ~> m or kg m-2].
  type(RealArray_t),  intent(in)  :: h_W_a !< West edge thickness in the
                                                                   !! reconstruction [H ~> m or kg m-2].
  type(RealArray_t),  intent(in)  :: h_E_a !< East edge thickness in the
                                                                   !! reconstruction [H ~> m or kg m-2].
  type(RealArray_t),  intent(inout):: h_u_a !< Effective thickness at zonal faces,
                                                                   !! scaled down to account for the effects of
                                                                   !! viscosity and the fractional open area
                                                                   !! [H ~> m or kg m-2].
  real,                                      intent(in)    :: dt   !< Time increment [T ~> s].
  type(RealArray_t),  intent(in)  :: dy_Cu_a  !< The grid cell's unblocked lengths of
                                              !! the u-faces of the h-cell [L ~> m].
  type(RealArray_t),  intent(in)  :: IareaT_a !< The grid cell's 1/areaT [L-2 ~> m-2].
  type(RealArray_t),  intent(in)  :: IdxT_a   !< The grid cell's 1/dxT [L-1 ~> m-1].
  logical,                                   intent(in)    :: vol_CFL !< If true, rescale the ratio
                          !! of face areas to the cell areas when estimating the CFL number.
  logical,                                   intent(in)    :: marginal !< If true, report the
                          !! marginal face thicknesses; otherwise report transport-averaged thicknesses.
  type(RealArray_t), &
                                   intent(in)    :: por_face_areaU_a !< fractional open area of
                                                                     !! U-faces [nondim]
  type(ocean_OBC_type),                      pointer       :: OBC !< Open boundaries control structure.
  type(RealArray_t), &
                                             intent(in)    :: visc_rem_u_a
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
  real, dimension(:,:,:), contiguous, pointer :: visc_rem_u
  real, dimension(:,:,:), contiguous, pointer :: u, h, h_W, h_E, h_u, por_face_areaU
  real, dimension(:,:), contiguous, pointer :: dy_Cu, IareaT, IdxT

  nullify(visc_rem_u)
  if (visc_rem_u_a%associated()) call visc_rem_u_a%view(visc_rem_u)
  call u_a%view(u)
  call h_a%view(h)
  call h_W_a%view(h_W)
  call h_E_a%view(h_E)
  call h_u_a%view(h_u)
  call por_face_areaU_a%view(por_face_areaU)
  call dy_Cu_a%view(dy_Cu)
  call IareaT_a%view(IareaT)
  call IdxT_a%view(IdxT)

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

end subroutine zonal_flux_thickness_fortran

!> Shim for zonal_flux_thickness -- dispatches via ZONAL_FLUX_THICKNESS_MODE env var.
subroutine zonal_flux_thickness(bxC, u_a, h_a, h_W_a, h_E_a, h_u_a, dt, &
                                dy_Cu_a, IareaT_a, IdxT_a, vol_CFL, &
                                marginal, OBC, por_face_areaU_a, visc_rem_u_a)
  type(box_t),                               intent(in)    :: bxC  !< Iteration box for continuity solver
  type(RealArray_t),  intent(in)  :: u_a   !< Zonal velocity [L T-1 ~> m s-1].
  type(RealArray_t),  intent(in)  :: h_a  !< Layer thickness used to
                                                                   !! calculate fluxes [H ~> m or kg m-2].
  type(RealArray_t),  intent(in)  :: h_W_a !< West edge thickness in the
                                                                   !! reconstruction [H ~> m or kg m-2].
  type(RealArray_t),  intent(in)  :: h_E_a !< East edge thickness in the
                                                                   !! reconstruction [H ~> m or kg m-2].
  type(RealArray_t),  intent(inout):: h_u_a !< Effective thickness at zonal faces,
                                                                   !! scaled down to account for the effects of
                                                                   !! viscosity and the fractional open area
                                                                   !! [H ~> m or kg m-2].
  real,                                      intent(in)    :: dt   !< Time increment [T ~> s].
  type(RealArray_t),  intent(in)  :: dy_Cu_a  !< The grid cell's unblocked lengths of
                                              !! the u-faces of the h-cell [L ~> m].
  type(RealArray_t),  intent(in)  :: IareaT_a !< The grid cell's 1/areaT [L-2 ~> m-2].
  type(RealArray_t),  intent(in)  :: IdxT_a   !< The grid cell's 1/dxT [L-1 ~> m-1].
  logical,                                   intent(in)    :: vol_CFL !< If true, rescale the ratio
                          !! of face areas to the cell areas when estimating the CFL number.
  logical,                                   intent(in)    :: marginal !< If true, report the
                          !! marginal face thicknesses; otherwise report transport-averaged thicknesses.
  type(RealArray_t), &
                                   intent(in)    :: por_face_areaU_a !< fractional open area of
                                                                     !! U-faces [nondim]
  type(ocean_OBC_type),                      pointer       :: OBC !< Open boundaries control structure.
  type(RealArray_t), &
                                             intent(in)    :: visc_rem_u_a
                          !< Both the fraction of the momentum originally in a layer that remains after
                          !! a time-step of viscosity, and the fraction of a time-step's worth of a
                          !! barotropic acceleration that a layer experiences after viscosity is applied [nondim].
                          !! Visc_rem_u is between 0 (at the bottom) and 1 (far above the bottom).

  integer :: mode, rc
  type(RealArray_C) :: u_c, h_c, h_W_c, h_E_c, h_u_c, dy_Cu_c, IareaT_c, IdxT_c
  type(RealArray_C) :: por_face_areaU_c, visc_rem_u_c
  type(Box_C)         :: bxC_c
  type(c_ptr)         :: OBC_c
  logical(c_bool)     :: vol_CFL_c, marginal_c
  type(io_recorder)   :: rec
  logical             :: capture
  character(len=80)   :: kernel
  character(len=100)  :: dir
  character(len=256)  :: binFile, metaFile

  kernel = "zonal_flux_thickness"

  mode = getenv_mode("ZONAL_FLUX_THICKNESS_MODE", default=TIMH_runFORTRAN)

  select case (mode)

    case (TIMH_capture)
      capture = (.not. already_recorded(trim(kernel))) .and. is_root_pe()
      if (capture) then
        dir = "capture"
        rc = mkdir_posix(trim(dir) // c_null_char, int(o'755', c_int))
        binFile  = trim(dir) // "/" // trim(kernel) // ".bin"
        metaFile = trim(dir) // "/" // trim(kernel) // ".meta"
        call rec%open_write(binFile, metaFile)
        call rec%add("_bxC",            bxC)
        call rec%add("_u",              u_a)
        call rec%add("_h",              h_a)
        call rec%add("_h_W",            h_W_a)
        call rec%add("_h_E",            h_E_a)
        call rec%add("_h_u_before",     h_u_a)
        call rec%add("_dt",             dt)
        call rec%add("_dy_Cu",          dy_Cu_a)
        call rec%add("_IareaT",         IareaT_a)
        call rec%add("_IdxT",           IdxT_a)
        call rec%add("_vol_CFL",        vol_CFL)
        call rec%add("_marginal",       marginal)
        call rec%add("_por_face_areaU", por_face_areaU_a)
        if (visc_rem_u_a%associated()) call rec%add("_visc_rem_u", visc_rem_u_a)
      endif

      call zonal_flux_thickness_fortran(bxC, u_a, h_a, h_W_a, h_E_a, h_u_a, dt, &
                                        dy_Cu_a, IareaT_a, IdxT_a, vol_CFL, &
                                        marginal, OBC, por_face_areaU_a, visc_rem_u_a)

      if (capture) then
        call rec%add("_h_u_after", h_u_a)
        call rec%close()
        call mark_recorded(trim(kernel))
      endif

#ifdef _TIM
    case (TIMH_runAMREX)
      bxC_c            = bxC%to_c()
      u_c              = u_a%to_c()
      h_c              = h_a%to_c()
      h_W_c            = h_W_a%to_c()
      h_E_c            = h_E_a%to_c()
      h_u_c            = h_u_a%to_c()
      dy_Cu_c          = dy_Cu_a%to_c()
      IareaT_c         = IareaT_a%to_c()
      IdxT_c           = IdxT_a%to_c()
      vol_CFL_c        = vol_CFL
      marginal_c       = marginal
      por_face_areaU_c = por_face_areaU_a%to_c()
      visc_rem_u_c     = visc_rem_u_a%to_c()
      if (associated(OBC)) then
        OBC_c = c_loc(OBC)
      else
        OBC_c = c_null_ptr
      endif
      call turbotmp_zonal_flux_thickness_bridge(bxC_c, u_c, h_c, h_W_c, h_E_c, h_u_c, dt, &
                                                dy_Cu_c, IareaT_c, IdxT_c, vol_CFL_c, &
                                                marginal_c, OBC_c, por_face_areaU_c, visc_rem_u_c)
#endif

    case default
      call zonal_flux_thickness_fortran(bxC, u_a, h_a, h_W_a, h_E_a, h_u_a, dt, &
                                        dy_Cu_a, IareaT_a, IdxT_a, vol_CFL, &
                                        marginal, OBC, por_face_areaU_a, visc_rem_u_a)

  end select

end subroutine zonal_flux_thickness

!> Returns the barotropic velocity adjustment that gives the
!! desired barotropic (layer-summed) transport.
!> Original Fortran implementation of zonal_flux_adjust (renamed). Takes containers.
subroutine zonal_flux_adjust_fortran(bxC, u_a, h_in_a, h_W_a, h_E_a, uh_tot_0_a, duhdu_tot_0_a, &
                             du_a, du_max_CFL_a, du_min_CFL_a, dt, dy_Cu_a, IareaT_a, IdxT_a, CS, &
                             visc_rem_a, do_I_in_a, por_face_areaU_a, uhbt_a, uh_3d_a, OBC)
  type(box_t),                                intent(in)    :: bxC  !< Iteration box for continuity solver
  type(RealArray_t),  intent(in)   :: u_a   !< Zonal velocity [L T-1 ~> m s-1].
  type(RealArray_t),  intent(in)   :: h_in_a !< Layer thickness used to
                                                                    !! calculate fluxes [H ~> m or kg m-2].
  type(RealArray_t),  intent(in)   :: h_W_a !< West edge thickness in the
                                                                    !! reconstruction [H ~> m or kg m-2].
  type(RealArray_t),  intent(in)   :: h_E_a !< East edge thickness in the
                                                                    !! reconstruction [H ~> m or kg m-2].
  type(RealArray_t),  intent(in)   :: visc_rem_a !< Both the fraction of the
                       !! momentum originally in a layer that remains after a time-step of viscosity, and
                       !! the fraction of a time-step's worth of a barotropic acceleration that a layer
                       !! experiences after viscosity is applied [nondim].
                       !! Visc_rem is between 0 (at the bottom) and 1 (far above the bottom).
  type(RealArray_t), intent(in)    :: uhbt_a !< The summed volume flux
                       !! through zonal faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(RealArray_t),  intent(in)   :: du_max_CFL_a  !< Maximum acceptable
                       !! value of du [L T-1 ~> m s-1].
  type(RealArray_t),  intent(in)   :: du_min_CFL_a  !< Minimum acceptable
                       !! value of du [L T-1 ~> m s-1].
  type(RealArray_t),  intent(in)   :: uh_tot_0_a    !< The summed transport
                       !! with 0 adjustment [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(RealArray_t),  intent(in)   :: duhdu_tot_0_a !< The partial derivative
                       !! of du_err with du at 0 adjustment [H L ~> m2 or kg m-1].
  type(RealArray_t),                          intent(inout) :: du_a !<
                       !! The barotropic velocity adjustment [L T-1 ~> m s-1].
  real,                                       intent(in)    :: dt  !< Time increment [T ~> s].
  type(RealArray_t),                          intent(in)    :: dy_Cu_a  !< The grid cell's unblocked
                       !! lengths of the u-faces of the h-cell [L ~> m].
  type(RealArray_t),                          intent(in)    :: IareaT_a !< The grid cell's 1/areaT
                       !! [L-2 ~> m-2].
  type(RealArray_t),                          intent(in)    :: IdxT_a   !< The grid cell's 1/dxT
                       !! [L-1 ~> m-1].
  type(transport_adjust_CS),           intent(in)    :: CS !< Options
                       !! controlling the transport adjustment and barotropic-consistency iteration.


  type(LogicalArray_t),  intent(in)  :: do_I_in_a !< A logical flag indicating
                                                                       !! which I values to work on.
  type(RealArray_t),  intent(in)  :: por_face_areaU_a !< fractional open area
                                                                              !! of U-faces [nondim].
  type(RealArray_t), &
                                              intent(inout) :: uh_3d_a !< Volume flux through zonal
                                                 !! faces = u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(ocean_OBC_type),             optional, pointer       :: OBC !< Open boundaries control structure.
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
  logical :: do_I(u_a%lb(1):u_a%ub(1))
  logical :: local_OBC, use_uhbt, use_uh_3d
  integer, parameter:: max_itts = 20
  real, dimension(:,:), contiguous, pointer :: uhbt
  real, dimension(:,:,:), contiguous, pointer :: uh_3d
  real, dimension(:,:,:), contiguous, pointer :: u, h_in, h_W, h_E, visc_rem, por_face_areaU
  real, dimension(:,:), contiguous, pointer :: du_max_CFL, du_min_CFL, uh_tot_0, duhdu_tot_0, du
  real, dimension(:,:), contiguous, pointer :: dy_Cu, IareaT, IdxT
  logical, dimension(:,:), contiguous, pointer :: do_I_in

  nullify(uhbt, uh_3d)
  if (uhbt_a%associated()) call uhbt_a%view(uhbt)
  if (uh_3d_a%associated()) call uh_3d_a%view(uh_3d)

  call u_a%view(u)
  call h_in_a%view(h_in)
  call h_W_a%view(h_W)
  call h_E_a%view(h_E)
  call visc_rem_a%view(visc_rem)
  call du_max_CFL_a%view(du_max_CFL)
  call du_min_CFL_a%view(du_min_CFL)
  call uh_tot_0_a%view(uh_tot_0)
  call duhdu_tot_0_a%view(duhdu_tot_0)
  call du_a%view(du)
  call do_I_in_a%view(do_I_in)
  call por_face_areaU_a%view(por_face_areaU)
  call dy_Cu_a%view(dy_Cu)
  call IareaT_a%view(IareaT)
  call IdxT_a%view(IdxT)

  local_OBC = .false.
  if (present(OBC)) then
    if (associated(OBC)) then
      local_OBC = OBC%open_u_BCs_exist_globally
    endif
  endif

  use_uhbt = uhbt_a%associated()
  use_uh_3d = uh_3d_a%associated()

  ish = bxC%idxS(1) ; ieh = bxC%idxE(1) ; jsh = bxC%idxS(2) ; jeh = bxC%idxE(2) ; nz  = bxC%idxE(3)

  tol_vel = CS%tol_vel

  ! NVIDIA needs private arrays to be alloc'ed to prevent data transfers.
  ! GCC doesn't understand map(alloc: ...) for variables also marked private
  !$omp target enter data map(alloc: do_I, du_max, du_min, duhdu_tot, uh_err, uh_err_best, uh_aux)

  ! NVIDIA do concurrent doesn't work with private arrays (private scalars OK)
  !$omp target teams loop &
  !$omp   private(uh_err, uh_err_best, duhdu_tot, du_min, du_max, do_I, uh_aux, itt, tol_eta)
  do j=jsh,jeh

    if (use_uh_3d) then
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
        case (:1) ; tol_eta = 1e-6 * CS%tol_eta
        case (2)  ; tol_eta = 1e-4 * CS%tol_eta
        case (3)  ; tol_eta = 1e-2 * CS%tol_eta
        case default ; tol_eta = CS%tol_eta
      end select

      do concurrent (I=ish-1:ieh, do_I(I)) &
          & DO_LOCALITY(local(ddu, du_prev))
        if (uh_err(I) > 0.0) then ; du_max(I) = du(I,j)
        elseif (uh_err(I) < 0.0) then ; du_min(I) = du(I,j)
        else ; do_I(I) = .false. ; endif
        if ((dt * min(IareaT(i,j),IareaT(i+1,j))*abs(uh_err(I)) > tol_eta) .or. &
            (CS%better_iter .and. &
             ((abs(uh_err(I)) > tol_vel * duhdu_tot(I)) .or. &
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

      if ((itt < max_itts) .or. use_uh_3d) then
        do concurrent (I=ish-1:ieh)
          uh_err(I) = 0.0 ; duhdu_tot(I) = 0.0
          if (use_uhbt) uh_err(I) = -uhbt(I,j)
        enddo
        do k=1,nz ; do concurrent (I=ish-1:ieh, do_I(I)) DO_LOCALITY(local(u_new, duhdu))
          u_new = u(I,j,k) + du(I,j) * visc_rem(I,j,k)
          call flux_elem(u_new, h_in(I,j,k), h_in(I+1,j,k), h_W(I,j,k), h_W(I+1,j,k), h_E(I,j,k), &
                         h_E(I+1,j,k), uh_aux(I,k), duhdu, visc_rem(I,j,k), dy_Cu(I,j), &
                         IareaT(I,j), IareaT(I+1,j), IdxT(I,j), IdxT(i+1,j), dt, &
                         CS%vol_CFL, por_face_areaU(I,j,k))
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
    if (use_uh_3d) then
      do concurrent (k=1:nz, I=ish-1:ieh)
        uh_3d(I,j,k) = uh_aux(I,k)
      enddo
    endif
  enddo ! j-loop
  ! If there are any faces which have not converged to within the tolerance,
  ! so-be-it, or else use a final upwind correction?
  ! This never seems to happen with 20 iterations as max_itt.

  !$omp target exit data map(release: do_I, du_max, du_min, duhdu_tot, uh_err, uh_err_best, uh_aux)

end subroutine zonal_flux_adjust_fortran

!> Shim for zonal_flux_adjust -- dispatches via ZONAL_FLUX_ADJUST_MODE env var.
subroutine zonal_flux_adjust(bxC, u_a, h_in_a, h_W_a, h_E_a, uh_tot_0_a, duhdu_tot_0_a, &
                             du_a, du_max_CFL_a, du_min_CFL_a, dt, dy_Cu_a, IareaT_a, IdxT_a, CS, &
                             visc_rem_a, do_I_in_a, por_face_areaU_a, uhbt_a, uh_3d_a, OBC)
  type(box_t),                                intent(in)    :: bxC  !< Iteration box for continuity solver
  type(RealArray_t),  intent(in)   :: u_a   !< Zonal velocity [L T-1 ~> m s-1].
  type(RealArray_t),  intent(in)   :: h_in_a !< Layer thickness used to
                                                                    !! calculate fluxes [H ~> m or kg m-2].
  type(RealArray_t),  intent(in)   :: h_W_a !< West edge thickness in the
                                                                    !! reconstruction [H ~> m or kg m-2].
  type(RealArray_t),  intent(in)   :: h_E_a !< East edge thickness in the
                                                                    !! reconstruction [H ~> m or kg m-2].
  type(RealArray_t),  intent(in)   :: visc_rem_a !< Both the fraction of the
                       !! momentum originally in a layer that remains after a time-step of viscosity, and
                       !! the fraction of a time-step's worth of a barotropic acceleration that a layer
                       !! experiences after viscosity is applied [nondim].
                       !! Visc_rem is between 0 (at the bottom) and 1 (far above the bottom).
  type(RealArray_t), intent(in)    :: uhbt_a !< The summed volume flux
                       !! through zonal faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(RealArray_t),  intent(in)   :: du_max_CFL_a  !< Maximum acceptable
                       !! value of du [L T-1 ~> m s-1].
  type(RealArray_t),  intent(in)   :: du_min_CFL_a  !< Minimum acceptable
                       !! value of du [L T-1 ~> m s-1].
  type(RealArray_t),  intent(in)   :: uh_tot_0_a    !< The summed transport
                       !! with 0 adjustment [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(RealArray_t),  intent(in)   :: duhdu_tot_0_a !< The partial derivative
                       !! of du_err with du at 0 adjustment [H L ~> m2 or kg m-1].
  type(RealArray_t),                          intent(inout) :: du_a !<
                       !! The barotropic velocity adjustment [L T-1 ~> m s-1].
  real,                                       intent(in)    :: dt  !< Time increment [T ~> s].
  type(RealArray_t),                          intent(in)    :: dy_Cu_a  !< The grid cell's unblocked
                       !! lengths of the u-faces of the h-cell [L ~> m].
  type(RealArray_t),                          intent(in)    :: IareaT_a !< The grid cell's 1/areaT
                       !! [L-2 ~> m-2].
  type(RealArray_t),                          intent(in)    :: IdxT_a   !< The grid cell's 1/dxT
                       !! [L-1 ~> m-1].
  type(transport_adjust_CS),           intent(in)    :: CS !< Options
                       !! controlling the transport adjustment and barotropic-consistency iteration.


  type(LogicalArray_t),  intent(in)  :: do_I_in_a !< A logical flag indicating
                                                                       !! which I values to work on.
  type(RealArray_t),  intent(in)  :: por_face_areaU_a !< fractional open area
                                                                              !! of U-faces [nondim].
  type(RealArray_t), &
                                              intent(inout) :: uh_3d_a !< Volume flux through zonal
                                                 !! faces = u*h*dy [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(ocean_OBC_type),             optional, pointer       :: OBC !< Open boundaries control structure.

  integer :: mode, rc
  type(RealArray_C) :: u_c, h_in_c, h_W_c, h_E_c, uh_tot_0_c, duhdu_tot_0_c, du_c
  type(RealArray_C) :: du_max_CFL_c, du_min_CFL_c, dy_Cu_c, IareaT_c, IdxT_c
  type(RealArray_C) :: visc_rem_c, por_face_areaU_c, uhbt_c, uh_3d_c
  type(LogicalArray_C)        :: do_I_in_c
  type(Box_C)                 :: bxC_c
  type(transport_adjust_CS_C) :: CS_c
  type(c_ptr)                 :: OBC_c
  type(io_recorder)  :: rec
  logical            :: capture
  character(len=80)  :: kernel
  character(len=100) :: dir
  character(len=256) :: binFile, metaFile

  kernel = "zonal_flux_adjust"

  mode = getenv_mode("ZONAL_FLUX_ADJUST_MODE", default=TIMH_runFORTRAN)

  select case (mode)

    case (TIMH_capture)
      capture = (.not. already_recorded(trim(kernel))) .and. is_root_pe()
      if (capture) then
        dir = "capture"
        rc = mkdir_posix(trim(dir) // c_null_char, int(o'755', c_int))
        binFile  = trim(dir) // "/" // trim(kernel) // ".bin"
        metaFile = trim(dir) // "/" // trim(kernel) // ".meta"
        call rec%open_write(binFile, metaFile)
        call rec%add("_bxC",            bxC)
        call rec%add("_u",              u_a)
        call rec%add("_h_in",           h_in_a)
        call rec%add("_h_W",            h_W_a)
        call rec%add("_h_E",            h_E_a)
        call rec%add("_uh_tot_0",       uh_tot_0_a)
        call rec%add("_duhdu_tot_0",    duhdu_tot_0_a)
        call rec%add("_du_before",      du_a)
        call rec%add("_du_max_CFL",     du_max_CFL_a)
        call rec%add("_du_min_CFL",     du_min_CFL_a)
        call rec%add("_dt",             dt)
        call rec%add("_dy_Cu",          dy_Cu_a)
        call rec%add("_IareaT",         IareaT_a)
        call rec%add("_IdxT",           IdxT_a)
        call rec%add("_tol_eta",        CS%tol_eta)
        call rec%add("_tol_vel",        CS%tol_vel)
        call rec%add("_better_iter",    CS%better_iter)
        call rec%add("_vol_CFL",        CS%vol_CFL)
        call rec%add("_visc_rem",       visc_rem_a)
        call rec%add("_do_I_in",        do_I_in_a)
        call rec%add("_por_face_areaU", por_face_areaU_a)
        if (uhbt_a%associated()) call rec%add("_uhbt", uhbt_a)
        if (uh_3d_a%associated()) call rec%add("_uh_3d_before", uh_3d_a)
      endif

      call zonal_flux_adjust_fortran(bxC, u_a, h_in_a, h_W_a, h_E_a, uh_tot_0_a, duhdu_tot_0_a, &
                           du_a, du_max_CFL_a, du_min_CFL_a, dt, dy_Cu_a, IareaT_a, IdxT_a, CS, &
                           visc_rem_a, do_I_in_a, por_face_areaU_a, uhbt_a, uh_3d_a, OBC)

      if (capture) then
        call rec%add("_du_after", du_a)
        if (uh_3d_a%associated()) call rec%add("_uh_3d_after", uh_3d_a)
        call rec%close()
        call mark_recorded(trim(kernel))
      endif

#ifdef _TIM
    case (TIMH_runAMREX)
      bxC_c            = bxC%to_c()
      u_c              = u_a%to_c()
      h_in_c           = h_in_a%to_c()
      h_W_c            = h_W_a%to_c()
      h_E_c            = h_E_a%to_c()
      uh_tot_0_c       = uh_tot_0_a%to_c()
      duhdu_tot_0_c    = duhdu_tot_0_a%to_c()
      du_c             = du_a%to_c()
      du_max_CFL_c     = du_max_CFL_a%to_c()
      du_min_CFL_c     = du_min_CFL_a%to_c()
      dy_Cu_c          = dy_Cu_a%to_c()
      IareaT_c         = IareaT_a%to_c()
      IdxT_c           = IdxT_a%to_c()
      CS_c             = transport_adjust_CS_to_c(CS)
      visc_rem_c       = visc_rem_a%to_c()
      do_I_in_c        = do_I_in_a%to_c()
      por_face_areaU_c = por_face_areaU_a%to_c()
      uhbt_c           = uhbt_a%to_c()
      uh_3d_c          = uh_3d_a%to_c()
      if (present(OBC)) then
        if (associated(OBC)) then
          OBC_c = c_loc(OBC)
        else
          OBC_c = c_null_ptr
        endif
      else
        OBC_c = c_null_ptr
      endif
      call turbotmp_zonal_flux_adjust_bridge(bxC_c, u_c, h_in_c, h_W_c, h_E_c, uh_tot_0_c, &
                                             duhdu_tot_0_c, du_c, du_max_CFL_c, du_min_CFL_c, dt, &
                                             dy_Cu_c, IareaT_c, IdxT_c, CS_c, visc_rem_c, &
                                             do_I_in_c, por_face_areaU_c, uhbt_c, uh_3d_c, OBC_c)
#endif

    case default
      call zonal_flux_adjust_fortran(bxC, u_a, h_in_a, h_W_a, h_E_a, uh_tot_0_a, duhdu_tot_0_a, &
                           du_a, du_max_CFL_a, du_min_CFL_a, dt, dy_Cu_a, IareaT_a, IdxT_a, CS, &
                           visc_rem_a, do_I_in_a, por_face_areaU_a, uhbt_a, uh_3d_a, OBC)

  end select

end subroutine zonal_flux_adjust


!> Sets a structure that describes the zonal barotropic volume or mass fluxes as a
!! function of barotropic flow to agree closely with the sum of the layer's transports.
!> Original Fortran implementation of set_zonal_BT_cont (renamed). Takes containers.
subroutine set_zonal_BT_cont_fortran(bxC, u_a, h_in_a, h_W_a, h_E_a, BT_cont, du0_a, uh_tot_0_a, &
                             duhdu_tot_0_a, du_max_CFL_a, du_min_CFL_a, dt, &
                             dxCu_a, dy_Cu_a, IareaT_a, IdxT_a, CS, &
                             visc_rem_a, visc_rem_max_a, do_I_a, por_face_areaU_a)
  type(box_t),             intent(in) :: bxC  !< Iteration box for continuity solver
  type(RealArray_t),       intent(in) :: u_a   !< Zonal velocity [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in) :: h_in_a !< Layer thickness used to calculate
                                                !! fluxes [H ~> m or kg m-2].
  type(RealArray_t),       intent(in) :: h_W_a !< West edge thickness in the reconstruction
                                                !! [H ~> m or kg m-2].
  type(RealArray_t),       intent(in) :: h_E_a !< East edge thickness in the reconstruction
                                                !! [H ~> m or kg m-2].
  type(BT_cont_type),   intent(inout) :: BT_cont !< A structure with elements
                       !! that describe the effective open face areas as a function of barotropic flow.
  type(RealArray_t),       intent(in) :: du0_a  !< The barotropic velocity increment that gives 0
                                                 !! transport [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in) :: uh_tot_0_a    !< The summed transport with 0 adjustment
                                                        !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(RealArray_t),       intent(in) :: duhdu_tot_0_a !< The partial derivative
                       !! of du_err with du at 0 adjustment [H L ~> m2 or kg m-1].
  type(RealArray_t),       intent(in) :: du_max_CFL_a  !< Maximum acceptable value of
                                                        !! du [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in) :: du_min_CFL_a  !< Minimum acceptable value of
                                                        !! du [L T-1 ~> m s-1].
  real,                    intent(in) :: dt   !< Time increment [T ~> s].
  type(RealArray_t),       intent(in) :: dxCu_a !< The grid cell's u-point x-extent [L ~> m].
  type(RealArray_t),       intent(in) :: dy_Cu_a !< The grid cell's unblocked lengths of the
                                                !! u-faces of the h-cell [L ~> m].
  type(RealArray_t),       intent(in) :: IareaT_a !< The grid cell's 1/areaT [L-2 ~> m-2].
  type(RealArray_t),       intent(in) :: IdxT_a   !< The grid cell's 1/dxT [L-1 ~> m-1].
  type(transport_adjust_CS), intent(in) :: CS !< Options controlling the
                       !! transport adjustment and barotropic-consistency iteration.
  type(RealArray_t),       intent(in) :: visc_rem_a !< Both the fraction of the
                       !! momentum originally in a layer that remains after a time-step of viscosity, and
                       !! the fraction of a time-step's worth of a barotropic acceleration that a layer
                       !! experiences after viscosity is applied [nondim].
                       !! Visc_rem is between 0 (at the bottom) and 1 (far above the bottom).
  type(RealArray_t),       intent(in) :: visc_rem_max_a !< Maximum allowable visc_rem [nondim].
  type(LogicalArray_t),    intent(in) :: do_I_a   !< A logical flag indicating
                                                        !! which I values to work on.
  type(RealArray_t),       intent(in) :: por_face_areaU_a !< fractional open area
                                                        !! of U-faces [nondim]
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
  real, dimension(:,:), contiguous, pointer :: du0, visc_rem_max
  real, dimension(:,:), contiguous, pointer :: dxCu, dy_Cu, IareaT, IdxT
  logical, dimension(:,:), contiguous, pointer :: do_I
  real, dimension(:,:), contiguous, pointer :: FA_u_W0, FA_u_E0, FA_u_WW, FA_u_EE, uBT_WW, uBT_EE

  call u_a%view(u)
  call h_in_a%view(h_in)
  call h_W_a%view(h_W)
  call h_E_a%view(h_E)
  call du0_a%view(du0)
  call visc_rem_a%view(visc_rem)
  call visc_rem_max_a%view(visc_rem_max)
  call do_I_a%view(do_I)
  call por_face_areaU_a%view(por_face_areaU)
  call dxCu_a%view(dxCu)
  call dy_Cu_a%view(dy_Cu)
  call IareaT_a%view(IareaT)
  call IdxT_a%view(IdxT)

  ish = bxC%idxS(1) ; ieh = bxC%idxE(1) ; jsh = bxC%idxS(2) ; jeh = bxC%idxE(2) ; nz  = bxC%idxE(3)
  Idt = 1.0 / dt
  min_visc_rem = 0.1 ; CFL_min = 1e-6

  call BT_cont%FA_u_W0%view(FA_u_W0) ; call BT_cont%FA_u_E0%view(FA_u_E0)
  call BT_cont%FA_u_WW%view(FA_u_WW) ; call BT_cont%FA_u_EE%view(FA_u_EE)
  call BT_cont%uBT_WW%view(uBT_WW)   ; call BT_cont%uBT_EE%view(uBT_EE)

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
                     CS%vol_CFL, por_face_areaU(I,j,k))
      call flux_elem(u_L, h_in(I,j,k), h_in(I+1,j,k), h_W(I,j,k), h_W(I+1,j,k), h_E(I,j,k), &
                     h_E(I+1,j,k), uh_L, duhdu_L, visc_rem(I,j,k), dy_Cu(I,j), &
                     IareaT(I,j), IareaT(I+1,j), IdxT(I,j), IdxT(i+1,j), dt, &
                     CS%vol_CFL, por_face_areaU(I,j,k))
      call flux_elem(u_R, h_in(I,j,k), h_in(I+1,j,k), h_W(I,j,k), h_W(I+1,j,k), h_E(I,j,k), &
                     h_E(I+1,j,k), uh_R, duhdu_R, visc_rem(I,j,k), dy_Cu(I,j), &
                     IareaT(I,j), IareaT(I+1,j), IdxT(I,j), IdxT(i+1,j), dt, &
                     CS%vol_CFL, por_face_areaU(I,j,k))
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

        FA_u_W0(I,j) = FA_0 ; FA_u_WW(I,j) = FAmt_L(I)
        if (abs(FA_0-FAmt_L(I)) <= 1e-12*FA_0) then ; uBT_WW(I,j) = 0.0 ; else
          uBT_WW(I,j) = (1.5 * (duL(I) - du0(I,j))) * &
                                ((FAmt_L(I) - FA_avg) / (FAmt_L(I) - FA_0))
        endif

        FA_0 = FAmt_0(I) ; FA_avg = FAmt_0(I)
        if ((duR(I) - du0(I,j)) /= 0.0) &
          FA_avg = uhtot_R(I) / (duR(I) - du0(I,j))
        if (FA_avg > max(FA_0, FAmt_R(I))) then ; FA_avg = max(FA_0, FAmt_R(I))
        elseif (FA_avg < min(FA_0, FAmt_R(I))) then ; FA_0 = FA_avg ; endif

        FA_u_E0(I,j) = FA_0 ; FA_u_EE(I,j) = FAmt_R(I)
        if (abs(FAmt_R(I) - FA_0) <= 1e-12*FA_0) then ; uBT_EE(I,j) = 0.0 ; else
          uBT_EE(I,j) = (1.5 * (duR(I) - du0(I,j))) * &
                                ((FAmt_R(I) - FA_avg) / (FAmt_R(I) - FA_0))
        endif
      else
        FA_u_W0(I,j) = 0.0 ; FA_u_WW(I,j) = 0.0
        FA_u_E0(I,j) = 0.0 ; FA_u_EE(I,j) = 0.0
        uBT_WW(I,j) = 0.0 ; uBT_EE(I,j) = 0.0
      endif
    enddo
  enddo

  !$omp target exit data map(release: duL, duR, du_CFL, FAmt_L, FAmT_R, FAmt_0, uhtot_L, uhtot_R)

end subroutine set_zonal_BT_cont_fortran

!> Shim for set_zonal_BT_cont -- dispatches via SET_ZONAL_BT_CONT_MODE env var.
subroutine set_zonal_BT_cont(bxC, u_a, h_in_a, h_W_a, h_E_a, BT_cont, du0_a, uh_tot_0_a, &
                             duhdu_tot_0_a, du_max_CFL_a, du_min_CFL_a, dt, &
                             dxCu_a, dy_Cu_a, IareaT_a, IdxT_a, CS, &
                             visc_rem_a, visc_rem_max_a, do_I_a, por_face_areaU_a)
  type(box_t),             intent(in) :: bxC  !< Iteration box for continuity solver
  type(RealArray_t),       intent(in) :: u_a   !< Zonal velocity [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in) :: h_in_a !< Layer thickness used to calculate
                                                !! fluxes [H ~> m or kg m-2].
  type(RealArray_t),       intent(in) :: h_W_a !< West edge thickness in the reconstruction
                                                !! [H ~> m or kg m-2].
  type(RealArray_t),       intent(in) :: h_E_a !< East edge thickness in the reconstruction
                                                !! [H ~> m or kg m-2].
  type(BT_cont_type),   intent(inout) :: BT_cont !< A structure with elements
                       !! that describe the effective open face areas as a function of barotropic flow.
  type(RealArray_t),       intent(in) :: du0_a  !< The barotropic velocity increment that gives 0
                                                 !! transport [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in) :: uh_tot_0_a    !< The summed transport with 0 adjustment
                                                        !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(RealArray_t),       intent(in) :: duhdu_tot_0_a !< The partial derivative
                       !! of du_err with du at 0 adjustment [H L ~> m2 or kg m-1].
  type(RealArray_t),       intent(in) :: du_max_CFL_a  !< Maximum acceptable value of
                                                        !! du [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in) :: du_min_CFL_a  !< Minimum acceptable value of
                                                        !! du [L T-1 ~> m s-1].
  real,                    intent(in) :: dt   !< Time increment [T ~> s].
  type(RealArray_t),       intent(in) :: dxCu_a !< The grid cell's u-point x-extent [L ~> m].
  type(RealArray_t),       intent(in) :: dy_Cu_a !< The grid cell's unblocked lengths of the
                                                !! u-faces of the h-cell [L ~> m].
  type(RealArray_t),       intent(in) :: IareaT_a !< The grid cell's 1/areaT [L-2 ~> m-2].
  type(RealArray_t),       intent(in) :: IdxT_a   !< The grid cell's 1/dxT [L-1 ~> m-1].
  type(transport_adjust_CS), intent(in) :: CS !< Options controlling the
                       !! transport adjustment and barotropic-consistency iteration.
  type(RealArray_t),       intent(in) :: visc_rem_a !< Both the fraction of the
                       !! momentum originally in a layer that remains after a time-step of viscosity, and
                       !! the fraction of a time-step's worth of a barotropic acceleration that a layer
                       !! experiences after viscosity is applied [nondim].
                       !! Visc_rem is between 0 (at the bottom) and 1 (far above the bottom).
  type(RealArray_t),       intent(in) :: visc_rem_max_a !< Maximum allowable visc_rem [nondim].
  type(LogicalArray_t),    intent(in) :: do_I_a   !< A logical flag indicating
                                                        !! which I values to work on.
  type(RealArray_t),       intent(in) :: por_face_areaU_a !< fractional open area
                                                        !! of U-faces [nondim]

  integer :: mode, rc
  type(RealArray_C) :: u_c, h_in_c, h_W_c, h_E_c
  type(RealArray_C) :: FA_u_W0_c, FA_u_E0_c, FA_u_WW_c, FA_u_EE_c, uBT_WW_c, uBT_EE_c
  type(RealArray_C) :: du0_c, uh_tot_0_c, duhdu_tot_0_c, du_max_CFL_c, du_min_CFL_c
  type(RealArray_C) :: dxCu_c, dy_Cu_c, IareaT_c, IdxT_c, visc_rem_c, visc_rem_max_c
  type(RealArray_C) :: por_face_areaU_c
  type(LogicalArray_C)        :: do_I_c
  type(Box_C)                 :: bxC_c
  type(transport_adjust_CS_C) :: CS_c
  type(io_recorder)  :: rec
  logical            :: capture
  character(len=80)  :: kernel
  character(len=100) :: dir
  character(len=256) :: binFile, metaFile

  kernel = "set_zonal_bt_cont"

  mode = getenv_mode("SET_ZONAL_BT_CONT_MODE", default=TIMH_runFORTRAN)

  select case (mode)

    case (TIMH_capture)
      capture = (.not. already_recorded(trim(kernel))) .and. is_root_pe()
      if (capture) then
        dir = "capture"
        rc = mkdir_posix(trim(dir) // c_null_char, int(o'755', c_int))
        binFile  = trim(dir) // "/" // trim(kernel) // ".bin"
        metaFile = trim(dir) // "/" // trim(kernel) // ".meta"
        call rec%open_write(binFile, metaFile)
        call rec%add("_bxC",            bxC)
        call rec%add("_u",              u_a)
        call rec%add("_h_in",           h_in_a)
        call rec%add("_h_W",            h_W_a)
        call rec%add("_h_E",            h_E_a)
        call rec%add("_FA_u_W0_before", BT_cont%FA_u_W0)
        call rec%add("_FA_u_E0_before", BT_cont%FA_u_E0)
        call rec%add("_FA_u_WW_before", BT_cont%FA_u_WW)
        call rec%add("_FA_u_EE_before", BT_cont%FA_u_EE)
        call rec%add("_uBT_WW_before",  BT_cont%uBT_WW)
        call rec%add("_uBT_EE_before",  BT_cont%uBT_EE)
        call rec%add("_du0",            du0_a)
        call rec%add("_uh_tot_0",       uh_tot_0_a)
        call rec%add("_duhdu_tot_0",    duhdu_tot_0_a)
        call rec%add("_du_max_CFL",     du_max_CFL_a)
        call rec%add("_du_min_CFL",     du_min_CFL_a)
        call rec%add("_dt",             dt)
        call rec%add("_dxCu",           dxCu_a)
        call rec%add("_dy_Cu",          dy_Cu_a)
        call rec%add("_IareaT",         IareaT_a)
        call rec%add("_IdxT",           IdxT_a)
        call rec%add("_vol_CFL",        CS%vol_CFL)
        call rec%add("_visc_rem",       visc_rem_a)
        call rec%add("_visc_rem_max",   visc_rem_max_a)
        call rec%add("_do_I",           do_I_a)
        call rec%add("_por_face_areaU", por_face_areaU_a)
      endif

      call set_zonal_BT_cont_fortran(bxC, u_a, h_in_a, h_W_a, h_E_a, BT_cont, du0_a, uh_tot_0_a, &
                             duhdu_tot_0_a, du_max_CFL_a, du_min_CFL_a, dt, &
                             dxCu_a, dy_Cu_a, IareaT_a, IdxT_a, CS, &
                             visc_rem_a, visc_rem_max_a, do_I_a, por_face_areaU_a)

      if (capture) then
        call rec%add("_FA_u_W0_after", BT_cont%FA_u_W0)
        call rec%add("_FA_u_E0_after", BT_cont%FA_u_E0)
        call rec%add("_FA_u_WW_after", BT_cont%FA_u_WW)
        call rec%add("_FA_u_EE_after", BT_cont%FA_u_EE)
        call rec%add("_uBT_WW_after",  BT_cont%uBT_WW)
        call rec%add("_uBT_EE_after",  BT_cont%uBT_EE)
        call rec%close()
        call mark_recorded(trim(kernel))
      endif

#ifdef _TIM
    case (TIMH_runAMREX)
      bxC_c            = bxC%to_c()
      u_c              = u_a%to_c()
      h_in_c           = h_in_a%to_c()
      h_W_c            = h_W_a%to_c()
      h_E_c            = h_E_a%to_c()
      FA_u_W0_c        = BT_cont%FA_u_W0%to_c()
      FA_u_E0_c        = BT_cont%FA_u_E0%to_c()
      FA_u_WW_c        = BT_cont%FA_u_WW%to_c()
      FA_u_EE_c        = BT_cont%FA_u_EE%to_c()
      uBT_WW_c         = BT_cont%uBT_WW%to_c()
      uBT_EE_c         = BT_cont%uBT_EE%to_c()
      du0_c            = du0_a%to_c()
      uh_tot_0_c       = uh_tot_0_a%to_c()
      duhdu_tot_0_c    = duhdu_tot_0_a%to_c()
      du_max_CFL_c     = du_max_CFL_a%to_c()
      du_min_CFL_c     = du_min_CFL_a%to_c()
      dxCu_c           = dxCu_a%to_c()
      dy_Cu_c          = dy_Cu_a%to_c()
      IareaT_c         = IareaT_a%to_c()
      IdxT_c           = IdxT_a%to_c()
      CS_c             = transport_adjust_CS_to_c(CS)
      visc_rem_c       = visc_rem_a%to_c()
      visc_rem_max_c   = visc_rem_max_a%to_c()
      do_I_c           = do_I_a%to_c()
      por_face_areaU_c = por_face_areaU_a%to_c()
      call turbotmp_set_zonal_bt_cont_bridge(bxC_c, u_c, h_in_c, h_W_c, h_E_c, &
                                             FA_u_W0_c, FA_u_E0_c, FA_u_WW_c, FA_u_EE_c, &
                                             uBT_WW_c, uBT_EE_c, du0_c, uh_tot_0_c, &
                                             duhdu_tot_0_c, du_max_CFL_c, du_min_CFL_c, dt, &
                                             dxCu_c, dy_Cu_c, IareaT_c, IdxT_c, CS_c, &
                                             visc_rem_c, visc_rem_max_c, do_I_c, por_face_areaU_c)
#endif

    case default
      call set_zonal_BT_cont_fortran(bxC, u_a, h_in_a, h_W_a, h_E_a, BT_cont, du0_a, uh_tot_0_a, &
                             duhdu_tot_0_a, du_max_CFL_a, du_min_CFL_a, dt, &
                             dxCu_a, dy_Cu_a, IareaT_a, IdxT_a, CS, &
                             visc_rem_a, visc_rem_max_a, do_I_a, por_face_areaU_a)

  end select

end subroutine set_zonal_BT_cont

!> Calculates the mass or volume fluxes through the meridional faces, and other related quantities.
subroutine meridional_mass_flux(bxC, v_a, h_in_a, h_S_a, h_N_a, vh_a, dt, &
                                dx_Cv_a, IareaT_a, IdyT_a, areaT_a, dyT_a, mask2dCv_a, dyCv_a, &
                                isd, ied, H_subroundoff, CS, &
                                OBC, por_face_areaV_a, vhbt_a, visc_rem_v_a, v_cor_a, BT_cont, &
                                dv_cor_a)
  type(Box_t),             intent(in)    :: bxC  !< Iteration box for continuity solver
  type(RealArray_t),  intent(in)  :: v_a    !< Meridional velocity [L T-1 ~> m s-1]
  type(RealArray_t),  intent(in)  :: h_in_a !< Layer thickness used to
                                                                  !! calculate fluxes [H ~> m or kg m-2]
  type(RealArray_t),  intent(in)  :: h_S_a !< South edge thickness in the
                                                                  !! reconstruction [H ~> m or kg m-2].
  type(RealArray_t),  intent(in)  :: h_N_a !< North edge thickness in the
                                                                  !! reconstruction [H ~> m or kg m-2].
  type(RealArray_t),  intent(inout) :: vh_a !< Volume flux through meridional
                                                                  !! faces = v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1]
  real,                                       intent(in)  :: dt   !< Time increment [T ~> s].
  type(RealArray_t),  intent(in)  :: dx_Cv_a !< The grid cell's unblocked lengths of the
                                                                  !! v-faces of the h-cell [L ~> m].
  type(RealArray_t),  intent(in)  :: IareaT_a !< The grid cell's 1/areaT [L-2 ~> m-2].
  type(RealArray_t),  intent(in)  :: IdyT_a   !< The grid cell's 1/dyT [L-1 ~> m-1].
  type(RealArray_t),  intent(in)  :: areaT_a  !< The area of the h-cell [L2 ~> m2].
  type(RealArray_t),  intent(in)  :: dyT_a    !< The y-extent of the h-cell [L ~> m].
  type(RealArray_t),  intent(in)  :: mask2dCv_a !< 0 for land points, 1 for ocean points at
                                                                  !! v-locations [nondim].
  type(RealArray_t),  intent(in)  :: dyCv_a   !< The grid cell's v-point y-extent [L ~> m].
  integer,                                    intent(in)  :: isd  !< The start i-index of
                                                                  !! the data domain.
  integer,                                    intent(in)  :: ied  !< The end i-index of
                                                                  !! the data domain.
  real,                                       intent(in)  :: H_subroundoff !< A negligibly small
                                                 !! thickness used to avoid division
                                                 !! by zero [H ~> m or kg m-2].
  type(transport_adjust_CS),                  intent(in)  :: CS !< Options
                       !! controlling the transport adjustment and barotropic-consistency iteration.
  type(ocean_OBC_type),                       pointer     :: OBC  !< Open boundary condition type
                                                                  !! specifies whether, where, and what
                                                                  !! open boundary conditions are used.
  type(RealArray_t), &
                 intent(in)  :: por_face_areaV_a !< fractional open area of V-faces [nondim]
  type(RealArray_t), intent(in) :: vhbt_a !< The summed volume flux through meridional
                                                                  !! faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(RealArray_t), &
                                              intent(in)  :: visc_rem_v_a !< Both the fraction of the momentum
                                   !! originally in a layer that remains after a time-step of viscosity,
                                   !! and the fraction of a time-step's worth of a barotropic acceleration
                                   !! that a layer experiences after viscosity is applied [nondim].
                                   !! Visc_rem_v is between 0 (at the bottom) and 1 (far above the bottom).
  type(RealArray_t), &
                                              intent(inout) :: v_cor_a
                                   !< The meridional velocities (v with a barotropic correction)
                                   !! that give vhbt as the depth-integrated transport [L T-1 ~> m s-1].
  type(BT_cont_type),               optional, pointer     :: BT_cont !< A structure with elements that describe
                                   !! the effective open face areas as a function of barotropic flow.
  type(RealArray_t), &
                                              intent(inout)   :: dv_cor_a !< The meridional velocity increments from v
                                                                  !! that give vhbt as the depth-integrated
                                                                  !! transports [L T-1 ~> m s-1].

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
  logical :: local_specified_BC, local_Flather_OBC, local_open_BC, any_simple_OBC  ! OBC-related logicals
  logical, dimension(v_a%lb(1):v_a%ub(1), v_a%lb(2):v_a%ub(2)) :: &
    simple_OBC_pt  ! Indicates points in a row with specified transport OBCs
  type(OBC_segment_type), pointer :: segment => NULL()
  real :: FAvi, FA_v    ! A sum of meridional face areas [H L ~> m2 or kg m-1].
  real, dimension(:,:), contiguous, pointer :: FA_v_S0, FA_v_N0, FA_v_SS, FA_v_NN, vBT_SS, vBT_NN
  real, dimension(:,:,:), contiguous, pointer :: visc_rem_v, v_cor
  real, dimension(:,:), contiguous, pointer :: dv_cor
  real, dimension(:,:,:), contiguous, pointer :: v, h_in, h_S, h_N, vh, por_face_areaV
  real, dimension(:,:), contiguous, pointer :: dx_Cv, IareaT, IdyT, areaT, dyT, mask2dCv, dyCv
  real, dimension(:,:,:), contiguous, pointer :: visc_rem_v_tmp
  real, dimension(:,:), contiguous, pointer :: dv_max_CFL, dv_min_CFL, dvhdv_tot_0, vh_tot_0
  real, dimension(:,:), contiguous, pointer :: visc_rem_max, dv
  logical, dimension(:,:), contiguous, pointer :: do_I
  type(RealArray_t) :: visc_rem_v_tmp_a
  type(RealArray_t) :: vh_tot_0_a, dvhdv_tot_0_a, dv_max_CFL_a, dv_min_CFL_a, dv_a
  type(RealArray_t) :: visc_rem_max_a
  type(LogicalArray_t) :: do_I_a
  ! Never allocated -- the zero-transport (dv0) correction never uses vhbt or reports vh_3d.
  type(RealArray_t) :: vhbt_none, vh_3d_none

  call cpu_clock_begin(id_clock_correct)

  call v_a%view(v)
  call h_in_a%view(h_in)
  call h_S_a%view(h_S)
  call h_N_a%view(h_N)
  call vh_a%view(vh)
  call por_face_areaV_a%view(por_face_areaV)
  call dx_Cv_a%view(dx_Cv)
  call IareaT_a%view(IareaT)
  call IdyT_a%view(IdyT)
  call areaT_a%view(areaT)
  call dyT_a%view(dyT)
  call mask2dCv_a%view(mask2dCv)
  call dyCv_a%view(dyCv)

  nullify(FA_v_S0, FA_v_N0, FA_v_SS, FA_v_NN, vBT_SS, vBT_NN)
  nullify(visc_rem_v, v_cor, dv_cor)
  if (visc_rem_v_a%associated()) call visc_rem_v_a%view(visc_rem_v)
  if (v_cor_a%associated()) call v_cor_a%view(v_cor)
  if (dv_cor_a%associated()) call dv_cor_a%view(dv_cor)

  use_visc_rem = visc_rem_v_a%associated()

  set_BT_cont = .false. ; if (present(BT_cont)) set_BT_cont = (associated(BT_cont))
  if (set_BT_cont) then
    call BT_cont%FA_v_S0%view(FA_v_S0) ; call BT_cont%FA_v_N0%view(FA_v_N0)
    call BT_cont%FA_v_SS%view(FA_v_SS) ; call BT_cont%FA_v_NN%view(FA_v_NN)
    call BT_cont%vBT_SS%view(vBT_SS)   ; call BT_cont%vBT_NN%view(vBT_NN)
  endif

  local_specified_BC = .false. ; local_Flather_OBC = .false. ; local_open_BC = .false.
  if (associated(OBC)) then ; if (OBC%OBC_pe) then
    local_specified_BC = OBC%specified_v_BCs_exist_globally
    local_Flather_OBC = OBC%Flather_v_BCs_exist_globally
    local_open_BC = OBC%open_v_BCs_exist_globally
  endif ; endif

  ish = bxC%idxS(1) ; ieh = bxC%idxE(1) ; jsh = bxC%idxS(2) ; jeh = bxC%idxE(2) ; nz  = bxC%idxE(3)

  CFL_dt = CS%CFL_limit_adjust / dt
  I_dt = 1.0 / dt
  if (CS%aggress_adjust) CFL_dt = I_dt

  call visc_rem_v_tmp_a%allocView(visc_rem_v_tmp, lb=[v_a%lb(1),v_a%lb(2),v_a%lb(3)], &
                                  ub=[v_a%ub(1),v_a%ub(2),v_a%ub(3)])
  if (vhbt_a%associated() .or. set_BT_cont) then
    call dv_max_CFL_a%allocView(dv_max_CFL, lb=[v_a%lb(1),v_a%lb(2)], ub=[v_a%ub(1),v_a%ub(2)])
    call dv_min_CFL_a%allocView(dv_min_CFL, lb=[v_a%lb(1),v_a%lb(2)], ub=[v_a%ub(1),v_a%ub(2)])
    call dvhdv_tot_0_a%allocView(dvhdv_tot_0, lb=[v_a%lb(1),v_a%lb(2)], ub=[v_a%ub(1),v_a%ub(2)])
    call vh_tot_0_a%allocView(vh_tot_0, lb=[v_a%lb(1),v_a%lb(2)], ub=[v_a%ub(1),v_a%ub(2)])
    call visc_rem_max_a%allocView(visc_rem_max, lb=[v_a%lb(1),v_a%lb(2)], ub=[v_a%ub(1),v_a%ub(2)])
  endif

  do concurrent (J=jsh-1:jeh)

    if (dv_cor_a%associated()) then
      do concurrent (i=ish:ieh)
        dv_cor(i,J) = 0.0
      enddo
    endif

    ! this is expensive
    if (.not.use_visc_rem) then
      do concurrent (k=1:nz, i=isd:ied)
        visc_rem_v_tmp(i,J,k) = 1.0
      enddo
    else
      do concurrent (k=1:nz, i=isd:ied)
        visc_rem_v_tmp(i,J,k) = visc_rem_v(i,J,k)
      enddo
    endif
    do concurrent (k=1:nz, i=ish:ieh)
      call flux_elem(v(i,J,k), h_in(i,J,k), h_in(i,J+1,k), h_S(i,J,k), h_S(i,J+1,k), h_N(i,J,k), &
                     h_N(i,J+1,k), vh(i,J,k), dvhdv(i,J,k), visc_rem_v_tmp(i,J,k), dx_Cv(i,J), &
                     IareaT(i,J), IareaT(i,J+1), IdyT(i,J), IdyT(i,J+1), dt, &
                     CS%vol_CFL, por_face_areaV(i,J,k))
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

    if (vhbt_a%associated() .or. set_BT_cont) then
      if (use_visc_rem .and. CS%use_visc_rem_max) then
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
        if (CS%vol_CFL) then
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
        if (CS%aggress_adjust) then
          ! untested
          do k=1,nz ; do concurrent (i=ish:ieh)
            if (CS%vol_CFL) then
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
            if (CS%vol_CFL) then
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
        if (CS%aggress_adjust) then
          ! untested
          do k=1,nz ; do concurrent (i=ish:ieh)
            if (CS%vol_CFL) then
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
            if (CS%vol_CFL) then
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
    endif ! vhbt_a%associated() .or. set_BT_cont

  enddo

  if (vhbt_a%associated() .or. set_BT_cont) then
    call do_I_a%allocView(do_I, lb=[v_a%lb(1),v_a%lb(2)], ub=[v_a%ub(1),v_a%ub(2)])
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

    call dv_a%allocView(dv, lb=[v_a%lb(1),v_a%lb(2)], ub=[v_a%ub(1),v_a%ub(2)], source=0.0)

    if (vhbt_a%associated()) then
      ! Find dv and vh.
      call meridional_flux_adjust(bxC, v_a, h_in_a, h_S_a, h_N_a, vh_tot_0_a, dvhdv_tot_0_a, dv_a, &
                             dv_max_CFL_a, dv_min_CFL_a, dt, dx_Cv_a, IareaT_a, IdyT_a, CS, &
                             visc_rem_v_tmp_a, &
                             do_I_a, por_face_areaV_a, vhbt_a=vhbt_a, vh_3d_a=vh_a, OBC=OBC)

      do concurrent (J=jsh-1:jeh)
        if (v_cor_a%associated()) then
          do concurrent (k=1:nz, i=ish:ieh)
            v_cor(i,J,k) = v(i,J,k) + dv(i,J) * visc_rem_v_tmp(i,J,k)
          enddo
          if (any_simple_OBC) then
            ! untested
            do concurrent (k=1:nz, i=ish:ieh, simple_OBC_pt(i,J))
              v_cor(i,J,k) = OBC%segment(abs(OBC%segnum_v(i,J)))%normal_vel(i,J,k)
            enddo
          endif
        endif ! v-corrected

        if (dv_cor_a%associated()) then
          do concurrent (i=ish:ieh)
            dv_cor(i,J) = dv(i,J)
          enddo
        endif ! dv-corrected
      enddo
    endif

    if (set_BT_cont) then
    ! Diagnose the zero-transport correction, dv0.
      call meridional_flux_adjust(bxC, v_a, h_in_a, h_S_a, h_N_a, vh_tot_0_a, dvhdv_tot_0_a, dv_a, &
                            dv_max_CFL_a, dv_min_CFL_a, dt, dx_Cv_a, IareaT_a, IdyT_a, CS, &
                            visc_rem_v_tmp_a, &
                            do_I_a, por_face_areaV_a, vhbt_a=vhbt_none, vh_3d_a=vh_3d_none)
      call set_merid_BT_cont(bxC, v_a, h_in_a, h_S_a, h_N_a, BT_cont, dv_a, vh_tot_0_a, &
                             dvhdv_tot_0_a, dv_max_CFL_a, dv_min_CFL_a, dt, &
                             dyCv_a, dx_Cv_a, IareaT_a, IdyT_a, CS, &
                             visc_rem_v_tmp_a, visc_rem_max_a, do_I_a, por_face_areaV_a)

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
          FA_v_S0(i,J) = FAvi ; FA_v_N0(i,J) = FAvi
          FA_v_SS(i,J) = FAvi ; FA_v_NN(i,J) = FAvi
          vBT_SS(i,J) = 0.0 ; vBT_NN(i,J) = 0.0
        enddo
      endif ! any_simple_OBC
    endif ! set_BT_cont
    call dv_a%free()
    call vh_tot_0_a%free() ; call dvhdv_tot_0_a%free() ; call dv_max_CFL_a%free()
    call dv_min_CFL_a%free() ; call visc_rem_max_a%free() ; call do_I_a%free()
  endif ! vhbt_a%associated() or set_BT_cont

  ! untested - probably needs to be refactored to be performant on GPU
  if (local_open_BC .and. set_BT_cont) then
    do n = 1, OBC%number_of_segments
      if (OBC%segment(n)%open .and. OBC%segment(n)%is_N_or_S) then
        J = OBC%segment(n)%HI%JsdB
        if (OBC%segment(n)%direction == OBC_DIRECTION_N) then
          do concurrent (i = OBC%segment(n)%HI%Isd:OBC%segment(n)%HI%Ied)
            FA_v = 0.0
            do k=1,nz ; FA_v = FA_v + h_in(i,j,k)*(dx_Cv(i,J)*por_face_areaV(i,J,k)) ; enddo
            FA_v_S0(i,J) = FA_v ; FA_v_N0(i,J) = FA_v
            FA_v_SS(i,J) = FA_v ; FA_v_NN(i,J) = FA_v
            vBT_SS(i,J) = 0.0 ; vBT_NN(i,J) = 0.0
          enddo
        else
          do concurrent (i = OBC%segment(n)%HI%Isd:OBC%segment(n)%HI%Ied)
            FA_v = 0.0
            do k=1,nz ; FA_v = FA_v + h_in(i,j+1,k)*(dx_Cv(i,J)*por_face_areaV(i,J,k)) ; enddo
            FA_v_S0(i,J) = FA_v ; FA_v_N0(i,J) = FA_v
            FA_v_SS(i,J) = FA_v ; FA_v_NN(i,J) = FA_v
            vBT_SS(i,J) = 0.0 ; vBT_NN(i,J) = 0.0
          enddo
        endif
      endif
    enddo
  endif

  if (set_BT_cont) then ; if (BT_cont%h_v%associated()) then
    if (v_cor_a%associated()) then
      call meridional_flux_thickness(bxC, v_cor_a, h_in_a, h_S_a, h_N_a, BT_cont%h_v, dt, &
                                    dx_Cv_a, IareaT_a, IdyT_a, CS%vol_CFL, &
                                    CS%marginal_faces, OBC, por_face_areaV_a, &
                                    visc_rem_v_tmp_a)
    else
      call meridional_flux_thickness(bxC, v_a, h_in_a, h_S_a, h_N_a, BT_cont%h_v, dt, &
                                    dx_Cv_a, IareaT_a, IdyT_a, &
                                    CS%vol_CFL, &
                                    CS%marginal_faces, OBC, por_face_areaV_a, &
                                    visc_rem_v_tmp_a)
    endif
  endif ; endif

  call visc_rem_v_tmp_a%free()

  call cpu_clock_end(id_clock_correct)

end subroutine meridional_mass_flux


!> Calculates the vertically integrated mass or volume fluxes through the meridional faces.
!> Original Fortran implementation of meridional_BT_mass_flux (renamed). Takes containers.
subroutine meridional_BT_mass_flux_fortran(bxC, v_a, h_in_a, h_S_a, h_N_a, vhbt_a, dt, &
                                   dx_Cv_a, IareaT_a, IdyT_a, CS, &
                                   OBC, por_face_areaV_a)

  type(box_t),                                intent(in)  :: bxC  !< Iteration box for continuity solver
  type(RealArray_t),  intent(in)  :: v_a    !< Meridional velocity [L T-1 ~> m s-1]
  type(RealArray_t),  intent(in)  :: h_in_a !< Layer thickness used to
                                                                  !! calculate fluxes [H ~> m or kg m-2]
  type(RealArray_t),  intent(in)  :: h_S_a !< Southern edge thickness in the PPM
                                                                  !! reconstruction [H ~> m or kg m-2].
  type(RealArray_t),  intent(in)  :: h_N_a !< Northern edge thickness in the PPM
                                                                  !! reconstruction [H ~> m or kg m-2].
  type(RealArray_t),  intent(inout) :: vhbt_a !< The summed volume flux through
                                                 !! meridional faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  real,                                       intent(in)  :: dt   !< Time increment [T ~> s].
  type(RealArray_t),  intent(in)  :: dx_Cv_a  !< The grid cell's unblocked lengths of
                                              !! the v-faces of the h-cell [L ~> m].
  type(RealArray_t),  intent(in)  :: IareaT_a !< The grid cell's 1/areaT [L-2 ~> m-2].
  type(RealArray_t),  intent(in)  :: IdyT_a   !< The grid cell's 1/dyT [L-1 ~> m-1].
  type(transport_adjust_CS),                  intent(in)  :: CS   !< Options controlling the
                                                                  !! transport adjustment and
                                                                  !! barotropic-consistency
                                                                  !! iteration.
  type(ocean_OBC_type),                       pointer     :: OBC  !< Open boundary condition type
                                                                  !! specifies whether, where, and what
                                                                  !! open boundary conditions are used.
  type(RealArray_t), &
                 intent(in)  :: por_face_areaV_a !< fractional open area of V-faces [nondim]

  ! Local variables
  real :: vh(v_a%lb(1):v_a%ub(1), v_a%lb(2):v_a%ub(2), v_a%lb(3):v_a%ub(3))
     ! Volume flux through meridional faces = v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1]
  real :: dvhdv(v_a%lb(1):v_a%ub(1), v_a%lb(2):v_a%ub(2), v_a%lb(3):v_a%ub(3))
     ! Partial derivative of vh with v [H L ~> m2 or kg m-1].
  integer :: i, j, k, ish, ieh, jsh, jeh, nz, l_seg
  logical :: local_specified_BC
  logical :: OBC_in_row(v_a%lb(2):v_a%ub(2))
  real, dimension(:,:,:), contiguous, pointer :: v, h_in, h_S, h_N, por_face_areaV
  real, dimension(:,:), contiguous, pointer :: vhbt
  real, dimension(:,:), contiguous, pointer :: dx_Cv, IareaT, IdyT

  call v_a%view(v)
  call h_in_a%view(h_in)
  call h_S_a%view(h_S)
  call h_N_a%view(h_N)
  call por_face_areaV_a%view(por_face_areaV)
  call vhbt_a%view(vhbt)
  call dx_Cv_a%view(dx_Cv)
  call IareaT_a%view(IareaT)
  call IdyT_a%view(IdyT)

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
                   CS%vol_CFL, por_face_areaV(i,J,k))
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

end subroutine meridional_BT_mass_flux_fortran

!> Shim for meridional_BT_mass_flux -- dispatches via MERIDIONAL_BT_MASS_FLUX_MODE env var.
subroutine meridional_BT_mass_flux(bxC, v_a, h_in_a, h_S_a, h_N_a, vhbt_a, dt, &
                                   dx_Cv_a, IareaT_a, IdyT_a, CS, &
                                   OBC, por_face_areaV_a)

  type(box_t),                                intent(in)  :: bxC  !< Iteration box for continuity solver
  type(RealArray_t),  intent(in)  :: v_a    !< Meridional velocity [L T-1 ~> m s-1]
  type(RealArray_t),  intent(in)  :: h_in_a !< Layer thickness used to
                                                                  !! calculate fluxes [H ~> m or kg m-2]
  type(RealArray_t),  intent(in)  :: h_S_a !< Southern edge thickness in the PPM
                                                                  !! reconstruction [H ~> m or kg m-2].
  type(RealArray_t),  intent(in)  :: h_N_a !< Northern edge thickness in the PPM
                                                                  !! reconstruction [H ~> m or kg m-2].
  type(RealArray_t),  intent(inout) :: vhbt_a !< The summed volume flux through
                                                 !! meridional faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  real,                                       intent(in)  :: dt   !< Time increment [T ~> s].
  type(RealArray_t),  intent(in)  :: dx_Cv_a  !< The grid cell's unblocked lengths of
                                              !! the v-faces of the h-cell [L ~> m].
  type(RealArray_t),  intent(in)  :: IareaT_a !< The grid cell's 1/areaT [L-2 ~> m-2].
  type(RealArray_t),  intent(in)  :: IdyT_a   !< The grid cell's 1/dyT [L-1 ~> m-1].
  type(transport_adjust_CS),                  intent(in)  :: CS   !< Options controlling the
                                                                  !! transport adjustment and
                                                                  !! barotropic-consistency
                                                                  !! iteration.
  type(ocean_OBC_type),                       pointer     :: OBC  !< Open boundary condition type
                                                                  !! specifies whether, where, and what
                                                                  !! open boundary conditions are used.
  type(RealArray_t), &
                 intent(in)  :: por_face_areaV_a !< fractional open area of V-faces [nondim]

  integer :: mode, rc
  type(RealArray_C) :: v_c, h_in_c, h_S_c, h_N_c, vhbt_c, dx_Cv_c, IareaT_c, IdyT_c
  type(RealArray_C) :: por_face_areaV_c
  type(Box_C)                 :: bxC_c
  type(transport_adjust_CS_C) :: CS_c
  type(c_ptr)                 :: OBC_c
  type(io_recorder)  :: rec
  logical            :: capture
  character(len=80)  :: kernel
  character(len=100) :: dir
  character(len=256) :: binFile, metaFile

  kernel = "meridional_bt_mass_flux"

  call cpu_clock_begin(id_clock_correct)

  mode = getenv_mode("MERIDIONAL_BT_MASS_FLUX_MODE", default=TIMH_runFORTRAN)

  select case (mode)

    case (TIMH_capture)
      capture = (.not. already_recorded(trim(kernel))) .and. is_root_pe()
      if (capture) then
        dir = "capture"
        rc = mkdir_posix(trim(dir) // c_null_char, int(o'755', c_int))
        binFile  = trim(dir) // "/" // trim(kernel) // ".bin"
        metaFile = trim(dir) // "/" // trim(kernel) // ".meta"
        call rec%open_write(binFile, metaFile)
        call rec%add("_bxC",            bxC)
        call rec%add("_v",              v_a)
        call rec%add("_h_in",           h_in_a)
        call rec%add("_h_S",            h_S_a)
        call rec%add("_h_N",            h_N_a)
        call rec%add("_vhbt_before",    vhbt_a)
        call rec%add("_dt",             dt)
        call rec%add("_dx_Cv",          dx_Cv_a)
        call rec%add("_IareaT",         IareaT_a)
        call rec%add("_IdyT",           IdyT_a)
        call rec%add("_vol_CFL",        CS%vol_CFL)
        call rec%add("_por_face_areaV", por_face_areaV_a)
      endif

      call meridional_BT_mass_flux_fortran(bxC, v_a, h_in_a, h_S_a, h_N_a, vhbt_a, dt, &
                                           dx_Cv_a, IareaT_a, IdyT_a, CS, &
                                           OBC, por_face_areaV_a)

      if (capture) then
        call rec%add("_vhbt_after", vhbt_a)
        call rec%close()
        call mark_recorded(trim(kernel))
      endif

#ifdef _TIM
    case (TIMH_runAMREX)
      bxC_c           = bxC%to_c()
      v_c             = v_a%to_c()
      h_in_c          = h_in_a%to_c()
      h_S_c           = h_S_a%to_c()
      h_N_c           = h_N_a%to_c()
      vhbt_c          = vhbt_a%to_c()
      dx_Cv_c         = dx_Cv_a%to_c()
      IareaT_c        = IareaT_a%to_c()
      IdyT_c          = IdyT_a%to_c()
      CS_c            = transport_adjust_CS_to_c(CS)
      por_face_areaV_c = por_face_areaV_a%to_c()
      if (associated(OBC)) then
        OBC_c = c_loc(OBC)
      else
        OBC_c = c_null_ptr
      endif
      call turbotmp_meridional_bt_mass_flux_bridge(bxC_c, v_c, h_in_c, h_S_c, h_N_c, vhbt_c, dt, &
                                                   dx_Cv_c, IareaT_c, IdyT_c, CS_c, OBC_c, &
                                                   por_face_areaV_c)
#endif

    case default
      call meridional_BT_mass_flux_fortran(bxC, v_a, h_in_a, h_S_a, h_N_a, vhbt_a, dt, &
                                           dx_Cv_a, IareaT_a, IdyT_a, CS, &
                                           OBC, por_face_areaV_a)

  end select

  call cpu_clock_end(id_clock_correct)

end subroutine meridional_BT_mass_flux


!> Sets the effective interface thickness associated with the fluxes at each meridional velocity point,
!! optionally scaling back these thicknesses to account for viscosity and fractional open areas.
!> Original Fortran implementation of meridional_flux_thickness (renamed). Takes containers.
subroutine meridional_flux_thickness_fortran(bxC, v_a, h_a, h_S_a, h_N_a, h_v_a, dt, &
                                     dx_Cv_a, IareaT_a, IdyT_a, vol_CFL, &
                                     marginal, OBC, por_face_areaV_a, visc_rem_v_a)
  type(box_t),                               intent(in)    :: bxC  !< Iteration box for continuity solver
  type(RealArray_t),  intent(in)  :: v_a   !< Meridional velocity [L T-1 ~> m s-1].
  type(RealArray_t),  intent(in)  :: h_a  !< Layer thickness used to
                                          !! calculate fluxes, [H ~> m or kg m-2].
  type(RealArray_t),  intent(in)  :: h_S_a !< South edge thickness in the
                                          !! reconstruction, [H ~> m or kg m-2].
  type(RealArray_t),  intent(in)  :: h_N_a !< North edge thickness in the
                                          !! reconstruction, [H ~> m or kg m-2].
  type(RealArray_t),  intent(inout):: h_v_a !< Effective thickness at meridional faces,
                                                                   !! scaled down to account for the effects of
                                                                   !! viscosity and the fractional open area
                                                                   !! [H ~> m or kg m-2].
  real,                                      intent(in)    :: dt   !< Time increment [T ~> s].
  type(RealArray_t),  intent(in)  :: dx_Cv_a  !< The grid cell's unblocked lengths of
                                              !! the v-faces of the h-cell [L ~> m].
  type(RealArray_t),  intent(in)  :: IareaT_a !< The grid cell's 1/areaT [L-2 ~> m-2].
  type(RealArray_t),  intent(in)  :: IdyT_a   !< The grid cell's 1/dyT [L-1 ~> m-1].
  logical,                                   intent(in)    :: vol_CFL !< If true, rescale the ratio
                          !! of face areas to the cell areas when estimating the CFL number.
  logical,                                   intent(in)    :: marginal !< If true, report the marginal
                          !! face thicknesses; otherwise report transport-averaged thicknesses.
  type(ocean_OBC_type),                      pointer       :: OBC !< Open boundaries control structure.
  type(RealArray_t), &
                                     intent(in) :: por_face_areaV_a  !< fractional open area of
                                                                     !! V-faces [nondim]
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
  real, dimension(:,:,:), contiguous, pointer :: visc_rem_v
  real, dimension(:,:,:), contiguous, pointer :: v, h, h_S, h_N, h_v, por_face_areaV
  real, dimension(:,:), contiguous, pointer :: dx_Cv, IareaT, IdyT
  type(box_t) :: bxV

  nullify(visc_rem_v)
  if (visc_rem_v_a%associated()) call visc_rem_v_a%view(visc_rem_v)
  call v_a%view(v)
  call h_a%view(h)
  call h_S_a%view(h_S)
  call h_N_a%view(h_N)
  call h_v_a%view(h_v)
  call por_face_areaV_a%view(por_face_areaV)
  call dx_Cv_a%view(dx_Cv)
  call IareaT_a%view(IareaT)
  call IdyT_a%view(IdyT)

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

end subroutine meridional_flux_thickness_fortran

!> Shim for meridional_flux_thickness -- dispatches via MERIDIONAL_FLUX_THICKNESS_MODE env var.
subroutine meridional_flux_thickness(bxC, v_a, h_a, h_S_a, h_N_a, h_v_a, dt, &
                                     dx_Cv_a, IareaT_a, IdyT_a, vol_CFL, &
                                     marginal, OBC, por_face_areaV_a, visc_rem_v_a)
  type(box_t),                               intent(in)    :: bxC  !< Iteration box for continuity solver
  type(RealArray_t),  intent(in)  :: v_a   !< Meridional velocity [L T-1 ~> m s-1].
  type(RealArray_t),  intent(in)  :: h_a  !< Layer thickness used to
                                          !! calculate fluxes, [H ~> m or kg m-2].
  type(RealArray_t),  intent(in)  :: h_S_a !< South edge thickness in the
                                          !! reconstruction, [H ~> m or kg m-2].
  type(RealArray_t),  intent(in)  :: h_N_a !< North edge thickness in the
                                          !! reconstruction, [H ~> m or kg m-2].
  type(RealArray_t),  intent(inout):: h_v_a !< Effective thickness at meridional faces,
                                                                   !! scaled down to account for the effects of
                                                                   !! viscosity and the fractional open area
                                                                   !! [H ~> m or kg m-2].
  real,                                      intent(in)    :: dt   !< Time increment [T ~> s].
  type(RealArray_t),  intent(in)  :: dx_Cv_a  !< The grid cell's unblocked lengths of
                                              !! the v-faces of the h-cell [L ~> m].
  type(RealArray_t),  intent(in)  :: IareaT_a !< The grid cell's 1/areaT [L-2 ~> m-2].
  type(RealArray_t),  intent(in)  :: IdyT_a   !< The grid cell's 1/dyT [L-1 ~> m-1].
  logical,                                   intent(in)    :: vol_CFL !< If true, rescale the ratio
                          !! of face areas to the cell areas when estimating the CFL number.
  logical,                                   intent(in)    :: marginal !< If true, report the marginal
                          !! face thicknesses; otherwise report transport-averaged thicknesses.
  type(ocean_OBC_type),                      pointer       :: OBC !< Open boundaries control structure.
  type(RealArray_t), &
                                     intent(in) :: por_face_areaV_a  !< fractional open area of
                                                                     !! V-faces [nondim]
  type(RealArray_t), intent(in) :: visc_rem_v_a !< Both the fraction
                          !! of the momentum originally in a layer that remains after a time-step of
                          !! viscosity, and the fraction of a time-step's worth of a barotropic
                          !! acceleration that a layer experiences after viscosity is applied [nondim].
                          !! Visc_rem_v is between 0 (at the bottom) and 1 (far above the bottom).

  integer :: mode, rc
  type(RealArray_C) :: v_c, h_c, h_S_c, h_N_c, h_v_c, dx_Cv_c, IareaT_c, IdyT_c
  type(RealArray_C) :: por_face_areaV_c, visc_rem_v_c
  type(Box_C)         :: bxC_c
  type(c_ptr)         :: OBC_c
  logical(c_bool)     :: vol_CFL_c, marginal_c
  type(io_recorder)   :: rec
  logical             :: capture
  character(len=80)   :: kernel
  character(len=100)  :: dir
  character(len=256)  :: binFile, metaFile

  kernel = "meridional_flux_thickness"

  mode = getenv_mode("MERIDIONAL_FLUX_THICKNESS_MODE", default=TIMH_runFORTRAN)

  select case (mode)

    case (TIMH_capture)
      capture = (.not. already_recorded(trim(kernel))) .and. is_root_pe()
      if (capture) then
        dir = "capture"
        rc = mkdir_posix(trim(dir) // c_null_char, int(o'755', c_int))
        binFile  = trim(dir) // "/" // trim(kernel) // ".bin"
        metaFile = trim(dir) // "/" // trim(kernel) // ".meta"
        call rec%open_write(binFile, metaFile)
        call rec%add("_bxC",            bxC)
        call rec%add("_v",              v_a)
        call rec%add("_h",              h_a)
        call rec%add("_h_S",            h_S_a)
        call rec%add("_h_N",            h_N_a)
        call rec%add("_h_v_before",     h_v_a)
        call rec%add("_dt",             dt)
        call rec%add("_dx_Cv",          dx_Cv_a)
        call rec%add("_IareaT",         IareaT_a)
        call rec%add("_IdyT",           IdyT_a)
        call rec%add("_vol_CFL",        vol_CFL)
        call rec%add("_marginal",       marginal)
        call rec%add("_por_face_areaV", por_face_areaV_a)
        if (visc_rem_v_a%associated()) call rec%add("_visc_rem_v", visc_rem_v_a)
      endif

      call meridional_flux_thickness_fortran(bxC, v_a, h_a, h_S_a, h_N_a, h_v_a, dt, &
                                             dx_Cv_a, IareaT_a, IdyT_a, vol_CFL, &
                                             marginal, OBC, por_face_areaV_a, visc_rem_v_a)

      if (capture) then
        call rec%add("_h_v_after", h_v_a)
        call rec%close()
        call mark_recorded(trim(kernel))
      endif

#ifdef _TIM
    case (TIMH_runAMREX)
      bxC_c            = bxC%to_c()
      v_c              = v_a%to_c()
      h_c              = h_a%to_c()
      h_S_c            = h_S_a%to_c()
      h_N_c            = h_N_a%to_c()
      h_v_c            = h_v_a%to_c()
      dx_Cv_c          = dx_Cv_a%to_c()
      IareaT_c         = IareaT_a%to_c()
      IdyT_c           = IdyT_a%to_c()
      vol_CFL_c        = vol_CFL
      marginal_c       = marginal
      por_face_areaV_c = por_face_areaV_a%to_c()
      visc_rem_v_c     = visc_rem_v_a%to_c()
      if (associated(OBC)) then
        OBC_c = c_loc(OBC)
      else
        OBC_c = c_null_ptr
      endif
      call turbotmp_meridional_flux_thickness_bridge(bxC_c, v_c, h_c, h_S_c, h_N_c, h_v_c, dt, &
                                                     dx_Cv_c, IareaT_c, IdyT_c, vol_CFL_c, &
                                                     marginal_c, OBC_c, por_face_areaV_c, &
                                                     visc_rem_v_c)
#endif

    case default
      call meridional_flux_thickness_fortran(bxC, v_a, h_a, h_S_a, h_N_a, h_v_a, dt, &
                                             dx_Cv_a, IareaT_a, IdyT_a, vol_CFL, &
                                             marginal, OBC, por_face_areaV_a, visc_rem_v_a)

  end select

end subroutine meridional_flux_thickness


!> Returns the barotropic velocity adjustment that gives the desired barotropic (layer-summed) transport.
!> Original Fortran implementation of meridional_flux_adjust (renamed). Takes containers.
subroutine meridional_flux_adjust_fortran(bxC, v_a, h_in_a, h_S_a, h_N_a, vh_tot_0_a, dvhdv_tot_0_a, &
                             dv_a, dv_max_CFL_a, dv_min_CFL_a, dt, dx_Cv_a, IareaT_a, IdyT_a, CS, &
                             visc_rem_a, do_I_in_a, por_face_areaV_a, vhbt_a, vh_3d_a, OBC)
  type(box_t),             intent(in)    :: bxC  !< Iteration box for continuity solver
  type(RealArray_t),       intent(in)    :: v_a   !< Meridional velocity [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in)    :: h_in_a !< Layer thickness used to calculate
                                                 !! fluxes [H ~> m or kg m-2].
  type(RealArray_t),       intent(in)    :: h_S_a !< South edge thickness in the
                                                 !! reconstruction [H ~> m or kg m-2].
  type(RealArray_t),       intent(in)    :: h_N_a !< North edge thickness in the
                                                 !! reconstruction [H ~> m or kg m-2].
  type(RealArray_t),       intent(in)    :: visc_rem_a
                             !< Both the fraction of the momentum originally
                             !! in a layer that remains after a time-step of viscosity, and the
                             !! fraction of a time-step's worth of a barotropic acceleration that
                             !! a layer experiences after viscosity is applied [nondim].
                             !! Visc_rem is between 0 (at the bottom) and 1 (far above the bottom).
  type(RealArray_t), &
                           intent(in)    :: vhbt_a !< The summed volume flux through
                                                 !! meridional faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(RealArray_t),       intent(in)    :: dv_max_CFL_a !< Maximum acceptable value of
                                                       !! dv [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in)    :: dv_min_CFL_a !< Minimum acceptable value of
                                                       !! dv [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in)    :: vh_tot_0_a   !< The summed transport with 0 adjustment
                                                       !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(RealArray_t),       intent(in)  :: dvhdv_tot_0_a !< The partial derivative of dv_err with
                                                      !! dv at 0 adjustment [H L ~> m2 or kg m-1].
  type(RealArray_t),       intent(inout) :: dv_a      !< The barotropic velocity
                                                       !! adjustment [L T-1 ~> m s-1].
  real,                    intent(in)  :: dt      !< Time increment [T ~> s].
  type(RealArray_t),       intent(in)  :: dx_Cv_a !< The grid cell's unblocked lengths of the
                                                 !! v-faces of the h-cell [L ~> m].
  type(RealArray_t),       intent(in)  :: IareaT_a !< The grid cell's 1/areaT [L-2 ~> m-2].
  type(RealArray_t),       intent(in)  :: IdyT_a   !< The grid cell's 1/dyT [L-1 ~> m-1].
  type(transport_adjust_CS), intent(in) :: CS !< Options controlling the
                       !! transport adjustment and barotropic-consistency iteration.
  type(LogicalArray_t),    intent(in)  :: do_I_in_a  !< A flag indicating which I values to work on.
  type(RealArray_t),       intent(in)  :: por_face_areaV_a !< fractional open area of
                                                       !! V-faces [nondim]
  type(RealArray_t), &
                           intent(inout) :: vh_3d_a !< Volume flux through meridional
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
  logical :: do_I(v_a%lb(1):v_a%ub(1))
  logical :: local_OBC, use_vhbt, use_vh_3d
  integer, parameter :: max_itts = 20
  integer :: ish     !< Start of i index range.
  integer :: ieh     !< End of i index range.
  integer :: jsh     !< Start of j index range.
  integer :: jeh     !< End of j index range.
  real, dimension(:,:), contiguous, pointer :: vhbt
  real, dimension(:,:,:), contiguous, pointer :: vh_3d
  real, dimension(:,:,:), contiguous, pointer :: v, h_in, h_S, h_N, visc_rem, por_face_areaV
  real, dimension(:,:), contiguous, pointer :: dv_max_CFL, dv_min_CFL, vh_tot_0, dvhdv_tot_0, dv
  real, dimension(:,:), contiguous, pointer :: dx_Cv, IareaT, IdyT
  logical, dimension(:,:), contiguous, pointer :: do_I_in

  nullify(vhbt, vh_3d)
  if (vhbt_a%associated()) call vhbt_a%view(vhbt)
  if (vh_3d_a%associated()) call vh_3d_a%view(vh_3d)

  call v_a%view(v)
  call h_in_a%view(h_in)
  call h_S_a%view(h_S)
  call h_N_a%view(h_N)
  call visc_rem_a%view(visc_rem)
  call dv_max_CFL_a%view(dv_max_CFL)
  call dv_min_CFL_a%view(dv_min_CFL)
  call vh_tot_0_a%view(vh_tot_0)
  call dvhdv_tot_0_a%view(dvhdv_tot_0)
  call dv_a%view(dv)
  call do_I_in_a%view(do_I_in)
  call por_face_areaV_a%view(por_face_areaV)
  call dx_Cv_a%view(dx_Cv)
  call IareaT_a%view(IareaT)
  call IdyT_a%view(IdyT)

  local_OBC = .false.
  if (present(OBC)) then
    if (associated(OBC)) then
      local_OBC = OBC%open_u_BCs_exist_globally
    endif
  endif

  use_vhbt = vhbt_a%associated()
  use_vh_3d = vh_3d_a%associated()

  ish = bxC%idxS(1) ; ieh = bxC%idxE(1) ; jsh = bxC%idxS(2) ; jeh = bxC%idxE(2) ; nz  = bxC%idxE(3)

  tol_vel = CS%tol_vel

  ! NVIDIA needs private arrays to be alloc'ed to prevent data transfers.
  ! GCC doesn't understand map(alloc: ...) for variables also marked private
  !$omp target enter data map(alloc: do_I, dv_max, dv_min, dvhdv_tot, vh_err, vh_err_best, vh_aux)

  !$omp target teams loop &
  !$omp   private(vh_err, vh_err_best, dvhdv_tot, dv_min, dv_max, do_I, vh_aux, itt, tol_eta)
  do J=jsh-1,jeh

    if (use_vh_3d) then
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
        case (:1) ; tol_eta = 1e-6 * CS%tol_eta
        case (2)  ; tol_eta = 1e-4 * CS%tol_eta
        case (3)  ; tol_eta = 1e-2 * CS%tol_eta
        case default ; tol_eta = CS%tol_eta
      end select

      do concurrent (i=ish:ieh)
        if (vh_err(i) > 0.0) then ; dv_max(i) = dv(i,j)
        elseif (vh_err(i) < 0.0) then ; dv_min(i) = dv(i,j)
        else ; do_I(i) = .false. ; endif
      enddo

      do concurrent (i=ish:ieh, do_I(i)) &
          & DO_LOCALITY(local(ddv, dv_prev))
        if ((dt * min(IareaT(i,j),IareaT(i,j+1))*abs(vh_err(i)) > tol_eta) .or. &
            (CS%better_iter .and. &
             ((abs(vh_err(i)) > tol_vel * dvhdv_tot(i)) .or. &
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

      if ((itt < max_itts) .or. use_vh_3d) then
        do concurrent (i=ish:ieh)
          vh_err(i) = 0.0 ; dvhdv_tot(i) = 0.0
          if (use_vhbt) vh_err(i) = -vhbt(i,J)
        enddo
        do k=1,nz ; do concurrent (i=ish:ieh, do_I(i)) DO_LOCALITY(local(v_new, dvhdv))
          v_new = v(i,J,k) + dv(i,j) * visc_rem(i,j,k)
          call flux_elem(v_new, h_in(i,J,k), h_in(i,J+1,k), h_S(i,J,k), h_S(i,J+1,k), &
                         h_N(i,J,k), h_N(i,J+1,k), vh_aux(i,k), dvhdv, visc_rem(i,J,k), &
                         dx_Cv(i,J), IareaT(i,J), IareaT(i,J+1), idyT(i,J), IdyT(i,J+1), &
                         dt, CS%vol_CFL, por_face_areaV(i,J,k))
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

    if (use_vh_3d) then
      do concurrent (k=1:nz, i=ish:ieh)
        vh_3d(i,J,k) = vh_aux(i,k)
      enddo
    endif
  enddo ! j-loop

  !$omp target exit data map(release: do_I, dv_max, dv_min, dvhdv_tot, vh_err, vh_err_best, vh_aux)

end subroutine meridional_flux_adjust_fortran

!> Shim for meridional_flux_adjust -- dispatches via MERIDIONAL_FLUX_ADJUST_MODE env var.
subroutine meridional_flux_adjust(bxC, v_a, h_in_a, h_S_a, h_N_a, vh_tot_0_a, dvhdv_tot_0_a, &
                             dv_a, dv_max_CFL_a, dv_min_CFL_a, dt, dx_Cv_a, IareaT_a, IdyT_a, CS, &
                             visc_rem_a, do_I_in_a, por_face_areaV_a, vhbt_a, vh_3d_a, OBC)
  type(box_t),             intent(in)    :: bxC  !< Iteration box for continuity solver
  type(RealArray_t),       intent(in)    :: v_a   !< Meridional velocity [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in)    :: h_in_a !< Layer thickness used to calculate
                                                 !! fluxes [H ~> m or kg m-2].
  type(RealArray_t),       intent(in)    :: h_S_a !< South edge thickness in the
                                                 !! reconstruction [H ~> m or kg m-2].
  type(RealArray_t),       intent(in)    :: h_N_a !< North edge thickness in the
                                                 !! reconstruction [H ~> m or kg m-2].
  type(RealArray_t),       intent(in)    :: visc_rem_a
                             !< Both the fraction of the momentum originally
                             !! in a layer that remains after a time-step of viscosity, and the
                             !! fraction of a time-step's worth of a barotropic acceleration that
                             !! a layer experiences after viscosity is applied [nondim].
                             !! Visc_rem is between 0 (at the bottom) and 1 (far above the bottom).
  type(RealArray_t), &
                           intent(in)    :: vhbt_a !< The summed volume flux through
                                                 !! meridional faces [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(RealArray_t),       intent(in)    :: dv_max_CFL_a !< Maximum acceptable value of
                                                       !! dv [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in)    :: dv_min_CFL_a !< Minimum acceptable value of
                                                       !! dv [L T-1 ~> m s-1].
  type(RealArray_t),       intent(in)    :: vh_tot_0_a   !< The summed transport with 0 adjustment
                                                       !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(RealArray_t),       intent(in)  :: dvhdv_tot_0_a !< The partial derivative of dv_err with
                                                      !! dv at 0 adjustment [H L ~> m2 or kg m-1].
  type(RealArray_t),       intent(inout) :: dv_a      !< The barotropic velocity
                                                       !! adjustment [L T-1 ~> m s-1].
  real,                    intent(in)  :: dt      !< Time increment [T ~> s].
  type(RealArray_t),       intent(in)  :: dx_Cv_a !< The grid cell's unblocked lengths of the
                                                 !! v-faces of the h-cell [L ~> m].
  type(RealArray_t),       intent(in)  :: IareaT_a !< The grid cell's 1/areaT [L-2 ~> m-2].
  type(RealArray_t),       intent(in)  :: IdyT_a   !< The grid cell's 1/dyT [L-1 ~> m-1].
  type(transport_adjust_CS), intent(in) :: CS !< Options controlling the
                       !! transport adjustment and barotropic-consistency iteration.
  type(LogicalArray_t),    intent(in)  :: do_I_in_a  !< A flag indicating which I values to work on.
  type(RealArray_t),       intent(in)  :: por_face_areaV_a !< fractional open area of
                                                       !! V-faces [nondim]
  type(RealArray_t), &
                           intent(inout) :: vh_3d_a !< Volume flux through meridional
                             !! faces = v*h*dx [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(ocean_OBC_type), optional, pointer :: OBC !< Open boundaries control structure.

  integer :: mode, rc
  type(RealArray_C) :: v_c, h_in_c, h_S_c, h_N_c, vh_tot_0_c, dvhdv_tot_0_c, dv_c
  type(RealArray_C) :: dv_max_CFL_c, dv_min_CFL_c, dx_Cv_c, IareaT_c, IdyT_c
  type(RealArray_C) :: visc_rem_c, por_face_areaV_c, vhbt_c, vh_3d_c
  type(LogicalArray_C)        :: do_I_in_c
  type(Box_C)                 :: bxC_c
  type(transport_adjust_CS_C) :: CS_c
  type(c_ptr)                 :: OBC_c
  type(io_recorder)  :: rec
  logical            :: capture
  character(len=80)  :: kernel
  character(len=100) :: dir
  character(len=256) :: binFile, metaFile

  kernel = "meridional_flux_adjust"

  mode = getenv_mode("MERIDIONAL_FLUX_ADJUST_MODE", default=TIMH_runFORTRAN)

  select case (mode)

    case (TIMH_capture)
      capture = (.not. already_recorded(trim(kernel))) .and. is_root_pe()
      if (capture) then
        dir = "capture"
        rc = mkdir_posix(trim(dir) // c_null_char, int(o'755', c_int))
        binFile  = trim(dir) // "/" // trim(kernel) // ".bin"
        metaFile = trim(dir) // "/" // trim(kernel) // ".meta"
        call rec%open_write(binFile, metaFile)
        call rec%add("_bxC",            bxC)
        call rec%add("_v",              v_a)
        call rec%add("_h_in",           h_in_a)
        call rec%add("_h_S",            h_S_a)
        call rec%add("_h_N",            h_N_a)
        call rec%add("_vh_tot_0",       vh_tot_0_a)
        call rec%add("_dvhdv_tot_0",    dvhdv_tot_0_a)
        call rec%add("_dv_before",      dv_a)
        call rec%add("_dv_max_CFL",     dv_max_CFL_a)
        call rec%add("_dv_min_CFL",     dv_min_CFL_a)
        call rec%add("_dt",             dt)
        call rec%add("_dx_Cv",          dx_Cv_a)
        call rec%add("_IareaT",         IareaT_a)
        call rec%add("_IdyT",           IdyT_a)
        call rec%add("_tol_eta",        CS%tol_eta)
        call rec%add("_tol_vel",        CS%tol_vel)
        call rec%add("_better_iter",    CS%better_iter)
        call rec%add("_vol_CFL",        CS%vol_CFL)
        call rec%add("_visc_rem",       visc_rem_a)
        call rec%add("_do_I_in",        do_I_in_a)
        call rec%add("_por_face_areaV", por_face_areaV_a)
        if (vhbt_a%associated()) call rec%add("_vhbt", vhbt_a)
        if (vh_3d_a%associated()) call rec%add("_vh_3d_before", vh_3d_a)
      endif

      call meridional_flux_adjust_fortran(bxC, v_a, h_in_a, h_S_a, h_N_a, vh_tot_0_a, dvhdv_tot_0_a, &
                           dv_a, dv_max_CFL_a, dv_min_CFL_a, dt, dx_Cv_a, IareaT_a, IdyT_a, CS, &
                           visc_rem_a, do_I_in_a, por_face_areaV_a, vhbt_a, vh_3d_a, OBC)

      if (capture) then
        call rec%add("_dv_after", dv_a)
        if (vh_3d_a%associated()) call rec%add("_vh_3d_after", vh_3d_a)
        call rec%close()
        call mark_recorded(trim(kernel))
      endif

#ifdef _TIM
    case (TIMH_runAMREX)
      bxC_c            = bxC%to_c()
      v_c              = v_a%to_c()
      h_in_c           = h_in_a%to_c()
      h_S_c            = h_S_a%to_c()
      h_N_c            = h_N_a%to_c()
      vh_tot_0_c       = vh_tot_0_a%to_c()
      dvhdv_tot_0_c    = dvhdv_tot_0_a%to_c()
      dv_c             = dv_a%to_c()
      dv_max_CFL_c     = dv_max_CFL_a%to_c()
      dv_min_CFL_c     = dv_min_CFL_a%to_c()
      dx_Cv_c          = dx_Cv_a%to_c()
      IareaT_c         = IareaT_a%to_c()
      IdyT_c           = IdyT_a%to_c()
      CS_c             = transport_adjust_CS_to_c(CS)
      visc_rem_c       = visc_rem_a%to_c()
      do_I_in_c        = do_I_in_a%to_c()
      por_face_areaV_c = por_face_areaV_a%to_c()
      vhbt_c           = vhbt_a%to_c()
      vh_3d_c          = vh_3d_a%to_c()
      if (present(OBC)) then
        if (associated(OBC)) then
          OBC_c = c_loc(OBC)
        else
          OBC_c = c_null_ptr
        endif
      else
        OBC_c = c_null_ptr
      endif
      call turbotmp_meridional_flux_adjust_bridge(bxC_c, v_c, h_in_c, h_S_c, h_N_c, vh_tot_0_c, &
                                                  dvhdv_tot_0_c, dv_c, dv_max_CFL_c, dv_min_CFL_c, &
                                                  dt, dx_Cv_c, IareaT_c, IdyT_c, CS_c, visc_rem_c, &
                                                  do_I_in_c, por_face_areaV_c, vhbt_c, vh_3d_c, OBC_c)
#endif

    case default
      call meridional_flux_adjust_fortran(bxC, v_a, h_in_a, h_S_a, h_N_a, vh_tot_0_a, dvhdv_tot_0_a, &
                           dv_a, dv_max_CFL_a, dv_min_CFL_a, dt, dx_Cv_a, IareaT_a, IdyT_a, CS, &
                           visc_rem_a, do_I_in_a, por_face_areaV_a, vhbt_a, vh_3d_a, OBC)

  end select

end subroutine meridional_flux_adjust


!> Sets of a structure that describes the meridional barotropic volume or mass fluxes as a
!! function of barotropic flow to agree closely with the sum of the layer's transports.
!> Original Fortran implementation of set_merid_BT_cont (renamed). Takes containers.
subroutine set_merid_BT_cont_fortran(bxC, v_a, h_in_a, h_S_a, h_N_a, BT_cont, dv0_a, vh_tot_0_a, &
                             dvhdv_tot_0_a, dv_max_CFL_a, dv_min_CFL_a, dt, &
                             dyCv_a, dx_Cv_a, IareaT_a, IdyT_a, CS, &
                             visc_rem_a, visc_rem_max_a, do_I_a, por_face_areaV_a)
  type(box_t),                                intent(in)    :: bxC  !< Iteration box for continuity solver
  type(RealArray_t),  intent(in)   :: v_a   !< Meridional velocity [L T-1 ~> m s-1].
  type(RealArray_t),  intent(in)   :: h_in_a !< Layer thickness used to
                                            !! calculate fluxes, [H ~> m or kg m-2].
  type(RealArray_t),  intent(in)   :: h_S_a !< South edge thickness in the
                                            !! reconstruction, [H ~> m or kg m-2].
  type(RealArray_t),  intent(in)   :: h_N_a !< North edge thickness in the
                                            !! reconstruction, [H ~> m or kg m-2].
  type(BT_cont_type),                         intent(inout) :: BT_cont !< A structure with elements
                       !! that describe the effective open face areas as a function of barotropic flow.
  type(RealArray_t),  intent(in)   :: dv0_a  !< The barotropic velocity increment
                                            !! that gives 0 transport [L T-1 ~> m s-1].
  type(RealArray_t),  intent(in)   :: vh_tot_0_a !< The summed transport
                       !! with 0 adjustment [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(RealArray_t),  intent(in)   :: dvhdv_tot_0_a !< The partial derivative
                       !! of du_err with dv at 0 adjustment [H L ~> m2 or kg m-1].
  type(RealArray_t),  intent(in)   :: dv_max_CFL_a !< Maximum acceptable value
                                                                          !!  of dv [L T-1 ~> m s-1].
  type(RealArray_t),  intent(in)   :: dv_min_CFL_a !< Minimum acceptable value
                                                                          !!  of dv [L T-1 ~> m s-1].
  real,                                       intent(in)    :: dt   !< Time increment [T ~> s].
  type(RealArray_t),  intent(in)   :: dyCv_a !< The grid cell's v-point y-extent [L ~> m].
  type(RealArray_t),  intent(in)   :: dx_Cv_a !< The grid cell's unblocked lengths of the
                                            !! v-faces of the h-cell [L ~> m].
  type(RealArray_t),  intent(in)   :: IareaT_a !< The grid cell's 1/areaT [L-2 ~> m-2].
  type(RealArray_t),  intent(in)   :: IdyT_a   !< The grid cell's 1/dyT [L-1 ~> m-1].
  type(transport_adjust_CS),           intent(in)    :: CS !< Options
                       !! controlling the transport adjustment and barotropic-consistency iteration.
  type(RealArray_t),  intent(in)   :: visc_rem_a !< Both the fraction of the
                       !! momentum originally in a layer that remains after a time-step
                       !! of viscosity, and the fraction of a time-step's worth of a barotropic
                       !! acceleration that a layer experiences after viscosity is applied [nondim].
                       !! Visc_rem is between 0 (at the bottom) and 1 (far above the bottom).
  type(RealArray_t), &
                 intent(in)    :: visc_rem_max_a !< Maximum allowable visc_rem [nondim]
  type(LogicalArray_t),  intent(in)  :: do_I_a !< A logical flag indicating
                                             !! which I values to work on.
  type(RealArray_t),  intent(in)   :: por_face_areaV_a !< fractional open
                                            !! area of V-faces [nondim]
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
  real, dimension(:,:), contiguous, pointer :: dv0, visc_rem_max
  real, dimension(:,:), contiguous, pointer :: dyCv, dx_Cv, IareaT, IdyT
  logical, dimension(:,:), contiguous, pointer :: do_I
  real, dimension(:,:), contiguous, pointer :: FA_v_S0, FA_v_N0, FA_v_SS, FA_v_NN, vBT_SS, vBT_NN

  call v_a%view(v)
  call h_in_a%view(h_in)
  call h_S_a%view(h_S)
  call h_N_a%view(h_N)
  call dv0_a%view(dv0)
  call visc_rem_a%view(visc_rem)
  call visc_rem_max_a%view(visc_rem_max)
  call do_I_a%view(do_I)
  call por_face_areaV_a%view(por_face_areaV)
  call dyCv_a%view(dyCv)
  call dx_Cv_a%view(dx_Cv)
  call IareaT_a%view(IareaT)
  call IdyT_a%view(IdyT)

  ish = bxC%idxS(1) ; ieh = bxC%idxE(1) ; jsh = bxC%idxS(2) ; jeh = bxC%idxE(2) ; nz  = bxC%idxE(3)
  Idt = 1.0 / dt
  min_visc_rem = 0.1 ; CFL_min = 1e-6

  call BT_cont%FA_v_S0%view(FA_v_S0) ; call BT_cont%FA_v_N0%view(FA_v_N0)
  call BT_cont%FA_v_SS%view(FA_v_SS) ; call BT_cont%FA_v_NN%view(FA_v_NN)
  call BT_cont%vBT_SS%view(vBT_SS)   ; call BT_cont%vBT_NN%view(vBT_NN)

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
                     IdyT(i,J+1), dt, CS%vol_CFL, &
                     por_face_areaV(i,J,k))
      call flux_elem(v_L, h_in(i,J,k), h_in(i,J+1,k), h_S(i,J,k), h_S(i,J+1,k), &
                     h_N(i,J,k), h_N(i,J+1,k), vh_L, dvhdv_L, visc_rem(i,J,k), &
                     dx_Cv(i,J), IareaT(i,J), IareaT(i,J+1), IdyT(i,J), &
                     IdyT(i,J+1), dt, CS%vol_CFL, &
                     por_face_areaV(i,J,k))
      call flux_elem(v_R, h_in(i,J,k), h_in(i,J+1,k), h_S(i,J,k), h_S(i,J+1,k), &
                     h_N(i,J,k), h_N(i,J+1,k), vh_R, dvhdv_R, visc_rem(i,J,k), &
                     dx_Cv(i,J), IareaT(i,J), IareaT(i,J+1), IdyT(i,J), &
                     IdyT(i,J+1), dt, CS%vol_CFL, &
                     por_face_areaV(i,J,k))
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
        FA_v_S0(i,J) = FA_0 ; FA_v_SS(i,J) = FAmt_L(i)
        if (abs(FA_0-FAmt_L(i)) <= 1e-12*FA_0) then ; vBT_SS(i,J) = 0.0 ; else
          vBT_SS(i,J) = (1.5 * (dvL(i) - dv0(i,J))) * &
                      ((FAmt_L(i) - FA_avg) / (FAmt_L(i) - FA_0))
        endif

        FA_0 = FAmt_0(i) ; FA_avg = FAmt_0(i)
        if ((dvR(i) - dv0(i,j)) /= 0.0) &
          FA_avg = vhtot_R(i) / (dvR(i) - dv0(i,j))
        if (FA_avg > max(FA_0, FAmt_R(i))) then ; FA_avg = max(FA_0, FAmt_R(i))
        elseif (FA_avg < min(FA_0, FAmt_R(i))) then ; FA_0 = FA_avg ; endif
        FA_v_N0(i,J) = FA_0 ; FA_v_NN(i,J) = FAmt_R(i)
        if (abs(FAmt_R(i) - FA_0) <= 1e-12*FA_0) then ; vBT_NN(i,J) = 0.0 ; else
          vBT_NN(i,J) = (1.5 * (dvR(i) - dv0(i,j))) * &
                      ((FAmt_R(i) - FA_avg) / (FAmt_R(i) - FA_0))
        endif
      else
        FA_v_S0(i,J) = 0.0 ; FA_v_SS(i,J) = 0.0
        FA_v_N0(i,J) = 0.0 ; FA_v_NN(i,J) = 0.0
        vBT_SS(i,J) = 0.0 ; vBT_NN(i,J) = 0.0
      endif
    enddo
  enddo

  !$omp target exit data map(release: dvL, dvR, dv_CFL, FAmt_L, FAmt_R, FAmt_0, vhtot_L, vhtot_R)

end subroutine set_merid_BT_cont_fortran

!> Shim for set_merid_BT_cont -- dispatches via SET_MERID_BT_CONT_MODE env var.
subroutine set_merid_BT_cont(bxC, v_a, h_in_a, h_S_a, h_N_a, BT_cont, dv0_a, vh_tot_0_a, &
                             dvhdv_tot_0_a, dv_max_CFL_a, dv_min_CFL_a, dt, &
                             dyCv_a, dx_Cv_a, IareaT_a, IdyT_a, CS, &
                             visc_rem_a, visc_rem_max_a, do_I_a, por_face_areaV_a)
  type(box_t),                                intent(in)    :: bxC  !< Iteration box for continuity solver
  type(RealArray_t),  intent(in)   :: v_a   !< Meridional velocity [L T-1 ~> m s-1].
  type(RealArray_t),  intent(in)   :: h_in_a !< Layer thickness used to
                                            !! calculate fluxes, [H ~> m or kg m-2].
  type(RealArray_t),  intent(in)   :: h_S_a !< South edge thickness in the
                                            !! reconstruction, [H ~> m or kg m-2].
  type(RealArray_t),  intent(in)   :: h_N_a !< North edge thickness in the
                                            !! reconstruction, [H ~> m or kg m-2].
  type(BT_cont_type),                         intent(inout) :: BT_cont !< A structure with elements
                       !! that describe the effective open face areas as a function of barotropic flow.
  type(RealArray_t),  intent(in)   :: dv0_a  !< The barotropic velocity increment
                                            !! that gives 0 transport [L T-1 ~> m s-1].
  type(RealArray_t),  intent(in)   :: vh_tot_0_a !< The summed transport
                       !! with 0 adjustment [H L2 T-1 ~> m3 s-1 or kg s-1].
  type(RealArray_t),  intent(in)   :: dvhdv_tot_0_a !< The partial derivative
                       !! of du_err with dv at 0 adjustment [H L ~> m2 or kg m-1].
  type(RealArray_t),  intent(in)   :: dv_max_CFL_a !< Maximum acceptable value
                                                                          !!  of dv [L T-1 ~> m s-1].
  type(RealArray_t),  intent(in)   :: dv_min_CFL_a !< Minimum acceptable value
                                                                          !!  of dv [L T-1 ~> m s-1].
  real,                                       intent(in)    :: dt   !< Time increment [T ~> s].
  type(RealArray_t),  intent(in)   :: dyCv_a !< The grid cell's v-point y-extent [L ~> m].
  type(RealArray_t),  intent(in)   :: dx_Cv_a !< The grid cell's unblocked lengths of the
                                            !! v-faces of the h-cell [L ~> m].
  type(RealArray_t),  intent(in)   :: IareaT_a !< The grid cell's 1/areaT [L-2 ~> m-2].
  type(RealArray_t),  intent(in)   :: IdyT_a   !< The grid cell's 1/dyT [L-1 ~> m-1].
  type(transport_adjust_CS),           intent(in)    :: CS !< Options
                       !! controlling the transport adjustment and barotropic-consistency iteration.
  type(RealArray_t),  intent(in)   :: visc_rem_a !< Both the fraction of the
                       !! momentum originally in a layer that remains after a time-step
                       !! of viscosity, and the fraction of a time-step's worth of a barotropic
                       !! acceleration that a layer experiences after viscosity is applied [nondim].
                       !! Visc_rem is between 0 (at the bottom) and 1 (far above the bottom).
  type(RealArray_t), &
                 intent(in)    :: visc_rem_max_a !< Maximum allowable visc_rem [nondim]
  type(LogicalArray_t),  intent(in)  :: do_I_a !< A logical flag indicating
                                             !! which I values to work on.
  type(RealArray_t),  intent(in)   :: por_face_areaV_a !< fractional open
                                            !! area of V-faces [nondim]

  integer :: mode, rc
  type(RealArray_C) :: v_c, h_in_c, h_S_c, h_N_c
  type(RealArray_C) :: FA_v_S0_c, FA_v_N0_c, FA_v_SS_c, FA_v_NN_c, vBT_SS_c, vBT_NN_c
  type(RealArray_C) :: dv0_c, vh_tot_0_c, dvhdv_tot_0_c, dv_max_CFL_c, dv_min_CFL_c
  type(RealArray_C) :: dyCv_c, dx_Cv_c, IareaT_c, IdyT_c, visc_rem_c, visc_rem_max_c
  type(RealArray_C) :: por_face_areaV_c
  type(LogicalArray_C)        :: do_I_c
  type(Box_C)                 :: bxC_c
  type(transport_adjust_CS_C) :: CS_c
  type(io_recorder)  :: rec
  logical            :: capture
  character(len=80)  :: kernel
  character(len=100) :: dir
  character(len=256) :: binFile, metaFile

  kernel = "set_merid_bt_cont"

  mode = getenv_mode("SET_MERID_BT_CONT_MODE", default=TIMH_runFORTRAN)

  select case (mode)

    case (TIMH_capture)
      capture = (.not. already_recorded(trim(kernel))) .and. is_root_pe()
      if (capture) then
        dir = "capture"
        rc = mkdir_posix(trim(dir) // c_null_char, int(o'755', c_int))
        binFile  = trim(dir) // "/" // trim(kernel) // ".bin"
        metaFile = trim(dir) // "/" // trim(kernel) // ".meta"
        call rec%open_write(binFile, metaFile)
        call rec%add("_bxC",            bxC)
        call rec%add("_v",              v_a)
        call rec%add("_h_in",           h_in_a)
        call rec%add("_h_S",            h_S_a)
        call rec%add("_h_N",            h_N_a)
        call rec%add("_FA_v_S0_before", BT_cont%FA_v_S0)
        call rec%add("_FA_v_N0_before", BT_cont%FA_v_N0)
        call rec%add("_FA_v_SS_before", BT_cont%FA_v_SS)
        call rec%add("_FA_v_NN_before", BT_cont%FA_v_NN)
        call rec%add("_vBT_SS_before",  BT_cont%vBT_SS)
        call rec%add("_vBT_NN_before",  BT_cont%vBT_NN)
        call rec%add("_dv0",            dv0_a)
        call rec%add("_vh_tot_0",       vh_tot_0_a)
        call rec%add("_dvhdv_tot_0",    dvhdv_tot_0_a)
        call rec%add("_dv_max_CFL",     dv_max_CFL_a)
        call rec%add("_dv_min_CFL",     dv_min_CFL_a)
        call rec%add("_dt",             dt)
        call rec%add("_dyCv",           dyCv_a)
        call rec%add("_dx_Cv",          dx_Cv_a)
        call rec%add("_IareaT",         IareaT_a)
        call rec%add("_IdyT",           IdyT_a)
        call rec%add("_vol_CFL",        CS%vol_CFL)
        call rec%add("_visc_rem",       visc_rem_a)
        call rec%add("_visc_rem_max",   visc_rem_max_a)
        call rec%add("_do_I",           do_I_a)
        call rec%add("_por_face_areaV", por_face_areaV_a)
      endif

      call set_merid_BT_cont_fortran(bxC, v_a, h_in_a, h_S_a, h_N_a, BT_cont, dv0_a, vh_tot_0_a, &
                             dvhdv_tot_0_a, dv_max_CFL_a, dv_min_CFL_a, dt, &
                             dyCv_a, dx_Cv_a, IareaT_a, IdyT_a, CS, &
                             visc_rem_a, visc_rem_max_a, do_I_a, por_face_areaV_a)

      if (capture) then
        call rec%add("_FA_v_S0_after", BT_cont%FA_v_S0)
        call rec%add("_FA_v_N0_after", BT_cont%FA_v_N0)
        call rec%add("_FA_v_SS_after", BT_cont%FA_v_SS)
        call rec%add("_FA_v_NN_after", BT_cont%FA_v_NN)
        call rec%add("_vBT_SS_after",  BT_cont%vBT_SS)
        call rec%add("_vBT_NN_after",  BT_cont%vBT_NN)
        call rec%close()
        call mark_recorded(trim(kernel))
      endif

#ifdef _TIM
    case (TIMH_runAMREX)
      bxC_c            = bxC%to_c()
      v_c              = v_a%to_c()
      h_in_c           = h_in_a%to_c()
      h_S_c            = h_S_a%to_c()
      h_N_c            = h_N_a%to_c()
      FA_v_S0_c        = BT_cont%FA_v_S0%to_c()
      FA_v_N0_c        = BT_cont%FA_v_N0%to_c()
      FA_v_SS_c        = BT_cont%FA_v_SS%to_c()
      FA_v_NN_c        = BT_cont%FA_v_NN%to_c()
      vBT_SS_c         = BT_cont%vBT_SS%to_c()
      vBT_NN_c         = BT_cont%vBT_NN%to_c()
      dv0_c            = dv0_a%to_c()
      vh_tot_0_c       = vh_tot_0_a%to_c()
      dvhdv_tot_0_c    = dvhdv_tot_0_a%to_c()
      dv_max_CFL_c     = dv_max_CFL_a%to_c()
      dv_min_CFL_c     = dv_min_CFL_a%to_c()
      dyCv_c           = dyCv_a%to_c()
      dx_Cv_c          = dx_Cv_a%to_c()
      IareaT_c         = IareaT_a%to_c()
      IdyT_c           = IdyT_a%to_c()
      CS_c             = transport_adjust_CS_to_c(CS)
      visc_rem_c       = visc_rem_a%to_c()
      visc_rem_max_c   = visc_rem_max_a%to_c()
      do_I_c           = do_I_a%to_c()
      por_face_areaV_c = por_face_areaV_a%to_c()
      call turbotmp_set_merid_bt_cont_bridge(bxC_c, v_c, h_in_c, h_S_c, h_N_c, &
                                             FA_v_S0_c, FA_v_N0_c, FA_v_SS_c, FA_v_NN_c, &
                                             vBT_SS_c, vBT_NN_c, dv0_c, vh_tot_0_c, &
                                             dvhdv_tot_0_c, dv_max_CFL_c, dv_min_CFL_c, dt, &
                                             dyCv_c, dx_Cv_c, IareaT_c, IdyT_c, CS_c, &
                                             visc_rem_c, visc_rem_max_c, do_I_c, por_face_areaV_c)
#endif

    case default
      call set_merid_BT_cont_fortran(bxC, v_a, h_in_a, h_S_a, h_N_a, BT_cont, dv0_a, vh_tot_0_a, &
                             dvhdv_tot_0_a, dv_max_CFL_a, dv_min_CFL_a, dt, &
                             dyCv_a, dx_Cv_a, IareaT_a, IdyT_a, CS, &
                             visc_rem_a, visc_rem_max_a, do_I_a, por_face_areaV_a)

  end select

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
  call get_param(param_file, mdl, "MONOTONIC_CONTINUITY", CS%reconstruction_CS%monotonic, &
                 "If true, CONTINUITY_PPM uses the Colella and Woodward "//&
                 "monotonic limiter.  The default (false) is to use a "//&
                 "simple positive definite limiter.", default=.false.)
  call get_param(param_file, mdl, "SIMPLE_2ND_PPM_CONTINUITY", CS%reconstruction_CS%simple_2nd, &
                 "If true, CONTINUITY_PPM uses a simple 2nd order "//&
                 "(arithmetic mean) interpolation of the edge values. "//&
                 "This may give better PV conservation properties. While "//&
                 "it formally reduces the accuracy of the continuity "//&
                 "solver itself in the strongly advective limit, it does "//&
                 "not reduce the overall order of accuracy of the dynamic "//&
                 "core.", default=.false.)
  call get_param(param_file, mdl, "UPWIND_1ST_CONTINUITY", CS%reconstruction_CS%upwind_1st, &
                 "If true, CONTINUITY_PPM becomes a 1st-order upwind "//&
                 "continuity solver.  This scheme is highly diffusive "//&
                 "but may be useful for debugging or in single-column "//&
                 "mode where its minimal stencil is useful.", default=.false.)
  call get_param(param_file, mdl, "ETA_TOLERANCE", CS%transport_adjust_CS%tol_eta, &
                 "The tolerance for the differences between the "//&
                 "barotropic and baroclinic estimates of the sea surface "//&
                 "height due to the fluxes through each face.  The total "//&
                 "tolerance for SSH is 4 times this value.  The default "//&
                 "is 0.5*NK*ANGSTROM, and this should not be set less "//&
                 "than about 10^-15*MAXIMUM_DEPTH.", units="m", scale=GV%m_to_H, &
                 default=0.5*GV%ke*GV%Angstrom_m)

  call get_param(param_file, mdl, "VELOCITY_TOLERANCE", CS%transport_adjust_CS%tol_vel, &
                 "The tolerance for barotropic velocity discrepancies "//&
                 "between the barotropic solution and  the sum of the "//&
                 "layer thicknesses.", units="m s-1", default=3.0e8, scale=US%m_s_to_L_T)
                 ! The speed of light is the default.

  call get_param(param_file, mdl, "CONT_PPM_AGGRESS_ADJUST", &
                 CS%transport_adjust_CS%aggress_adjust,&
                 "If true, allow the adjusted velocities to have a "//&
                 "relative CFL change up to 0.5.", default=.false.)
  CS%transport_adjust_CS%vol_CFL = CS%transport_adjust_CS%aggress_adjust
  call get_param(param_file, mdl, "CONT_PPM_VOLUME_BASED_CFL", CS%transport_adjust_CS%vol_CFL, &
                 "If true, use the ratio of the open face lengths to the "//&
                 "tracer cell areas when estimating CFL numbers.  The "//&
                 "default is set by CONT_PPM_AGGRESS_ADJUST.", &
                 default=CS%transport_adjust_CS%aggress_adjust, &
                 do_not_read=CS%transport_adjust_CS%aggress_adjust)
  call get_param(param_file, mdl, "CONTINUITY_CFL_LIMIT", &
                 CS%transport_adjust_CS%CFL_limit_adjust, &
                 "The maximum CFL of the adjusted velocities.", units="nondim", &
                 default=0.5)
  call get_param(param_file, mdl, "CONT_PPM_BETTER_ITER", CS%transport_adjust_CS%better_iter, &
                 "If true, stop corrective iterations using a velocity "//&
                 "based criterion and only stop if the iteration is "//&
                 "better than all predecessors.", default=.true.)
  call get_param(param_file, mdl, "CONT_PPM_USE_VISC_REM_MAX", &
                 CS%transport_adjust_CS%use_visc_rem_max, &
                 "If true, use more appropriate limiting bounds for "//&
                 "corrections in strongly viscous columns.", default=.true.)
  call get_param(param_file, mdl, "CONT_PPM_MARGINAL_FACE_AREAS", &
                 CS%transport_adjust_CS%marginal_faces, &
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

  stencil = 3 ; if (CS%reconstruction_CS%simple_2nd) stencil = 2
  if (CS%reconstruction_CS%upwind_1st) stencil = 1

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
