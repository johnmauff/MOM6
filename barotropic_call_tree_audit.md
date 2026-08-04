# `MOM_barotropic.F90` (`btstep`) — call tree / isolation audit

Context: evaluating this file as a candidate for the same raw-array-to-container
conversion campaign already completed for `src/core/MOM_continuity_PPM.F90`.
Entry point requested: the barotropic solver, precise subroutine unknown up
front -- identified below as `btstep`.

## Summary

Good candidate, arguably the best of the three evaluated so far (continuity
solver aside). The call tree is genuinely tree-shaped -- a top-level driver
orchestrating a sequence of well-defined, single-purpose helper subroutines,
the same shape as `continuity_PPM` -- rather than one giant monolithic
subroutine (contrast with `horizontal_viscosity`, see the comparison table).
The one real difference from the other three files audited: `btstep` is called
from only 2 of the 4 dynamics-core drivers, not all 4.

## Entry point

`btstep`, `src/core/MOM_barotropic.F90`, lines 479-2370 (~1890 lines -- large,
but a dispatcher/orchestrator, not a monolith; see call tree below).

## Callers

`btstep` is called from:

- `src/core/MOM_dynamics_split_RK2.F90`
- `src/core/MOM_dynamics_split_RK2b.F90`

**Not** called from `MOM_dynamics_unsplit.F90` or `MOM_dynamics_unsplit_RK2.F90`.
This is architecturally expected -- the unsplit time-stepping schemes don't
separate barotropic and baroclinic time-stepping, so they have no need for a
separate barotropic solve -- but it does mean this module serves a narrower
slice of the dynamics core than the continuity solver, `MOM_vert_friction.F90`,
or `MOM_hor_visc.F90`, all three of which are called from all four driver
files.

## Call tree

Verified with two independent passes: an initial regex pass matching
`call X(...)` at the start of a line, then a second, more permissive pass
matching `call X(...)` anywhere on the line (this caught one inline
`if (cond) call X(...)` pattern the first pass missed -- `truncate_velocities`,
gated by `CS%clip_velocity`). The depth conclusion below held up after the
correction.

```
btstep
  -> BT_cont_to_face_areas       (same file, LEAF)
  -> adjust_local_BT_cont_types  (same file, LEAF)
  -> btstep_find_Cor             (same file, LEAF)
  -> btstep_layer_accel          (same file, LEAF)
  -> btstep_ubt_from_layer       (same file, LEAF)
  -> find_face_areas             (same file, LEAF)
  -> set_local_BT_cont_types     (same file, LEAF)
  -> set_up_BT_OBC               (same file, LEAF)
  -> btstep_timeloop             (same file)
       -> apply_u_velocity_OBCs  (same file, LEAF)
       -> apply_v_velocity_OBCs (same file, LEAF)
       -> btloop_add_dyn_PF      (same file, LEAF)
       -> btloop_eta_predictor   (same file, LEAF)
       -> btloop_find_PF         (same file, LEAF)
       -> btloop_update_u        (same file, LEAF)
       -> btloop_update_v        (same file, LEAF)
       -> truncate_velocities    (same file, LEAF; gated by CS%clip_velocity)
```

**3 levels deep, but wide**: `btstep` directly orchestrates 8 helpers, and its
one non-leaf callee (`btstep_timeloop`) orchestrates 8 more, every one of
which is confirmed to make no further subroutine calls. Same shape as
`continuity_PPM` (top-level driver -> sequence of single-purpose helpers), not
`horizontal_viscosity`'s single 2182-line subroutine.

### Pure-function helpers (not subroutine calls, but a recurring pattern)

A family of small `pure` functions is invoked via expression syntax
(e.g. `uhbt0(I,j) = uhbt(I,j) - find_uhbt(dt*ubt(I,j), BTCL_u(I,j)) * Idt`)
rather than `call`, extensively, inside several of the Level-2 subroutines
above (`btloop_eta_predictor`, `btloop_update_u`, `btloop_update_v`,
`set_up_BT_OBC`, and others):

- `find_uhbt`, `find_duhbt_du`, `uhbt_to_ubt`
- `find_vhbt`, `find_dvhbt_dv`, `vhbt_to_vbt`
- `swap`

All are leaves themselves (confirmed no further calls), so they don't add
depth -- but they're a recurring, heavily-used pattern that will need
classifying when this is converted, and look like natural `elemental`-kernel
candidates, similar to `flux_elem` in the continuity solver.

## Infra / diagnostics calls (not counted as "depth")

`Bchksum`, `hchksum`, `uvchksum`, `chksum0` (debug checksums); `MOM_error`,
`MOM_mesg` (error/messaging); `create_group_pass`, `do_group_pass`,
`start_group_pass`, `complete_group_pass`, `pass_var`, `pass_vector` (halo
exchange); `post_data` (diagnostics posting); `cpu_clock_begin`/`cpu_clock_end`
(timing); `enable_averages`/`enable_averaging`; `min_across_PEs` (MPI
reduction). Same flavor of infra as the continuity and vert-friction work.

## External (non-infra) physics dependencies

Four found, and **all four are gated behind optional-feature flags** -- same
pattern as `horizontal_viscosity`'s `use_GME`/`use_QG_Leith_visc`/`use_ZB2020`
gates:

| Call | Source module | Gate |
|---|---|---|
| `Filt_accum` | `MOM_streaming_filter.F90` | `if (CS%use_filter)` |
| `wave_drag_calc` | `MOM_wave_drag.F90` | `if (CS%use_filter .and. CS%linear_freq_drag)` |
| `HA_accum` | `MOM_harmonic_analysis.F90` | `if (associated(CS%HA_CSp) .and. find_etaav)` |
| `scalar_SAL_sensitivity` | `MOM_self_attr_load.F90` | `if (CS%calculate_SAL)` |

So the default/core execution path of `btstep` likely never leaves this file
at all, aside from ordinary infra.

`HA_accum` is also called once from `MOM.F90` (for a different field, `'ssh'`)
-- a simple accumulator call, not a deep entanglement.

## Isolation -- one nuance not present in the other three files

`MOM_barotropic.F90` is not purely single-audience the way continuity and
vert-friction are. Besides `btstep`, it also defines `barotropic_get_tav`,
which is called by `horizontal_viscosity` (`MOM_hor_visc.F90`), behind that
file's own `CS%use_GME` gate, to read back the current barotropic velocity
state. `barotropic_get_tav` is itself a leaf (no further calls), so it doesn't
add depth to `btstep`'s tree -- but it means this file has two distinct
external consumers reaching into two different entry points, not one.

Other subroutines in the file not reached by `btstep` at all:
`barotropic_init`, `barotropic_end`, `register_barotropic_restarts` (init/
teardown, expected -- same relationship `continuity_PPM_init` has to the
continuity solver), and `barotropic_get_tav` (covered above).

## Not yet container-converted

Zero occurrences of `RealArray_t`, `IntArray_t`, or `Box_t` in this file --- a
clean, fresh target. File is 6865 lines total (larger than continuity's ~4200
or vert-friction's ~3900), but that's a consequence of having ~25
subroutines/functions in a wide, shallow tree, not one oversized subroutine.

## Comparison table

| | Continuity solver (`continuity_PPM`) | Pressure solver (`Set_pbce_Bouss`) | Vertical friction (`vertvisc`/`vertvisc_coef`/`vertvisc_remnant`) | Lateral viscosity (`horizontal_viscosity`) | Barotropic solver (`btstep`) |
|---|---|---|---|---|---|
| Shape | Tree of subroutines | Shallow tree into a shared subsystem | Tree of subroutines | One 2182-line monolith | Tree of subroutines |
| Max physics call depth | 4-5 levels | 3-4 levels (TEOS10 path) | 2 levels | N/A -- inline branching, not depth | 3 levels |
| External callers | 4 dynamics-core files | 1 dispatcher (`MOM_PressureForce.F90`) | Same 4 dynamics-core files | Same 4 dynamics-core files | 2 of the 4 dynamics-core files (split schemes only) |
| Shared subsystem? | No -- self-contained | Yes -- EOS layer used by 48 files | No -- no other file calls in | No for the default path; yes if optional schemes enabled | No for the default path; `barotropic_get_tav` shared with `horizontal_viscosity`'s GME path |
| Runtime polymorphism? | No | Yes -- 9 possible EOS implementations | No | No | No |
| Vendored external library? | No | Yes -- GSW/TEOS10 toolbox (TEOS10 path only) | Yes -- CVMix, but only via `vertFPmix`, not reached by the 3 real entry points | No, but pulls in `MOM_barotropic.F90`, `MOM_thickness_diffuse.F90`, `MOM_lateral_mixing_coeffs.F90`, `MOM_Zanna_Bolton.F90` -- all gated behind optional-scheme flags | No |
| Already container-converted? | Fully (this campaign) | No | No | No | No |
