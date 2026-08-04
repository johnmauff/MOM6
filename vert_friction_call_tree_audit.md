# `MOM_vert_friction.F90` — call tree / isolation audit

Context: evaluating this file as a candidate for the same raw-array-to-container
conversion campaign already completed for `src/core/MOM_continuity_PPM.F90`.
Four subroutines were named as potential entry points: `vertvisc`,
`vertvisc_limit_vel`, `vertvisc_remnant`, `vertvisc_coef`.

## Summary

Good candidate. Call tree is shallow (2 levels of real physics depth, at most),
and the module is called from exactly the same four dynamics-core driver files
that called into the continuity solver — no other file in the repo depends on
it. This is a much closer match to the continuity solver's profile than
`MOM_PressureForce_Montgomery.F90`'s `Set_pbce_Bouss` was (that one immediately
fans out into the equation-of-state layer, which is shared by 48 files across
the model and involves runtime-polymorphic dispatch — see the comparison notes
below).

## One correction to the entry-point list

`vertvisc_limit_vel` is **not** an independent entry point. Nothing outside
`MOM_vert_friction.F90` calls it directly — it's called internally by
`vertvisc` itself (same relationship `zonal_flux_adjust` had to
`present_uhbt_or_set_BT_cont` before that was converted). So there are really
three external entry points here, not four: `vertvisc`, `vertvisc_coef`,
`vertvisc_remnant`, with `vertvisc_limit_vel` one level down the tree from
`vertvisc`.

## Call tree (physics-relevant calls only; infra omitted below, see next section)

```
vertvisc            -> vertvisc_limit_vel        (same file, internal helper)

vertvisc_coef       -> find_coupling_coef         (same file, LEAF -- no further calls)
                     -> find_coupling_coef_gl90    (same file, LEAF -- no further calls)
                     -> find_coupling_coef_k       (same file, LEAF -- no further calls)
                     -> find_ustar                 (MOM_forcing_type.F90; itself a LEAF,
                                                     only calls MOM_error for validation)

vertvisc_remnant    -> (no physics-relevant calls; only a debug checksum utility)
```

Maximum physics-relevant depth from any of the three real entry points is
**2 levels** (`vertvisc_coef` -> `find_ustar`, which is itself a leaf). This is
shallower than the continuity solver's tree, which ran 4-5 levels deep through
`flux_elem`.

## Infra / diagnostics calls (not counted as "depth", same treatment as
`cpu_clock_begin`/`MOM_error` in the continuity work)

- `vertvisc`: `post_data`, `post_product_u`, `post_product_v`,
  `post_product_sum_u`, `post_product_sum_v` (diagnostics posting),
  `create_group_pass`, `do_group_pass` (halo exchange), `MOM_error`.
- `vertvisc_coef`: `post_data`, `uvchksum` (debug checksum).
- `vertvisc_remnant`: `uvchksum` only.
- `vertvisc_limit_vel`: `write_u_accel`, `write_v_accel`
  (`MOM_PointAccel.F90` diagnostic writers).

Worth flagging for whoever picks this up: `vertvisc`'s own callees
(`post_product_u/v`, `post_product_sum_u/v`) suggest there's real
diagnostics-posting logic interleaved with the physics, which will need the
same "what's data vs. what's infra" classification work that `G`/`GV`/`CS`/`OBC`
needed in the continuity file. Not a red flag, just expected process overhead.

## Isolation

- `vertvisc`, `vertvisc_coef`, `vertvisc_remnant` (and a fourth subroutine in
  this file, `vertFPmix`, not one of the four named) are called **only** from:
  - `src/core/MOM_dynamics_unsplit_RK2.F90`
  - `src/core/MOM_dynamics_split_RK2b.F90`
  - `src/core/MOM_dynamics_unsplit.F90`
  - `src/core/MOM_dynamics_split_RK2.F90`

  These are the exact same four files that called into `continuity_PPM`,
  `continuity_3d_fluxes`, `continuity_2d_fluxes`, and `continuity_adjust_vel`.
  No other file in the repo calls any of `MOM_vert_friction`'s public entry
  points.
- No runtime polymorphism anywhere in this tree -- every call is a fixed,
  statically-determined subroutine call, same as the continuity solver. (Contrast
  with the equation-of-state layer, where `EOS%type%density_elem` dispatches at
  runtime to one of 9 possible concrete implementations.)
- There is one external, vendored-third-party-library dependency in this file:
  `cvmix_kpp_composite_Gshape`, imported from `CVMix_kpp` (the CVMix mixing
  library, same flavor of dependency as the GSW/TEOS10 toolbox used by the
  equation-of-state layer). **It is used only by `vertFPmix`**, which is not
  reached by any of the three real entry points (`vertvisc`, `vertvisc_coef`,
  `vertvisc_remnant`) -- confirmed `vertFPmix` is called from a separate site
  in `MOM_dynamics_split_RK2.F90` only, not from within `vertvisc`/`vertvisc_coef`/
  `vertvisc_remnant`. So the CVMix dependency does not touch the call tree in
  question.
- The file is currently **untouched by any container conversion** -- zero
  occurrences of `RealArray_t`, `IntArray_t`, or `Box_t` -- a clean, fresh
  target. It is 3902 lines, comparable in size to
  `MOM_continuity_PPM.F90` (~4200 lines after the continuity conversion work).

## Comparison table

| | Continuity solver (`continuity_PPM`) | Pressure solver (`Set_pbce_Bouss`) | Vertical friction (`vertvisc`/`vertvisc_coef`/`vertvisc_remnant`) |
|---|---|---|---|
| Max physics call depth | 4-5 levels | 3-4 levels (TEOS10 path) | 2 levels |
| External callers | 4 dynamics-core files | 1 dispatcher (`MOM_PressureForce.F90`) | Same 4 dynamics-core files |
| Shared subsystem? | No -- self-contained | Yes -- EOS layer used by 48 files | No -- no other file calls in |
| Runtime polymorphism? | No | Yes -- 9 possible EOS implementations | No |
| Vendored external library? | No | Yes -- GSW/TEOS10 toolbox (TEOS10 path only) | Yes -- CVMix, but only via `vertFPmix`, not reached by the 3 real entry points |
| Already container-converted? | Fully (this campaign) | No | No |
