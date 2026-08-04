# `MOM_hor_visc.F90` (`horizontal_viscosity`) — call tree / isolation audit

Context: evaluating this file as a candidate for the same raw-array-to-container
conversion campaign already completed for `src/core/MOM_continuity_PPM.F90`.

## Summary

**Correction (see below): this is not purely a monolith.** An earlier pass
called `horizontal_viscosity` "essentially one enormous subroutine" with a
shallow, contained call tree. That was wrong for two of its three
optional-scheme branches: `calc_QG_slopes` (QG-Leith) and
`ZB2020_lateral_stress` (ZB2020) each fan out into real, multi-level call
trees of their own, one of which reaches all the way into the same
runtime-polymorphic equation-of-state dispatch layer flagged in the
pressure-solver audit. See "Call trees hidden behind the optional-scheme
calls" below for the corrected picture.

The top-level body itself is still one 2182-line subroutine rather than a
tree of small same-file helpers (contrast `continuity_PPM`/`btstep`), and
that size/branching-complexity problem (9+ config-gated code paths in one
body) is real and stands as originally written. But "not a call-tree
problem" was an overstatement -- the QG-Leith and ZB2020 branches mean a
container conversion here cannot stop at this file's boundary the way the
default/GME path suggested.

## Entry point

`horizontal_viscosity`, `src/parameterizations/lateral/MOM_hor_visc.F90`,
lines 270-2451 -- **2182 lines in a single subroutine**. For scale, that is
roughly half the entire 4042-line file, and far larger than any individual
subroutine encountered in the continuity, vertical-friction, or barotropic
audits.

## Shape: one monolith, not a tree

The rest of the module besides `horizontal_viscosity` itself:

- `hor_visc_init` / `hor_visc_end` -- init/teardown, not part of the runtime
  physics call (same relationship `continuity_PPM_init` has to the continuity
  solver, `barotropic_init`/`_end` to `btstep`).
- `align_aniso_tensor_to_grid` -- confirmed called only from `hor_visc_init`,
  not from `horizontal_viscosity`. Init-time only.
- `smooth_x9_h`, `smooth_x9_uv`, `smooth_GME` -- three small same-file helper
  subroutines, called from *within* `horizontal_viscosity` (confirmed via the
  call list below). These are the only same-file "tree" `horizontal_viscosity`
  actually has.

Rather than delegating most of the work to a call tree of smaller
subroutines (the shape seen in `continuity_PPM` and `btstep`), essentially
everything is inlined into the one subroutine body, which internally branches
on at least 9 distinct optional-scheme flags found by direct search:
`use_GME`, `use_QG_Leith_visc`, `use_ZB2020`, `use_Leithy`, `use_circulation`,
`use_beta_in_Leith`, `use_cont_thick`, `use_cont_thick_bug`, `use_land_mask`
(likely more that don't match a `use_*` naming pattern). It implements
Smagorinsky viscosity, Leith viscosity, biharmonic viscosity, MEKE-based
backscatter, GME backscatter, QG-Leith, ZB2020 (a newer ML-derived closure),
and an anisotropic-tensor variant, all inline in one body, each behind its own
`if (CS%use_...)` gate.

## Full call list from `horizontal_viscosity`

Verified with a permissive regex (matching `call X(...)` anywhere on the
line, not just at line start, the same check that caught a missed inline call
in the `btstep` audit) -- no additional external dependencies turned up beyond
what a first pass found:

```
Bchksum, MOM_error, ZB2020_copy_gradient_and_thickness, ZB2020_lateral_stress,
barotropic_get_tav, calc_QG_Leith_viscosity, calc_QG_slopes, hchksum,
pass_var, pass_vector, post_data, post_product_sum_u, post_product_sum_v,
post_product_u, post_product_v, smooth_GME, smooth_x9_h, smooth_x9_uv,
thickness_diffuse_get_KH, thickness_to_dz, uvchksum
```

## Infra / diagnostics calls (not counted as "depth")

`post_data` (**33 separate call sites** -- a lot of diagnostic output points
scattered through the subroutine), `pass_var`/`pass_vector` (halo exchange),
`MOM_error`, `post_product_u`/`post_product_v`/`post_product_sum_u`/
`post_product_sum_v` (diagnostics posting), `hchksum`/`Bchksum`/`uvchksum`
(debug checksums). Same flavor of infra as the other three files audited.
Worth flagging for whoever picks this up: the density of `post_data` calls
and the `post_product_*` calls suggest real diagnostics-posting logic
tightly interleaved with the physics throughout the subroutine, not
segregated to one section.

## External (non-infra) physics dependencies -- all gated behind optional-scheme flags

| Call | Source module | Gate | Leaf, or its own call tree? |
|---|---|---|---|
| `barotropic_get_tav` | `MOM_barotropic.F90` | `if (CS%use_GME)` | **Leaf** -- confirmed empty body of calls |
| `thickness_diffuse_get_KH` | `MOM_thickness_diffuse.F90` | `if (CS%use_GME)` | **Leaf** -- confirmed empty body of calls |
| `thickness_to_dz` | `MOM_interface_heights.F90` | `if (CS%use_QG_Leith_visc .and. (CS%Leith_Kh .or. CS%Leith_Ah))` | **Leaf** -- only calls `MOM_error` |
| `calc_QG_slopes` | `MOM_lateral_mixing_coeffs.F90` | same QG-Leith gate | **Not a leaf -- see below** |
| `calc_QG_Leith_viscosity` | `MOM_lateral_mixing_coeffs.F90` | same QG-Leith gate | Leaf apart from a `post_data` diagnostic |
| `ZB2020_lateral_stress` | `MOM_Zanna_Bolton.F90` | `if (CS%use_ZB2020)` | **Not a leaf -- see below** |
| `ZB2020_copy_gradient_and_thickness` | `MOM_Zanna_Bolton.F90` | `if (CS%use_ZB2020)` | Not yet traced further |

So the **default/core execution path (Smagorinsky/Leith/biharmonic, no
optional add-ons enabled) likely never leaves this file** apart from ordinary
infra, and the **GME branch is genuinely leaf-level** (both of its calls are
confirmed leaves). The QG-Leith and ZB2020 branches are a different story --
see the next section, which corrects the original "shallow, contained"
characterization of this file's call tree.

Separately, `thickness_to_dz` is otherwise a **widely-shared utility used by
~30 other files** across the model (tracer code, every vertical mixing
scheme, core dynamics, diagnostics, user initial conditions).
`horizontal_viscosity`'s own use of it is still conditional on
`use_QG_Leith_visc` and is itself a leaf call, so this doesn't add depth --
it's a breadth/shared-ownership concern, not a depth one.

## Call trees hidden behind the optional-scheme calls (correction)

Two of the three optional-scheme dependencies are not leaves -- each fans out
into a real, multi-level call tree, discovered by reading each callee's own
body rather than stopping at the first level (the same check applied to
`find_uhbt`/etc. in the barotropic audit, but missed on the first pass here):

**QG-Leith, via `calc_QG_slopes`** (`MOM_lateral_mixing_coeffs.F90:1360`) --
now traced to full depth, every node's body read directly:
```
horizontal_viscosity  [if CS%use_QG_Leith_visc]
  -> calc_QG_slopes                          (MOM_lateral_mixing_coeffs.F90)
       -> find_eta                           (generic; resolves to find_eta_3d)
            -> find_dz_for_eta               (MOM_interface_heights.F90)
                 -> int_specific_vol_dp      (MOM_density_integrals.F90)
                      -> int_spec_vol_dp_generic_pcm   (same file, LEAF apart
                      -> analytic_int_specific_vol_dp   from calling into...)
                           -> calculate_spec_vol        (MOM_EOS.F90 --
                              runtime-polymorphic EOS dispatch, same
                              mechanism as calculate_density)
                           -> int_spec_vol_dp_linear/wright/wright_full/
                              wright_red  (MOM_EOS.F90, concrete per-EOS-form
                              analytic leaves)
       -> calc_isoneutral_slopes             (MOM_isopycnal_slopes.F90)
            -> vert_fill_TS                  (same file, LEAF)
            -> calculate_density_derivs        (MOM_EOS.F90 -- confirmed
            -> calculate_density_second_derivs  dispatch via
                                                 EOS%type%calculate_density_derivs_2d,
                                                 the same runtime-polymorphic,
                                                 9-implementation dispatch
                                                 documented for the pressure
                                                 solver)
```
**At least 6 levels deep**, and it reaches the equation-of-state dispatch
layer through *two separate paths* (`find_eta` -> ... -> `int_specific_vol_dp`,
and `calc_isoneutral_slopes` -> `calculate_density_derivs`), not just one.
Every node above was read directly to confirm it either dispatches further or
is a genuine leaf -- this is a real, substantial call tree, gated behind
`CS%use_QG_Leith_visc` so it doesn't affect the default path, but not a leaf
by any measure.

**ZB2020, via `ZB2020_lateral_stress`** (`MOM_Zanna_Bolton.F90`) -- also now
traced to full depth:
```
horizontal_viscosity  [if CS%use_ZB2020]
  -> ZB2020_lateral_stress                   (MOM_Zanna_Bolton.F90)
       -> compute_c_diss                     (same file, LEAF -- cpu_clock only)
       -> filter_velocity_gradients          (same file)
            -> filter_hq                    (same file)
                 -> filter_3D                (same file, LEAF)
       -> compute_stress_ANN_collocated      (same file)
            -> ANN_apply_array_sio           (MOM_framework/MOM_ANN.F90)
                 -> layer_apply_sio          (same file, internal to
                                              ANN_apply_array_sio, LEAF)
          / compute_stress                   (same file, LEAF -- cpu_clock only)
       -> filter_stress                      (same file)
            -> filter_hq                    (same file, LEAF, as above)
       -> compute_stress_divergence          (same file)
            -> compute_energy_source        (same file, LEAF -- infra/post_data only)
  -> ZB2020_copy_gradient_and_thickness      (same file, LEAF -- cpu_clock only)
```
**4 levels deep**, entirely self-contained within `MOM_Zanna_Bolton.F90` plus
one hop out to `MOM_ANN.F90` (a small neural-net inference module) for the
`compute_stress_ANN_collocated` variant. Every node confirmed by reading its
body. Shallower than the QG-Leith branch and does not touch the EOS layer,
but still a genuine multi-level tree, not a leaf call.

**Net correction, now fully verified:** the depth/isolation picture for
`horizontal_viscosity` is "one large monolith at the top, whose GME branch is
genuinely leaf-level, but whose QG-Leith branch (6+ levels, reaching the
EOS dispatch layer via two separate paths) and ZB2020 branch (4 levels,
self-contained plus `MOM_ANN.F90`) are real, fully-traced, multi-level call
trees." Any container conversion here would need to scope in both subtrees
in full, not treat them as single-call dependencies.

## Isolation

- Called from the same four dynamics-core files seen throughout this
  campaign: `MOM_dynamics_unsplit.F90`, `MOM_dynamics_unsplit_RK2.F90`,
  `MOM_dynamics_split_RK2b.F90` (2 call sites), `MOM_dynamics_split_RK2.F90`
  (2 call sites). No other file calls it.
- Every external, non-infra call target is uniquely consumed by
  `horizontal_viscosity` alone at the call-site level (each of the 7 calls in
  the table above has exactly one call site in the whole repo, all inside this
  subroutine) -- but three of those targets live in substantial modules of
  their own (`MOM_barotropic.F90`, `MOM_thickness_diffuse.F90`,
  `MOM_Zanna_Bolton.F90`) that would need separate handling if those optional
  paths are ever converted.
- No runtime polymorphism.
- Not yet touched by any container conversion (zero `RealArray_t`/
  `IntArray_t`/`Box_t` references).

## Comparison table

| | Continuity solver (`continuity_PPM`) | Pressure solver (`Set_pbce_Bouss`) | Vertical friction (`vertvisc`/`vertvisc_coef`/`vertvisc_remnant`) | Lateral viscosity (`horizontal_viscosity`) | Barotropic solver (`btstep`) |
|---|---|---|---|---|---|
| Shape | Tree of subroutines | Shallow tree into a shared subsystem | Tree of subroutines | One 2182-line top-level body; but QG-Leith and ZB2020 branches open into real multi-level subtrees | Tree of subroutines |
| Max physics call depth | 4-5 levels | 3-4 levels (TEOS10 path) | 2 levels | Default/GME path: 1 level (leaf calls only). QG-Leith path: 6+ levels, fully traced (reaches EOS dispatch via two separate paths). ZB2020 path: 4 levels, fully traced (self-contained + one hop to `MOM_ANN.F90`) | 3 levels |
| External callers | 4 dynamics-core files | 1 dispatcher (`MOM_PressureForce.F90`) | Same 4 dynamics-core files | Same 4 dynamics-core files | 2 of the 4 dynamics-core files (split schemes only) |
| Shared subsystem? | No -- self-contained | Yes -- EOS layer used by 48 files | No -- no other file calls in | No for the default/GME path; yes for QG-Leith (reaches the same EOS layer as the pressure solver) and ZB2020 (own module, `MOM_Zanna_Bolton.F90`) | No for the default path; `barotropic_get_tav` shared with `horizontal_viscosity`'s GME path |
| Runtime polymorphism? | No | Yes -- 9 possible EOS implementations | No | No for default/GME; **yes for QG-Leith** (bottoms out in the same `EOS%type%density_elem` dispatch) | No |
| Vendored external library? | No | Yes -- GSW/TEOS10 toolbox (TEOS10 path only) | Yes -- CVMix, but only via `vertFPmix`, not reached by the 3 real entry points | No, but pulls in `MOM_barotropic.F90` (GME, leaf), `MOM_thickness_diffuse.F90` (GME, leaf), `MOM_lateral_mixing_coeffs.F90`/`MOM_isopycnal_slopes.F90`/`MOM_EOS.F90` (QG-Leith, multi-level), `MOM_Zanna_Bolton.F90` (ZB2020, multi-level) -- all gated behind optional-scheme flags | No |
| Already container-converted? | Fully (this campaign) | No | No | No | No |

## Recommendation

Still likely the lowest priority of the three parameterization files audited
(vert-friction, barotropic, this one) -- but for a different reason than
originally stated. The obstacle isn't a shallow-but-oversized single
subroutine; it's that a real conversion would need to scope in the QG-Leith
and ZB2020 subtrees as well (each a separate module of meaningful depth), on
top of the 2182-line top-level body's own 9+-way internal branching. Worth
considering splitting the optional-scheme branches (GME, QG-Leith, ZB2020,
anisotropic) out into separate subroutines *before* attempting a container
conversion, both for the top-level body's own sake and because it would let
QG-Leith's and ZB2020's subtrees be scoped and converted independently rather
than as part of one combined effort.
