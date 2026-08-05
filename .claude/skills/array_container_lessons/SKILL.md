---
name: array_container_lessons
description: Reference material for converting MOM6 subroutines from raw Fortran arrays (real, dimension(SZI_(G),SZJ_(G),SZK_(GV))) to the RealArray_t / IntArray_t / Box_t container types -- the in-place conversion pattern, where the marshalling lives at the call site, the complete alloc/allocView/view/copy2F/copy2Array/free API with exact signatures, intent-driven copy rules, and recurring pitfalls, organized in numbered sections §1-§10. Companion to the convert_array_containers skill, which cites these sections by number throughout its procedure. Invoke once near the start of a session that will run convert_array_containers one or more times; not needed again per-subroutine within that session.
user-invocable: true
---

# Array-container conversion: reference

> Invoke this skill once, near the start of a session, before running
> `convert_array_containers`. Its content then stays available for the
> rest of the session — you do not need to invoke it again before each
> subroutine. Sections are numbered §1–§10 and are referenced by that
> numbering from the calling skill; do not renumber them.

Worked reference in the tree, `src/core/MOM_continuity_PPM.F90`:

- **A converted callee:** `zonal_edge_thickness_fortran` (line ~603) and
  its mirror `meridional_edge_thickness_fortran` (~747) — container
  dummies, `%view` pointers, `do concurrent` over a `Box_t`.
- **A call site doing the marshalling:** `zonal_edge_thickness` (~636)
  and `meridional_edge_thickness` (~780) — `alloc` / `copy2F` / `free`
  around the call.
- **A fully container-native subtree:** `PPM_reconstruction_x/y`,
  `PPM_limit_pos`, `PPM_limit_cw84` and their `_fortran` halves — every
  one of these already speaks containers end to end.

**Scope of this reference — read this before anything else.** This
covers the Fortran-side *data-structure* conversion only:

- **In scope:** changing a subroutine's array dummies from
  `real, dimension(SZI_(G),…)` to `type(RealArray_t)`, rewriting its body
  to obtain pointers via `%view`, and updating its call sites to build,
  write back, and free the containers.
- **Explicitly out of scope:** renaming anything to `*_fortran`,
  introducing a dispatcher shim, `select case (mode)` / `getenv_mode`,
  capture mode / `io_recorder`, `#ifdef _TIM` arms, `bind(C)` interfaces,
  `%to_c`, and anything C++. Those belong to `generate_cpp_bridge` and
  `generate_amrex_code` and are layered on *afterwards*. See §10.

The `*_fortran` routines named above are useful examples of what a
container-native subroutine *looks like*, but the `_fortran` suffix and
the shim wrapped around them came from the bridge skill, not from a
conversion. **A conversion never renames a subroutine.**

---

## 1. The big picture — moving the boundary, one routine at a time

The refactor drives a boundary down through the call tree. Above it,
routines pass native Fortran arrays shaped by the grid macros
(`SZI_(G)`, `SZJ_(G)`, `SZK_(GV)`, …). Below it, everything speaks
containers.

Converting subroutine `FOO` moves that boundary up past `FOO`:

1. **`FOO` keeps its name.** Its array dummies become
   `type(RealArray_t)` / `type(IntArray_t)`, and its body obtains
   ordinary Fortran pointers from them via `%view`. The math is
   untouched.
2. **Every caller of `FOO` takes on the marshalling** it used to get for
   free: `alloc` a container per array argument, call `FOO`, `copy2F`
   the outputs back, `free` everything.

So the marshalling does not disappear — it moves up one level. Convert
the callers next and it moves up again. The boundary comes to rest
wherever the module's genuine public entry points are (in
`MOM_continuity_PPM.F90`, that is `continuity_PPM` and friends, which
have no `Box_t` and are probably meant to stay raw).

**Direction of travel.** One subroutine is converted per invocation;
nothing recurses. Both orders work, and they differ only in churn:

- **Top-down** — start at the highest routine you intend to convert. Its
  caller is the routine that will stay raw permanently, so the
  marshalling block is written once and never moves. Callees that are
  still raw receive `%view` pointers in the meantime (§9 #10 shows the
  same technique for `elemental` kernels), and each later conversion just
  swaps a pointer argument for the container it came from. **Preferred.**
- **Bottom-up** — leaves first. Each conversion adds a marshalling block
  to its caller, which is then deleted when that caller is itself
  converted. Correct, but every level's scaffolding is written and thrown
  away.

Either way, a converted routine calling a still-raw callee is a valid,
buildable intermediate state — pass the callee a `%view` pointer and note
it as the next conversion candidate.

**This preference is specifically about marshalling churn on the
container-*dummy* conversion — it is not the right direction for
deciding whether `G`/`GV`/`US` can be dropped from a signature
entirely.** Those are two different questions. Containerizing a
subroutine's own array dummies is a local decision — it needs nothing
from its callees, so top-down avoids rewriting scaffolding. Whether `G`
can be *dropped* from a subroutine's signature (`convert_array_containers`
"Settle these decisions" #2, and its Step 2b) is not local: a subroutine
that still forwards `G`/`GV` wholesale to even one still-raw callee
cannot drop it, no matter how many of its own `G%`/`GV%` field references
have been promoted to containers. That fact is only knowable once every
callee below it has itself finished the same check — which is inherently
bottom-up information. Converting container dummies top-down and then
deciding `G`/`GV` elimination bottom-up (leaves first, working back up)
are not in conflict; the first is about where marshalling scaffolding
lands, the second is about what a signature is still forced to carry.

Because no signature *semantics* change — same values, same order, same
bounds — a conversion is **numerically inert**. Any numerical difference
afterwards is a bug, not an expected consequence.

---

## 2. Template — the converted subroutine (the callee)

Before:

```fortran
subroutine zonal_edge_thickness(bxC, h_in, h_W, h_E, G, GV, US, CS, OBC)
  type(box_t),             intent(in)    :: bxC
  type(ocean_grid_type),   intent(in)    :: G
  real,  dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h_in !< Tracer cell layer thickness [H ~> m or kg m-2].
  real,  dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out)   :: h_W  !< Western edge layer thickness [H ~> m or kg m-2].
  ...
  do concurrent (k=..., j=..., i=...)
    h_W(i,j,k) = h_in(i,j,k)
  enddo
```

After:

```fortran
subroutine zonal_edge_thickness(bxC, h_in_a, h_W_a, h_E_a, mask2dT_a, &
                                h_min, upwind_1st, monotonic, simple_2nd, OBC)
  type(Box_t),          intent(in)    :: bxC        !< Iteration box for continuity solver
  type(RealArray_t),    intent(in)    :: h_in_a     !< Tracer cell layer thickness [H ~> m or kg m-2]
  type(RealArray_t),    intent(inout) :: h_W_a      !< Western edge layer thickness [H ~> m or kg m-2]
  type(RealArray_t),    intent(inout) :: h_E_a      !< Eastern edge layer thickness [H ~> m or kg m-2]
  type(RealArray_t),    intent(in)    :: mask2dT_a  !< Cell land/ocean mask [nondim]
  real,                 intent(in)    :: h_min      !< Minimum layer thickness (2*Angstrom_H) [H ~> m or kg m-2]
  logical,              intent(in)    :: upwind_1st !< If true, use 1st-order upwind reconstruction
  type(ocean_OBC_type), pointer       :: OBC        !< Open boundaries control structure

  integer :: i, j, k
  type(Box_t) :: bx
  real, dimension(:,:,:), contiguous, pointer :: h_in, h_W, h_E

  bx = bxC%grow(dim=1, n=1)
  call h_in_a%view(h_in)
  call h_W_a%view(h_W)
  call h_E_a%view(h_E)
  do concurrent (k=bx%idxS(3):bx%idxE(3), j=bx%idxS(2):bx%idxE(2), i=bx%idxS(1):bx%idxE(1))
    h_W(i,j,k) = h_in(i,j,k) ; h_E(i,j,k) = h_in(i,j,k)
  enddo
  call bx%free()
end subroutine zonal_edge_thickness
```

Key points:

- **The name is unchanged.** No `_fortran`, no wrapper.
- **The loop body is unchanged.** `%view` restores the container's
  `lb:ub` onto the pointer, so `h_W(i,j,k)` means exactly what it meant
  before. If a conversion diff shows edits inside the math, something is
  wrong.
- Container dummies are named `<orig>_a`; the `%view` pointer keeps the
  **original** name, which is what makes the body edit-free (§7).
- Scalars, `logical`s and `pointer` dummies (`OBC`) are **not**
  containerised — they pass through unchanged.
- `G`, `GV`, `US`, `CS` frequently disappear from the signature, replaced
  by the specific arrays and scalars the body actually used (§8).
- Written containers are `intent(inout)`, **never** `intent(out)` (§9 #6).

---

## 3. Template — the call site (the caller's new responsibility)

Every caller of a converted routine grows an alloc/call/write-back/free
block. Strip the dispatch machinery out of `zonal_edge_thickness`
(~636) and this is what remains:

```fortran
  real :: h_min
  type(RealArray_t) :: h_in_a, h_W_a, h_E_a, mask2dT_a

  h_min = 2.0 * GV%Angstrom_H

  ! Inputs: allocate AND copy in.
  call h_in_a%alloc(lb=LBOUND(h_in), ub=UBOUND(h_in), source=h_in)
  call mask2dT_a%alloc(lb=LBOUND(G%mask2dT), ub=UBOUND(G%mask2dT), source=G%mask2dT)
  ! Pure outputs: copy in too -- the write below only covers the Box_t's
  ! range, not h_W/h_E's full halo-inclusive extent (§6).
  call h_W_a%alloc(lb=LBOUND(h_W), ub=UBOUND(h_W), source=h_W)
  call h_E_a%alloc(lb=LBOUND(h_E), ub=UBOUND(h_E), source=h_E)

  call zonal_edge_thickness(bxC, h_in_a, h_W_a, h_E_a, mask2dT_a, &
                            h_min, CS%upwind_1st, CS%monotonic, CS%simple_2nd, OBC)

  ! Outputs and inouts: copy back out.
  call h_W_a%copy2F(h_W)
  call h_E_a%copy2F(h_E)

  call h_in_a%free()
  call h_W_a%free()
  call h_E_a%free()
  call mask2dT_a%free()
```

Fixed ordering: **derive scalars → alloc → call → `copy2F` → free.** If
the call site sits inside a `cpu_clock_begin`/`cpu_clock_end` pair, keep
the clock boundaries where they already are; do not move timing around a
conversion.

**This block is only for a caller that still takes raw arrays.** If the
caller is *already* container-based, it needs no marshalling at all — it
passes its own containers straight through:

```fortran
  ! Caller already holds h_in_a / h_W_a / h_E_a as its own dummies.
  call zonal_edge_thickness(bxC, h_in_a, h_W_a, h_E_a, mask2dT_a, &
                            h_min, upwind_1st, monotonic, simple_2nd, OBC)
```

Only an array the caller has no container for — a local work array, say —
still needs its own `alloc`/`free` pair.

So the block above is either **permanent** (the caller is a module entry
point that will stay raw) or **temporary scaffolding** (the caller is
itself queued for conversion, at which point the block is deleted and the
call becomes the pass-through form). Converting top-down avoids ever
writing the temporary kind — see §1, "Direction of travel".

> `h_W`/`h_E` are `intent(out)` in the original signature, yet the
> template above still copies them in with `source=`. That's deliberate,
> not an oversight — see §6 and §9 #2 for why "pure output" doesn't mean
> "safe to skip `source=`" in this codebase.

---

## 4. Container API — exact signatures

`RealArray_t` / `IntArray_t` live in `config_src/infra/FMS2/array_mod.F90`
and `config_src/infra/TIM/array_mod.F90`. The two infra layers are
**argument-for-argument identical** for everything in this section.

### 4.1 `alloc` (generic)

Five specifics, disambiguated by the **rank of `source`**:

```fortran
! Scalar / no-source form -- `source` is OPTIONAL and SCALAR
subroutine allocReal(this, dims, lb, ub, source)
  class(RealArray_t), intent(inout) :: this
  integer, intent(in), optional :: dims(:)
  integer, intent(in), optional :: lb(:)
  integer, intent(in), optional :: ub(:)
  real(kind=real64), intent(in), optional :: source   ! fills every element

! Rank-N forms (1D/2D/3D/4D) -- `source` is REQUIRED and is an ARRAY
subroutine allocReal3D(this, dims, lb, ub, source)
  class(RealArray_t), intent(inout) :: this
  integer, intent(in), optional :: dims(:)
  integer, intent(in), optional :: lb(:)
  integer, intent(in), optional :: ub(:)
  real(kind=real64), intent(in) :: source(:,:,:)      ! deep copy-in
```

Four usable call forms:

| Call | Effect |
|---|---|
| `call a%alloc(lb=LBOUND(x), ub=UBOUND(x), source=x)` | allocate to `x`'s bounds **and deep-copy `x` in** |
| `call a%alloc(lb=LBOUND(x), ub=UBOUND(x))` | allocate only; contents **undefined** |
| `call a%alloc(lb=LBOUND(x), ub=UBOUND(x), source=0.0)` | allocate and fill with a scalar |
| `call a%alloc(dims=[n1,n2,n3])` | allocate with `lb=1`, `ub=dims` |

Rules enforced at runtime:

- Supply **either** `lb`+`ub` **or** `dims`, never both, never neither →
  otherwise `MOM_err(FATAL, "allocReal: Must specify either ub and lb or dims")`.
- `size(lb) /= size(ub)` → `FATAL "allocReal: size of lb and ub must match"`.
- `alloc` **deallocates any existing payload first**, so re-allocating a
  live container is safe and does not leak.

The array-`source` forms are what an older API spelled `%dup` followed by
`%copy2Array`; that functionality was folded into `alloc`. See §4.6.

### 4.2 `allocView` (generic) — allocate and view in one call

```fortran
subroutine allocViewReal3D(this, a, dims, lb, ub, source)
   class(RealArray_t), intent(inout) :: this
   real(kind=real64), intent(inout), pointer :: a(:,:,:)   ! 2nd POSITIONAL arg
   integer, intent(in), optional :: dims(:)
   integer, intent(in), optional :: lb(:)
   integer, intent(in), optional :: ub(:)
   real(kind=real64), intent(in), optional :: source       ! SCALAR, optional
```

The pointer is the **second positional argument**, before
`dims`/`lb`/`ub`, and `source` here is scalar-and-optional even on the
rank-N specifics. This is the cleanest primitive for a pure output:
allocate storage and get a writable pointer without ever reading the
caller's undefined array.

### 4.3 `view` (generic) — zero-copy pointer remap

```fortran
subroutine viewReal3D(this, a)
   class(RealArray_t), intent(in) :: this
   real(kind=real64), intent(inout), pointer :: a(:,:,:)
```

Body is literally `a(lb(1):ub(1), lb(2):ub(2), lb(3):ub(3)) => this%data`.
**No allocation, no data movement.** Consequences:

- The pointer comes back carrying the container's original `lb:ub`, so
  index expressions in the ported loop body are unchanged.
- `this` is `intent(in)`, so `view` may be called on an `intent(in)`
  container and still yields a **writable** pointer. The in-tree code
  relies on this.
- Overloaded on the **rank of `a`** (1-D…4-D, real and integer). Rank is
  checked at runtime against `this%rank` and FATALs on mismatch.
- The dummy `a` is not declared `CONTIGUOUS`; passing a
  `contiguous, pointer` actual is legal and is what the in-tree code does.

### 4.4 `copy2F` and `copy2Array` — directions

```fortran
subroutine copy2FReal3D(this, var)          ! container --> native Fortran array
  class(RealArray_t), intent(in) :: this
  real, dimension(:,:,:), intent(inout) :: var

subroutine copy2AReal3D(this, var)          ! native Fortran array --> container
  class(RealArray_t), intent(inout) :: this
  real, dimension(:,:,:), intent(in) :: var
```

Mnemonic: `copy2F` = *copy to Fortran* (out of the container);
`copy2Array` = *copy to Array container* (into it).

- `copy2F` exists for ranks 1–4. There is **no 0-D `copy2F`**.
- `copy2Array` exists for ranks **0**–4; the 0-D form broadcasts a scalar.
- `copy2Array` does **not** allocate — the container must already be
  allocated, or it will dereference a null `shape`.

### 4.5 `free`

```fortran
subroutine freeReal(this)
  class(RealArray_t), intent(inout) :: this
```

Deallocates payload plus all three metadata arrays and resets `rank = 0`.
Every deallocate is guarded by `associated`, so **`free` is safe on an
unallocated or already-freed container**, and calling it twice is safe.

### 4.6 `dup` no longer exists

`%dup` was restructured into `%alloc`: what used to be
`call a%dup(x) ; call a%copy2Array(x)` is now the single call
`call a%alloc(lb=LBOUND(x), ub=UBOUND(x), source=x)`. There is no `dup`
type-bound procedure in either infra layer — `grep -rn "%dup"` over the
repo's `.F90` returns zero hits. Any documentation still showing `%dup`
predates that change and will not compile.

---

## 5. `Box_t` — iteration domain

Defined in `config_src/infra/{FMS2,TIM}/box_mod.F90`, which are
**byte-identical**. Components `idxS(:)` / `idxE(:)` are `allocatable`
and public. All bindings are plain `procedure ::` — no generics, so each
name has exactly one signature.

| Call | Semantics |
|---|---|
| `call bx%safe_alloc(ndims)` | Deallocate if needed, then allocate `idxS`/`idxE` to `ndims`, zero-filled. Idempotent. |
| `call bx%set(idxS=[...], idxE=[...])` | Element-wise assign. **FATALs if not already allocated** — `safe_alloc` must come first. |
| `call bx%free()` | Deallocate both components. Safe on an unallocated box. |
| `new = bx%grow(dim, n)` | **Returns a NEW `Box_t` by value.** `idxS(dim) -= n`, `idxE(dim) += n`. Does not mutate `bx`. |
| `new = bx%growLo(dim, n)` / `growHi(dim, n)` | As `grow`, but only the start / only the end. |
| `new = bx%shrink(dim, n)` | New box with `idxS(dim) += n`, `idxE(dim) -= n`. |
| `cdesc = bx%to_c()` | `Box_C` of `c_ptr`s. Null-safe. Available under **both** infra layers. |

`dim` is 1-based and is **not** bounds-checked against `SIZE(idxS)`.

Axis convention in this module: **`dim=1` is the zonal / i / x axis,
`dim=2` is the meridional / j / y axis.** This is the single most
important difference between a zonal routine and its meridional mirror.

Standard iteration idiom:

```fortran
do concurrent (k=bx%idxS(3):bx%idxE(3), j=bx%idxS(2):bx%idxE(2), i=bx%idxS(1):bx%idxE(1))
```

Because `grow`/`shrink` **return a new box**, the result is a fresh
allocation you own:

```fortran
bx = bxC%grow(dim=1, n=1)
...
call bx%free()          ! required -- see §9 #3
```

---

## 6. Intent-driven copy rules

What a dummy's original `intent` implies for the call-site block (§3):

| Original intent | Allocate | Copy in | Copy back (`copy2F`) | Container dummy intent |
|---|---|---|---|---|
| `intent(in)` | yes | **yes** — `source=x` | no | `intent(in)` |
| `intent(out)` | yes | **default yes** — `source=x` (see below) | **yes** | `intent(inout)` |
| `intent(inout)` | yes | **yes** — `source=x` | **yes** | `intent(inout)` |

**`intent(out)` does not mean "the callee writes every element."** It
only means Fortran treats the array as having no meaningful value on
entry — it says nothing about how much of the array actually gets
written, and for a plain (non-allocatable, non-pointer) array, `intent(out)`
does **not** clear memory. Before conversion, the raw array was the
caller's actual memory, passed by reference: whatever the callee didn't
write, the caller kept from before (a previous halo exchange, a prior
timestep, whatever was already there).

**This matters concretely in this codebase because almost every array is
allocated wider than any single call writes it.** Arrays shaped by
`SZI_`/`SZIB_`/`SZJ_`/`SZK_` are sized over the full halo-inclusive data
domain (`G%isd:G%ied`, `G%IsdB:G%IedB`, ...), while a write scoped to a
`Box_t` built from `set_continuity_box` only covers the narrower
*computational* domain (`G%isc:G%iec`, optionally stencil-widened) —
that gap between compute domain and data domain is the halo, and it
exists precisely so most calls *don't* touch it. So a container
allocated without `source=` for a `Box_t`-scoped write starts as
uninitialized memory over the *entire* declared extent, gets real values
written only inside the box, and then `copy2F` copies the whole thing —
including that never-written, garbage-filled halo — back over the
caller's array, silently destroying whatever meaningful halo content was
there before. This is not hypothetical: it is the root cause of a real
CI-breaking numerical failure in an already-converted, already-merged
routine (`zonal_mass_flux`'s `uh`) that took several sessions to trace,
specifically *because* it doesn't crash — the corrupted halo cells sit
quietly wrong until something downstream reads them, often much later.

**The corrected default: copy in with `source=x`, even for `intent(out)`,
unless you can positively confirm the callee writes the array's entire
declared extent** (not just the `Box_t`'s range — the actual `LBOUND`/
`UBOUND`). For any write scoped by a `Box_t`, that confirmation essentially
never holds in this codebase, so treat `intent(out)` the same as
`intent(inout)` for copy-in purposes by default. The allocate-only forms
below are real, but they are the exception, not the default — use them
only when you have checked, specifically, that the write loop's bounds
equal the array's full `LBOUND`/`UBOUND` in every dimension:

```fortran
call h_W_a%alloc(lb=LBOUND(h_W), ub=UBOUND(h_W))                  ! only if h_W is written in full
call h_W_a%allocView(h_W_p, lb=LBOUND(h_W), ub=UBOUND(h_W))       ! same caveat -- no copy-in either
```

If in doubt, copy in. The wasted copy of a few halo cells is cheap; a
silently corrupted halo is not.

### Optional dummies

An `optional` array dummy keeps its `optional` attribute when it becomes
a container, and the `present()` machinery around it has to move with it.

**In the callee**, three things change together:

```fortran
  type(RealArray_t), optional, intent(in) :: hin_a  !< Initial layer thickness [H ~> m or kg m-2].
                                                    !! If hin is absent, h is also the initial thickness.
  real, dimension(:,:,:), contiguous, pointer :: hin

  if (present(hin_a)) call hin_a%view(hin)      ! guard the view
  ...
  if (present(hin_a)) then                      ! existing logic now tests the _a name
```

1. Keep `optional` on the container declaration.
2. **Guard the `%view`.** Calling `%view` on an absent dummy is illegal.
3. **Rewrite every existing `present(<orig>)` test to `present(<orig>_a)`.**
   This is the easiest thing to miss: the body's control flow already
   references the old name, and a stale `present(hin)` will not compile
   once `hin` is a local pointer rather than a dummy.

The pointer stays unassociated when the dummy is absent, so any use of it
must sit inside the same `present` guard.

**At the call site**, it depends on whether the caller's own array is
optional:

- **Caller's array is always present** (the common case — it is one of the
  caller's ordinary dummies or a local): allocate and pass it normally.
  Nothing special.

- **Caller's array is itself `optional`**: a container local is *always*
  present as an actual argument, so passing an unallocated container
  would make `present()` return `.true.` in the callee with meaningless
  data behind it. Allocate conditionally and branch the call:

  ```fortran
  if (present(hin)) then
    call hin_a%alloc(lb=LBOUND(hin), ub=UBOUND(hin), source=hin)
    call continuity_zonal_convergence(bxC, h_a, uh_a, dt, IareaT_a, hin_a=hin_a)
  else
    call continuity_zonal_convergence(bxC, h_a, uh_a, dt, IareaT_a)
  endif
  call hin_a%free()      ! safe even if never allocated (§4.5)
  ```

  With one or two optional arrays this is fine. If a call site has more
  than two, the branch count explodes combinatorially — **stop and ask
  the user** rather than generating a thicket of nested conditionals.

Optional *scalars* (`hmin`, `h_min`) are unaffected — they are not
containerised and pass through with their `optional` attribute and
`present()` tests unchanged.

---

## 7. Naming conventions

| Thing | Convention | Example |
|---|---|---|
| Subroutine name | **unchanged** | `zonal_mass_flux` stays `zonal_mass_flux` |
| Container dummy | original name + `_a` | `h_in` → `h_in_a` |
| `%view` pointer inside the body | the **original** name | `h_in` |
| Container local at a call site | original name + `_a` | `h_W_a` |
| Grid-derived container | field name + `_a` | `G%mask2dT` → `mask2dT_a` |

Giving the `%view` pointer the original name is what lets the loop body
stay edit-free. A clean conversion diff touches declarations and adds
`%view` calls — nothing inside the math.

### Doc comments

Every dummy in this codebase carries a `!<` Doxygen comment, and **every
one of them must survive a conversion verbatim** — same wording, same
unit annotation (`[H ~> m or kg m-2]`, `[nondim]`), same trailing
punctuation.

The common way this breaks: in the original, a grid-shaped array
declaration usually spans two lines with the comment on the
*continuation* line, whereas the converted declaration fits on one. The
comment is easy to lose in that collapse.

```fortran
! BEFORE -- two lines, comment on the continuation
  real,  dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h_in !< Tracer cell layer thickness [H ~> m or kg m-2].

! AFTER -- one line, comment carried across
  type(RealArray_t),    intent(in)    :: h_in_a     !< Tracer cell layer thickness [H ~> m or kg m-2]
```

Record each comment before rewriting the declaration, then re-attach it.
Dummies **added** by the conversion — grid-derived containers like
`mask2dT_a`, scalars lifted out of `CS` or `GV` — need a new `!<`
comment in the same house style; a missing one is a Doxygen build
warning, and the repo's CI treats doc coverage as part of the build.

Comments that are not attached to a dummy — the `!>` header above the
subroutine, blocks between declarations, anything in the body — stay
exactly where and as they are.

### No commentary about the conversion itself

This codebase documents ocean physics and numerics — never the fact
that a subroutine was migrated to containers, or why a particular
Fortran mechanism is being used here. Do not add comments like
`! Containers for zonal_mass_flux`, `! %view pointers for forwarding
this routine's optional arguments to the still-raw child`, or
`! Must be nullified -- see below` anywhere in converted code. The
container API and the `_a` naming convention already make the pattern
legible to anyone who has read this reference; a maintainer reading the
physics code should not need to re-derive the migration's own rationale
from prose scattered through it. This applies to declaration-block
labels, multi-line explanatory blocks, and anything in between — a
conversion diff should read as pure mechanism (types, `%view`/`%alloc`/
`%free` calls, renamed dummies), not narration.

The one exception, already covered above: a genuinely **new** dummy
still needs a `!<` doc comment describing its physical meaning, in house
style — exactly the same requirement as before, not an exception to this
rule.

---

## 8. Grid-derived arrays and shrinking argument lists

Some arrays a routine needs are not dummies at all — they are reached
through a derived type, most commonly `G%mask2dT`. These become **new**
container dummies, built by the caller:

```fortran
call mask2dT_a%alloc(lb=LBOUND(G%mask2dT), ub=UBOUND(G%mask2dT), source=G%mask2dT)
```

The converted routine then takes `mask2dT_a` and may no longer need `G`
at all. So a converted argument list is **not** simply the old list with
types swapped — it gains container dummies and often sheds `G` / `GV` /
`US` / `CS`.

Work through the original body and classify every derived-type
reference:

- **Array** (`G%mask2dT`, `G%areaT`) → new container dummy, allocated by
  the caller.
- **Scalar** (`GV%Angstrom_H`, `CS%upwind_1st`, `CS%monotonic`) → plain
  scalar dummy, computed or read at the call site
  (`h_min = 2.0 * GV%Angstrom_H`).
- **Still genuinely needed as a whole** (passed on to an unconverted
  callee, used for `US` unit scaling in a way that resists splitting) →
  keep the derived type in the signature. Do not force it.

Only drop `G`/`GV`/`US`/`CS` from the signature once nothing in the body
references them.

---

## 9. Recurring pitfalls

1. **Do not rename the subroutine.** No `_fortran` suffix, no wrapper
   pair. That transformation belongs to `generate_cpp_bridge` (§10).
2. **Skipping `source=` on a `Box_t`-scoped "pure output" silently
   corrupts its halo.** `intent(out)` does not mean "written in full" —
   it only means undefined on entry. Any array shaped by `SZI_`/`SZIB_`/
   etc. is allocated over the full halo-inclusive domain, while a write
   scoped to a `Box_t` covers only the narrower computational range in
   between. Allocate-without-`source=` starts the untouched halo as
   garbage; `copy2F` then copies that garbage back over the caller's
   real array, overwriting whatever meaningful halo content was there
   before. Default to `source=x` even for `intent(out)` (§6); this
   caused a real, hard-to-trace CI failure (`zonal_mass_flux`'s `uh`) —
   only skip the copy-in when you have confirmed the write covers the
   array's full `LBOUND`/`UBOUND`, not just the `Box_t`'s range.
3. **`grow`/`shrink`/`growLo`/`growHi` return a NEW box that you must
   `free`.** Forgetting `call bx%free()` leaks on every call — and these
   routines run inside timestep loops.
4. **`set` before `safe_alloc` is FATAL.** Always
   `call bx%safe_alloc(ndims=3)` before `call bx%set(...)`.
5. **`copy2Array` does not allocate.** The container must already exist.
   `alloc(..., source=x)` is the one-step form; `copy2Array` is for
   re-filling an existing container.
6. **Container dummies that are written must be `intent(inout)`, never
   `intent(out)`.** A container passed `intent(out)` arrives with its
   pointer components reset, discarding the allocation the caller just
   made. This inverts the raw-array convention and is easy to get
   backwards.
7. **Container extents must match the source array's extents.** The copy
   loops are driven by `this%shape`, not by `SHAPE(source)`; a mismatch
   reads or writes out of bounds instead of failing cleanly.
8. **`view` FATALs on rank mismatch at run time, not compile time.** A
   3-D container viewed through a 2-D pointer compiles fine and dies at
   run time.
8a. **`optional` dummies: rename the `present()` tests too.** When
   `hin` becomes `hin_a`, every `present(hin)` in the body must become
   `present(hin_a)`, and the `%view` must be guarded by it (§6). A
   missed rename does not compile once `hin` is a local pointer — but a
   missed *guard* compiles and fails at run time.
8b. **Never pass an unallocated container to an `optional` dummy.** A
   container local is always present as an actual argument, so
   `present()` in the callee returns `.true.` regardless. If the
   caller's source array is optional, branch the call (§6).
9. **Every call site must be updated in the same change.** Converting a
   subroutine's dummies without updating all its callers breaks the
   build. Find them with
   `grep -rn "call[[:space:]]\+<name>(" src/ config_src/` and check for
   continuation lines.
10. **`elemental` callees cannot take containers.** `flux_elem` /
    `flux_elem_OBC` are `elemental` scalar kernels invoked inside
    `do concurrent` loops. They keep their scalar signatures; the
    converted caller passes `%view` pointer elements to them, unchanged.
11. **Derived types with array members are not containers.**
    `BT_cont_type` (`BT_cont%h_u`, `%h_v`) holds raw arrays. Converting a
    routine that passes `BT_cont` around does not convert its members;
    leave the derived type alone unless explicitly asked.
12. **`%to_c` on the array types is TIM-only.** The FMS2 implementation
    is commented out and there is no binding, so every array `%to_c()`
    must sit inside `#ifdef _TIM`. (`Box_t%to_c` is available under
    both.) A conversion emits no `to_c` at all — this matters only when a
    bridge is added later.
13. **`IntArray_C` is not public** in either infra layer, so a caller
    cannot declare the result type of `to_c_Int`.
14. **Default `real` is assumed to be `real64`.** `copy2*` declare their
    `var` as default `real` while `alloc*` declare `source` as
    `real(kind=real64)`. This compiles only because MOM6 builds with
    default real = real64. Never emit code that deliberately mixes kinds
    across this boundary.
15. **Free every container on every path.** Early returns at a call site
    must still reach the `free` block.
16. **Forwarding one of the routine's own optional container dummies to
    a still-raw callee needs a nullified `%view` pointer, not a branch —
    and the `nullify` is mandatory, not defensive.** Once a dummy is
    `optional, type(RealArray_t)`, you cannot pass it directly to a
    callee whose matching dummy is still a raw optional array (type
    mismatch); the fix is a local pointer, guarded `%view`
    (`if (present(x_a)) call x_a%view(x)`), forwarded to the callee
    *unconditionally* by keyword. Fortran reads a disassociated pointer
    passed to a non-pointer, non-allocatable `optional` dummy as
    genuinely absent, so this is branch-free in both directions — but
    only if the pointer is nullified first. A plain local pointer's
    association status on subroutine entry is **undefined, not
    disassociated**; skip the `nullify` and the callee's `present()` can
    read `.true.` by accident off stack garbage whenever the container
    is absent, silently corrupting results several call levels
    downstream with no compile error and no crash — confirmed as a real,
    CI-caught bug in this codebase. Nullify every such pointer once,
    before the first guarded `%view` that might set it, regardless of
    which branch of the routine actually runs. Do **not** nullify inside
    the declaration (`pointer :: x => null()`) — that implicitly adds
    `SAVE`, so it only runs once ever and the pointer then carries
    whatever association it had at the end of the *previous* call,
    which is worse than not nullifying at all for any routine invoked
    more than once. This also requires the callee's own dummy to be a
    plain optional array — never pointer or allocatable; check every
    callee individually. See the `convert_optional_args_to_containers`
    skill for the full mechanism and worked example, including the
    caller side of this same forwarding problem.

---

## 10. Verification, and composing with the other skills

A conversion is **numerically inert by construction** (§1), so
verification is mostly about proving nothing moved:

1. **Diff review.** The math inside the converted routine should be
   unchanged. Any edit inside a loop body is a red flag. The diff should
   consist of: dummy declarations, added `%view` locals and calls, and
   the call-site alloc/copy2F/free blocks.
2. **Call-site completeness.** Every caller found by grep is updated;
   the argument count and order at each site matches the new signature.
3. **Both infra layers.** The conversion must compile under FMS2 *and*
   TIM. The container API is identical across them (§4), so anything
   that compiles under only one has reached outside the shared surface.
4. **Build.** Do not assume a Fortran compiler is present — many
   developer machines on this project do not have one. Build where one
   exists (e.g. Derecho), or state plainly that the build was not run and
   is pending. Never report a conversion as verified on the strength of
   inspection alone.
5. **Numerics.** Results must be bit-identical to the pre-conversion
   build.

**Composing with the bridge skills.** A converted routine is the natural
input to `generate_cpp_bridge`, which is what adds the `*_fortran`
rename, the dispatcher shim, `getenv_mode`, capture mode, the
`#ifdef _TIM` AMReX arm, and the `bind(C)` interface —
`generate_amrex_code` then writes the C++ side. Keeping the conversion
separate keeps each diff small and reviewable, and means a routine can be
containerised without any commitment to ever bridging it.

**A `_fortran` sibling doesn't change what this skill does.** A
subroutine with raw dummies is a conversion target regardless of its
history — whether it has never been touched, or is a bridge shim whose
own signature hasn't been migrated yet (a leftover from before the skill
split, when one combined skill produced the container worker and the
shim in the same pass). Convert it exactly like any other raw-array
subroutine (§2–§8). The one constraint: some of its content belongs to
`generate_cpp_bridge`, not this skill, and must be left alone —

- Never edit anything between `select case (mode)` and `end select`,
  except renaming a `CS%x`/`GV%x` reference if that derived type is
  dropped from the signature (§8). No other change inside that region —
  not the case structure, not `getenv_mode`, not the capture block, not
  any `%to_c()` call, not the `#ifdef _TIM` guard.
- Never touch a `bind(C) interface … end interface` declaration. It
  takes `RealArray_C`/`Box_C` regardless of where the `RealArray_t` came
  from, so a conversion never has a reason to reach it.

If you find yourself rewriting either, stop — that's a sign the task has
drifted into bridge territory rather than a data-structure conversion.
Container-building/writeback bookkeeping *around* that content (a
"build containers for all dispatch paths" block before it, a trailing
`copy2F`/`free` block after it) is ordinary conversion scaffolding and
gets replaced the same way it would at any call site (§3) — it's simply
inside this subroutine rather than in one of its callers.

`meridional_edge_thickness` (not yet migrated) and `zonal_edge_thickness`
(migrated) are an exact real-world before/after pair in this tree, if you
want to see it worked out in full.

The reliable signal that a raw-array subroutine has bridge content to
preserve is **not** its dummy types — a `<name>_fortran` sibling's
existence, and dispatch machinery in the body (`select case (mode)`,
`getenv_mode`, `TIMH_*`, `io_recorder`, `#ifdef _TIM`), which
`convert_array_containers`'s Step 0 already checks for.
