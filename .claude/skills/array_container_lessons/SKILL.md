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
  ! Pure outputs: allocate ONLY -- do not read undefined memory (§6).
  call h_W_a%alloc(lb=LBOUND(h_W), ub=UBOUND(h_W))
  call h_E_a%alloc(lb=LBOUND(h_E), ub=UBOUND(h_E))

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

> The in-tree exemplar calls `alloc(..., source=h_W)` on its
> `intent(out)` arrays. That reads undefined memory. Prefer the
> allocate-only form above for pure outputs — see §6 and §9 #2.

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
| `intent(out)` | yes | **no** — omit `source` | **yes** | `intent(inout)` |
| `intent(inout)` | yes | **yes** — `source=x` | **yes** | `intent(inout)` |

The `intent(out)` row is the one the in-tree exemplars get wrong. An
`intent(out)` dummy is undefined on entry, so `source=x` copies garbage.
It is numerically harmless *when the callee fully overwrites the array*,
but it is wasted work on every call. Two correct spellings:

```fortran
call h_W_a%alloc(lb=LBOUND(h_W), ub=UBOUND(h_W))                  ! allocate only
call h_W_a%allocView(h_W_p, lb=LBOUND(h_W), ub=UBOUND(h_W))       ! allocate + get pointer
```

**Caution:** only omit the copy-in when the callee genuinely writes every
element the caller will later read. If the callee writes an interior
region and the caller reads the halo, the array is *effectively*
`intent(inout)` regardless of how it is declared — copy it in, and
consider whether the declared intent is itself wrong.

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

Doc comments (`!<`) move with their dummy and are preserved verbatim,
including unit annotations like `[H ~> m or kg m-2]`.

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
2. **`source=` on a pure output reads undefined memory.** Use the
   allocate-only form (§6). Both in-tree exemplars have this issue; do
   not propagate it just because the examples do.
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
