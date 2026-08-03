---
name: hoist_container_marshalling
description: Reduce redundant %alloc/%free/%copy2F/%copy2Array churn in a MOM6 subroutine that already marshals RealArray_t/IntArray_t containers around calls to converted descendants (a Case-A caller, in convert_array_containers terms), and group the surviving alloc calls before -- and free/copy-back calls after -- the block of descendant calls, so that block reads and profiles as pure computation rather than memory bookkeeping. Use this on a subroutine that already has containers, and only when every descendant it calls is already container-native AND the subroutine's own dummy list has no optional raw array -- it refuses to run if even one call site still passes a raw array to a still-raw callee, or if the subroutine itself still takes an optional raw array (Step 0), since both cases pin a runtime-conditional alloc/free or copy2F that this skill cannot hoist away. A still-optional-and-raw dummy needs convert_optional_args_to_containers run on it first (often eliminating the branch entirely, not just relocating it, since a converted optional container forwards straight through to a same-typed optional dummy on the callee with no present() check needed). It never converts a new raw array, and it never touches a callee's signature or body. Best applied once every one of the subroutine's container conversions -- its own dummies, its own optional dummies, and every one of its descendants' -- are settled and every call site inside it is visible in one pass -- exactly the moment convert_array_containers' own Step 8/10 "hoisting candidate" deferral points to.
user-invocable: true
argument-hint: <work-directory> <function-name> [--enable_git_commit] [--disable_git_commit]
---

# Hoist container marshalling out of a Case-A caller's call sequence

## Why this exists

`convert_array_containers` deliberately does not do this work. Its Step 8
tells a Case-A caller to emit the straightforward per-call-site
`alloc`/`copy2F`/`free` block, even when a container's source is
loop-invariant grid metadata, and its Step 10 output lists these as
**hoisting candidates** instead of hoisting them on the spot -- because
hoisting correctly requires seeing every call site in the subroutine at
once, and a conversion pass only ever looks at one callee at a time.

The result, once a subroutine has had several callees converted
underneath it (the common case for a module entry point like
`continuity_PPM`, which calls out to `zonal_edge_thickness`,
`zonal_mass_flux`, `continuity_zonal_convergence`, and their meridional
mirrors), is a caller full of container churn that was each individually
correct but collectively wasteful: the same container freed and
reallocated from a source that never changed, or copied out to a raw
array and immediately copied into a *different* container from that same
raw array, purely to hand a value from one call to the next.

This skill is the deferred pass: given a subroutine whose container
conversions are done and stable, find every one of these patterns, and
either eliminate the operation entirely, replace a free+realloc with a
`%copy2Array` in place, or recognize that two containers are actually one
continuous value and can be merged into a single variable passed directly
into both callees' dummy slots. It also groups what survives so that the
subroutine's block of descendant calls (the actual physics) is
uninterrupted by marshalling, `alloc`ing before and `free`ing (and any
truly necessary copy-back) after.

This skill assumes `array_container_lessons` has already been invoked
this session, for the container API (`%alloc`/`%view`/`%copy2F`/
`%copy2Array`/`%free` signatures, in particular §4.4 and §4.6) and the
intent-mapping rules. It is not reproduced here.

## Scope

**In scope:** the marshalling code in one subroutine that already speaks
containers, and whose every descendant call already speaks containers
too -- reordering, merging, or eliminating `%alloc`, `%free`, `%copy2F`,
and `%copy2Array` calls, and the plain scalar/loop-invariant computations
that feed them (e.g. hoisting `edge_h_min = 2.0 * GV%Angstrom_H` out of a
branch that computes it twice).

**Hard precondition, checked in Step 0: every descendant this subroutine
calls must already be container-native.** If even one call site passes a
raw array to a callee that still takes raw arrays (a `zonal_BT_mass_flux`,
`meridional_BT_mass_flux`, or anything else `convert_array_containers`
hasn't reached yet), stop -- do not run on this subroutine. A still-raw
callee needs its raw array populated (via `%copy2F`) immediately before
the call that consumes it, which pins that copy2F in the middle of the
call sequence and makes the "uninterrupted descendant-call block" goal
(§4) permanently unreachable for that call, no matter how the rest of the
subroutine's containers are hoisted. Convert the raw callee first
(`convert_array_containers`), then come back.

**Second hard precondition, also checked in Step 0: the subroutine's own
dummy list must have no `optional` raw array.** A subroutine like
`continuity_adjust_vel`, whose own `optional` dummy (`visc_rem_u`) is
still a raw array, is forced to build its container for that argument
conditionally (`if (present(visc_rem_u)) then ... alloc ... free ...
endif`, lessons §6) -- and that conditional alloc/free can never be
hoisted to the unconditional top/bottom blocks, because doing so would
make the container always "present" to the callee regardless of whether
the caller actually supplied `visc_rem_u` (lessons §9 #8b). Convert that
optional dummy first with `convert_optional_args_to_containers`. This
usually does more than relocate the problem: once the subroutine's own
`visc_rem_u` is `visc_rem_u_a` (an optional container), it forwards
straight through to the callee's own optional `visc_rem_u_a` dummy with
no `present()` check at all -- the same-type optional-to-optional
forwarding rule preserves absence for free -- so the branch this skill
couldn't hoist typically disappears entirely rather than surviving as an
unavoidable exception.

**Explicitly out of scope:**
- Converting any raw array dummy or local to a container. That is
  `convert_array_containers`'s job.
- Touching any callee's signature or body. A callee is only ever a call
  target here.
- Touching a `%view`-guarded optional-forwarding block (lessons §9 #16)
  beyond what this skill's own hoisting requires -- its `nullify` and
  guard structure stays exactly as correct as it already is.
- Rewriting loop bodies, `do concurrent` bounds, or numerical logic of
  any kind. A hoist is numerically inert by construction (lessons §1);
  any difference afterward is a bug.
- Moving `cpu_clock_begin`/`cpu_clock_end` boundaries.
- Adding a new `cpu_clock` pair, or any other timing instrumentation,
  unless explicitly asked for separately. This skill's "grouping" step is
  about source-level readability and manual profiling, not about adding
  measurement code.

## The core technique -- one litmus test, applied three ways

Every eliminable pattern in this class reduces to a single question,
asked at each `%free()` that is followed (anywhere later in the same
subroutine, not just the next line) by a re-`%alloc()` of the same
container, and at each `%copy2F(X)` that is followed by a *different*
container's `%alloc(source=X)`:

> **Does anything else read or write the raw array `X` between these two
> operations?**

Walk the subroutine's actual body to answer this -- never assume from the
variable's name or its position in the source. If the answer is no, one
of three fixes applies:

### 1. Source never changes -- hoist to a single alloc/free

If a container's source (`G%mask2dT`, a mandatory dummy like `u` or `v`
that this subroutine never mutates, `hin` if the subroutine never writes
it, ...) is the same at every one of its use sites, there is no reason to
free and reallocate it between them at all. Allocate it once, at the top
of the subroutine (or the top of whichever branch uses it, if it is used
in only one branch and both branches are mutually exclusive at runtime --
check this explicitly, e.g. a `logical` set once near the top and never
reassigned), pass the same container to every call site that needs it,
and free it once, at the bottom.

```fortran
! BEFORE -- freed and rebuilt from the same unchanging source, twice
call mask2dT_a%alloc(lb=LBOUND(G%mask2dT), ub=UBOUND(G%mask2dT), source=G%mask2dT)
call zonal_edge_thickness(bxC, h_in_a, h_W_a, h_E_a, mask2dT_a, ...)
call mask2dT_a%free()
...
call mask2dT_a%alloc(lb=LBOUND(G%mask2dT), ub=UBOUND(G%mask2dT), source=G%mask2dT)
call meridional_edge_thickness(bxC, h_in_a, h_S_a, h_N_a, mask2dT_a, ...)
call mask2dT_a%free()

! AFTER -- one alloc, one free, for the whole subroutine
call mask2dT_a%alloc(lb=LBOUND(G%mask2dT), ub=UBOUND(G%mask2dT), source=G%mask2dT)
...
call zonal_edge_thickness(bxC, h_in_a, h_W_a, h_E_a, mask2dT_a, ...)
...
call meridional_edge_thickness(bxC, h_in_a, h_S_a, h_N_a, mask2dT_a, ...)
...
call mask2dT_a%free()
```

A container written by a call (`intent(inout)`) can be hoisted the same
way, as long as nothing *else* reads or writes its raw source in between
its uses -- see §2 and §3 below for what to do when it does need to
change.

### 2. Source genuinely changes once -- `%copy2Array` in place, not free+realloc

If a container's value must change partway through the subroutine (a
later call needs the *updated* value of something an earlier call wrote
back to a raw array), that is not a reason to `%free()` and `%alloc()`
again -- `%copy2Array` (lessons §4.4) overwrites an *already-allocated*
container's payload from a raw array of the same shape, with no
deallocate/reallocate cycle:

```fortran
! BEFORE -- free, then realloc from a source that DID change once
call h_in_a%free()
call h_in_a%alloc(lb=LBOUND(h), ub=UBOUND(h), source=h)

! AFTER -- same container, updated in place
call h_in_a%copy2Array(h)
```

This requires the container's shape to be identical across every use --
true whenever the two source arrays share the same declared dimensions
(e.g. `hin` and `h` are both `(SZI_(G),SZJ_(G),SZK_(GV))` in
`MOM_continuity_PPM.F90`). Confirm this from the declarations, not by
assumption.

### 3. Two containers are one continuous value -- merge them

The deepest form of this pattern, and the easiest to miss: if container
`A` is written by call 1 and then, soon after, container `B` is built
fresh from `A`'s own copied-out value purely to feed call 2 -- ask
whether call 2 could simply take `A` directly. **Fortran does not require
an actual argument's name to match the callee's dummy name; only type,
rank, and shape have to match.** If they do, there is no reason for `B`
to exist as a separate container at all for that use:

```fortran
! BEFORE -- A written, copied out to raw h, a DIFFERENT container B
! rebuilt from that same raw h, purely to feed the next call
call continuity_zonal_convergence(bxC, h_a, uh_a, dt, IareaT_a, hin_a=hin_a)
call h_a%copy2F(h)
call h_in_a%copy2Array(h)
call meridional_edge_thickness(bxC, h_in_a, h_S_a, h_N_a, mask2dT_a, ...)

! AFTER -- h_a IS the value meridional_edge_thickness needs as its h_in;
! pass it directly, positionally, into that slot. No round trip at all.
call continuity_zonal_convergence(bxC, h_a, uh_a, dt, IareaT_a, hin_a=hin_a)
call meridional_edge_thickness(bxC, h_a, h_S_a, h_N_a, mask2dT_a, ...)
```

`B` (`h_in_a` above) may still be needed under its own name for a
*different* use elsewhere in the subroutine (e.g. the first call in the
sequence, sourced from a genuinely different raw array such as `hin`
rather than `h`) -- merging only removes the uses where `A` and `B`
would otherwise hold byte-identical data. Do not delete `B`'s declaration
or its other, still-needed alloc/free if this is the case.

This fix requires tracing the actual data lineage through the
subroutine's body, not pattern-matching on variable names -- `h_a` and
`h_in_a` looking similar is a coincidence of this file's naming
convention (lessons §7), not a signal either way.

## Steps

### 0. Validate inputs

`$0` = work-directory, `$1` = function-name.

1. If `$0` or `$1` is empty, or `$1` is `help`/`--help`/`-h` -> print
   usage (mirroring `convert_array_containers --help`'s format) and stop.
2. If `$0` is not an existing directory -> stop:
   `Error: work directory "$0" does not exist.`
3. Parse `--enable_git_commit` / `--disable_git_commit`; stop if both are
   passed (mutually exclusive), same wording as `convert_array_containers`.
4. If `array_container_lessons` has not been invoked this session (no
   memory of §4.4 `copy2Array`, §4.6, or the naming conventions in §7),
   invoke it now before proceeding.
5. Locate `$1`'s declaration under `$0/src` and `$0/config_src`. If it has
   no `type(RealArray_t)`/`type(IntArray_t)` dummies or locals at all, stop
   and say so -- this skill has nothing to do on a still-raw subroutine;
   that is `convert_array_containers`'s job.
6. **Every descendant must already be container-native.** Find every
   `call <callee>(...)` inside `$1` (same grep as
   `convert_array_containers` Step 3) and check each callee's own
   declaration. If any actual argument passed to any callee is a raw,
   grid-shaped array (not a `type(RealArray_t)`/`type(IntArray_t)`
   container) -> stop:
   `Error: $1 calls <callee> with a raw array argument (<arg>). hoist_container_marshalling requires every descendant to already be container-native -- run convert_array_containers on <callee> first.`
   This is a hard requirement, not a judgment call: a still-raw callee
   pins a `%copy2F` immediately before the call that feeds it (§ Scope),
   which this skill cannot hoist or defer around, and mixing "fully
   hoistable" containers with "must stay pinned" ones in the same pass is
   exactly the confusion this precondition exists to avoid.
7. **The subroutine's own dummy list must have no `optional` raw array.**
   Read `$1`'s own declaration. If any `optional` dummy is still a raw,
   grid-shaped array (not a `type(RealArray_t)`/`type(IntArray_t)`
   container) -> stop:
   `Error: $1 has an optional raw array dummy (<arg>). hoist_container_marshalling requires every optional dummy to already be a container -- run convert_optional_args_to_containers on <arg> first.`
   Same reasoning as item 6: an optional raw dummy forces a runtime
   `if (present(...))`-guarded alloc/free (lessons §6) that can never be
   hoisted to the unconditional top/bottom blocks without corrupting the
   callee's own `present()` check (lessons §9 #8b). Converting it first
   usually removes the branch entirely rather than just relocating it --
   see § Scope for why.
8. Read the full subroutine. Build the table described in Step 1 below
   *before* changing anything, and if it shows every container already at
   exactly one `%alloc` and one `%free`, stop and report there is nothing
   to hoist.

### 1. Build the "before" table

For every container variable declared in the subroutine (`type(RealArray_t)
::` / `type(IntArray_t) ::` locals), walk the body in program order and
record: each `%alloc` call and its `source=` expression, each `%free`,
each `%copy2F`/`%copy2Array` and its target/source expression, and which
call(s) it is passed to in between. A short script (grep the subroutine's
extracted text for `%alloc(`, `%free()`, `%copy2F(`, `%copy2Array(`,
tally by container name) is the fastest way to build this and to verify
it afterward -- there is no need to read the whole body by eye once the
table exists.

This table is what makes every later step mechanical rather than
judgment-based: once you know container `X` is alloc'd 4 times, all with
`source=Y`, the fix in §1 above applies without re-reading the
surrounding code each time.

### 2. Classify every branch that matters

If the subroutine has an `if`/`else` (or `select case`) that runs exactly
one of several mutually-exclusive paths, determined by a value fixed
before any of the paths execute (e.g. `x_first = (MOD(G%first_direction,2)
== 0)` decided once near the top, never reassigned) -- confirm this
explicitly, then treat the two paths as symmetric for hoisting purposes:
a container used in exactly one of the two paths, with a provably
unchanging source, is still hoistable to a single subroutine-wide
alloc/free, because only one path ever runs per call. Do not hoist across
a branch boundary if the branch condition could, in principle, be
re-evaluated or the paths are not truly mutually exclusive -- check the
actual control flow, not just the presence of an `if`.

### 3. Apply the three fixes from the table

For each container in the Step 1 table, in order:
- **Unchanging source, multiple alloc/free pairs** -> hoist per §1 above.
- **Source changes exactly once, no other container built from the same
  raw array afterward** -> `%copy2Array` in place per §2 above.
- **A second container is built purely to receive what the first one
  just wrote out** -> merge per §3 above; trace the actual data lineage,
  do not assume from names.

After every fix, a container's remaining `%alloc`/`%free` should each
appear exactly once, unless it is genuinely used in two disjoint,
mutually-exclusive-at-runtime scopes with different shapes or you have a
specific, stated reason it cannot be merged (e.g. two callees expect
incompatible shapes for what looks like "the same" quantity) -- state that
reason in the Step 5 report rather than leaving it unexplained.

### 4. Group the survivors

Move every surviving `%alloc` to a single block before the subroutine's
main branch (or before the sequence of descendant calls, if there is no
branch). Move every surviving `%free`, plus any `%copy2F` that is not
itself required *before* a later step in the same call sequence (i.e. it
has no reader between where it currently sits and the end of the
subroutine), to a single block after that sequence, in this fixed order
within the block: all `%copy2F` calls first, then all `%free` calls.

A `%copy2F` that *does* feed a later `%copy2Array`/`%alloc(source=...)` in
the same sequence -- the Step 3 §3 merge case, before the merge is fully
applied, or any pairing a merge turned out not to be possible for --
stays exactly where it is needed; do not move a copy-back earlier or
later than the point that makes it correct. Since Step 0 already refused
to run on any subroutine with a still-raw descendant, this is the *only*
reason a `%copy2F` should remain mid-sequence -- there is no other
legitimate case for one, and the goal below is always fully achievable.

The goal is that the subroutine's block of descendant calls -- the
`bxC = set_continuity_box(...)` / `call zonal_edge_thickness(...)` /
`call zonal_mass_flux(...)` / `call continuity_zonal_convergence(...)`
style sequence -- reads as pure computation, with no `%alloc`, `%free`,
`%copy2F`, or `%copy2Array` calls interleaved, once this step is done.

### 5. Verify

- Rebuild the Step 1 table against the edited code. Every container's
  `%alloc` count must equal its `%free` count.
- For every remaining `%copy2F`/`%copy2Array` call, state in one line why
  it is necessary (feeds a specific later call or dummy write-back) --
  an unexplained one is a sign a fix from Step 3 was missed.
- Diff review: no edits inside any loop body; no callee's signature or
  body touched; only the calling subroutine's marshalling and branch
  structure changed.
- Line length <= 100 on every added line (lessons §7, same check as
  `convert_array_containers` Step 9).
- No doc comment (`!<`) added, removed, or reworded -- this skill never
  touches a dummy declaration, so the `!<` count before and after must be
  identical, not just balanced.
- **Build:** do not assume a Fortran compiler is available on this
  machine (this project's environment routinely has none); if one is
  present, build under both infra layers; if not, say plainly that the
  build was not run. A hoist is numerically inert by construction, but
  that claim is only as good as the trace in Step 3 -- flag it as worth a
  real build once one is available, same as any other container change.

## Hard rules

- Never run on a subroutine that calls even one still-raw descendant
  with a raw array argument (Step 0 item 6). Convert that descendant
  first, or pick a different target.
- Never run on a subroutine that itself still has an `optional` raw
  array dummy (Step 0 item 7). Run `convert_optional_args_to_containers`
  on it first.
- Never convert a raw array to a container here -- that is
  `convert_array_containers`'s job, run separately if something in scope
  turns out to still be raw.
- Never touch a callee's signature or body.
- Never move a `cpu_clock_begin`/`cpu_clock_end` boundary, and never add
  a new clock or other timing instrumentation unless explicitly asked.
- Never hoist a container's alloc/free across a point where something
  else -- not itself -- reads or writes its raw source, without first
  confirming from the actual body that no such reader/writer exists.
  "It doesn't look like anything else uses it" is not a substitute for
  checking.
- Never merge two containers (§3) based on their names looking related
  -- trace the value lineage through the body first.
- Never leave a `%copy2F`/`%copy2Array` call whose purpose cannot be
  stated in one sentence; either justify it in the Step 5 report or
  remove it.
- Do not add comments about the hoist itself -- no `! Hoisted to avoid
  re-allocating`, no `! Merged with h_a -- see chat`, anywhere in the
  code. This codebase documents ocean physics, not its own refactoring
  history (lessons §7, "No commentary about the conversion itself"). The
  Step 5/Output report is where this reasoning belongs, not the source.
- Do not attempt to install a Fortran compiler, and do not claim a build
  passed that was never run.

## Commit gating

Same as `convert_array_containers`: `--enable_git_commit` /
`--disable_git_commit` override; otherwise follow
`~/.claude/preferences.json`'s `git_commit_and_push` key (`"auto"` runs
the commit; `"manual"`, a missing file, a missing key, or unparseable
JSON all mean skip it and report modified files for manual commit).
Branch name: `claude_<lowercased_function-name>_hoist`.

## Output to the user on success

Report:
1. The subroutine, and the before/after `%alloc` count total (sum across
   all containers) as the headline number.
2. Per container: before count -> after count, and which of the three
   fixes (§1/§2/§3) applied, in one line each.
3. Any container left with more than one alloc/free pair, and why.
4. Confirmation that the descendant-call sequence is now uninterrupted by
   marshalling (Step 0's precondition guarantees this is always fully
   achievable; if a `%copy2F` still had to stay mid-sequence, say what it
   feeds -- this should only ever be a Step 3 §3 merge case, never a
   still-raw callee).
5. Doc-comment count unchanged, line-length check, diff-scope confirmation
   (only this subroutine's body changed).
6. Build status.
7. Whether committed, or the list of modified files for manual commit.
