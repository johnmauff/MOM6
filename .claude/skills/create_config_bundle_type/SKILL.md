---
name: create_config_bundle_type
version: "0.3"
description: Cluster a private (module-local) control-structure/config derived type's fields — scalars, and already-containerized arrays — into small, purpose-built bind(C)-ready bundle types, nested back into the original type, based on which fields actually travel together across the call tree. Worked example in this codebase: continuity_PPM_CS's 20 individually-flattened fields were clustered into reconstruction_opts_type (upwind_1st/monotonic/simple_2nd) and transport_adjust_opts_type (CFL_limit_adjust/aggress_adjust/vol_CFL/better_iter/use_visc_rem_max/marginal_faces/tol_eta/tol_vel). Use this instead of leaving convert_array_containers's flattening as the end state, and instead of create_shadow_container_type, whenever the type in question is private to the module/call tree being worked on (no shadow-and-copy-back dance needed) rather than shared with code outside it.
user-invocable: true
argument-hint: <work-directory> <type-name> [<call-tree-entry-point>] [--enable_git_commit] [--disable_git_commit]
---

# Cluster a private control structure's fields into bind(C)-ready bundles

## Why this exists

`convert_array_containers`'s argument-list-reduction step (its own
"Settle these decisions" item 2) flattens a private control structure's
fields into individual dummies wherever a leaf kernel needs them — this
is necessary groundwork whenever the type itself can never be `bind(C)`
(e.g. it holds a `type(diag_ctrl), pointer`), so a leaf destined for a
bridge can't take the whole struct. Left at that end state, the same
handful of fields end up reappearing individually, repeatedly, across
every signature and every `bind(C)` bridge interface down the call
tree — the flattening obscures exactly the logical relationships between
arguments that the original struct expressed for free.

**The fix:** trace which fields actually travel together across the
call tree, and cluster each such group into its own small bundle type —
plain Fortran fields, a `bind(C)` mirror, and a converter — nested back
into the original struct in place of the flat fields. Threading one
bundle dummy through a dozen signatures reads the same as threading the
original struct did, without reintroducing the struct's own
bridge-incompatibility (the bundle holds only bridge-representable
fields; the un-bindable ones, like the `diag_ctrl` pointer, stay on the
original struct, untouched).

**Not `create_shadow_container_type`.** That skill solves a different
problem: a type *shared* with code outside the call tree being worked
on, where you can't touch the real fields at all and must build a
temporary parallel shadow, copy in, forward, copy back, free. A private,
module-local control structure has no such blast radius — restructure
its real fields in place, once, with no shadow and no copy-back. If the
type you're about to cluster turns out to be used outside this call
tree (grep `use <its module>, only : <TypeName>` repo-wide), stop and
use `create_shadow_container_type` instead.

**Runs on an already-flattened struct or an unflattened one.** Starting
unflattened is mechanically simpler (straight to clustered bundles, no
flattened intermediate state), but usage-tracing costs the same either
way, and every leaf still needs Step 6's explicit per-leaf call — keep
the whole struct, or stop taking it (Step 7) — never both.

This skill assumes `array_container_lessons` (container API,
`%to_c()`/`bind(C)` conventions) and, if a bridge already exists
anywhere in the call tree, `cpp_bridge_lessons`, have already been
invoked this session. Neither is reproduced here.

## Scope — field eligibility

Walk every field of the target type and classify it before clustering
anything:

- **Scalar** (any intrinsic type) — always eligible, no precondition.
- **Array** — eligible *only if it is already a container*
  (`RealArray_t`/`IntArray_t`/`LogicalArray_t`). Never convert a raw
  array as part of this skill — that is `convert_array_containers`'s
  job. If a field you want to cluster is still a raw array, stop and run
  `convert_array_containers` on the subroutines that need it first,
  then come back.
- **Pointer to another derived type, or anything else that can never be
  `bind(C)`-representable** (e.g. `continuity_PPM_CS`'s
  `type(diag_ctrl), pointer :: diag`) — never bundle. Leave it exactly
  where it is on the original struct, alongside any other field that
  ends up in no cluster at all.

A cluster may freely mix eligible scalars and eligible container fields
— the converter (below) handles each member's own kind. This is not a
new mechanism invented for this skill: `reconstruction_opts_to_c`/
`transport_adjust_opts_to_c` already show the pure-scalar case
(unconditional, direct assignment per field) and `BT_cont_container_to_c`
already shows the pure-container case (`#ifdef _TIM`-guarded, delegates
to each field's own `%to_c()`); a mixed cluster's converter is just both
rules applied per field in the same function.

## Hard precondition

Confirm the target type is genuinely private to the call tree you're
about to touch: grep `use <module>, only : <TypeName>` (and any
qualified `<module>%<TypeName>` reference) across the whole repo, not
just the module that defines it. Any hit outside the call tree means
this skill does not apply — use `create_shadow_container_type` instead.
This is exactly the check that made `continuity_PPM_CS` eligible here
while `BT_cont_type` (used directly by `MOM_barotropic.F90` and the
RK2/RK2b dynamics files) was not.

## Settle these decisions

1. **Bundle uniformly, with one documented exception.** A subroutine
   that needs only 1–2 fields of a cluster still takes the whole bundle,
   rather than inventing a narrower type. The one exception: a leaf in
   the hot inner-loop path (`elemental`, called per grid cell) may keep
   a bare scalar, unpacked from the bundle by its caller at that one
   call site — mirror `flux_elem`'s treatment of `vol_CFL` in the
   `continuity_PPM_CS` work. Call out any such leaf explicitly; don't
   apply the exception silently or to more than the one hottest call.
2. **How many clusters** comes from tracing actual usage (Step 2 below),
   never from guessing by field name or proximity in the struct
   definition. Two disjoint clusters is common but not a rule.
3. **Don't bundle a field with no repeated companion.** If a field is
   read or forwarded at only one or two call sites and never travels
   with another field, leave it as a bare individual dummy. Bundling
   only pays off when a cluster of fields recurs together across many
   signatures — confirmed negatively in this codebase by `GV`
   (`verticalGrid_type`) and `US` (`unit_scale_type`) in the continuity
   call tree: each of their touched fields is consumed at 1–2 call sites
   with no recurring companion, so neither is a bundling candidate
   there despite both being flattened structs in the same call tree.
4. **Naming**: `<cluster>_opts_type` / `_C` / `_to_c`, matching the
   `reconstruction_opts`/`transport_adjust_opts` precedent, unless this
   codebase's existing naming for the domain suggests otherwise.
5. **Stage multi-cluster rollouts.** If more than one cluster is found,
   roll each out as its own atomic pass (every call boundary for that
   cluster's fields must agree at once, so a partial rewrite doesn't
   compile) — mirror the Stage A/Stage B split used for
   `reconstruction_opts_type`/`transport_adjust_opts_type`.

If the user already specified any of these, take their values as-is.

## Steps

### 0. Validate inputs

`<work-directory>` and `<type-name>` are required; `<call-tree-entry-point>`
defaults to every public subroutine in the type's own module. If either
required argument is empty, or equals `help`/`--help`/`-h`, print usage
(mirroring this file's frontmatter) and stop. If `<work-directory>` is
not an existing directory, stop. Locate `<type-name>`'s declaration
(`grep -irn "^[[:space:]]*type.*::[[:space:]]*<type-name>\b"`); zero
matches is an error, more than one asks the user which to use.

### 1. Classify every field

Read the type's full declaration and classify each field per Scope
above. Stop and report if any field you'd otherwise want to cluster is
still a raw array (refer to `convert_array_containers`); do not convert
it yourself.

### 2. Trace usage across the call tree

For every subroutine reachable from `<call-tree-entry-point>` down to
every leaf, record which eligible fields it reads, writes, or forwards
(whether currently reached via individually flattened dummies, direct
`<dummy>%<field>` access, or the whole struct passed wholesale). This is
the expensive step — for anything beyond a handful of subroutines, use
an Explore/fork agent to build the usage map rather than guessing from
memory or from field names.

### 3. Cluster by co-occurrence

Group fields that are read or forwarded together across most or all of
the same subroutines. A field with no recurring companion is not
bundled (Settle-these-decisions item 3) — leave it exactly as found.

### 4. Define each cluster's bundle

For each cluster, add (in the same module as the original type):

```fortran
!> <One-line purpose of this cluster of options>.
type, public :: <cluster>_opts_type
  <field's own declared type> :: field1 !< <doc, copied from the original struct's field doc>
  ...
end type <cluster>_opts_type

!> bind(C) mirror of <cluster>_opts_type, field-for-field, same order.
type, bind(C) :: <cluster>_opts_C
  <c_kind or _C container-mirror type> :: field1
  ...
end type <cluster>_opts_C
```

Then the converter. If every member is a scalar, it is unconditional —
no infra dependency, matching `reconstruction_opts_to_c`/
`transport_adjust_opts_to_c`:

```fortran
!> Converts a <cluster>_opts_type to its bind(C) mirror.
function <cluster>_opts_to_c(opts) result(cdesc)
  type(<cluster>_opts_type), intent(in) :: opts
  type(<cluster>_opts_C) :: cdesc
  cdesc%field1 = opts%field1   ! scalar: direct assignment
  ...
end function <cluster>_opts_to_c
```

If the cluster has at least one container member, guard the whole
function with `#ifdef _TIM` (it now depends on `%to_c()`, which only
exists under the TIM infra layer) and delegate that member, matching
`BT_cont_container_to_c`:

```fortran
cdesc%field1 = opts%field1        ! scalar: direct assignment
cdesc%field2 = opts%field2_a%to_c() ! container: delegate
```

### 5. Restructure the original type

Replace each clustered field with one
`type(<cluster>_opts_type), public :: <cluster>_opts` member per
cluster. Leave every non-clustered field (pointers, lone scalars from
Settle-these-decisions item 3) exactly where it is.

### 6. Per leaf, decide once: keeps the whole struct, or stops taking it

Every subroutine that currently touches the target type's fields —
anywhere in the call tree, not just its own module — falls into exactly
one of these; a leaf cannot straddle both for the same type:

- **Keeps taking the whole struct.** True whenever the leaf has no
  reason to stop (no bind(C) bridge on it now or planned, and nothing
  else forcing individual-field access). This needs **no signature
  change at all** — only a body rename, `<dummy>%<field>` →
  `<dummy>%<cluster>_opts%<field>` for every clustered field it touches.
  This is Step 5's restructuring alone rippling through; it applies
  identically whether the leaf is the type's own init/accessor routine
  or an arbitrary subroutine three levels down that happened to still
  take the whole struct.
- **Stops taking the whole struct.** True for a leaf that needs (now or
  per its own bridging work) individual-field access instead — already
  flattened by `convert_array_containers`, or about to be. Its dummy
  list changes to exactly: this cluster's bundle (Step 7) for every
  cluster it uses fields from, plus any standalone (never-bundled)
  field it uses promoted to its own bare dummy via
  `convert_array_containers`'s Step 2b, in the same pass — never a mix
  of "still takes the whole struct" *and* "also takes a bundle
  extracted from that same struct," which would be redundant by
  construction (the bundle is nested inside the struct after Step 5).

Decide this per leaf before touching anything, so Step 7 below has an
unambiguous target signature for every leaf it touches.

### 7. Roll out one cluster at a time

Per the staging decision, for one cluster at a time, touch only the
leaves classified above as "stops taking the whole struct" that use
this cluster's fields: replace whatever currently carries them —
individually flattened dummies, or direct struct-field access if this
is the first cluster reaching that leaf — with the one bundle dummy,
alongside whatever other bundles/standalone fields that same leaf also
needs (Step 6). For an already-bridged leaf that means its
`_fortran`/shim/`bind(C)` triple together (all three must agree); for a
leaf with no bridge yet, just its own plain signature. Apply the
hot-loop exception (item 1) only at the one leaf it was scoped to.

Fix any naming inconsistency uncovered for free during the rewrite (a
dummy renamed away from the field's real name at some leaf) by
restoring the bundle field's real name at that call site.

## Verify

- Programmatic argument count/order check across every `_fortran`/shim/
  `bind(C)` triple touched (or plain signature, if unbridged) — a script
  comparing each, not a by-eye read.
- Every cluster's `_to_c` converter: scalar members assigned directly,
  container members delegate to `%to_c()`, and the function is
  `#ifdef _TIM`-guarded if and only if it has at least one container
  member.
- Zero remaining bare references to any bundled field name outside of
  `%field` access, anywhere in the touched call tree.
- Every field left un-bundled (Settle-these-decisions item 3) is named
  explicitly in the report, with the call-site count that justified
  leaving it alone.
- Every touched leaf is classified as exactly one of "keeps the whole
  struct" (body rename only) or "stops taking it" (Step 7 dummy-list
  change) — none mixes both for the same type.
- Line length ≤ 100; `#ifdef`/`#ifndef`/`#if` vs `#endif` depth-tracked
  across the whole file (must never go negative, must end at 0) rather
  than a raw token count.
- No Fortran compiler locally — report the build as unrun.

## Versioning marker

Every Fortran file this skill creates or modifies gets a `!!SKILLS: 0.3`
marker line — the shared version number for this whole skill family,
not just this one skill (bump every skill file's `version:` field and
this marker in lockstep when any of them changes in a way that affects
generated code, never per-skill). If the file doesn't already have one,
add it as its own line immediately after the file's existing license/
header comment block, before the `module` statement. If it already has
one, update it to the current version rather than adding a second line.
Deliberately grep-able (`grep -rn "!!SKILLS:"`) and meant to be stripped
later, once these markers are no longer useful.

## Hard rules

- Never skip the `!!SKILLS: 0.3` marker on a file this skill touches,
  and never add a second marker line if one already exists — update it
  in place instead.
- Never touch a type that's used outside the call tree being worked on
  — that's `create_shadow_container_type`'s job; check first.
- Never convert a raw array to a container here — `convert_array_containers`,
  run separately, first.
- Never bundle a pointer-to-derived-type or other field that can never
  be `bind(C)`-representable.
- Never guess a cluster from field names or declaration order — trace
  actual usage (Step 2) first.
- Never bundle a field with no recurring companion across the call
  tree; leave it as a bare individual dummy instead.
- Never apply the hot-loop bare-scalar exception beyond the one leaf it
  was scoped to.
- Never leave a leaf both still taking the whole struct *and* also
  taking a bundle dummy extracted from that same struct — decide once,
  per leaf (Step 6), which of the two it is.
- Do not add commentary about the refactor itself in the source — this
  codebase documents ocean physics, not its own refactoring history.
- Do not claim a build passed that was never run; do not attempt to
  install a Fortran compiler.

## Commit gating

Same as `convert_array_containers`: `--enable_git_commit`/
`--disable_git_commit` override; otherwise follow
`~/.claude/preferences.json`'s `git_commit_and_push` key (`"auto"` runs
the commit step; `"manual"`, a missing file, a missing key, or
unparseable JSON all mean skip it and report modified files for manual
commit). Branch name: `claude_<lowercased_type-name>_config_bundle`.

## Output to the user on success

1. The clusters found, and each one's field list.
2. Every field classified but left un-bundled, and why (no recurring
   companion, or excluded as never-bind(C)-representable).
3. Whether the hot-loop bare-scalar exception was applied, and where.
4. Which clusters have been rolled out vs remain, if staged.
5. Any naming inconsistencies fixed along the way.
6. Build status — explicitly unrun if no compiler is available.
7. Whether committed, or the list of modified files for manual commit.
