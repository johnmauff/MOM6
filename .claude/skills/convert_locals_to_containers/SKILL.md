---
name: convert_locals_to_containers
description: Convert a subroutine's own internal scratch-local arrays (never its dummy arguments) to RealArray_t/IntArray_t containers, wherever every callee that consumes a given local already accepts a container for that argument -- eliminating the alloc-source-copy-free round trip currently needed to hand a freshly-computed local off to an already-converted callee. Accepts a comma-separated list of subroutine names so a whole call tree (or however much of it has had convert_array_containers run on its interfaces already) can be swept in one invocation, each subroutine verified independently. Distinct from convert_array_containers (which only ever touches dummy arguments, never locals) and from hoist_container_marshalling (which only reorders/merges alloc, free and copy2F calls for containers that already exist -- it never changes a variable's fundamental type). Use this once a callee a subroutine calls has already been converted and that subroutine still computes the matching argument into a raw local purely to copy it in and immediately throw it away.
user-invocable: true
argument-hint: <work-directory> <function-name>[,<function-name>...] [--enable_git_commit] [--disable_git_commit]
---

# Convert a subroutine's scratch locals to containers, not just its interface

## Why this exists

`convert_array_containers` converts a subroutine's *dummy* arguments and
stops there -- it explicitly never looks at local variables. Once a
callee has been converted, though, a caller that computes that callee's
argument into a plain local (rather than receiving it as one of its own
dummies) is stuck doing `alloc(source=local)` / call / `free` at every
call site, purely to hand off a value that was just computed and will
never be read again in raw form -- a real, observed pattern in
`zonal_mass_flux`/`meridional_mass_flux`, whose scratch locals (`du`,
`dv`, the CFL-limit arrays, the resolved viscosity-remnant copy, the
summed-transport accumulators) get exactly this treatment at their calls
into `present_uhbt_or_set_BT_cont`/`present_vhbt_or_set_BT_cont`.

This is a third, distinct kind of work from the other two skills in this
family:

- `convert_array_containers` decides raw-vs-container for **dummy
  arguments** and rewrites the subroutine's public interface.
- `hoist_container_marshalling` reduces `%alloc`/`%free`/`%copy2F`/
  `%copy2Array` **churn for containers that already exist** -- it never
  changes what type a variable is.
- This skill decides raw-vs-container for a **local, non-dummy
  variable**, specifically because retyping it lets one or more of its
  consuming call sites skip the alloc-and-copy-in step entirely. It
  doesn't touch the subroutine's own signature at all.

Composing all three, the natural order down a call tree is:
`convert_array_containers` (make the interface a container) ->
`convert_optional_args_to_containers` (narrower interface work, if the
subroutine has optional arguments worth pushing up) -> **this skill**
(stop rebuilding containers from scratch for values the subroutine
already computed itself) -> `hoist_container_marshalling` (group and
minimize whatever alloc/free calls survive). Running this skill can
create new containers for `hoist_container_marshalling` to then group;
running it does **not**, however, change the subroutine's own dummy
list, so it never makes that subroutine newly eligible as *its own*
caller's target -- the effect is strictly downward, one level at a time,
same as the other two.

This skill assumes `array_container_lessons` has already been invoked
this session, for the container API (`%alloc`/`%view`/`%free`
signatures, lessons §4) and the naming convention (lessons §7). It is
not reproduced here.

## Scope

**In scope:** local array variables declared inside one or more named
subroutines (never a dummy argument), where at least one callee they are
passed to already takes a container for that exact argument.

**The eligibility test is per-local, not per-subroutine** -- the same
granularity `hoist_container_marshalling` uses for its own per-container
entanglement check, just applied to a different question (whether to
convert a local at all, versus whether an existing container's
marshalling can be hoisted). A subroutine can have some locals convert
and others stay raw in the very same pass -- eligibility depends only on
where *that particular* local is used, not on the state of the rest of
the subroutine. See "The core technique" below for the exact test,
including the case where a local feeds *both* a converted and a still-raw
callee (common, and not a blocker).

**Explicitly out of scope:**
- Converting a subroutine's own dummy arguments. That is
  `convert_array_containers`'s job; this skill assumes it has already
  been run on whatever interfaces it depends on.
- Reordering, merging, or grouping `%alloc`/`%free` calls for containers
  that already exist independent of this specific local-to-container
  promotion. Run `hoist_container_marshalling` after this skill for that
  cleanup -- it's a natural next step, not part of this one.
- Touching any callee's signature or body.
- Rewriting loop bodies or changing any numerical result -- exactly as
  inert as the other two skills; the `%view` pointer keeps the local's
  original name so the body needs zero edits beyond the declaration
  itself.
- Silently reinterpreting an `!$omp target`/`map(...)` directive that
  names a candidate local. Whether a `%view`-derived pointer maps to a
  device the same way a plain automatic array does is compiler- and
  offload-target-specific and cannot be verified by inspection alone on
  a machine with no compiler (lessons §10, same standing caveat as every
  other skill in this family) -- see Hard rules.

## The core technique

For every local array declared in a target subroutine (skip anything
that's a dummy argument -- that's out of scope entirely):

1. **Never called anywhere -> not in scope.** If the local is pure
   internal scratch that never appears as an actual argument to any
   `call`, there is nothing to gain from containerizing it. Leave it
   alone.

2. **Called into at least one already-converted callee -> eligible.**
   Check every call site the local is passed to, one at a time, against
   that callee's *current* dummy declaration at the matching position:
   - If the callee's dummy there is already `type(RealArray_t)`/
     `type(IntArray_t)`, this call site currently needs the
     alloc-source-copy-free round trip -- convert the local so this call
     site can instead pass the container directly, with nothing built
     fresh.
   - If the callee's dummy there is still a raw array, this call site
     already just passes the raw local (or a `%view` pointer, if the
     *caller's own* dummy that fed this local was itself already
     converted) -- nothing changes here, and this does **not** block
     conversion of the local elsewhere in the subroutine.
   A local converts as long as *at least one* of its call sites benefits;
   it is completely normal for the same local to be passed as a
   container to one callee and as a `%view` pointer to another,
   still-raw one in the same subroutine -- this needs no special
   handling beyond the ordinary guard-free `%view` pattern (lessons §4.3),
   since a local is a plain mandatory value, never `optional`, so none of
   the disassociated-pointer/`nullify` machinery from lessons §9 #16
   applies here at all.

3. **Referenced by an `!$omp target`/`map(...)`/similar directive by
   name -> stop and ask, don't guess.** Converting the local's storage
   from a plain automatic array to a container's backing store changes
   what that directive is actually mapping, and this skill cannot verify
   on a machine with no compiler whether a `%view`-derived pointer
   behaves identically under an existing device-data-mapping directive.
   Flag every such local explicitly and let the user decide, the same
   way `convert_array_containers` Step 2 stops and asks about every
   `optional` array dummy rather than assuming an answer. **This is not
   a hypothetical edge case for the flagship use of this skill** -- the
   locals named above (`zonal_mass_flux`/`meridional_mass_flux`'s `du`,
   `dv`, ...) are all currently named in exactly this kind of
   `!$omp target enter data map(alloc: ...)` directive, so expect this
   question on the very first real run, not as a rare exception.

4. **Otherwise -> convert.** Rename the local's declaration `<name>` to
   `<name>_a` (container, same rank), preserving its existing plain `!`
   comment verbatim if it has one -- locals in this codebase carry plain
   trailing comments, not the `!<` Doxygen form required on dummies, so
   the doc-comment check here is "the comment text survives," not a
   `!<` count (lessons §7 still applies to the *wording*, just not to
   the marker). Declare a pointer with the local's *original* name and
   matching rank, `%view` it once, and the rest of the subroutine's body
   -- every loop, every element reference -- needs zero edits, exactly
   like converting a dummy (lessons §4.3, §7).

   **Allocate with shape only, no `source=`.** Unlike a call-site
   container built from a real, externally-meaningful array (lessons §6,
   which defaults to copy-in precisely because the caller's raw array
   already held real data worth preserving), a local that is being
   containerized *is* the only home this data has ever had -- there is no
   separate "old" array whose contents would otherwise be lost. Confirm
   the local is genuinely write-before-read within this same invocation
   (true of ordinary Fortran locals by default, since they carry no
   initial value either) before skipping `source=`; if you find a
   genuine read-before-write on some path, treat it like any other
   `intent(inout)`-shaped case and copy in instead (lessons §6).

5. **Delete the old marshalling at every call site that now receives the
   container directly.** The `%alloc(source=...)` / call / `%free()`
   block that used to build a *separate*, throwaway container from this
   local disappears entirely -- the local's own new container is passed.

## Steps

### 0. Validate inputs

`$0` = work-directory, `$1` = comma-separated function-name list.

1. Help/empty-argument handling, work-directory existence, and
   `--enable_git_commit`/`--disable_git_commit` parsing, identical to
   `convert_array_containers` Step 0 items 1-3.
2. If `array_container_lessons` has not been invoked this session,
   invoke it now.
3. Split `$1` on commas. For each name: locate its declaration under
   `$0/src` and `$0/config_src` (same grep as `convert_array_containers`
   Step 0 item 6); stop on 0 or >1 matches for any single name, same
   wording. A target does **not** need to already be container-native
   itself -- a fully raw subroutine can still have a local that happens
   to feed an already-converted callee, and that local is still in
   scope.
4. For each name, if it has no local array declarations at all besides
   its dummies, drop it from the list and report "nothing to do" for
   that one; do not error the whole run over it.
5. **Name-collision check.** Before renaming any local `<name>` to
   `<name>_a`, confirm no dummy or other local already named `<name>_a`
   exists in that subroutine. If one does, stop and ask -- do not
   silently pick a different suffix.

### 1. Build the "before" table

Per target subroutine, per local array: its rank, every call site it
feeds (callee name, argument position), whether that position is
currently a container in the callee's live declaration, and any
`!$omp`/similar directive that names it. Classify each local as
**eligible**, **not called anywhere** (skip), or **flagged** (an OMP
directive names it -- ask before proceeding).

### 2. Convert each eligible local

Per "The core technique" step 4: rename the declaration, add the `%view`
pointer under the original name, alloc with shape only (verify
write-before-read first). Do this for every eligible local in every
target subroutine before touching any call site, so the full set of
renames is known when Step 3 runs.

### 3. Update every call site the local feeds

For each (eligible local) x (call site) pair: if that callee's dummy is
already a container, delete the old alloc/copy-in/free block and pass
`<name>_a` directly; if the callee's dummy is still raw, leave that call
site passing the `%view` pointer `<name>`, unchanged.

### 4. Verify

- Rebuild the Step 1 table against the edited code: every eligible local
  now has exactly one `%alloc` and one `%free`; every call site into an
  already-container-native callee passes the container, not a rebuilt
  one.
- Plain-comment text preserved verbatim for every renamed local (not a
  `!<` count -- see "The core technique" step 4).
- Line length <= 100 on every added line (lessons §7).
- Diff review: no edits inside any loop body; no callee's signature or
  body touched; only the target subroutines' own local declarations and
  call-site marshalling changed.
- Every flagged (OMP-directive) local reported with the user's decision,
  not silently resolved either way.
- **Build:** do not assume a compiler is available; if one is present,
  build under both infra layers; if not, say so plainly and flag the
  OMP-directive interaction specifically as the highest-value thing to
  confirm once a build is possible -- more so than the ordinary
  disassociated-pointer caveat the other skills carry, since this is a
  genuinely new mechanism (a `%view` pointer standing in for a plain
  automatic array under an existing device-mapping directive).

## Hard rules

- Never convert a local that is never passed to any call -- no benefit,
  don't add the complexity.
- Never let one still-raw callee block conversion of a local that also
  feeds an already-converted one -- convert it, and leave the still-raw
  call site passing the plain `%view` pointer, unchanged.
- Never silently convert a local named in an `!$omp target`/`map(...)`/
  similar directive. Stop and ask every time, even if this is the
  fifth local in a row that hits this -- it is not automatable on a
  machine without a compiler to verify the interaction.
- Never add a `source=` copy-in for a newly-containerized local unless
  you have confirmed a genuine read-before-write on some path -- default
  to shape-only allocation, the opposite default from call-site
  container marshalling (lessons §6), and state why in the Step 4
  report.
- Never touch a subroutine's own dummy argument declarations -- that is
  `convert_array_containers`'s job, run separately.
- Never touch a callee's signature or body.
- Never rewrite loop math; the `%view` pointer keeps the local's
  original name specifically so the body needs zero edits.
- Do not add comments about the conversion itself -- same standing rule
  as every other skill in this family (lessons §7, "No commentary about
  the conversion itself").
- Do not attempt to install a Fortran compiler, and do not claim a build
  passed that was never run.

## Commit gating

Same as the other skills in this family: `--enable_git_commit` /
`--disable_git_commit` override; otherwise follow
`~/.claude/preferences.json`'s `git_commit_and_push` key. Branch name:
`claude_<lowercased_first_function_name>_local_containers` if one name
was given, or `claude_<lowercased_first_function_name>_and_N_more_local_containers`
for a multi-name run (`N` = the count of additional names).

## Output to the user on success

Report, per target subroutine, then a combined total:
1. Which locals were converted, and which call sites for each now pass
   the container directly versus still pass a `%view` pointer to a
   still-raw callee.
2. Which locals were left alone and why (never called, or flagged for
   an OMP directive -- with the user's decision recorded for the
   flagged ones).
3. Plain-comment preservation, line-length, and diff-scope confirmation.
4. Build status, with the OMP-directive caveat called out explicitly if
   any local was converted under a directive the user approved.
5. Whether committed, or the list of modified files for manual commit.
6. A one-line pointer to `hoist_container_marshalling` as the natural
   next step if any new containers were created.
