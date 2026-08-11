---
name: convert_optional_args_to_containers
version: "0.3"
description: Convert only the optional array arguments of a MOM6 subroutine that sits in the middle of a call tree to RealArray_t/IntArray_t containers, leaving its mandatory arguments and every still-raw child subroutine it calls completely untouched. Uses a guarded %view producing a disassociated pointer, forwarded unconditionally to the child's still-raw optional dummy, to avoid the combinatorial if(present(...)) branch trees that a naive Case-A conversion produces. Use this when a full convert_array_containers pass on the subroutine (or its neighbors) isn't wanted yet, or when the branching cost of converting an optional argument bottom-up has already bitten once and you want to fix it by moving the container boundary up past the argument's real point of use instead.
user-invocable: true
argument-hint: <work-directory> <function-name> <optional-dummy-name>[,<optional-dummy-name>...] [--enable_git_commit] [--disable_git_commit]
---

# Convert a subroutine's optional arguments to containers, without converting its children

## Why this exists

`convert_array_containers` converts every array dummy of one subroutine
— mandatory and optional alike — and assumes the callee ends up fully
container-based. That's wrong for a narrower, common case: a **pivot**
subroutine simply forwards a few optional array arguments, unmodified,
to still-raw **child** subroutines further down — and converting the
child first, bottom-up, produces a branch-tree mess at the pivot.

**Concretely:** converting a child's optional arguments while its pivot
caller stays raw forces the pivot to build a *fresh* container from its
own raw optional dummy and conditionally include it — but a
freshly-built local container isn't the same actual argument as the
pivot's own dummy, so Fortran's "absent optional forwarded to another
optional stays absent" rule no longer applies, and the pivot ends up
enumerating every present/absent combination by hand. A real case in
this codebase produced a 16-way branch tree at each of two call sites,
plus an unrelated 8-way tree inside the child, just to forward the same
values one level further down.

**The fix:** convert the *pivot's* optional dummies instead, and thread
containers through its external callers. The pivot can then forward
into its still-raw children with zero branching in either direction
(mechanism below). The children stay untouched; if they need containers
too, that's a separate `convert_array_containers` invocation.

This skill assumes `array_container_lessons` has already been invoked
this session, for the container API and intent-mapping rules.

## The mechanism (read this before doing anything)

Two independent directions of data flow through the pivot, each with
its own branch-free trick. Get either wrong and you either reintroduce
combinatorial branching or corrupt "absent" into "present with garbage."

### Direction 1 — forwarding OUT of the pivot, into a still-raw child

A `%view` on an `optional` container yields an ordinary `pointer`
local. Guard the view, then forward the pointer unconditionally:

```fortran
! pivot's own dummy is now: type(RealArray_t), optional, intent(in) :: uhbt_a
real, dimension(:,:), contiguous, pointer :: uhbt   ! local, same name as the child's raw dummy

nullify(uhbt)
if (present(uhbt_a)) call uhbt_a%view(uhbt)
call some_still_raw_child(..., uhbt=uhbt, ...)
```

(These three lines are exactly what belongs in the subroutine — no
trailing comments; see Hard rules.)

This relies on a specific Fortran rule: **a disassociated pointer,
passed as the actual argument for an `optional` dummy that is not
itself `POINTER`/`ALLOCATABLE`, reads as not-present.** Different
mechanism from the "container passed by value is always present" trap
(`array_container_lessons` §9 #8b) — that's about derived-type values,
not pointers, and doesn't conflict here.

**The `nullify` is mandatory, not defensive — omitting it is a live,
confirmed silent-corruption bug.** A local pointer's association status
on entry is **undefined, not disassociated**, unless nullified as an
*executable statement* each call. (`pointer :: uhbt => null()` at
declaration is a different, also-wrong fix: it implicitly adds `SAVE`,
so the pointer is null only on the first-ever call and keeps whatever
it was left as afterward — wrong for anything called more than once.)
Skip the `nullify` and the pointer holds garbage on entry; `present()`
in the child can read `.true.` by accident and read that garbage as
real data — no crash, just silently wrong numbers. `nullify` every
Direction-1 pointer local, once, before the first guard that might view
it, regardless of which branch a given test happens to exercise.

A pointer only ever *dereferenced* inside the same `if (present(..._a))`
branch that set it needs no `nullify` (ordinary single-subroutine
container code) — the risk here is specific to *forwarding* a
possibly-untouched pointer into another procedure's `present()` check.

**Hard precondition, checked per child:** its dummy must be a plain
optional array — no `pointer`, no `allocatable`. A `pointer`/
`allocatable` child dummy reads a disassociated pointer as present but
disassociated, not absent — the wrong outcome. Read the declaration;
don't assume.

### Direction 2 — forwarding INTO the pivot, from its external callers

Whether a call site needs `if (present(...))` guarding depends on the
caller's actual argument for that slot, not on the pivot's new signature:

- **Fixed, compile-time-known source** (a plain array or expression,
  whether or not this particular call passes the keyword — that choice
  is baked into the line) → no branching. Build the container
  unconditionally (`source=` copy-in), pass by keyword, `free` after.
  The common case — true at all 15 external call sites in this
  codebase's worked example, since MOM6 dynamics code routinely
  hand-picks which optional outputs a call needs from persistent
  working arrays allocated unconditionally at init.
- **The caller's source is itself one of *its own* optional dummies**
  (a pass-through link) → the ordinary guarded-`alloc`-and-branch
  problem (`array_container_lessons` §6/§9 #8b), unchanged — a
  freshly-built container has no disassociated-value equivalent to
  Direction 1's pointer trick. Three or more such optionals live at one
  call site reintroduces combinatorial branching — stop and ask rather
  than nesting conditionals.

Check every caller individually; a tree usually has both kinds, and
getting the determination wrong is a real bug no linter or type-checker
catches — only a test that exercises the absent-argument path.

## Steps

### 1. Scope the conversion

Name the pivot and exactly which optional dummies convert (not
necessarily all of them). List: its mandatory dummies and any `G`/
`GV`/`US`/`CS`-style derived types staying untouched (mandatory
conversion is `convert_array_containers`'s job, separately); every
child each named argument forwards to, with the Direction 1 precondition
checked per child (stop and flag any failure); every external caller
(grep `call <pivot>(`, including continuations) — the Step 4 work queue.

### 2. Convert the named optional dummies on the pivot's signature

Same intent-mapping as `convert_array_containers` (lessons §6):
`optional, intent(out)` → `type(RealArray_t), optional, intent(inout)`
(never `intent(out)` on a container); `intent(in)` stays. Preserve
every `!<` doc comment verbatim; rename with the `_a` suffix.

Rename every `present(<name>)` test within the pivot's own body to
`present(<name>_a)`, including grouped checks (e.g.
`present(x) .neqv. present(y)`) — if a grouped partner isn't being
converted in this pass, stop and reconsider scope; a mixed group is a
correctness hazard, same as in `convert_present_to_associated`.

### 3. Rewire internal forwarding into still-raw children (Direction 1)

Declare one pointer local per converted argument, named like the
child's own raw dummy, once per pivot even across several call sites.
**Immediately `nullify` all of them in one statement, before anything
else runs** — not optional (see mechanism). For each (child call site)
× (converted argument): add the guarded `%view`, forward the local
pointer by keyword instead of the now-nonexistent raw dummy. Cover
every internal call site — a pivot commonly calls the same child once
per branch or per zonal/meridional direction; one top-level `nullify`
covers all of them regardless of which branch runs.

The child's own signature and body are untouched — verify this
explicitly in Step 5.

### 4. Update every external call site (Direction 2)

Apply the mechanism's determination per call site per converted
argument: fixed source → unconditional build (`source=` copy-in even
for `intent(out)`/`inout`, per lessons §6/§9#2), keyword pass, `copy2F`
back if written, `free`, no guard; caller's-own-optional source →
guarded `alloc` and branched call, flagging ≥3 live optionals at one
site before nesting conditionals. Arguments outside this pass's scope
pass through unchanged everywhere.

### 5. Verify

- Doc comments: `!<` count/text check on converted dummies (same
  mechanical diff as `convert_array_containers` Step 9).
- Line length ≤ 100.
- **Every child subroutine byte-identical to before** — the single most
  important check here, since staying out of scope is the whole point.
- Mechanical argument-count/keyword audit at every internal and
  external call site.
- Repeat `convert_array_containers` Step 9's diff-review/repo-grep
  checks (no loop-body edits, no missed call site).
- **Nullify audit:** every Direction-1 pointer-local block has a
  `nullify` covering all its members before any guard that views one —
  mechanical, always-do-this; omitting it compiles clean and silently
  corrupts results several call levels downstream, confirmed in this
  codebase.
- **This mechanism can't be verified by inspection alone** — it relies
  on a specific standard rule about disassociated-pointer arguments.
  State explicitly whether a build ran; if none is available (don't
  install one), say so and flag this as the top thing to confirm once a
  build is possible, ideally with a test exercising both the present
  and absent paths.

## Versioning marker

Every Fortran file this skill creates or modifies gets a `!!SKILLS: 0.3`
marker line — the shared version for this whole skill family. If
missing, add it right after the file's license/header block, before
`module`; if present, update it in place. Grep-able
(`grep -rn "!!SKILLS:"`), meant to be stripped later.

## Hard rules

- Never skip or duplicate the `!!SKILLS: 0.3` marker — add once, update
  in place thereafter.
- Never touch a child's signature or body — a separate
  `convert_array_containers` invocation if it needs converting too.
- Never apply the Direction-1 trick to a child dummy declared `pointer`
  or `allocatable` — check every child individually.
- Never skip the Direction-1 `nullify`, and never put it in the
  declaration instead of as an executable statement — the former is a
  live corruption bug, the latter implicitly adds `SAVE` (wrong for
  anything called more than once).
- Never assume a caller's source is safe to build unconditionally
  without checking whether it's the caller's own optional dummy —
  conflating the two reintroduces the "always present" bug, caught only
  by an absent-argument test.
- Do not convert the pivot's mandatory arguments here — that's
  `convert_array_containers`, separately.
- Do not convert `present()` to `associated()` here — that's
  `convert_present_to_associated`'s job, strictly afterward.
- Do not drop a `!<` comment, exceed 100 characters, or leave a call
  site unupdated.
- Do not claim an unrun build passed; never install a compiler.
- **No commentary about the conversion itself** — the `nullify`,
  guarded `%view`, and unconditional forward speak for themselves to
  anyone who's read this skill; don't annotate them in the source. A
  genuinely new dummy still gets a physical-meaning `!<` comment.

## Commit gating

Same as `convert_array_containers`: `--enable_git_commit`/
`--disable_git_commit` override; otherwise
`~/.claude/preferences.json`'s `git_commit_and_push` key. Branch:
`claude_<lowercased_function-name>_optional_containers`.

## Output to the user on success

1. The pivot and which optional dummies converted (and which were left
   out of scope, and why).
2. Every internal forwarding call site touched, with each child's
   pointer/allocatable precondition confirmed.
3. Every external call site touched, and which Direction-2 scenario
   applied at each.
4. Diff confirmation every child is untouched.
5. Build status, with an explicit callout if the disassociated-pointer
   mechanism is unverified by inspection alone.
6. Whether committed, or modified files for manual commit.
