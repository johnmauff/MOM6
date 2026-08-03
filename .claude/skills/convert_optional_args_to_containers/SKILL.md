---
name: convert_optional_args_to_containers
description: Convert only the optional array arguments of a MOM6 subroutine that sits in the middle of a call tree to RealArray_t/IntArray_t containers, leaving its mandatory arguments and every still-raw child subroutine it calls completely untouched. Uses a guarded %view producing a disassociated pointer, forwarded unconditionally to the child's still-raw optional dummy, to avoid the combinatorial if(present(...)) branch trees that a naive Case-A conversion produces. Use this when a full convert_array_containers pass on the subroutine (or its neighbors) isn't wanted yet, or when the branching cost of converting an optional argument bottom-up has already bitten once and you want to fix it by moving the container boundary up past the argument's real point of use instead.
user-invocable: true
argument-hint: <work-directory> <function-name> <optional-dummy-name>[,<optional-dummy-name>...] [--enable_git_commit] [--disable_git_commit]
---

# Convert a subroutine's optional arguments to containers, without converting its children

## Why this exists

`convert_array_containers` converts every array dummy of one subroutine —
mandatory and optional alike — and its call-site template (lessons §3)
assumes the callee ends up container-based. That's the right tool when
the whole subroutine is moving to containers. It is the *wrong* tool for
a narrower, common situation: a subroutine (call it the **pivot**) has a
few optional array arguments that it simply forwards on, unmodified, to
one or more still-raw **child** subroutines further down the call tree —
and only those optional arguments need containers, typically because
converting them in isolation, bottom-up, at the child instead of the
pivot, produced a branch-tree mess.

**The mess this fixes, concretely:** converting a child's optional
arguments to containers while its pivot caller stays raw forces the
pivot to build a *fresh* container from its own raw optional dummy and
then conditionally include it in the call — and a freshly-built local
container is not the same actual argument as the pivot's own dummy, so
Fortran's "absent optional forwarded to another optional stays absent"
rule no longer applies. The pivot ends up enumerating every
present/absent combination of every optional argument by hand: a real
case in this codebase produced a 16-way branch tree at each of a
pivot's two call sites into the child, plus an *unrelated* 8-way branch
tree inside the child itself, for forwarding those same values one
level further down.

**The fix:** convert the *pivot's* optional dummies to containers
instead (or in addition), and thread containers through its own
external callers. The pivot can then forward into its still-raw
children with **zero branching in either direction** — see mechanism
below. The children are never touched; if they also need containers
eventually, that's a separate `convert_array_containers` invocation.

This skill assumes the `array_container_lessons` skill has already been
invoked this session for the container API (`%alloc`/`%view`/`%copy2F`/
`%free` signatures, intent-mapping rules, doc-comment preservation). It
is not reproduced here.

## The mechanism (read this before doing anything)

There are two independent directions of data flow through the pivot,
and each has its own branch-free trick. Get either one wrong and you
either reintroduce combinatorial branching or corrupt "absent" into
"present with garbage."

### Direction 1 — forwarding OUT of the pivot, into a still-raw child

A `%view` on an `optional` container yields an ordinary Fortran
`pointer` local. Guard the view itself, then forward the pointer
**unconditionally**:

```fortran
! pivot's own dummy is now: type(RealArray_t), optional, intent(in) :: uhbt_a
real, dimension(:,:), contiguous, pointer :: uhbt   ! local, same name as the child's raw dummy

nullify(uhbt)
if (present(uhbt_a)) call uhbt_a%view(uhbt)
call some_still_raw_child(..., uhbt=uhbt, ...)
```

(The three lines above are exactly what belongs in the actual subroutine
— no trailing comments. The explanation of *why* each line is there
lives in this skill document, not in the code; see "No commentary about
the conversion itself" under Hard rules below.)

This works because of a specific, standard Fortran rule: **a
disassociated pointer, passed as the actual argument for a dummy that is
`optional` but does *not* itself have the `POINTER` or `ALLOCATABLE`
attribute, causes that dummy to read as not-present (`present()` is
`.false.`) in the callee.** This is a different mechanism from the
"a container passed by value is always present" trap described in
`array_container_lessons` §9 #8b — that trap is about derived-type
*values*; it does not apply to pointers, and this trick does not
contradict it.

**The `nullify` is not optional and not merely defensive — omitting it
is a live, silent-corruption bug, confirmed in this codebase.** A plain
local pointer's association status on subroutine entry is **undefined**,
not disassociated, unless the pointer is explicitly nullified as an
*executable statement* each time the subroutine runs. (Initializing it
in the declaration itself, `pointer :: uhbt => null()`, is a different
and *also wrong* fix here — that gives the pointer the `SAVE` attribute
implicitly, so it is null only on the very first call ever and then
retains whatever it was left as after every subsequent call, which is
not what a subroutine invoked once per timestep needs.) Skip the
`nullify` and the pointer is left in whatever undefined state the stack
happened to hold whenever the corresponding `_a` container is absent —
`present()` in the child can then read `.true.` by accident, and the
child proceeds to read garbage floating-point values as if they were
real data. This does not reliably crash; it silently corrupts results,
exactly the failure mode this skill exists to avoid introducing. `nullify`
every Direction-1 pointer local, for every converted argument, once,
before the first `if (present(..._a))` guard that might view it — not
just the ones a particular test happens to exercise as absent.

This is a materially different situation from a pointer local that is
only ever *dereferenced* within the same `if (present(..._a))` branch
that produced it (as in ordinary, single-subroutine container code) —
that pattern is safe even without `nullify`, because the undefined
pointer is never read. The risk here is specific to *forwarding* a
possibly-untouched pointer as an actual argument to another procedure's
`present()` check, which is exactly what Direction 1 does.

**Hard precondition — check this for every child before applying the
trick:** the child's own dummy must be a plain optional array (no
`pointer`, no `allocatable`). If the child's dummy is itself declared
`pointer` or `allocatable`, passing a disassociated pointer makes the
child's dummy *disassociated but still present* — `present()` returns
`.true.` there, which is the wrong outcome. Read the child's actual
declaration; do not assume.

### Direction 2 — forwarding INTO the pivot, from its external callers

At every call site that calls the pivot, whether a container needs
`if(present(...))` guarding depends on what the caller's actual argument
for that slot *is*, not on the pivot's new signature:

- **The caller's source is a plain, unconditionally-available raw array
  (or expression), and whether this specific call includes the keyword
  is a fixed, compile-time choice already baked into that line of
  source** (a different call site elsewhere in the same file might omit
  it — that's fine, it's a different, separately-authored line) — then
  **no branching is needed at all.** Build the container unconditionally
  right before the call (`source=` copy-in), pass it by keyword, `free`
  it after. This is the common case, and it was true at all 15 external
  call sites in the worked example below — MOM6 dynamics-core code
  routinely hand-picks which optional outputs a given call needs, as a
  compile-time decision, from persistent working arrays that are
  allocated unconditionally at init.
- **The caller's source is itself one of *the caller's own* optional
  dummies** (i.e. the caller is a pass-through link in a longer optional
  chain) — this is the ordinary problem from `array_container_lessons`
  §6 / §9 #8b, unchanged: guard the `alloc` with `if (present(...))`,
  and branch the call itself, because a freshly-built container is not
  the same forwarding scenario as Direction 1's pointer trick (the
  pivot's dummy is a container *value*, not a pointer — there is no
  disassociated-value equivalent). With three or more such caller-side
  optionals live at one call site, this reintroduces combinatorial
  branching at that *one* site — stop and ask the user rather than
  generating a thicket of nested conditionals, exactly as
  `convert_array_containers` already advises for this situation.

Check every caller individually — a real call tree usually has some of
each kind, and getting this determination wrong in either direction is a
real bug that neither a linter nor a type-checker catches, only a test
that actually exercises the absent-argument path.

## Steps

### 1. Scope the conversion

Name the pivot subroutine and exactly which of its optional array
dummies are being converted (some or all — this skill does not require
converting every optional argument the pivot has). List:
- Its mandatory dummies and any derived types it still needs whole
  (`G`/`GV`/`US`/`CS`-style) — confirm these are staying untouched. If
  the user also wants the mandatory dummies converted, that's
  `convert_array_containers`'s job, run separately or first.
- Every child subroutine the pivot forwards each named optional argument
  to. For each, read the child's own dummy declaration and confirm the
  Direction 1 precondition (plain optional array, not pointer/
  allocatable). If a child fails this precondition, stop and flag it —
  this skill cannot proceed for that argument until resolved.
- Every external caller of the pivot (grep for `call <pivot>(`,
  including continuation-line forms). This is the Step 4 work queue.

### 2. Convert the named optional dummies on the pivot's own signature

Same intent-mapping as `convert_array_containers` (lessons §6): a
dummy that was `optional, intent(out)` becomes
`type(RealArray_t), optional, intent(inout)` (never `intent(out)` on a
container — lessons §9 #6); `intent(in)` stays `intent(in)`. Preserve
every `!<` doc comment verbatim, rewrapping only for the 100-character
line limit (lessons §7, same procedure as `convert_array_containers`
Step 4). Rename the dummy with the established `_a` suffix.

Rename every `present(<name>)` test on this argument **within the
pivot's own body** to `present(<name>_a)` — including any grouped
presence check (e.g. `present(x) .neqv. present(y)`) if `y` is also
being converted in this pass; if a grouped partner is *not* being
converted in this pass, stop and reconsider scope — a mixed group is a
correctness hazard the same way it is in `convert_present_to_associated`.

### 3. Rewire the internal forwarding calls into still-raw children (Direction 1)

Declare one pointer local per converted argument, named the same as the
child's own raw dummy (so the call itself reads unchanged and no other
body edits are needed) — once per pivot, even if used across several
internal call sites. **Immediately after declaring them, `nullify` all
of them in a single executable statement, before any other code in the
subroutine runs.** This is not optional — see the mechanism section
above for why an un-nullified pointer is a live correctness bug, not a
defensive nicety.

For each (child call site) × (converted argument) pair after that: add
the guarded `%view`, and change the call from directly forwarding the
now-nonexistent raw dummy to forwarding the local pointer by keyword.
Do this for every internal call site — a pivot commonly calls the same
child more than once (e.g. once per branch of an `if`/`else`, or once
per direction of a zonal/meridional pair). One `nullify` at the top
covers every branch; do not scatter per-branch `nullify` calls or skip
one because "that branch always sets it anyway" — the whole point is
not to depend on which branch runs.

The child's own signature and body are not touched at all. Verify this
explicitly in Step 6 — it's the cheapest, most mechanical confirmation
that scope was actually kept narrow.

### 4. Update every external call site (Direction 2)

For each call site found in Step 1, and each converted optional
argument it passes: apply the caller-source determination from the
mechanism section above.

- **Fixed/compile-time source:** build the container unconditionally
  (`source=` copy-in — the argument may be `intent(out)`/`inout` in the
  pivot and Box_t-scoped, so copy-in is required per lessons §6/§9#2
  even here), pass by keyword, `copy2F` back afterward if the pivot
  writes to it, `free`. No guard anywhere.
- **Caller's own optional source:** guard the `alloc` and branch the
  call (lessons §6 pattern). Flag to the user if ≥3 such optionals are
  live at one call site before generating nested conditionals.

Arguments not being converted in this pass (mandatory, or optional
arguments intentionally left out of scope) pass through completely
unchanged at every call site.

### 5. Verify

- Doc comments on the converted dummies: `!<` count check, verbatim-text
  confirmation (same mechanical diff check as `convert_array_containers`
  Step 9).
- Line length ≤ 100 on every added line.
- **Confirm every child subroutine's signature and body are
  byte-identical to before** — `git diff` should show zero hunks inside
  any child's own subroutine boundaries. This is the single most
  important verification for this skill specifically, since the whole
  point is that children are out of scope.
- Mechanical argument-count/keyword audit (a short Python script over
  the extracted call statements, as used for `convert_array_containers`)
  at every internal forwarding call and every external call site.
- Repeat the diff-review and repo-wide-grep checks from
  `convert_array_containers` Step 9 (no loop-body edits, no call site
  missed).
- **Grep every Direction-1 pointer-local block for a `nullify` covering
  every one of them, before any `if (present(..._a))` guard that views
  one.** This is a mechanical, always-do-this check — do not skip it
  because a particular subroutine "always takes the present branch in
  practice." Confirmed in this codebase: omitting it compiles cleanly,
  passes casual inspection, and silently corrupts numerical results
  several call levels downstream, only surfacing as unexplained
  numerical drift in CI rather than a build failure.
- **This skill's Direction-1 mechanism cannot be verified by inspection
  alone even with the `nullify` in place** — it depends on a specific
  Fortran standard rule about disassociated-pointer arguments to
  non-pointer optional dummies. State explicitly whether a compiler
  build was run; if none is available on this machine (do not attempt
  to install one), say so plainly and flag this as the highest-value
  thing to confirm once a build is possible — ideally with a test that
  actually exercises both the present and absent paths through the
  pivot, not just a clean compile.

## Hard rules

- Never touch a child subroutine's own signature or body. If a child
  also needs converting, that's a separate `convert_array_containers`
  invocation, not part of this one.
- Never apply the Direction-1 `%view` + unconditional-forward trick to a
  child dummy declared `pointer` or `allocatable` — check every child
  individually; the trick's correctness depends on this precondition and
  it is not automatically true.
- **Never skip the `nullify` of a Direction-1 pointer local, and never
  put it inside the declaration (`=> null()`) instead of as an
  executable statement** — the former is a live silent-corruption bug
  (undefined initial association status read as spuriously present
  downstream); the latter implicitly gives the pointer the `SAVE`
  attribute, which is wrong for anything called more than once (it
  would carry a stale association from a previous call instead of
  resetting each time).
- Never assume a caller's source is "safe to build unconditionally"
  without first checking whether it is itself the caller's own optional
  dummy — conflating the two reintroduces the "container is always
  present" bug, and nothing but an absent-argument test run catches it.
- Do not convert the pivot's mandatory arguments as part of this skill —
  that's `convert_array_containers`, run separately if wanted.
- Do not convert `present()` to `associated()` here — that is
  `convert_present_to_associated`'s job, and is about bridge-readiness,
  not about eliminating branching; keep this skill's changes strictly to
  introducing containers under the existing `optional`/`present()`
  idiom.
- Do not leave a `!<` doc comment dropped, an added line over 100
  characters, or a call site unupdated — same standing rules as
  `convert_array_containers`.
- Do not claim a build passed that was never run, and do not attempt to
  install a Fortran compiler.
- **No commentary about the conversion itself.** Do not add comments to
  the pivot's code explaining that a dummy became a container, that a
  `nullify` is mandatory, why Direction 1 works, or anything else about
  the migration — this codebase documents ocean physics, not its own
  refactoring history (lessons §7). The `nullify` statement, the guarded
  `%view`, and the unconditional forward speak for themselves to anyone
  who has read this skill and the lessons reference; do not annotate
  them in the source. A genuinely new dummy still gets a `!<` doc
  comment describing its physical meaning — unchanged, and not an
  exception to this rule.

## Commit gating

Same as `convert_array_containers`: `--enable_git_commit` /
`--disable_git_commit` override; otherwise follow
`~/.claude/preferences.json`'s `git_commit_and_push` key (`"auto"` runs
Step 6's commit; `"manual"`, missing file, missing key, or unparseable
JSON all mean skip it and report modified files for manual commit).
Branch name: `claude_<lowercased_function-name>_optional_containers`.

## Output to the user on success

Report, mirroring `convert_array_containers`'s style:
1. The pivot subroutine and which optional dummies were converted (and
   which, if any, were explicitly left out of scope, and why).
2. Every internal forwarding call site touched, and confirmation that
   each target child's pointer/allocatable precondition was checked.
3. Every external call site touched, and which of the two Direction-2
   scenarios applied at each (fixed-source vs. caller's-own-optional).
4. Confirmation (via diff) that every child subroutine is untouched.
5. Build status — and an explicit callout that the disassociated-pointer
   mechanism is unverified by inspection alone if no compiler was run.
6. Whether committed, or the list of modified files for manual commit.
