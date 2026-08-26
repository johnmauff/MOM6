---
name: convert_calltree
version: "0.3"
description: Given a call-tree entry point (a subroutine whose external signature must stay fixed, like continuity()), survey every derived type and array/optional argument reachable beneath it, classify each against the existing container/bridge-readiness skills, get a human decision wherever a fixed rule can't resolve one, record the plan, then execute it end to end. Use this instead of re-deriving the survey-and-decide process by hand for each new entry point.
user-invocable: true
argument-hint: <work-directory> <entry-point-name> [--enable_git_commit] [--disable_git_commit]
---

# Convert a call-tree entry point to be bridge-ready, end to end

## Why this exists

Making a call tree bridge-ready takes several different targeted
transformations — containerizing raw arrays, converting optional
dummies, shadowing or bundling a derived type, bridging each leaf — and
none of the skills that perform them decide *which* one applies to
*which* target, or in what order. That survey-and-decide work was done
by hand for `continuity()`; it doesn't scale to the ~10 more entry
points expected to need the same treatment.

**Why the entry point is the unit of work.** Its external signature is
the one thing that can't change — everything below it is free to be
restructured, because nothing below it is visible to the rest of MOM6.
This skill adds no new transformation; it adds the missing layer above
the existing ones: survey, classify, decide once with a human, record
the decisions, then execute them.

## The three phases

**Phase 1 (interactive)** — survey the tree, classify every target, ask
about the ones a fixed rule can't resolve, write a durable plan.
**Phase 2 (checkpointed)** — build every container, bundle, and shadow;
stabilize the tree. **Phase 3 (checkpointed)** — bridge every
subroutine, leaf to root, once the tree is stable. Both execution
phases share the same discipline: no decisions left to make, but a
commit/push/CI gate between every stage, not one unattended run.

Phase 1 alone is a complete, valid use of this skill (e.g. to review
before committing to the campaign); Phase 2 and Phase 3 are separate,
later invocations that read the plan back. Phase 3 assumes Phase 2 is
done — it bridges a tree that's already fully container/bundle-ready,
never a still-raw one.

## Scope

An **entry point** is a subroutine whose external signature must not
change; its descendants, transitively, are fair game. A **target** is
anything Phase 1 must classify: a derived type referenced anywhere in
the tree, a raw array dummy, an optional dummy.

**Wrapper / TreeRoot convention.** The entry point's own body needs
somewhere to attach entry-level logic (a shadow's build/copy-back, a
bundle's construction) without touching the tree below it. Step 1d
classifies which of three cases applies; Phase 2 (Step 5a) acts on it
before anything else:

- **Case 1 — a wrapper already exists**, any shape or name. Nothing to do.
- **Case 2 — the entry point is a bodiless alias for a distinctly-named
  implementation** (e.g. a `use ..., only : entry=>impl` rename). Keep
  the implementation's name; author a new subroutine under the entry
  point's own name whose body just calls the implementation, and fix up
  the `use`/`public` lines.
- **Case 3 — no separate name exists at all** (one subroutine, called
  directly). Rename it `<entry>_TR` ("TreeRoot"), then author the
  wrapper the same way as Case 2, calling `_TR`.

Worked example: pre-campaign `continuity()` was a zero-body
`use MOM_continuity_PPM, only : continuity=>continuity_PPM` alias over
`continuity_PPM` — Case 2. `_TR` is a distinct convention from
`generate_cpp_bridge`'s own `_fortran` rename: `_TR` marks the tree's
root, once; `_fortran` marks a bridged leaf's original implementation,
per leaf, later.

Whichever case, the wrapper's body starts as a pure pass-through and
never gets bridged — every later stage that needs entry-level logic
edits it, never the implementation/`_TR` subroutine.

## Hard precondition

1. Confirm the entry point's signature really is externally fixed — if
   it isn't, there may be no reason to scope the campaign this tightly.
2. Check every descendant for callers *outside* this tree, repo-wide. A
   descendant reached from more than one entry point can't have its
   containerization decided in isolation — flag it in the plan (Step 4)
   as "shared with `<other entry point>`, decide jointly."

## Phase 1 — survey, classify, decide, record

### Step 1. Inventory

1a. The entry point's own current public dummy list — the contract
    Phase 2 must never change.
1b. Every derived type referenced (directly, or via `%field`), every
    raw array dummy, and every optional dummy, down to the leaves. Use
    an Explore/fork agent for anything beyond a handful of subroutines
    — don't re-derive a survey this size from memory.
1c. Every descendant's caller list, repo-wide (the shared-descendant
    check, above) — do this for the whole tree now, not per-target later.
1d. Which Scope-section case applies (1/2/3). Record the tree's actual
    root name, and whether Phase 2 must author a wrapper (Cases 2/3)
    and/or a rename (Case 3 only).

### Step 2. Classify each target against fixed rules

- Raw array dummy → `convert_array_containers`.
- Optional array/scalar dummy → `convert_present_to_associated`, after
  containerizing if it's an array.
- A struct/pointer only ever forwarded opaquely (confirm by grep, never
  assume from its name) → leave alone; record the confirmation.
- A private, scalar-or-container-field control structure whose fields
  recur together across signatures → `create_config_bundle_type`.
- A struct used elsewhere in the repo outside this tree →
  `create_shadow_container_type`, not a wholesale conversion.

### Step 3. Quantify and decide what Step 2 couldn't resolve

Typically a large multi-field struct (`G`/`GV`/`US`-shaped) too big to
shadow or bundle without measuring first, or anything whose shared-vs-
private status isn't obvious. For each: measure total fields touched
vs. total fields, and whether the touched ones recur together (a
bundling/shadowing candidate) or are each used at only one or two sites
(not worth it). Present the measurement and the real tradeoff via
`AskUserQuestion` — one topic at a time, never bundled into one
multi-question call. Never pick a default silently.

### Step 4. Record the plan

Write a durable manifest to
`<work-directory>/.claude/calltree-plans/<entry-point-name>.md` — not
the ephemeral Plan-mode file. For every target: which skill handles it,
its settings (field lists, "leave alone" plus its grep confirmation,
"shared with `<X>`" flags), and the execution order below. This file is
what makes the campaign resumable without re-running Step 3, and what
Phase 2 and Phase 3 read.

**Phase 2 execution order.** Each position is forced by some sub-
skill's own precondition or direction rule, not chosen arbitrarily;
record any deviation the survey justifies rather than reordering
silently.

1. **TreeRoot split** (Scope) — skip for Case 1.
2. **`create_shadow_container_type`**, per target classified "shared
   outside the tree" — built/copied-back inside the wrapper.
3. **`create_config_bundle_type`**, per private clustering candidate —
   skip if none.
4. **Optional-array containerization**, per subroutine: default to
   `convert_array_containers`; use `convert_optional_args_to_containers`
   instead, scoped to one pivot, when converting its children first
   would branch combinatorially on `present()`. Never both, for the
   same argument. `convert_present_to_associated` (item 8) is separate
   and later — not an alternative here. **Never leave an optional array
   dummy raw on a target this plan also schedules for Phase 3
   bridging** — `generate_cpp_bridge` requires every array dummy to
   already be a container (its own Step 0 precondition), so "leave as
   raw" is only a valid choice here for a target that will never be
   bridged. Check Phase 3's list (below) before answering this question.
5. **`convert_array_containers` — downward pass**, root to leaves,
   mandatory dummies (that skill's preferred direction for dummy
   conversion).
6. **`convert_array_containers` — upward pass**, leaves to root,
   `G`/`GV`/`US`-drop decisions and Step 2b promotions (that skill's
   *required* direction for these — undecidable top-down).
7. **`convert_locals_to_containers`**, per subroutine once its dummies
   and its callees' are stable — scratch locals feeding an
   already-converted callee.
8. **`convert_present_to_associated`**, per dummy about to cross a
   `bind(C)` boundary, ahead of Phase 3 bridging that leaf. Precondition:
   already a container (item 4).
9. **`hoist_container_marshalling`**, once, at the tree's root, last.

**Phase 3 execution order.** Not a fixed list — computed from Step 1b's
call graph, restricted to every subroutine in the tree except the
wrapper (which never bridges): a subroutine is ready to bridge once
every in-tree callee it still calls is already bridged, so leaves go
first and the tree's root goes last. Record the computed wave order
(the ready-together groups) in the plan; Phase 3 (Step 7) reads it.

## Phase 2 — checkpointed execution

### Step 5. Execute one stage at a time, gated by commit/push/CI

One branch for the whole run, `claude_<lowercased-entry-point>_calltree`,
checked out once before Stage 1 (confirm the tree is clean first). Pass
`--disable_git_commit` to every invoked sibling skill, unconditionally
— this skill's Step 5 is the only thing that commits, so the run lands
on one branch, not one per sibling.

For each of the 9 stages (Step 4's Phase 2 order), in sequence:

1. **Do the work** — author code directly (Stage 1) or invoke every
   sibling skill the stage's targets are recorded against, each with
   `--disable_git_commit`. A stage often means several invocations; the
   checkpoint is per stage, not per invocation.
2. **Verify** — every invoked skill's own Verify section, plus this
   skill's Step 6 checks scoped to just this stage's files (the
   full-tree sweep is Stage 9's job). If verification finds a problem,
   fix it before committing.
3. **Commit** onto the branch, message naming the stage and what ran.
4. **Push**, then check CI (confirm how this repo's CI is wired; don't
   assume). No local compiler exists here, so this is the first real
   build check the stage gets.
5. **Stop and report** — stage, commit, CI status — then wait. Never
   start the next stage in the same turn, and never proceed past a
   failed stage or an unresolved verification problem; surface it and
   let the user decide.

If a stage surfaces a target Phase 1 didn't anticipate, stop and return
to Phase 1 for it rather than guessing.

**Step 5a — Stage 1's work**, per Step 1d's case: Case 1 → nothing,
commit skipped. Case 2 → author the wrapper under the entry point's
name, calling the implementation's existing name. Case 3 → rename the
implementation to `_TR` first, then author the wrapper. (See Scope for
what each case means.) This is the one piece of code this skill authors
itself rather than dispatching to a sibling.

### Step 6. Whole-tree verification (Stage 9)

One full sweep, before Stage 9's commit, alongside
`hoist_container_marshalling` itself: argument count/order across the
whole tree, zero stray references anywhere to a field that should have
moved, doc-comment/line-length/`#ifdef` balance for every touched file,
and the entry point's external signature (1a) still byte-identical.
Also confirm Step 1d's case played out correctly (Case 2: implementation
kept its name; Case 3: `_TR` exists, pre-split body unchanged), and that
the wrapper's own signature matches 1a and carries no `bind(C)`
interface.

## Phase 3 — leaf-to-root bridging

### Step 7. Bridge one wave at a time, leaves first, gated by commit/push/CI

Same branch as Phase 2 (`claude_<lowercased-entry-point>_calltree`,
already checked out); do not create a new one. Pass
`--disable_git_commit` to every `generate_cpp_bridge` invocation,
unconditionally — it has its own commit step, and this skill's Step 7
is what commits instead, same reasoning as Phase 2's Step 5.

**Precondition check, per target, before invoking `generate_cpp_bridge`
on it — do not assume Phase 2 left every target ready.** Confirm every
array dummy is already `type(RealArray_t)`/`type(IntArray_t)`, and the
iteration domain already a `type(Box_t)` (`generate_cpp_bridge`'s own
Step 0 item 5). If a target fails this — most likely a Step 4-item-4
decision that left an optional array dummy deliberately raw before this
target was scheduled for bridging — stop for that target: either
containerize it now (`convert_array_containers`, a Phase-2-shaped fix
applied here because Phase 2 already ran) or exclude the target from
this campaign's bridging and record why. Never force a container
conversion through unreviewed, and never skip the target silently.

For each wave (Step 4's Phase 3 order — every target whose in-tree
callees are already bridged), in sequence:

1. **Do the work** — invoke `generate_cpp_bridge` on every target in the
   wave, each with `--disable_git_commit`, after its precondition check
   passes.
2. **Verify** — each invocation's own Step 9 (the three-mode matrix,
   `cpp_bridge_lessons` §17); CAPTURE mode is the real bar here, since
   no AMReX C++ implementation exists yet (that's a separate, later
   skill's job — see Step 8).
3. **Commit** onto the branch, message naming the wave and which
   subroutines it bridged.
4. **Push**, then check CI.
5. **Stop and report** — wave, commit, CI status — then wait. Same
   rules as Phase 2's Step 5: never start the next wave in the same
   turn, never proceed past a failed wave or a precondition-check
   failure.

The tree's root (Step 1d's recorded root — `continuity_PPM`-equivalent)
is the last wave; the wrapper is never bridged, in any wave.

### Step 8. Whole-tree bridging verification

After the last wave: confirm every non-wrapper subroutine in the tree
has a `_fortran`/shim/`bind(C)` triple, the wrapper still has none, and
every shim's public signature is unchanged from what Phase 2 left it
(bridging must not alter a signature — only rename and wrap). Report
the AMReX side (the C++ implementation behind each `bind(C)` interface)
as the explicit next deliverable this skill does not produce.

## Versioning marker

Every Fortran file this skill creates or modifies — directly (Step 5a)
or via an invoked sibling, in Phase 2 or Phase 3 — gets a
`!!SKILLS: 0.3` marker line, the shared version for this whole skill
family. Deliberately grep-able (`grep -rn "!!SKILLS:"`) and meant to be
stripped later.

## Hard rules

- Never skip the `!!SKILLS: 0.3` marker on a file Step 5a edits directly.
- Never skip the shared-descendant check (1c).
- Never let Phase 2 or Phase 3 resolve a judgment call Phase 1 didn't
  record — return to Phase 1 instead.
- Never bypass an invoked sibling's own hard rules or preconditions —
  this skill only decides *which* skill and *when*, never *how*.
- Never default an ambiguous classification (Step 3) silently.
- Never deviate from Step 1d's case: no cosmetic rename of an existing
  (Case 1) wrapper, no renaming a Case 2 implementation, no skipping
  Step 5a for Case 2/3.
- Never bridge the wrapper, and never confuse its `_TR` rename with
  `generate_cpp_bridge`'s `_fortran` rename — different marks, different
  points in the pipeline.
- Never let an invoked sibling commit or branch during Phase 2 or
  Phase 3 — always pass `--disable_git_commit`.
- Never advance to the next stage or wave in the same turn, past one
  with an unresolved verification problem, or past a red CI run.
- Never invoke `generate_cpp_bridge` on a target whose precondition
  check (Step 7) hasn't passed — fix or exclude the target first, never
  force a container conversion through unreviewed.
- Never bridge a target before every in-tree callee it still calls is
  already bridged (Step 7's wave order) — leaves first, root last.
- Do not claim an unrun build passed; never attempt to install a
  compiler.

## Commit gating

Phase 2 and Phase 3 each commit on one branch — the same branch,
created once before Phase 2's Stage 1, never per-sibling branches.
Whether either phase commits/pushes at all follows the usual chain:
`--enable_git_commit`/`--disable_git_commit` on `convert_calltree`
itself, else `~/.claude/preferences.json`'s `git_commit_and_push` key.
Once committing, every invoked sibling still gets `--disable_git_commit`
regardless (Hard rules). Applies to Phase 2/3 only — Phase 1 produces a
plan file, not code changes.

## Output to the user on success

**Phase 1:** the plan file path, and a summary table (target → skill →
decision) for review before Phase 2 runs.
**Phase 2, per stage:** the stage number, what ran, the commit hash, CI
status — then stop.
**Phase 2, after Stage 9:** Step 6's verification results, confirmation
the entry point's signature is unchanged, and a rollup of all 9 stages'
commits.
**Phase 3, per wave:** the wave number, which subroutines it bridged,
the commit hash, CI status — then stop.
**Phase 3, after the last wave:** Step 8's verification results, and the
explicit list of `bind(C)` interfaces still needing an AMReX-side
implementation.
