---
name: convert_present_to_associated
version: "0.3"
description: Convert an already-array-container-ified optional dummy argument from Fortran's `optional`/`present()` idiom to a non-optional container checked via `x%associated()`, for bind(C)/C++ bridge-readiness (C/C++ has no `optional`, only nullable pointers). The unit of conversion is an ARGUMENT NAME across its whole by-value forwarding chain, not a single subroutine. Precondition: the argument must already be a RealArray_t/IntArray_t/LogicalArray_t at every subroutine in that chain (run convert_array_containers first where it is still raw). Does not touch the container API beyond the %associated() method, does not add a bind(C) bridge -- that is generate_cpp_bridge's job.
user-invocable: true
argument-hint: <work-directory> <dummy-argument-name>[,<dummy-argument-name>...] [<seed-function-name>]
---

# Convert `present(x)` to `x%associated()`

## Why this exists

`bind(C)` can't express "this argument may be omitted." The fix is to
make the argument **always present as a container**, whose internal
pointer may be unassociated — "absent" becomes "associated with
nothing," not "not passed." Query this via `%associated()` (wraps
`associated(this%data)`, identical across container types and infra
layers) rather than `x%data` inline, which would leak the field name
outside `array_mod`.

Separate from `convert_array_containers`: converting a raw array to a
container, and converting its presence-test idiom, are independent
transformations that happen to often apply to the same variable.

## Step 0 — scope: one subroutine, or the whole forwarding chain

**The unit here is an argument name, not a subroutine** — the one
respect this skill differs from the rest of the family. A dummy's
presence-idiom must match at both ends of every place it's forwarded
**by value** (not via a guarded `%view`, which carries its own presence
signal). Naming one subroutine is only safe if the argument is
self-contained there; otherwise find the full set first:

1. **Seed.** Start from `<seed-function-name>` if given (confirm it
   declares the argument as an optional container there); otherwise
   grep for every subroutine declaring an optional dummy of that name
   and keep the ones that are the same logical quantity. A known root
   of the chain is a good seed — its own upward check should terminate
   almost immediately, confirming it's really the top.

2. **Close the set.** From every member, check both directions and
   repeat until a pass adds nothing:
   - *Downward:* does it forward the argument by value into another
     subroutine's own same-named optional dummy? Add that callee.
   - *Upward:* grep the whole tree for every caller. If a caller's
     actual argument is *its own* same-named optional dummy, add it.
   Apply both checks to every member, including ones found by the
   *other* direction — a node reached only upward can have its own
   downward callees, and checking only from the seed can miss a
   sibling caller that shares a callee further down but never appears
   in the seed's own tree.
   A branch ends at a subroutine that consumes the argument directly
   (`%view`/`%associated()`), or hands it to a *mandatory* dummy under
   a different name (a separately-decided variable). A caller that
   passes a fixed value or omits the argument outright (not its own
   same-named optional dummy) is a **boundary**: its call site still
   needs Step 4, but it isn't a set member.

3. The closed set is one consistent boundary: every member gets Steps
   1–3 applied; every edge call site (internal, external, or
   non-`optional`) gets Step 4. An empty closure — one subroutine — is
   still the common outcome, just no longer assumed without checking.

## Hard precondition

Every argument in the set must already be a container (run
`convert_array_containers` first where still raw). Every caller
currently omitting the argument needs Step 4 in the same pass — never
leave one end expecting a container while the other still omits it.

## Grouped arguments

If this argument's presence is checked together with another's (e.g.
`present(a) .and. present(b)`, or the routine treats a set as
all-or-nothing), convert the whole group together — converting one
member while a sibling stays on `present()` produces mixed semantics
for what was meant to be one logical state.

## Steps

1. **List every `present(<name>)` occurrence**, for every subroutine in
   the set, before changing anything.

2. **Drop `optional`** from each declaration:
   ```fortran
   ! before
   type(RealArray_t), optional, intent(in) :: uhbt
   ! after
   type(RealArray_t), intent(in) :: uhbt
   ```

3. **Replace `present(x)` with `x%associated()`**, everywhere in the
   set — confirm each occurrence actually tests the container's
   presence, not an unrelated similarly-named flag.
   ```fortran
   ! before
   if (present(uhbt)) then
   ! after
   if (uhbt%associated()) then
   ```
   `%associated()` must be `pure` — F2018 forbids referencing an impure
   procedure inside `do concurrent`, where these checks often sit.
   Confirm once per container type, not per call site.

   Two independent, confirmed-real-bug checks apply to every guard
   rewritten here, regardless of when this subroutine last saw
   `convert_array_containers`:
   - If the guard also conditionally sets a `pointer` local via `%view`,
     that pointer needs a `nullify` before the first such guard.
   - If the guard sits inside an `!$omp target`/`teams`/`parallel`/
     `loop`-family construct, never leave `x%associated()` itself as
     the in-loop condition — precompute a scalar before the construct.

4. **Fix every call site at the edge of the set.** One that already
   builds the container unconditionally needs nothing. One that omits
   the argument needs a dedicated, never-`%alloc`'d placeholder local
   (its `data` stays `null()`) declared *within* the subroutine that
   actually contains that call — confirm by that subroutine's own line
   range, never by proximity to other edits or a same-named "real"
   container elsewhere in the file.

   **Never let one placeholder serve two dummies in the same call.**
   Fortran forbids the same actual argument associated with two dummies
   in one invocation when either is definable (`intent(out)`/`inout`) —
   the restriction is on declared intent, not on which branch runs, and
   a compiler may assume non-aliased memory. Check the intents in play
   before sharing a variable across two argument positions.

5. **Leave `%to_c()` alone** — already null-safe for an unassociated
   container; this skill only changes the Fortran-side presence idiom.

## Verify

Programmatically, not by eye — a wrong-scope placeholder or a missing
`pure` both compile clean past a casual read.

- Zero remaining `present(<name>)`/`optional ... <name>`, repo-wide,
  for every argument in the set.
- Every call site's placeholder declared within its actual enclosing
  subroutine's line range.
- No actual argument repeated within one call's full list (Step 4's
  aliasing hazard).
- **Nullify audit:** every `pointer` local set by a guarded
  `%view`, anywhere in each set subroutine, has a `nullify` before the
  first such guard.
- **OpenMP-construct audit:** no `%associated()` guard, anywhere in
  each set subroutine, is the condition of an `!$omp target`/`teams`/
  `parallel`/`loop`-family construct.
- Line length ≤ 100; doc comments preserved verbatim.

## Explicitly not performance motivated

`%associated()` and `present()` are both O(1); Fortran resolves neither
at compile time. The only reason to do this is bridge-readiness.

## Versioning marker

Every Fortran file this skill creates or modifies gets a `!!SKILLS: 0.3`
marker line — the shared version for this whole skill family. If
missing, add it right after the file's license/header block, before
`module`; if present, update it in place. Grep-able
(`grep -rn "!!SKILLS:"`), meant to be stripped later.

## Hard rules

- Never skip or duplicate the `!!SKILLS: 0.3` marker.
- Never convert an argument still declared as a raw array.
- Never convert one member of a `present()`-linked group alone.
- Never leave a caller omitting the argument after its callee drops
  `optional` — same pass, always.
- Never scope this to one subroutine without running Step 0's closure.
- Never assume a call site's enclosing subroutine from proximity —
  check its line range.
- Never let two dummies in one call share an actual argument when
  either is definable.
- Never declare `isAssociated<Type>` without `pure`.
- Never rewrite a `present(x)`/`x%associated()` guard without checking
  both the `nullify` and OpenMP-construct hazards independently — two
  distinct, confirmed bugs.
- Do not add a `bind(C)` interface, dispatcher, or capture mode — that's
  `generate_cpp_bridge`'s job, strictly afterward.
