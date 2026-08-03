---
name: convert_array_containers
description: Convert a MOM6 Fortran subroutine from raw grid-shaped arrays (real, dimension(SZI_(G),SZJ_(G),SZK_(GV))) to the RealArray_t / IntArray_t container types in place -- array dummies become containers, the body obtains ordinary Fortran pointers via %view so the math is untouched, and every call site is updated to alloc / copy2F / free around the call. Use when pushing the raw-array boundary further up the MOM6 call tree. Does NOT rename anything, add a dispatcher shim, capture mode, or a bind(C) bridge -- those belong to generate_cpp_bridge.
user-invocable: true
argument-hint: <work-directory> <function-name> [--enable_src_validate] [--enable_git_commit] [--disable_git_commit]
---

# Convert a MOM6 subroutine to array containers

This skill is the **execution checklist**. All templates, the container
API reference, intent rules, and pitfalls live in the
`array_container_lessons` skill (numbered sections §1–§10) — invoke
`array_container_lessons` once near the start of a session, before
running this skill, and refer back to its numbered sections from each
step below (Step 0 checks that it has been loaded). Do not reproduce
templates here.

**This skill converts data structures only.** It does not rename the
subroutine, does not add a `*_fortran` variant, and does not introduce a
dispatcher, `getenv_mode`, capture mode, `#ifdef _TIM`, `%to_c`, or a
`bind(C)` interface. Those belong to `generate_cpp_bridge` and are
layered on afterwards (lessons §10). A routine may be converted with no
commitment to ever bridging it.

## Help message

If `$ARGUMENTS` is empty, or equals `help`, or equals `--help`, or
equals `-h`, do NOT run any steps. Print the following help message
verbatim and stop:

```
Usage: /convert_array_containers <work-directory> <function-name> [--enable_src_validate] [--enable_git_commit] [--disable_git_commit]

Convert one MOM6 subroutine from raw grid-shaped Fortran arrays to
RealArray_t / IntArray_t containers, in place. The subroutine keeps its
name; its array dummies become containers; its body obtains pointers via
%view so the math is unchanged; and every call site is updated to
alloc / copy2F / free around the call.

Does NOT rename to *_fortran, add a dispatcher shim, capture mode, or a
bind(C) bridge -- use /generate_cpp_bridge for that, afterwards.

Arguments:
  <work-directory>   Absolute path to an existing TURBO-ESM/MOM6 checkout
                     (must contain src/ and config_src/). Cloning is not
                     performed.
  <function-name>    Name of the Fortran subroutine to convert
                     (case-insensitive match against its declaration).
  --enable_src_validate  (optional) Run Step 1: verify the work directory is a
                         TURBO-ESM/MOM6 checkout with a clean tree, and that
                         the subroutine exists. Off by default.
  --enable_git_commit    (optional) Force Step 9 to run this invocation
                         only, overriding the global git_commit_and_push
                         preference in ~/.claude/preferences.json.
                         Mutually exclusive with --disable_git_commit.
  --disable_git_commit   (optional) Force Step 9 to be skipped this
                         invocation only, overriding the global
                         preference. Mutually exclusive with
                         --enable_git_commit.

When neither --enable_git_commit nor --disable_git_commit is passed, Step 9
follows ~/.claude/preferences.json's "git_commit_and_push" key ("auto" or
"manual"; treated as "manual" when the file is missing, the key is absent,
or the JSON fails to parse).

Example:
  /convert_array_containers /glade/derecho/scratch/sunjian/MOM6 zonal_mass_flux
```

## Step 0 — validate inputs

`$0` = work-directory, `$1` = function-name. Run these checks **before any
other step** and stop on the first failure with a one-line, actionable
error. Do not retry, do not assume defaults, do not create anything.

1. **Argument count.** If `$0` empty OR `$1` empty → stop:
   `Error: missing arguments. Run "/convert_array_containers --help" for usage.`
2. **Work directory is an existing MOM6 checkout.** If `$0` is not an
   existing directory → stop:
   `Error: work directory "<value>" does not exist.`
   The directory must already contain a TURBO-ESM/MOM6 checkout —
   cloning is not performed by this skill.
3. **Parse optional flags.** Scan remaining arguments for
   `--enable_src_validate`, `--enable_git_commit`, and
   `--disable_git_commit`. Store as boolean flags (default: off). If both
   `--enable_git_commit` and `--disable_git_commit` are passed → stop:
   `Error: --enable_git_commit and --disable_git_commit are mutually exclusive.`
   Any unrecognised argument that starts with `--` → stop:
   `Error: unknown option "<value>". Run "/convert_array_containers --help" for usage.`
4. **Reference material primed.** This skill assumes the
   `array_container_lessons` skill has already been invoked earlier in
   this session. If you have no memory of its numbered sections (§2, §3,
   §4, §5, §6, §7, §8, §9 are referenced throughout this procedure),
   invoke the `array_container_lessons` skill now, before proceeding to
   Step 1 or any step below. Do not guess at template content or API
   signatures from memory of a similar prior task — invoke it and read
   its content first.
5. **A `<function-name>_fortran` sibling doesn't mean stop.** Run
   `grep -irn "^[[:space:]]*subroutine[[:space:]]\+<function-name>_fortran\b" $0/src $0/config_src`.
   If it matches, `<function-name>` is a bridge shim over an
   already-converted worker — convert it exactly like any other
   raw-array subroutine (Steps 2–9), with one constraint (lessons §10):
   never edit anything between `select case (mode)` and `end select`
   except renaming a `CS%x`/`GV%x` reference when Step 2 drops that
   derived type, and never touch a `bind(C)` interface declaration. If
   either seems necessary, stop — that means the task has drifted into
   `generate_cpp_bridge`'s territory.
6. **Target is not already converted.** Locate `<function-name>`'s own
   declaration with
   `grep -irn "^[[:space:]]*subroutine[[:space:]]\+<function-name>\b" $0/src $0/config_src`
   and read its dummy list. If every array dummy is already
   `type(RealArray_t)` / `type(IntArray_t)`, stop and report that the
   subroutine is already container-native; there is nothing to do.

Items 4, 5, and 6 always run, regardless of flags — a bad target is
worth catching before any git/branch work happens. Step 1 runs only when
`--enable_src_validate` is set, and adds deeper checkout/branch
validation on top of the same subroutine lookup. Step 9's behavior is
decided by `--enable_git_commit` / `--disable_git_commit` if passed, or
otherwise by the global preference described in Step 9. No other
conversion work executes until Step 0 validation passes.

## Settle these decisions (ask if not obvious from the tree)

1. **Iteration box** — if the subroutine already takes a `type(Box_t)`
   dummy (commonly `bxC`), use it. If it does not, the loop bounds are
   currently derived from `G`/`GV` index fields; decide with the user
   whether to add a `Box_t` dummy or keep the existing bound expressions.
   Do not invent a `Box_t` silently.
2. **Argument-list reduction** — whether `G` / `GV` / `US` / `CS` can be
   dropped from the signature once their array members become containers
   and their scalar members become plain dummies (lessons §8). Default:
   drop only what the body genuinely no longer references.
3. **Scope of call-site updates** — default is to update every call site
   found in Step 3 in the same change, so the build stays green. Each
   site is handled per its own caller's state (Step 8, Case A or B); a
   raw caller gets a marshalling block, an already-converted caller
   passes its containers straight through.

This skill converts **one subroutine per invocation** and does not
recurse. A multi-level conversion is a sequence of invocations, and
either direction works — but they differ in cost:

- **Top-down** (start at the highest routine you intend to convert, work
  toward the leaves): the marshalling block lands once, in the caller
  that will remain raw permanently, and never moves. Each subsequent
  conversion is mostly Step 8 Case B — swap a `%view` pointer argument
  for the container it came from. **Preferred.**
- **Bottom-up** (leaves first): every level's marshalling is written when
  its callee is converted and deleted again when the caller itself is
  converted. Correct, but churns the same code repeatedly.

Report the follow-up candidates found in Step 7 so the user can pick the
next invocation deliberately, and say which direction the current
invocation implies.

If the user already specified any of these, take their values as-is.

## Procedure

Each step is one action with a pointer to the `array_container_lessons`
section that holds the template or rationale.

### 1. Validate the checkout *(runs only when `--enable_src_validate` is passed; skip otherwise and proceed to Step 2)*
   If `--enable_src_validate` was not supplied, skip this entire step and go to Step 2.

   Verify `$0` is a TURBO-ESM/MOM6 checkout by running
   `git -C $0 remote -v` and confirming the output contains
   `TURBO-ESM/MOM6`. If not → stop:
   `Error: "<value>" is not a TURBO-ESM/MOM6 checkout.`

   - If the working tree is dirty (`git -C $0 status --porcelain` is
     non-empty) → stop and surface to the user; do not stash or discard.
   - Report the current branch (`git -C $0 branch --show-current`) in the
     plan confirmation below. This skill does **not** require a specific
     base branch — a conversion is branch-agnostic — but the user should
     see which branch is about to be modified. Do not switch branches.

   Then validate the tree — stop on the first failure with a one-line,
   actionable error:
   - **Looks like MOM6.** If `$0/src` or `$0/config_src` is missing → stop:
     `Error: "<value>" is not a MOM6 source tree (missing src/ or config_src/).`
   - **Subroutine present.** Run
     `grep -irn "^[[:space:]]*subroutine[[:space:]]\+<function-name>\b" $0/src $0/config_src`.
     - 0 matches → stop: `Error: subroutine "<function-name>" not found under <work-directory>/{src,config_src}.`
     - >1 match → list candidates and ask the user which to convert.
     (Step 0 items 5 and 6 already ruled out an already-bridged or
     already-converted target before this step runs; no need to
     re-check either here.)

### 2. Classify every dummy argument
   Read the full subroutine. For each dummy record: name, declared type,
   `intent`, rank, and whether it is a grid-shaped array
   (`real, dimension(SZI_(G),…)`), a scalar, a `logical`, a `pointer`, or
   a derived type.

   Decide per dummy (lessons §6, §7):
   - Grid-shaped `real` array → `type(RealArray_t)`, renamed `<name>_a`.
   - Grid-shaped `integer` array → `type(IntArray_t)`, renamed `<name>_a`.
   - Scalar / `logical` / `pointer` (`OBC`) → unchanged, passes through.
   - Derived type (`BT_cont_type`, …) → unchanged (lessons §9 #11).

   **Stop and ask about every `optional` array dummy — do not decide
   either way silently.** This skill has produced both outcomes in
   practice on real subroutines (e.g. `zonal_mass_flux`'s `uhbt`,
   `visc_rem_u`, `u_cor`, and `du_cor` were left as raw optional arrays;
   other conversions have containerized theirs) and neither is a safe
   default to assume.

   For each optional array dummy found, present it to the user by name
   and ask: convert it now, or leave it as a raw array, passed through
   unconverted.

   - **Convert now.** It becomes a container and keeps its `optional`
     attribute; its `%view` must be guarded by `present(<name>_a)`, and
     every existing `present(<name>)` test in the body must be renamed
     to `present(<name>_a)` (lessons §6, §9 #8a). Note it now so Steps 4,
     5 and 8 all handle it; a missed `present` guard compiles and fails
     at run time.
   - **Leave as raw.** Make no change to its declaration, rank, or any
     `present(<name>)` test in the body. At call sites (Step 8) this
     argument is simply forwarded as-is — it was never a container, so
     there is nothing to marshal for it.

   Record each answer before proceeding to Step 3; every optional array
   dummy must have an explicit answer, not an inferred one.

   Then scan the **body** for derived-type array references
   (`G%mask2dT`, `G%areaT`, …) and scalar references (`GV%Angstrom_H`,
   `CS%upwind_1st`, …). Array references become **new** container
   dummies; scalar references become plain scalar dummies computed at
   the call site (lessons §8).

   Print the proposed new argument list and pause for user confirmation
   before writing anything. Name every added and removed dummy
   explicitly.

### 3. Find every call site
   Run `grep -rn "call[[:space:]]\+<function-name>(" $0/src $0/config_src`
   and also check for continuation-line call forms
   (`grep -rn -A2 "call[[:space:]]\+<function-name>[[:space:]]*&"`).
   List every file and line. This list is the Step 8 work queue — a
   conversion that misses a call site breaks the build (lessons §9 #9).

   If any call site is inside a routine you were not asked to touch,
   that is expected: the marshalling block lands there (lessons §1).

### 4. Rewrite the dummy declarations
   Apply the classification from Step 2 to the subroutine's declaration
   block. Written containers take `intent(inout)`, never `intent(out)`
   (lessons §9 #6).

   **Every `!<` doc comment must survive verbatim.** This is the single
   most commonly broken part of a conversion, because the original
   declarations are usually *two* lines with the comment on the
   continuation line, while the converted declaration is one line — so
   collapsing the two silently drops the comment:

   ```fortran
   ! BEFORE -- comment lives on the continuation line
     real,  dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(in)    :: h_in !< Tracer cell layer thickness [H ~> m or kg m-2].

   ! AFTER -- WRONG, comment dropped during the collapse
     type(RealArray_t),    intent(in)    :: h_in_a

   ! AFTER -- CORRECT, comment carried onto the new single line
     type(RealArray_t),    intent(in)    :: h_in_a     !< Tracer cell layer thickness [H ~> m or kg m-2]
   ```

   Procedure to guarantee this:
   1. Before editing, record each array dummy's doc comment text
      verbatim, including unit annotations like `[H ~> m or kg m-2]` and
      `[nondim]`. Scan the **whole** declaration — the comment may be on
      the `real, dimension(...)` line or on any continuation line.
   2. Write the new declaration and re-attach the recorded comment to it.
      Keep the wording exactly as it was; do not reword, abbreviate, or
      "improve" it, and never shorten or drop words to make a line fit —
      wrap it instead (item 5). A trailing period stays or goes exactly
      as it was.
   3. For a **new** dummy that had no prior comment (a grid-derived
      container such as `mask2dT_a`, or a scalar lifted out of `CS`),
      write a fresh `!<` comment in the same house style — every dummy in
      this codebase carries one, and Doxygen builds warn on any that
      does not.
   4. Align the `!<` column with its neighbours for readability, but
      never at the cost of the text.
   5. **Check the resulting line length — MOM6 enforces a 100-character
      limit.** Collapsing a two-line raw-array declaration onto one
      container-declaration line very often pushes it past 100
      characters, especially when the comment carries a unit annotation
      (`[H ~> m or kg m-2]`, `[nondim]`). If the new line exceeds 100
      characters, wrap the comment onto a `!!` continuation line at a
      word boundary — do not truncate, shorten, or drop any of the
      comment to make it fit. Match the continuation style already used
      by neighbouring dummies in the *same* subroutine (indent the `!!`
      under the `!<` column above it); this file already wraps long
      comments this way in many places, so there is always a local
      example to match. Recheck the length of both resulting lines
      before moving on.
   6. **The `subroutine <name>(...)` header line is subject to the same
      100-character limit — check it separately from the dummy
      declarations.** Adding new grid-derived container dummies (Step 2)
      lengthens the argument list on this line, and it is easy to check
      every declaration's length while forgetting the one line that
      names them all. If the header's first physical line exceeds 100
      characters, move the `&` continuation earlier — wrap the argument
      list across more lines, breaking at a comma — the same mechanism
      already used to wrap it at all; do not remove or rename an
      argument to make it fit. Recheck the length of every physical line
      of the header, not just the first.

   Comments that are **not** on a dummy declaration — the `!>` Doxygen
   header above the subroutine, comment blocks between declarations, and
   anything inside the body — must be left exactly where they are and
   exactly as they read.

   Do **not** rename the subroutine (lessons §9 #1).

### 5. Add the `%view` locals and calls
   Declare one pointer per container, using the **original** array name:
   `real, dimension(:,:,:), contiguous, pointer :: h_in, h_W, h_E`
   (matching the rank of each array). Immediately before the first use,
   `call <name>_a%view(<name>)` for each (lessons §2, §4.3).

   For an `optional` dummy flagged in Step 2, guard the view and rename
   its `present` tests (lessons §6):
   `if (present(hin_a)) call hin_a%view(hin)`, and every existing
   `present(hin)` in the body becomes `present(hin_a)`.

   The loop body must then require **zero edits** *apart from those
   `present` renames* — the pointers carry the container's original
   `lb:ub` bounds. If you find yourself editing an index expression or a
   formula, stop and re-check: something is wrong.

### 6. Wire up the iteration box
   If the routine already takes a `Box_t`, use it as-is. If the body
   needs a grown or shrunk domain, derive it with
   `bx = bxC%grow(dim=N, n=M)` — `dim=1` zonal/i, `dim=2` meridional/j
   (lessons §5) — and **`call bx%free()`** when done; `grow` returns a
   new box you own (lessons §9 #3).

   Leave existing loop-bound expressions alone unless the user chose to
   add a `Box_t` in "Settle these decisions" #1.

### 7. Handle `elemental` callees and unconverted descendants
   If the body calls an `elemental` scalar kernel (`flux_elem`,
   `flux_elem_OBC`) inside a loop, leave that kernel's signature alone
   and pass `%view` pointer elements to it exactly as before
   (lessons §9 #10).

   If the body calls a subroutine that still takes raw arrays, that
   callee is *below* the boundary and is out of scope for this
   invocation: pass it the `%view` pointer, not the container. Note it
   in the Step 10 report as a follow-up conversion candidate.

### 8. Update every call site from Step 3
   **Classify each call site first.** For every site, look at the
   *calling* subroutine's own dummy declarations and decide which of two
   cases applies. Do not assume — a conversion campaign will hit both,
   often in the same file.

   **Case A — the caller still takes raw arrays.** Insert the marshalling
   block (lessons §3): declare a container per array argument, `alloc`
   every array **with** `source=` — including `intent(out)` ones, unless
   you have confirmed the callee's write covers the array's full
   `LBOUND`/`UBOUND` rather than just a `Box_t`-scoped sub-region
   (lessons §6) — call the converted routine, `copy2F` outputs and
   inouts back, then `free` every container. Derive any new scalar
   dummies here (`h_min = 2.0 * GV%Angstrom_H`).

   This is the common case when converting bottom-up, and it is also the
   permanent end state when the caller is a routine that will never be
   converted (a module entry point such as `continuity_PPM`).

   *If the caller's own source array is `optional`*, a container local is
   always present as an actual argument, so an unallocated one would make
   `present()` true in the callee with garbage behind it. Allocate
   conditionally and branch the call (lessons §6, §9 #8b). With more than
   two optional arguments at one site the branches multiply — stop and ask
   the user rather than nesting conditionals.

   *Loop-invariant sources:* if a container's source never changes —
   grid metadata such as `G%IareaT`, `G%mask2dT`, `G%areaT` — do **not**
   restructure the caller to hoist it out of a loop or share it across
   call sites as part of this conversion. Emit the straightforward
   per-site `alloc`/`free` and **report it in Step 10** as a hoisting
   candidate instead. Hoisting is a whole-subroutine optimisation that is
   safer to do once, with every call site visible, than incrementally.

   **Case B — the caller is already container-based.** Do **not** add a
   marshalling block. The caller already holds containers, so pass them
   straight through:

   - If the caller has the matching container as one of its own dummies,
     pass it directly. No `alloc`, no `copy2F`, no `free`.
   - If the caller was passing a `%view` pointer down to this routine
     because it used to be raw (Step 7 of that caller's own conversion),
     replace the pointer argument with the container it came from, and
     delete the `%view` pointer local if nothing else in the body uses it.
   - If the caller holds a marshalling block for *this* call from an
     earlier bottom-up pass, that block is now dead scaffolding
     (lessons §3): delete the `alloc` / `copy2F` / `free` calls and the
     container locals, and pass the caller's own containers instead.
   - Only a genuinely new array — one the caller has no container for,
     e.g. a local work array — still needs its own `alloc`/`free` pair.

   Case B is the common case when converting top-down, and it is what
   makes that direction cheaper: the marshalling is written once at the
   permanent boundary rather than added and removed at every level.

   **Both cases:** do not move existing `cpu_clock_begin`/`cpu_clock_end`
   boundaries, and ensure early returns still reach any `free` block
   (lessons §9 #15). After editing, re-check that the argument count and
   order at every site match the new signature.

   **Raw-vs-container check, argument by argument.** Before moving to the
   next call site, re-read Step 2's list of dummies that stayed raw for
   *this* subroutine (an `optional` array left unconverted, a `logical`
   array with no container type, or anything else) and confirm the actual
   argument at that position is the raw pointer/array/value, not a
   container. This is easiest to get wrong precisely when a subroutine
   mixes converted and raw array dummies: a call site is usually built by
   mechanically appending `_a` across most of the argument list, and a
   raw dummy's caller-side actual often has a same-stem container sibling
   sitting right next to it in scope — e.g. a raw `uh_3d` dummy fed by the
   caller's own `uh`, when the caller also holds `uh_a` for a different,
   already-converted argument. Passing `uh_a` there instead of `uh` is a
   real type mismatch (`type(RealArray_t)` vs `real`), and only a Fortran
   compiler catches it — which this skill cannot assume is available
   (Step 9). This caused a real build failure in an already-reviewed
   conversion (`zonal_flux_adjust`'s `uh_3d`), caught only once a CI build
   ran. The more raw dummies a subroutine has, the higher the odds of
   this slip — check it explicitly rather than trusting the mechanical
   edit.

### 9. Verify
   Run the checks in lessons §10:
   - **Doc comments intact.** Count `!<` occurrences in each edited
     subroutine before and after; the count must not drop. Any dummy
     whose comment went missing is a defect — restore it from the
     original before continuing. A quick mechanical check:
     `git -C $0 diff -U0 -- <file> | grep '^-' | grep -c '!<'` versus
     `git -C $0 diff -U0 -- <file> | grep '^+' | grep -c '!<'`; the added
     count must be greater than or equal to the removed count (greater
     when new dummies were introduced in Step 4). Then read the removed
     lines and confirm each comment text reappears somewhere in the
     added lines.
   - **Line length.** MOM6 enforces a 100-character limit (Step 4 items 5
     and 6) — this applies to every added line, not just dummy
     declarations, including the `subroutine <name>(...)` header line
     itself. Check every added line: `git -C $0 diff -U0 -- <file> |
     grep '^+' | grep -v '^+++' | awk '{ print length($0)-1 }' | sort -rn
     | head -1` — if the top value exceeds 100, find and wrap the
     offending line(s) (`git -C $0 diff -U0 -- <file> | grep '^+' | grep
     -v '^+++' | awk 'length($0)-1 > 100'`). A collapsed declaration line
     is the usual cause (wrap the comment onto a continuation line); a
     lengthened argument list on the header line is the other common
     cause (move the `&` continuation earlier, Step 4 item 6). Never
     shorten text or drop an argument to make a line fit.
   - Diff review: no edits inside any loop body; the subroutine name is
     unchanged; only declarations, `%view` calls, and call-site blocks
     changed.
   - Every call site from Step 3 updated; argument count and order match
     the new signature at each.
   - **Raw-vs-container argument audit.** Walk every call site's actual
     arguments position by position against the callee's current
     signature — matching count and order (above) does not catch a
     same-position type swap (Step 8).
   - Nothing used outside the FMS2 ∩ TIM shared API surface
     (lessons §4) — in particular, no `%to_c` (lessons §9 #12).

   **Build:** do not assume a Fortran compiler is available on this
   machine; many developer machines on this project do not have one, and
   this skill must never try to install one. If a compiler is present,
   build under **both** infra layers. If not, say plainly that the build
   was not run and is pending, and never describe the conversion as
   verified on inspection alone.

### 10. Commit and push to `claude_<function-name>_containers` *(gated by the global git-commit preference, overridable per run — see below)*
   **Decide whether to run this step:**
   1. If `--disable_git_commit` was passed, skip this entire step and
      report the files that were modified so the user can commit
      manually — this run's explicit override wins.
   2. Else if `--enable_git_commit` was passed, run this step — this
      run's explicit override wins.
   3. Else, read the global preference: run
      `cat ~/.claude/preferences.json 2>/dev/null` and inspect the
      `git_commit_and_push` key.
      - If it is `"auto"`, run this step.
      - If it is `"manual"`, or the file does not exist, or the key is
        absent, or the JSON fails to parse — skip this entire step and
        report the files that were modified so the user can commit
        manually. Absent or unreadable configuration means "manual";
        never push on an unconfigured machine.

   From the work directory, create (or check out, if it already exists)
   the branch `claude_<lowercased_$1>_containers`, based on whatever
   branch the checkout is currently on. Stage every modified file and
   commit with a message naming the converted subroutine, the arguments
   that became containers, and the call sites updated. Use a HEREDOC so
   the trailer is preserved verbatim:

   ```
   BRANCH="claude_$(echo "$1" | tr '[:upper:]' '[:lower:]')_containers"
   git -C "$0" checkout -B "$BRANCH"
   git -C "$0" add -A
   git -C "$0" commit -m "$(cat <<'EOF'
   <one-line summary of the container conversion for $1>

   <1-3 line body: which dummies became containers, which grid-derived
   arrays were added, which call sites were updated, build status>

   Co-authored-by: Claude <noreply@anthropic.com>
   EOF
   )"
   git -C "$0" push -u origin "$BRANCH"
   ```

   If the push is rejected because the branch already exists upstream
   with unrelated history, stop and surface the conflict to the user
   rather than force-pushing.

## Hard rules

- Do not rename the subroutine, and do not create a `*_fortran` variant
  or a wrapper pair. That is `generate_cpp_bridge`'s job.
- Do not add a dispatcher, `select case (mode)`, `getenv_mode`, capture
  mode, `io_recorder`, `#ifdef _TIM`, `%to_c`, or a `bind(C)` interface
  where none already exists — that is always `generate_cpp_bridge`'s job.
  If a target already has this content (a `_fortran` sibling exists),
  convert it normally (item 5) but never edit inside
  `select case … end select` except renaming a `CS%x`/`GV%x` reference,
  and never touch a `bind(C)` interface declaration (lessons §10).
- Do not drop a `!<` doc comment. Every dummy that had one keeps it,
  word for word, when its declaration is rewritten — including when a
  two-line declaration collapses to one line, which is where they are
  usually lost (Step 4). Every dummy added by the conversion gets a new
  one.
- Do not leave a collapsed declaration line, or the subroutine's own
  `subroutine <name>(...)` header line, over MOM6's 100-character limit.
  Never fix an overlong line by shortening or dropping comment text or
  an argument — wrap a declaration's comment onto a `!!` continuation
  line, or move the header's `&` continuation earlier, instead (Step 4
  items 5–6, Step 9).
- Do not edit anything inside a loop body. A conversion changes
  declarations and adds `%view` calls; the math is untouched.
- Do not change the numerical result. A conversion is inert by
  construction; any difference is a bug.
- Do not pass a written container as `intent(out)` — always
  `intent(inout)` (lessons §9 #6).
- Do not skip `source=` for an `intent(out)` array just because it's
  "pure output" — that only reads garbage if the write covers the array
  in full, which a `Box_t`-scoped write almost never does (lessons §6,
  §9 #2). Default to copying in; skip it only when confirmed safe.
- Do not leave a `grow`/`shrink` result unfreed (lessons §9 #3).
- Do not use `%dup` — it no longer exists; its functionality is folded
  into `%alloc` (lessons §4.6).
- Do not convert the members of a derived type (`BT_cont_type`) unless
  explicitly asked.
- Do not decide unilaterally whether an `optional` array dummy becomes a
  container or stays raw — Step 2 stops and asks per dummy; both
  outcomes are common and neither is a default.
- Do not leave a call site unupdated.
- Do not pass a container actual to a dummy Step 2 classified as raw, or
  vice versa (Step 8, Step 9).
- Do not attempt to install a Fortran compiler, and do not claim a build
  passed that was never run.

If something not covered here comes up, consult lessons §9 (recurring
pitfalls) before improvising.

## Output to the user on success

Report:

1. The subroutine converted, and its new argument list.
2. Which dummies became containers, and which grid-derived arrays were
   added as new container dummies.
3. Which of `G` / `GV` / `US` / `CS` were dropped from the signature, if
   any.
4. Every call site updated (file and line).
5. Any callee still taking raw arrays, as a follow-up conversion
   candidate (Step 7).
6. **Hoisting candidates:** every container whose source is
   loop-invariant grid metadata (`G%IareaT`, `G%mask2dT`, …), listed with
   the call sites that now rebuild it. These are deliberately left
   per-site (Step 8) and are worth a single cleanup pass once the
   enclosing subroutine's conversions are complete.
7. Every `optional` array dummy found, the user's decision for each
   (converted vs. left raw), and whether call sites needed conditional
   allocation for any converted ones.
8. Build status — which infra layers were built, or explicitly that no
   build was run and why.
9. Whether the change was committed, or the list of modified files for
   manual commit.
