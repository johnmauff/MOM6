---
name: generate_cpp_bridge
description: Wrap an already-container-based MOM6 Fortran subroutine in a runtime-dispatched shim that selects between (a) the original Fortran code, (b) a binary capture mode that records inputs+outputs to disk for offline validation, and (c) a C++/AMReX bridge invoked through bind(C). Use when porting any MOM6 kernel to AMReX while keeping the Fortran truth available as a numerical reference. Requires the target subroutine's array dummies to already be RealArray_t / IntArray_t -- run convert_array_containers first. Mirrors the pattern established in TURBO-ESM/MOM6 PR #15.
user-invocable: true
argument-hint: <work-directory> <function-name> [--enable_src_validate] [--enable_git_commit] [--disable_git_commit]
---

# Generate C++ bridge for a MOM6 Fortran subroutine

This skill is the **execution checklist**. All templates, rationale,
type-mapping tables, and pitfalls live in the `cpp_bridge_lessons` skill
(same numbered sections, §1–§17) — invoke `cpp_bridge_lessons` once near
the start of a session, before running this skill, and refer back to its
numbered sections from each step below (Step 0 checks that it has been
loaded). Do not reproduce templates here.

**Precondition — data structures are already converted.** This skill
operates on a subroutine whose array dummies are *already*
`type(RealArray_t)` / `type(IntArray_t)` and whose iteration domain is
already a `type(Box_t)`. Converting raw grid-shaped arrays
(`real, dimension(SZI_(G),…)`) into containers is the job of the
**`convert_array_containers`** skill; run it first. Step 0 checks this
and stops if the target is not yet converted.

The two skills are **orthogonal** and must stay that way:

| Owned by `convert_array_containers` | Owned by this skill |
|---|---|
| Array dummies → containers | The `*_fortran` rename |
| `%view` pointers in the body | The dispatcher shim + `select case (mode)` |
| Loops → `do concurrent` over a `Box_t` | `getenv_mode` and the env-var contract |
| Introducing a `Box_t` dummy | Capture mode / `io_recorder` |
| Call-site `alloc` / `copy2F` / `free` | `#ifdef _TIM` AMReX arm, `%to_c` |
| | The `bind(C)` interface declaration |
| | Halo-check and cpu-clock relocation |

Never perform work from the left-hand column here. If the target needs
it, stop and tell the user to run `convert_array_containers` first.

## Help message

If `$ARGUMENTS` is empty, or equals `help`, or equals `--help`, or equals
`-h`, do NOT run any steps. Print the following help message verbatim
and stop:

```
Usage: /generate_cpp_bridge <work-directory> <function-name> [--enable_src_validate] [--enable_git_commit] [--disable_git_commit]

Wrap a MOM6 Fortran subroutine in a runtime-dispatched shim that selects
between the original Fortran code, a capture mode for offline validation,
and a C++/AMReX bridge. Mirrors TURBO-ESM/MOM6 PR #15.

PRECONDITION: the target subroutine's array dummies must ALREADY be
RealArray_t / IntArray_t containers, with a Box_t iteration domain. If it
still takes raw arrays (real, dimension(SZI_(G),...)), run
/convert_array_containers on it first -- this skill does not convert data
structures.

Arguments:
  <work-directory>   Absolute path to an existing TURBO-ESM/MOM6 checkout
                     (must contain src/ and config_src/, and should be on
                     the dev/turbo-debug branch). Cloning is not performed.
  <function-name>    Name of the Fortran subroutine to wrap (case-insensitive
                     match against its declaration in the tree). Must
                     already be container-based.
  --enable_src_validate  (optional) Run Step 1: verify the work directory is a
                         TURBO-ESM/MOM6 checkout on dev/turbo-debug and that
                         the subroutine exists. Off by default.
  --enable_git_commit    (optional) Force Step 10 to run this invocation
                         only, overriding the global git_commit_and_push
                         preference in ~/.claude/preferences.json.
                         Mutually exclusive with --disable_git_commit.
  --disable_git_commit   (optional) Force Step 10 to be skipped this
                         invocation only, overriding the global
                         preference. Mutually exclusive with
                         --enable_git_commit.

When neither --enable_git_commit nor --disable_git_commit is passed, Step 10
follows ~/.claude/preferences.json's "git_commit_and_push" key ("auto" or
"manual"; treated as "manual" when the file is missing, the key is absent,
or the JSON fails to parse).

Example:
  /generate_cpp_bridge /glade/derecho/scratch/sunjian/MOM6 PPM_limit_pos
```

## Step 0 — validate inputs

`$0` = work-directory, `$1` = function-name. Run these checks **before any
other step** and stop on the first failure with a one-line, actionable
error. Do not retry, do not assume defaults, do not create anything.

1. **Argument count.** If `$0` empty OR `$1` empty → stop:
   `Error: missing arguments. Run "/generate_cpp_bridge --help" for usage.`
2. **Work directory is an existing MOM6 checkout.** If `$0` is not an
   existing directory → stop:
   `Error: work directory "<value>" does not exist.`
   The directory must already contain a TURBO-ESM/MOM6 checkout —
   cloning is not performed by this skill. Step 1 verifies the checkout
   identity and branch state.
3. **Parse optional flags.** Scan remaining arguments for `--enable_src_validate`,
   `--enable_git_commit`, and `--disable_git_commit`.
   Store as boolean flags (default: off). If both `--enable_git_commit` and
   `--disable_git_commit` are passed → stop: `Error: --enable_git_commit and
   --disable_git_commit are mutually exclusive.` Any unrecognised argument
   that starts with `--` → stop:
   `Error: unknown option "<value>". Run "/generate_cpp_bridge --help" for usage.`
4. **Reference material primed.** This skill assumes the `cpp_bridge_lessons`
   skill has already been invoked earlier in this session. If you have no
   memory of its numbered sections (§3.1, §3.2, §5, §12, §13, §14, §15,
   §16, §17 are referenced throughout this procedure), invoke the
   `cpp_bridge_lessons` skill now, before proceeding to Step 1 or any step
   below. Do not guess at template content from memory of a similar prior
   task — invoke it and read its content first.
5. **Target is already container-based.** Locate the subroutine's
   declaration and inspect its dummy arguments. Every grid-shaped array
   dummy must already be `type(RealArray_t)` / `type(IntArray_t)`, and the
   iteration domain must already be a `type(Box_t)`. If any dummy is still
   a raw grid-shaped array (`real, dimension(SZI_(G),…)` or similar) →
   stop:
   `Error: "<function-name>" still takes raw arrays. Run "/convert_array_containers <work-directory> <function-name>" first, then re-run this skill.`
   Do not convert the data structures yourself — that is
   `convert_array_containers`' job, and doing it here would duplicate
   that skill's work and produce an unreviewable combined diff.

Step 1 runs only when `--enable_src_validate` is set. Items 4 (reference-
material priming) and 5 (container precondition) always run, regardless of
flags. Step 10's behavior is decided by `--enable_git_commit` /
`--disable_git_commit` if passed, or otherwise by the global preference
described in Step 10. No other wrapping work executes until Step 0
validation passes.

## Settle these decisions (ask if not obvious from the tree)

1. **Bridge prefix** — default to whatever existing `_bridge) bind(C)`
   declarations in the tree use; otherwise `turbotmp_`.
2. **Env-var name** — default `<UPPERCASE_$1>_MODE`.
3. **AMReX gating** — default `#ifdef _TIM` around the `case (TIMH_runAMREX)`
   arm only; capture mode is never gated.

The shim's public argument list is **not** a decision: it is inherited
verbatim from the converted subroutine (containers in, containers out),
so existing call sites need no changes.

If the user already specified any of these, take their values as-is.

## Procedure

Each step is one action with a pointer to the lessons.md section that
holds the template or rationale.

### 1. Validate the existing checkout *(runs only when `--enable_src_validate` is passed; skip otherwise and proceed to Step 2)*
   If `--enable_src_validate` was not supplied, skip this entire step and go to Step 2.

   **Verify the work directory.** All wrapping work in this skill is
   based on the `dev/turbo-debug` branch — that is the agreed base for
   bridge work and what the C++ side expects to merge against.

   Verify `$0` is a TURBO-ESM/MOM6 checkout by running
   `git -C $0 remote -v` and confirming the output contains
   `TURBO-ESM/MOM6`. If not → stop:
   `Error: "<value>" is not a TURBO-ESM/MOM6 checkout.`

   Then run `git -C $0 fetch origin dev/turbo-debug` and check the
   branch state:
   - If the working tree is dirty (`git -C $0 status --porcelain` is
     non-empty) → stop and surface to the user; do not stash or discard.
   - If HEAD is not already at `origin/dev/turbo-debug` (compare
     `git -C $0 rev-parse HEAD` and
     `git -C $0 rev-parse origin/dev/turbo-debug`) → stop and ask the
     user whether to `git checkout dev/turbo-debug` before continuing.
     Do not switch branches silently — the user may have unrelated work
     on the current branch.

   **Then validate the tree** before proceeding to Step 2 — stop on the
   first failure with a one-line, actionable error:
   - **Looks like MOM6.** If `$0/src` or `$0/config_src` is missing → stop:
     `Error: "<value>" is not a MOM6 source tree (missing src/ or config_src/).`
   - **Subroutine present.** Run `grep -irn "^[[:space:]]*subroutine[[:space:]]\+<function-name>\b" $0/src $0/config_src`.
     - 0 matches → stop: `Error: subroutine "<function-name>" not found under <work-directory>/{src,config_src}.`
     - >1 match → list candidates and ask the user which file to wrap.
   - **`lessons.md` present.** If `$0/.claude/skills/generate_cpp_bridge/lessons.md` is missing → stop:
     `Error: lessons.md not found at <work-directory>/.claude/skills/generate_cpp_bridge/lessons.md.`
   - **Confirm the plan.** Print one paragraph naming: the resolved
     subroutine file path, the proposed env-var (`<UPPERCASE_$1>_MODE`),
     the proposed bridge symbol (`<prefix>_$1_bridge`; pick `<prefix>` by
     `grep -h "_bridge) bind(C)" $0/src $0/config_src` — fall back to
     `turbotmp_`), and the caller files found by
     `grep -irl "call[[:space:]]\+<function-name>(" $0/src $0/config_src`.
     Then proceed to Step 2.

### 2. Classify each dummy argument for the C boundary
   The dummies are already containers (Step 0 item 5 verified this). For
   every arg, record intent / kind / rank / optional, and map it to its
   `bind(C)` counterpart using lessons.md §3.1–§3.2:
   `type(RealArray_t)` → `type(RealArray_C)`, `type(Box_t)` → `type(Box_C)`,
   scalars → `real(c_double)` / `integer(c_int)` / `logical(c_bool)` by
   value, pointers → `type(c_ptr)`.

   Do **not** introduce a `Box_t`, and do not alter any dummy's type. If
   the iteration domain is not already a `Box_t`, stop and refer the user
   to `convert_array_containers` — Step 0 item 5 should have caught this.

### 3. Rename the original implementation
   Rename `$1` → `$1_fortran` in place. **That is the only change to this
   subroutine.** Its body, dummy declarations, `%view` calls, and loop
   structure are already container-based and must be left byte-for-byte
   unchanged apart from the name on the `subroutine` / `end subroutine`
   lines. Template for what the result should look like: lessons.md §12.

   Do not convert array dummies (already done), do not rewrite loops
   (already `do concurrent` over the `Box_t`), and do not touch the math.

### 4. Add the `bind(C)` interface block
   At the top of the host module, declare `<prefix>_$1_bridge` per the
   template in lessons.md §3.3. Doc-comment every dummy with `!<`.

### 5. Write the shim subroutine `$1`
   Takes the public name the original had, and **inherits that
   subroutine's container-based dummy list verbatim** — same names, same
   types, same intents, same `!<` doc comments. Because the signature is
   unchanged, existing call sites keep working untouched.

   No array↔container wrapping happens anywhere in this shim: the
   containers arrive ready to use. Each dispatch arm either passes them
   straight to `$1_fortran` or converts them to their `_C` mirrors with
   `%to_c()` for the bridge call (inside `#ifdef _TIM` — see
   lessons.md §9 on why `%to_c` cannot appear unguarded).

   Body is one `select case (mode)` over
   `getenv_mode("<ENV_VAR>", default=TIMH_runFORTRAN)` with three arms:
   `TIMH_capture`, `TIMH_runAMREX` (inside `#ifdef _TIM`), and `case
   default` falling through to `$1_fortran`. Template: lessons.md §2.
   Capture-record naming: lessons.md §13.

### 6. Update host-module `use` statements
   Add the imports listed in lessons.md §14, deduplicating with what the
   module already has.

### 7. Confirm callers need no changes
   The shim's signature is identical to the pre-rename subroutine's
   (Step 5), so every existing call site remains valid by construction.
   Verify this rather than assuming: run
   `grep -rn "call[[:space:]]\+$1(" $0/src $0/config_src` and confirm each
   site's argument count and order still match. Do not edit any call site.

   If a call site *does* need changing, something has gone wrong — most
   likely the shim's signature drifted from the original in Step 5. Stop
   and fix the shim instead of rewriting callers.

   Any call-site work involving `alloc` / `copy2F` / `free` belongs to
   `convert_array_containers`, not here.

### 8. Relocate halo checks and CPU clocks
   Move any halo-sufficiency `MOM_error(FATAL,...)` from inside the
   kernel to the caller or the shim entry (lessons.md §15). Leave any
   pre-existing `cpu_clock_begin/end` at the caller, around the shim
   call — never inside the shim (lessons.md §16).

### 9. Verify
   Run the three-mode matrix in lessons.md §17. If only the Fortran
   shim + capture are being delivered, stop after CAPTURE verification
   and report that the C++ side of `<prefix>_$1_bridge` is the next
   deliverable.

### 10. Commit and push to `claude_<function-name>_bridge` *(gated by the global git-commit preference, overridable per run — see below)*
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
   the branch `claude_$1_bridge` (use the lowercased function name so
   different functions land on different branches and can be committed
   in parallel). Step 1 has already guaranteed HEAD is at
   `dev/turbo-debug`, so the new branch will be based on it. Stage every
   modified and newly created file, and commit with a message that
   briefly summarizes the changes (the wrapped subroutine name, the new
   bridge symbol, and the affected caller files) and explicitly notes
   that the work is co-authored by Claude. Use a HEREDOC so the trailer
   is preserved verbatim:

   ```
   BRANCH="claude_$(echo "$1" | tr '[:upper:]' '[:lower:]')_bridge"
   git checkout -B "$BRANCH"
   git add -A
   git commit -m "$(cat <<'EOF'
   <one-line summary of the wrapping work for $1>

   <optional 1–3 line body describing the bridge symbol, callers touched,
   and verification mode reached>

   Co-authored-by: Claude <noreply@anthropic.com>
   EOF
   )"
   git push -u origin "$BRANCH"
   ```

   If the push is rejected because the branch already exists upstream
   with unrelated history, stop and surface the conflict to the user
   rather than force-pushing.

## Hard rules

- Do not change the public name of `$1`, and do not change its argument
  list at all — the shim inherits the converted subroutine's container
  signature verbatim and must drop into existing call sites unchanged.
- Do not convert data structures. Array dummies are already containers
  before this skill runs; if they are not, stop and refer the user to
  `convert_array_containers`. Specifically: never change a dummy's type,
  never introduce a `Box_t`, never rewrite a loop nest, and never add
  `alloc` / `%view` / `copy2F` / `free` calls at a call site.
- Do not introduce a global mode variable; per-kernel env vars only.
- Do not gate capture mode behind `#ifdef _TIM`; only the AMReX arm is gated.
- Do not omit the `case default` arm.
- Do not pass default-kind `logical`/`integer` to a `bind(C)` routine —
  cast to `c_bool` / `c_int` (lessons.md §3.2).
- Do not call `c_loc(this%data(1))` on an unallocated container; allocate
  first.
- Do not re-implement `getenv_mode`, `already_recorded`, `mark_recorded`,
  or `io_recorder` — `use` them from `turbotmp_helperF`.
- Do not reuse a `kernel` string across kernels.

If something not covered here comes up, consult lessons.md §10
(recurring pitfalls) before improvising.

## Output to the user on success

Report:

1. The shim's env-var name and the three accepted values.
2. The `kernel` string used for capture filenames.
3. The bridge symbol name(s) the C++ side must implement.
4. The list of caller files that were rewritten.
5. What still needs to happen on the C++ side (if anything).
