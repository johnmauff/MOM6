---
name: generate_cpp_bridge
description: Wrap an existing MOM6 Fortran subroutine in a runtime-dispatched shim that selects between (a) the original Fortran code, (b) a binary capture mode that records inputs+outputs to disk for offline validation, and (c) a C++/AMReX bridge invoked through bind(C). Use when porting any MOM6 kernel to AMReX while keeping the Fortran caller unchanged and the Fortran truth available as a numerical reference. Mirrors the pattern established in TURBO-ESM/MOM6 PR #15.
user-invocable: true
argument-hint: <work-directory> <function-name> [--enable_src_validate] [--enable_git_commit] [--disable_git_commit] [--change-shim-interface]
---

# Generate C++ bridge for a MOM6 Fortran subroutine

This skill is the **execution checklist**. All templates, rationale,
type-mapping tables, and pitfalls live in the `cpp_bridge_lessons` skill
(same numbered sections, §1–§17) — invoke `cpp_bridge_lessons` once near
the start of a session, before running this skill, and refer back to its
numbered sections from each step below (Step 0 checks that it has been
loaded). Do not reproduce templates here.

## Help message

If `$ARGUMENTS` is empty, or equals `help`, or equals `--help`, or equals
`-h`, do NOT run any steps. Print the following help message verbatim
and stop:

```
Usage: /generate_cpp_bridge <work-directory> <function-name> [--enable_src_validate] [--enable_git_commit] [--disable_git_commit] [--change-shim-interface]

Wrap a MOM6 Fortran subroutine in a runtime-dispatched shim that selects
between the original Fortran code, a capture mode for offline validation,
and a C++/AMReX bridge. Mirrors TURBO-ESM/MOM6 PR #15.

Arguments:
  <work-directory>   Absolute path to an existing TURBO-ESM/MOM6 checkout
                     (must contain src/ and config_src/, and should be on
                     the dev/turbo-debug branch). Cloning is not performed.
  <function-name>    Name of the Fortran subroutine to wrap (case-insensitive
                     match against its declaration in the tree).
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
  --change-shim-interface  (optional) Allow the shim subroutine's dummy argument
                           list to change: array dummies become RealArray_t /
                           IntArray_t containers in the public interface.
                           When set, callers are rewritten to build containers
                           and call the new interface. When omitted (default),
                           the shim preserves the original argument list and
                           performs the array↔container conversion internally,
                           leaving call sites unchanged.

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
   `--enable_git_commit`, `--disable_git_commit`, and `--change-shim-interface`.
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

Step 1 runs only when `--enable_src_validate` is set. Item 4 (reference-
material priming) always runs, regardless of flags — this is what replaces
the old flag-gated read of lessons.md. Step 10's behavior is decided by
`--enable_git_commit` / `--disable_git_commit` if passed, or otherwise by
the global preference described in Step 10. No other wrapping work executes
until Step 0 validation passes.

## Settle these decisions (ask if not obvious from the tree)

1. **Bridge prefix** — default to whatever existing `_bridge) bind(C)`
   declarations in the tree use; otherwise `turbotmp_`.
2. **Env-var name** — default `<UPPERCASE_$1>_MODE`.
3. **AMReX gating** — default `#ifdef _TIM` around the `case (TIMH_runAMREX)`
   arm only; capture mode is never gated.
4. **Shim interface mode** — if `--change-shim-interface` was passed, the shim's
   public argument list will use `RealArray_t`/`IntArray_t` containers and
   callers will be rewritten. Otherwise (default) the shim preserves the
   original array-based argument list and converts internally.

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

### 2. Classify each dummy argument
   For every arg, record intent / kind / rank / optional. Map to bridge
   types using lessons.md §3.1–§3.2. If the original takes loop-bound
   integers, collapse them into a single `Box_t` (lessons.md §5).

### 3. Rename the original implementation
   Rename `$1` → `$1_fortran` in place; convert array dummies to
   `RealArray_t` / `IntArray_t`; rewrite loops as `do concurrent` over
   `bx%idxS`/`bx%idxE`. Template: lessons.md §12.
   `$1_fortran` always uses the container-based signature regardless of the
   shim interface mode — only the public-facing shim (`$1`) differs.

### 4. Add the `bind(C)` interface block
   At the top of the host module, declare `<prefix>_$1_bridge` per the
   template in lessons.md §3.3. Doc-comment every dummy with `!<`.

### 5. Write the shim subroutine `$1`
   Always uses the same public name as the original. The dummy argument
   list depends on the shim interface mode:

   - **Default (interface preserved):** Keep the original array-based dummy
     list. In each dispatch arm, wrap raw arrays into `RealArray_t` /
     `IntArray_t` containers before calling `$1_fortran` or the bridge, and
     unwrap them on return. Call sites in existing callers require no changes.

   - **`--change-shim-interface` (interface changed):** Replace array dummies
     with `RealArray_t` / `IntArray_t` containers in the public signature.
     No internal wrapping is needed in the dispatch arms. Callers must be
     rewritten (Step 7).

   Body is one `select case (mode)` over
   `getenv_mode("<ENV_VAR>", default=TIMH_runFORTRAN)` with three arms:
   `TIMH_capture`, `TIMH_runAMREX` (inside `#ifdef _TIM`), and `case
   default` falling through to `$1_fortran`. Template: lessons.md §2.
   Capture-record naming: lessons.md §13.

### 6. Update host-module `use` statements
   Add the imports listed in lessons.md §14, deduplicating with what the
   module already has.

### 7. Rewrite each caller
   - **Default (interface preserved):** The shim accepts the same raw arrays
     as the original; call sites are already compatible. No caller rewriting
     is required. Proceed to Step 8.
   - **`--change-shim-interface` (interface changed):** Wrap raw arrays into
     `RealArray_t` containers, build the iteration `Box_t`, call the shim,
     copy results back, free. Recipe: lessons.md §4. For `intent(out)` args,
     skip the inbound `copy2Array`; for pure-input args, skip the outbound
     `copy2F`. Reuse adjacent-kernel containers if the caller already has them.

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

- Do not change the public name of `$1`. Do not change its argument list
  unless `--change-shim-interface` was explicitly passed — if it was not,
  the shim must drop into existing call sites unchanged.
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
