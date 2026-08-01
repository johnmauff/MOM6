<!-- This file is intentionally a pointer, not a copy. Its content lives in
     ../cpp_bridge_lessons/SKILL.md. Do not re-inline that content here:
     the two copies drifted twice while both were maintained, and nothing
     reads this file's body any more. -->

# Lessons from PR #15 — moved

**The content that used to live here is now the `cpp_bridge_lessons`
skill:** [`../cpp_bridge_lessons/SKILL.md`](../cpp_bridge_lessons/SKILL.md).
Section numbering (§1–§17) is unchanged, so every "lessons.md §N"
citation in `generate_cpp_bridge/SKILL.md` still resolves — read §N of
`cpp_bridge_lessons` instead.

Invoke it once near the start of a session:

```
/cpp_bridge_lessons
```

## Why this file still exists

`generate_cpp_bridge/SKILL.md` Step 1 (under `--enable_src_validate`)
checks that
`<work-directory>/.claude/skills/generate_cpp_bridge/lessons.md` is
present, as a cheap proxy for "this checkout has the skill installed."
The check tests existence only, never content. Keeping a stub satisfies
it without maintaining a second copy of a 500-line document.

## If you are reading an older copy of this file

Some checkouts still carry the full pre-move text. Two things in it are
now wrong:

1. **`%dup` no longer exists.** It was restructured into `%alloc`. The
   old two-call recipe

   ```fortran
   call h_in_a%dup(h_in) ; call h_in_a%copy2Array(h_in)     ! WILL NOT COMPILE
   ```

   is now a single call:

   ```fortran
   call h_in_a%alloc(lb=LBOUND(h_in), ub=UBOUND(h_in), source=h_in)
   ```

   For a pure output, omit `source` rather than reading undefined
   memory: `call h_W_a%alloc(lb=LBOUND(h_W), ub=UBOUND(h_W))`.

2. **Converting data structures is no longer this skill's job.** Raw
   arrays → containers, `%view` in the body, and call-site
   `alloc`/`copy2F`/`free` all belong to the **`convert_array_containers`**
   skill, with its own reference in **`array_container_lessons`**. Run
   that first; `generate_cpp_bridge` then adds only the `*_fortran`
   rename, the dispatcher shim, capture mode, the `#ifdef _TIM` AMReX
   arm, and the `bind(C)` interface.

## Pre-condition (unchanged)

`generate_cpp_bridge` operates on a pre-existing TURBO-ESM/MOM6
checkout. The work directory must already contain the source tree (i.e.
have `src/` and `config_src/`) and be on (or rebased onto) the
`dev/turbo-debug` branch. Cloning is not performed — use

```
git clone -b dev/turbo-debug git@github.com:TURBO-ESM/MOM6.git <dir>
```

once to set up the directory, then pass it as `<work-directory>` on
every subsequent skill invocation.
