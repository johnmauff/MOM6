---
name: convert_present_to_associated
description: Convert an already-array-container-ified optional dummy argument from Fortran's `optional`/`present()` idiom to a non-optional container checked via `associated(x%data)`. Use this to make an optional argument bridge-ready for eventual bind(C)/C++ dispatch, since C/C++ has no `optional` concept -- only nullable pointers. Hard precondition: the argument's declared type at THIS call site must already be RealArray_t/IntArray_t (run convert_array_containers first if it is still a raw array). Does not touch the container API itself, does not add a bind(C) bridge -- that is generate_cpp_bridge's job.
user-invocable: true
argument-hint: <work-directory> <function-name> <dummy-argument-name>[,<dummy-argument-name>...]
---

# Convert `present(x)` to `associated(x%data)`

## Why this exists

Fortran's `optional` attribute has no C/C++ equivalent — `bind(C)`
interfaces can't express "this argument may be omitted." The bridging
alternative is to make the argument **always present as a container**,
whose internal pointer (`RealArray_t%data` / `IntArray_t%data`) may be
unassociated. "Absent" becomes "associated with nothing" instead of
"not passed." This skill performs that transformation on the Fortran
side only, one call level at a time.

This is intentionally a **separate skill from `convert_array_containers`**.
Converting a raw array to a container and converting an optional
argument's presence-test idiom are independent transformations that
happen to often apply to the same variable at different times.

## Hard precondition — check this before doing anything else

The dummy argument named in `<dummy-argument-name>` must **already** be
declared as `type(RealArray_t)` or `type(IntArray_t)` (with `optional`)
at this call site. If it is still a raw array
(`real, dimension(...), optional`), **stop** — run
`convert_array_containers` on this subroutine first. This skill does
not perform the raw-array-to-container conversion itself.

Also check the immediate caller(s): every call site that currently
omits this argument (relying on it being `optional`) will, after this
change, need to pass an explicit container instead (see Step 3). If a
caller's own corresponding argument is itself still `optional`/raw,
that caller is out of scope for this pass — convert bottom-up, one
fully-consistent call-level boundary at a time. Do not leave a
half-converted boundary where the callee expects a container but a
caller is still omitting the argument.

## Grouped arguments

If the target argument's presence is checked together with one or more
other optional arguments (e.g. `if (present(a) .and. present(b))`, or
the routine's logic treats a set of optional arguments as
all-present-or-all-absent), convert the whole group in the same pass.
Converting one member of such a group while leaving siblings on the old
`present()` idiom produces inconsistent semantics — the routine can end
up with `associated(a%data)` true and `present(b)` false for what was
meant to be a single logical state.

## Steps

1. **Confirm the precondition and grouping** above. List every
   `present(<name>)` occurrence for the target argument(s) in the
   subroutine body.

2. **Drop `optional` from the dummy argument declaration.** The
   container's own `data` pointer (default `=> null()`) now carries the
   "is this here" signal — no separate `optional` flag is needed.

   ```fortran
   ! before
   type(RealArray_t), optional, intent(in) :: uhbt
   ! after
   type(RealArray_t), intent(in) :: uhbt
   ```

3. **Replace every `present(x)` with `associated(x%data)`** in the
   subroutine body. Do not just wrap the old checks — read each one to
   confirm it is testing the same thing `associated()` would test
   (presence of the container's data), not some unrelated flag that
   happens to share a name.

   ```fortran
   ! before
   if (present(uhbt)) then
   ! after
   if (associated(uhbt%data)) then
   ```

4. **Fix every call site.** A call site that previously omitted this
   argument (legal only because it was `optional`) must now pass an
   explicit container — typically one that was never `%alloc`'d, so its
   `data` pointer stays `null()` and `associated()` correctly reads as
   absent downstream. A call site that already builds the container
   unconditionally needs no change here.

5. **Leave `%to_c()` alone.** TIM's `to_c_Real`/`to_c_Int` are already
   null-safe (an unassociated container converts to an all-null
   `cdesc`, `rank=0` — see `cpp_bridge_lessons` §3.1/§10). This skill
   does not need to touch bridge code; it only changes how "presence"
   is represented on the Fortran side before any bridge exists.

## Explicitly NOT performance motivated

`associated()` and `present()` are both O(1) runtime checks — Fortran
does not resolve either at compile time. Do not justify this
conversion, or decline it, on performance grounds in either direction;
the only reason to do it is bridge-readiness.

## Hard rules

- Never convert an argument whose declared type at this call site is
  still a raw array — that is `convert_array_containers`' job, run it
  first.
- Never convert one member of a present()-linked group without the
  rest of the group in the same pass.
- Never leave a caller passing nothing (relying on `optional`) after
  the callee has dropped `optional` — every call site must be updated
  in the same pass as the callee.
- Do not add a bind(C) interface, dispatcher, or capture mode as part
  of this skill — that is `generate_cpp_bridge`'s job, and strictly
  comes after this conversion, not as part of it.
