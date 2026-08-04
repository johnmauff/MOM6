# Converting a MOM6 subroutine call tree to AMReX-style Fortran

This document describes how to apply a series of Claude skills to convert a
subroutine call tree from MOM6-style Fortran to AMReX-style Fortran. The
AMReX-style Fortran is constructed so that it is easy for another Claude
skill to convert into AMReX code.

## 1. Choose a root subroutine

To convert a subroutine call tree, first identify a target root. Look for a
root that has a non-trivial amount of work in its descendant subroutines. It
is also useful if the root is fundamentally isolated from the rest of the
model, rather than fanning out into a widely shared subsystem.

Descriptions of candidate roots can be found in `barotropic_call_tree_audit.md`,
`hor_visc_call_tree_audit.md`, and `vert_friction_call_tree_audit.md`. The
rest of this document uses the barotropic solver (`btstep`) as a worked
example.

## 2. Load the conversion context

Start Claude in agentic mode and invoke the `array_container_lessons` skill.
This sets the context for the rest of the transformation and only needs to
be done once per session.

```
/array_container_lessons
```

## 3. Convert optional arguments at the root

Next, run `convert_optional_args_to_containers` at the top of the call tree:

```
/convert_optional_args_to_containers . btstep
```

This converts all optional array arguments into optional containers. It is
critical to do this first: converting the call tree bottom-up instead would
force enumerating every present/absent combination of optional arguments at
each call site into the next level down — a 2^n-way branch, where n is the
number of optional array arguments involved.

## 4. Convert each descendant, top-down

Run `convert_array_containers` on each descendant of the root, one level at
a time. For the barotropic solver:

```
/convert_array_containers . BT_cont_to_face_areas
/convert_array_containers . adjust_local_BT_cont_types
...
/convert_array_containers . btstep_timeloop
```

Only after `btstep_timeloop` has been converted should you progress to its
own descendants. Convert the tree one level at a time, top-down.

## 5. Clean up

Once the entire call tree has been converted, apply the cleanup skills.

**Cleanup of the call site** minimizes the number of `%alloc` and `%free`
calls used:

```
/hoist_container_marshalling . btstep
```

**Cleanup of the descendants:**

```
/convert_locals_to_containers . BT_cont_to_face_areas
```

Note that this last skill has not yet been tested, so it is unclear whether
it will actually perform the necessary cleanup.

## Known issues

- There is currently no logical-array container. If the subroutine call tree
  contains logical arrays, the raw arrays must be passed through to the next
  level instead of being containerized.

- One of the goals of this transformation pass is to eliminate any presence
  of the `G` and `GV` data structures from converted subroutine interfaces,
  replacing them with the primitive scalars and containers a subroutine
  actually uses. In the continuity solver, this could not be fully achieved:
  a lower-level leaf subroutine still declared `G` and `GV` dummies it never
  actually referenced internally. Removing those two unused dummies from the
  leaf would be a real, if small, win on its own — but it would *not*
  eliminate `G`/`GV` further up the call tree, since callers still need them
  for sizing local arrays (via the `SZI_`/`SZJ_`/`SZK_` macros) and for
  extracting per-point grid metrics that feed into the leaf's own scalar
  arguments.
