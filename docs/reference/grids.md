# Grids

Each grid comes in three types, by whether its values and vectors are real
(`d`) or complex (`c`): `dd`, `dc` and `cc`. They have the same methods, so
only the `dc` type, the usual one for phonons, is shown in full.
[Interpolation grids](../explanation/grids.md) explains how the grids differ.

## Structured mesh

::: brille._brille.BZMeshQdc

::: brille._brille.BZMeshQdd
    options:
      members: false

::: brille._brille.BZMeshQcc
    options:
      members: false

## Hybrid trellis

::: brille._brille.BZTrellisQdc

::: brille._brille.BZTrellisQdd
    options:
      members: false

::: brille._brille.BZTrellisQcc
    options:
      members: false

## Nested mesh

::: brille._brille.BZNestQdc

::: brille._brille.BZNestQdd
    options:
      members: false

::: brille._brille.BZNestQcc
    options:
      members: false

## Helpers

::: brille._brille.RotatesLike

::: brille._brille.SortingStatus
