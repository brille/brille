# Symmetry

brille holds each symmetry operation as an integer matrix, the identity, a
rotation or a rotoinversion in units of the real-space basis vectors, and a
vector, its translation as fractions of those basis vectors. A
[`Symmetry`][brille._brille.Symmetry] can also be made from operations written
as [CIF xyz](https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Ispace_group_symop_operation_xyz.html)
strings, and a [`HallSymbol`][brille._brille.HallSymbol] decodes a
[Hall symbol](https://doi.org/10.1107/S0567739481001228) into the generators of
its space group. [Give a lattice its symmetry](../how-to/symmetry.md) shows each
in use.

::: brille._brille.Symmetry

::: brille._brille.PointSymmetry

::: brille._brille.HallSymbol

::: brille._brille.Spacegroup

::: brille._brille.Pointgroup
