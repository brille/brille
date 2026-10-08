# Give a lattice its symmetry

A lattice's symmetry decides its irreducible Brillouin zone, and so how much
of reciprocal space a grid must cover. [`brille.Lattice`][brille.lattice.Lattice]
takes it in any of four forms. The code here comes from
[`symmetry.py`](symmetry.py), which the tests run.

## By name

Give a Hall symbol or a Hermann-Mauguin symbol as `spacegroup`:

```python
--8<-- "how-to/symmetry.py:name"
```

Names are matched ignoring spaces and subscript marks, so `P 2_1/c`,
`P21/c` and `P 21/c` are the same. For a space group with more than one
setting, give the choice as well, `spacegroup=("P 2/m", "b")`, or the full
Hermann-Mauguin symbol that names it, `P 1 2/m 1`: here, the unique axis
along **b**. A string that is none of these is refused with a `ValueError`.

## By generators

A [`HallSymbol`][brille._brille.HallSymbol] decodes a Hall symbol into the
generators of its space group, and a lattice accepts those as `symmetry`:

```python
--8<-- "how-to/symmetry.py:generators"
```

```text
11 generators make 192 operations
```

## By explicit operations

Give every operation as a rotation matrix, in units of the real-space basis
vectors, and a translation, as fractions of them. This is the form that
[spglib](https://spglib.readthedocs.io) returns for a cell, so it is the
route for a primitive cell from a phonon calculation:

```python
--8<-- "how-to/symmetry.py:explicit"
```

## From a CIF file

A lattice, or a [`Symmetry`][brille._brille.Symmetry], reads operations written
as CIF `x, y, z` strings, separated by `;`; generators are enough:

```python
--8<-- "how-to/symmetry.py:cif"
```

## Check the result

All four routes give the same zone:

```python
--8<-- "how-to/symmetry.py:check"
```

```text
Hall symbol       48 point operations, irreducible zone 0.112207 Å⁻³
Hermann-Mauguin   48 point operations, irreducible zone 0.112207 Å⁻³
generators        48 point operations, irreducible zone 0.112207 Å⁻³
operations        48 point operations, irreducible zone 0.112207 Å⁻³
```

## Eigenvectors and centred cells

A grid's eigenvectors describe the atoms of one primitive cell,
[`primitive_basis`][brille._brille.Lattice.primitive_basis]. For a lattice given
by a centred conventional cell (A, B, C, F, I or R), that is the first atom of
each group of centring copies in the basis you give, a half, a third or a
quarter of it:

```python
--8<-- "how-to/symmetry.py:basis"
```

```text
8 atoms in the cell given, 2 in the primitive cell a grid's eigenvectors describe
```

Give the lattice the conventional cell's full basis, and the grid the
eigenvectors of the primitive cell, for those atoms in that order, in the
cell phase convention (see [the phase convention](../explanation/phase-convention.md)).

If your code computed the conventional cell instead, the grid also takes its
eigenvectors as they are, for every atom of the basis, and interpolates all of
that cell's modes. These are the primitive cell's modes at q folded together
with those at the wavevectors the larger cell cannot tell from q. To keep only
the primitive cell's modes, convert them first: the primitive cell's modes at q folded together with
[`conventional_to_primitive`][brille.utils.conventional_to_primitive] picks out
the modes at q and the primitive cell's atoms:

```python
from brille.utils import conventional_to_primitive

values, vectors = conventional_to_primitive(lattice, grid.rlu, values, vectors)
grid.fill(values, value_elements, vectors, vector_elements)
```

Filling a grid with the conventional cell's eigenvectors is refused when
interpolating, with an error that says so.
