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

If you will fill a grid with eigenvectors, the cell you describe must be the
cell your eigenvectors describe:

- **Conventional cell, by name:** for a centred lattice (A, B, C, F, I or R),
  brille builds the zone from the primitive cell but keeps the conventional
  cell's operations and atoms. Eigenvectors must then cover every atom of the
  conventional cell, two to four times as many as the primitive cell has.
- **Primitive cell, by explicit operations:** phonon codes give eigenvectors
  for the primitive cell. Give brille the primitive lattice vectors and their
  operations, from spglib for example.

Eigenvalues alone (energies, for instance) work either way.
