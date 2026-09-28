# Release notes

## Unreleased

### Eigenvectors of centred crystals describe the primitive cell

!!! warning "Breaking change"
    If you give brille a centred crystal (A, B, C, F, I or R lattice) by its
    conventional cell and fill a grid with eigenvectors, those eigenvectors
    must now describe the atoms of the **primitive** cell, not the
    conventional cell.

**Who is affected.** You are affected if you both:

- give a lattice by its centred conventional cell, whether by space group
  name, Hall symbol, or symmetry operations that include centring
  translations; and
- fill a grid with eigenvectors (`RotatesLike.Gamma`) computed for that
  conventional cell.

Primitive lattices, and grids whose vectors are not phonon eigenvectors, are
unaffected. Euphonic's `BrilleInterpolator` gives brille the cell of its force
constants: most force constants are for a primitive cell, and are unaffected;
force constants for a centred conventional cell are affected.

**Why.** brille built a centred crystal's Brillouin zone from its primitive
cell, but described its eigenvectors by the conventional cell. The
conventional cell has two to four times the atoms, and its modes at each q
include those of other wavevectors folded in. Phonon codes give eigenvectors
of the primitive cell, which could not be used on this route before.

**What changes.**

- A grid's eigenvectors describe the atoms of
  [`Lattice.primitive_basis`][brille._brille.Lattice.primitive_basis]: the
  first atom given of each set of centring copies, in the order given. For a
  primitive lattice this is the whole basis, as before.
- The lattice still takes the conventional cell's full basis.
- Filling a grid with the conventional cell's eigenvectors now fails when
  interpolating, with

    ```text
    RuntimeError: The eigenvectors have 8 atoms per mode, but the primitive cell has 2. Give eigenvectors of the primitive cell; for those of a centred conventional cell, convert them with brille.utils.conventional_to_primitive.
    ```

    Before, such eigenvectors were accepted.

**What to do.** Either compute eigenvectors of the primitive cell, for the
atoms of `lattice.primitive_basis`, or convert those of the conventional
cell with [`conventional_to_primitive`][brille.utils.conventional_to_primitive]
before filling the grid:

```python
from brille.utils import conventional_to_primitive

values, vectors = conventional_to_primitive(lattice, grid.rlu, values, vectors)
grid.fill(values, value_elements, vectors, vector_elements)
```

with `vector_elements` giving three components for each atom of the
primitive cell. The eigenvectors must be in the cell phase convention (see
[the phase convention](explanation/phase-convention.md)). The interpolated
results then have the primitive cell's modes: a half, a third or a quarter
as many as before.

**New.**

- [`Lattice.primitive_basis`][brille._brille.Lattice.primitive_basis]: the
  atoms a grid's eigenvectors describe.
- [`Lattice.centring_vectors`][brille._brille.Lattice.centring_vectors]: the
  cell's centring vectors.
- [`brille.utils.conventional_to_primitive`][brille.utils.conventional_to_primitive]:
  converts a conventional cell's eigen-solutions to its primitive cell.

[Give a lattice its symmetry](how-to/symmetry.md#eigenvectors-and-centred-cells)
describes the primitive basis and the conversion.
