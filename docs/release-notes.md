# Release notes

## Unreleased

Changes since v0.8.3. The largest are a new mesh for `BZMeshQ`, which fills the
irreducible zone exactly and can be refined; eigenvectors of centred crystals
now describing the primitive cell; and many lattices that failed to mesh,
or to find their irreducible zone, now working.

### Changes that can affect existing code

#### Eigenvectors of centred crystals describe the primitive cell

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

**New for this.**

- [`Lattice.primitive_basis`][brille._brille.Lattice.primitive_basis]: the
  atoms a grid's eigenvectors describe.
- [`Lattice.centring_vectors`][brille._brille.Lattice.centring_vectors]: the
  cell's centring vectors.
- [`brille.utils.conventional_to_primitive`][brille.utils.conventional_to_primitive]:
  converts a conventional cell's eigen-solutions to its primitive cell.

[Give a lattice its symmetry](how-to/symmetry.md#eigenvectors-and-centred-cells)
describes the primitive basis and the conversion.

#### Eigenvectors in Cartesian units are normalized by default

Linear interpolation between unit eigenvectors gives vectors shorter than one
wherever neighbouring eigenvectors differ, so structure factors computed from
them came out too small. Interpolated eigenvectors stored in Cartesian units
(`LengthUnit` `angstrom` or `inverse_angstrom`, as Euphonic stores them) are
now scaled to unit norm by default; those in lattice units are not.
Structure factors computed from interpolated eigenvectors therefore change, for
the better. [`set_vector_normalization`][brille._brille.BZMeshQdc.set_vector_normalization]
turns it on or off, and takes a metric, such as η for spin-wave (Bogoliubov)
vectors.

#### `BZMeshQ` is a new, structured mesh

[`BZMeshQdc`][brille._brille.BZMeshQdc] and its siblings no longer use TetGen.
The mesh is now a grid of the reciprocal lattice clipped exactly to the
irreducible zone ([Interpolation grids](explanation/grids.md#the-structured-mesh-bzmeshq)).
The constructor's arguments are unchanged, and `max_size` is still the largest
tetrahedron volume in Å⁻³, but:

- the vertices differ, so interpolated values differ by about the interpolation
  error;
- `num_levels` is ignored;
- a mesh saved by an earlier version still loads, but cannot be refined.

The mesh matches itself across equivalent zone faces, so interpolation is
continuous there; interpolating with it is 2 to 5 times faster than with
[`BZTrellisQdc`][brille._brille.BZTrellisQdc]. `BZTrellisQ` and `BZNestQ` still
use TetGen. [Switch from BZTrellisQ to BZMeshQ](how-to/trellis-to-mesh.md)
compares the two.

#### `BRILLE_NUM_THREADS` replaces `OMP_NUM_THREADS`

brille no longer uses OpenMP; its parallel sections run on a pool of native
threads. `OMP_NUM_THREADS` no longer limits them: set `BRILLE_NUM_THREADS`
instead, for example when several processes use brille at once
([Control the number of threads](how-to/threads.md)).

#### Space group names and symmetry operations are checked

- A space group name is matched ignoring spaces and subscript marks, so
  `P 2/m`, `P21/c` and `P 21/c` now name the groups they should. Before, a
  name not spelt exactly as in brille's table was misread as a different
  group. A string that is not a name, a Hall symbol or CIF xyz operations
  now raises `ValueError`
  ([Give a lattice its symmetry](how-to/symmetry.md#by-name)).
- Symmetry operations that do not map the lattice onto itself (such as the
  Hall symbol `R 3 -2`, where `R 3 -2"` was meant) are refused with
  `ValueError` when the lattice is made, instead of failing later with no
  irreducible zone.

#### Other changes

- Python 3.11 or later is required.
- The printing functions `emit` and `emit_datetime` take their argument as
  `status`.

### New

- **Mesh refinement.** [`refinement_points`][brille._brille.BZMeshQdc.refinement_points]
  and [`refine`][brille._brille.BZMeshQdc.refine] add vertices where
  interpolation needs them, keeping existing vertices and their data, with a
  resolution limit; [`release_triangulation`][brille._brille.BZMeshQdc.release_triangulation]
  frees the memory refinement holds
  ([Building and refining a mesh](tutorials/mesh.md)).
- **The GIL is released** during zone and grid construction, interpolation,
  moving points into the zone, and mode sorting, so other Python threads run
  meanwhile, and several threads can use brille at once.
- **`NearSymmetryWarning`**: a lattice within 10⁻⁴ of a more symmetric one
  has a zone with features far smaller than itself; brille keeps the exact
  zone and warns (`BrillouinZone(..., warn_near_symmetry=False)` silences it).
- **A type stub**, `brille/_brille.pyi`, and `py.typed` give type checkers
  and editors the compiled module's signatures and docstrings.
- `Lattice.primitive_basis`, `Lattice.centring_vectors` and
  `brille.utils.conventional_to_primitive` (see above).

### Fixed

- With `time_reversal_symmetry=True`, crystals without inversion symmetry
  failed with "All atoms in the basis *must* be mapped…". Time reversal is now
  applied as complex conjugation of the eigenvectors.
- About 2 % of real lattices, mostly centred or nearly more symmetric, could
  not be meshed ("Duplicate intersection point").
- Some centred lattices given by a primitive cell with explicit symmetry
  operations, as phonon codes give them, had no irreducible zone
  (ITA 22, 197, 199, 203), and some meshes ignored `max_size` (ITA 127, 200,
  201, 215, 221). The irreducible wedge is now found exactly.
- Meshing a lattice within 6 ppm of face-centred cubic never finished.
- Interpolating at points that round-off leaves just outside the mesh, such as
  zone-surface points, raised "N of M points not found in tetrahedral mesh"
  (148 of 581 real lattices, at their own mesh vertices).
- An error inside a parallel section aborted Python; it is now raised as a
  Python exception.
- Parallel calls into TetGen could crash, or build different meshes; they are
  now serialised.
- `Lattice(..., symmetry="x,y,z;...")` failed, though documented.
- `import brille.vis` failed on Python 3.11 and later.
- `to_file` printed "Provided flags ... is translated to ...".

### Performance

- Moving points into the first and irreducible zones works on plain numbers
  and runs in parallel: interpolation on a small mesh went from 44 to 10 µs
  per point on one thread, and now gains from more.
- The structured mesh locates points in constant time and builds far faster
  than TetGen for large meshes.

### Documentation

The documentation is now a [Zensical](https://zensical.org) site arranged as
tutorials, how-to guides, explanation and reference, with new pages on the
structured mesh, switching grids, giving a lattice its symmetry, and the
eigenvector phase convention. Every page's example code is tested.
