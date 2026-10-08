# Release notes

## 0.9.1 (unreleased)

### Fixes

#### Eigenvectors of a centred conventional cell are accepted again

0.9.0 refused eigenvectors of a centred crystal's conventional cell with
`RuntimeError: The eigenvectors have 8 atoms per mode, but the primitive cell
has 2`, though brille up to 0.8.3 rotated them correctly. This broke Euphonic's
`BrilleInterpolator` for force constants of a conventional cell, such as
NaCl's 8-atom cubic cell. A grid now takes eigenvectors for either the atoms of
[`Lattice.primitive_basis`][brille._brille.Lattice.primitive_basis] or every
atom of the conventional cell's basis, and interpolates the modes of the cell
they describe. Other atom counts are still refused.

#### Compact irreducible zones for point groups without mirrors

0.9.0 found the irreducible wedge as the Dirichlet cone about one fixed point.
For a point group without mirrors (such as 32, 422, 222, -3 or 432) its planes
tilt with that point, so the irreducible zone had more vertices than needed:
13 for quartz (P3₂21), not the 6 of the Γ–K–K′ prism found by 0.8.3. Meshes
of such zones had more vertices and were less accurate for their size. The
wedge now comes from whichever of several points gives the zone with the
fewest vertices, which on the AFLOW lattices is never more than in 0.9.0 or
0.8.3. Point groups with mirrors keep their 0.9.0 wedge. Grids of the other
point groups have different vertices, so interpolated values change by about
the interpolation error.

## 0.9.0

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

#### `create_grid` makes a `BZMeshQ` by default

[`create_grid`][brille.utils.create_grid] now makes a `BZMeshQ` unless you ask
for another grid: with `trellis=True` (new), `nest=True`, or a trellis's own
arguments (`node_volume_fraction`, `always_triangulate`), so calls that size a
trellis still get one. Without a size, the mesh has `max_size = 1e-5 / 6`,
about as many vertices as the trellis's default; `mesh=True` without a size
used to give the coarsest mesh.

Euphonic makes its own grids, and its `BrilleInterpolator` still defaults to
`grid_type="trellis"`. The mesh is now the recommended grid: pass
`grid_type="mesh"`, and size it with
`grid_kwargs={"max_size": volume / (6 * grid_npts)}`, since for the same
`grid_npts` Euphonic 2.1 asks the mesh for fewer points than the trellis
([Switch from BZTrellisQ to BZMeshQ](how-to/trellis-to-mesh.md#through-euphonic)).

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

#### Hall numbers are gone from the API

Hall numbers only number the rows of spglib's table of settings; symbols name
space groups. The constructors that took them are removed:
`Symmetry(hall_number)`, `PointSymmetry(hall_number, time_reversal)` (already
deprecated), `Spacegroup(hall_number)` and its `hall_number` property, and
`PrimitiveTransform(hall_number)`. Instead:

- [`Spacegroup`][brille._brille.Spacegroup] takes a Hall symbol, or a
  Hermann-Mauguin symbol or International Tables name with an optional
  setting choice, as `Lattice` does, and `Spacegroup.all()` lists every setting;
- `PointSymmetry(Symmetry(...))` or `Lattice.pointgroup` give point groups;
- [`PrimitiveTransform`][brille._brille.PrimitiveTransform] takes a
  [`Bravais`][brille._brille.Bravais] centring type.

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
- `Spacegroup(symbol, choice)` and `Spacegroup.all()`, and `create_grid`'s
  `trellis` flag (see above).

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
- Vector or matrix data in units whose rotation is not implemented (such as
  Cartesian vectors, or Γ eigenvectors in reciprocal-lattice units) were
  accepted by `fill` and failed only when interpolating; `fill` now raises
  `ValueError`, naming the supported combinations. Scalars need no units.
- Building a `BZTrellisQ` could corrupt it: a worker adding a vertex grew the
  shared vertex array while others read it. Its polyhedron nodes are now
  triangulated by one task, which also makes it about twice as fast on many
  threads.
- The API reference listed keyword arguments, and a few names, as parameters
  the functions do not have.

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
