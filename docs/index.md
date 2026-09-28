# brille

brille finds the Brillouin zone of a crystal, and the smaller irreducible part
of it that symmetry makes unique, and interpolates your calculations over that
part. Evaluate an expensive model, such as phonon frequencies and
eigenvectors, only at the points of a grid over the irreducible zone, and
brille estimates it anywhere in reciprocal space, using the crystal's
symmetry to rotate eigenvectors where needed.

You can use brille to:

- get the space group and point group symmetry operations of a lattice;
- find the first Brillouin zone and an irreducible Brillouin zone of any
  lattice;
- build a grid over either zone, and refine it where your model needs more
  points;
- linearly interpolate data you provide at the grid points, at any point in
  reciprocal space.

```bash
python -m pip install brille
```

## Where to go next

<div class="grid cards" markdown>

- **Tutorials**

    ---

    Learn by doing: [interpolate with brille](tutorials/interpolation.md),
    then [build and refine a mesh](tutorials/mesh.md).

- **How-to guides**

    ---

    Get a task done: [install](how-to/install.md), [give a lattice its
    symmetry](how-to/symmetry.md), [switch from the trellis to the
    mesh](how-to/trellis-to-mesh.md), [control threads](how-to/threads.md).

- **Explanation**

    ---

    Understand the ideas: [lattices and Brillouin
    zones](explanation/background.md), [interpolation
    grids](explanation/grids.md), [the eigenvector phase
    convention](explanation/phase-convention.md), [memory
    use](explanation/memory-use.md).

- **Reference**

    ---

    Look up the API: [lattices](reference/lattice.md),
    [symmetry](reference/symmetry.md), [zones](reference/brillouin-zone.md),
    [grids](reference/grids.md), and more.

</div>

Worked examples that interpolate phonons from force constants, for NaCl with
Euphonic, are in the documentation of
[brilleu](https://github.com/brille/brilleu), which joins brille to
Euphonic.
