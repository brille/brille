.. _howto_trellis_to_mesh:

=================================
Switch from BZTrellisQ to BZMeshQ
=================================

:py:class:`~brille._brille.BZMeshQdc` and its siblings can replace
:py:class:`~brille._brille.BZTrellisQdc` in most code: they are filled and
interpolated in the same way. This page shows what to change, what to watch
for, and how to check whether the mesh suits your system.

The code on this page comes from :download:`trellis_to_mesh.py`, which the
tests run.

Change the constructor
----------------------

A trellis like this:

.. literalinclude:: trellis_to_mesh.py
  :start-after: # [before]
  :end-before: # [/before]

becomes:

.. literalinclude:: trellis_to_mesh.py
  :start-after: # [after]
  :end-before: # [/after]

Both size parameters are volumes in Å⁻³, not fractions:

- The trellis's ``node_volume_fraction`` is, despite its name, the volume of
  one cube of the trellis.
- The mesh's ``max_size`` is the largest volume of one tetrahedron.

A cell of the mesh's grid holds six tetrahedra, so dividing by six gives a mesh
with about as many vertices as the trellis:

.. list-table:: Vertices for ``node_volume_fraction = V/N`` and ``max_size = V/(6N)``, where ``V`` is the irreducible zone's volume
  :header-rows: 1

  * - Lattice
    - N
    - Trellis
    - Mesh
  * - triclinic ``-P 1``
    - 10,000
    - 14,139
    - 18,215
  * - hexagonal ``-P 6``
    - 10,000
    - 13,248
    - 14,641
  * - fcc ``-F 4 2 3``
    - 10,000
    - 14,572
    - 14,840

All other arguments, and :py:meth:`~brille._brille.BZMeshQdc.fill`,
:py:meth:`~brille._brille.BZMeshQdc.ir_interpolate_at`,
:py:meth:`~brille._brille.BZMeshQdc.sort`, and saving and loading, work as
before.

If you create grids with :py:func:`brille.utils.create_grid`, pass
``mesh=True`` and ``max_size`` in place of ``node_volume_fraction``:

.. code-block:: python

  grid = create_grid(bz, complex_vectors=True, mesh=True, max_size=volume / (6 * points))

Through Euphonic
----------------

Euphonic's ``BrilleInterpolator.from_force_constants`` takes
``grid_type="mesh"``. For the same ``grid_npts``, it asks the mesh for
between a fifth and a half as many points as the trellis (in Euphonic 2.1). To keep the
number of points, give the mesh's size directly:

.. code-block:: python

  from brille import BrillouinZone
  from euphonic.brille import BrilleInterpolator

  volume = BrillouinZone(lattice).ir_polyhedron.volume   # your crystal's lattice
  interpolator = BrilleInterpolator.from_force_constants(
      force_constants, grid_type="mesh", grid_kwargs={"max_size": volume / (6 * grid_npts)})

What else changes
-----------------

.. list-table::
  :header-rows: 1

  * - Only on :py:class:`~brille._brille.BZTrellisQdc`
    - Instead, on :py:class:`~brille._brille.BZMeshQdc`
  * - ``always_triangulate``
    - Not needed: the mesh is always tetrahedra.
  * - ``interpolate_at``, which moves points into the first zone only
    - :py:meth:`~brille._brille.BZMeshQdc.ir_interpolate_at`
  * - ``inner_rlu``, ``outer_rlu``, ``node_at`` and other node methods
    - None: the mesh has no cubic nodes.

The mesh adds:

- :py:meth:`~brille._brille.BZMeshQdc.refinement_points` and
  :py:meth:`~brille._brille.BZMeshQdc.refine`, to add vertices where the
  interpolation needs them (see :ref:`tutorial_mesh`);
- ``max_points``, which caps the number of vertices.

Interpolated values change slightly, because the vertices are different: by
about the interpolation error, not more. On the zone boundary, the mesh
matches itself across equivalent faces, so interpolation is continuous there.

A saved trellis can only be read back as a trellis, so rebuild your grids as
meshes rather than converting files. Meshes saved by earlier versions of brille,
which built them with TetGen, still load, but cannot be refined.

Trade-offs
----------

These figures compare the two grids for three lattices. They were measured on
one machine, with one thread, interpolating a smooth symmetric function at
100,000 random points. Accuracy is given as rms error × vertices\ :sup:`2/3`:
for linear interpolation of a smooth function this hardly depends on the
number of vertices, so it compares grids of different sizes, and lower is
better.

.. list-table::
  :header-rows: 1

  * - Lattice
    - Grid
    - Vertices
    - Build / s
    - Interpolation / µs per point
    - Error × N\ :sup:`2/3`
  * - triclinic
    - trellis
    - 13,555
    - 7.56
    - 7.0
    - 5.91
  * -
    - mesh
    - 9,947
    - 0.44
    - 2.2
    - 5.96
  * - hexagonal 6/m
    - trellis
    - 4,008
    - 1.74
    - 5.9
    - 1.96
  * -
    - mesh
    - 8,141
    - 0.36
    - 2.2
    - 1.93
  * - fcc
    - trellis
    - 926
    - 0.46
    - 6.6
    - 1.51
  * -
    - mesh
    - 8,160
    - 0.66
    - 3.2
    - 1.19

In short:

- **Interpolation** with the mesh is 2–3 times faster: it finds the tetrahedron
  holding a point in constant time.
- **Building** the mesh is faster, by a factor that grows with its size:
  17 times for the triclinic lattice above.
- **Accuracy per vertex** is the same for fine grids, or better for the mesh
  (fcc above). For coarse grids of low-symmetry lattices, a trellis of a few
  hundred vertices can be up to twice as accurate per vertex as a mesh of the
  same size.
- **Memory** held by the mesh is similar to the trellis's, or less for the
  same number of vertices.

Stay with the trellis if you rely on its node methods or on ``interpolate_at``,
if you need results identical to earlier runs, or if you use very coarse grids
of low-symmetry lattices. Otherwise, the mesh is faster, and refinement lets it
put vertices where your model needs them.

Measure for your own system
---------------------------

The best test is your own model. This builds both grids at three sizes, then
reports build time, interpolation time and rms error; replace ``lattice`` and
``model`` in the script with yours.

.. literalinclude:: trellis_to_mesh.py
  :start-after: # [compare]
  :end-before: # [/compare]

For the example's triclinic lattice:

.. code-block:: text

  grid     vertices  build/s  µs/point  rms error
  trellis      1941     0.71      6.19   2.61e-02
  mesh         2886     0.48      1.55   2.53e-02
  trellis      7787     3.10      6.61   9.11e-03
  mesh         9947     1.12      2.15   9.06e-03
  trellis     26432    12.76      7.30   3.86e-03
  mesh        32505     2.76      3.11   3.59e-03
