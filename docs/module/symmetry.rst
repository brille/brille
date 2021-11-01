Symmetry
--------
A set of symmetry operations can be represented in various different ways.
Within :py:mod:`brille` each operation is represented as an integer-valued matrix and a floating-point-valued vector.
The matrix must be the identity matrix, a rotation matrix, or a rotoinversion matrix expressed as integer multiples of the real-space basis vectors.
The vector is the translation part of the operation expressed in fractional coordinates of the same basis vectors.

Other representations of a symmetry operation are convenient for human interpretation and communication.
One such which can be used to specify symmetry operations for :py:class:`~brille._brille.Symmetry` is `CIF xyz`_,
which allows for specifying either a list of matrices and list of vectors or a single CIF xyz string:

>>> s0 = brille._brille.Symmetry([[[1, 0, 0], [0, -1, 0], [0, 0, 1]]], [[0, 1/2, 1/2]])
>>> s1 = brille._brille.Symmetry('x, 1/2-y, 1/2+z')
>>> s0 == s1
True

.. _`CIF xyz`: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Ispace_group_symop_operation_xyz.html


With a CIF xyz encoded string from, e.g., a CIF file, one can directly create a :py:class:`~brille._brille.Direct` lattice object
with the appropriate symmetry operations without also specifying, e.g., the spacegroup::

  basis_lengths, basis_angles, xyz_string = hypothetical_CIF_parser('some_system.cif')
  lattice = brille._brille.Direct(basis_lengths, basis_angles)
  lattice.spacegroup = brille._brille.Symmetry(xyz_string)

in setting the lattice spacegroup, the :py:meth:`~brille._brille.Symmetry.generate()` method is called,
so the provided CIF xyz operators need only be the generators of the spacegroup.


.. autoclass:: brille._brille.Symmetry
  :members:

Hall Symbol
===========
A more compact representation of a full space group is the Hall_ symbol_
which has the advantage of representing all unique crystallographic lattices as an encoded set of generators.

.. _Hall: https://doi.org/10.1107/S0567739481001228
.. _symbol: http://cci.lbl.gov/sginfo/hall_symbols.html

Within :py:mod:`brille`, :py:class:`~brille._brille.HallSymbol` can be used to decode an arbitrary Hall symbol.
And then its generators can be extracted to create the spacegroup for lattice::

  basis_lengths, basis_angles, hall_symbol = hypothetical_CIF_parser('some_system.cif')
  lattice = brille._brille.Direct(basis_lengths, basis_angles)
  lattice.spacegroup = brille._brille.HallSymbol(hall_symbol).generators


.. autoclass:: brille._brille.HallSymbol
  :members:


Other classes
=============

.. autoclass:: brille._brille.PointSymmetry
  :members:
