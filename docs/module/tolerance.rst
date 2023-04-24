Floating point tolerance
------------------------
Internally :py:mod:`brille` uses a mixture of 'exact' geometric predicates
from Jonathan Richard Shewchuk's `Fast Robust Predicates for Computational Geometry <http://www.cs.cmu.edu/~quake/robust.html>`_
via `TetGen <https://wias-berlin.de/software/index.jsp?id=TetGen>`_, and approximate floating point comparisons.

The latter comparisons were are intended to allow for errors in the floating point
representation of user-provided input *and* compounded rounding errors in geometric calculations.
The internal floating point tolerance check compares the difference between two numbers to *both*
the minimum floating point representation *difference*, so-called machine epsilon, :math:`\epsilon`,
scaled up by the sum of the two values times a constant,
*and* an absolute floating point tolerance.

It has been found through experience that a single sensible default floating point tolerance does not exist
for general user-provided input.
And for particularly egregious cases the precision to which, e.g., the positions of atoms within a unit cell
is provided may be significantly different than the same-lattice basis vectors.
This difference can lead to a need for different absolute floating point tolerances for different operations within
:py:mod:`brille._brille`.
In the specific case of the atom positions being specified from, e.g., 32-bit floating point values written to
a text file, it may be sufficient to specify a reduced absolute floating point tolerance only for real space calculations.
For this reason, it is possible to change the absolute floating point tolerance for real and reciprocal space
comparisons independently.

At present the absolute floating point tolerances for real and reciprocal space comparisons has been set to
:math:`10^{-12}` angstrom and inverse angstrom, respectively.
If these values are found lacking, they can be modified at runtime via, e.g., as a global setting for the module:

>>> import brille
>>> brille.real_space_tolerance(1e-10)
>>> brille.reciprocal_space_tolerance(1e-14)

or as a local setting, to be passed to a specific function,

>>> ac = brille.ApproxConfig()
>>> ac.real_space_tolerance = 1e-10
>>> ac.reciprocal_space_tolerance = 1e-14


.. note::
  It is hoped that these functions and methods are **no longer necessary** for typical input and typical users.
  Strange errors may occur when they are modified.
  Attempt to make use of, e.g., the `snap_to_symmetry`, functionality in :py:class:`brille._brille.Lattice` instead.


.. autofunction:: brille._brille.real_space_tolerance

.. autofunction:: brille._brille.reciprocal_space_tolerance

.. autoclass:: brille._brille.ApproxConfig
  :members:
