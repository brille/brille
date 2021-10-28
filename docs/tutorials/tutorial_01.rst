===============================
Specifying a lattice Spacegroup
===============================

The symmetry information for a lattice can be specified by providing a string at
construction which represents a Hall symbol or an International Tables
spacegroup (typically in Hermann-Mauguin notation).

If this is undesirable for any reason, the generators of a Spacegroup can be
provided after the lattice has been created. One must only specify the generators
of the Spacegroup since they are a subset of the full group which can be
combined to reproduce all group elements. If one prefers, the entire Spacegroup
motions can be provided instead.

Sodium Chloride symmetry
------------------------

The Hall symbol for cubic :math:`F m \bar{3} m` is :math:`-F 4 2 3`.
The Hall symbol *is* the encoded generators of the spacegroup, and can be used
to define a :py:class:`~brille._brille.HallSymbol` which has the functionality
to decode the generators. Those decoded generators can then be used to generate
the full spacegroup:

.. literalinclude:: tutorial_01.py
  :lines: 201-203


This spacegroup is comprised of :math:`192` motions,

.. raw:: html

   <details>
   <summary><a>(show the 192 lines of motions)</a></summary>

.. literalinclude:: tutorial_01.py
  :lines: 4-195

.. raw:: html

   </details>


In this compact form each :math:`4 \times 3` array is the :math:`3 \times 3`
rotation-like part of the motion while the remaining three values are twice the
translation part of the motion.
They can be separated and used to define a :py:class:`~brille._brille.Symmetry`
object plus a set of generators for the group:

.. literalinclude:: tutorial_01.py
  :lines: 196-199

The symmetry information for NaCl is available in :download:`a script <tutorial_01.py>`
which also contains different methods of constructing a real space lattice
which contains the spacegroup symmetry information.

As a teaser, the two lattices `lat0` and `lat1` are identically the same:

.. code-block:: python

  a = 5.69  # angstrom, the approximate lattice constant for NaCl

  lat0 = brille.Direct((a, a, a), (90, 90, 90), '-F 4 2 3')

  lat1 = brille.Direct((a, a, a), (90, 90, 90))
  lat1.spacegroup = brille.Symmetry('-F 4 2 3').generators
