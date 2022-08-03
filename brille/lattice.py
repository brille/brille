def _t(args):
    return args if isinstance(args, tuple) else tuple(args)


def _make(to_type, args):
    return args if isinstance(args, to_type) else to_type(*_t(args))


# def Lattice(values, /, *, spacegroup=None, symmetry=None, basis=None, **kwargs):
def Lattice(values, spacegroup=None, symmetry=None, basis=None, **kwargs):
    """  Construct a space-spanning lattice in three dimensions

    A space-spanning lattice in :math:`N` dimensions has :math:`N` basis vectors
    which can be described fully by their :math:`N` lengths and the
    :math:`\\frac{1}{2} N (N-1)` angles between each set of basis vectors, or
    :math:`\\frac{1}{2}N(N+1)` scalars in total.
    This class stores the basis vectors of the lattice described in an orthonormal space,
    plus the metric of the space, and the equivalent information for the dual of the lattice.

    Examples
    --------

    >>> from brille import Lattice
    >>> a, c = 3.95, 12.9
    >>> avec, bvec, cvec = [a, 0, 0], [0, a, 0], [0, 0, c]
    >>> basis = ([[0, 0, 0], [0.5, 0.5, 0.5], [0.25, 0.25, 0.3]], [0, 0, 1]))
    >>> symmetry = (rotations, translations)
    >>> lat = Lattice(([a, a, c], [90, 90, 90]), spacegroup='I4/mmm', basis=basis)
    >>> lat = Lattice(([avec, bvec, cvec],), symmetry=symmetry, basis=basis)
    >>> lat = Lattice([avec, bvec, cvec], spacegroup='I4/mmm', basis=basis)

    Note
    ----
    The parameters packed into ``values`` are position-based and should be *either* ``(lengths, angles)``
    *or* ``(vectors,)`` [which also requires setting ``row_vectors=False`` if the matrix represents column vectors].
    The lengths or vector components are interpreted in angstrom or inverse angstrom for real or reciprocal lattice
    parameters, respectively, and are assumed to be angstrom if the keyword ``real_space`` is missing or True.


    Parameters
    ----------

    values: (lengths, angles), vectors
        See the note above and the documenation for ``lengths``, ``angles``, and ``vectors`` below.
    lengths : list, tuple, numpy.ndarray
        The three basis vector lengths in angstrom for real lattice, or inverse angstrom for reciprocal lattices
    angles : list, tuple, numpy.ndarray
        The three angles between the basis vectors in degrees or radians. The angles are interpreted as radians
        if none are greater than pi and are otherwise assumed to be degrees.
    vectors: list[list,...], tuple(tuple,...), numpy.ndarray
        The basis vectors in angstrom for real lattices, or inverse angstrom for reciprocal lattices,
        expressed in an orthorhombic coordinate system. An optional keyword argument, ``row_vectors``, identifies
        if the provided basis vectors are row vectors [`row_vectors=True`, default] or column vectors
        [``row_vectors=False``].
    spacegroup : str, tuple(str, str)
        The International Tables name, Hermann-Mauguin symbol with optional choice, or Hall symbol for the
        spacegroup of the lattice. The spacegroup may be provided as positional argument(s) or by keyword.
        If present, the ``symmetry`` keyword must not be used.
        Valid syntax for ``(Hermann-Mauguin, choice)`` input depends on the spacegroup but is generally one of:

        * a single letter ('a', 'b', or 'c') with possible prepended '-', denoting unique-axis choice
        * a single digit ('1', '2', or '3'), denoting origin choice
        * a letter and digit, denoting unique-axis and origin choice
        * a permutation of 'abc' with possible '-' before one of the letters, denoting axis permutation
        * or 'R' or 'H' for trigonal systems with Rhombohedral or Hexagonal lattice settings, respectively.
        
        Acceptable values are contained in the C++ source code in the seventh column of
        `this table <https://github.com/brille/brille/blob/eecb4cb28227665908793abc47e88c69518c09fc/src/spg_database.cpp#L63-L610>`_
        with each line representing one spacegroup with values, in order, defined by the
        `class signature <https://github.com/brille/brille/blob/eecb4cb28227665908793abc47e88c69518c09fc/src/spg_database.hpp#L80-L91>`_
        These values come from spglib, which likely obtained them from
        `Seto's Home Page <https://web.archive.org/web/20210621195003/http://pmsl.planet.sci.kobe-u.ac.jp/~seto/?page_id=37&lang=en>`_.
    symmetry: brille.Symmetry, tuple(matrices, vectors), str
        The spacegroup symmetry operations as an object, a tuple of the (pseudo)rotation matrices and translation
        vectors, or a CIF xyz encoded string. The symmetry information must be provided by keyword.
        If present, the 'spacegroup' positional argument(s) or keyword must not be used.
    basis: brille.Basis, tuple(positions, types)
        The atom basis information of the lattice, expressed in units of the real space basis vectors. The types
        must be integer and are used only to identify equivalent atoms -- they should probably be contiguous from
        zero to 1-N where N is the number of unique atoms in the atom basis.
        If present, either spacegroup or symmetry information must be provided.
    kwargs:
        Keyword arguments are passed to the :py:class:`brille._brille.Lattice` constructor,
        see its documentation for details.

    """
    from ._brille import Lattice as _Lattice, Symmetry as _Symmetry, Basis as _Basis

    if isinstance(spacegroup, str):
        spacegroup = spacegroup,
    if spacegroup is None:
        spacegroup = tuple()
    if len(spacegroup) and not all([isinstance(x, str) for x in spacegroup]):
        raise ValueError("Keyword argument spacegroup should be one of: str, (str,), or (str, str)")

    if basis is not None and symmetry is None and len(spacegroup) == 0:
        raise ValueError("Providing basis without symmetry or spacegroup is not allowed")
    if symmetry is not None and len(spacegroup):
        raise ValueError("Providing both spacegroup and symmetry is not allowed")

    # The pybind11 wrapper has overloads which accept:
    #   3-tuple-like, 3-tuple-like : 3-element arrays of lengths and angles
    #   3-tuple(3-tuple-like)-like : 3 3-element arrays of basis vectors
    #
    # with the basis vectors ((ax, ay, az), (bx, by, bz), (cx, cy, cz)) or ((ax, bx, cx), (ay, by, cy), (az, bz, cz))
    #
    # If the user provided a 2-tuple then it must be (lengths, angles), otherwise it must be vectors
    # We use less-than to allow a user to provide (vectors,) as the first input
    args = values if len(values) < 3 else (values, )

    if symmetry is None and basis is None:
        return _Lattice(*args, *spacegroup, **kwargs)
    elif symmetry is None:
        return _Lattice(*args, *spacegroup, _make(_Basis, basis), **kwargs)

    if basis is None:
        return _Lattice(*args, _make(_Symmetry, symmetry), **kwargs)
    else:
        return _Lattice(*args, _make(_Symmetry, symmetry), _make(_Basis, basis), **kwargs)
