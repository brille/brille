def _t(args):
    return args if isinstance(args, tuple) else tuple(args)


def _make(to_type, args):
    return args if isinstance(args, to_type) else to_type(*_t(args))


def wrap_lattice(lattice_type, symmetry_type, basis_type):
    def wrapped(*args, spacegroup=None, symmetry=None, basis=None, **kwargs):
        """  Construct a space-spanning lattice in three dimensions

        A space-spanning lattice in :math:`N` dimensions has :math:`N` basis vectors
        which can be described fully by their :math:`N` lengths and the
        :math:`\\sum_1^{N-1} 1` angles between each set of basis vectors, or
        :math:`\\sum_1^N 1 = \\frac{1}{2}N(N+1)` scalars in total.
        This class stores the basis vectors of the lattice described in an orthonormal space,
        plus the metric of the space, and the equivalent information for the dual of the lattice.

        Examples
        --------

            lat = brille.Lattice([3.95, 3.95, 12.9], [90, 90, 90], 'I4/mmm',
                                 basis=([[0, 0, 0], [0.5, 0.5, 0.5], [0.25, 0.25, 0.3]], [0, 0, 1]))
            lat = brille.Lattice([3.95, 3.95, 12.9], [90, 90, 90], spacegroup='I4/mmm',
                                 basis=([[0, 0, 0], [0.5, 0.5, 0.5], [0.25, 0.25, 0.3]], [0, 0, 1]))
            lat = brille.Lattice([3.95, 3.95, 12.9], [90, 90, 90], symmetry=(rotations, translations),
                                 basis=([[0, 0, 0], [0.5, 0.5, 0.5], [0.25, 0.25, 0.3]], [0, 0, 1]))
            lat = brille.Lattice([[3.95, 0, 0], [0, 3.95, 0], [0, 0, 12.9]], 'I4/mmm',
                                 basis=([[0, 0, 0], [0.5, 0.5, 0.5], [0.25, 0.25, 0.3]], [0, 0, 1]))

        Note
        ----
        The parameters packed into `args` are position-based and should be *either* `(lengths, angles)`
        *or* `(vectors)` [which also requires setting `row_vectors=False` if the matrix represents column vectors].


        Parameters
        ----------
        args:
            lengths : list, tuple, numpy.ndarray
                The three basis vector lengths in angstrom for real lattice, or inverse angstrom for reciprocal lattices
            angles : list, tuple, numpy.ndarray
                The three angles between the basis vectors in degrees or radians. The angles are interpreted as radians
                if none are greater than pi and are otherwise assumed to be degrees.
            vectors: list[list,...], tuple(tuple,...), numpy.ndarray
                The basis vectors in angstrom for real lattices, or inverse angstrom for reciprocal lattices,
                expressed in an orthorhombic coordinate system. An optional keyword argument, `row_vectors`, identifies
                if the provided basis vectors are row vectors [`row_vectors=True`, default] or column vectors
                [`row_vectors=False`].
        spacegroup : str, tuple(str, str)
            The International Tables name, Hermann-Mauguin symbol with optional choice, or Hall symbol for the
            spacegroup of the lattice. The spacegroup may be provided as positional argument(s) or by keyword.
            If present, the `symmetry` keyword must not be used.
        symmetry: brille.Symmetry, tuple(matrices, vectors), str
            The spacegroup symmetry operations as an object, a tuple of the (pseudo)rotation matrices and translation
            vectors, or a CIF xyz encoded string. The symmetry information must be provided by keyword.
            If present, the 'spacegroup' positional argument(s) or keyword must not be used.
        basis: brille.Basis, tuple(positions, types)
            The atom basis information of the lattice, expressed in units of the real space basis vectors. The types
            must be integer and are used only to identify equivalent atoms -- they should probably be contiguous from
            zero to 1-N where N is the number of unique atoms in the atom basis.
            If present, either spacegroup or symmetry information must be provided.
        """
        if spacegroup is None:
            # look for a string (or two!) on the end of the position arguments
            if len(args) > 2 and isinstance(args[-1], str) and isinstance(args[-2], str):
                spacegroup = args[-2], args[-1]
                args = args[:-2]
            elif len(args) > 1 and isinstance(args[-1], str):
                spacegroup = args[-1],
                args = args[:-1]
            else:
                spacegroup = tuple()

        if symmetry is None and basis is None:
            return lattice_type(*args, *spacegroup, **kwargs)
        elif symmetry is None:
            if len(spacegroup) == 0:
                raise ValueError("Providing a Basis without Symmetry or Spacegroup is not allowed")
            return lattice_type(*args, *spacegroup, _make(basis_type, basis), **kwargs)

        if len(spacegroup):
            raise ValueError("Providing both Spacegroup and Symmetry is not allowed")
        if basis is None:
            return lattice_type(*args, _make(symmetry_type, symmetry), **kwargs)
        else:
            return lattice_type(*args, _make(symmetry_type, symmetry), _make(basis_type, basis), **kwargs)

    return wrapped
