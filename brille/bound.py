# The compiled module is installed at the base of the package directory structure such that its properties
# can be included here. If the module is loaded, however, without being installed first these import statements will
# raise a ModuleNotFoundError, which prevents automated pre-installation testing from taking place.
try:
    from ._brille import __version__, version, AngleUnit, LengthUnit, Bravais, RotatesLike, NodeType
    from ._brille import real_space_tolerance, reciprocal_space_tolerance, ApproxConfig
    from ._brille import Basis, HallSymbol, Spacegroup, Pointgroup, Symmetry, PointSymmetry, PrimitiveTransform
    from ._brille import BrillouinZone, Polyhedron, LPolyhedron, SortingStatus
    from ._brille import BZTrellisQcc, BZTrellisQdc, BZTrellisQdd
    from ._brille import BZMeshQcc, BZMeshQdc, BZMeshQdd
    from ._brille import BZNestQcc, BZNestQdc, BZNestQdd

    # Store the grid types for use in, e.g., the plotting routines
    __grid_types__ = (
        BZTrellisQcc, BZTrellisQdc, BZTrellisQdd,
        BZMeshQcc, BZMeshQdc, BZMeshQdd,
        BZNestQcc, BZNestQdc, BZNestQdd,
    )
except ModuleNotFoundError:
    pass
except ImportError as imp_error:
    import platform
    if platform.system == 'Windows' and 'DLL' in imp_error.msg:
        msg = 'You may be missing the latest Visual C++ redistributable package'
        msg += ' install it from Microsoft @ https://support.microsoft.com/en-us/help/2977003'
        msg += ' before trying to import brille again'
        print(msg)
    raise imp_error
