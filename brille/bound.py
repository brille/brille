# The compiled module is installed at the base of the package directory structure such that its properties
# can be included here. If the module is loaded, however, without being installed first these import statements will
# raise a ModuleNotFoundError, which prevents automated pre-installation testing from taking place.
from importlib.util import find_spec

bound_spec = find_spec("brille._brille")  # won't this cause a circular import?
if bound_spec is not None:
    try:
        from ._brille import (
            __version__,
            version,
            AngleUnit,
            LengthUnit,
            Bravais,
            RotatesLike,
            NodeType,
            real_space_tolerance,
            reciprocal_space_tolerance,
            ApproxConfig,
            Basis,
            HallSymbol,
            Spacegroup,
            Pointgroup,
            Symmetry,
            PointSymmetry,
            PrimitiveTransform,
            BrillouinZone,
            Polyhedron,
            LPolyhedron,
            SortingStatus,
            BZTrellisQcc,
            BZTrellisQdc,
            BZTrellisQdd,
            BZMeshQcc,
            BZMeshQdc,
            BZMeshQdd,
            BZNestQcc,
            BZNestQdc,
            BZNestQdd,
        )

        # Store the grid types for use in, e.g., the plotting routines
        __grid_types__ = (
            BZTrellisQcc,
            BZTrellisQdc,
            BZTrellisQdd,
            BZMeshQcc,
            BZMeshQdc,
            BZMeshQdd,
            BZNestQcc,
            BZNestQdc,
            BZNestQdd,
        )

        __all__ = [
            "__version__",
            "version",
            "AngleUnit",
            "LengthUnit",
            "Bravais",
            "RotatesLike",
            "NodeType",
            "real_space_tolerance",
            "reciprocal_space_tolerance",
            "ApproxConfig",
            "Basis",
            "HallSymbol",
            "Spacegroup",
            "Pointgroup",
            "Symmetry",
            "PointSymmetry",
            "PrimitiveTransform",
            "BrillouinZone",
            "Polyhedron",
            "LPolyhedron",
            "SortingStatus",
            "BZTrellisQcc",
            "BZTrellisQdc",
            "BZTrellisQdd",
            "BZMeshQcc",
            "BZMeshQdc",
            "BZMeshQdd",
            "BZNestQcc",
            "BZNestQdc",
            "BZNestQdd",
            "__grid_types__",
        ]
    except ImportError as imp_error:
        import platform

        if platform.system == "Windows" and "DLL" in imp_error.msg:
            msg = "You may be missing the latest Visual C++ redistributable package"
            msg += " install it from Microsoft @ https://support.microsoft.com/en-us/help/2977003"
            msg += " before trying to import brille again"
            print(msg)
        raise imp_error
else:
    __version__ = "0.0.0"
    __all__ = ["__version__"]
