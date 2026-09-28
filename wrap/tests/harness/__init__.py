"""Reference-model test harness for brille.

The harness compares brille's symmetry-unfolded interpolation results against
a model that can be evaluated exactly at any q. The model is an interface
(:py:class:`ReferenceModel`) so that phonon models can later be joined by
spin-wave (LSWT) and SU(N) models without changing the comparison machinery.

Modules
-------
crystals
    Crystals on real (AFLOW) lattices with a prescribed space group, and their primitive reductions.
models
    The :py:class:`ReferenceModel` interface and a Born–von Kármán phonon model.
compare
    Degeneracy-aware comparison of modes, parameterised by an inner-product metric.
adapter
    Construction of brille objects from a model, and the calls into brille.
queries
    Choosing q points that brille can be asked about safely.
"""
from .crystals import Crystal, aflow_crystal, aflow_lattices, hall_numbers_by_ita, primitive
from .models import Modes, ReferenceModel, BornVonKarman
from .compare import ModeComparison, compare_modes, degenerate_groups
from .adapter import BrilleSetup, brille_lattice, build, interpolate
from .queries import interior_vertex_indices, point_group_images

__all__ = [
    "Crystal",
    "aflow_crystal",
    "aflow_lattices",
    "primitive",
    "hall_numbers_by_ita",
    "Modes",
    "ReferenceModel",
    "BornVonKarman",
    "ModeComparison",
    "compare_modes",
    "degenerate_groups",
    "BrilleSetup",
    "brille_lattice",
    "build",
    "interpolate",
    "interior_vertex_indices",
    "point_group_images",
]
