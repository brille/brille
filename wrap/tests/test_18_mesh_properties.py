"""Properties any mesh of the irreducible zone must have (see mesh_properties.py).

BZMeshQdc's structured mesh fills the zone, is conforming, and its boundary
matches itself under the zone's face pairings. The TetGen mesh it replaced
failed the last test except where the zone is bounded by mirror planes only
(GitHub #114).
"""
import warnings

import pytest
from brille import BrillouinZone, BZMeshQdc, Lattice
from mesh_properties import check

CASES = {
    "triclinic -P 1": (((4.02, 4.90, 3.29), (98.97, 86.24, 88.47)), "-P 1", False),
    "monoclinic C": (((5.1, 8.8, 5.2), (90, 104.5, 90)), "-C 2y", False),
    "tetragonal I": (((4.0, 4.0, 9.0), (90, 90, 90)), "-I 4 2", False),
    "trigonal P 3, time reversal": (((3.1, 3.1, 4.7), (90, 90, 120)), "P 3", True),
    "hexagonal 6/m": (((3.0, 3.0, 5.0), (90, 90, 120)), "-P 6", False),
    "hexagonal 6/mmm": (((3.0, 3.0, 5.0), (90, 90, 120)), "-P 6 2", False),
    "rhombohedral -3m": (((4.9, 4.9, 13.8), (90, 90, 120)), '-R 3 2"', False),
    "fcc": (((5.64,) * 3, (90,) * 3), "-F 4 2 3", False),
    "cubic P": (((4.0,) * 3, (90,) * 3), "-P 4 2 3", False),
}


def properties(name):
    lengths_angles, symmetry, time_reversal = CASES[name]
    lattice = Lattice(lengths_angles, spacegroup=symmetry)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        bz = BrillouinZone(lattice, time_reversal_symmetry=time_reversal)
        mesh = BZMeshQdc(bz, max_size=bz.ir_polyhedron.volume / 50)
    return check(bz, mesh, lattice, time_reversal=time_reversal)


@pytest.mark.parametrize("name", CASES)
def test_mesh_fills_the_irreducible_zone(name):
    assert properties(name)["volume"] == pytest.approx(1, abs=1e-9)


@pytest.mark.parametrize("name", CASES)
def test_mesh_is_conforming(name):
    assert properties(name)["open_faces"] == 0


@pytest.mark.parametrize("name", CASES)
def test_mesh_boundary_matches_itself(name):
    result = properties(name)
    assert result["boundary_missing"] == 0, f"{result['boundary_missing']} of {result['boundary_checked']} boundary images missing"
