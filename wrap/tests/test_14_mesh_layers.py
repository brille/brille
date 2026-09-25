"""Building a mesh links its coarse and fine tetrahedral layers.

The link test once cut one tetrahedron by the other's faces and compared the
computed points, which failed on round-off ("Duplicate intersection point") for
about 2 % of real lattices, mostly centred or near-pseudo-symmetric ones. It now
decides overlap with exact predicates. Lattice values must be exact: one bit
different can avoid the error.
"""
import json
from pathlib import Path

import pytest
from brille import BrillouinZone, BZMeshQdc, Lattice

AFLOW = json.loads((Path(__file__).parent / "aflow_lattices.json").read_text())
# entries of aflow_lattices.json that failed; the last is rhombohedral graphite
FAILED = (9, 14, 152, 166, 179, 195, 303, 313, 344, 359, 371, 381)


@pytest.mark.parametrize("index", FAILED)
def test_aflow_mesh_builds(index):
    _, lengths, angles, hall = AFLOW[index]
    bz = BrillouinZone(Lattice((lengths, angles), spacegroup=hall))
    mesh = BZMeshQdc(bz, max_size=bz.ir_polyhedron.volume / 100)
    assert mesh.rlu.shape[0] > 4


def test_near_monoclinic_mesh_builds():
    """The lattice from the original report, beta = 90.04 degrees"""
    lat = Lattice(((3.580975428, 3.5819754211999997, 3.5869753869), (90.0, 90.04, 90.0)), spacegroup="P 1 m 1")
    bz = BrillouinZone(lat)
    mesh = BZMeshQdc(bz, max_size=bz.ir_polyhedron.volume / 100)
    assert mesh.rlu.shape[0] > 4
