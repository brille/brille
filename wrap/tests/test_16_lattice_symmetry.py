"""The lattice's own symmetry: operations that do not fit it, and near-symmetry.

Symmetry operations that are not symmetries of the lattice are refused when the
Lattice is made, instead of failing later with no irreducible zone.

A lattice close to a more symmetric one has a Brillouin zone with features far
smaller than itself: this R lattice is within 6 ppm of face-centred cubic, and
one zone edge is 1e-6 of the zone's size. brille keeps the exact zone, warns,
and limits mesh refinement, which otherwise never ended.
"""
import warnings

import pytest
from brille import BrillouinZone, BZMeshQdc, Lattice, NearSymmetryWarning


def near_fcc():
    """An R lattice (AFLOW, Hall 433) within 6 ppm of face-centred cubic"""
    return Lattice(([6.885854708699999, 6.885854708699999, 8.4334670046], [90, 90, 120]), spacegroup="R 3")


def test_operations_must_map_the_centred_lattice_onto_itself():
    # 'R 3 -2' lost the quote of 'R 3 -2"' (R3m); its mirrors swap the obverse and reverse settings
    with pytest.raises(ValueError, match="rhombohedrally centred lattice onto itself"):
        Lattice(((9.619, 9.619, 3.1499), (90, 90, 120)), spacegroup="R 3 -2")
    Lattice(((9.619, 9.619, 3.1499), (90, 90, 120)), spacegroup='R 3 -2"')


def test_operations_must_preserve_the_metric():
    with pytest.raises(ValueError, match="changes the lattice metric"):
        Lattice(((3.5, 3.5, 12.9), (90, 90, 120)), spacegroup="P 4")


def test_near_symmetry_warns():
    with pytest.warns(NearSymmetryWarning, match="within 1e-4 of m-3m"):
        bz = BrillouinZone(near_fcc())
    assert bz.lattice_symmetry_counts() == (12, 48)


def test_near_symmetry_warning_can_be_silenced():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        BrillouinZone(near_fcc(), warn_near_symmetry=False)


@pytest.mark.parametrize("lattice, symmetry, count", [
    (((4.0, 4.0, 4.0), (90, 90, 90)), "P 1", 48),
    (((3.0, 3.0, 5.0), (90, 90, 120)), "-P 6 2", 24),
    (((5.1, 8.8, 5.2), (90, 104.5, 90)), "-C 2y", 4),
])
def test_no_warning_away_from_higher_symmetry(lattice, symmetry, count):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        bz = BrillouinZone(Lattice(lattice, spacegroup=symmetry))
    assert bz.lattice_symmetry_counts() == (count, count)


def test_refinement_stops_near_higher_symmetry():
    """Quality refinement of this zone never ended; now it stops at a default limit"""
    bz = BrillouinZone(near_fcc(), warn_near_symmetry=False)
    with pytest.warns(RuntimeWarning, match="refinement stopped"):
        mesh = BZMeshQdc(bz, num_levels=1)
    assert mesh.refinement_limited


def test_refinement_limit_is_reported():
    bz = BrillouinZone(near_fcc(), warn_near_symmetry=False)
    with pytest.warns(RuntimeWarning, match="refinement stopped"):
        mesh = BZMeshQdc(bz, max_size=bz.ir_polyhedron.volume / 100, max_points=300)
    assert mesh.refinement_limited and mesh.rlu.shape[0] < 400


def test_ordinary_mesh_is_not_limited():
    bz = BrillouinZone(Lattice(((5.64,) * 3, (90,) * 3), spacegroup="-F 4 2 3"))
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        mesh = BZMeshQdc(bz, max_size=bz.ir_polyhedron.volume / 1000)
    assert not mesh.refinement_limited
