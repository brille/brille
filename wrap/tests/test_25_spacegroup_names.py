"""A space group given by name is the one named, however it is spaced, or refused.

Names are matched ignoring spaces and subscript marks, which no two names in the
table differ by alone. A string that is not a name, a Hall symbol or CIF xyz
operations is refused: before, the Hall symbol parser skipped the characters it
did not know, so 'P 2/m' became the Hall symbol 'P 2', a different group.
"""
import pytest
from brille import Lattice, Symmetry

MONOCLINIC = ([4.0, 5.0, 6.0], [90, 100, 90])      # unique axis b
ORTHORHOMBIC = ([4.0, 5.0, 6.0], [90, 90, 90])


def group(cell, spacegroup):
    return Lattice(cell, spacegroup=spacegroup).spacegroup


@pytest.mark.parametrize("name", ["P2/m", "P 2/m", ("P 2/m", "b"), ("P2/m", "b"), "P 1 2/m 1", "-P 2y"])
def test_p2m_unique_axis_b(name):
    assert group(MONOCLINIC, name) == group(MONOCLINIC, "P 1 2/m 1")


@pytest.mark.parametrize("name", ["P2_1/c", "P 2_1/c", "P21/c", "P 21/c", "P 1 21/c 1", ("P21/c", "b1")])
def test_p21c_spellings(name):
    reference = group(MONOCLINIC, "-P 2ybc")
    assert group(MONOCLINIC, name) == reference
    assert len(reference.W) == 4


def test_choice_selects_the_setting():
    # unique axis c needs a lattice with its non-right angle between a and b
    c_setting = group(([4.0, 5.0, 6.0], [90, 90, 100]), ("P 2/m", "c"))
    assert c_setting == group(([4.0, 5.0, 6.0], [90, 90, 100]), "P 1 1 2/m")


@pytest.mark.parametrize("name", ["Pmmm", "P m m m", "P 2/m 2/m 2/m"])
def test_orthorhombic_spellings(name):
    assert len(group(ORTHORHOMBIC, name).W) == 8


def test_hall_symbols_with_a_change_of_basis():
    hexagonal = ([4.0, 4.0, 6.0], [90, 90, 120])
    assert len(group(hexagonal, "P 31 2c (0 0 1)").W) == 6


@pytest.mark.parametrize("name", ["P 2/q", "Pnonsense", ("P2/m", "q"), "P 2 garbage"])
def test_unknown_names_are_refused(name):
    with pytest.raises(ValueError, match="not a space group"):
        Lattice(MONOCLINIC, spacegroup=name)


def test_cif_xyz_operations_still_work():
    assert len(Lattice(ORTHORHOMBIC, symmetry=Symmetry("x,y,z;-x,-y,-z")).spacegroup.W) == 2
