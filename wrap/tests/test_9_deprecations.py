"""Deprecated parts of the Python API still work, and warn."""
import numpy as np
import pytest
from brille import PointSymmetry, Symmetry


def test_point_symmetry_from_hall_number_warns():
    with pytest.warns(DeprecationWarning, match="PointSymmetry\\(Symmetry"):
        old = PointSymmetry(1)
    assert old.size == PointSymmetry(Symmetry(1)).size


def test_point_symmetry_from_hall_number_can_be_an_error():
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        with pytest.raises(DeprecationWarning):
            PointSymmetry(1)


def test_replacement_matches_for_every_hall_number():
    key = lambda ps: sorted(map(tuple, np.asarray(ps.W).reshape(-1, 9)))
    for hall in range(1, 531):
        with pytest.warns(DeprecationWarning):
            old = PointSymmetry(hall)
        assert key(old) == key(PointSymmetry(Symmetry(hall))), hall
