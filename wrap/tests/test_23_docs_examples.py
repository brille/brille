"""The scripts behind the documentation's tutorials and how-to guides run, and
still show what the pages say they do.
"""
import runpy
from pathlib import Path

import numpy as np
import pytest

DOCS = Path(__file__).resolve().parents[2] / "docs"


def run(relative):
    """Run a documentation script, without its __main__ section, and return its globals"""
    path = DOCS / relative
    if not path.exists():
        pytest.skip(f"{path} is not available")
    return runpy.run_path(str(path), run_name="docs_example")


def test_mesh_tutorial():
    tutorial = run("tutorials/tutorial_03.py")
    mesh, uniform, test_points = tutorial["mesh"], tutorial["uniform"], tutorial["test_points"]
    initial_rlu = tutorial["initial"][0]
    # refinement kept the initial vertices, and made the interpolation better
    # than uniform refinement with more vertices
    assert np.array_equal(np.asarray(mesh.rlu)[:len(initial_rlu)], initial_rlu)
    assert len(uniform.rlu) > len(mesh.rlu)
    assert tutorial["rms_error"](mesh, test_points) < tutorial["rms_error"](uniform, test_points)
    assert tutorial["tetrahedron_errors"](mesh).max() < 0.2


def test_trellis_to_mesh(capsys):
    howto = run("howto/trellis_to_mesh.py")
    howto["compare"]((500,))
    rows = [line.split() for line in capsys.readouterr().out.splitlines()[1:]]
    assert [row[0] for row in rows] == ["trellis", "mesh"]
    # matching sizes give similar vertex counts and errors
    trellis, mesh = (int(row[1]) for row in rows)
    assert 0.5 < mesh / trellis < 2
    trellis, mesh = (float(row[4]) for row in rows)
    assert mesh < 2 * trellis
