#!/usr/bin/env python3
"""Run tests of Brillouin zone creation."""
from importlib import util
import os
import sys
import unittest
import numpy as np
import pymatgen
try:
    from importlib.util import find_spec
except ImportError:
    raise Exception("This test requires importlib.util (and Python3)")

# We need to tell Python where it can find the brille module.
ADDPATH = os.getcwd()
# It's either in the working directory where python was called or in
# a sub-directory (called Debug using Visual Studio under Windows)
if os.path.exists('Debug'):
    ADDPATH += "\\Debug"
sys.path.append(ADDPATH)

if util.find_spec('brille') is not None and util.find_spec('brille._brille') is not None:
    import brille as s
elif util.find_spec('_brille') is not None:
    # pylint: disable=e0401
    import _brille as s
else:
    raise Exception("brille module not found!")

def cif2direct(filename):
    print(filename)
    cif = pymatgen.Structure.from_file(filename)
    sym = pymatgen.symmetry.analyzer.SpacegroupAnalyzer(cif) # uses spglib?
    l = cif.lattice
    d = s.Direct(l.abc, l.angles, sym.get_hall())
    return d

def findcifdir():
    """Locate the cifs directory in the repository tests folder

    Note that this only works if the package is installed via
        python setup.py develop [--user]
    which builds the C++ library and installs it inside of the repository
    structure with simlinks in the usual package location(s).
    """
    test_spec = find_spec('brille')
    if test_spec is None:
        brilleroot = '.' #punt. maybe we're running this from a build directory?
    else:
        brilleroot = test_spec.submodule_search_locations[0]
    cifsdir = os.path.join(brilleroot,'..','tests','cifs')
    return cifsdir

def test_aflow_crystaldatabase():
    tested = 0
    failed = 0
    errored = 0
    failed_file = []
    failed_lat = []
    failed_ratio = []
    errored_file = []
    errored_lat = []
    errored_arg = []
    hall_groups_passed = np.zeros(530,dtype='int')
    hall_groups_failed = np.zeros(530,dtype='int')
    for dirpath, dirs, files in os.walk(findcifdir()):
        for filename in files:
            dlat = cif2direct(os.path.join(dirpath, filename))
            i = dlat.hall
            try:
                bz = s.BrillouinZone(dlat.star)
                vol_bz = bz.polyhedron.volume
                vol_ir = bz.ir_polyhedron.volume
                tested += 1
                if not np.isclose(vol_ir, vol_bz/s.PointSymmetry(i).size):
                    failed += 1
                    failed_file.append(filename)
                    failed_lat.append(dlat)
                    failed_ratio.append(vol_ir/vol_bz*s.PointSymmetry(i).size)
                    hall_groups_failed[i-1] += 1
                else:
                    hall_groups_passed[i-1] += 1
            except Exception as err:
                errored += 1
                errored_file.append(filename)
                errored_lat.append(dlat)
                errored_arg.append(err.args)
    if failed > 0:
        print("\nFailed to find correct irreducible Brillouin zone for",failed,"out of",tested," Hall groups")
        for file, lat, rat in zip(failed_file, failed_lat, failed_ratio):
            print(file,lat,rat)
    if errored > 0:
        print("\nException raised for",errored,"out of",tested,"(max 530) Hall Groups")
        for file, lat, arg in zip(errored_file, errored_lat, errored_arg):
            print(file,lat,arg)
    print("Hall groups passed\n",hall_groups_passed)
    print("Hall groups failed\n",hall_groups_failed)


if __name__ == '__main__':
    test_aflow_crystaldatabase()
