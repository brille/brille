#!/usr/bin/env python3
"""Run tests of Brillouin zone creation."""
from importlib import util
import os
import sys
import unittest
import json
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

def find_tests_dir():
    """Locate the tests directory in the repository folder

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
    cifsdir = os.path.join(brilleroot,'..','tests')
    return cifsdir

def get_aflow_lattices():
    jsonfile = os.path.join(find_tests_dir(),'aflow_lattices.json')
    with open(jsonfile,'r') as f:
        data = json.load(f)
    return data

def n2chr(n):
    return chr(48+n) if n<10 else chr(65-10+n)

def test_aflow_crystaldatabase():
    tested = 0
    failed = 0
    errored = 0
    failed_afl = []
    failed_ratio = []
    errored_afl = []
    errored_arg = []
    hall_groups_passed = np.zeros(530,dtype='int')
    hall_groups_failed = np.zeros(530,dtype='int')
    for afl in get_aflow_lattices():
        dlat = s.Direct(afl[1], afl[2], afl[0]) # use the pre-determined Hall number
        i = dlat.hall
        try:
            bz = s.BrillouinZone(dlat.star)
            vol_bz = bz.polyhedron.volume
            vol_ir = bz.ir_polyhedron.volume
            tested += 1
            if not np.isclose(vol_ir, vol_bz/s.PointSymmetry(i).size):
                failed += 1
                failed_afl.append(afl)
                failed_ratio.append(vol_ir/vol_bz*s.PointSymmetry(i).size)
                hall_groups_failed[i-1] += 1
            else:
                hall_groups_passed[i-1] += 1
        except Exception as err:
            errored += 1
            errored_afl.append(afl)
            errored_arg.append(err.args)
    if failed > 0:
        print("\nFailed to find correct irreducible Brillouin zone for",failed,"out of",tested,"lattices")
        for file,  rat in zip(failed_afl, failed_ratio):
            print(file,rat)
    if errored > 0:
        print("\nException raised for",errored,"out of",tested,"lattices")
        for file,  arg in zip(errored_afl,  errored_arg):
            print(file,arg)
    print("\nHall groups passed (total =",hall_groups_passed.sum(),"of",tested,"tested)")
    encoded_hgp = [n2chr(x) for x in hall_groups_passed]
    for x in [encoded_hgp[i*53:(i+1)*53] for i in range(10)]:
        print(''.join(x))

if __name__ == '__main__':
    test_aflow_crystaldatabase()
