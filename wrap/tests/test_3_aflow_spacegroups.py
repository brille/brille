#!/usr/bin/env python3
"""Run tests of Brillouin zone creation."""
import os
import sys
import unittest
import numpy as np
from pathlib import Path
from importlib.util import find_spec
addpaths = [Path(), Path('..')] # we might be in .../build/wrap/
config = os.environ.get('CMAKE_CONFIG_TYPE') # set by ctest -C <cfg>
if config:
  for path in addpaths:
    if Path(path, config).exists():
      addpaths.append(Path(path,config))
sys.path[:0] = [str(path.absolute()) for path in addpaths]

if find_spec('_brille') is not None:
    import _brille as s
elif find_spec('brille') is not None and find_spec('brille._brille') is not None:
    import brille as s
else:
    abspaths = [str(path.absolute()) for path in addpaths]
    raise Exception("brille module not found in {}!".format(abspaths))

def get_local_JSON(filename):
    import os
    import pathlib
    import json
    dir = os.path.dirname(os.path.abspath(__file__))
    with open(str(pathlib.Path(dir, filename)), 'r') as f:
        data = json.load(f)
    return data

def get_aflow_lattices():
    return get_local_JSON('aflow_lattices.json')

def n2chr(n):
    return chr(48+n) if n<10 else chr(65-10+n)

class AflowTest (unittest.TestCase):
    def test_aflow_crystaldatabase(self):
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
            # afl == [hall_number, basis_vector_lengths, basis_vector_angles, Hall_symbol]
            dlat = s.Direct(afl[1], afl[2], afl[3])
            tested += 1
            try:
                bz = s.BrillouinZone(dlat.star)
                vol_bz = bz.polyhedron.volume
                vol_ir = bz.ir_polyhedron.volume
                if not np.isclose(vol_ir, vol_bz/s.PointSymmetry(afl[0]).size):
                    failed += 1
                    failed_afl.append(afl)
                    failed_ratio.append(vol_ir/vol_bz*s.PointSymmetry(afl[0]).size)
                    hall_groups_failed[afl[0]-1] += 1
                else:
                    hall_groups_passed[afl[0]-1] += 1
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

        self.assertTrue(failed == 0)
        self.assertTrue(errored == 0)

if __name__ == '__main__':
    unittest.main()
