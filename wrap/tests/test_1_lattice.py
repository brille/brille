#!/usr/bin/env python3
import os
import sys
import unittest
import numpy as np
from pathlib import Path
from importlib.util import find_spec
addpath = Path() # the relative working directory path
config = os.environ.get('CMAKE_CONFIG_TYPE') # set by ctest -C <cfg>
if config:
    # print('CMAKE_BUILD_TYPE = {}'.format(config))
    if Path(addpath, config).exists():
        addpath = Path(addpath, config)
    elif Path(addpath, '..', config).exists():
        addpath = Path(addpath, '..', config)
# print('adding {} to Python search path'.format(addpath))
sys.path.append(str(addpath.absolute()))

if find_spec('_brille') is not None:
    import _brille as s
    # print(find_spec('_brille'))
elif find_spec('brille') is not None and find_spec('brille._brille') is not None:
    import brille as s
    # print(find_spec('brille._brille'))
else:
    raise Exception("brille module not found!")

class Lattice (unittest.TestCase):
    def test_a_init2(self):
        a=1.0
        b=np.pi/2
        v = 1.0
        l = s.Lattice(a,a,a,b,b,b) # compute the volume on instantiation
        self.assertEqual(l.a,a)
        self.assertEqual(l.b,a)
        self.assertEqual(l.c,a)
        self.assertEqual(l.alpha,b)
        self.assertEqual(l.beta,b)
        self.assertEqual(l.gamma,b)
        self.assertAlmostEqual(l.volume,v) # since 0.9999999999999999 != 1.0
    def test_a_init3(self):
        a = (1,1,1)
        b = np.array([1,1,1])*np.pi/2
        v = 1.0
        l = s.Lattice(a,b)
        self.assertEqual(l.a,a[0])
        self.assertEqual(l.b,a[1])
        self.assertEqual(l.c,a[2])
        self.assertEqual(l.alpha,b[0])
        self.assertEqual(l.beta ,b[1])
        self.assertEqual(l.gamma,b[2])
        self.assertAlmostEqual(l.volume,v) # since 0.9999999999999999 != 1.0
    def test_b_tensors(self):
        a = (6,6,9)
        b = np.array([90,90,120])*np.pi/180
        l = s.Lattice(a,b)
        cmt = np.array([[a[0]**2, a[0]*a[1]*np.cos(b[2]), 0],
                        [a[0]*a[1]*np.cos(b[2]), a[1]**2, 0],
                        [0, 0, a[2]**2]])
        self.assertAlmostEqual( (l.get_covariant_metric_tensor()-cmt).sum(), 0)
        cnt = np.linalg.inv(cmt)
        self.assertAlmostEqual( (l.get_contravariant_metric_tensor()-cnt).sum(), 0)
    def test_c_star(self):
        a = (6,6,9)
        b = np.array([90,90,120])*np.pi/180
        l = s.Lattice(a,b)
        with self.assertRaises(RuntimeError):
            l.star
        with self.assertRaises(RuntimeError):
            l.isstar(l)
    def test_d_equality(self):
        a=1.0
        b=np.pi/2
        v = 1.0
        l1 = s.Lattice(a,a,a,b,b,b)
        l2 = s.Lattice(a,a,a,b,b,b)
        self.assertTrue(l1.issame(l2))
        self.assertEqual(l1,l2)
    def test_e_subclass_star(self):
        d = s.Direct((1,1,1),np.array([1,1,1])*np.pi/2)
        r = s.Reciprocal(np.array([1,1,1])*np.pi*2, np.array([1,1,1])*np.pi/2)
        self.assertTrue( d.isstar(r) )
        self.assertTrue( r.isstar(d) )
        self.assertEqual( d, r.star )
        self.assertEqual( r, d.star )

if __name__ == '__main__':
  unittest.main()
