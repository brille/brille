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

HASMPL = find_spec('matplotlib') is not None
HASMPL &= find_spec('mpl_toolkits') is not None
# protect against trying to load a submodule of a non-existant module
if HASMPL:
    HASMPL &= find_spec('mpl_toolkits.mplot3d') is not None
if HASMPL:
    import matplotlib.pyplot as pp
    from mpl_toolkits.mplot3d.axes3d import Axes3D

#pylint: disable=c0103
def plot_points(x, title=''):
    """Plot the 3D points with a given title."""
    if HASMPL:
        fig = pp.figure()
        ax = Axes3D(fig)
        ax.scatter(x[:, 0], x[:, 1], x[:, 2], s=10)
        ax.set_title(title)
        pp.show()

def plot_points_with_lines(x, y, title=''):
    """Plot 3D points and lines with a given title."""
    if HASMPL:
        fig = pp.figure()
        ax = Axes3D(fig)
        ax.plot(y[:, 0], y[:, 1], y[:, 2])
        ax.scatter(x[:, 0], x[:, 1], x[:, 2], s=10)
        ax.set_title(title)
        pp.show()

def plot_2d_points(x, title=''):
    """Plot 2D points and lines with a given title."""
    if HASMPL:
        ax = pp.axes()
        ax.scatter(x[:, 0], x[:, 1], s=10)
        ax.set_title(title)
        pp.show()

def make_dr(a, b, c, al=np.pi/2, be=np.pi/2, ga=np.pi/2, hall=1):
    """Make a Direct and Reciprocal lattice from Direct lattice parameters."""
    d = s.Direct(a, b, c, al, be, ga, hall)
    r = d.star
    return (d, r)

def norm(x):
    return np.sqrt(np.dot(x, x))

def vector_lists_match(A, B):
    pA, pB = np.where(np.isclose(np.abs(A[:,np.newaxis,:] - B).sum(axis=2),0.))
    if pA.size != A.shape[0] or pB.size != B.shape[0]:
        return False
    return np.allclose(A[pA], B[pB])

class BrillouinZone (unittest.TestCase):
    def test_a_init_unit_cube(self):
        # instantiate a cubic lattice with unit-length vectors
        # and its reciprocal lattice, still cubic with 2π-length vectors
        d, r = make_dr(1, 1, 1)
        # creating a BrillouinZone objects requires that we pass the reciprocal
        # lattice object, and passing a Direct object is a TypeError:
        with self.assertRaises(TypeError):
            s.BrillouinZone(d)
        # its first Brillouin zone is a cube in reciprocal space
        bz = s.BrillouinZone(r)
        # with six faces, defined by: (00̄1),(0̄10),(̄100),(100),(010),(001)
        p = bz.points
        self.assertEqual(p.ndim, 2)
        self.assertEqual(p.shape[0], 6) # and there is one for each of the six faces
        self.assertEqual(p.shape[1], 3) # the face vectors are 3-vectors
        n = bz.normals
        self.assertEqual(n.ndim, 2)
        self.assertEqual(n.shape[0], 6) # and there is one for each of the six faces
        self.assertEqual(n.shape[1], 3) # the face normals are 3-vectors
        n_dot_p = np.array([np.dot(x/norm(x),y/norm(y)) for x,y in zip(n,p)])
        self.assertTrue((n_dot_p == 1.0).all())

        expected = np.array([[-1,0,0],[0,-1,0],[0,0,-1],[0,0,1],[0,1,0],[1,0,0]])
        self.assertTrue(vector_lists_match(p, expected/2))
        self.assertTrue(vector_lists_match(np.array([x/norm(x) for x in n]), expected))
        #
        # the vertices of the first Brillouin zone are the 8 corners of the 1/2-unit cube:
        expected = np.array([[-1,-1,-1], [-1,-1, 1], [-1, 1,-1], [-1, 1, 1], [1,-1,-1], [1,-1, 1], [1, 1,-1], [1, 1, 1]])/2
        verts = bz.vertices
        self.assertEqual(verts.ndim, 2)
        self.assertEqual(verts.shape[0], 8)
        self.assertEqual(verts.shape[1], 3)
        self.assertTrue(vector_lists_match(verts, expected))

    def test_b_isinside_unit_cube(self):
        d, r = make_dr(1, 1, 1)
        bz = s.BrillouinZone(r)
        # the first Brillouin zone of this unit-cube is bounded by
        # {(00̄1),(0̄10),(̄100),(100),(010),(001)}/2
        face_centres = np.array([ [-1, 0, 0], [0,-1, 0], [0, 0,-1], [0, 0, 1], [0, 1, 0], [1, 0, 0] ], dtype='double')/2
        self.assertTrue( np.all(bz.isinside(face_centres)) )
        # including the vertices (since they are the corners of the zone)
        corners = np.array([[-1,-1,-1], [-1,-1, 1], [-1, 1,-1], [-1, 1, 1], [1,-1,-1], [1,-1, 1], [1, 1,-1], [1, 1, 1]], dtype='double')/2
        self.assertTrue( np.all(bz.isinside(corners)) )
        # so all points with h, k, and l in the range [-0.5, 0.5] are in the zone
        Q = np.random.rand(100, 3) - 0.5 # this is a uniform distribution over [-0.5, 0.5) -- close enough
        self.assertTrue( np.all(bz.isinside(Q)) )
        self.assertFalse( np.all(bz.isinside(Q+5)) )

    def test_c_moveinto_unit_cube(self):
        d, r = make_dr(1, 1, 1)
        bz = s.BrillouinZone(r)
        Q = (np.random.rand(100, 3)-0.5) * 10 # this is a uniform distribution over [-5, 5)
        if np.all(bz.isinside(Q)): # this is vanishingly-unlikely
            Q += 100.0
        self.assertFalse( np.all(bz.isinside(Q)) )
        (q, tau) = bz.moveinto(Q)
        self.assertTrue( np.all(bz.isinside(q)) )
        self.assertAlmostEqual( np.abs(Q-q-tau).sum(), 0)
#
#
#
    def test_a_init_hexagonal(self):
        # instantiate a hexagonal lattice and its reciprocal lattice, still hexagonal
        d, r = make_dr(3, 3, 9, np.pi/2, np.pi/2, np.pi*2/3)
        # creating a BrillouinZone objects requires that we pass the reciprocal
        # lattice object, and passing a Direct object is a TypeError:
        with self.assertRaises(TypeError):
            s.BrillouinZone(d)
        # its first Brillouin zone is a cube in reciprocal space
        bz = s.BrillouinZone(r)
        # with eight faces, defined by: (̄100),(̄110),(0̄10),(001),(001),(010),(1̄10),(100)
        p = bz.points
        self.assertEqual(p.ndim, 2)
        self.assertEqual(p.shape[0], 8) # and there is one for each of the eight faces
        self.assertEqual(p.shape[1], 3) # the face vectors are 3-vectors
        n = bz.normals
        self.assertEqual(n.ndim, 2)
        self.assertEqual(n.shape[0], 8) # and there is one for each of the six faces
        self.assertEqual(n.shape[1], 3) # the face normals are 3-vectors
        n_dot_p = np.array([np.dot(x/norm(x),y/norm(y)) for x,y in zip(n,p)])
        self.assertTrue(np.allclose(n_dot_p, 1.))

        expected= np.array([ [-1, 0, 0],[-1, 1, 0],[0,-1, 0],[0, 0,-1],[0, 0, 1],[0, 1, 0],[1,-1, 0],[1, 0, 0] ])
        self.assertTrue(vector_lists_match(p, expected/2))
        n_compare = np.array([x/np.max(np.abs(x)) for x in n])
        self.assertTrue(vector_lists_match(n_compare, expected))
        self.assertTrue( np.all(bz.isinside(expected/2)) )
        #
        # the vertices of the first Brillouin zone are the 12 corners of the hexagonal-prism:
        expected = np.array([[-4, 2,-3],[-2,-2,-3],[-4, 2, 3],[-2,-2, 3],[-2, 4,-3],[-2, 4, 3],
                             [ 2,-4,-3],[ 2,-4, 3],[ 2, 2,-3],[ 4,-2,-3],[ 2, 2, 3],[ 4,-2, 3]])/6
        verts = bz.vertices
        self.assertEqual(verts.ndim, 2)
        self.assertEqual(verts.shape[0], 12)
        self.assertEqual(verts.shape[1], 3)
        self.assertTrue(vector_lists_match(verts, expected))
        self.assertTrue( np.all(bz.isinside(expected)) )

    def test_b_isinside_hexagonal(self):
        d, r = make_dr(3, 3, 9, np.pi/2, np.pi/2, np.pi*2/3)
        bz = s.BrillouinZone(r)
        # Q = (np.random.rand(1000, 3)-0.5) * 2 # this is a uniform distribution over [-1, 1)
        x=np.linspace(-1, 1, 100)
        X, Y, Z=np.meshgrid(x, x, 0)
        Q = np.stack( (X.flatten(), Y.flatten(), Z.flatten()), axis=-1)
        Qin = bz.isinside(Q)
        B = r.B
        X = np.stack( [ np.matmul(B, v) for v in Q[Qin,:] ] )
        # plot_2d_points(X)


    def test_c_moveinto_hexagonal(self):
        d, r = make_dr(3, 3, 9, np.pi/2, np.pi/2, np.pi*2/3)
        bz = s.BrillouinZone(r)
        Q = (np.random.rand(100, 3)-0.5) * 10 # this is a uniform distribution over [-5, 5)
        if np.all(bz.isinside(Q)): # this is vanishingly-unlikely
            Q += 100.0
        self.assertFalse( np.all(bz.isinside(Q)) )
        (q, tau) = bz.moveinto(Q)
        self.assertTrue( np.all(bz.isinside(q)) )
        self.assertAlmostEqual( np.abs(Q-q-tau).sum(), 0)

    def test_d_all_hallgroups(self):
        tested = 0
        failed = 0
        errored = 0
        failed_spg = []
        failed_ptg = []
        failed_lat = []
        failed_ratio = []
        errored_spg = []
        errored_ptg = []
        errored_lat = []
        errored_arg = []
        print()
        for i in range(1,531):
            spacegroup = s.Spacegroup(i)
            pointgroup = s.Pointgroup(spacegroup.pointgroup_number)
            a = 5; b = 5; c = 5; al = np.pi/2; be = np.pi/2; ga = np.pi/2
            # nothing to do for cubic spacegroups
            if 'hexa' in pointgroup.holohedry:
                ga = 2*np.pi/3
            elif 'trig' in pointgroup.holohedry:
                if 'R' in spacegroup.choice:
                    al = be = ga = np.pi/3
                else: # 'H' setting or normally hexagonal
                    c = 10;
                    ga = 2*np.pi/3
            elif 'tetr' in pointgroup.holohedry:
                c = 10
            elif 'orth' in pointgroup.holohedry:
                axperm = spacegroup.choice.replace('-','')
                if 'cab' in axperm:
                    c = 5; a = 10; b = 15;
                elif 'cba' in axperm:
                    c = 5; b = 10; a = 15;
                elif 'bca' in axperm:
                    b = 5; c = 10; a = 15;
                elif 'bac' in axperm:
                    b = 5; a = 10; c = 15;
                elif 'acb' in axperm:
                    a = 5; c = 10; b = 15;
                else:
                    a = 5; b = 10; c = 15;
            elif 'mono' in pointgroup.holohedry:
                # continue # skip all monoclinic pointgroups for now
                if 'a' in spacegroup.choice:
                    a = 10; al = np.pi/180*(91+19*np.random.rand())
                elif 'b' in spacegroup.choice:
                    b = 10; be = np.pi/180*(91+19*np.random.rand())
                elif 'c' in spacegroup.choice:
                    c = 10; ga = np.pi/180*(91+19*np.random.rand())
                else:
                    print("Monoclinic without 'a', 'b', or 'c' choice?? ", spacegroup.choice)
                    continue
            elif 'tric' in pointgroup.holohedry:
                a = 5; b = 10; c = 15;
                al = np.pi/3*(1 + np.random.rand())
                be = np.pi/3*(1 + np.random.rand())
                ga = np.pi/3*(1 + np.random.rand())

            dlat, rlat = make_dr(a, b, c, al, be, ga, i)

            # print("Hall ", i, " ", dlat)
            # print(spacegroup,pointgroup)
            try:
                bz = s.BrillouinZone(rlat)
                vol_bz = bz.polyhedron.volume
                vol_ir = bz.ir_polyhedron.volume
                tested += 1
                if not np.isclose(vol_ir, vol_bz/s.PointSymmetry(i).size):
                    # print(dlat,": ",vol_ir," != ",vol_bz/s.PointSymmetry(i).size)
                    failed += 1
                    failed_spg.append(spacegroup)
                    failed_ptg.append(pointgroup)
                    failed_lat.append(dlat)
                    failed_ratio.append(vol_ir/vol_bz*s.PointSymmetry(i).size)
            except Exception as err:
                errored += 1
                errored_spg.append(spacegroup)
                errored_ptg.append(pointgroup)
                errored_lat.append(dlat)
                errored_arg.append(err.args)


        if failed > 0:
            print("\nFailed to find irreducible Brillouin zone for",failed,"out of",tested,"(max 530) Hall groups")
            for spg, ptg, lat, rat in zip(failed_spg, failed_ptg, failed_lat, failed_ratio):
                print(spg,ptg,lat,rat)
        if errored > 0:
            print("\nException raised for",errored,"out of",tested,"(max 530) Hall Groups")
            for spg, ptg, lat, arg in zip(errored_spg, errored_ptg, errored_lat, errored_arg):
                print(arg)
                print(spg,ptg,lat)
        self.assertTrue(errored == 0)


if __name__ == '__main__':
  unittest.main()
