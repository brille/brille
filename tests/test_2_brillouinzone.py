#!/usr/bin/env python3
import os, sys, unittest
import numpy as np

# We need to tell Python where it can find the symbz module.
addpath = os.getcwd()
# It's either in the working directory where python was called or in
# a sub-directory (called Debug using Visual Studio under Windows)
if os.path.exists('Debug'):
    addpath += "\\Debug"
sys.path.append(addpath)

import symbz as s

from importlib import util
hasmpl  = util.find_spec('matplotlib') is not None
hasmpl &= util.find_spec('mpl_toolkits') is not None
# protect against trying to load a submodule of a non-existant module
if hasmpl:
    hasmpl &= util.find_spec('mpl_toolkits.mplot3d') is not None
if hasmpl:
    import matplotlib.pyplot as pp
    from mpl_toolkits.mplot3d.axes3d import Axes3D

def plot_points(x,title=''):
    if hasmpl:
        fig = pp.figure()
        ax = Axes3D(fig)
        ax.scatter(x[:,0],x[:,1],x[:,2],s=10)
        ax.set_title(title)
        pp.show()
def plot_points_with_lines(x,y,title=''):
    if hasmpl:
        fig = pp.figure()
        ax = Axes3D(fig)
        ax.plot(y[:,0],y[:,1],y[:,2])
        ax.scatter(x[:,0],x[:,1],x[:,2],s=10)
        ax.set_title(title)
        pp.show()
def plot_2d_points_with_lines(x,y,title=''):
    if hasmpl:
        ax = pp.axes()
        ax.plot(y[:,0],y[:,1])
        ax.scatter(x[:,0],x[:,1],s=10)
        ax.set_title(title)
        pp.show()


def make_dr(a,b,c,al=np.pi/2,be=np.pi/2,ga=np.pi/2):
    d = s.Direct(a,b,c,al,be,ga)
    r = d.star()
    return (d,r)


class BrillouinZone (unittest.TestCase):
    def test_a_init_unit_cube(self):
        # instantiate a cubic lattice with unit-length vectors
        # and its reciprocal lattice, still cubic with 2π-length vectors
        d,r = make_dr(1,1,1)
        # creating a BrillouinZone objects requires that we pass the reciprocal
        # lattice object, and passing a Direct object is a TypeError:
        with self.assertRaises(TypeError):
            s.BrillouinZone(d)
        # its first Brillouin zone is a cube in reciprocal space
        bz = s.BrillouinZone(r)
        # with six faces, defined by: (00̄1),(0̄10),(̄100),(100),(010),(001)
        expected_faces = np.array([ [-1,0,0], [0,-1,0], [0,0,-1], [0,0,1], [0,1,0], [1,0,0] ], dtype='int32')
        faces = bz.faces
        self.assertEqual(faces.ndim,2)
        self.assertEqual(faces.shape[0],6) # and there is one for each of the six faces
        self.assertEqual(faces.shape[1],3) # the face vectors are 3-vectors
        self.assertTrue( (faces==expected_faces).all() )
        #
        # the vertices of the first Brillouin zone are the 8 corners of the 1/2-unit cube:
        expected_verts = np.array([[-1,-1,-1], [-1,-1,1], [-1,1,-1], [-1,1,1], [1,-1,-1], [1,-1,1], [1,1,-1], [1,1,1]], dtype='double')/2
        verts = bz.vertices
        self.assertEqual(verts.ndim,2)
        self.assertEqual(verts.shape[0],8)
        self.assertEqual(verts.shape[1],3)
        self.assertAlmostEqual( np.abs(expected_verts-verts).sum(), 0)
    def test_b_isinside_unit_cube(self):
        d,r = make_dr(1,1,1)
        bz = s.BrillouinZone(r)
        # the first Brillouin zone of this unit-cube is bounded by
        # {(00̄1),(0̄10),(̄100),(100),(010),(001)}/2
        face_centres = np.array([ [-1,0,0], [0,-1,0], [0,0,-1], [0,0,1], [0,1,0], [1,0,0] ], dtype='double')/2
        self.assertTrue( (bz.isinside(face_centres)).all() )
        # including the vertices (since they are the corners of the zone)
        corners = np.array([[-1,-1,-1], [-1,-1,1], [-1,1,-1], [-1,1,1], [1,-1,-1], [1,-1,1], [1,1,-1], [1,1,1]], dtype='double')/2
        self.assertTrue( (bz.isinside(corners)).all() )
        # so all points with h, k, and l in the range [-0.5,0.5] are in the zone
        Q = np.random.rand(100,3) - 0.5 # this is a uniform distribution over [-0.5,0.5) -- close enough
        self.assertTrue( bz.isinside(Q).all() )
        self.assertFalse( bz.isinside(Q+5).all() )
    def test_c_moveinto_unit_cube(self):
        d,r = make_dr(1,1,1)
        bz = s.BrillouinZone(r)
        Q = (np.random.rand(100,3)-0.5) * 10 # this is a uniform distribution over [-5,5)
        if bz.isinside(Q).all(): # this is vanishingly-unlikely
            Q += 100.0
        self.assertFalse( bz.isinside(Q).all() )
        (q,tau) = bz.moveinto(Q)
        self.assertTrue( bz.isinside(q).all() )
        self.assertAlmostEqual( np.abs(Q-q-tau).sum(), 0)
#
#
#
    def test_a_init_hexagonal(self):
        # instantiate a hexagonal lattice and its reciprocal lattice, still hexagonal
        d,r = make_dr(3,3,9,np.pi/2,np.pi/2,np.pi*2/3)
        # creating a BrillouinZone objects requires that we pass the reciprocal
        # lattice object, and passing a Direct object is a TypeError:
        with self.assertRaises(TypeError):
            s.BrillouinZone(d)
        # its first Brillouin zone is a cube in reciprocal space
        bz = s.BrillouinZone(r)
        # with eight faces, defined by: (̄100),(̄110),(0̄10),(001),(001),(010),(1̄10),(100)
        expected_faces = np.array([ [-1,0,0],[-1,1,0],[0,-1,0],[0,0,-1],[0,0,1],[0,1,0],[1,-1,0],[1,0,0] ], dtype='int32')
        faces = bz.faces
        self.assertEqual(faces.ndim,2)
        self.assertEqual(faces.shape[0],8) # there is one for each of the eight faces
        self.assertEqual(faces.shape[1],3) # the face vectors are 3-vectors
        self.assertTrue( (faces==expected_faces).all() )
        #
        # the vertices of the first Brillouin zone are the 12 corners of the hexagonal-prism:
        expected_verts = np.array([[-4, 2,-3],[-2,-2,-3],[-4, 2, 3],[-2,-2, 3],[-2, 4,-3],[-2, 4, 3],
                                   [ 2,-4,-3],[ 2,-4, 3],[ 2, 2,-3],[ 4,-2,-3],[ 2, 2, 3],[ 4,-2, 3]], dtype='double')/6
        verts = bz.vertices
        self.assertEqual(verts.ndim,2)
        self.assertEqual(verts.shape[0],12)
        self.assertEqual(verts.shape[1],3)
        self.assertAlmostEqual( np.abs(expected_verts-verts).sum(), 0)
        self.assertTrue( (bz.isinside(expected_faces/2)).all() )
        self.assertTrue( (bz.isinside(expected_verts)).all() )

        B = r.get_B_matrix()
        print(B)
        X = np.stack( [np.matmul(B,v) for v in expected_verts])
        plot_points(X)

    def test_b_isinside_hexagonal(self):
        d,r = make_dr(3,3,9,np.pi/2,np.pi/2,np.pi*2/3)
        bz = s.BrillouinZone(r)
        # Q = (np.random.rand(1000,3)-0.5) * 2 # this is a uniform distribution over [-1,1)
        x=np.linspace(-1,1,100)
        X,Y,Z=np.meshgrid(x,x,0)
        Q = np.stack( (X.flatten(),Y.flatten(),Z.flatten()),axis=-1)
        Qin = bz.isinside(Q)
        B = r.get_B_matrix()
        X = np.stack( [ np.matmul(B,v) for v in Q[Qin,:] ] )
        plot_2d_points_with_lines(X,bz.vertices_invA)


    def test_c_moveinto_hexagonal(self):
        d,r = make_dr(3,3,9,np.pi/2,np.pi/2,np.pi*2/3)
        bz = s.BrillouinZone(r)
        Q = (np.random.rand(100,3)-0.5) * 10 # this is a uniform distribution over [-5,5)
        if bz.isinside(Q).all(): # this is vanishingly-unlikely
            Q += 100.0
        self.assertFalse( bz.isinside(Q).all() )
        (q,tau) = bz.moveinto(Q)
        self.assertTrue( bz.isinside(q).all() )
        self.assertAlmostEqual( np.abs(Q-q-tau).sum(), 0)


if __name__ == '__main__':
  unittest.main()
