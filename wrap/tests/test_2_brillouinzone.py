#!/usr/bin/env python3
"""Run tests of Brillouin zone creation."""
import unittest
import numpy as np
from importlib.util import find_spec
import brille as br_py


HASMPL = find_spec('matplotlib') is not None
HASMPL &= find_spec('mpl_toolkits') is not None
# protect against trying to load a submodule of a non-existant module
if HASMPL:
    HASMPL &= find_spec('mpl_toolkits.mplot3d') is not None
if HASMPL:
    import matplotlib.pyplot as pp
    from mpl_toolkits.mplot3d.axes3d import Axes3D


# pylint: disable=c0103
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


def make_dr(a, b, c, al=np.pi / 2, be=np.pi / 2, ga=np.pi / 2):
    """Make a Direct and Reciprocal lattice from Direct lattice parameters."""
    return br_py.Lattice(([a, b, c], [al, be, ga]))


def norm(x):
    return np.sqrt(np.dot(x, x))


def vector_lists_match(A, B):
    pA, pB = np.where(np.isclose(np.abs(A[:, np.newaxis, :] - B).sum(axis=2), 0.))
    if pA.size != A.shape[0] or pB.size != B.shape[0]:
        return False
    return np.allclose(A[pA], B[pB])


class BrillouinZone(unittest.TestCase):
    def test_a_init_unit_cube(self):
        # instantiate a cubic lattice with unit-length vectors
        # and its reciprocal lattice, still cubic with 2π-length vectors
        lat = make_dr(1, 1, 1)
        # creating a BrillouinZone objects requires that we pass the lattice
        # its first Brillouin zone is a cube in reciprocal space
        bz = br_py.BrillouinZone(lat)
        # with six faces, defined by: (00̄1),(0̄10),(̄100),(100),(010),(001)
        p = bz.points
        self.assertEqual(p.ndim, 2)
        self.assertEqual(p.shape[0], 6)  # and there is one for each of the six faces
        self.assertEqual(p.shape[1], 3)  # the face vectors are 3-vectors
        n = bz.normals
        self.assertEqual(n.ndim, 2)
        self.assertEqual(n.shape[0], 6)  # and there is one for each of the six faces
        self.assertEqual(n.shape[1], 3)  # the face normals are 3-vectors
        n_dot_p = np.array([np.dot(x / norm(x), y / norm(y)) for x, y in zip(n, p)])
        self.assertTrue((n_dot_p == 1.0).all())

        expected = np.array([[-1, 0, 0], [0, -1, 0], [0, 0, -1], [0, 0, 1], [0, 1, 0], [1, 0, 0]])
        self.assertTrue(vector_lists_match(p, expected / 2))
        self.assertTrue(vector_lists_match(np.array([x / norm(x) for x in n]), expected))
        #
        # the vertices of the first Brillouin zone are the 8 corners of the 1/2-unit cube:
        expected = np.array([[-1, -1, -1], [-1, -1, 1], [-1, 1, -1], [-1, 1, 1],
                             [1, -1, -1], [1, -1, 1], [1, 1, -1], [1, 1, 1]]) / 2
        verts = bz.vertices
        self.assertEqual(verts.ndim, 2)
        self.assertEqual(verts.shape[0], 8)
        self.assertEqual(verts.shape[1], 3)
        self.assertTrue(vector_lists_match(verts, expected))

    def test_b_isinside_unit_cube(self):
        lat = make_dr(1, 1, 1)
        bz = br_py.BrillouinZone(lat)
        # the first Brillouin zone of this unit-cube is bounded by
        # {(00̄1),(0̄10),(̄100),(100),(010),(001)}/2
        face_centres = np.array([[-1, 0, 0], [0, -1, 0], [0, 0, -1], [0, 0, 1], [0, 1, 0], [1, 0, 0]],
                                dtype='double') / 2
        self.assertTrue(np.all(bz.isinside(face_centres)))
        # including the vertices (since they are the corners of the zone)
        corners = np.array([[-1, -1, -1], [-1, -1, 1], [-1, 1, -1], [-1, 1, 1],
                            [1, -1, -1], [1, -1, 1], [1, 1, -1], [1, 1, 1]], dtype='double') / 2
        self.assertTrue(np.all(bz.isinside(corners)))
        # so all points with h, k, and l in the range [-0.5, 0.5] are in the zone
        Q = np.random.rand(100, 3) - 0.5  # this is a uniform distribution over [-0.5, 0.5) -- close enough
        self.assertTrue(np.all(bz.isinside(Q)))
        self.assertFalse(np.all(bz.isinside(Q + 5)))

    def test_c_moveinto_unit_cube(self):
        lat = make_dr(1, 1, 1)
        bz = br_py.BrillouinZone(lat)
        Q = (np.random.rand(100, 3) - 0.5) * 10  # this is a uniform distribution over [-5, 5)
        if np.all(bz.isinside(Q)):  # this is vanishingly-unlikely
            Q += 100.0
        self.assertFalse(np.all(bz.isinside(Q)))
        (q, tau) = bz.moveinto(Q)
        self.assertTrue(np.all(bz.isinside(q)))
        self.assertAlmostEqual(np.abs(Q - q - tau).sum(), 0)

    #
    #
    #
    def test_a_init_hexagonal(self):
        # instantiate a hexagonal lattice and its reciprocal lattice, still hexagonal
        lat = make_dr(3, 3, 9, np.pi / 2, np.pi / 2, np.pi * 2 / 3)
        # creating a BrillouinZone objects requires that we pass the lattice object
        bz = br_py.BrillouinZone(lat)
        # with eight faces, defined by: (̄100),(̄110),(0̄10),(001),(001),(010),(1̄10),(100)
        p = bz.points
        self.assertEqual(p.ndim, 2)
        self.assertEqual(p.shape[0], 8)  # and there is one for each of the eight faces
        self.assertEqual(p.shape[1], 3)  # the face vectors are 3-vectors
        n = bz.normals
        self.assertEqual(n.ndim, 2)
        self.assertEqual(n.shape[0], 8)  # and there is one for each of the six faces
        self.assertEqual(n.shape[1], 3)  # the face normals are 3-vectors
        n_dot_p = np.array([np.dot(x / norm(x), y / norm(y)) for x, y in zip(n, p)])
        self.assertTrue(np.allclose(n_dot_p, 1.))

        expected = np.array(
            [[-1, 0, 0], [-1, 1, 0], [0, -1, 0], [0, 0, -1], [0, 0, 1], [0, 1, 0], [1, -1, 0], [1, 0, 0]])
        self.assertTrue(vector_lists_match(p, expected / 2))
        n_compare = np.array([x / np.max(np.abs(x)) for x in n])
        self.assertTrue(vector_lists_match(n_compare, expected))
        self.assertTrue(np.all(bz.isinside(expected / 2)))
        #
        # the vertices of the first Brillouin zone are the 12 corners of the hexagonal-prism:
        expected = np.array([[-4, 2, -3], [-2, -2, -3], [-4, 2, 3], [-2, -2, 3], [-2, 4, -3], [-2, 4, 3],
                             [2, -4, -3], [2, -4, 3], [2, 2, -3], [4, -2, -3], [2, 2, 3], [4, -2, 3]]) / 6
        verts = bz.vertices
        self.assertEqual(verts.ndim, 2)
        self.assertEqual(verts.shape[0], 12)
        self.assertEqual(verts.shape[1], 3)
        self.assertTrue(vector_lists_match(verts, expected))
        self.assertTrue(np.all(bz.isinside(expected)))

    def test_b_isinside_hexagonal(self):
        lat = make_dr(3, 3, 9, np.pi / 2, np.pi / 2, np.pi * 2 / 3)
        bz = br_py.BrillouinZone(lat)
        # Q = (np.random.rand(1000, 3)-0.5) * 2 # this is a uniform distribution over [-1, 1)
        x = np.linspace(-1, 1, 100)
        X, Y, Z = np.meshgrid(x, x, 0)
        Q = np.stack((X.flatten(), Y.flatten(), Z.flatten()), axis=-1)
        Qin = bz.isinside(Q)
        B = lat.reciprocal_vectors
        X = np.stack([np.matmul(B, v) for v in Q[Qin, :]])
        # plot_2d_points(X)


    def test_c_moveinto_hexagonal(self):
        lat = make_dr(3, 3, 9, np.pi / 2, np.pi / 2, np.pi * 2 / 3)
        bz = br_py.BrillouinZone(lat)
        Q = (np.random.rand(100, 3) - 0.5) * 10  # this is a uniform distribution over [-5, 5)
        if np.all(bz.isinside(Q)):  # this is vanishingly-unlikely
            Q += 100.0
        self.assertFalse(np.all(bz.isinside(Q)))
        (q, tau) = bz.moveinto(Q)
        self.assertTrue(np.all(bz.isinside(q)))
        self.assertAlmostEqual(np.abs(Q - q - tau).sum(), 0)


    # This test can not work since I've removed the Hall-number symmetry specification constructor from Lattice
    # construction. ... but not from Spacegroup!
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
        for i in range(1, 531):
            spacegroup = br_py.Spacegroup(i)
            pointgroup = br_py.Pointgroup(spacegroup.pointgroup_number)
            a, b, c, al, be, ga = 5, 5, 5, np.pi/2, np.pi/2, np.pi/2
            # nothing to do for cubic spacegroups
            if 'hexa' in pointgroup.holohedry:
                ga = 2 * np.pi / 3
            elif 'trig' in pointgroup.holohedry:
                if 'R' in spacegroup.choice:
                    al = be = ga = np.pi / 3
                else:  # 'H' setting or normally hexagonal
                    c, ga = 10, 2 * np.pi / 3
            elif 'tetr' in pointgroup.holohedry:
                c = 10
            elif 'orth' in pointgroup.holohedry:
                axperm = spacegroup.choice.replace('-', '')
                if 'cab' in axperm:
                    c, a, b = 5, 10, 15
                elif 'cba' in axperm:
                    c, b, a = 5, 10, 15
                elif 'bca' in axperm:
                    b, c, a = 5, 10, 15
                elif 'bac' in axperm:
                    b, a, c = 5, 10, 15
                elif 'acb' in axperm:
                    a, c, b = 5, 10, 15
                else:
                    a, b, c = 5, 10, 15
            elif 'mono' in pointgroup.holohedry:
                # continue # skip all monoclinic pointgroups for now
                if 'a' in spacegroup.choice:
                    a, al = 10, np.pi / 180 * (91 + 19 * np.random.rand())
                elif 'b' in spacegroup.choice:
                    b, be = 10, np.pi / 180 * (91 + 19 * np.random.rand())
                elif 'c' in spacegroup.choice:
                    c, ga = 10, np.pi / 180 * (91 + 19 * np.random.rand())
                else:
                    print("Monoclinic without 'a', 'b', or 'c' choice?? ", spacegroup.choice)
                    continue
            elif 'tric' in pointgroup.holohedry:
                ang = lambda: np.pi / 3 * (1 + np.random.rand())
                a, b, c, al, be, ga = 5, 10, 15, ang(), ang(), ang()

            lat = br_py.Lattice(([a, b, c], [al, be, ga]), spacegroup=spacegroup.hall_symbol)

            # print("Hall ", i, " ", dlat)
            # print(spacegroup,pointgroup)
            try:
                bz = br_py.BrillouinZone(lat)
                vol_bz = bz.polyhedron.volume
                vol_ir = bz.ir_polyhedron.volume
                tested += 1
                if not np.isclose(vol_ir, vol_bz / br_py.PointSymmetry(i).size):
                    # print(dlat,": ",vol_ir," != ",vol_bz/br_py.PointSymmetry(i).size)
                    failed += 1
                    failed_spg.append(spacegroup)
                    failed_ptg.append(pointgroup)
                    failed_lat.append(lat)
                    failed_ratio.append(vol_ir / vol_bz * br_py.PointSymmetry(i).size)
            except Exception as err:
                errored += 1
                errored_spg.append(spacegroup)
                errored_ptg.append(pointgroup)
                errored_lat.append(lat)
                errored_arg.append(err.args[0])

        if failed > 0:
            print("\nFailed to find irreducible Brillouin zone for", failed, "out of", tested, "(max 530) Hall groups")
            for spg, ptg, lat, rat in zip(failed_spg, failed_ptg, failed_lat, failed_ratio):
                print(spg, ptg, lat, rat)
        if errored > 0:
            print("\nException raised for", errored, "out of", tested, "(max 530) Hall Groups")
            for spg, ptg, lat, arg in zip(errored_spg, errored_ptg, errored_lat, errored_arg):
                print(arg)
                print(spg, ptg, lat)
        self.assertTrue(errored == 0)


if __name__ == '__main__':
    unittest.main()
