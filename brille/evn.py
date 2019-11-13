# Copyright 2019 Greg Tucker
#
# This file is part of brille.
#
# brille is free software: you can redistribute it and/or modify it under the
# terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# brille is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.
#
# See the GNU Affero General Public License for more details.
# You should have received a copy of the GNU Affero General Public License
# along with brille. If not, see <https://www.gnu.org/licenses/>.

"""Functions for eigenvector normalisation."""
import numpy as np


def __r_z(theta):
    c_t = np.cos(theta)
    s_t = np.sin(theta)
    return np.array([[c_t, -s_t, 0], [s_t, c_t, 0], [0, 0, 1]])


def __r_x(theta):
    c_t = np.cos(theta)
    s_t = np.sin(theta)
    return np.array([[1, 0, 0], [0, c_t, -s_t], [0, s_t, c_t]])


def __r_y(theta):
    c_t = np.cos(theta)
    s_t = np.sin(theta)
    return np.array([[c_t, 0, s_t], [0, 1, 0], [-s_t, 0, c_t]])


def __spherical_r_theta_phi(vec):
    rho = np.sqrt(np.dot(vec, vec))
    if rho > 0:
        theta = np.arccos(vec[2]/rho)
        phi = np.arctan2(vec[1], vec[0])
    else:
        # If the length of the vector is zero, there's nothing we can do
        theta = 0
        phi = 0
    return (rho, theta, phi)


def local_xyz(vec):
    """Return a local coordinate system based on the vector provided."""
    e_z = vec / np.sqrt(np.dot(vec, vec))
    e_x, e_y = __local_xy(e_z)
    return (e_x, e_y, e_z)


def __local_xy(vec):
    # Using a spherical coordinate system has the undesirable
    # quality that the local orthonormal coordinate system is discontinuous
    # at the poles, which is along the z-axis normally.
    # Since any of the axes is likely a high-symmetry direction with a high
    # chance of hosting degenerate modes, rotate the z-axis away to
    # an arbitrary (but constant) direction before determining theta and phi
    r_zx = np.matmul(__r_z(np.pi/3), __r_x(np.pi/5))
    inv_r_zx = np.matmul(__r_x(-np.pi/5), __r_z(-np.pi/3))
    _, theta, phi = __spherical_r_theta_phi(np.matmul(r_zx, vec))
    # _, theta, phi = __spherical_r_theta_phi(vec)
    sin_theta = np.sin(theta)
    cos_theta = np.cos(theta)
    sin_phi = np.sin(phi)
    cos_phi = np.cos(phi)
    e_x = np.array((-sin_phi, cos_phi, 0))
    e_y = np.array((-cos_theta*cos_phi, -cos_theta*sin_phi, sin_theta))
    return (np.matmul(inv_r_zx, e_x), np.matmul(inv_r_zx, e_y))
    # return (e_x, e_y)


def __arbitrary_xy(q_pt, e_vecs, primary_ion):
    e_x, _ = __local_xy(q_pt)
    _, n_ions, n_dims = e_vecs.shape
    if n_dims != 3:
        raise Exception("Only three dimesnional vectors supported.")
    assert primary_ion < n_ions
    e_0 = e_vecs[-1, primary_ion, :]
    v0_x = np.dot(e_x, e_0)
    if np.real(np.conj(v0_x)*v0_x) == 0:
        raise Exception("Primary ion ϵ⋅x̂ is zero!")
    phase_angle = np.arctan2(np.imag(v0_x), np.real(v0_x))
    if np.real(v0_x)*np.cos(phase_angle)+np.imag(v0_x)*np.sin(phase_angle) < 0:
        phase_angle = np.pi - phase_angle
    phase = np.cos(phase_angle) - 1j*np.sin(phase_angle)
    e_vecs *= phase
    # The only thing we could check at this point is that ϵ₀⋅x̂ is purely real.
    # Which unfortunately does not imply anything about ϵ₀⋅ŷ, ϵ₁⋅x̂, or ϵ₁⋅ŷ.
    return e_vecs


def __find_degenerate(e_vals):
    degenerate = []
    for val in e_vals:
        equiv = np.isclose(e_vals, val)
        if equiv.sum() > 1:
            degenerate.append(np.flatnonzero(equiv))
    return degenerate


def __check_and_fix(q_pt, e_vals, e_vecs, primary):
    degenerate = __find_degenerate(e_vals)
    while degenerate:
        d_set = degenerate.pop()
        if d_set.size == 2:
            e_vecs[d_set, :] = __arbitrary_xy(q_pt, e_vecs[d_set, :], primary)
        # if there are more than two degenerate modes we are likely at the
        # zone centre or boundary, in either case we need more than one point
        # to try and figure things out, so skip over this for now
    return e_vecs


def __find_primary_ion(q_pts, e_vecs):
    n_pts, n_modes, n_ions, n_dims = e_vecs.shape
    q_dot_v = q_pts.reshape(n_pts, 1, 1, n_dims) * e_vecs
    q_v_norm = np.sqrt(np.sum(q_dot_v*np.conj(q_dot_v), axis=3))
    ion_q_v = np.real(np.mean(q_v_norm.reshape(n_pts*n_modes, n_ions), axis=0))
    # The ion with smallest <|q⋅ϵ|> is most likely to have q⟂ϵ for any
    # given grid point (maybe?)
    return np.argmin(ion_q_v)


def degenerate_check(q_pts, e_vals, e_vecs, primary=None):
    """Check for degenerate eigenvalues and standardise their eigenvectors."""
    n_pts, n_modes, n_ions, n_dims = e_vecs.shape
    assert (n_pts, n_dims) == q_pts.shape
    assert (n_pts, n_modes) == e_vals.shape
    if n_ions > 1 and primary is None:
        primary = __find_primary_ion(q_pts, e_vecs)
    else:
        primary = 0
    for i in range(n_pts):
        e_vecs[i, :] = __check_and_fix(q_pts[i, :],
                                       e_vals[i, :],
                                       e_vecs[i, :],
                                       primary)
    return e_vecs


def align_eigenvectors(q_pts, e_vals, e_vecs, primary=None):
    """Align all eigenvectors such that ℑ(ϵₚ⋅x̂)≡0 at each point.

    Pick a smothly-varying local coordinate system based on q for each point
    provided and apply an arbitrary phase such that the primary ion eigenvector
    is purely real and positive along the local x direction.
    If a primary ion is not specified, pick one automatically which has the
    smallest mean displacement along q̂, and therefore is most likely to have
    a significant |ϵₚ⋅x̂|² for all q.
    """
    n_pts, n_modes, n_ions, n_dims = e_vecs.shape
    assert (n_pts, n_dims) == q_pts.shape
    assert (n_pts, n_modes) == e_vals.shape
    if n_ions > 1 and primary is None:
        primary = __find_primary_ion(q_pts, e_vecs)
    else:
        primary = 0
    for i in range(n_pts):
        e_vecs[i, :] = __arbitrary_xy(q_pts[i, :], e_vecs[i, :], primary)
    return e_vecs
