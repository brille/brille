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
        # Using a spherical coordinate system has the undesirable
        # quality that the local orthonormal coordinate system is discontinuous
        # at the poles, which is along the z-axis normally.
        # Since any of the axes is likely a high-symmetry direction with a high
        # chance of hosting degenerate modes, rotate the z-axis away to
        # an arbitrary (but constant) direction
        vec = np.matmul(__r_z(np.pi/3), np.matmul(__r_x(np.pi/5), vec))
        theta = np.arccos(vec[2]/rho)
        phi = np.arctan2(vec[1], vec[0])
    else:
        # If the length of the vector is zero, there's nothing we can do
        theta = 0
        phi = 0
    return (rho, theta, phi)


def __local_xy(vec):
    _, theta, phi = __spherical_r_theta_phi(vec)
    sin_theta = np.sin(theta)
    cos_theta = np.cos(theta)
    sin_phi = np.sin(phi)
    cos_phi = np.cos(phi)
    e_x = np.array((-sin_phi, cos_phi, 0))
    e_y = np.array((cos_theta*cos_phi, cos_theta*sin_phi, -sin_theta))
    return (e_x, e_y)


def __set_arbitrary_xy(q_pt, e_vecs, primary_ion):
    e_x, e_y = __local_xy(q_pt)
    n_vecs, n_ions, n_dims = e_vecs.shape
    if n_vecs != 2:
        raise Exception("Only two degenerate modes allowed.")
    if n_dims != 3:
        raise Exception("Only three dimesnional vectors supported.")
    assert primary_ion < n_ions
    v0_x = np.dot(e_x, e_vecs[0, primary_ion, :])
    if np.real(np.conj(v0_x)*v0_x) == 0:
        raise Exception("Primary ion ϵ⋅x̂ is zero!")
    phase_angle = np.arctan2(np.imag(v0_x), np.real(v0_x))
    phase = np.cos(phase_angle) - 1j*np.sin(phase_angle)
    e_vecs *= phase
    v0_y = np.dot(e_y, e_vecs[0, primary_ion, :])
    is_ok = np.isclose(np.real(np.conj(v0_y)*v0_y), 0)
    v1_x = np.dot(e_x, e_vecs[1, primary_ion, :])
    is_ok &= np.isclose(np.real(np.conj(v1_x)*v1_x), 0)
    v1_y = np.dot(e_y, e_vecs[1, primary_ion, :])
    is_ok &= not np.isclose(np.real(np.conj(v1_y)*v1_y), 0)
    if not is_ok:
        raise Exception("Something wrong with applying phase factor")
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
        if d_set.size > 2:
            raise Exception("More than two degenerate modes. Extend my code!")
        e_vecs[d_set, :] = __set_arbitrary_xy(q_pt, e_vecs[d_set, :], primary)
    return e_vecs


def __find_primary_ion(q_pts, e_vecs):
    n_pts, n_modes, n_ions, n_dims = e_vecs.shape
    q_dot_v = q_pts.reshape(n_pts, 1, 1, n_dims) * e_vecs
    q_v_norm = np.sqrt(np.sum(q_dot_v*np.conj(q_dot_v), axis=3))
    ion_q_v = np.mean(q_v_norm.reshape(n_modes*n_dims, n_ions), axis=0)
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
