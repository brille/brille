"""
Define a class SymSim to act as the interface between SymBZ and SimPhony.

Thus enabling efficient interpolation of CASTEP-derived phonons at arbitrary
Q points.
"""
import numpy as np
import spglib

# from simphony.data.bands import BandsData
# from simphony.data.phonon import PhononData
from simphony.data.interpolation import InterpolationData
from simphony.calculate.scattering import structure_factor
from simphony import ureg    # avoid creating a second Pint UnitRegistry

from scipy import special
from scipy.stats import norm, cauchy

import symbz as sbz
from symbz.evn import degenerate_check


class SymSim:
    """
    Efficient interpolation of phonon intensity at arbitrary Q points.

    The SimPhony data classes can be used to interpolate CASTEP dynamical
    force matrices at arbitary Q points. It can then use that information to
    determine eigen values (excitation energy) and eigen vectors
    (ion displacements) for the 3×[number of ions] phonon branches. Finally
    Q and the eigen vectors can be used to determine the structure factors of
    the phonon branches, which are proportional to the doubly-differential
    cross section measured by neutron scattering.

    SymBZ provides classes to hold arbitrarily-shaped data at the points of a
    grid filling the first Brillouin zone of the primitive lattice reciprocal
    unit cell. The SymBZ classes can then use linear interpolation to estimate
    their held-data at points between grid-points, and the translational
    symmetry of the primitive lattice to find an equivalent
    first-Brillouin-zone point, q, for any arbitrary reciprocal space point, Q.

    The SymSim object uses lattice information from the SimPhony object and the
    package spglib to determine the conventional unit cell equivalent to that
    used in the CASTEP calculations. This unit cell information is then used
    to construct the primitive first Brillouin zone and a SymBZ grid. The
    gridded-points are then used by the SimPhony object to calculate ωᵢ(q)
    and ϵᵢⱼ(q), which are placed in the SymBZ object at their respective grid
    points. When an external request is made to the SymSim object to calculate
    Sᵢ(Q) it first uses the filled SymBZ object to interpolate ωᵢ(Q) and ϵᵢⱼ(Q)
    and then the SimPhony object to convert Q and ϵᵢⱼ(Q) into Sᵢ(Q).
    The more-likely use for SymSim, however, will be in calculating S(Q,ω) in
    which case ωᵢ(Q) and Sᵢ(Q) are calculated as above and then a simple
    distribution (selected by the caller) is used to broaden each of the i
    phonon branches before combining their intensities.
    """

    # pylint: disable=r0913,r0914
    def __init__(self, SPData,
                 scattering_lengths=None, cell_is_primitive=None,
                 hall_number=None, parallel=False, **kwds):
        """Initialize a new SymSim object from an existing SimPhony object."""
        if not isinstance(SPData, InterpolationData):
            msg = "Unexpected data type {}, expect failures."
            print(msg.format(type(SPData)))
        if not isinstance(scattering_lengths, dict):
            scattering_lengths = {k: 1 for k in np.unique(SPData.ion_type)}
        self.data = SPData
        self.scattering_lengths = scattering_lengths
        self.hall_number = hall_number
        self.data_cell_is_primitive = cell_is_primitive
        self.__check_if_cell_is_primitive()
        # Construct the BZGrid, by default using the conventional unit cell
        grid_q = self.__make_grid(**kwds)
        # Calculate ωᵢ(Q) and ⃗ϵᵢⱼ(Q), and fill the BZGrid:
        # Select only those keyword arguments which SimPhony expects:
        cfp_keywords = ('asr', 'precondition', 'set_attrs', 'dipole',
                        'eta_scale', 'splitting')
        cfp_dict = {k: kwds[k] for k in cfp_keywords if k in kwds}
        freq, vecs = self.data.calculate_fine_phonons(grid_q, **cfp_dict)
        n_pt = grid_q.shape[0]
        n_br = self.data.n_branches
        n_io = self.data.n_ions
        vecs = degenerate_check(grid_q, freq.magnitude, vecs)
        frqs_vecs = np.concatenate(
            ((freq.magnitude).reshape((n_pt, n_br, 1)),
             vecs.reshape(n_pt, n_br, 3*n_io)), axis=2)
        # fill requires input shaped (n_pt, n_br, [anything])
        # or (n_pt, [anything]) if n_br==1. So frqs_vecs is fine as is.
        self.grid.fill(frqs_vecs,
                       scalar_elements=1,
                       eigenvector_elements=3*n_io)
        # self.sort_branches()
        self.parallel = parallel

    def sort_branches(self, energy_weight=1.0, eigenvector_weight=1.0, weight_function=0):
        """Sort the phonon branches stored at all mapped grip points.

        By comparing the difference in phonon branch energy and the angle
        between the branch eigenvectors it is possible to determine a cost
        matrix for assigning the branches on one grid point to those on a
        neighbouring grid point. The Munkres' Assignment algorithm is then used
        to determine a local branch permutation, which is ultimately used in
        determining a global branch permutation for each grid point.

        The cost for each branch-branch assignment is the weighted sum of the
        difference in eigen energies and the angle between eigenvectors:

            Cᵢⱼ = [energy_weight]*√(ωᵢ-ωⱼ)²
                + [angle_weight]*acos(<ϵᵢ,ϵⱼ>/|ϵᵢ||ϵⱼ|)

        The weights are both one by default but can be modified as necessary.
        """
        # The input to sort_perm indicates what weight should be given to
        # each part of the resultant cost matrix. In this case, each phonon
        # branch consists of one energy, n_ions three-vectors, and no matrix;
        # perm = self.grid.centre_sort_perm(
        perm = self.grid.multi_sort_perm(
            scalar_cost_weight=energy_weight,
            eigenvector_cost_weight=eigenvector_weight,
            vector_cost_weight=0,
            matrix_cost_weight=0,
            eigenvector_weight_function=weight_function)
        frqs_vecs = np.array([x[y, :] for (x, y) in zip(self.grid.data, perm)])
        self.grid.fill(frqs_vecs,
                       scalar_elements=1,
                       eigenvector_elements=3*self.data.n_ions)
        return frqs_vecs

    def __check_if_cell_is_primitive(self):
        # CASTEP can calculate in a standard, non-standard, primitive, or
        # non-primitive cell. Try to be clever about how its stored cell
        # relates to the one a calling function will use.
        _, ion_indexes = np.unique(self.data.ion_type, return_inverse=True)
        lattice_vectors = (self.data.cell_vec.to('angstrom')).magnitude
        cell = (lattice_vectors, self.data.ion_r, ion_indexes)
        # spglib can come to our rescue here a bit
        prim_lv, _, _ = spglib.find_primitive(cell)
        # if (prim_lv == lattice_vectors).all():
        # pylint: disable=no-member
        if sbz.Direct(lattice_vectors).isapprox(sbz.Direct(prim_lv), 1e-10):
            # check for approximate equivalence of lattice parameters instead
            # of equivalent lattice vectors as spglib may introduce a rotation
            if self.data_cell_is_primitive is None:
                self.data_cell_is_primitive = True
            if not self.data_cell_is_primitive:
                msg = "Cell vectors and positions seem to represent a "
                msg += "primitive lattice but a set flag indicates otherwise."
                raise Exception(msg)
        else:
            # eventually we could try to get *really* clever and look for
            # rotation/transformation matrices which SymBZ doesn't know about
            if self.data_cell_is_primitive is None:
                self.data_cell_is_primitive = False
            if self.data_cell_is_primitive:
                msg = "Cell vectors and positions seem not to represent a "
                msg += "primitive lattice but a set flag indicates otherwise."
                raise Exception(msg)
        # while we're here with cell defined, check for its Hall number:
        if self.hall_number is None:
            symmetry_data = spglib.get_symmetry_dataset(cell)
            self.hall_number = symmetry_data['hall_number']

    # pylint: disable=no-member
    def __get_primitive_transform(self):
        return sbz.PrimitiveTransform(self.hall_number)

    # pylint: disable=c0103,w0613,no-member
    def __make_grid(self, halfN=None, step=None, units='rlu', **kwds):
        prim_tran = self.__get_primitive_transform()
        lattice_vectors = (self.data.cell_vec.to('angstrom')).magnitude
        # And we can check whether there's anything to do using SymBZ
        if prim_tran.does_anything and self.data_cell_is_primitive:
            # this multiplication is backwards compared to, e.g.,
            # https://atztogo.github.io/spglib/definition.html#space-group-operation-and-change-of-basis
            # where the lattice_vectors form the columns of a matrix and here
            # they form the rows
            lattice_vectors = np.matmul(prim_tran.invP.transpose(1, 0), lattice_vectors)
        #
        dlat = sbz.Direct(lattice_vectors, self.hall_number)
        rlat = dlat.star()
        brillouin_zone = sbz.BrillouinZone(rlat)
        if halfN is not None:
            if isinstance(halfN, (tuple, list)):
                halfN = np.array(halfN)
            if not isinstance(halfN, np.ndarray):
                raise Exception("halfN must be a tuple, list, or ndarray")
            halfN = halfN.astype('uint64')
            self.grid = sbz.BZGridQcomplex(brillouin_zone, halfN)
        elif step is not None:
            if isinstance(step, ureg.Quantity):
                step = (step.to('angstrom')).magnitude
                isrlu = False
            else:
                isrlu = units.lower() == 'rlu'
            self.grid = sbz.BZGridQcomplex(brillouin_zone, step, isrlu)
        else:
            raise Exception("You must provide a halfN or step keyword")
        # We need to make sure that we pass gridded Q points in the primitive
        # lattice, since that is what SimPhony expects:
        grid_q = self.grid.mapped_rlu
        if prim_tran.does_anything and self.data_cell_is_primitive:
            grid_q = np.array([np.matmul(prim_tran.P, x) for x in grid_q])
        return grid_q

    def s_q(self, q_hkl, **kwargs):
        """Calculate Sᵢ(Q) where Q = (q_h,q_k,q_l)."""
        self.w_q(q_hkl, **kwargs)
        # Finally calculate Sᵢ(Q)
        # using SymPhony.calculate.scattering.structure_factor
        # which only allows a limited number of keyword arguments
        sf_keywords = ('T', 'scale', 'dw_seed', 'dw_grid', 'calc_bose')
        sf_dict = {k: kwargs[k] for k in sf_keywords if k in kwargs}
        return structure_factor(self.data, self.scattering_lengths, **sf_dict)

    def w_q(self, q_pt, interpolate=True, moveinto=True, **kwargs):
        """Calculate ωᵢ(Q) where Q = (q_h,q_k,q_l)."""
        prim_tran = self.__get_primitive_transform()
        if interpolate:
            # Interpolate the previously-stored eigen values for each Q
            # each grid point has a (n_br, 1+3*n_io) array and interpolate_at
            # returns an (n_pt, n_br, 1+3*n_io) array.
            frqs_vecs = self.grid.interpolate_at(q_pt, moveinto, self.parallel)
            # Separate the frequencies and eigenvectors
            # pylint: disable=w0632
            frqs, vecs = np.split(frqs_vecs, (1,), axis=2)
            # And reshape to what SimPhony expects
            n_pt = q_pt.shape[0]
            n_br = self.data.n_branches
            n_io = self.data.n_ions
            frqs = np.ascontiguousarray(frqs.reshape(n_pt, n_br))
            vecs = np.ascontiguousarray(vecs.reshape(n_pt, n_br, n_io, 3))
            # Store all information in the SimPhony Data object
            self.data.n_qpts = n_pt
            if prim_tran.does_anything and self.data_cell_is_primitive:
                q_pt = np.array([np.matmul(prim_tran.P, x) for x in q_pt])
            self.data.qpts = q_pt
            self.data.freqs = frqs*self.data.freqs.units
            self.data.eigenvecs = vecs
        else:
            if prim_tran.does_anything and self.data_cell_is_primitive:
                q_pt = np.array([np.matmul(prim_tran.P, x) for x in q_pt])
            cfp_keywords = ('asr', 'precondition', 'set_attrs', 'dipole',
                            'eta_scale', 'splitting')
            cfp_dict = {k: kwargs[k] for k in cfp_keywords if k in kwargs}
            self.data.calculate_fine_phonons(q_pt, **cfp_dict)
        return self.data.freqs

    def frqs_vecs(self, q_hkl, **kwargs):
        """Calculate and return ωᵢ(Q) and ϵᵢ(Q)."""
        self.w_q(q_hkl, **kwargs)
        return (self.data.freqs, self.data.eigenvecs)

    def s_qw(self, q_hkl, energy, p_dict):
        """Calculate S(Q,E) for Q = (q_h, q_k, q_l) and E=energy.

        The last input, p_dict, should be a dict with keys 'resfun' and 'param'
        controlling the phonon-linewidth.

        ========================== ===================== ================
        Linewidth function           'resfun' (one of)       'param'
        ========================== ===================== ================
        Simple Harmonic Oscillator 'sho', 's'                 fwhm
        Gaussian                   'gauss', 'g'               fwhm
        Lorentzian                 'lorentz', 'lor', 'l'      fwhm
        Voigt                      'voi', 'v'            [g_fwhm, l_fwhm]
        ========================== ===================== ================

        For each linewidth function, the full name is also a valid value for
        'resfun', e.g., 'resfun':'Simple Harmonic Oscillator'.
        Functions taking a single 'param' value will use the first element in
        any non-scalar value.
        The Simple Harmonic Oscillator function looks for an additional key,
        'T', in p_dict to optionally include the temperature.

        Additional keys in p_dict are allowed and are passed on to SimPhony
        as keyword arguments to simphony.calculate.scattering.structure_factor.

        """
        res_par_tem = ('delta',)
        if 'resfun' in p_dict and 'param' in p_dict:
            res_par_tem = (p_dict['resfun'].replace(' ', '').lower(),
                           p_dict['param'])
        n_pt = energy.size
        n_br = self.data.n_branches
        # Check if we might perform the Bose factor correction twice:
        # Replicate SimPhony's behaviour of calc_bose=True by default,
        # and T=5 by default.
        if res_par_tem[0] in ('s', 'sho', 'simpleharmonicoscillator'):
            # pull out T, or 5 if it's not present
            temp_k = p_dict['T'] if 'T' in p_dict else 5
            # keep T if calc_bose is present and True
            # discard T if calc_bose is present and False
            if 'calc_bose' in p_dict:
                temp_k = temp_k if p_dict['calc_bose'] else 0
            # Prevent SimPhony from performing the Bose correction twice
            p_dict['calc_bose'] = False
            res_par_tem = (*res_par_tem, temp_k)
        # Calculate Sᵢ(Q) after interpolating ωᵢ(Q) and ⃗ϵᵢⱼ(Q)
        if 'unique_q' in p_dict and p_dict['unique_q']:
            # Avoid repeated Q entries for, e.g., (Q,E) maps
            # Finding unique points is 𝒪(q_hkl.shape[0])
            q_hkl, u_inv = np.unique(q_hkl, return_inverse=True, axis=0)
            s_i = self.s_q(q_hkl, **p_dict)[u_inv]
            omega = (self.data.freqs.to('millielectron_volt')).magnitude[u_inv]
        else:
            s_i = self.s_q(q_hkl, **p_dict)
            omega = (self.data.freqs.to('millielectron_volt')).magnitude
        # The resulting array should be (n_pt,n_br)
        if s_i.shape[0] != n_pt or s_i.shape[1] != n_br:
            msg = "Expected S(Q) shape ({}, {}) but got {}."
            msg = msg.format(n_pt, n_br, s_i.shape)
            raise Exception(msg)

        shapein = energy.shape
        energy = energy.flatten()[:, None]
        # ωᵢ(Q)  is (n_pt,n_br)
        # Sᵢ(Q)  is (n_pt,n_br)
        # energy is (n_pt,1) [instead of (n_pt,)]

        # Broaden and then sum over the n_br branches
        s_q_e = broaden_modes(energy, omega, s_i, res_par_tem).sum(1)
        if s_q_e.shape != shapein:
            s_q_e = s_q_e.reshape(shapein)
        return s_q_e


def broaden_modes(energy, omega, s_i, res_par_tem):
    """Compute S(Q,E) for a number of dispersion relations and intensities.

    Given any number of dispersion relations, ω(Q), and the intensities of the
    modes which they represent, S(Q), plus energy-broadening information in
    the form of a function name plus parameters (if required), calculate S(Q,E)
    at the provided energy positions.

    The energy positions must have shape (Npoints,).
    The dispersion and intensities must have been precalculated and should have
    shape similar to (Npoints, Nmodes). This function calls one of five
    available broadening functions, a simple harmonic oscillator, gaussian,
    lorentzian, voigt, or delta function.
    The retuned S(Q,E) array will have shape (Npoints, Nmodes).
    """
    if res_par_tem[0] in ('s', 'sho', 'simpleharmonicoscillator'):
        s_q_e = sho(energy, omega, s_i, res_par_tem[1], res_par_tem[2])
    elif res_par_tem[0] in ('g', 'gauss', 'gaussian'):
        s_q_e = gaussian(energy, omega, s_i, res_par_tem[1])
    elif res_par_tem[0] in ('l', 'lor', 'lorentz', 'lorentzian'):
        s_q_e = lorentzian(energy, omega, s_i, res_par_tem[1])
    elif res_par_tem[0] in ('v', 'voi', 'voigt'):
        s_q_e = voigt(energy, omega, s_i, res_par_tem[1])
    elif res_par_tem[0] in ('d', 'del', 'delta'):
        s_q_e = delta(energy, omega, s_i)
    else:
        print("Unknown function {}".format(res_par_tem[0]))
        s_q_e = s_i
    return s_q_e


def delta(x_0, x_i, y_i):
    """
    Compute the δ-function.

    y₀ = yᵢ×δ(x₀-xᵢ)
    """
    y_0 = np.zeros(y_i.shape, dtype=y_i.dtype)
    # y_0 = np.zeros_like(y_i)
    y_0[x_0 == x_i] = y_i[x_0 == x_i]
    return y_0


def gaussian(x_0, x_i, y_i, fwhm):
    """Compute the normal distribution with full-width-at-half-maximum fwhm."""
    if not np.isscalar(fwhm):
        fwhm = fwhm[0]
    sigma = fwhm/np.sqrt(np.log(256))
    z_0 = (x_0-x_i)/sigma
    y_0 = norm.pdf(z_0) * y_i
    return y_0


def lorentzian(x_0, x_i, y_i, fwhm):
    """Compute the Cauchy distribution with full-width-at-half-maximum fwhm."""
    if not np.isscalar(fwhm):
        fwhm = fwhm[0]
    gamma = fwhm/2
    z_0 = (x_0-x_i)/gamma
    y_0 = cauchy.pdf(z_0) * y_i
    return y_0


def voigt(x_0, x_i, y_i, params):
    """Compute the convolution of a normal and Cauchy distribution.

    The Voigt function is the exact convolution of a normal distribution (a
    Gaussian) with full-width-at-half-max gᶠʷʰᵐ and a Cauchy distribution
    (a Lorentzian) with full-with-at-half-max lᶠʷʰᵐ. Computing the Voigt
    function exactly is computationally expensive, but it can be approximated
    to (almost always nearly) machine precision quickly using the [Faddeeva
    distribution](http://ab-initio.mit.edu/wiki/index.php/Faddeeva_Package).

    The Voigt distribution is the real part of the Faddeeva distribution,
    given an appropriate rescaling of the parameters. See, e.g.,
    https://en.wikipedia.org/wiki/Voigt_profile.
    """
    if np.isscalar(params):
        g_fwhm = params
        l_fwhm = 0
    else:
        g_fwhm = params[0]
        l_fwhm = params[1]
    if l_fwhm == 0:
        return gaussian(x_0, x_i, y_i, g_fwhm)
    if g_fwhm == 0:
        return lorentzian(x_0, x_i, y_i, l_fwhm)

    area = np.sqrt(np.log(2)/np.pi)
    gamma = g_fwhm/2
    real_z = np.sqrt(np.log(2))*(x_0-x_i)/gamma
    imag_z = np.sqrt(np.log(2))*np.abs(l_fwhm/g_fwhm)
    # pylint: disable=no-member
    y_0 = area*np.real(special.wofz(real_z + 1j*imag_z))/gamma
    return y_0


def sho(x_0, x_i, y_i, fwhm, t_k):
    """Compute the Simple-Harmonic-Oscillator distribution."""
    # (partly) ensure that all inputs have the same shape:
    if np.isscalar(fwhm):
        fwhm = fwhm * np.ones(y_i.shape)
    if np.isscalar(t_k):
        t_k = t_k * np.ones(y_i.shape)
    if x_0.ndim < x_i.ndim or (x_0.shape[1] == 1 and x_i.shape[1] > 1):
        x_0 = np.repeat(x_0, x_i.shape[1], 1)
    # include the Bose factor if the temperature is non-zero
    bose = x_0 / (1-np.exp(-11.602*x_0/t_k))
    bose[t_k == 0] = 1.0
    # We need x₀² the same shape as xᵢ
    x_02 = x_0**2
    # and to ensure that only valid (finite) modes are included
    flag = (x_i != 0) * np.isfinite(x_i)
    # create an output array
    y_0 = np.zeros(y_i.shape)
    # flatten everything so that we can use logical indexing
    # keeping the original output shape
    outshape = y_0.shape
    bose = bose.flatten()
    fwhm = fwhm.flatten()
    y_0 = y_0.flatten()
    x_i = x_i.flatten()
    y_i = y_i.flatten()
    x_02 = x_02.flatten()
    flag = flag.flatten()
    # and actually calculate the distribution
    part1 = bose[flag]*(4/np.pi)*fwhm[flag]*x_i[flag]*y_i[flag]
    part2 = ((x_02[flag]-x_i[flag]**2)**2 + 4*fwhm[flag]**2*x_02[flag])
    y_0[flag] = part1/part2
    return y_0.reshape(outshape)
