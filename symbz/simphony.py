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

from scipy.special import wofz
from scipy.stats import norm, cauchy

import symbz as sbz


def delta(x_0, x_i, y_i):
    """
    Compute the δ-function.

    y₀ = yᵢ×δ(x₀-xᵢ)
    """
    y_0 = np.zeros_like(y_i)
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
    y_0 = area*np.real(wofz(real_z + 1j*imag_z))/gamma
    return y_0


def sho(x_0, x_i, y_i, fwhm, t_k):
    """Compute the Simple-Harmonic-Oscillator distribution."""
    if not np.isscalar(fwhm):
        fwhm = fwhm[0]
    bose = x_0 / (1-np.exp(-11.602*x_0/t_k)) if np.abs(t_k) > 0 else 1.0
    x_02 = np.repeat(x_0**2, x_i.shape[1], 1)
    v_fl = (x_i > 0) * np.isfinite(x_i)
    y_0 = np.zeros_like(y_i)
    y_0[v_fl] = bose*(4/np.pi)*fwhm*x_i[v_fl]*y_i[v_fl]
    y_0[v_fl] /= ((x_02[v_fl]-x_i[v_fl]**2)**2 + 4*fwhm**2*x_02[v_fl])
    return y_0


class SymSim(object):
    """
    An object to enable efficient interpolation of CASTEP-derived phonons at
    arbitrary Q points.

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
    package spglib to determine the conventional unit cell used in the CASTEP
    calculations. This unit cell information is then used to construct the
    primitive first Brillouin zone and a SymBZ grid. The gridded-points are
    then used by the SimPhony object to calculate ωᵢ(q) and ϵᵢⱼ(q), which are
    placed in the SymBZ object at their respective grid points.
    When an external request is made to the SymSim object to calculate Sᵢ(Q)
    it first uses the filled SymBZ object to interpolate ωᵢ(Q) and ϵᵢⱼ(Q) and
    then the SimPhony object to convert Q and ϵᵢⱼ(Q) into Sᵢ(Q).
    The more-likely use for SymSim, however, will be in calculating S(Q,ω) in
    which case ωᵢ(Q) and Sᵢ(Q) are calculated as above and then a simple
    distribution (selected by the caller) is used to broaden each of the i
    phonon branches before combining their intensities.

    """
    def __init__(self, SPData,
                 scattering_lengths=None,
                 parallel=False, **kwds):
        """Initialize a new SymSim object from an existing SimPhony object"""
        if not isinstance(SPData, InterpolationData):
            print(f"Unexpected data type {type(SPData)}, expect failures.")
        if not isinstance(scattering_lengths, dict):
            scattering_lengths = {k: 1 for k in np.unique(SPData.ion_type)}
        self.data = SPData
        self.scattering_lengths = scattering_lengths

        self.__make_grid(**kwds)

        # Calculate ωᵢ(Q) and ⃗ϵᵢⱼ(Q), and fill the BZGrid:
        grid_q = self.grid.mapped_rlu
        # Select only those keyword arguments which SimPhony expects:
        cfp_keywords = ('asr', 'precondition', 'set_attrs', 'dipole',
                        'eta_scale', 'splitting')
        cfp_dict = {k: kwds[k] for k in cfp_keywords}
        freq, vecs = self.data.calculate_fine_phonons(grid_q, **cfp_dict)
        n_pt = grid_q.shape[0]
        n_br = self.data.n_branches
        n_io = self.data.n_ions
        frqs_vecs = np.concatenate(
            (
                (freq.magnitude).reshape((n_pt, n_br, 1)),
                vecs.reshape(n_pt, n_br, 3*n_io)
            ),
            ax_is=2)
        self.grid.fill(frqs_vecs)
        self.parallel = parallel

    def __make_grid(self, n_half=None, step=None, units='rlu', **kwds):
        _, ion_indexes = np.unique(self.data.ion_type, return_inverse=True)
        lattice_vectors = (self.data.cell_vec.to('angstrom')).magnitude
        cell = (lattice_vectors, self.data.ion_r, ion_indexes)
        symmetry_data = spglib.get_symmetry_dataset(cell)
        #
        dlat = sbz.Direct(lattice_vectors, symmetry_data['hall_number'])
        rlat = dlat.star()
        brillouin_zone = sbz.BrillouinZone(rlat)
        if n_half is not None:
            self.grid = sbz.BZGridQcomplex(brillouin_zone, n_half)
        elif step is not None:
            if isinstance(step, ureg.Quantity):
                step = (step.to('angstrom')).magnitude
                isrlu = False
            else:
                isrlu = units.lower() == 'rlu'
            self.grid = sbz.BZGridQcomplex(brillouin_zone, step, isrlu)

    def s_q(self, q_hkl, **kwargs):
        """Calculate Sᵢ(Q) where Q = (q_h,q_k,q_l)."""
        # Interpolate the previously-stored eigen values and vectors for each Q
        frqs_vecs = self.grid.interpolate_at(q_hkl, True, self.parallel)
        # Separate them
        frqs, vecs = np.split(frqs_vecs, np.array([1]), axis=2)
        # And reshape to what SimPhony expects
        n_pt = q_hkl.shape[0]
        n_br = self.data.n_branches
        n_io = self.data.n_ions
        frqs = frqs.reshape((n_pt, n_br))
        vecs = vecs.reshape((n_pt, n_br, n_io, 3))
        # Store all information in the SimPhony Data object
        self.data.n_qpts = n_pt
        self.data.qpts = q_hkl
        self.data.freqs = frqs*self.data.freqs.units
        self.data.eigenvecs = vecs
        # Finally calculate Sᵢ(Q)
        # using SymPhony.calculate.scattering.structure_factor
        # which only allows a limited number of keyword arguments
        sf_keywords = ('T', 'scale', 'dw_seed', 'dw_grid', 'calc_bose')
        sf_dict = {k: kwargs[k] for k in sf_keywords}
        return structure_factor(self.data, self.scattering_lengths, **sf_dict)

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
        if 'resfun' in p_dict and 'param' in p_dict:
            resfun = p_dict['resfun'].replace(' ', '').lower()
            params = p_dict['param']
        else:
            resfun = 'delta'
        n_pt = energy.size
        n_br = self.data.n_branches
        # Check if we might perform the Bose factor correction twice:
        # Replicate SimPhony's behaviour of calc_bose=True by default,
        # and T=5 by default.
        if resfun in ('s', 'sho', 'simpleharmonicoscillator'):
            # pull out T, or 5 if it's not present
            temp_k = p_dict['T'] if 'T' in p_dict else 5
            # keep T if calc_bose is present and True
            # discard T if calc_bose is present and False
            if 'calc_bose' in p_dict:
                temp_k = temp_k if p_dict['calc_bose'] else 0
            # Prevent SimPhony from performing the Bose correction twice
            p_dict['calc_bose'] = False
        # Calculate Sᵢ(Q) after interpolating ωᵢ(Q) and ⃗ϵᵢⱼ(Q)
        s_i = self.s_q(q_hkl, **p_dict)
        # The resulting array should be (n_pt,n_br)
        msg = f"Expected S(Q) shape ({n_pt},{n_br}) but got {s_i.shape}."
        assert s_i.shape[0] is n_pt and s_i.shape[1] is n_br, msg

        omega = (self.data.freqs.to('millielectron_volt')).magnitude
        shapein = energy.shape
        energy = energy.flatten()[:, None]
        # ωᵢ(Q)  is (n_pt,n_br)
        # Sᵢ(Q)  is (n_pt,n_br)
        # energy is (n_pt,1) [instead of (n_pt,)]

        if resfun in ('s', 'sho', 'simpleharmonicoscillator'):
            s_q_e = sho(energy, omega, s_i, params, temp_k)
        elif resfun in ('g', 'gauss', 'gaussian'):
            s_q_e = gaussian(energy, omega, s_i, params)
        elif resfun in ('l', 'lor', 'lorentz', 'lorentzian'):
            s_q_e = lorentzian(energy, omega, s_i, params)
        elif resfun in ('v', 'voi', 'voigt'):
            s_q_e = voigt(energy, omega, s_i, params)
        elif resfun in ('d', 'del', 'delta'):
            s_q_e = delta(energy, omega, s_i)
        else:
            print(f'Unknown function {resfun}')
            s_q_e = s_i
        # Sum over the n_br branches
        s_q_e = s_q_e.sum(1)
        if s_q_e.shape != shapein:
            s_q_e = s_q_e.reshape(shapein)

        return s_q_e
