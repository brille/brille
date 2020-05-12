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
"""
Define a class BrEu to act as the interface between brille and Euphonic.

Thus enabling efficient interpolation of CASTEP-derived phonons at arbitrary
Q points.
"""
import numpy as np
import spglib

from scipy import special
from scipy.stats import norm, cauchy

from euphonic.data.interpolation import InterpolationData
from euphonic import ureg    # avoid creating a second Pint UnitRegistry

import brille as br
from brille.spglib import BrSpgl
from brille.evn import degenerate_check


class BrEu:
    """
    Efficient interpolation of phonon intensity at arbitrary Q points.

    The Euphonic data classes can be used to interpolate CASTEP dynamical
    force matrices at arbitary Q points. It can then use that information to
    determine eigen values (excitation energy) and eigen vectors
    (ion displacements) for the 3×[number of ions] phonon branches. Finally
    Q and the eigen vectors can be used to determine the structure factors of
    the phonon branches, which are proportional to the doubly-differential
    cross section measured by neutron scattering.

    brille provides classes to hold arbitrarily-shaped data at the points of a
    grid filling the first irreducible Brillouin zone of the primitive lattice
    reciprocal unit cell. The brille classes can then use linear interpolation
    to estimate their held-data at points between grid-points, the translational
    and rotational symmetry of the primitive lattice to find an equivalent
    first-irreducible-Brillouin-zone point, q, for any arbitrary reciprocal
    space point, Q.

    The BrEu object uses lattice information from the Euphonic object and the
    package spglib to determine the conventional unit cell equivalent to that
    used in, e.g., CASTEP calculations. This unit cell information is then used
    to construct the primitive first irreducible Brillouin zone and a brille
    grid. The gridded-points are then used by the Euphonic object to calculate
    ωᵢ(q) and ϵᵢⱼ(q), which are placed in the brille object at their respective
    grid points. When an external request is made to the BrEu object to
    calculate Sᵢ(Q) it first uses the filled brille object to interpolate ωᵢ(Q)
    and ϵᵢⱼ(Q) and then the Euphonic object to convert Q and ϵᵢⱼ(Q) into Sᵢ(Q).
    The more-likely use for BrEu, however, will be in calculating S(Q,ω) in
    which case ωᵢ(Q) and Sᵢ(Q) are calculated as above and then a simple
    distribution (selected by the caller) is used to broaden each of the i
    phonon branches before combining their intensities.
    """

    # pylint: disable=r0913,r0914
    def __init__(self, SPData,
                 scattering_lengths=None, cell_is_primitive=None,
                 sort=True, hall=None, parallel=False, **kwds):
        """Initialize a new BrEu object from an existing Euphonic object."""
        if not isinstance(SPData, InterpolationData):
            msg = "Unexpected data type {}, expect failures."
            print(msg.format(type(SPData)))
        if not isinstance(scattering_lengths, dict):
            scattering_lengths = {k: 1 for k in np.unique(SPData.ion_type)}
        self.data = SPData
        self.scattering_lengths = scattering_lengths
        self.brspgl = BrSpgl(SPData.cell_vec.to('angstrom').magnitude, SPData.ion_r, SPData.ion_type, hall=hall)
        # Construct the BZGrid, by default using the conventional unit cell
        grid_q = self.__define_grid_or_mesh(**kwds)
        # Calculate ωᵢ(Q) and ⃗ϵᵢⱼ(Q), and fill the BZGrid:
        # Select only those keyword arguments which Euphonic expects:
        cfp_keywords = ('asr', 'precondition', 'set_attrs', 'dipole',
                        'eta_scale', 'splitting')
        cfp_dict = {k: kwds[k] for k in cfp_keywords if k in kwds}
        if parallel:
            cfp_dict['use_c'] = True
            cfp_dict['n_threads'] = half_cpu_count()
        # calculate_fine_phonos returns the frequencies and eigenvectors
        # equivalent to the properties .freqs and .eigenvecs
        # but we need to make sure we grab _freqs (or _reduced_freqs)
        # since Euphonic no longer attempts to handle varying units
        # internally
        self.data.calculate_fine_phonons(grid_q, **cfp_dict)
        freq = self.data._freqs    # (n_pt, n_br) # can this be replaced by _reduced_freqs?
        vecs = self.data.eigenvecs # (n_pt, n_br, n_io, 3) # can this be replaced by _reduced_eigenvecs?
        vecs = degenerate_check(grid_q, freq, vecs)
        self.sort_branches(sort, freq, vecs)
        self.parallel = parallel

    def _fill_grid(self, freq, vecs):
        n_pt = self.grid.invA.shape[0]
        n_br = self.data.n_branches
        n_io = self.data.n_ions
        if freq.shape == (n_pt, n_br):
            freq = freq.reshape(n_pt, n_br, 1)
        if freq.shape != (n_pt, n_br, 1):
            raise Exception('freqiencies have wrong shape')
        if vecs.shape != (n_pt, n_br, n_io, 3):
            raise Exception('eigenvectors have wrong shape')
        # We must provide extra content information to enable efficient
        # interpolation and possible rotations. All arrays provided to fill
        # must be 3+ dimensional -- the first dimension is over the points of
        # the grid, the second dimension is over the phonon branches, and any
        # additional dimensions are collapsed from highest to lowest dimension
        # (indexed as row-ordered linear indexing)
        # So, the vectors in (n_pt, n_br, n_io, 3) end up as
        # (n_pt, n_br, 3*n_io) with each three-vector's elements contiguous in
        # memory. The extra information provided details *what* the last
        # collapsed dimension contains, as a tuple, list, or 1-D array
        #   ( the number of scalar elements,
        #     the total number of vectors elements [must be divisible by 3],
        #     the total number of matrix elements [must be divisible by 9],
        #     how the vectors/matrices behave under rotation [0, 1, or 2] )
        # The last dimension elements must be ordered as scalars, then vectors
        # then matrices. The three rotation behaviour values are 0 ≡ Realspace
        # vectors, 1 ≡ Reciprocal space vectors, 2 ≡ (Realspace) Axial vectors.
        # Any missing entries from the tuple/list/array are treated as zeros.
        #
        # For the case of phonon eigenvalues and eigenvectors: there is one
        # eigenvalue per branch (scalars do not change under rotation), and
        # each branch eigenvector is comprised of n_ions displacement 3-vectors
        # which transform via the phonon Γ function,
        # so [1,0,0,0] ≡ (1,) and [0,n_ions*3,0,3] ≡ (0,3*n_io,0,3)
        self.grid.fill(freq, (1,), self.brspgl.orthogonal_to_conventional_eigenvectors(vecs), (0,3*n_io,0,3))


    def sort_branches(self, sort, frqs, vecs, energy_weight=1.0, eigenvector_weight=1.0, weight_function=0):
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
        self._fill_grid(frqs, vecs)
        if sort:
            # The input to sort_perm indicates what weight should be given to
            # each part of the resultant cost matrix. In this case, each phonon
            # branch consists of one energy, n_ions three-vectors, and no matrix;
            # perm = self.grid.centre_sort_perm(
            perm = self.grid.multi_sort_perm(
                scalar_cost_weight=energy_weight,
                vector_cost_weight=eigenvector_weight,
                matrix_cost_weight=0,
                vector_weight_function=weight_function)
            # use the permutations to sort the eigenvalues and vectors
            frqs = np.array([x[y] for (x, y) in zip(frqs, perm)])
            vecs = np.array([x[y,:,:] for (x,y) in zip(vecs, perm)])
            # # The gridded frequencies are (n_pt, n_br, 1)
            # frqs = np.array([x[y,:] for (x, y) in zip(self.grid.values, perm)])
            # # The gridded eigenvectors are (n_pt, n_br, n_io, 3)
            # vecs = np.array([x[y,:,:] for (x, y) in zip(self.grid.vectors, perm)])
            self._fill_grid(frqs, vecs)
        return frqs, vecs

    # pylint: disable=c0103,w0613,no-member
    def __define_grid_or_mesh(self, mesh=False, trellis=False, nest=False, **kwds):
        brillouin_zone = self.brspgl.get_conventional_BrillouinZone()
        if mesh:
            self.__make_mesh(brillouin_zone, **kwds)
        elif trellis:
            self.__make_trellis(brillouin_zone, **kwds)
        elif nest:
            self.__make_nest(brillouin_zone, **kwds)
        else:
            self.__make_grid(brillouin_zone, **kwds)
        # We need to make sure that we pass gridded Q points in the primitive
        # lattice, since that is what Euphonic expects:
        return self.brspgl.conventional_to_input_Q(self.grid.rlu)

    def __make_grid(self, bz, halfN=None, step=None, units='rlu', **kwds):
        if halfN is not None:
            # pybind is supposed to be able to convert to the C-expected type
            # automatically, but it can't distinguish between the call requiring
            # 'size_t' and 'double' so we need to help it out.
            if isinstance(halfN, (tuple, list)):
                halfN = np.array(halfN)
            if not isinstance(halfN, np.ndarray):
                raise Exception("'halfN' must be a tuple, list, or ndarray")
            halfN = halfN.astype('uint64')
            self.grid = br.BZGridQdc(bz, halfN)
        elif step is not None:
            if isinstance(step, ureg.Quantity):
                step = (step.to('angstrom')).magnitude
                isrlu = False
            else:
                isrlu = units.lower() == 'rlu'
            self.grid = br.BZGridQdc(bz, step, isrlu)
        else:
            raise Exception("keyword 'halfN' or 'step' required")

    def __make_mesh(self, bz, max_size=-1, max_points=-1, num_levels=3, **kwds):
        self.grid = br.BZMeshQdc(bz, max_size, num_levels, max_points)

    def __make_trellis(self, bz, max_volume=None, number_density=None, **kwds):
        if max_volume is not None:
            self.grid = br.BZTrellisQdc(bz, max_volume)
        elif number_density is not None:
            self.grid = br.BZTrellisQdc(bz, number_density)
        else:
            raise Exception("keyword 'max_volume' or 'number_density' required")

    def __make_nest(self, bz, max_branchings=5, max_volume=None, number_density=None, **kwds):
        if max_volume is not None:
            self.grid = br.BZNestQdc(bz, max_volume, max_branchings)
        elif number_density is not None:
            self.grid = br.BZNestQdc(bz, number_density, max_branchings)
        else:
            raise Exception("keyword 'max_volume' or 'number_density' required")

    def s_q(self, q_hkl, interpolate=True,  T=5.0, scale=1.0, dw_data=None, **kwargs):
        """Calculate Sᵢ(Q) where Q = (q_h,q_k,q_l)."""
        self.w_q(q_hkl, interpolate=interpolate, **kwargs)
        # Finally calculate Sᵢ(Q)
        if interpolate:
            sf = self._calculate_structure_factor(T=T, scale=scale, dw_data=dw_data)
        else:
            # using InterpolationData.calculate_structure_factor
            # which only allows a limited number of keyword arguments
            sf_keywords = ('calc_bose',)
            sf_dict = {k: kwargs[k] for k in sf_keywords if k in kwargs}
            sf = self.data.calculate_structure_factor(self.scattering_lengths, **sf_dict)
        return sf

    def dw(self, q_hkl, T=0):
        """Calculates the Debye-Waller factor using the Brillouin zone grid."""
        DWfactor = self.grid.debye_waller(q_hkl,
                                          self.data.ion_mass.to('meV*s**2/angstrom**2').magnitude,
                                          T)
        return DWfactor

    def _calculate_structure_factor(self, T=5.0, scale=1.0, dw_data=None):
            """
            Calculate the one phonon inelastic scattering at each q-point
            See M. Dove Structure and Dynamics Pg. 226

            Parameters
            ----------
            scattering_lengths : dictionary
                Dictionary of spin and isotope averaged coherent scattering legnths
                for each element in the structure in fm e.g.
                {'O': 5.803, 'Zn': 5.680}
            T : float, optional, default 5.0
                The temperature in Kelvin to use when calculating the Bose and
                Debye-Waller factors
            scale : float, optional, default 1.0
                Apply a multiplicative factor to the final structure factor.
            calc_bose : boolean, optional, default True
                Whether to calculate and apply the Bose factor
            dw_data : InterpolationData or PhononData object
                A PhononData or InterpolationData object with
                frequencies/eigenvectors calculated on a q-grid over which the
                Debye-Waller factor will be calculated

            Returns
            -------
            sf : (n_qpts, n_branches) float ndarray
                The structure factor for each q-point and phonon branch
            """
            sl = [self.scattering_lengths[x] for x in self.data.ion_type]


            freqs = self.data._freqs # abuse Euphonic to get the right units
            ion_mass = self.data._ion_mass # ditto
            sl = (sl*ureg('fm').to('bohr')).magnitude

            # Calculate normalisation factor
            norm_factor = sl/np.sqrt(ion_mass)

            # Calculate the exponential factor for all ions and q-points
            # ion_r in fractional coords, so Qdotr = 2pi*qh*rx + 2pi*qk*ry...
            exp_factor = np.exp(1J*2*np.pi*np.einsum('ij,kj->ik',

                                                       self.data.qpts, self.data.ion_r))

            # brille prefers eigenvectors in units of the lattice, so we don't need
            # to modify the units of Q

            # Calculate dot product of Q and eigenvectors for all branches, ions
            # and q-points
            eigenv_dot_q = np.einsum('ijkl,il->ijk', np.conj(self.data.eigenvecs), self.data.qpts)

            # Calculate Debye-Waller factors
            if dw_data:
                if dw_data.n_ions != self.data.n_ions:

                    raise Exception((
                        'The Data object used as dw_data is not compatible with the'
                        ' object that calculate_structure_factor has been called on'
                        ' (they have a different number of ions). Is dw_data '
                        'correct?'))
                dw = dw_data._dw_coeff(T)
                dw_factor = np.exp(-np.einsum('jkl,ik,il->ij', dw, self.data.qpts, self.data.qpts)/2)
                exp_factor *= dw_factor

            # Multiply Q.eigenvector, exp factor and normalisation factor
            term = np.einsum('ijk,ik,k->ij', eigenv_dot_q, exp_factor, norm_factor)

            # Take mod squared and divide by frequency to get intensity
            sf = np.absolute(term*np.conj(term))/np.absolute(freqs)

            sf = np.real(sf*scale)

            return sf

    def w_q(self, q_pt, moveinto=True, interpolate=True, threads=-1, **kwargs):
        """Calculate ωᵢ(Q) where Q = (q_h,q_k,q_l)."""
        if interpolate:
            # Interpolate the previously-stored eigen values/vectors for each Q
            # each grid point has a (n_br, 1) values array and a (n_br, n_io, 3)
            # eigenvectors array and interpolate_at returns a tuple with the
            # first entry a (n_pt, n_br, 1) values array and the second a
            # (n_pt, n_br, n_io, 3) eigenvectors array
            frqs, vecs = self.grid.ir_interpolate_at(q_pt, self.parallel, threads, not moveinto)
            n_pt = q_pt.shape[0]
            # Store all information in the Euphonic Data object
            self.data.n_qpts = n_pt
            self.data.qpts = q_pt # no lattice conversion since we use internal Structure Factor calculation
            self.data._reduced_freqs = np.squeeze(frqs) # Euphonic expects this to be (n_pt, n_br)
            self.data._reduced_eigenvecs = vecs # already (n_pt, n_br, n_io, 3)
            self.data._qpts_i = np.arange(n_pt, dtype=np.int32)
        else:
            cfp_keywords = ('asr', 'precondition', 'dipole', 'eta_scale'
                            'splitting', 'reduce_qpts', 'use_c', 'n_threads')
            cfp_dict = {k: kwargs[k] for k in cfp_keywords if k in kwargs}
            self.data.calculate_fine_phonons(self.brspgl.conventional_to_input_Q(q_pt), **cfp_dict)
        return self.data.freqs # now a property, handles unit conversion

    def __call__(self, *args, **kwargs):
        """Calculate and return Sᵢ(Q) and ωᵢ(Q) or S(Q,ω) depending on input.

        If one positional argument is provided it is assumed to be Q in which
        case both the intensity, Sᵢ(Q), and eigen-energy, ωᵢ(Q), for all phonon
        branches are returned as a tuple.
        If two positional arguments are provided they are assumed to be Q and ω
        in which case Sᵢ(Q) and ωᵢ(Q) are used in conjunction with keyword
        arguments 'resfun' and 'param' to calculate a convolved S(Q,ω).
        """
        if len(args) < 1:
            raise RuntimeError('At least one argument, Q, is required')
        elif len(args) is 1:
            s_i_of_q = self.s_q(*args, **kwargs)
            return (s_i_of_q, self.data.freqs)
        elif len(args) is 2:
            return self.s_qw(*args, kwargs) # keep kwargs as a dictionary
        else:
            raise RuntimeError('Only one or two arguments expected, (Q,) or (Q,ω), expected')

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

        Additional keys in p_dict are allowed and are passed on to Euphonic
        as keyword arguments to the calculate_structure_factor method.

        """
        res_par_tem = ('delta',)
        if 'resfun' in p_dict and 'param' in p_dict:
            res_par_tem = (p_dict['resfun'].replace(' ', '').lower(),
                           p_dict['param'])
        n_pt = energy.size
        n_br = self.data.n_branches
        # Check if we might perform the Bose factor correction twice:
        # Replicate Euphonic's behaviour of calc_bose=True by default,
        # and T=5 by default.
        if res_par_tem[0] in ('s', 'sho', 'simpleharmonicoscillator'):
            # pull out T, or 5 if it's not present
            temp_k = p_dict.get('T', 5)
            # If calc_bose is present
            if 'calc_bose' in p_dict:
                # keep T if it's True, discard T if it's False
                temp_k = temp_k if p_dict['calc_bose'] else 0
            # Prevent Euphonic from performing the Bose correction twice
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
    # if the brille object is holding complex values (it is) then its returned
    # interpolated values are all complex too, even, e.g., energies which are
    # purely real with identically zero imaginary components.
    # The division of two purely-real complex numbers in Python will annoyingly
    # raise a warning about discarding the imaginary part. So preempt it here.
    if not np.isclose(np.sum(np.abs(np.imag(part1))+np.abs(np.imag(part2))), 0.):
        raise RuntimeError('Unexpected imaginary components.')
    y_0[flag] = np.real(part1)/np.real(part2)
    return y_0.reshape(outshape)

def half_cpu_count():
    import os
    count = os.cpu_count()
    if 'sched_get_affinity' in dir(os):
        count = len(os.sched_getaffinity(0))
    return count//2
