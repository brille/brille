import symbz
import numpy as np
import spglib

from simphony.data.bands import BandsData
from simphony.data.phonon import PhononData
from simphony.data.interpolation import InterpolationData
from simphony.calculate.scattering import structure_factor
from simphony import ureg # avoid creating a second Pint UnitRegistry

from scipy.special import wofz
from scipy.stats import norm, cauchy
def delta(x,xi,yi):
  y = np.zeros_like(yi)
  y[x==xi] = yi[x==xi]
  return y.sum(1)

def gaussian(x,xi,yi,fwhm):
  if not np.isscalar(fwhm):
    fwhm=fwhm[0]
  sigma = fwhm/np.sqrt(np.log(256))
  z = (x-xi)/sigma
  y = norm.pdf(z) * yi
  return y.sum(1)

def lorentzian(x,xi,yi,fwhm):
  if not np.isscalar(fwhm):
    fwhm=fwhm[0]
  gamma = fwhm/2
  z = (x-xi)/gamma
  y = cauchy.pdf(z) * yi
  return y.sum(1)

def voigt(x,xi,yi,params):
  if np.isscalar(params):
    g = params
    l = 0
  else:
    g = params[0]
    l = params[1]
  if l is 0:
    return gaussian(x,xi,yi,g)
  elif g is 0:
    return lorentzian(x,xi,yi,l)

  A = np.sqrt(np.log(2)/np.pi)
  gamma = g/2
  rz = np.sqrt(np.log(2))*(x-xi)/gamma
  iz = np.sqrt(np.log(2))*np.abs(l/g)
  y = A*np.real(wofz(rz + 1j*iz))/gamma
  return y.sum(1)

def sho(x,xi,yi,fwhm,T):
  if not np.isscalar(fwhm):
    fwhm=fwhm[0]
  if np.abs(T) > 0:
    bose = x / (1-np.exp(-11.602*x/T))
  else:
    bose = 1.0
  x2 = np.repeat( x**2, xi.shape[1], 1)
  ok = (xi>0)*np.isfinite(xi)
  y = np.zeros_like(yi)
  y[ok] = (4/np.pi)*fwhm*xi[ok]*yi[ok]/( (x2[ok]-xi[ok]**2)**2 + 4*fwhm**2*x2[ok] )
  y = y*bose
  return y.sum(1)

class SymSim(object):
  def __init__(self, SPData, scattering_lengths=None, N=None, d=None, units='rlu', parallel=False, **kwds):
    dt = type(SPData)
    if dt is not BandsData or dt is not PhononData or dt is not InterpolationData:
      print( f"Unexpected data type {dt}, expect failures.")
    if scattering_lengths is None or type(scattering_lengths) is not dict:
      scattering_lengths = {k:1 for k in np.unique(SPData.ion_type)}
    self.data = SPData
    self.scattering_lengths = scattering_lengths

    lattice_vectors = (self.data.cell_vec.to('angstrom')).magnitude
    ion_positions = self.data.ion_r
    _, ion_indexes = np.unique(self.data.ion_type, return_inverse=True)
    cell = (lattice_vectors, ion_positions, ion_indexes)
    symmetry_data = spglib.get_symmetry_dataset(cell)

    dlat = symbz.Direct(lattice_vectors, symmetry_dat['hall_number'])
    rlat = dlat.star()
    brillouinZone = symbz.BrillouinZone(rlat)
    if N is not None:
      self.grid = symbz.BZGridQcomplex(brillouinZone, N)
    elif d is not None:
      if isinstance(d,ureg.Quantity):
        d = (d.to('angstrom')).magnitude
        isrlu = False
      else:
        isrlu = units.lower()=='rlu'
      self.grid = symbz.BZGridQcomplex(brillouinZone, d, isrlu)

    # Calculate ωᵢ(Q) and ⃗ϵᵢⱼ(Q), and fill the BZGrid:
    Q = self.grid.mapped_rlu
    freq,vecs = self.data.calculate_fine_phonons(Q,**kwds)
    nQ = Q.shape[0]
    nB = self.data.n_branches
    nI = self.data.n_ions
    fv = np.concatenate(( (freq.magnitude).reshape((nQ,nB,1)), vecs.reshape(nQ,nB,3*nI) ), axis=2)
    self.grid.fill(fv)

    self.parallel=parallel

  def sq(self, qh, qk, ql, **kwargs):
    # Construct the (nQ,3) Q array
    Q = np.stack( (qh.flatten(), qk.flatten(), ql.flatten()) ).transpose()
    # Interpolate the previously-stored eigen values and vectors for each Q
    fv = self.grid.interpolate_at(Q,True,self.parallel)
    # Separate them
    f,v = np.split(fv, np.array([1]), axis=2)
    # And reshape to what SimPhony expects
    nQ = Q.shape[0]
    nB = self.data.n_branches
    nI = self.data.n_ions
    f = f.reshape((nQ,nB))
    v = v.reshape((nQ,nB,nI,3))
    # Store all information in the SimPhony Data object
    self.data.n_qpts = nQ
    self.data.qpts = Q
    self.data.freqs = f*self.data.freqs.units
    self.data.eigenvecs = v
    # Finally calculate Sᵢ(Q) using SymPhony.calculate.scattering.structure_factor
    return structure_factor(self.data,self.scattering_lengths,**kwargs)

  def sqw(self, qh,qk,ql,en, p):
    nQ = en.size
    nB = self.data.n_branches
    # Calculate S(Q) after interpolating ωᵢ(Q) and ⃗ϵᵢⱼ(Q)
    sf = self.sq(qh,qk,ql, **p)
    # The resulting array should be (nQ,nB)
    assert sf.shape[0] is nQ and sf.shape[1] is nB, f"Expected S(Q) shape ({nQ},{nB}) but got {sf.shape}."

    omega = (self.data.freqs.to('millielectron_volt')).magnitude
    shapein = en.shape
    en = en.flatten()[:,None]
    # ωᵢ(Q) is (nQ,nB)
    # Sᵢ(Q) is (nQ,nB)
    # en    is (nQ,1)

    # The dict p should contain parameters *only* for non-resolution broadening effects
    if 'resfun' in p and 'param' in p:
      resfun = p['resfun'].lower()
      params = p['param']
      if resfun is 'g' or resfun is 'gauss' or resfun is 'gaussian':
        S = gaussian(en,omega,sf,params)
      elif resfun is 'l' or resfun is 'lor' or resfun is 'lorentz' or resfun is 'lorentzian':
        S = lorentzian(en,omega,sf,params)
      elif resfun is 'v' or resfun is 'voi' or resfun is 'voigt':
        S = voigt(en,omega,sf,params)
      elif resfun is 's' or resfun is 'sho' or resfun is 'simple harmonic oscillator':
        if 'T' in p:
          T = p['T']
        else:
          T = 0;
        S = sho(en,omega,sf,params,T)
      else:
        print(f'Unknown function {resfun}')
    else:
      S = delta(en,omega,sf)

    if S.shape != shapein:
      S = S.reshape(shapein)

    return S
