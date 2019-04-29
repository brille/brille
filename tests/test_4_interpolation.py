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

from importlib import util

if util.find_spec('symbz') is not None and util.find_spec('symbz._symbz') is not None:
    import symbz as s
elif util.find_spec('_symbz') is not None:
    import _symbz as s
else:
    raise Exception("symbz module not found!")

def sqwfunc_ones(Q):
  return 1.0+0*(Q[:,0]+Q[:,1]+Q[:,2])
def sqwfunc_x(Q):
  return Q[:,0]
def sqwfunc_y(Q):
  return Q[:,1]
def sqwfunc_z(Q):
  return Q[:,2]
def sqwfunc_xy(Q):
  return Q[:,0]+Q[:,1]
def sqwfunc_xyz(Q):
  return Q[:,0]+Q[:,1]+Q[:,2]

def vecfun_ident(Q):
  return Q;
def vecfun_rotz(Q,θ=0):
  x=Q[:,[0]]
  y=Q[:,[1]]
  z=Q[:,[2]]
  c = np.cos(θ)
  s = np.sin(θ)
  return np.concatenate( ( x*c-y*s, x*s+y*c, z ) ,axis=1)
def vecfun_rotx(Q,θ=0):
  x=Q[:,[0]]
  y=Q[:,[1]]
  z=Q[:,[2]]
  c = np.cos(θ)
  s = np.sin(θ)
  return np.concatenate( ( x, y*c-z*s, y*s+z*c ) ,axis=1)

def matfun_ident(Q):
  sh = Q.shape
  z = np.ndarray( (sh[0],sh[1],sh[1]), dtype=Q.dtype)
  for i in range(sh[0]):
    z[i,:,:] = np.diag(Q[i,:])
  return z

def setup_grid():
  rlat = s.Reciprocal( (1,1,1), np.array([1,1,1])*np.pi/2 )
  bz = s.BrillouinZone(rlat)
  bzg = s.BZGrid(bz, N=(2,2,2))
  return bzg

def define_Q_points():
  return np.array( [[0,0,0],[0.1,0,0],[0,0.1,0],[0,0,0.1]],dtype='double')

class Interpolate (unittest.TestCase):
  def test_a_norm(self):
    bzg = setup_grid()
    Qi = define_Q_points()
    bzg.fill( sqwfunc_ones( bzg.mapped_rlu ) )
    interpolated_ones = bzg.interpolate_at(Qi)
    self.assertTrue( (interpolated_ones==1).all() )
  def test_b_x(self):
    bzg = setup_grid()
    Qi = define_Q_points()
    bzg.fill( sqwfunc_x(bzg.mapped_rlu) )
    intres = bzg.interpolate_at(Qi)
    self.assertTrue( np.isclose(intres,Qi[:,0]).all() )
  def test_b_y(self):
    bzg = setup_grid()
    Qi = define_Q_points()
    bzg.fill( sqwfunc_y(bzg.mapped_rlu) )
    intres = bzg.interpolate_at(Qi)
    self.assertTrue( np.isclose(intres,Qi[:,1]).all() )
  def test_b_z(self):
    bzg = setup_grid()
    Qi = define_Q_points()
    bzg.fill( sqwfunc_z(bzg.mapped_rlu) )
    intres = bzg.interpolate_at(Qi)
    self.assertTrue( np.isclose(intres,Qi[:,2]).all() )
  def test_c_xy(self):
    bzg = setup_grid()
    Qi = define_Q_points()
    bzg.fill( sqwfunc_xy(bzg.mapped_rlu) )
    intres = bzg.interpolate_at(Qi)
    self.assertTrue( np.isclose(intres,Qi[:,0]+Qi[:,1]).all() )
  def test_c_xyz(self):
    bzg = setup_grid()
    Qi = define_Q_points()
    bzg.fill( sqwfunc_xyz(bzg.mapped_rlu) )
    intres = bzg.interpolate_at(Qi)
    self.assertTrue( np.isclose(intres,Qi[:,0]+Qi[:,1]+Qi[:,2]).all() )
  def test_d_vec_ident(self):
    bzg = setup_grid();
    Qi = define_Q_points()
    bzg.fill( vecfun_ident(bzg.mapped_rlu) )
    intres = bzg.interpolate_at(Qi)
    self.assertTrue(np.isclose(intres,Qi).all())
  def test_d_vec_rotx(self):
    bzg = setup_grid();
    Qi = define_Q_points()
    ang = np.pi/3;
    bzg.fill( vecfun_rotx(bzg.mapped_rlu,ang) )
    intres = bzg.interpolate_at(Qi)
    antres = vecfun_rotx(Qi,ang)
    self.assertTrue(np.isclose(intres,antres).all())
  def test_d_vec_rotz(self):
    bzg = setup_grid();
    Qi = define_Q_points()
    ang = 3*np.pi/5;
    bzg.fill( vecfun_rotz(bzg.mapped_rlu,ang) )
    intres = bzg.interpolate_at(Qi)
    antres = vecfun_rotz(Qi,ang)
    self.assertTrue(np.isclose(intres,antres).all())
  def test_e_mat_ident(self):
    bzg = setup_grid()
    Qi = define_Q_points()
    bzg.fill( matfun_ident(bzg.mapped_rlu) )
    intres = bzg.interpolate_at(Qi)
    antres = matfun_ident(Qi)
    self.assertTrue(np.isclose(intres,antres).all())

if __name__ == '__main__':
  unittest.main()
