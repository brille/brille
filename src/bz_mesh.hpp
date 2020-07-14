/* Copyright 2019 Greg Tucker
//
// This file is part of brille.
//
// brille is free software: you can redistribute it and/or modify it under the
// terms of the GNU Affero General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// brille is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with brille. If not, see <https://www.gnu.org/licenses/>.            */

typedef long slong;

#ifndef _BZ_MESH_
#define _BZ_MESH_

#include "bz.hpp"
#include "mesh.hpp"

template<class T, class S> class BrillouinZoneMesh3: public Mesh3<T,S>{
protected:
  BrillouinZone brillouinzone;
public:
  /*! Construct using a maximum tetrahedron volume -- makes a tetrahedron mesh
      instead of a orthogonal grid.
      @param bz The BrillouinZone object
      @param vol The maximum tetrahedron volume
      @param isrlu A flag to indicate if vol is in units of rlu (isrlu=1) or inverse angstrom (isrlu=0)
      @note If vol is in relative lattice units an absolute volume will be
            determined using the unit cell volume of the underlying lattice.
  */
  // BrillouinZoneMesh3(const BrillouinZone& bz, const double max_size_invA=-1., const double min_angle=20.0, const double max_angle=-1.0, const double max_ratio=-1., const int max_points=-1):
  //   Mesh3<T>(bz.get_ir_vertices().get_xyz(), bz.get_ir_vertices_per_face(), max_size_invA, min_angle, max_angle, max_ratio, max_points),
  //   brillouinzone(bz) {}
  // BrillouinZoneMesh3(const BrillouinZone& bz) Mesh3<T>(bz.get_ir_vertices().get_xyz(), bz.get_ir_vertices_per_face());
  template<typename... A>
  BrillouinZoneMesh3(const BrillouinZone& bz, A... args):
    Mesh3<T,S>(bz.get_ir_vertices().get_xyz(), bz.get_ir_vertices_per_face(), args...),
    brillouinzone(bz) {}
  // get the BrillouinZone object
  BrillouinZone get_brillouinzone(void) const {return this->brillouinzone;}
  // get the mesh vertices in relative lattice units
  ArrayVector<double> get_mesh_hkl(void) const {
    ArrayVector<double> xyz = this->get_mesh_xyz();
    double toxyz[9], fromxyz[9];
    const BrillouinZone bz = this->get_brillouinzone();
    bz.get_lattice().get_xyz_transform(toxyz);
    if (!matrix_inverse(fromxyz,toxyz)) throw std::runtime_error("transform matrix toxyz has zero determinant");
    ArrayVector<double> hkl(3,xyz.size());
    for (size_t i=0; i<xyz.size(); i++) multiply_matrix_vector<double,double,double,3>(hkl.data(i), fromxyz, xyz.data(i));
    return hkl;
  }

  template<typename R>
  std::tuple<ArrayVector<T>,ArrayVector<S>>
  ir_interpolate_at(const LQVec<R>& x, const int nthreads, const bool no_move=false) const {
    LQVec<R> ir_q(x.get_lattice(), x.size());
    LQVec<int> tau(x.get_lattice(), x.size());
    std::vector<size_t> rot(x.size(),0u), invrot(x.size(),0u);
    if (no_move){
      ir_q = x;
    } else if (!brillouinzone.ir_moveinto(x, ir_q, tau, rot, invrot, nthreads)){
      std::string msg;
      msg = "Moving all points into the irreducible Brillouin zone failed.";
      throw std::runtime_error(msg);
    }
    // perform the interpolation within the irreducible Brillouin zone
    ArrayVector<T> vals;
    ArrayVector<S> vecs;
    std::tie(vals,vecs) = (nthreads > 1)
        ? this->Mesh3<T,S>::parallel_interpolate_at(ir_q.get_xyz(), nthreads)
        : this->Mesh3<T,S>::interpolate_at(ir_q.get_xyz());
    // we always need the pointgroup operations to 'rotate'
    PointSymmetry psym = brillouinzone.get_pointgroup_symmetry();
    // and might need the Phonon Gamma table
    GammaTable pgt{GammaTable()};
    if (RotatesLike::Gamma == this->data().vectors().rotateslike()){
      pgt.construct(brillouinzone.get_lattice().star(), brillouinzone.add_time_reversal());
    }
    // actually perform the rotation to Q
    this->data().values() .rotate_in_place(vals, ir_q, pgt, psym, rot, invrot, nthreads);
    this->data().vectors().rotate_in_place(vecs, ir_q, pgt, psym, rot, invrot, nthreads);
    // we're done so bundle the output
    return std::make_tuple(vals, vecs);
  }
};

#endif // _BZ_MESH_
