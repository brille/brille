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

template<class T> class BrillouinZoneMesh3: public Mesh3<T>{
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
    Mesh3<T>(bz.get_ir_vertices().get_xyz(), bz.get_ir_vertices_per_face(), args...),
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
  template<typename R> ArrayVector<T> interpolate_at(const LQVec<R>& x, const int nthreads, const bool no_move=false) const{
    LQVec<R> ir_q(x.get_lattice(), x.size());
    LQVec<int> tau(x.get_lattice(), x.size());
    std::vector<std::array<int,9>> rots(x.size());
    BrillouinZone bz = this->get_brillouinzone();

    std::string msg;
    if (no_move){
      // Special mode for testing where no specified points are moved
      // IT IS IMPERITIVE THAT THE PROVIDED POINTS ARE *INSIDE* THE IRREDUCIBLE
      // POLYHEDRON otherwise the interpolation will fail or give garbage back.
      ir_q = x;
      for (size_t i=0; i<x.size(); ++i) rots[i] = {1,0,0, 0,1,0, 0,0,1};
    } else if (!bz.ir_moveinto(x, ir_q, tau, rots, nthreads)){
      msg = "Moving all points into the irreducible Brillouin zone failed.";
      throw std::runtime_error(msg);
    }
    ArrayVector<T> ir_result = (nthreads < 2)
      ? this->Mesh3<T>::interpolate_at(ir_q.get_xyz())
      : this->Mesh3<T>::parallel_interpolate_at(ir_q.get_xyz(), nthreads);

    if (nthreads < 2)
      this->data().rotate_in_place(ir_result, rots);
    else
      this->data().rotate_in_place(ir_result, rots, nthreads);
      
    return ir_result;
  }
};

#endif // _BZ_MESH_
