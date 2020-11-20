/* This file is part of brille.

Copyright © 2019,2020 Greg Tucker <greg.tucker@stfc.ac.uk>

brille is free software: you can redistribute it and/or modify it under the
terms of the GNU Affero General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option)
any later version.

brille is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with brille. If not, see <https://www.gnu.org/licenses/>.            */

#ifndef BRILLE_BZ_MESH_
#define BRILLE_BZ_MESH_
#include "bz.hpp"
#include "mesh.hpp"
namespace brille {

template<class T, class S>
class BrillouinZoneMesh3: public Mesh3<T,S>{
  using SuperClass = Mesh3<T,S>;
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
    SuperClass(bz.get_ir_vertices().get_xyz(), bz.get_ir_vertices_per_face(), args...),
    brillouinzone(bz) {}
  // get the BrillouinZone object
  BrillouinZone get_brillouinzone(void) const {return this->brillouinzone;}
  // get the mesh vertices in relative lattice units
  bArray<double> get_mesh_hkl(void) const {
    auto xyz = this->get_mesh_xyz();
    double toxyz[9], fromxyz[9];
    const BrillouinZone bz = this->get_brillouinzone();
    bz.get_lattice().get_xyz_transform(toxyz);
    if (!brille::utils::matrix_inverse(fromxyz,toxyz)) throw std::runtime_error("transform matrix toxyz has zero determinant");
    auto shape = xyz.shape();
    bArray<double> hkl(shape);
    std::vector<double> tmp(3);
    for (size_t i=0; i<shape[0]; ++i){
      auto vxyz = xyz.view(i).to_std();
      brille::utils::multiply_matrix_vector<double,double,double,3>(tmp.data(), fromxyz, vxyz.data());
      hkl.set(i, tmp);
    }
    return hkl;
  }

  template<class R>
  std::tuple<brille::Array<T>,brille::Array<S>>
  ir_interpolate_at(const LQVec<R>& x, const int nthreads, const bool no_move=false) const {
    LQVec<R> ir_q(x.get_lattice(), x.size(0));
    LQVec<int> tau(x.get_lattice(), x.size(0));
    std::vector<size_t> rot(x.size(0),0u), invrot(x.size(0),0u);
    if (no_move){
      ir_q = x;
    } else if (!brillouinzone.ir_moveinto(x, ir_q, tau, rot, invrot, nthreads)){
      std::string msg;
      msg = "Moving all points into the irreducible Brillouin zone failed.";
      throw std::runtime_error(msg);
    }
    // perform the interpolation within the irreducible Brillouin zone
    auto [vals, vecs] = (nthreads > 1)
        ? this->SuperClass::parallel_interpolate_at(ir_q.get_xyz(), nthreads)
        : this->SuperClass::interpolate_at(ir_q.get_xyz());
    // we always need the pointgroup operations to 'rotate'
    PointSymmetry psym = brillouinzone.get_pointgroup_symmetry();
    // and might need the Phonon Gamma table
    GammaTable pgt{GammaTable()};
    if (RotatesLike::Gamma == this->data().vectors().rotateslike()){
      pgt.construct(brillouinzone.get_lattice().star(), brillouinzone.add_time_reversal());
    }
    brille::Array2<T> vals2(vals);
    brille::Array2<S> vecs2(vecs);
    // actually perform the rotation to Q
    this->data().values() .rotate_in_place(vals2, ir_q, pgt, psym, rot, invrot, nthreads);
    this->data().vectors().rotate_in_place(vecs2, ir_q, pgt, psym, rot, invrot, nthreads);
    // we're done so bundle the output
    return std::make_tuple(vals, vecs);
  }
};

} // end namespace brille
#endif // _BZ_MESH_
