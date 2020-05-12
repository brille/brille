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

#ifndef _BZ_NEST_
#define _BZ_NEST_

#include <tuple>
#include "bz.hpp"
#include "nest.hpp"

template<class T, class S> class BrillouinZoneNest3: public Nest<T,S>{
  BrillouinZone brillouinzone;
public:
  template<typename... A>
  BrillouinZoneNest3(const BrillouinZone& bz, A... args):
    Nest<T,S>(bz.get_ir_polyhedron(), args...),
    brillouinzone(bz) {}
  //! get the BrillouinZone object
  BrillouinZone get_brillouinzone(void) const {return this->brillouinzone;}
  //! get the vertices of the leaf vertices in inverse Angstrom
  ArrayVector<double> get_xyz(void) const {return this->vertices();}
  //! get the vertices of all vertices in absolute units
  const ArrayVector<double>& get_all_xyz(void) const {return this->all_vertices(); }
  //! get the vertices of the leaf vertices in relative lattice units
  ArrayVector<double> get_hkl(void) const { return xyz_to_hkl(brillouinzone.get_lattice(),this->vertices());}
  //! get the vertices of the inner (cubic) nodes in relative lattice units
  ArrayVector<double> get_all_hkl(void) const {return xyz_to_hkl(brillouinzone.get_lattice(),this->all_vertices()); }
  // //! get the indices forming the faces of the tetrahedra
  // std::vector<std::array<size_t,4>> get_vertices_per_tetrahedron(void) const {return this->tetrahedra();}

  template<typename R>
  std::tuple<ArrayVector<T>,ArrayVector<S>>
  ir_interpolate_at(const LQVec<R>& x, const int nth, const bool no_move=false) const{
    LQVec<R> ir_q(x.get_lattice(), x.size());
    LQVec<int> tau(x.get_lattice(), x.size());
    std::vector<size_t> rot(x.size(),0u), invrot(x.size(),0u);
    if (no_move){
      ir_q = x;
    } else if (!brillouinzone.ir_moveinto(x, ir_q, tau, rot, invrot, nth)){
      std::string msg;
      msg = "Moving all points into the irreducible Brillouin zone failed.";
      throw std::runtime_error(msg);
    }
    // perform the interpolation within the irreducible Brillouin zone
    ArrayVector<T> vals;
    ArrayVector<S> vecs;
    std::tie(vals,vecs) = (nth > 1)
        ? this->Nest<T,S>::interpolate_at(ir_q.get_xyz(), nth)
        : this->Nest<T,S>::interpolate_at(ir_q.get_xyz());
    // we always need the pointgroup operations to 'rotate'
    PointSymmetry psym = brillouinzone.get_pointgroup_symmetry();
    // and might need the Phonon Gamma table
    GammaTable pgt{GammaTable()};
    if (RotatesLike::Gamma == this->data().vectors().rotateslike()){
      pgt.construct(brillouinzone.get_lattice().star(), brillouinzone.add_time_reversal());
    }
    // actually perform the rotation to Q
    this->data().values().rotate_in_place(vals, ir_q, pgt, psym, rot, invrot, nth);
    this->data().vectors().rotate_in_place(vecs, ir_q, pgt, psym, rot, invrot, nth);
    // we're done so bundle the output
    return std::make_tuple(vals, vecs);
  }
};

#endif
