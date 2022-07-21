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
/*! \file
    \author Greg Tucker
    \brief Defines a class to extend `Nest` with `BrillouinZone` information.
*/
#ifndef BRILLE_BZ_NEST_
#define BRILLE_BZ_NEST_
// #include <tuple>
#include "bz.hpp"
#include "nest.hpp"
namespace brille {
  /*! \brief A Nest in a BrillouinZone

  The first or irreducible Brillouin zone Polyhedron contained in a BrillouinZone
  object can be used to define the domain of a Nest triangulation.
  The symmetries of the Brillouin zone can then be used to interpolate at any
  point in reciprocal space by finding an equivalent point within the triangulated
  domain.
  */
template<class T, class S, class A>
class BrillouinZoneNest3: public Nest<T,S,A,Array2>{
  using class_t = BrillouinZoneNest3<T,S,A>;
  using base_t = Nest<T,S,A,Array2>;
  template<class V> using lv_t = lattice::LVec<V>;
  template<class V> using bv_t = Array2<V>;
  BrillouinZone bz_;
public:
  BrillouinZoneNest3(const base_t& pt, BrillouinZone bz): base_t(pt), bz_(std::move(bz)) {}
  BrillouinZoneNest3(base_t&& pt, BrillouinZone&& bz): base_t(std::move(pt)), bz_(std::move(bz)) {}
  /*! \brief Construct a `BrillouinZoneNest3` from a `BrillouinZone` and variable arguments

  All arguments beyond the `BrillouinZone` are passed to the `Nest` constructor.
  \param bz the `BrillouinZone` used to define the boundaries of the `Mesh3`
  \param args the construction arguments for `Nest`
  */
  template<typename... Args>
  BrillouinZoneNest3(const BrillouinZone& bz, Args... args):
    base_t({bz.get_ir_polyhedron().vertices().xyz(), bz.get_ir_polyhedron().faces()}, args...),
    bz_(bz) {}
  //! get the BrillouinZone object
  BrillouinZone get_brillouinzone(void) const {return this->bz_;}
  //! get the vertices of the leaf vertices in inverse Angstrom
  bv_t<A> get_xyz(void) const {return this->vertices();}
  //! get the vertices of all vertices in absolute units
  const bv_t<A>& get_all_xyz(void) const {return this->all_vertices(); }
  //! get the vertices of the leaf vertices in relative lattice units
  bv_t<A> get_hkl(void) const { return from_xyz_like(LengthUnit::inverse_angstrom, bz_.get_lattice(), this->vertices()).hkl();}
  //! get the vertices of the inner (cubic) nodes in relative lattice units
  bv_t<A> get_all_hkl(void) const {return from_xyz_like(LengthUnit::inverse_angstrom, bz_.get_lattice(), this->all_vertices()).hkl(); }
  // //! get the indices forming the faces of the tetrahedra
  // std::vector<std::array<size_t,4>> get_vertices_per_tetrahedron(void) const {return this->tetrahedra();}

  /*! \brief Interpolate at an equivalent irreducible reciprocal lattice point

  \param x        One or more points expressed in the same reciprocal lattice as
                  the stored `BrillouinZone`
  \param nth      the number of parallel OpenMP workers to utilize
  \param no_move  If all provided points are *already* within the irreducible
                  Brillouin zone this optional parameter can be used to skip a
                  call to `BrillouinZone::ir_moveinto`.
  \return a tuple of the interpolated eigenvalues and eigenvectors

  The interpolation is performed by `Nest::interpolate_at` and then
  corrected for the pointgroup operation by `Interpolator::rotate_in_place`.
  If the stored data has the same behaviour under application of the pointgroup
  operation as Phonon eigenvectors, then the appropriate `GammaTable` is
  constructed and used as well.

  \warning The last parameter should only be used with extreme caution as no
           check is performed to ensure that the points are actually in the
           irreducible Brillouin zone. If this condition is not true and the
           parameter is set to true, the subsequent interpolation call may raise
           an error or access unassigned memory and will produce garbage output.
  */
  template<bool NO_MOVE=false, class... Args>
  std::tuple<brille::Array<T>,brille::Array<S>>
  ir_interpolate_at(const lv_t<A>& x, Args... args) const {
    lv_t<A> ir_q(x.type(), x.lattice(), x.size(0));
    lv_t<int> tau(x.type(), x.lattice(), x.size(0));
    std::vector<size_t> rot(x.size(0),0u), invrot(x.size(0),0u);
    if constexpr (NO_MOVE){
      ir_q = x;
    } else if (!bz_.ir_moveinto(x, ir_q, tau, rot, invrot, args...)){
      std::string msg;
      msg = "Moving all points into the irreducible Brillouin zone failed.";
      throw std::runtime_error(msg);
    }
    // perform the interpolation within the irreducible Brillouin zone
    auto [vals, vecs] = this->base_t::interpolate_at(brille::get_xyz(ir_q), args...);
    // we always need the pointgroup operations to 'rotate'
    PointSymmetry psym = bz_.get_pointgroup_symmetry();
    // and might need the Phonon Gamma table
    auto cfg = this->approx_config();
    auto s_tol = cfg.template direct<double>();
    auto n_tol = cfg.digit();
    bool needed = RotatesLike::Gamma == this->data().vectors().rotateslike();
    auto pgt = GammaTable(needed, bz_.get_lattice(), bz_.add_time_reversal(), s_tol, n_tol);
    //
    brille::Array2<T> vals2(vals);
    brille::Array2<S> vecs2(vecs);
    // actually perform the rotation to Q
    this->data().values().rotate_in_place(vals2, ir_q, pgt, psym, rot, invrot, args...);
    this->data().vectors().rotate_in_place(vecs2, ir_q, pgt, psym, rot, invrot, args...);
    // we're done so bundle the output
    return std::make_tuple(vals, vecs);
  }

#ifdef USE_HIGHFIVE
    template<class HF>
    std::enable_if_t<std::is_base_of_v<HighFive::Object, HF>, bool>
    to_hdf(HF& obj, const std::string& entry) const{
        auto group = overwrite_group(obj, entry);
        bool ok{true};
        ok &= base_t::to_hdf(group, "nest");
        ok &= bz_.to_hdf(group, "bz_");
        return ok;
    }
    // Implementing this requires Nest3::to/from_hdf and therefore NestNode and NestLeaf to/from_hdf, which is problematic
//    // Input from HDF5 file/object
//    template<class HF>
//    static std::enable_if_t<std::is_base_of_v<HighFive::Object, HF>, BrillouinZoneNest3<T,S>>
//    from_hdf(HF& obj, const std::string& entry){
//        auto group = obj.getGroup(entry);
//        auto nest = super_t::from_hdf(group, "nest");
//        auto bz = BrillouinZone::from_hdf(group, "bz_");
//        return {bz, nest};
//    }

    bool to_hdf(const std::string&, const std::string&, unsigned) const {
//    bool to_hdf(const std::string& filename, const std::string& entry, const unsigned perm=HighFive::File::OpenOrCreate) const {
      throw std::logic_error("to_hdf not implemented yet due to NestNode and NestLeaf");
//        HighFive::File file(filename, perm);
//        return this->to_hdf(file, entry);
      return false;
    }
    static class_t from_hdf(const std::string&, const std::string&) {
      throw std::logic_error("from_hdf not implemented yet due to NestNode and NestLeaf");
//    static BrillouinZoneNest3<T,S> from_hdf(const std::string& filename, const std::string& entry){
//        HighFive::File file(filename, HighFive::File::ReadOnly);
//        return BrillouinZoneNest3<T,S>::from_hdf(file, entry);
    }
#endif // USE_HIGHFIVE
};

}
#endif
