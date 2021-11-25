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
template<class T, class S>
class BrillouinZoneNest3: public Nest<T,S>{
  using SuperClass = Nest<T,S>;
  BrillouinZone brillouinzone;
public:
  BrillouinZoneNest3(const SuperClass& pt, BrillouinZone bz): SuperClass(pt), brillouinzone(std::move(bz)) {}
  BrillouinZoneNest3(SuperClass&& pt, BrillouinZone&& bz): SuperClass(std::move(pt)), brillouinzone(std::move(bz)) {}
  /*! \brief Construct a `BrillouinZoneNest3` from a `BrillouinZone` and variable arguments

  All arguments beyond the `BrillouinZone` are passed to the `Nest` constructor.
  \param bz the `BrillouinZone` used to define the boundaries of the `Mesh3`
  \param args the construction arguments for `Nest`
  */
  template<typename... A>
  BrillouinZoneNest3(const BrillouinZone& bz, A... args):
    SuperClass(bz.get_ir_polyhedron(), args...),
    brillouinzone(bz) {}
  //! get the BrillouinZone object
  BrillouinZone get_brillouinzone(void) const {return this->brillouinzone;}
  //! get the vertices of the leaf vertices in inverse Angstrom
  bArray<double> get_xyz(void) const {return this->vertices();}
  //! get the vertices of all vertices in absolute units
  const bArray<double>& get_all_xyz(void) const {return this->all_vertices(); }
  //! get the vertices of the leaf vertices in relative lattice units
  bArray<double> get_hkl(void) const { return xyz_to_hkl(brillouinzone.get_lattice(),this->vertices());}
  //! get the vertices of the inner (cubic) nodes in relative lattice units
  bArray<double> get_all_hkl(void) const {return xyz_to_hkl(brillouinzone.get_lattice(),this->all_vertices()); }
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
  template<class R>
  std::tuple<brille::Array<T>,brille::Array<S>>
  ir_interpolate_at(const LQVec<R>& x, const int nth, const bool no_move=false) const {
    LQVec<R> ir_q(x.get_lattice(), x.size(0));
    LQVec<int> tau(x.get_lattice(), x.size(0));
    std::vector<size_t> rot(x.size(0),0u), invrot(x.size(0),0u);
    if (no_move){
      ir_q = x;
    } else if (!brillouinzone.ir_moveinto(x, ir_q, tau, rot, invrot, nth)){
      std::string msg;
      msg = "Moving all points into the irreducible Brillouin zone failed.";
      throw std::runtime_error(msg);
    }
    auto ir_q_invA = ir_q.get_xyz();
    // perform the interpolation within the irreducible Brillouin zone
    auto [vals, vecs] = (nth > 1)
        ? this->SuperClass::interpolate_at(ir_q_invA, nth)
        : this->SuperClass::interpolate_at(ir_q_invA);
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
    this->data().values().rotate_in_place(vals2, ir_q, pgt, psym, rot, invrot, nth);
    this->data().vectors().rotate_in_place(vecs2, ir_q, pgt, psym, rot, invrot, nth);
    // we're done so bundle the output
    return std::make_tuple(vals, vecs);
  }

#ifdef USE_HIGHFIVE
    template<class HF>
    std::enable_if_t<std::is_base_of_v<HighFive::Object, HF>, bool>
    to_hdf(HF& obj, const std::string& entry) const{
        auto group = overwrite_group(obj, entry);
        bool ok{true};
        ok &= SuperClass::to_hdf(group, "nest");
        ok &= brillouinzone.to_hdf(group, "brillouinzone");
        return ok;
    }
    // Implementing this requires Nest3::to/from_hdf and therefore NestNode and NestLeaf to/from_hdf, which is problematic
//    // Input from HDF5 file/object
//    template<class HF>
//    static std::enable_if_t<std::is_base_of_v<HighFive::Object, HF>, BrillouinZoneNest3<T,S>>
//    from_hdf(HF& obj, const std::string& entry){
//        auto group = obj.getGroup(entry);
//        auto nest = SuperClass::from_hdf(group, "nest");
//        auto bz = BrillouinZone::from_hdf(group, "brillouinzone");
//        return {bz, nest};
//    }

    bool to_hdf(const std::string&, const std::string&, unsigned) const {
//    bool to_hdf(const std::string& filename, const std::string& entry, const unsigned perm=HighFive::File::OpenOrCreate) const {
      throw std::logic_error("to_hdf not implemented yet due to NestNode and NestLeaf");
//        HighFive::File file(filename, perm);
//        return this->to_hdf(file, entry);
      return false;
    }
    static BrillouinZoneNest3<T,S> from_hdf(const std::string&, const std::string&) {
      throw std::logic_error("from_hdf not implemented yet due to NestNode and NestLeaf");
//    static BrillouinZoneNest3<T,S> from_hdf(const std::string& filename, const std::string& entry){
//        HighFive::File file(filename, HighFive::File::ReadOnly);
//        return BrillouinZoneNest3<T,S>::from_hdf(file, entry);
    }
#endif // USE_HIGHFIVE
};

}
#endif
