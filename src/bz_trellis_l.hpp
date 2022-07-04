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
    \brief Defines a class to extend `PolyhedronTrellis` with `BrillouinZone` information.
*/
#ifndef BRILLE_BZ_TRELLIS_
#define BRILLE_BZ_TRELLIS_
#include <utility>

#include "bz.hpp"
#include "trellis_poly.hpp"
// #include "phonon.hpp"
namespace brille {
  /*! \brief A PolyhedronTrellis in a BrillouinZone

  The first or irreducible Brillouin zone Polyhedron contained in a BrillouinZone
  object can be used to define the domain of a PolyhedronTrellis hybrid regular
  parallelpiped and triangulation grid.
  The symmetries of the Brillouin zone can then be used to interpolate at any
  point in reciprocal space by finding an equivalent point within the domain of
  the PolyhedronTrellis.
  */
template<class T, class R, class S>
class BrillouinZoneTrellis3: public polytrellis::PolyTrellis<T,R,S,lattice::LVec>{
  using child_t = BrillouinZoneTrellis3<T,R,S>;
  using super_t = polytrellis::PolyTrellis<T,R,S,lattice::LVec>;
  template<class V> using VectorClass = lattice::LVec<V>;
  using vector_t = VectorClass<S>;
  BrillouinZone bz_;
public:
  BrillouinZoneTrellis3(const super_t & pt, BrillouinZone bz): super_t(pt), bz_(std::move(bz)) {}
  BrillouinZoneTrellis3(super_t && pt, BrillouinZone&& bz): super_t(std::move(pt)), bz_(std::move(bz)) {}
  /*! \brief Construct a `BrillouinZoneTrellis3` from a `BrillouinZone` and variable arguments

  All arguments beyond the `BrillouinZone` are passed to the `PolyhedronTrellis` constructor.
  \param bz the `BrillouinZone` used to define the boundaries of the `Mesh3`
  \param args the construction arguments for `PolyhedronTrellis`
  */
  template<typename... A>
  explicit BrillouinZoneTrellis3(const BrillouinZone& bz, A... args): super_t(bz.get_ir_polyhedron(), args...), bz_(bz) {}
  //! get the BrillouinZone object
  [[nodiscard]] BrillouinZone get_brillouinzone() const {return this->bz_;}
  //! get the vertices of the trellis in absolute units
  [[nodiscard]] bArray<S> get_xyz() const {return this->vertices().xyz();}
  //! get the vertices of the inner (cubic) nodes in absolute units
  [[nodiscard]] bArray<S> get_inner_xyz() const {return this->cube_vertices().xyz(); }
  //! get the vertices of the outer (polyhedron) nodes in absolute units
  [[nodiscard]] bArray<S> get_outer_xyz() const {return this->poly_vertices().xyz(); }
  //! get the vertices of the trellis in relative lattice units
  [[nodiscard]] bArray<S> get_hkl() const { return this->vertices().hkl();}
  //! get the vertices of the inner (cubic) nodes in relative lattice units
  [[nodiscard]] bArray<S> get_inner_hkl() const {return this->cube_vertices().hkl(); }
  //! get the vertices of the outer (polyhedron) nodes in relative lattice units
  [[nodiscard]] bArray<S> get_outer_hkl() const {return this->poly_vertices().hkl(); }
  //! get the indices forming the faces of the tetrahedra
  [[nodiscard]] std::vector<std::array<brille::ind_t,4>> get_vertices_per_tetrahedron() const {return this->vertices_per_tetrahedron();}

  /*! \brief Interpolate at an equivalent first Brillouin zone reciprocal lattice point

  \param x        One or more points expressed in the same reciprocal lattice as
                  the stored `BrillouinZone`
  \param nth      the number of parallel OpenMP workers to utilize
  \param no_move  If all provided points are *already* within the first
                  Brillouin zone this optional parameter can be used to skip a
                  call to `BrillouinZone::moveinto`.
  \return a tuple of the interpolated eigenvalue
   s and eigenvectors

  The interpolation is performed by `PolyhedronTrellis::interpolate_at` and no
  post processing is performed.

  \warning The last parameter should only be used with extreme caution as no
           check is performed to ensure that the points are actually in the
           first Brillouin zone. If this condition is not true and the
           parameter is set to true, the subsequent interpolation call may raise
           an error or access unassigned memory and will produce garbage output.
  */
  std::tuple<brille::Array<T>,brille::Array<R>>
  interpolate_at(const vector_t& x, const int nth, const bool no_move=false) const {
    profile_update("Start BrillouinZoneTrellis3::interpolate_at");
    vector_t q(x.type(), x.lattice(), x.size(0));
    VectorClass<int> tau(x.type(), x.lattice(), x.size(0));
    if (no_move){
      // Special mode for testing where no specified points are moved
      // IT IS IMPERATIVE THAT THE PROVIDED POINTS ARE *INSIDE* THE IRREDUCIBLE
      // POLYHEDRON otherwise the interpolation will fail or give garbage back.
      q = x;
    } else if (!bz_.moveinto(x, q, tau, nth)){
      std::string msg;
      msg = "Moving all points into the first Brillouin zone failed.";
      throw std::runtime_error(msg);
    }
    if (nth < 2)
      return this->super_t::interpolate_at(q);
    return this->super_t::interpolate_at(q, nth);
  }

  /*! \brief Interpolate at an equivalent irreducible reciprocal lattice point

  \param x        One or more points expressed in the same reciprocal lattice as
                  the stored `BrillouinZone`
  \param nth      the number of parallel OpenMP workers to utilize
  \param no_move  If all provided points are *already* within the irreducible
                  Brillouin zone this optional parameter can be used to skip a
                  call to `BrillouinZone::ir_moveinto`.
  \return a tuple of the interpolated eigenvalues and eigenvectors

  The interpolation is performed by `PolyhedronTrellis::interpolate_at` and then
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
  std::tuple<brille::Array<T>,brille::Array<R>>
  ir_interpolate_at(const vector_t& x, const int nth, const bool no_move=false) const {
    verbose_update("BZTrellisQ::ir_interpolate_at called with ",nth," threads");
    profile_update("Start BrillouinZoneTrellis3::ir_interpolate_at");
    VectorClass<S> ir_q(x.type(), x.lattice(), x.size(0));
    VectorClass<int> tau(x.type(), x.lattice(), x.size(0));
    std::vector<size_t> rot(x.size(0),0u), invrot(x.size(0),0u);
    if (no_move){
      // Special mode for testing where no specified points are moved
      // IT IS IMPERATIVE THAT THE PROVIDED POINTS ARE *INSIDE* THE IRREDUCIBLE
      // POLYHEDRON otherwise the interpolation will fail or give garbage back.
      ir_q = x;
    } else if (!bz_.ir_moveinto(x, ir_q, tau, rot, invrot, nth)){
      std::string msg;
      msg = "Moving all points into the irreducible Brillouin zone failed.";
      throw std::runtime_error(msg);
    }
    auto [vals, vecs] = (nth>1)
        ? this->super_t::interpolate_at(ir_q, nth)
        : this->super_t::interpolate_at(ir_q);
    profile_update("Apply rotations/permutations to interpolated results");
    // we always need the pointgroup operations to 'rotate'
    PointSymmetry psym = bz_.get_pointgroup_symmetry();
    // and might need the Phonon Gamma table
    GammaTable pgt{GammaTable()};
    if (RotatesLike::Gamma == this->data().vectors().rotateslike()){
      auto cfg = this->approx_config();
      auto s_tol = cfg.template direct<double>();
      auto n_tol = cfg.digit();
      pgt.construct(bz_.get_lattice(), bz_.add_time_reversal(), s_tol, n_tol);
    }
    // make data-sharing Array2 objects for the rotation functions
    brille::Array2<T> vals2(vals);
    brille::Array2<R> vecs2(vecs);
    // actually perform the rotation to Q
    this->data().values() .rotate_in_place(vals2, ir_q, pgt, psym, rot, invrot, nth);
    this->data().vectors().rotate_in_place(vecs2, ir_q, pgt, psym, rot, invrot, nth);
    // we're done so bundle the output
    profile_update("  End BrillouinZoneTrellis3::ir_interpolate_at");
    return std::make_tuple(vals, vecs);
  }

#ifdef USE_HIGHFIVE
    template<class HF>
    std::enable_if_t<std::is_base_of_v<HighFive::Object, HF>, bool>
    to_hdf(HF& obj, const std::string& entry) const{
        auto group = overwrite_group(obj, entry);
        bool ok{true};
        ok &= super_t::to_hdf(group, "trellis");
        ok &= bz_.to_hdf(group, "bz_");
        return ok;
    }
    // Input from HDF5 file/object
    template<class HF>
    static std::enable_if_t<std::is_base_of_v<HighFive::Object, HF>,child_t>
    from_hdf(HF& obj, const std::string& entry){
        auto group = obj.getGroup(entry);
        auto trellis = super_t::from_hdf(group, "trellis");
        auto bz = BrillouinZone::from_hdf(group, "bz_");
        return child_t(trellis, bz);
    }

    [[nodiscard]] bool to_hdf(const std::string& filename, const std::string& entry, const unsigned perm=HighFive::File::OpenOrCreate) const {
        HighFive::File file(filename, perm);
        return this->to_hdf(file, entry);
    }
    static child_t from_hdf(const std::string& filename, const std::string& entry){
        HighFive::File file(filename, HighFive::File::ReadOnly);
        return child_t::from_hdf(file, entry);
    }
#endif // USE_HIGHFIVE

    bool operator!=(const child_t& other) const {
      if (bz_ != other.bz_) return true;
      return this->super_t::operator!=(other);
    }
    bool operator==(const child_t& other) const {
      return !this->operator!=(other);
    }
};

} // end namespace brille
#endif // BRILLE_BZ_TRELLIS_
