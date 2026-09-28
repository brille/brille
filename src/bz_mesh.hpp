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
    \brief Defines a class to extend `Mesh3` with `BrillouinZone` information.
*/
#ifndef BRILLE_BZ_MESH_
#define BRILLE_BZ_MESH_
#include "bz.hpp"
#include "mesh.hpp"

#include <utility>
namespace brille {
/*! \brief A Mesh3 in a BrillouinZone

The first or irreducible Brillouin zone Polyhedron contained in a BrillouinZone
object can be used to define the domain of a Mesh3 triangulation.
The symmetries of the Brillouin zone can then be used to interpolate at any
point in reciprocal space by finding an equivalent point within the triangulated
domain.
*/
template<class T, class S, class V>
class BrillouinZoneMesh3: public Mesh3<T,S,V,Array2>{
  using class_t = BrillouinZoneMesh3<T,S,V>;
  using base_t = Mesh3<T,S,V,Array2>;
  template<class A> using lv_t = lattice::LVec<A>;
  template<class A> using bv_t = Array2<A>;
protected:
  BrillouinZone bz_;
public:
  BrillouinZoneMesh3(const base_t& pt, BrillouinZone bz): base_t(pt), bz_(std::move(bz)) {}
  BrillouinZoneMesh3(base_t&& pt, BrillouinZone&& bz): base_t(std::move(pt)), bz_(std::move(bz)) {}
  /*! \brief The structured mesh of a `BrillouinZone`'s irreducible zone

  \param bz the `BrillouinZone`, whose irreducible polyhedron the mesh fills exactly
  \param max_size the largest grid tetrahedron volume, in Å⁻³, which sets the grid
         spacing; not positive for the coarsest grid
  \param num_levels unused; the TetGen mesh it controlled is gone, and points are
         located in constant time without layers
  \param max_points if positive, the grid is coarsened until its estimated vertex
         count is at most this (see `refinement_limited`)
  */
  explicit BrillouinZoneMesh3(const BrillouinZone& bz, const double max_size=-1., [[maybe_unused]] const int num_levels=3,
                              const int max_points=-1):
    base_t(LatticeMesh::from_zone(bz, max_size, max_points)), bz_(bz) {}
  //! The vertices that refining the tetrahedra `tets` would add, in relative lattice units
  [[nodiscard]] bv_t<V> refinement_points_hkl(const std::vector<ind_t>& tets, const double min_edge) const {
    return from_xyz_like(LengthUnit::inverse_angstrom, bz_.get_lattice(), this->refinement_points(tets, min_edge)).hkl();
  }
  //! Apply a refinement (see Mesh3::apply_refinement), returning the added vertices in relative lattice units
  template<class... A>
  bv_t<V> apply_refinement_hkl(typename base_t::refinement_t && plan, A... args) {
    return from_xyz_like(LengthUnit::inverse_angstrom, bz_.get_lattice(), this->apply_refinement(std::move(plan), args...)).hkl();
  }
  //! \brief Return the BrillouinZone object
  [[nodiscard]] BrillouinZone get_brillouinzone() const {return this->bz_;}
  //! Return the mesh vertices in relative lattice units
  [[nodiscard]] bv_t<V> get_mesh_hkl() const {
    return from_xyz_like(LengthUnit::inverse_angstrom, bz_.get_lattice(), this->get_mesh_xyz()).hkl();
  }

  /*! \brief Interpolate at an equivalent irreducible reciprocal lattice point

  \param x        One or more points expressed in the same reciprocal lattice as
                  the stored `BrillouinZone`
  \param nthreads the number of parallel threads to utilize
  \param no_move  If all provided points are *already* within the irreducible
                  Brillouin zone this optional parameter can be used to skip a
                  call to `BrillouinZone::ir_moveinto`.
  \return a tuple of the interpolated eigenvalues and eigenvectors

  The interpolation is performed by `Mesh3::interpolate_at` and then
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
  ir_interpolate_at(const lv_t<V>& x, Args... args) const {
    lv_t<V> ir_q(x.type(), x.lattice(), x.size(0));
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
    if constexpr (NO_MOVE) {
      // set rot and invrot to the identity symmetry operation (which is not necessarily the 0th one)
      auto identity = psym.find_identity_index();
      std::fill(rot.begin(), rot.end(), identity);
      std::fill(invrot.begin(), invrot.end(), identity);
    }
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
    this->data().values() .rotate_in_place(vals2, ir_q, pgt, psym, rot, invrot, args...);
    this->data().vectors().rotate_in_place(vecs2, ir_q, pgt, psym, rot, invrot, args...);
    // we're done so bundle the output
    return std::make_tuple(vals, vecs);
  }

    template<class HF>
    std::enable_if_t<std::is_base_of_v<HighFive::Object, HF>, bool>
    to_hdf(HF& obj, const std::string& entry) const{
        auto group = overwrite_group(obj, entry);
        bool ok{true};
        ok &= base_t::to_hdf(group, "mesh");
        ok &= bz_.to_hdf(group, "bz_");
        return ok;
    }
    // Input from HDF5 file/object
    template<class HF>
    static std::enable_if_t<std::is_base_of_v<HighFive::Object, HF>, class_t>
    from_hdf(HF& obj, const std::string& entry){
        auto group = obj.getGroup(entry);
        auto mesh = base_t::from_hdf(group, "mesh");
        auto bz = BrillouinZone::from_hdf(group, "bz_");
        return {mesh, bz};
    }

    [[nodiscard]] bool to_hdf(const std::string& filename, const std::string& entry, const unsigned perm=HighFive::File::OpenOrCreate) const {
        HighFive::File file(filename, perm);
        return this->to_hdf(file, entry);
    }
    static class_t from_hdf(const std::string& filename, const std::string& entry){
        HighFive::File file(filename, HighFive::File::ReadOnly);
        return class_t::from_hdf(file, entry);
    }

};

} // end namespace brille
#endif // _BZ_MESH_
