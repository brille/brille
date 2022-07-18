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
  /* Construct using a maximum tetrahedron volume -- makes a tetrahedron mesh
      instead of a orthogonal grid.
      @param bz The BrillouinZone object
      @param vol The maximum tetrahedron volume
      @param isrlu A flag to indicate if vol is in units of rlu (isrlu=1) or inverse angstrom (isrlu=0)
      @note If vol is in relative lattice units an absolute volume will be
            determined using the unit cell volume of the underlying lattice.
  */
  // BrillouinZoneMesh3(const BrillouinZone& bz, const double max_size_invA=-1., const double min_angle=20.0, const double max_angle=-1.0, const double max_ratio=-1., const int max_points=-1):
  //   Mesh3<T>(bz.get_ir_vertices().get_xyz(), bz.get_ir_vertices_per_face(), max_size_invA, min_angle, max_angle, max_ratio, max_points),
  //   bz_(bz) {}
  // BrillouinZoneMesh3(const BrillouinZone& bz) Mesh3<T>(bz.get_ir_vertices().get_xyz(), bz.get_ir_vertices_per_face());
  /*! \brief Construct a `BrillouinZoneMesh3` from a `BrillouinZone` and variable arguments

  All arguments beyond the `BrillouinZone` are passed to the `Mesh3` constructor.
  \param bz the `BrillouinZone` used to define the boundaries of the `Mesh3`
  \param args the construction arguments for `Mesh3`
  */
  template<typename... A>
  explicit BrillouinZoneMesh3(const BrillouinZone& bz, A... args):
    base_t(bz.get_ir_vertices().xyz(), bz.get_ir_vertices_per_face(), args...),
    bz_(bz) {}
  //! \brief Return the BrillouinZone object
  [[nodiscard]] BrillouinZone get_brillouinzone() const {return this->bz_;}
  //! Return the mesh vertices in relative lattice units
  [[nodiscard]] bv_t<V> get_mesh_hkl() const {
    return from_xyz_like(LengthUnit::inverse_angstrom, bz_.get_lattice(), this->get_mesh_xyz()).hkl();
  }

  /*! \brief Interpolate at an equivalent irreducible reciprocal lattice point

  \param x        One or more points expressed in the same reciprocal lattice as
                  the stored `BrillouinZone`
  \param nthreads the number of parallel OpenMP workers to utilize
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

#ifdef USE_HIGHFIVE
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
#endif //USE_HIGHFIVE
};

} // end namespace brille
#endif // _BZ_MESH_
