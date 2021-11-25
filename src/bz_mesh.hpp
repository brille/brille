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
template<class T, class S>
class BrillouinZoneMesh3: public Mesh3<T,S>{
  using SuperClass = Mesh3<T,S>;
protected:
  BrillouinZone brillouinzone;
public:
  BrillouinZoneMesh3(const SuperClass& pt, BrillouinZone bz): SuperClass(pt), brillouinzone(std::move(bz)) {}
  BrillouinZoneMesh3(SuperClass&& pt, BrillouinZone&& bz): SuperClass(std::move(pt)), brillouinzone(std::move(bz)) {}
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
  //   brillouinzone(bz) {}
  // BrillouinZoneMesh3(const BrillouinZone& bz) Mesh3<T>(bz.get_ir_vertices().get_xyz(), bz.get_ir_vertices_per_face());
  /*! \brief Construct a `BrillouinZoneMesh3` from a `BrillouinZone` and variable arguments

  All arguments beyond the `BrillouinZone` are passed to the `Mesh3` constructor.
  \param bz the `BrillouinZone` used to define the boundaries of the `Mesh3`
  \param args the construction arguments for `Mesh3`
  */
  template<typename... A>
  explicit BrillouinZoneMesh3(const BrillouinZone& bz, A... args):
    SuperClass(bz.get_ir_vertices().get_xyz(), bz.get_ir_vertices_per_face(), args...),
    brillouinzone(bz) {}
  //! \brief Return the BrillouinZone object
  [[nodiscard]] BrillouinZone get_brillouinzone() const {return this->brillouinzone;}
  //! Return the mesh vertices in relative lattice units
  [[nodiscard]] bArray<double> get_mesh_hkl() const {
    auto xyz = this->get_mesh_xyz();
    double toxyz[9], fromxyz[9];
    const BrillouinZone bz = this->get_brillouinzone();
    bz.get_lattice().get_xyz_transform(toxyz);
    if (!brille::utils::matrix_inverse(fromxyz,toxyz)) throw std::runtime_error("transform matrix toxyz has zero determinant");
    auto shape = xyz.shape();
    bArray<double> hkl(shape);
    std::vector<double> tmp(3);
    for (ind_t i=0; i<shape[0]; ++i){
      auto vxyz = xyz.view(i).to_std();
      brille::utils::multiply_matrix_vector<double,double,double,3>(tmp.data(), fromxyz, vxyz.data());
      hkl.set(i, tmp);
    }
    return hkl;
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

#ifdef USE_HIGHFIVE
    template<class HF>
    std::enable_if_t<std::is_base_of_v<HighFive::Object, HF>, bool>
    to_hdf(HF& obj, const std::string& entry) const{
        auto group = overwrite_group(obj, entry);
        bool ok{true};
        ok &= SuperClass::to_hdf(group, "mesh");
        ok &= brillouinzone.to_hdf(group, "brillouinzone");
        return ok;
    }
    // Input from HDF5 file/object
    template<class HF>
    static std::enable_if_t<std::is_base_of_v<HighFive::Object, HF>, BrillouinZoneMesh3<T,S>>
    from_hdf(HF& obj, const std::string& entry){
        auto group = obj.getGroup(entry);
        auto mesh = SuperClass::from_hdf(group, "mesh");
        auto bz = BrillouinZone::from_hdf(group, "brillouinzone");
        return {mesh, bz};
    }

    [[nodiscard]] bool to_hdf(const std::string& filename, const std::string& entry, const unsigned perm=HighFive::File::OpenOrCreate) const {
        HighFive::File file(filename, perm);
        return this->to_hdf(file, entry);
    }
    static BrillouinZoneMesh3<T,S> from_hdf(const std::string& filename, const std::string& entry){
        HighFive::File file(filename, HighFive::File::ReadOnly);
        return BrillouinZoneMesh3<T,S>::from_hdf(file, entry);
    }
#endif //USE_HIGHFIVE
};

} // end namespace brille
#endif // _BZ_MESH_
