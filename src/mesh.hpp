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
#ifndef BRILLE_MESH_H_
#define BRILLE_MESH_H_
/*! \file
    \author Greg Tucker
    \brief A class holding a triangulated tetrahedral mesh and data for interpolation
*/
#include "interpolatordual.hpp"
#include <atomic>
#include <queue>
#include <optional>
#include <utility>
#include "lattice_mesh.hpp"
#include "approx_config.hpp"
namespace brille {

/*!
\brief A tetrahedral mesh with eigenvalue and eigenvector data

The mesh is the structured mesh of the irreducible Brillouin zone
(`LatticeMesh`): a lattice grid clipped exactly to the zone, whose boundary
matches itself under the zone's face pairings. If one or more values are defined
for every vertex then it can be used to perform linear interpolation for any
point within the zone. A bucket grid finds the tetrahedron holding a point in
constant time.
*/
template<class DataValues, class DataVectors, class VertexComponents, template<class> class VertexType>
class Mesh3{
  using class_t = Mesh3<DataValues, DataVectors, VertexComponents, VertexType>;
  using mesh_t = LatticeMesh;
  using data_t = DualInterpolator<DataValues, DataVectors>;
  using vert_t = VertexType<VertexComponents>;
  using approx_t = approx_float::Config;
protected:
  mesh_t mesh;
  data_t data_;
  approx_t approx_;
public:
  explicit Mesh3(mesh_t m, approx_t a = approx_float::config): mesh(std::move(m)), approx_(a) {
    data_.initialize_permutation_table(this->size(), this->mesh.collect_keys());
  }
  Mesh3(const class_t& other){
    this->mesh = other.mesh;
    this->data_ = other.data_;
  }
  //Mesh3(TetTri  m, const data_t& d): mesh(std::move(m)), data_(d) {}
  Mesh3(mesh_t m, data_t d, approx_t a): mesh(std::move(m)), data_(std::move(d)), approx_(a) {}
  Mesh3(mesh_t&& m, data_t&& d, approx_t a): mesh(m), data_(d), approx_(a) {}
  //
  auto& operator=(const class_t& other){
    this->mesh = other.mesh;
    this->data_ = other.data_;
    return *this;
  }
  [[nodiscard]] approx_t approx_config() const {return approx_;}
  //! Return the number of mesh vertices
  [[nodiscard]] ind_t size() const { return this->mesh.number_of_vertices(); }
  //! Return the number of mesh vertices
  [[nodiscard]] ind_t vertex_count() const { return this->mesh.number_of_vertices(); }
  //! Return the positions of all vertices in the mesh
  [[nodiscard]] const vert_t& get_mesh_xyz() const{ return this->mesh.get_vertex_positions(); }
  //! Return the tetrahedron indices of the mesh
  [[nodiscard]] const bArray<ind_t>& get_mesh_tetrehedra() const{ return this->mesh.get_vertices_per_tetrahedron();}
  // Get a constant reference to the stored data
  const data_t& data() const {return data_;}
  // Replace the data stored in the object
  template<typename... A> void replace_data(A... args) { data_.replace_data(args...); }
  template<typename... A> void replace_value_data(A... args) { data_.replace_value_data(args...); }
  template<typename... A> void replace_vector_data(A... args) { data_.replace_vector_data(args...); }
  template<typename... A> void set_value_cost_info(A... args) { data_.set_value_cost_info(args...); }
  template<typename... A> void set_vector_cost_info(A... args) {data_.set_vector_cost_info(args...);}
  template<typename... A> void set_vector_normalization(A... args) {data_.set_vector_normalization(args...);}
  //! Return the number of bytes used per Q point
  [[nodiscard]] size_t bytes_per_point() const {return data_.bytes_per_point(); }

  //! Perform sanity checks before attempting to interpolate
  template<class R, template<class> class V> unsigned int check_before_interpolating(const V<R>& x) const {
    unsigned int mask = 0u;
    if (data_.size()==0)
      throw std::runtime_error("The interpolation data must be filled before interpolating.");
    if (x.ndim()!=2 || x.size(1)!=3u)
      throw std::runtime_error("Only (n,3) two-dimensional Q vectors supported in interpolating.");
    if (x.stride().back()!=1)
      throw std::runtime_error("Contiguous vectors required for interpolation.");
    return mask;
  }
  //! Perform linear interpolation at the specified Reciprocal lattice points
//  template<class R>
//  std::tuple<brille::Array<DataValues>,brille::Array<DataVectors>>
//  interpolate_at(const lattice::LVec<R>& x) const {return this->interpolate_at(x.xyz());}
  //! Perform linear interpolating at the specified points in the mesh's orthonormal frame
  std::tuple<brille::Array<DataValues>,brille::Array<DataVectors>>
  interpolate_at(const vert_t& x) const {
    this->check_before_interpolating(x);
    auto valsh = data_.values().shape();
    auto vecsh = data_.vectors().shape();
    valsh[0] = x.size(0);
    vecsh[0] = x.size(0);
    brille::Array<DataValues> vals(valsh);
    brille::Array<DataVectors> vecs(vecsh);
    // vals and vecs are row-ordered contiguous by default, so we can create
    // mutable data-sharing Array2 objects for use with
    // Interpolator2::interpolate_at through the constructor:
    brille::Array2<DataValues> vals2(vals);
    brille::Array2<DataVectors> vecs2(vecs);
    for (ind_t i=0; i<x.size(0); ++i){
      verbose_update("Locating ",x.to_string(i));
      auto verts_weights = this->mesh.locate(x.view(i));
      if (verts_weights.size()<1){
        debug_update("Point ",x.to_string(i)," not found in tetrahedra!");
        throw std::runtime_error("Point not found in tetrahedral mesh");
      }
      data_.interpolate_at(verts_weights, vals2, vecs2, i);
    }
    return std::make_tuple(vals, vecs);
  }
  std::tuple<brille::Array<DataValues>,brille::Array<DataVectors>>
  interpolate_at(const vert_t& x, const int threads) const {
    this->check_before_interpolating(x);
    // not used in parallel region
    auto valsh = data_.values().shape();
    auto vecsh = data_.vectors().shape();
    valsh[0] = x.size(0);
    vecsh[0] = x.size(0);
    // shared between threads
    Array<DataValues> vals(valsh);
    Array<DataVectors> vecs(vecsh);
    // vals and vecs are row-ordered contiguous by default, so we can create
    // mutable data-sharing Array2 objects for use with
    // Interpolator2::interpolate_at through the constructor:
    Array2<DataValues> vals2(vals);
    Array2<DataVectors> vecs2(vecs);

    const auto pool = ThreadPool::getInstance();
    if (threads > 0) pool->resize(threads); else pool->resize();
    const auto workers = pool->size();
    std::atomic<size_t> missing{0};
    auto task = [&](const size_t worker) {
      auto [f, l] = thread_slice(x.size(0), workers, worker);
      return [&,first=f,last=l]() {
        for (size_t i=first; i<last; ++i) {
          // round-off can put a point on the mesh surface just outside every tetrahedron
          if (auto verts_weights = mesh.locate(x.view(i)); verts_weights.size()) {
            data_.interpolate_at(verts_weights, vals2, vecs2, i);
          } else {
            ++missing;
          }
        }
      };
    };
    for (size_t i=0; i<workers; ++i) pool->enqueue(task(i));
    pool->wait();

    if (missing) {
      throw std::runtime_error(std::to_string(missing.load()) + " of " + std::to_string(x.size(0))
                               + " points not found in tetrahedral mesh");
    }
    return std::make_tuple(vals, vecs);
  }
  //! Return the neighbours for which a passed boolean array holds true
  // template<typename R> std::vector<ind_t> which_neighbours(const std::vector<R>& t, const R value, const ind_t idx) const;
  [[nodiscard]] std::string to_string() const {
    std::string str= data_.to_string();
    str += " for the points of a " + mesh.to_string();
    return str;
  }
  void sort() {data_.sort();}

  template<class HF>
  std::enable_if_t<std::is_base_of_v<HighFive::Object, HF>, bool>
  to_hdf(HF& obj, const std::string& entry) const{
    auto group = overwrite_group(obj, entry);
    bool ok{true};
    ok &= mesh.to_hdf(group, "triangulation");
    ok &= data_.to_hdf(group, "data");
    ok &= approx_.to_hdf(group,"approx");
    return ok;
  }
  // Input from HDF5 file/object
  template<class HF>
  static std::enable_if_t<std::is_base_of_v<HighFive::Object, HF>, class_t>
  from_hdf(HF& obj, const std::string& entry){
    auto group = obj.getGroup(entry);
    auto m = mesh_t::from_hdf(group, "triangulation");
    auto d = data_t::from_hdf(group, "data");
    auto a = approx_t::from_hdf(group, "approx");
    return class_t(m, d, a);
  }

  //! Whether the mesh can be refined (not one read from a file without what built it)
  [[nodiscard]] bool refinable() const { return mesh.refinable(); }
  //! Whether the triangulation refinement works on is held; see `release_triangulation`
  [[nodiscard]] bool holds_triangulation() const { return mesh.holds_triangulation(); }
  //! Free the triangulation refinement works on; it is rebuilt when next needed
  void release_triangulation() { mesh.release_triangulation(); }
  /*! \brief The vertices (Cartesian) that refining the tetrahedra `tets` would add

  Nothing is changed. `refine` with the same arguments adds exactly these, in this
  order, after the existing vertices.
  */
  [[nodiscard]] vert_t refinement_points(const std::vector<ind_t>& tets, const double min_edge) const {
    return mesh.plan(tets, min_edge).points;
  }
  /*! \brief Bisect the tetrahedra `tets`, with closure, adding vertices after the existing ones

  \param tets the tetrahedra to split, by index
  \param min_edge tetrahedra whose longest edge is at most twice this (Å⁻¹) are not split
  \param values,vectors data for the new vertices, laid out like the data held;
         required if the mesh holds data, and then for exactly the vertices added
  \return the new vertices (Cartesian)

  Existing vertices keep their indices and data. The mode permutations found by
  `sort` are reset: sort again after refining.
  */
  vert_t refine(const std::vector<ind_t>& tets, const double min_edge,
                const std::optional<brille::Array<DataValues>>& values = std::nullopt,
                const std::optional<brille::Array<DataVectors>>& vectors = std::nullopt){
    return apply_refinement(plan_refinement(tets, min_edge), values, vectors);
  }
  using refinement_t = typename mesh_t::Refinement;
  //! The expensive half of `refine`, which changes nothing
  [[nodiscard]] refinement_t plan_refinement(const std::vector<ind_t>& tets, const double min_edge) const {
    return mesh.plan(tets, min_edge);
  }
  /*! The other half of `refine`. It replaces the held data, which may be the last
  reference to a Python buffer, so a binding must hold the GIL while calling it. */
  vert_t apply_refinement(refinement_t && plan,
                          const std::optional<brille::Array<DataValues>>& values = std::nullopt,
                          const std::optional<brille::Array<DataVectors>>& vectors = std::nullopt){
    const ind_t added = plan.points.size(0);
    const bool filled = data_.size() > 0;
    if (filled) {
      if (!values || !vectors)
        throw std::runtime_error("The mesh holds data: give the values and vectors for the new vertices too");
      if (values->size(0) != added || vectors->size(0) != added)
        throw std::runtime_error("Refining adds " + std::to_string(added) + " vertices, but data for "
                                 + std::to_string(values->size(0)) + " values and " + std::to_string(vectors->size(0)) + " vectors was given");
    } else if (values || vectors) {
      throw std::runtime_error("The mesh holds no data yet: fill it before giving data for new vertices");
    }
    vert_t points = plan.points;
    mesh.commit(std::move(plan));
    if (filled) data_.append(*values, *vectors);
    data_.initialize_permutation_table(this->size(), mesh.collect_keys());
    return points;
  }

  /*! \brief Whether `max_points` made the mesh coarser than `max_size` asked for

  The mesh is valid, but its tetrahedra are larger than requested.
  */
  [[nodiscard]] bool refinement_limited() const {return mesh.refinement_limited();}
};

} // namespace brille
#endif // BRILLE_MESH_H_
