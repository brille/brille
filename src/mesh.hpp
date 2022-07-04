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
// #include <set>
// #include "array.hpp"
// #include "array2.hpp"
// #include <vector>
// #include <array>
// #include <omp.h>
#include "interpolatordual.hpp"
// #include "utilities.hpp"
// #include "permutation.hpp"
#include <queue>
#include <utility>
#include "triangulation_layers.hpp"
#include "polyhedron_flex.hpp"
#include "approx_config.hpp"
namespace brille {

/*!
\brief A triangulated tetrahedral mesh with eigenvalue and eigenvector data

One way of dividing three dimensional space is to fill it with a tiling of
tetrehedra. Each tetrahedra has four vertices and four triangular faces.
As each face can be shared with one other tetrahedra, each tetrahedron has up
to four neighbours (tetrahedra on the surface of a space will have one fewer
neighbour per surface). There is no limit to how many tetrahedra any vertex
can contribute to.

If one or more values are defined for every vertex in the tetrahedral mesh then
it can be used to perform linear interpolation for any arbitrary point within
the bound space.

With no guaranteed ordering of the tethrahedra finding which tetrahedra contains
the interpolation point can require testing all tetrahedra for inclusion.
As such a simple tetrahedral mesh is not well suited for *fast* interpolation.
In order to overcome this limitation, this class uses a hierarchy of overlapping
tetrahedra arranged in layers to limit the number of inclusion tests required
during interpolation.
If a point is within a tetrahedron at a given layer then a list of tetrahedra
that it might be in at the next lower layer is available. These next-lower
tetrahedra all have the property that they have an non-null intersection with
the connected higher-level tetrahedra.
*/
template<class DataValues, class DataVectors, class VertexComponents, template<class> class VertexType>
class Mesh3{
  using class_t = Mesh3<DataValues, DataVectors, VertexComponents, VertexType>;
  using mesh_t = TetTri;
  using data_t = DualInterpolator<DataValues, DataVectors>;
  using vert_t = VertexType<VertexComponents>;
  using approx_t = approx_float::Config;
  using poly_t = polyhedron::Poly<VertexComponents, VertexType>;
protected:
  mesh_t mesh;
  data_t data_;
  approx_t approx_;
public:
  template<class... Args>
  Mesh3(const vert_t& vertices, Args... args){
    this->construct(vertices, args...);
  }
  template<class... Args>
  Mesh3(const poly_t& poly, Args... args){
    this->construct(poly.vertices(), poly.facets().facets(), args...);
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
    omp_set_num_threads( (threads > 0) ? threads : omp_get_max_threads() );
    this->check_before_interpolating(x);
    // not used in parallel region
    auto valsh = data_.values().shape();
    auto vecsh = data_.vectors().shape();
    valsh[0] = x.size(0);
    vecsh[0] = x.size(0);
    // shared between threads
    brille::Array<DataValues> vals(valsh);
    brille::Array<DataVectors> vecs(vecsh);
    // vals and vecs are row-ordered contiguous by default, so we can create
    // mutable data-sharing Array2 objects for use with
    // Interpolator2::interpolate_at through the constructor:
    brille::Array2<DataValues> vals2(vals);
    brille::Array2<DataVectors> vecs2(vecs);
    // OpenMP < v3.0 (VS uses v2.0) requires signed indexes for omp parallel
    long xsize = brille::utils::u2s<long, ind_t>(x.size(0));
#pragma omp parallel for default(none) shared(x, vals2, vecs2, xsize) schedule(dynamic)
    for (long si=0; si<xsize; ++si){
      auto i = brille::utils::s2u<ind_t, long>(si);
      auto verts_weights = this->mesh.locate(x.view(i));
      data_.interpolate_at(verts_weights, vals2, vecs2, i);
    }
    return std::make_tuple(vals, vecs);
  }
  //! Return the neighbours for which a passed boolean array holds true
  // template<typename R> std::vector<ind_t> which_neighbours(const std::vector<R>& t, const R value, const ind_t idx) const;
  [[nodiscard]] std::string to_string() const {
    std::string str= data_.to_string();
    str += " for the points of a TetTri[" + mesh.to_string() + "]";
    return str;
  }
  void sort() {data_.sort();}

#ifdef USE_HIGHFIVE
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
#endif // USE_HIGHFIVE

private:
  template<class I>
  void construct(const vert_t& vertices,
                 const std::vector<std::vector<I>>& facets,
                 const double max_volume,
                 const int num_levels,
                 const int max_points,
                 approx_t a){
    this->mesh = triangulate(vertices, facets, max_volume, num_levels, max_points);
    data_.initialize_permutation_table(this->size(), this->mesh.collect_keys());
    approx_ = a;
  }
  template<class I>
  void construct(const vert_t& vertices,
                 const std::vector<std::vector<I>>& facets,
                 const double max_volume,
                 const int num_levels,
                 const int max_points
                 ){
    this->construct(vertices, facets, max_volume, num_levels, max_points, approx_float::config);
  }
  template<class I>
  void construct(const vert_t& vertices,
                 const std::vector<std::vector<I>>& facets,
                 const double max_volume,
                 const int num_levels
  ){
    this->construct(vertices, facets, max_volume, num_levels, -1);
  }
  template<class I>
  void construct(const vert_t& vertices,
                 const std::vector<std::vector<I>>& facets,
                 const double max_volume
  ){
    this->construct(vertices, facets, max_volume, 5);
  }
  template<class I>
  void construct(const vert_t& vertices,
                 const std::vector<std::vector<I>>& facets
  ){
    this->construct(vertices, facets, -1.);
  }
};

} // namespace brille
#endif // BRILLE_MESH_H_
