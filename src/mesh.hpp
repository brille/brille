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
// #include "array_latvec.hpp" // defines bArray
// #include <vector>
// #include <array>
// #include <omp.h>
#include "interpolatordual.hpp"
// #include "utilities.hpp"
// #include "permutation.hpp"
#include <queue>
#include <utility>
#include "triangulation_layers.hpp"
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
template<class T, class S>
class Mesh3{
  using data_t = DualInterpolator<T,S>;
protected:
  TetTri mesh;
  data_t data_;
public:
  /* \brief Triangulate a space and build a TetTri within it

  \param verts the vertices bounding the space
  \param facets the vertex lists of the surfaces of the boundary
  \param max_volume The largest volume that a leaf tetrahedron can have
  \param num_levels The number of branchings
  \param max_points The largest number of triangulated vertices
  */
  Mesh3(const bArray<double>& verts,
        const std::vector<std::vector<int>>& facets,
        const double max_volume=-1.0,
        const int num_levels=3,
        const int max_points=-1
      ){
    // this->mesh = triangulate(verts, facets, max_volume, min_angle, max_angle, min_ratio, max_points, trellis_fraction);
    this->mesh = triangulate(verts, facets, max_volume, num_levels, max_points);
    data_.initialize_permutation_table(this->size(), this->mesh.collect_keys());
  }
  Mesh3(const Mesh3<T,S>& other){
    this->mesh = other.mesh;
    this->data_ = other.data_;
  }
  //Mesh3(TetTri  m, const data_t& d): mesh(std::move(m)), data_(d) {}
  Mesh3(TetTri m, data_t d): mesh(std::move(m)), data_(std::move(d)) {}
  Mesh3(TetTri&& m, data_t&& d): mesh(m), data_(d) {}
  //
  Mesh3<T,S>& operator=(const Mesh3<T,S>& other){
    this->mesh = other.mesh;
    this->data_ = other.data_;
    return *this;
  }
  //! Return the number of mesh vertices
  ind_t size() const { return this->mesh.number_of_vertices(); }
  //! Return the number of mesh vertices
  ind_t vertex_count() const { return this->mesh.number_of_vertices(); }
  //! Return the positions of all vertices in the mesh
  const bArray<double>& get_mesh_xyz() const{ return this->mesh.get_vertex_positions(); }
  //! Return the tetrahedron indices of the mesh
  const bArray<ind_t>& get_mesh_tetrehedra() const{ return this->mesh.get_vertices_per_tetrahedron();}
  // Get a constant reference to the stored data
  const data_t& data(void) const {return data_;}
  // Replace the data stored in the object
  template<typename... A> void replace_data(A... args) { data_.replace_data(args...); }
  template<typename... A> void replace_value_data(A... args) { data_.replace_value_data(args...); }
  template<typename... A> void replace_vector_data(A... args) { data_.replace_vector_data(args...); }
  template<typename... A> void set_value_cost_info(A... args) { data_.set_value_cost_info(args...); }
  template<typename... A> void set_vector_cost_info(A... args) {data_.set_vector_cost_info(args...);}
  //! Return the number of bytes used per Q point
  size_t bytes_per_point() const {return data_.bytes_per_point(); }
  // Calculate the Debye-Waller factor for the provided Q points and ion masses
  template<template<class> class A>
  brille::Array<double> debye_waller(const A<double>& Qpts, const std::vector<double>& Masses, const double t_K) const{
    return data_.debye_waller(Qpts,Masses,t_K);
  }

  //! Perform sanity checks before attempting to interpolate
  template<class R> unsigned int check_before_interpolating(const bArray<R>& x) const;
  //! Perform linear interpolation at the specified Reciprocal lattice points
  template<class R>
  std::tuple<brille::Array<T>,brille::Array<S>>
  interpolate_at(const LQVec<R>& x) const {return this->interpolate_at(x.get_xyz());}
  //! Perform linear interpolating at the specified points in the mesh's orthonormal frame
  template<class R>
  std::tuple<brille::Array<T>,brille::Array<S>>
  interpolate_at(const bArray<R>& x) const;
  template<class R>
  std::tuple<brille::Array<T>,brille::Array<S>>
  parallel_interpolate_at(const bArray<R>& x, const int nthreads) const;
  //! Return the neighbours for which a passed boolean array holds true
  // template<typename R> std::vector<ind_t> which_neighbours(const std::vector<R>& t, const R value, const ind_t idx) const;
  std::string to_string(void) const {
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
    return ok;
  }
  // Input from HDF5 file/object
  template<class HF>
  static std::enable_if_t<std::is_base_of_v<HighFive::Object, HF>, Mesh3<T,S>>
  from_hdf(HF& obj, const std::string& entry){
    auto group = obj.getGroup(entry);
    auto m = TetTri::from_hdf(group, "triangulation");
    auto d = data_t::from_hdf(group, "data");
    return Mesh3(m, d);
  }
#endif // USE_HIGHFIVE
};

#include "mesh.tpp"

} // namespace brille
#endif // _MESH_H_
