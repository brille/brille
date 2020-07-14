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

#ifndef _MESH_H_
#define _MESH_H_

#include <set>
#include "latvec.hpp"
#include <vector>
#include <array>
#include <omp.h>
#include "interpolation.hpp"
#include "interpolation_data.hpp"
#include "utilities.hpp"
#include "permutation.hpp"
#include <queue>
#include "triangulation_layers.hpp"


template<class T, class S> class Mesh3{
protected:
  TetTri mesh;
  InterpolationData<T,S> data_;
public:
  //constructors
  Mesh3(const ArrayVector<double>& verts,
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
  Mesh3<T,S>& operator=(const Mesh3<T,S>& other){
    this->mesh = other.mesh;
    this->data_ = other.data_;
    return *this;
  }
  //! Return the number of mesh vertices
  size_t size() const { return this->mesh.number_of_vertices(); }
  //! Return the positions of all vertices in the mesh
  const ArrayVector<double>& get_mesh_xyz() const{ return this->mesh.get_vertex_positions(); }
  //! Return the tetrahedron indices of the mesh
  const ArrayVector<size_t>& get_mesh_tetrehedra() const{ return this->mesh.get_vertices_per_tetrahedron();}
  // Get a constant reference to the stored data
  const InterpolationData<T,S>& data(void) const {return data_;}
  // Replace the data stored in the object
  // template<typename... A> void replace_data(A... args) { data_.replace_data(args...); }
  template<typename... A> void replace_value_data(A... args) { data_.replace_value_data(args...); }
  template<typename... A> void replace_vector_data(A... args) { data_.replace_vector_data(args...); }
  template<typename... A> void set_value_cost_info(A... args) { data_.set_value_cost_info(args...); }
  template<typename... A> void set_vector_cost_info(A... args) {data_.set_vector_cost_info(args...);}
  // Calculate the Debye-Waller factor for the provided Q points and ion masses
  template<template<class> class A>
  ArrayVector<double> debye_waller(const A<double>& Q, const std::vector<double>& M, const double t_K) const{
    return data_.debye_waller(Q,M,t_K);
  }

  //! Perform sanity checks before attempting to interpolate
  template<typename R> unsigned int check_before_interpolating(const ArrayVector<R>& x) const;
  //! Perform linear interpolation at the specified Reciprocal lattice points
  template<typename R>
  std::tuple<ArrayVector<T>,ArrayVector<S>>
  interpolate_at(const LQVec<R>& x) const {return this->interpolate_at(x.get_xyz());}
  //! Perform linear interpolating at the specified points in the mesh's orthonormal frame
  template<typename R>
  std::tuple<ArrayVector<T>,ArrayVector<S>>
  interpolate_at(const ArrayVector<R>& x) const;
  template<typename R>
  std::tuple<ArrayVector<T>,ArrayVector<S>>
  parallel_interpolate_at(const ArrayVector<R>& x, const int nthreads) const;
  //! Return the neighbours for which a passed boolean array holds true
  template<typename R> std::vector<size_t> which_neighbours(const std::vector<R>& t, const R value, const size_t idx) const;
  std::string to_string(void) const {
    std::string str= data_.to_string();
    str += " for the points of a TetTri[" + mesh.to_string() + "]";
    return str;
  }
  void sort() {data_.sort();}
private:
};

#include "mesh.tpp"

#endif // _MESH_H_
