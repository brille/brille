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


template<class T> class Mesh3{
protected:
  TetTri mesh;
  InterpolationData<T> data_;
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
  }
  Mesh3(const Mesh3<T>& other){
    this->mesh = other.mesh;
    this->data_ = other.data_;
  }
  Mesh3<T>& operator=(const Mesh3<T>& other){
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
  const InterpolationData<T>& data(void) const {return data_;}
  // Replace the data stored in the object
  template<typename... A> void replace_data(A... args) { data_.replace_data(args...); }
  // Calculate the Debye-Waller factor for the provided Q points and ion masses
  template<template<class> class A>
  ArrayVector<double> debye_waller(const A<double>& Q, const std::vector<double>& M, const double t_K) const{
    return data_.debye_waller(Q,M,t_K);
  }

  //! Perform sanity checks before attempting to interpolate
  template<typename R> unsigned int check_before_interpolating(const ArrayVector<R>& x) const;
  //! Perform linear interpolation at the specified Reciprocal lattice points
  template<typename R> ArrayVector<T> interpolate_at(const LQVec<R>& x) const {return this->interpolate_at(x.get_xyz());}
  //! Perform linear interpolating at the specified points in the mesh's orthonormal frame
  template<typename R> ArrayVector<T> interpolate_at(const ArrayVector<R>& x) const;
  template<typename R> ArrayVector<T> parallel_interpolate_at(const ArrayVector<R>& x, const int nthreads) const;
  //! Return the neighbours for which a passed boolean array holds true
  template<typename R> std::vector<size_t> which_neighbours(const std::vector<R>& t, const R value, const size_t idx) const;
  //! Sort the values stored in the mesh by already-sorted neighbour consensus
  template<typename R=double>
  ArrayVector<size_t> multi_sort_perm(const R s=R(1), const R v=R(1), const R m=R(1), const int vf=0) const;
  template<typename R> size_t consensus_sort_from(const size_t s, const R w[3],
    const int f, const size_t so[2], ArrayVector<size_t>& p,
    std::vector<bool>& d, std::vector<bool>& l, std::vector<size_t>& v) const;
  template<typename R> bool consensus_sort_difference(const R w[3],
    const int f, const size_t so[2], ArrayVector<size_t>& p,
    std::vector<bool>& d, const size_t i, const std::vector<size_t> n) const;
  std::string to_string(void) const {
    std::string str= data_.to_string();
    str += " for the points of a TetTri[" + mesh.to_string() + "]";
    return str;
  }
private:
};

#include "mesh.tpp"

#endif // _MESH_H_
