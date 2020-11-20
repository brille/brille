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
/*! \file */
#ifndef BRILLE_INTERPOLATION_H_
#define BRILLE_INTERPOLATION_H_
#include <iostream>
#include <vector>
#include <cmath>
/*! \brief Find the linear index and weights of interpolation points for a given position

Linear interpolation in D dimensions requires 2ᴰ points surrounding the position
where an interpolated value is to be determined. The value at each of the 2ᴰ
points contributes in proportion to how close it is to the interpolation point.
This function finds the linear indices for any object which has a `sub2map`
method with signature `(const size_t* subscripted_index, size_t* linear_index)`
and determines the appropriate normalized weights for each of the 2ᴰ points.

@param that A pointer to an object with a `sub2map(const size_t*, size_t&)`
            method: MapGrid3, MapGrid4, InterpolateGrid3, InterpolateGrid4,
                    BrillouinZoneGrid3, BrillouinZoneGrid4, etc.
@param      zero  Pointer to the zero-point values of the grid, with N elements
@param      step  Pointer to the step-size values of the grid, with N elements
@param      ijk   Pointer to the nearest grid subscripted index to the
                  interpolation point, with N elements
@param      x     Pointer to the coordinates of the interpolation point,
                  with N elements
@param[out] c     Pointer where the linear indices will be stored,
                  with 2ᴺ elements
@param[out] w     Pointer where the weights will be stored, with 2ᴺ elements
@param      N     Dimensionality of the gridded space
@param      dirs  Directions in the gridded space over which interpolation is
                  to be performed; D=1≤`dirs.size()`≤N
@returns An unsigned integer with bits acting as flags to indicate whether a
         point required for the interpolation is out-of-bounds. The order of
         bits depends on `dirs.size()`, see the source if you need to decode
         this information.
@note This function supports linear, bilinear, trilinear, and quadralinear
      interpolation, e.g., 1≤D≤4, in an arbitrary number of dimensions with D≤N;
      extending it to higher-dimensional interpolation should be possible
      if ever desired.
*/
template <class R, template<class> class T>
int corners_and_weights(const T<R>* that, const double* zero, const double* step, const size_t *ijk, const double *x, size_t *c, double *w, const size_t N, const std::vector<size_t>& dirs){
    size_t ndims = dirs.size();
    std::vector<double> p(ndims), m(ndims);
    std::vector<int> d(ndims);
    for (size_t i=0; i<ndims; ++i){
      // tmp = interpolation_direction_and_distance(zero[dirs[i]],step[dirs[i]],ijk[dirs[i]],x[dirs[i]]);
      double tmp = (x[dirs[i]]-(zero[dirs[i]]+ijk[dirs[i]]*step[dirs[i]]))/step[dirs[i]];
      d[i] = tmp < 0 ? -1 : 1;
      p[i] = std::abs(tmp);
      m[i] = 1.0 - p[i];
    }
    std::vector<size_t> t(N);
    for (size_t i=0; i<N; ++i) t[i] = ijk[i];
    int oob=0;
    switch (ndims){
      case 4:
                          oob +=       that->sub2map(t.data(),c[ 0u]); w[0] = m[0]*m[1]*m[2]*m[3]; // (0000)
      t[dirs[0]] += d[0]; oob +=     2*that->sub2map(t.data(),c[ 1u]); w[1] = p[0]*m[1]*m[2]*m[3]; // (1000)
      t[dirs[1]] += d[1]; oob +=     4*that->sub2map(t.data(),c[ 2u]); w[2] = p[0]*p[1]*m[2]*m[3]; // (1100)
      t[dirs[0]] -= d[0]; oob +=     8*that->sub2map(t.data(),c[ 3u]); w[3] = m[0]*p[1]*m[2]*m[3]; // (0100)
      t[dirs[2]] += d[2]; oob +=    16*that->sub2map(t.data(),c[ 4u]); w[4] = m[0]*p[1]*p[2]*m[3]; // (0110)
      t[dirs[0]] += d[0]; oob +=    32*that->sub2map(t.data(),c[ 5u]); w[5] = p[0]*p[1]*p[2]*m[3]; // (1110)
      t[dirs[1]] -= d[1]; oob +=    64*that->sub2map(t.data(),c[ 6u]); w[6] = p[0]*m[1]*p[2]*m[3]; // (1010)
      t[dirs[0]] -= d[0]; oob +=   128*that->sub2map(t.data(),c[ 7u]); w[7] = m[0]*m[1]*p[2]*m[3]; // (0010)
      t[dirs[3]] += d[3]; oob +=   256*that->sub2map(t.data(),c[ 8u]); w[0] = m[0]*m[1]*p[2]*p[3]; // (0011)
      t[dirs[0]] += d[0]; oob +=   512*that->sub2map(t.data(),c[ 9u]); w[1] = p[0]*m[1]*p[2]*p[3]; // (1011)
      t[dirs[1]] += d[1]; oob +=  1024*that->sub2map(t.data(),c[10u]); w[2] = p[0]*p[1]*p[2]*p[3]; // (1111)
      t[dirs[0]] -= d[0]; oob +=  2048*that->sub2map(t.data(),c[11u]); w[3] = m[0]*p[1]*p[2]*p[3]; // (0111)
      t[dirs[2]] -= d[2]; oob +=  4096*that->sub2map(t.data(),c[12u]); w[4] = m[0]*p[1]*m[2]*p[3]; // (0101)
      t[dirs[0]] += d[0]; oob +=  8192*that->sub2map(t.data(),c[13u]); w[5] = p[0]*p[1]*m[2]*p[3]; // (1101)
      t[dirs[1]] -= d[1]; oob += 16384*that->sub2map(t.data(),c[14u]); w[6] = p[0]*m[1]*m[2]*p[3]; // (1001)
      t[dirs[0]] -= d[0]; oob += 32768*that->sub2map(t.data(),c[15u]); w[7] = m[0]*m[1]*m[2]*p[3]; // (0001)
      break;
      case 3:
                          oob +=       that->sub2map(t.data(),c[ 0u]); w[0] = m[0]*m[1]*m[2]; // (000)
      t[dirs[0]] += d[0]; oob +=     2*that->sub2map(t.data(),c[ 1u]); w[1] = p[0]*m[1]*m[2]; // (100)
      t[dirs[1]] += d[1]; oob +=     4*that->sub2map(t.data(),c[ 2u]); w[2] = p[0]*p[1]*m[2]; // (110)
      t[dirs[0]] -= d[0]; oob +=     8*that->sub2map(t.data(),c[ 3u]); w[3] = m[0]*p[1]*m[2]; // (010)
      t[dirs[2]] += d[2]; oob +=    16*that->sub2map(t.data(),c[ 4u]); w[4] = m[0]*p[1]*p[2]; // (011)
      t[dirs[0]] += d[0]; oob +=    32*that->sub2map(t.data(),c[ 5u]); w[5] = p[0]*p[1]*p[2]; // (111)
      t[dirs[1]] -= d[1]; oob +=    64*that->sub2map(t.data(),c[ 6u]); w[6] = p[0]*m[1]*p[2]; // (101)
      t[dirs[0]] -= d[0]; oob +=   128*that->sub2map(t.data(),c[ 7u]); w[7] = m[0]*m[1]*p[2]; // (001)
      break;
      case 2:
                          oob +=       that->sub2map(t.data(),c[ 0u]); w[0] = m[0]*m[1]; // (00)
      t[dirs[0]] += d[0]; oob +=     2*that->sub2map(t.data(),c[ 1u]); w[1] = p[0]*m[1]; // (10)
      t[dirs[1]] += d[1]; oob +=     4*that->sub2map(t.data(),c[ 2u]); w[2] = p[0]*p[1]; // (11)
      t[dirs[0]] -= d[0]; oob +=     8*that->sub2map(t.data(),c[ 3u]); w[3] = m[0]*p[1]; // (01)
      break;
      case 1:
                          oob +=       that->sub2map(t.data(),c[ 0u]); w[0] = m[0]; // (0)
      t[dirs[0]] += d[0]; oob +=     2*that->sub2map(t.data(),c[ 1u]); w[1] = p[0]; // (1)
      break;
      default:
      std::string msg = "Can not interpolate in " + std::to_string(ndims) + " dimensions.";
      throw std::runtime_error(msg);
    }
    return oob;
}

/*! \brief Find the linear index and weights of interpolation points for a given position

Linear interpolation in D dimensions requires 2ᴰ points surrounding the position
where an interpolated value is to be determined. The value at each of the 2ᴰ
points contributes in proportion to how close it is to the interpolation point.
This function finds the linear indices for any object which has a `sub2map`
method with signature `(const size_t* subscripted_index, size_t* linear_index)`
and determines the appropriate normalized weights for each of the 2ᴰ points.

@param that A pointer to an object with a `sub2map(const size_t*, size_t&)`
            method: MapGrid3, MapGrid4, InterpolateGrid3, InterpolateGrid4,
                    BrillouinZoneGrid3, BrillouinZoneGrid4, etc.
@param      zero  Pointer to the zero-point values of the grid, with N elements
@param      step  Pointer to the step-size values of the grid, with N elements
@param      ijk   Pointer to the nearest smaller grid subscripted index to the
                  interpolation point, with N elements
@param      x     Pointer to the coordinates of the interpolation point,
                  with N elements
@param[out] c     Pointer where the linear indices will be stored,
                  with 2ᴺ elements
@param[out] w     Pointer where the weights will be stored, with 2ᴺ elements
@param      N     Dimensionality of the gridded space
@param      dirs  Directions in the gridded space over which interpolation is
                  to be performed; D=1≤`dirs.size()`≤N
@returns An unsigned integer with bits acting as flags to indicate whether a
         point required for the interpolation is out-of-bounds. The order of
         bits depends on `dirs.size()`, see the source if you need to decode
         this information.
@note This function supports linear, bilinear, trilinear, and quadralinear
      interpolation, e.g., 1≤D≤4, in an arbitrary number of dimensions with D≤N;
      extending it to higher-dimensional interpolation should be possible
      if ever desired.
      For this function to work as anticipated the position must have the
      property
        zero + ijk*step <= x < zero+(ijk+1)*step
      for all indices in dirs.
*/
template <class R, template<class> class T>
int floor_corners_and_weights(const T<R>* that, const double* zero, const double* step, const size_t *ijk, const double *x, size_t *c, double *w, const size_t N, const std::vector<size_t>& dirs){
    size_t ndims = dirs.size();
    std::vector<double> p(ndims), m(ndims);
    for (size_t i=0; i<ndims; ++i){
      double tmp = zero[dirs[i]] + ((double)ijk[dirs[i]]) * step[dirs[i]];
      p[i] = (x[dirs[i]]-tmp)/step[dirs[i]];
      m[i] = 1.0 - p[i];
    }
    std::vector<size_t> t(N);
    for (size_t i=0; i<N; ++i) t[i] = ijk[i];
    int oob=0;
    switch (ndims){
      case 4:
                       oob +=       that->sub2map(t.data(),c[ 0u]); w[0] = m[0]*m[1]*m[2]*m[3]; // (0000)
      t[dirs[0]] += 1; oob +=     2*that->sub2map(t.data(),c[ 1u]); w[1] = p[0]*m[1]*m[2]*m[3]; // (1000)
      t[dirs[1]] += 1; oob +=     4*that->sub2map(t.data(),c[ 2u]); w[2] = p[0]*p[1]*m[2]*m[3]; // (1100)
      t[dirs[0]] -= 1; oob +=     8*that->sub2map(t.data(),c[ 3u]); w[3] = m[0]*p[1]*m[2]*m[3]; // (0100)
      t[dirs[2]] += 1; oob +=    16*that->sub2map(t.data(),c[ 4u]); w[4] = m[0]*p[1]*p[2]*m[3]; // (0110)
      t[dirs[0]] += 1; oob +=    32*that->sub2map(t.data(),c[ 5u]); w[5] = p[0]*p[1]*p[2]*m[3]; // (1110)
      t[dirs[1]] -= 1; oob +=    64*that->sub2map(t.data(),c[ 6u]); w[6] = p[0]*m[1]*p[2]*m[3]; // (1010)
      t[dirs[0]] -= 1; oob +=   128*that->sub2map(t.data(),c[ 7u]); w[7] = m[0]*m[1]*p[2]*m[3]; // (0010)
      t[dirs[3]] += 1; oob +=   256*that->sub2map(t.data(),c[ 8u]); w[0] = m[0]*m[1]*p[2]*p[3]; // (0011)
      t[dirs[0]] += 1; oob +=   512*that->sub2map(t.data(),c[ 9u]); w[1] = p[0]*m[1]*p[2]*p[3]; // (1011)
      t[dirs[1]] += 1; oob +=  1024*that->sub2map(t.data(),c[10u]); w[2] = p[0]*p[1]*p[2]*p[3]; // (1111)
      t[dirs[0]] -= 1; oob +=  2048*that->sub2map(t.data(),c[11u]); w[3] = m[0]*p[1]*p[2]*p[3]; // (0111)
      t[dirs[2]] -= 1; oob +=  4096*that->sub2map(t.data(),c[12u]); w[4] = m[0]*p[1]*m[2]*p[3]; // (0101)
      t[dirs[0]] += 1; oob +=  8192*that->sub2map(t.data(),c[13u]); w[5] = p[0]*p[1]*m[2]*p[3]; // (1101)
      t[dirs[1]] -= 1; oob += 16384*that->sub2map(t.data(),c[14u]); w[6] = p[0]*m[1]*m[2]*p[3]; // (1001)
      t[dirs[0]] -= 1; oob += 32768*that->sub2map(t.data(),c[15u]); w[7] = m[0]*m[1]*m[2]*p[3]; // (0001)
      break;
      case 3:
                       oob +=       that->sub2map(t.data(),c[ 0u]); w[0] = m[0]*m[1]*m[2]; // (000)
      t[dirs[0]] += 1; oob +=     2*that->sub2map(t.data(),c[ 1u]); w[1] = p[0]*m[1]*m[2]; // (100)
      t[dirs[1]] += 1; oob +=     4*that->sub2map(t.data(),c[ 2u]); w[2] = p[0]*p[1]*m[2]; // (110)
      t[dirs[0]] -= 1; oob +=     8*that->sub2map(t.data(),c[ 3u]); w[3] = m[0]*p[1]*m[2]; // (010)
      t[dirs[2]] += 1; oob +=    16*that->sub2map(t.data(),c[ 4u]); w[4] = m[0]*p[1]*p[2]; // (011)
      t[dirs[0]] += 1; oob +=    32*that->sub2map(t.data(),c[ 5u]); w[5] = p[0]*p[1]*p[2]; // (111)
      t[dirs[1]] -= 1; oob +=    64*that->sub2map(t.data(),c[ 6u]); w[6] = p[0]*m[1]*p[2]; // (101)
      t[dirs[0]] -= 1; oob +=   128*that->sub2map(t.data(),c[ 7u]); w[7] = m[0]*m[1]*p[2]; // (001)
      break;
      case 2:
                       oob +=       that->sub2map(t.data(),c[ 0u]); w[0] = m[0]*m[1]; // (00)
      t[dirs[0]] += 1; oob +=     2*that->sub2map(t.data(),c[ 1u]); w[1] = p[0]*m[1]; // (10)
      t[dirs[1]] += 1; oob +=     4*that->sub2map(t.data(),c[ 2u]); w[2] = p[0]*p[1]; // (11)
      t[dirs[0]] -= 1; oob +=     8*that->sub2map(t.data(),c[ 3u]); w[3] = m[0]*p[1]; // (01)
      break;
      case 1:
                       oob +=       that->sub2map(t.data(),c[ 0u]); w[0] = m[0]; // (0)
      t[dirs[0]] += 1; oob +=     2*that->sub2map(t.data(),c[ 1u]); w[1] = p[0]; // (1)
      break;
      default:
      std::string msg = "Can not interpolate in " + std::to_string(ndims) + " dimensions.";
      throw std::runtime_error(msg);
    }
    return oob;
}

template<class T, template<class> class A>
T
triangle_area(const A<T>& a, const A<T>& b, const A<T>& c){
  T ab, bc, ac, s;
  ab = (a-b).norm(0);
  bc = (b-c).norm(0);
  ac = (a-c).norm(0);
  s = (ab+bc+ac)/2.0;
  return std::sqrt(s*(s-ab)*(s-bc)*(s-ac)); // Heron's formula
}
template<class T, template<class> class A>
T
tetrahedron_volume(const A<T>& a, const A<T>& b, const A<T>& c, const A<T>& d){
  A<T> dumb({4,3});
  dumb.set(0, b-a);
  dumb.set(1, c-a);
  dumb.set(2, d-a);
  brille::utils::vector_cross(dumb.ptr(3), dumb.ptr(1), dumb.ptr(2));
  return std::abs(dumb.dot(0,3)/6.0); // abs so that we don't have to worry about the permutation
}
/*
  For between 1 and 4 vertices of a tetrahedron, v, determine the linear
  inerpolation weight with with each vertex contributes at a point, p.
  The required relationship between p and the vertices depends on the number
  of supplied vertices:
    v.size()  requisite relationship
    --------  -----------------------------------------------------------------
       1        p-v ≡ ⃗0
       2        p must lie on the line between v[0] and v[1]
       3        p must lie within the triangle formed by v[0], v[1], and v[2]
       4        p must lie within the volume of the tetrahedron
*/
template<class T, template<class> class A>
std::vector<double>
tetrahedron_weights(const A<T>& p, const A<T>& v){
  std::vector<double> weights(v.size(0));
  switch (v.size(0)){
    case 1:{
      info_update("interpolation weight for one point is 1.");
      weights[0] = 1.0;
      break;
    }
    case 2:{
      info_update("finding interpolation weights for point along a line.");
      T len = (v.view(1)-v.view(0)).norm(0);
      weights[0] = (v.view(1)-p).norm(0)/len;
      weights[1] = (v.view(0)-p).norm(0)/len;
      break;
    }
    case 3:{
      info_update("finding interpolation weights for point in a triangle.");
      T area = triangle_area(v.view(0), v.view(1), v.view(2));
      weights[0] = triangle_area(p, v.view(1), v.view(2))/area;
      weights[1] = triangle_area(v.view(0), p, v.view(2))/area;
      weights[2] = triangle_area(v.view(0), v.view(1), p)/area;
      break;
    }
    case 4:{
      info_update("finding interpolation weights for point in a tetrahedron.");
      T vol = tetrahedron_volume(v.view(0), v.view(1), v.view(2), v.view(3));
      weights[0] = tetrahedron_volume(p, v.view(1), v.view(2), v.view(3))/vol;
      weights[1] = tetrahedron_volume(v.view(0), p, v.view(2), v.view(3))/vol;
      weights[2] = tetrahedron_volume(v.view(0), v.view(1), p, v.view(3))/vol;
      weights[3] = tetrahedron_volume(v.view(0), v.view(1), v.view(2), p)/vol;
      break;
    }
    default:
      throw std::runtime_error("interpolation.tetrahedron_weights: this should be impossible.");
  }
  return weights;
}

#endif
