/*! \file */
#ifndef _INTERPOLATION_H_
#define _INTERPOLATION_H_
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

@param that A pointer to an object with a `sub2map(const size_t*, size_t*)`
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
int corners_and_weights( T<R>* that, const double* zero, const double* step, const size_t *ijk, const double *x, size_t *c, double *w, const size_t N, const std::vector<size_t> dirs){
    size_t ndims = dirs.size();
    std::vector<double> p(ndims), m(ndims);
    std::vector<int> d(ndims);
    size_t ti;
    double tmp;
    for (size_t i=0; i<ndims; ++i){
      // tmp = interpolation_direction_and_distance(zero[dirs[i]],step[dirs[i]],ijk[dirs[i]],x[dirs[i]]);
      tmp = (x[dirs[i]]-(zero[dirs[i]]+ijk[dirs[i]]*step[dirs[i]]))/step[dirs[i]];
      d[i] = tmp < 0 ? -1 : 1;
      p[i] = std::abs(tmp);
      m[i] = 1.0 - p[i];
    }
    std::vector<size_t> t(N);
    for (int i=0; i<N; ++i) t[i] = ijk[i];
    int oob=0;
    switch (ndims){
      case 4:
                          oob +=       that->sub2map(t.data(),c    ); w[0] = m[0]*m[1]*m[2]*m[3]; // (0000)
      t[dirs[0]] += d[0]; oob +=     2*that->sub2map(t.data(),c+ 1u); w[1] = p[0]*m[1]*m[2]*m[3]; // (1000)
      t[dirs[1]] += d[1]; oob +=     4*that->sub2map(t.data(),c+ 2u); w[2] = p[0]*p[1]*m[2]*m[3]; // (1100)
      t[dirs[0]] -= d[0]; oob +=     8*that->sub2map(t.data(),c+ 3u); w[3] = m[0]*p[1]*m[2]*m[3]; // (0100)
      t[dirs[2]] += d[2]; oob +=    16*that->sub2map(t.data(),c+ 4u); w[4] = m[0]*p[1]*p[2]*m[3]; // (0110)
      t[dirs[0]] += d[0]; oob +=    32*that->sub2map(t.data(),c+ 5u); w[5] = p[0]*p[1]*p[2]*m[3]; // (1110)
      t[dirs[1]] -= d[1]; oob +=    64*that->sub2map(t.data(),c+ 6u); w[6] = p[0]*m[1]*p[2]*m[3]; // (1010)
      t[dirs[0]] -= d[0]; oob +=   128*that->sub2map(t.data(),c+ 7u); w[7] = m[0]*m[1]*p[2]*m[3]; // (0010)
      t[dirs[3]] += d[3]; oob +=   256*that->sub2map(t.data(),c+ 8u); w[0] = m[0]*m[1]*p[2]*p[3]; // (0011)
      t[dirs[0]] += d[0]; oob +=   512*that->sub2map(t.data(),c+ 9u); w[1] = p[0]*m[1]*p[2]*p[3]; // (1011)
      t[dirs[1]] += d[1]; oob +=  1024*that->sub2map(t.data(),c+10u); w[2] = p[0]*p[1]*p[2]*p[3]; // (1111)
      t[dirs[0]] -= d[0]; oob +=  2048*that->sub2map(t.data(),c+11u); w[3] = m[0]*p[1]*p[2]*p[3]; // (0111)
      t[dirs[2]] -= d[2]; oob +=  4096*that->sub2map(t.data(),c+12u); w[4] = m[0]*p[1]*m[2]*p[3]; // (0101)
      t[dirs[0]] += d[0]; oob +=  8192*that->sub2map(t.data(),c+13u); w[5] = p[0]*p[1]*m[2]*p[3]; // (1101)
      t[dirs[1]] -= d[1]; oob += 16384*that->sub2map(t.data(),c+14u); w[6] = p[0]*m[1]*m[2]*p[3]; // (1001)
      t[dirs[0]] -= d[0]; oob += 32768*that->sub2map(t.data(),c+15u); w[7] = m[0]*m[1]*m[2]*p[3]; // (0001)
      break;
      case 3:
                          oob +=       that->sub2map(t.data(),c    ); w[0] = m[0]*m[1]*m[2]; // (000)
      t[dirs[0]] += d[0]; oob +=     2*that->sub2map(t.data(),c+ 1u); w[1] = p[0]*m[1]*m[2]; // (100)
      t[dirs[1]] += d[1]; oob +=     4*that->sub2map(t.data(),c+ 2u); w[2] = p[0]*p[1]*m[2]; // (110)
      t[dirs[0]] -= d[0]; oob +=     8*that->sub2map(t.data(),c+ 3u); w[3] = m[0]*p[1]*m[2]; // (010)
      t[dirs[2]] += d[2]; oob +=    16*that->sub2map(t.data(),c+ 4u); w[4] = m[0]*p[1]*p[2]; // (011)
      t[dirs[0]] += d[0]; oob +=    32*that->sub2map(t.data(),c+ 5u); w[5] = p[0]*p[1]*p[2]; // (111)
      t[dirs[1]] -= d[1]; oob +=    64*that->sub2map(t.data(),c+ 6u); w[6] = p[0]*m[1]*p[2]; // (101)
      t[dirs[0]] -= d[0]; oob +=   128*that->sub2map(t.data(),c+ 7u); w[7] = m[0]*m[1]*p[2]; // (001)
      break;
      case 2:
                          oob +=       that->sub2map(t.data(),c    ); w[0] = m[0]*m[1]; // (00)
      t[dirs[0]] += d[0]; oob +=     2*that->sub2map(t.data(),c+ 1u); w[1] = p[0]*m[1]; // (10)
      t[dirs[1]] += d[1]; oob +=     4*that->sub2map(t.data(),c+ 2u); w[2] = p[0]*p[1]; // (11)
      t[dirs[0]] -= d[0]; oob +=     8*that->sub2map(t.data(),c+ 3u); w[3] = m[0]*p[1]; // (01)
      break;
      case 1:
                          oob +=       that->sub2map(t.data(),c    ); w[0] = m[0]; // (0)
      t[dirs[0]] += d[0]; oob +=     2*that->sub2map(t.data(),c+ 1u); w[1] = p[0]; // (1)
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

@param that A pointer to an object with a `sub2map(const size_t*, size_t*)`
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
int floor_corners_and_weights( T<R>* that, const double* zero, const double* step, const size_t *ijk, const double *x, size_t *c, double *w, const size_t N, const std::vector<size_t> dirs){
    size_t ndims = dirs.size();
    std::vector<double> p(ndims), m(ndims);
    size_t ti;
    double tmp;
    for (size_t i=0; i<ndims; ++i){
      tmp = zero[dirs[i]] + ((double)ijk[dirs[i]]) * step[dirs[i]];
      p[i] = (x[dirs[i]]-tmp)/step[dirs[i]];
      m[i] = 1.0 - p[i];
    }
    std::vector<size_t> t(N);
    for (int i=0; i<N; ++i) t[i] = ijk[i];
    int oob=0;
    switch (ndims){
      case 4:
                       oob +=       that->sub2map(t.data(),c    ); w[0] = m[0]*m[1]*m[2]*m[3]; // (0000)
      t[dirs[0]] += 1; oob +=     2*that->sub2map(t.data(),c+ 1u); w[1] = p[0]*m[1]*m[2]*m[3]; // (1000)
      t[dirs[1]] += 1; oob +=     4*that->sub2map(t.data(),c+ 2u); w[2] = p[0]*p[1]*m[2]*m[3]; // (1100)
      t[dirs[0]] -= 1; oob +=     8*that->sub2map(t.data(),c+ 3u); w[3] = m[0]*p[1]*m[2]*m[3]; // (0100)
      t[dirs[2]] += 1; oob +=    16*that->sub2map(t.data(),c+ 4u); w[4] = m[0]*p[1]*p[2]*m[3]; // (0110)
      t[dirs[0]] += 1; oob +=    32*that->sub2map(t.data(),c+ 5u); w[5] = p[0]*p[1]*p[2]*m[3]; // (1110)
      t[dirs[1]] -= 1; oob +=    64*that->sub2map(t.data(),c+ 6u); w[6] = p[0]*m[1]*p[2]*m[3]; // (1010)
      t[dirs[0]] -= 1; oob +=   128*that->sub2map(t.data(),c+ 7u); w[7] = m[0]*m[1]*p[2]*m[3]; // (0010)
      t[dirs[3]] += 1; oob +=   256*that->sub2map(t.data(),c+ 8u); w[0] = m[0]*m[1]*p[2]*p[3]; // (0011)
      t[dirs[0]] += 1; oob +=   512*that->sub2map(t.data(),c+ 9u); w[1] = p[0]*m[1]*p[2]*p[3]; // (1011)
      t[dirs[1]] += 1; oob +=  1024*that->sub2map(t.data(),c+10u); w[2] = p[0]*p[1]*p[2]*p[3]; // (1111)
      t[dirs[0]] -= 1; oob +=  2048*that->sub2map(t.data(),c+11u); w[3] = m[0]*p[1]*p[2]*p[3]; // (0111)
      t[dirs[2]] -= 1; oob +=  4096*that->sub2map(t.data(),c+12u); w[4] = m[0]*p[1]*m[2]*p[3]; // (0101)
      t[dirs[0]] += 1; oob +=  8192*that->sub2map(t.data(),c+13u); w[5] = p[0]*p[1]*m[2]*p[3]; // (1101)
      t[dirs[1]] -= 1; oob += 16384*that->sub2map(t.data(),c+14u); w[6] = p[0]*m[1]*m[2]*p[3]; // (1001)
      t[dirs[0]] -= 1; oob += 32768*that->sub2map(t.data(),c+15u); w[7] = m[0]*m[1]*m[2]*p[3]; // (0001)
      break;
      case 3:
                       oob +=       that->sub2map(t.data(),c    ); w[0] = m[0]*m[1]*m[2]; // (000)
      t[dirs[0]] += 1; oob +=     2*that->sub2map(t.data(),c+ 1u); w[1] = p[0]*m[1]*m[2]; // (100)
      t[dirs[1]] += 1; oob +=     4*that->sub2map(t.data(),c+ 2u); w[2] = p[0]*p[1]*m[2]; // (110)
      t[dirs[0]] -= 1; oob +=     8*that->sub2map(t.data(),c+ 3u); w[3] = m[0]*p[1]*m[2]; // (010)
      t[dirs[2]] += 1; oob +=    16*that->sub2map(t.data(),c+ 4u); w[4] = m[0]*p[1]*p[2]; // (011)
      t[dirs[0]] += 1; oob +=    32*that->sub2map(t.data(),c+ 5u); w[5] = p[0]*p[1]*p[2]; // (111)
      t[dirs[1]] -= 1; oob +=    64*that->sub2map(t.data(),c+ 6u); w[6] = p[0]*m[1]*p[2]; // (101)
      t[dirs[0]] -= 1; oob +=   128*that->sub2map(t.data(),c+ 7u); w[7] = m[0]*m[1]*p[2]; // (001)
      break;
      case 2:
                       oob +=       that->sub2map(t.data(),c    ); w[0] = m[0]*m[1]; // (00)
      t[dirs[0]] += 1; oob +=     2*that->sub2map(t.data(),c+ 1u); w[1] = p[0]*m[1]; // (10)
      t[dirs[1]] += 1; oob +=     4*that->sub2map(t.data(),c+ 2u); w[2] = p[0]*p[1]; // (11)
      t[dirs[0]] -= 1; oob +=     8*that->sub2map(t.data(),c+ 3u); w[3] = m[0]*p[1]; // (01)
      break;
      case 1:
                       oob +=       that->sub2map(t.data(),c    ); w[0] = m[0]; // (0)
      t[dirs[0]] += 1; oob +=     2*that->sub2map(t.data(),c+ 1u); w[1] = p[0]; // (1)
      break;
      default:
      std::string msg = "Can not interpolate in " + std::to_string(ndims) + " dimensions.";
      throw std::runtime_error(msg);
    }
    return oob;
}

#endif
