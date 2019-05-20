#ifndef _INTERPOLATION_H_
#define _INTERPOLATION_H_
#include <iostream>
// #include <functional>
#include <vector>

double interpolation_direction_and_distance(const double, const double, const size_t, const double);

template <class R, template<class> class T>
int corners_and_weights( T<R>* that, const double* zero, const double* step, const size_t *ijk, const double *x, size_t *c, double *w, const size_t N, const std::vector<size_t> dirs){
    size_t ndims = dirs.size();
    std::vector<double> p(ndims), m(ndims);
    std::vector<int> d(ndims);
    size_t ti;
    for (size_t i=0; i<ndims; ++i){
      p[i] = interpolation_direction_and_distance(zero[dirs[i]],step[dirs[i]],ijk[dirs[i]],x[dirs[i]]);
      d[i] = p[i] < 0 ? -1 : 1;
      p[i] = std::abs(p[i]);
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

#endif
