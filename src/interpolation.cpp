#include "interpolation.h"
// #include <functional>
#include <vector>
// using namespace std::placeholders;

double interpolation_direction_and_distance(const double zero, const double step, const size_t nrstidx, const double x){
  double tmp = x - (zero+nrstidx*step);
  double d_p = tmp/(double)(step);
  return d_p;
}

// template <int N>
// int corners_and_weights(std::function<int(const size_t*,size_t*)> crnrfun, const double* zero, const double* step, const size_t *ijk, const double *x, size_t *c, double *w, const std::vector<size_t> dirs){
