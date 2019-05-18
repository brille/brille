#include "interpolation.h"

double interpolation_direction_and_distance(const double zero, const double step, const size_t nrstidx, const double x){
  double tmp = x - (zero+nrstidx*step);
  double d_p = tmp/(double)(step);
  return d_p;
}
