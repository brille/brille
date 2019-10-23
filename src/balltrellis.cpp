#include "balltrellis.h"
#include "debug.h"

Trellis construct_trellis(const std::vector<TrellisLeaf>& leaves, const size_t Nxyz){
  std::array<size_t,3> Nxyz_array{Nxyz, Nxyz, Nxyz};
  return construct_trellis(leaves, Nxyz_array);
}
Trellis construct_trellis(const std::vector<TrellisLeaf>& leaves, const std::array<size_t,3> Nxyz){
  std::array<double,3> p, c{0.,0.,0.};
  std::array<double,9> xyz;
  // find the centroid of the leaf positions
  for (size_t i=0; i<3u; ++i) c[i] = 0.;
  for (auto leaf: leaves){
    p = leaf.centre();
    for (size_t i=0; i<3u; ++i) c[i] += p[i];
  }
  for (size_t i=0; i<3u; ++i) c[i] /= static_cast<double>(leaves.size());
  // plus the vector pointing to the farthest leaf from the centroid
  // choose this direction to be our x-axis
  double d, dmax = 0.;
  for (auto leaf: leaves){
    p = leaf.centre();
    d = 0.;
    for (size_t i=0; i<3u; ++i) d += (p[i]-c[i])*(p[i]-c[i]);
    if (d>dmax){
      dmax = d;
      for (size_t i=0; i<3u; ++i) xyz[i] = p[i] - c[i];
    }
  }
  // normalise the x-axis:
  d = 0.;
  for (size_t i=0; i<3u; ++i) d += xyz[i]*xyz[i];
  d = std::sqrt(d);
  for (size_t i=0; i<3u; ++i) xyz[i] /= d;
  // find the perpendicular direction with the farthest leaf from the centroid
  // choose this to be our y-axis
  dmax = 0.;
  std::array<double,3> pperp;
  for (auto leaf: leaves){
    p = leaf.centre();
    d = 0.; // p⋅x
    for (size_t i=0; i<3u; ++i) d += (p[i]-c[i])*xyz[i];
    for (size_t i=0; i<3u; ++i) pperp[i] = (p[i]-c[i])*(1-d*xyz[i]);
    d = 0.;
    for (size_t i=0; i<3u; ++i) d += pperp[i]*pperp[i];
    if (d > dmax){
      dmax = d;
      for (size_t i=0; i<3u; ++i) xyz[i+3] = pperp[i];
    }
  }
  // normalise the y-axis;
  d = 0;
  for (size_t i=0; i<3u; ++i) d += xyz[i+3]*xyz[i+3];
  d = std::sqrt(d);
  for (size_t i=0; i<3u; ++i) xyz[i+3] /= d;
  // pick the z-axis as x × y
  vector_cross(xyz.data()+6, xyz.data(), xyz.data()+3);
  // find the extents of the data
  std::array<std::array<double,2>,3> minmax;
  for (size_t i=0; i<3u; ++i){
    minmax[i][0] = (std::numeric_limits<double>::max)();
    minmax[i][1] = std::numeric_limits<double>::lowest();
  }
  // minimum does not need to include the radius but maximum does:
  for (auto leaf: leaves){
    p = leaf.centre();
    double r = leaf.radius();
    for (size_t i=0; i<3u; ++i){
      d = 0.;
      for (size_t j=0; j<3u; ++j) d += p[j]*xyz[i*3u+j];
      if (d - r < minmax[i][0]) minmax[i][0] = d - r;
      if (d + r > minmax[i][1]) minmax[i][1] = d + r;
    }
  }
  // and construct the boundaries
  std::array<std::vector<double>,3> boundaries;
  for (size_t i=0; i<3u; ++i){
    // find the boundary step size
    d = (minmax[i][1]-minmax[i][0]) / static_cast<double>(Nxyz[i]);
    // and construct the boundaries:
    boundaries[i].push_back(minmax[i][0]);
    while( boundaries[i].back() < minmax[i][1] ) boundaries[i].push_back(boundaries[i].back()+d);
    // replace the first and last boundaries to -∞ and +∞, respectively
    boundaries[i].front() = std::numeric_limits<double>::lowest();
    boundaries[i].back() = (std::numeric_limits<double>::max)();
  }

  // Construct the Trellis object assigning the leaves to nodes
  Trellis trellis(xyz, boundaries, leaves);
  return trellis;
}
