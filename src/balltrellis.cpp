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

#include "balltrellis.hpp"
#include "debug.hpp"

using namespace brille;

/* A note about the directions chosen for binning the leaves:
  A better way of doing this would be to compute the covariance and/or
  correlation matrix for the leaf centres, e.g.,
    | ∑|xᵢ-̄x|² ∑|yᵢ-̄x|² ∑|zᵢ-̄x|² |
    | ∑|xᵢ-̄y|² ∑|yᵢ-̄y|² ∑|zᵢ-̄y|² |
    | ∑|xᵢ-̄z|² ∑|yᵢ-̄z|² ∑|zᵢ-̄z|² |
  and then find its eigenvectors. If we ever decide to include a real
  linear algebra library then we could easily update this function to make
  use of its eigenvalue solver as well.
  For now, an additional library just to address this issue is probably overkill.
*/
Trellis construct_trellis(const std::vector<TrellisLeaf>& leaves, const double fraction){
  std::array<double,3> p, c{0.,0.,0.};
  std::array<double,9> xyz{1.,0.,0., 0.,1.,0., 0.,0.,1.};
  double d;
  // we need more than one leaf to be able to calculate x, y, and z
  if (leaves.size()>1){
    // find the centroid of the leaf positions
    for (size_t i=0; i<3u; ++i) c[i] = 0.;
    for (auto leaf: leaves){
      p = leaf.centre();
      for (size_t i=0; i<3u; ++i) c[i] += p[i];
    }
    for (size_t i=0; i<3u; ++i) c[i] /= static_cast<double>(leaves.size());
    // plus the vector pointing to the farthest leaf from the centroid
    // choose this direction to be our x-axis
    double dmax = 0.;
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
      for (size_t i=0; i<3u; ++i) d += p[i]*xyz[i]; // dot(p, ̂x)
      for (size_t i=0; i<3u; ++i) pperp[i] = p[i] - d*xyz[i]; // p - ̂x dot(p, ̂x)
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
    brille::utils::vector_cross(xyz.data()+6, xyz.data(), xyz.data()+3);
  }
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
  // >>>>>>>>>>>>>>>>>>>>>>>>>>> optimisation testing
  double max_radius = 0.;
  for (auto leaf: leaves) if (leaf.radius() > max_radius) max_radius = leaf.radius();
  debug_update("The maximum leaf radius is ",max_radius);
  debug_update("And the leaves extend over");
  debug_update(" x = ",xyz[0]," ",xyz[1]," ",xyz[2]," : ",minmax[0]);
  debug_update(" y = ",xyz[3]," ",xyz[4]," ",xyz[5]," : ",minmax[1]);
  debug_update(" z = ",xyz[6]," ",xyz[7]," ",xyz[8]," : ",minmax[2]);
  d = max_radius * fraction;
  debug_update("Using a boundary step size of max_radius*",fraction," which is ",d);
  // <<<<<<<<<<<<<<<<<<<<<<<<<<<
  // and construct the boundaries
  std::array<std::vector<double>,3> boundaries;
  for (size_t i=0; i<3u; ++i){
    // // find the boundary step size
    // d = (minmax[i][1]-minmax[i][0]) / static_cast<double>(Nxyz[i]);
    // and construct the boundaries:
    boundaries[i].push_back(minmax[i][0]);
    while( boundaries[i].back() < minmax[i][1] ) boundaries[i].push_back(boundaries[i].back()+d);
    // replace the first and last boundaries to -∞ and +∞, respectively
    // boundaries[i].front() = -std::numeric_limits<double>::infinity();
    // boundaries[i].back()  =  std::numeric_limits<double>::infinity();
    debug_update("Step size along axis ",i,", ",d,", yields ",boundaries[i].size()-1," bins with boundaries ",boundaries[i]);
  }

  // Construct the Trellis object assigning the leaves to nodes
  Trellis trellis(xyz, boundaries, leaves);
  return trellis;
}
