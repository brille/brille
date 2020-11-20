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

#include "neighbours.hpp"
#include <array>
#include <vector>

using namespace brille;

static std::vector<int> all_from_extent(const int extent){
  int min = -extent, max = extent+1;
  std::vector<int> vec;
  for (int i=min; i<max; ++i) vec.push_back(i);
  return vec;
}

bArray<int> brille::make_relative_neighbour_indices(const int extent){
  auto vec = all_from_extent(extent);
  size_t num = vec.size();
  size_t n{0};
  std::vector<std::array<int,3>> o(num*num*num-1);  //-1 because we're skipping the origin)
  for (auto i: vec) for (auto j: vec) for (auto k: vec) if (i||j||k)
    o[n++] = {{i,j,k}};
  return bArray<int>::from_std(o);
}

bArray<int> brille::make_relative_neighbour_indices_prime(const int extent){
  auto vec = all_from_extent(extent);
  size_t num = vec.size();
  size_t n{0};
  std::vector<std::array<int,3>> o(num*num*num-1);
  for (size_t i=0; i<num*num*num-1; ++i) o[i] = {{0,0,0}};
  // first follow each axis:
  for (int e=0; e<3; ++e) for (auto j: vec) if (j) o[n++][e] = j;
  // next the ab, bc, and ac planes:
  for (int e0=0; e0<2; ++e0) for (int e1=e0+1; e1<3; ++e1)
  for (auto k: vec) for (auto j: vec) if (k&&j){
    o[n][e0] = k; o[n++][e1] = j;
  }
  // then fill in everything else
  for (auto i: vec) for (auto j: vec) for (auto k: vec) if (i&&j&&k)
    o[n++] = {{i,j,k}};
  return bArray<int>::from_std(o);
}

bArray<int> brille::make_relative_neighbour_indices4(const int extent){
  auto v = all_from_extent(extent);
  std::vector<std::array<int,4>> o;
  for (auto i: v) for (auto j: v) for (auto k: v) for (auto l: v)
  if (i||j||k||l) o.push_back({{i,j,k,l}});
  return bArray<int>::from_std(o);
}
