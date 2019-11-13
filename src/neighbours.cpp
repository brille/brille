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

#include "neighbours.hpp"

// int make_all_indices(LQVec<int> *ijk, const int extent){
//   if (ijk->numel()!=3u) throw "expected three vectors";
//   int min = -extent, max = extent+1;
//   int num = max-min;
//   ijk->resize(num*num*num-1); //-1 because we're skipping the origin
//   int *tmp = new int[3];
//   int n=0;
//   for (int i=min; i<max; i++){
//     tmp[0] = i;
//     for (int j=min; j<max; j++){
//       tmp[1] = j;
//       for (int k=min; k<max; k++){
//         tmp[2] = k;
//         if (i||j||k) ijk->set(n++,tmp); // as long as any are non-zero
//       }
//     }
//   }
//   delete[] tmp;
//   return n;
// }

ArrayVector<int> make_relative_neighbour_indices(const int extent){
  int min = -extent, max = extent+1;
  int num = max-min;
  ArrayVector<int> out(3u,num*num*num-1);  //-1 because we're skipping the origin
  int *tmp = new int[3];
  int n=0;
  for (int i=min; i<max; i++){
    tmp[0] = i;
    for (int j=min; j<max; j++){
      tmp[1] = j;
      for (int k=min; k<max; k++){
        tmp[2] = k;
        if (i||j||k) out.set(n++,tmp); // as long as any are non-zero
      }
    }
  }
  delete[] tmp;
  return out;
}

ArrayVector<int> make_relative_neighbour_indices_prime(const int extent){
  int min = -extent, max = extent+1;
  int num = max-min;
  ArrayVector<int> out(3u,num*num*num-1);  //-1 because we're skipping the origin
  int *tmp = nullptr;
  tmp = new int[3];
  int n=0;
  for (int e=0; e<3; ++e){
    for (int j=min; j<max; ++j){
      for (int i=0; i<3; ++i) tmp[i]=0;
      tmp[e] = j;
      if (j) out.set(n++, tmp);
    }
  }
  for (int e0=0; e0<2; ++e0){
    for (int e1=e0+1; e1<3; ++e1){
      for (int k=min; k<max; ++k){
        for (int j=min; j<max; ++j){
          for (int i=0; i<3; ++i) tmp[i]=0;
          tmp[e0]=k;
          tmp[e1]=j;
          if (k&&j) out.set(n++,tmp);
        }
      }
    }
  }
  for (int i=min; i<max; i++){
    tmp[0] = i;
    for (int j=min; j<max; j++){
      tmp[1] = j;
      for (int k=min; k<max; k++){
        tmp[2] = k;
        if (i&&j&&k) out.set(n++,tmp); // as long as any are non-zero
      }
    }
  }
  delete[] tmp;
  return out;
}

ArrayVector<int> make_relative_neighbour_indices4(const int extent){
  int min = -extent, max = extent+1;
  int num = max-min;
  ArrayVector<int> out(4u,num*num*num*num-1);  //-1 because we're skipping the origin
  int *tmp = nullptr;
  tmp = new int[4];
  int n=0;
  for (int i=min; i<max; i++){
    tmp[0] = i;
    for (int j=min; j<max; j++){
      tmp[1] = j;
      for (int k=min; k<max; k++){
        tmp[2] = k;
        for (int l=min; l<max; l++){
          tmp[3] = l;
          if (i||j||k||l) out.set(n++,tmp); // as long as any are non-zero
        }
      }
    }
  }
  delete[] tmp;
  return out;
}
