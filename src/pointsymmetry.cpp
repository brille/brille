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
#include <cstring>
#include <string>
#include <sstream>
#include "pointsymmetry.hpp"
#include <iostream>
#include "pointgroup.hpp"
#include <algorithm>
#include "approx.hpp"
using namespace brille;
/*****************************************************************************\
| PointSymmetry class Member functions:                                       |
\*****************************************************************************/
size_t PointSymmetry::add(const int *r){
  Matrix<int> newR;
  for (size_t i=0; i<9; ++i) newR[i] = r[i];
  return this->add(newR);
}
size_t PointSymmetry::add(const Matrix<int>& r){
  this->R.push_back(r);
  return this->R.size();
}
size_t PointSymmetry::add(const std::string& s){
  Matrix<int> newR{0,0,0,0,0,0,0,0,0};
  // a valid encoded Motion will be of the form
  // ∑ₙiₙαₙ + f where each iₙ is an integer, αₙ∈{x,y,z}
  // and f one of 0, ±1/12, ±1/6, ±1/4, ±1/3, ±1/2
  // but we should support any rational or floating point value for f
  int n{0}, i{1};
  char c{' '};
  std::istringstream stream(s);
  size_t ret{0};
  while (stream.good()){
    c = char(stream.get());
    switch (c){
      case ';': // end of one motion
        ret = this->add(newR);
        newR = {0,0,0,0,0,0,0,0,0};
        n=0;
        i=1;
        break;
      case ',': // end of motion component
        i=1;
        break;
      case 'x':
        newR[n*3] = i;
        i=1;
        break;
      case 'y':
        newR[n*3+1] = i;
        i=1;
        break;
      case 'z':
        newR[n*3+2] = i;
        i=1;
        break;
      default: // not ;,xyz
        std::stringstream tmp;
        tmp << c;
        bool haspoint{false}, hasslash{false};
        while (std::strchr(";,xyz-+",stream.peek())==nullptr){
          hasslash |= '/'==stream.peek();
          haspoint |= '.'==stream.peek();
          tmp << stream.get();
        }
        if (haspoint && hasslash) throw std::logic_error("Number has both point and slash?!");
        if (!(haspoint^hasslash)) tmp >> i; // TT (error) or FF (integer)
        // either floating point or rational numbers represent part of the
        // translation which strictly shouldn't be here for pointgroup
        // operations, but we should just skip for now
    }
  }
  if (';'!=c) ret = this->add(newR);
  return ret;
}
bool PointSymmetry::set(const size_t i, const int *r){
  if ( i>=this->size() ) return false;
  for(size_t j=0; j<9; j++) this->R[i][j] = r[j];
  return true;
}
bool PointSymmetry::get(const size_t i, int *r) const {
  if ( i>=this->size() ) return false;
  for(size_t j=0; j<9; j++) r[j]=this->R[i][j];
  return true;
}
int * PointSymmetry::data(const size_t i) {
  return (i<this->size()) ? this->R[i].data() : nullptr;
}
const int * PointSymmetry::data(const size_t i) const {
  return (i<this->size()) ? this->R[i].data() : nullptr;
}
Matrix<int> PointSymmetry::get(const size_t i) const {
  return (i<this->size()) ? this->R[i] : Matrix<int>({0,0,0,0,0,0,0,0,0});
}
Matrix<int> PointSymmetry::get_proper(const size_t i) const {
  Matrix<int> out{0,0,0, 0,0,0, 0,0,0};
  if (i<this->size()){
    out = this->R[i];
    if (this->isometry(i) < 0) for (size_t j=0; j<9u; ++j) out[j] *= -1;
  }
  return out;
}
Matrix<int> PointSymmetry::get_inverse(const size_t i) const {
  Matrix<int> identity{1,0,0, 0,1,0, 0,0,1}, result;
  if (i>=this->size()) return identity;
  // We *could* calculate R⁻¹ but might run into rounding/type issues.
  // We know that there must be an element of R for which RᵢRⱼ=E,
  // the identity element. Let's look for Rⱼ instead.
  for (size_t j=0; j<this->size(); ++j){
    brille::utils::multiply_matrix_matrix(result.data(), this->R[j].data(), this->R[i].data());
    if (brille::approx::matrix(3, identity.data(), result.data())) return this->R[j];
  }
  throw std::runtime_error("Incomplete pointgroup. Missing inverse of operation.");
}
size_t PointSymmetry::get_inverse_index(const size_t i) const {
  if (i>=this->size()) return i;
  Matrix<int> identity{1,0,0, 0,1,0, 0,0,1}, result;
  // We *could* calculate R⁻¹ but might run into rounding/type issues.
  // We know that there must be an element of R for which RᵢRⱼ=E,
  // the identity element. Let's look for Rⱼ instead.
  for (size_t j=0; j<this->size(); ++j){
    brille::utils::multiply_matrix_matrix(result.data(), this->R[j].data(), this->R[i].data());
    if (brille::approx::matrix(3, identity.data(), result.data())) return j;
  }
  throw std::runtime_error("Incomplete pointgroup. Missing inverse of operation.");
}
size_t PointSymmetry::find_index(const Matrix<int>& a) const {
  auto iseq = [ap=a.data()](const Matrix<int>& r){ return brille::approx::matrix(3,ap,r.data());};
  auto itr = std::find_if(this->R.begin(), this->R.end(), iseq);
  return std::distance(this->R.begin(), itr);
}
// const Matrix<int>& PointSymmetry::get(const size_t i) const {
//   if (i>=this->size())
//     throw std::out_of_range("The requested symmetry operation is out of range");
//   return this->R[i];
// }
size_t PointSymmetry::resize(const size_t newsize){
  this->R.resize(newsize);
  return this->R.size();
}
void PointSymmetry::print(const size_t i) const {
  for (size_t j=0; j<3u; ++j){
    for (size_t k=0; k<3u; ++k) std::cout << " " << R[i][j*3u+k];
    std::cout << std::endl;
  }
}
void PointSymmetry::print(void) const {
  for (size_t i=0; i<this->R.size(); ++i){
    this->print(i);
    std::cout << std::endl;
  }
}
void PointSymmetry::sort(const int mode) {
  auto oup = [](Matrix<int> a, Matrix<int> b){ return rotation_order(a.data()) < rotation_order(b.data());};
  auto odw = [](Matrix<int> a, Matrix<int> b){ return rotation_order(a.data()) > rotation_order(b.data());};
  auto iup = [](Matrix<int> a, Matrix<int> b){ return isometry_value(a.data()) < isometry_value(b.data());};
  auto idw = [](Matrix<int> a, Matrix<int> b){ return isometry_value(a.data()) > isometry_value(b.data());};
  switch(mode){
    // sort by isometry, decreasing
    case -2: std::sort(this->R.begin(),this->R.end(),idw); break;
    // sort by order, decreasing
    case -1: std::sort(this->R.begin(),this->R.end(),odw); break;
    // sort by isometry, increasing
    case  2: std::sort(this->R.begin(),this->R.end(),iup); break;
    // sort by order, increasing
    default: std::sort(this->R.begin(),this->R.end(),oup);
  }
}
void PointSymmetry::permute(const std::vector<size_t>& p){
  std::string msg;
  std::vector<size_t> o(this->size());
  std::iota(o.begin(), o.end(), 0u); // 0, 1, 2, …, size()-1
  std::vector<size_t> s=p;
  std::sort(s.begin(), s.end()); // includes requires *sorted* ranges
  if (!std::includes(o.begin(), o.end(), s.begin(), s.end()) || p.size()!=this->size()){
    msg="The provided permutation vector [";
    for (auto x: p) msg += " " + std::to_string(x);
    msg += " ] is invalid.";
    msg += " A permutation of [";
    for (auto x: o) msg += " " + std::to_string(x);
    msg += " ] was expected.";
    throw std::runtime_error(msg);
  }
  // swapping requires we have the inverse permutation
  for (size_t i=0; i<p.size(); ++i) s[p[i]] = i;
  // perform the actual element swapping
  for (size_t i=0; i<this->size();){
    if (s[i]==i) {
      ++i;
    } else {
      std::swap(this->R[i], this->R[s[i]]);
      std::swap(s[i], s[s[i]]);
    }
  }
  // confirm that the permuting worked:
  if (!std::is_sorted(s.begin(), s.end())){
    msg="Undoing the permutation [";
    for (auto x: p) msg += " " + std::to_string(x);
    msg +=" ] failed. End result is [";
    for (auto x: s) msg += " " + std::to_string(x);
    msg +=" ]";
    throw std::runtime_error(msg);
  }
}

std::vector<int> PointSymmetry::orders(void) const {
  std::vector<int> o;
  for (auto r: this->R) o.push_back(rotation_order(r.data()));
  return o;
}
std::vector<int> PointSymmetry::isometries(void) const {
  std::vector<int> i;
  for (auto r: this->R) i.push_back(isometry_value(r.data()));
  return i;
}
Vectors<int> PointSymmetry::axes(void) const {
  Vectors<int> ax;
  for (auto r: this->R) ax.push_back(brille::rotation_axis_and_perpendicular_vectors(r.data())[0]);
  return ax;
}

bool PointSymmetry::has_space_inversion() const {
  // space inversion means that a pointgroup operation has isometry -1.
  for (auto op: this->R) if (-1 == isometry_value(op.data())) return true;
  return false;
}

bool PointSymmetry::has(const Matrix<int>& m) const {
  for (auto r: R) if (brille::approx::matrix(m.data(), r.data())) return true;
  return false;
}
Matrix<int> PointSymmetry::pop(const size_t i) {
  Matrix<int> out = this->get(i);
  this->erase(i);
  return out;
}
size_t PointSymmetry::erase(const size_t i){
  if (i>=this->size())
    throw std::out_of_range("The requested symmetry operation is out of range");
  this->R.erase(this->R.begin()+i);
  return this->size();
}
int PointSymmetry::order(const size_t i) const {
  if (i>=this->size())
    throw std::out_of_range("The requested symmetry operation is out of range");
  return rotation_order(this->R[i].data());
}
int PointSymmetry::isometry(const size_t i) const {
  if (i>=this->size())
    throw std::out_of_range("The requested symmetry operation is out of range");
  return isometry_value(this->R[i].data());
}
Vector<int> PointSymmetry::axis(const size_t i) const {
  if (i>=this->size())
    throw std::out_of_range("The requested symmetry operation is out of range");
  return brille::rotation_axis_and_perpendicular_vectors(this->R[i].data())[0];
}
Vector<int> PointSymmetry::perpendicular_axis(const size_t i) const {
  if (i>=this->size())
    throw std::out_of_range("The requested symmetry operation is out of range");
  return brille::rotation_axis_and_perpendicular_vectors(this->R[i].data())[1];
}
Vectors<int> PointSymmetry::perpendicular_axes(void) const {
  Vectors<int> ax;
  for (auto r: this->R) ax.push_back(brille::rotation_axis_and_perpendicular_vectors(r.data())[1]);
  return ax;
}
PointSymmetry PointSymmetry::generate(void) const {
  Matrix<int> x, y, E{1,0,0, 0,1,0, 0,0,1}; // E, the identity operation
  PointSymmetry gen;
  gen.add(E); // any generated operations always need the identity operation
  size_t gensize;
  for (size_t i=0; i<this->size(); ++i){
    x = E;
    // for each of R, R¹,…, Rⁿ⁻¹
    for (int count=this->order(i); count--;){ // post-decrimate to run the loop count times
      brille::utils::multiply_matrix_matrix(y.data(), this->data(i), x.data());
      x = y;
      // multiply Rⁱ by the operators in gen
      gensize = gen.size(); // only multiply against pre-existing operations
      for (size_t j=0; j<gensize; ++j){
        brille::utils::multiply_matrix_matrix(y.data(), x.data(), gen.data(j));
        // and add them to gen if they are not already present
        if (!gen.has(y)) gen.add(y);
      }
    }
  }
  return gen;
}
PointSymmetry PointSymmetry::generators(void) const {
  PointSymmetry listA(this->getall()), listB, listC;
  while (listA.size()){
    // move the first element of A to B
    listB.add(listA.get(0));
    // use the operations in B as generators
    listC = listB.generate();
    // remove from A all generated operations (in C)
    for (size_t i=0; i<listA.size(); listC.has(listA.get(i))?listA.erase(i):++i);
  }
  // listA is empty, listB contains at least the generators
  // for each element of B, use all other elements to generate operations C
  // if the skipped element is in C, remove it from B as it is not a generator
  for (size_t i=0; i<listB.size(); listC.has(listB.get(i))?listB.erase(i):++i){
    listA.resize(0); // empty-out A again
    // copy every element but i from B to A
    for (size_t j=0; j<listB.size(); ++j) if (i!=j) listA.add(listB.get(j));
    // use the reduced list as generators
    listC = listA.generate();
  }
  // listB now contains a complete set of generators (which might not be unique for the pointgroup)
  return listB;
}
PointSymmetry PointSymmetry::nfolds(const int min_order) const {
  Vectors<int> all_axis = this->axes();
  std::vector<size_t> eqv_idx(this->size());
  std::vector<size_t> unq_idx;
  Vectors<int> unq_axis;
  for (size_t i=0; i<this->size(); ++i){
    bool isunique = true;
    for (auto j: unq_idx) if (brille::approx::vector(all_axis[j].data(), all_axis[i].data())){
      eqv_idx[i] = j;
      isunique = false;
      break;
    }
    if (isunique){
      unq_idx.push_back(i);
      eqv_idx[i]=i;
      unq_axis.push_back(all_axis[i]);
    }
  }
  // equal-valued entries of eqv_idx have the same stationary axis
  // but the order of the selected unique rotation may not be what we want
  for (size_t i=0; i<unq_idx.size(); ++i){
    size_t idx = unq_idx[i];
    for (size_t j=0; j<this->size(); ++j)
    if (idx == eqv_idx[j] && this->order(j) > this->order(unq_idx[i]))
    unq_idx[i] = j;
  }
  /* Now each unique index points to a symmetry operation of highest order along
     its unique axis. What is not determined at this point is whether each
     symmetry operation is N or ̄N -- both of which are of order N.
     In order to ensure we return just the rotation part of any symmetry
     operation, we must multiply any rotoinversion operations by the inversion
     operation.*/
  int onebar[9]{-1,0,0, 0,-1,0, 0,0,-1};
  PointSymmetry out;
  for (auto i: unq_idx)
  if (this->order(i) > min_order){
    out.add(this->get(i));
    if (this->isometry(i) < 0)
      brille::utils::multiply_matrix_matrix(out.data(out.size()-1), onebar, this->data(i));
  }
  out.sort(); // double-check that the orders are sorted
  return out;
}
PointSymmetry PointSymmetry::higher(const int min_order) const {
  PointSymmetry out;
  for (size_t i=0; i<this->size(); ++i)
  if (this->order(i) > min_order)
  out.add(this->get(i));
  out.sort();
  return out;
}
