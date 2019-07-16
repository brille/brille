#include "symmetry.h"
#include <iostream>
#include "pointgroup.h"
#include <algorithm>

/*****************************************************************************\
| Symmetry class Member functions:                                            |
|-----------------------------------------------------------------------------|
|   set, setrot, settran   copy a given matrix and/or vector into R and/or T  |
|                          for the motion i.                                  |
|   get, getrot, gettran   copy the rotation and/or translation of motion i   |
|                          into the provided matrix and/or vector.            |
|   getrot, gettran        return a pointer to the first element of the       |
|                          rotation or translation of motion i.               |
|   resize                 grow or shrink the number of motions that the      |
|                          object can hold -- this necessitates a memory copy |
|                          if the old and new sizes are non-zero.             |
|   size                   return the number of motions the object can store. |
\*****************************************************************************/
size_t Symmetry::add(const int *r, const double *t){
  std::array<int,9> newR;
  std::array<double,3> newT;
  for (size_t i=0; i<9; ++i) newR[i]=r[i];
  for (size_t i=0; i<3; ++i) newT[i]=t[i];
  return this->add(newR,newT);
}
size_t Symmetry::add(const std::array<int,9>& r, const std::array<double,3>& t){
  this->R.push_back(r);
  this->T.push_back(t);
  return this->size();
}
bool Symmetry::set(const size_t i, const int *r, const double *t){
  if ( i>=this->size() ) return false;
  for(size_t j=0; j<9; j++) this->R[i][j] = r[j];
  for(size_t j=0; j<3; j++) this->T[i][j] = t[j];
  return true;
}
bool Symmetry::setrot(const size_t i, const int *r){
  if ( i>=this->size() ) return false;
  for(size_t j=0; j<9; j++) this->R[i][j] = r[j];
  return true;
}
bool Symmetry::settran(const size_t i, const double *t){
  if ( i>=this->size() ) return false;
  for(size_t j=0; j<3; j++) this->T[i][j] = t[j];
  return true;
}
bool Symmetry::get(const size_t i, int *r, double *t) const {
  if ( i>=this->size() ) return false;
  for(size_t j=0; j<9; j++) r[j]=this->R[i][j];
  for(size_t j=0; j<3; j++) t[j]=this->T[i][j];
  return true;
}
bool Symmetry::getrot(const size_t i, int *r) const {
  if ( i>=this->size() ) return false;
  for(size_t j=0; j<9; j++) r[j]=this->R[i][j];
  return true;
}
bool Symmetry::gettran(const size_t i, double *t) const {
  if ( i>=this->size() ) return false;
  for(size_t j=0; j<3; j++) t[j]=this->T[i][j];
  return true;
}
int    * Symmetry::getrot (const size_t i) {
  return (i<this->size()) ? this->R[i].data() : nullptr;
}
double * Symmetry::gettran(const size_t i) {
  return (i<this->size()) ? this->T[i].data() : nullptr;
}
const int    * Symmetry::getrot (const size_t i) const {
  return (i<this->size()) ? this->R[i].data() : nullptr;
}
const double * Symmetry::gettran(const size_t i) const {
  return (i<this->size()) ? this->T[i].data() : nullptr;
}
std::array<int,9> Symmetry::getrotarray(const size_t i) const {
  return (i<this->size()) ? this->R[i] : std::array<int,9>({0,0,0,0,0,0,0,0,0});
}
std::array<double,3> Symmetry::gettranarray(const size_t i) const {
  return (i<this->size()) ? this->T[i] : std::array<double,3>({0,0,0});
}
size_t Symmetry::resize(const size_t newsize){
  this->R.resize(newsize);
  this->T.resize(newsize);
  return newsize;
}

/*****************************************************************************\
| PointSymmetry class Member functions:                                       |
|-----------------------------------------------------------------------------|
|   set      copy a given matrix into R at index i.                           |
|   get      copy the rotation at index i into the provided matrix.           |
|   resize   grow or shrink the number of rotations that the object can hold  |
|            -- this causes a memory copy if the old and new sizes are finite.|
|   size     return the number of rotations the object can/does store.        |
\*****************************************************************************/
size_t PointSymmetry::add(const int *r){
  std::array<int,9> newR;
  for (size_t i=0; i<9; ++i) newR[i] = r[i];
  return this->add(newR);
}
size_t PointSymmetry::add(const std::array<int,9>& r){
  this->R.push_back(r);
  return this->R.size();
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
int * PointSymmetry::get(const size_t i) {
  return (i<this->size()) ? this->R[i].data() : nullptr;
}
const int * PointSymmetry::get(const size_t i) const {
  return (i<this->size()) ? this->R[i].data() : nullptr;
}
std::array<int,9> PointSymmetry::getarray(const size_t i) const {
  return (i<this->size()) ? this->R[i] : std::array<int,9>({0,0,0,0,0,0,0,0,0});
}
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
void PointSymmetry::sort(const int ad) {
  if (ad<0)
  std::sort(this->R.begin(), this->R.end(),
    [](std::array<int,9> a, std::array<int,9> b){
      return rotation_order(a.data()) > rotation_order(b.data());
  });
  else
  std::sort(this->R.begin(), this->R.end(),
    [](std::array<int,9> a, std::array<int,9> b){
      return rotation_order(a.data()) < rotation_order(b.data());
  });
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
std::vector<std::array<int,3>> PointSymmetry::axes(void) const {
  std::vector<std::array<int,3>> ax;
  for (auto r: this->R) ax.push_back(rotation_axis_and_perpendicular_vectors(r.data())[0]);
  return ax;
}
