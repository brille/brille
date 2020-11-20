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
#include "symmetry.hpp"
#include <iostream>
#include "pointgroup.hpp"
#include <algorithm>
using namespace brille;
/*****************************************************************************\
| Symmetry class Member functions:                                            |
|-----------------------------------------------------------------------------|
|   set             copy a given matrix and vector into R and T for motion i. |
|   get             copy the rotation and translation of motion i into the    |
|                   provided matrix and/or vector.                            |
|   rdata, tdata    return a pointer to the first element of the rotation or  |
|                   translation of motion i.                                  |
|   resize          grow or shrink the number of motions that the object can  |
|                   hold -- this necessitates a memory copy if the old and    |
|                   new sizes are non-zero.                                   |
|   size            return the number of motions the object can store.        |
\*****************************************************************************/
Matrices<int>   Symmetry::getallr() const{
  Matrices<int> r;
  for (auto m: this->M) r.push_back(m.getr());
  return r;
}
Vectors<double> Symmetry::getallt() const{
  Vectors<double> t;
  for (auto m: this->M) t.push_back(m.gett());
  return t;
}
size_t Symmetry::add(const int *r, const double *t){
  Matrix<int> newR;
  Vector<double> newT;
  for (size_t i=0; i<9; ++i) newR[i]=r[i];
  for (size_t i=0; i<3; ++i) newT[i]=t[i];
  return this->add(newR,newT);
}
size_t Symmetry::add(const Matrix<int>& r, const Vector<double>& t){
  return this->add(Motion<int,double>(r,t));
}
size_t Symmetry::add(const Motion<int,double>& m){
  this->M.push_back(m);
  return this->size();
}
bool Symmetry::from_ascii(const std::string& s){
  // split the string on ';' character(s) which denote individual motions
  std::vector<std::string> motions;
  std::istringstream stream(s);
  for (std::string m; std::getline(stream, m, ';'); ) motions.push_back(m);
  this->M.resize(motions.size());
  for (size_t i=0; i<motions.size(); ++i) this->M[i].from_ascii(motions[i]);
  return true;
}
size_t Symmetry::add(const std::string& s){
  Motion<int,double> m;
  m.from_ascii(s);
  return this->add(m);
}
Matrix<int> Symmetry::getr(const size_t i) const {
  return (i<this->size()) ? this->M[i].getr() : Matrix<int>({0,0,0,0,0,0,0,0,0});
}
Vector<double> Symmetry::gett(const size_t i) const {
  return (i<this->size()) ? this->M[i].gett() : Vector<double>({0,0,0});
}
Motion<int,double> Symmetry::getm(const size_t i) const {
  return (i<this->size()) ? this->M[i] : Motion<int,double>();
}
size_t Symmetry::resize(const size_t newsize){
  this->M.resize(newsize);
  return newsize;
}

bool Symmetry::has(const Motion<int,double>& om) const {
  for (auto m: this->M) if (m == om) return true;
  return false;
}
Motion<int,double> Symmetry::pop(const size_t i) {
  Motion<int,double> out = this->getm(i);
  this->erase(i);
  return out;
}
size_t  Symmetry::erase(const size_t i){
  if (i>=this->size())
    throw std::out_of_range("The requested symmetry operation is out of range");
  this->M.erase(this->M.begin()+i);
  return this->size();
}

int Symmetry::order(const size_t i) const {
  if (i>=this->size())
    throw std::out_of_range("The requested symmetry operation is out of range");
  Motion<int,double> e, t(this->getm(i)); // x initalised to {E,0}
  int ord{0};
  while (++ord < 10 && t != e) { // if M[i]==e then order is 1 and this is done
    // info_update("M^",ord," = {");
    // info_update(t.getr());
    // info_update("|",t.gett(),"}");
    t = this->getm(i)*t;
  }
  info_update_if(t!=e,"Order not found for Motion {\n",this->getr(i),"|\n",this->gett(i),"\n}");
  return ord;
}
std::vector<int> Symmetry::orders() const {
  std::vector<int> o;
  for (size_t i=0; i<this->size(); ++i) o.push_back(this->order(i));
  return o;
}
Symmetry Symmetry::generate() const {
  Motion<int,double> x, y, e;
  Symmetry gen;
  gen.add(e);
  size_t gensize;
  for (size_t i=0; i<this->size(); ++i){
    x = e;
    // for each of {W,w}, {W²,Ww+w}, {W³,W²w+Ww+w}, …
    for (int count=this->order(i); count--;) /* loop runs count times */ {
      x = this->getm(i)*x;
      // multiply {Wⁱ,Wⁱ⁻¹w+…} by the motions in gen
      gensize = gen.size(); // only multiply against pre-existing motions
      for (size_t j=0; j<gensize; ++j){
        y = x * gen.getm(j);
        // add missing motions to gen
        if (!gen.has(y)) gen.add(y);
      }
    }
  }
  return gen;
}
Symmetry Symmetry::generators(void) const{
  Symmetry lA(this->getallm()), lB, lC;
  while (lA.size()){
    lB.add(lA.getm(0)); // move the first element of lA to lB
    lC = lB.generate(); // generate all possible Motions from lB
    for (size_t i=0; i<lA.size(); lC.has(lA.getm(i))?lA.erase(i):++i); // remove generated Motions still in lA
  }
  // lA is empty. lB contains the generators and possibly extraneous Motions.
  // Go through lB, skipping one at a time and removing any which still show up in lC
  for (size_t i=0; i<lB.size(); lC.has(lB.getm(i))?lB.erase(i):++i){ // lC.has is evaluated *after* lA.generate below!
    lA.resize(0);
    // copy each of lB to lA, skipping the chosen entry
    for (size_t j=0; j<lB.size(); ++j) if (i!=j) lA.add(lB.getm(j));
    // generate all motions for the shortened list
    lC = lA.generate();
  }
  // what remains in lB can be used to generate all Motions in this->getallm()
  // this list may or may not be unique.
  return lB;
}

bool Symmetry::operator==(const Symmetry& other) const {
  // both must have the same number of Motions
  if (this->size() != other.size()) return false;
  // all motions in one must be in the other (order does not matter)
  for (const auto & m: this->M) if (!other.has(m)) return false;
  // if both conditions hold then they are equivalent
  return true;
}

bool Symmetry::has_space_inversion() const {
  Motion<int,double> space_inversion({{-1,0,0, 0,-1,0, 0,0,-1}},{{0.,0.,0.}});
  for (auto op: this->M) if (op == space_inversion) return true;
  return false;
}

size_t Symmetry::find_matrix_index(const Matrix<int>& m) const {
  auto itr = std::find_if(this->M.begin(), this->M.end(), [m](const Motion<int,double>& Mot){return Mot.equal_matrix(m);});
  return std::distance(this->M.begin(), itr);
}

Bravais Symmetry::getcentring() const {
  double ot{1./3.}, tt{2./3.};
  std::array<int,9> ii{{1,0,0, 0,1,0, 0,0,1}};
  std::array<double,3> wa{{0,0.5,0.5}}, wb{{0.5,0,0.5}}, wc{{0.5,0.5,0}},
                 wi{{0.5,0.5,0.5}}, wr0{{ot,tt,tt}}, wr1{{tt,ot,ot}};
  Motion<int,double> Ma(ii,wa), Mb(ii,wb), Mc(ii,wc), Mi(ii,wi), Mr0(ii,wr0), Mr1(ii,wr1);
  bool hasa = this->has(Ma);
  bool hasb = this->has(Mb);
  bool hasc = this->has(Mc);
  if (hasa && hasb && hasc) return Bravais::F; // three face centred
  if (hasa) return Bravais::A; // A face centred
  if (hasb) return Bravais::B; // B face centred
  if (hasc) return Bravais::C; // C face centred
  if (this->has(Mi)) return Bravais::I; // body centred
  if (this->has(Mr0) && this->has(Mr1)) return Bravais::R;
  // With none of the pure translations it *must* be primitive (or have a problem)
  return Bravais::P;
}
