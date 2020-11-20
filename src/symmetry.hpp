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

/*! \file */

#ifndef BRILLE_SYMMETRY_H_
#define BRILLE_SYMMETRY_H_
#include <algorithm>
#include "utilities.hpp"
#include "bravais.hpp"
#include "approx.hpp"
#include "symmetry_common.hpp"
namespace brille {

// template<class T> using Matrix = std::array<T,9>;
// template<class T> using Vector = std::array<T,3>;
// template<class T> using Matrices = std::vector<Matrix<T>>;
// template<class T> using Vectors = std::vector<Vector<T>>;

template<class R, class T> class Motion{
  Matrix<R> W;
  Vector<T> w;
public:
  Motion(): W({1,0,0, 0,1,0, 0,0,1}), w({0,0,0}) {}
  Motion(const Matrix<R>& X): W(X), w({0,0,0}) {}
  Motion(const Vector<T>& x): W({1,0,0, 0,1,0, 0,0,1}), w(x) {}
  Motion(const Matrix<R>& X, const Vector<T>& x): W(X), w(x) {}
  Motion(const std::string& x, bool change_of_basis=false) {this->from_ascii(x,change_of_basis);}
  Motion<R,T> operator*(const Motion<R,T>& m) const {
    // {W0, w0} * {W1, w1} == {W0*W1, (W0*w1 + w0)%1}
    Matrix<R> outW;
    Vector<T> outw;
    // output rotation-part is W0*W1
    brille::utils::multiply_matrix_matrix(outW.data(), this->W.data(), m.getr().data());
    // ouput translation-part is W0*w1
    brille::utils::multiply_matrix_vector(outw.data(), this->W.data(), m.gett().data());
    // *plus* w0
    for (int i=0; i<3; ++i) outw[i] += this->w[i];
    // we want the output translation elements to be ∈ [0,1)
    for (int i=0; i<3; ++i) outw[i] -= std::floor(outw[i]);
    return Motion(outW, outw);
  }
  Motion<R,T> inverse() const {
    Matrix<R> outW;
    Vector<T> outw;
    // output rotation part is W⁻¹
    brille::utils::matrix_inverse(outW.data(), this->W.data());
    // output translation part is -W⁻¹*w
    brille::utils::multiply_matrix_vector(outw.data(), outW.data(), this->w.data());
    for (int i=0; i<3; ++i) outw[i] *= -1;
    // we want the output translation elements to be ∈ [0,1)
    for (int i=0; i<3; ++i) outw[i] -= std::floor(outw[i]);
    return Motion(outW, outw);
  }
  template<class S>
  Vector<S> move_point(const Vector<S>& point) const {
    Vector<S> outp;
    // (W,w) * ⃗p  = W * ⃗p + w
    brille::utils::multiply_matrix_vector(outp.data(), this->W.data(), point.data());
    for (int i=0; i<3; ++i) outp[i] += this->w[i];
    return outp;
  }
  template<class S>
  Vector<S> move_vector(const Vector<S>& v) const {
    Vector<S> outv;
    brille::utils::multiply_matrix_vector(outv.data(), this->W.data(), v.data());
    return outv;
  }
  template<class S>
  Vector<S> move_axial(const Vector<S>& a) const {
    Vector<S> outa = this->move_vector(a);
    R det = brille::utils::matrix_determinant(this->W.data());
    if (1 != det) for (int i=0; i<3; ++i) outa[i] *= det;
    return outa;
  }
  // template<class T, typename S/*=promotion stuff*/>
  // Vector<S> operator*(const Vector<T>& x) const;
  Matrix<R> getr(void) const {return W;}
  Vector<T> gett(void) const {return w;}
  bool operator==(const Motion<R,T>& m) const {
    Vector<T> z{{0,0,0}}, d{m.gett()};
    Matrix<R> mW{m.getr()};
    // construct the difference vector
    for (int i=0; i<3; ++i) d[i] = d[i]-w[i];
    // limit its elements to the range [0,1) and transform such that 0 and 1 are near zero
    for (int i=0; i<3; ++i) d[i] = T(0.5) - std::abs(d[i] - std::floor(d[i]) - T(0.5));
    // so if the difference vector is ~ ⃗0 or ~ ⃗1 and the matrix parts match
    // then the Motions are equivalent
    return brille::approx::vector(d.data(), z.data()) && brille::approx::matrix(W.data(), mW.data());
  }
  template<class S>
  bool equal_matrix(const Matrix<S>& m) const {
    return brille::approx::matrix(W.data(), m.data());
  }
  // Some compiler documentation (GCC,+?) claims that operator!= will be
  // automatically constructed if operator== is defined. This is apparently not
  // true for MSVC. This statement might need to be precompiler filtered if GCC
  // complains about it.
  bool operator!=(const Motion<R,T>& m) const { return !this->operator==(m); }
  bool has_identity_rotation() const {
    Matrix<R> i{{1,0,0, 0,1,0, 0,0,1}};
    return brille::approx::matrix(W.data(), i.data());
  }
  bool has_identity_translation() const {
    Vector<T> i{{0,0,0}};
    return brille::approx::vector(w.data(), i.data());
  }
  size_t from_ascii(const std::string& x, const bool cob=false);
};

template<class R, class T>
size_t Motion<R,T>::from_ascii(const std::string& s, const bool cob){
  // if R==T and the string has no commas, it represents a special version of a
  // change of basis matrix where W≡𝟙 and 12×w is specified as space-separated integers
  bool special = cob && std::is_same<R,T>::value && s.find(',')==std::string::npos;
  this->W ={0,0,0, 0,0,0, 0,0,0};
  if (special) this->W[0] = this->W[4] = this->W[8] = R(1); // set to 𝟙
  this->w = {0,0,0};
  // a valid encoded Motion will be of the form
  // ∑ₙiₙαₙ + f where each iₙ is an integer, αₙ∈{x,y,z}
  // and f one of 0, ±1/12, ±1/6, ±1/4, ±1/3, ±1/2
  // but we should support any rational or floating point value for f
  int n{0};
  R i{1};
  T f{0}; // always floating point
  char c{' '};
  size_t ret{0}; // number of characters read
  std::stringstream stream(s);
  std::string nosearch = special ? "-+ " : ";,xyz-+";
  std::string digits = "0123456789";
  while (stream.good()){
    c = char(stream.get());
    debug_update_if(special, "next character ", c);
    ++ret; // we've read another character
    switch (c){
      case ' ': // end of motion component in the special case, otherwise nothing
        if (special) {
          this->w[n++] = static_cast<T>(i)/T(12);
          debug_update("set w[",n-1,"] to ", this->w[n-1]);
        }
        break;
      case ';': // end of one motion -- must fall through for third component end
      case ',': // end of motion component
        i=1;
        this->w[n++] = f;
        f=0;
        break;
      case 'x':
        this->W[n*3] = i;
        i=1;
        break;
      case 'y':
        this->W[n*3+1] = i;
        i=1;
        break;
      case 'z':
        this->W[n*3+2] = i;
        i=1;
        break;
      default: // not ;,xyz
        debug_update_if(special, "which is not one of ';,xyz ' so search until one of '",nosearch,"' is found next");
        std::stringstream tmp;
        tmp << c; // if a number, tmp might be the sign.
        bool haspoint{false}, hasslash{false}, hasdigit{digits.find(c)!=std::string::npos};
        while (nosearch.find(char(stream.peek()))==std::string::npos && stream.good()){
          hasdigit |= digits.find(char(stream.peek()))!=std::string::npos;
          hasslash |= '/'==stream.peek();
          haspoint |= '.'==stream.peek();
          tmp << char(stream.get()); // without char() the output of stream.get() is interpreted as an integer (under MSVC, at least)
        }
        // Make sure the end-of-motion character is next even if we're at the end of the string
        if (!stream.good()) stream.putback(special ? ' ' : ';');
        debug_update_if(special, "Finished searching. Result ",tmp.str()," which");
        debug_update_if(special, " ", (haspoint ? "has" : "does not have"), " a decimal point");
        debug_update_if(special, " ", (hasslash ? "has" : "does not have"), " a slash");
        debug_update_if(special, " ", (hasdigit ? "has" : "does not have"), " a digit");
        if (haspoint && hasslash) throw std::logic_error("Number has both point and slash?!");
        if (!(haspoint^hasslash)){ // TT (error) or FF (integer)
          if (!hasdigit) tmp << '1'; // in case tmp == "+" or "-"
          tmp >> i;
        }
        if (haspoint) tmp >> f; // floating point
        if (hasslash){ // rational
          std::vector<std::string> numden;
          for (std::string a; std::getline(tmp, a, '/');) numden.push_back(a);
          if (numden.size() != 2u) throw std::runtime_error("String representation of fractional number does not have two parts.");
          int pre{std::stoi(numden[0])}, post{std::stoi(numden[1])};
          f = static_cast<T>(pre)/static_cast<T>(post);
        }
    }
  }
  return ret;
}
// template<class T>
// size_t Motion<T,T>::from_ascii(const std::string& s){
//   // this overload is for reading a general transformation 4x4 matrix, i.e., Motion<T,T>
//   this->W = {0,0,0, 0,0,0, 0,0,0}
//   this->w = {0,0,0};
//   int n{0}, i{1};
//   T f{0}; // always floating point
//   char c{' '};
//   size_t ret{0}; // number of characters read
//   std::istringstream stream(s);
//   } else {
//     this->W = {{1,0,0, 0,1,0, 0,0,1}};
//     while (stream.good()) {
//       c = stream.get();
//       ++ret;
//       switch (c){
//         case ' ': // end of a translation component
//         this->w[n++] = f;
//         f = 0;
//         break;
//         default:
//         std::stringstream tmp;
//         tmp << c;
//         bool haspoint{false}, hasslash{false};
//         while (std::strchr(";,xyz-+",stream.peek())==nullptr){
//           hasslash |= '/'==stream.peek();
//           haspoint |= '.'==stream.peek();
//           tmp << stream.get();
//         }
//         if (haspoint && hasslash) throw std::logic_error("Number has both point and slash?!");
//         if (!(haspoint^hasslash)){ // TT (error) or FF (integer)
//           tmp >> i;
//           f = static_cast<T>(pre)/T(12); // the normal case
//         }
//         if (haspoint) tmp >> f; // floating point
//         if (hasslash){ // rational
//           std::vector<std::string> numden;
//           for (std::string a; std::getline(tmp, a, '/');) numden.push_back(a);
//           if (numden.size() != 2u) throw std::runtime_error("String representation of fractional number does not have two parts.");
//           int pre{std::stoi(numden[0])}, post{std::stoi(numden[1])};
//           f = static_cast<T>(pre)/static_cast<T>(post);
//         }
//       }
//     }
//   }
//   return ret;
// }

/*
  The following is probably a very bad idea, but will work so long as it is only
  ever used when the rotation part of b is directly convertable to R, as is
  the case when it is 𝟙.
  For 'safety' this will enforced at runtime.
*/
template<class R, class T> Motion<R,T> operator*(const Motion<R,T>& a, const Motion<T,T>& b){
  if (!b.has_identity_rotation())
    throw std::runtime_error("Differing type Motion multiplication requires the identity rotation");
  Motion<R,T> newb(b.gett());
  return a*newb;
}
template<class R, class T> Motion<R,T> operator*(const Motion<T,T>& a, const Motion<R,T>& b){
  if (!a.has_identity_rotation())
    throw std::runtime_error("Differing type Motion multiplication requires the identity rotation");
  Motion<R,T> newa(a.gett());
  return newa*b;
}

/*****************************************************************************\
| Symmetry class                                                              |
|-----------------------------------------------------------------------------|
|   Holds N 3x3 matrices R and N 3 vectors T which are the motions of a       |
|   space group. A motion is the combination of a rotation and a translation. |
|-----------------------------------------------------------------------------|
| Member functions:                                                           |
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
class Symmetry{
public:
  using Motions = std::vector<Motion<int,double>>;
private:
  Motions M;
public:
  Symmetry(size_t n=0) { M.resize(n); }
  Symmetry(const Motions& m): M(m) {};
  Symmetry(const std::string& strrep){ this->from_ascii(strrep); };
  Bravais                getcentring() const;
  const Motions&         getallm() const {return this->M;}
  Matrices<int>          getallr() const;
  Vectors<double>        getallt() const;
  size_t                 size()    const { return M.size(); }
  size_t                 resize(const size_t i)                                ;
  size_t                 add(const int *r, const double *t)                    ;
  size_t                 add(const Matrix<int>&, const Vector<double>&)        ;
  size_t                 add(const Motion<int,double>&)                        ;
  size_t                 add(const std::string&)                               ;
  bool                   from_ascii(const std::string& s)                      ;
  size_t                 erase(const size_t i)                                 ;
  Motion<int,double>     pop(const size_t i)                                   ;
  bool                   has(const Motion<int,double>&)                   const;
  Matrix<int>            getr(const size_t i)                             const;
  Vector<double>         gett(const size_t i)                             const;
  Motion<int,double>     getm(const size_t i)                             const;
  int order(const size_t i) const;
  std::vector<int> orders(void) const;
  Symmetry generate(void) const;
  Symmetry generators(void) const;
  bool operator==(const Symmetry& other) const;
  bool operator!=(const Symmetry& other) const { return !this->operator==(other);}
  bool has_space_inversion() const;
  size_t  find_matrix_index(const Matrix<int>&) const;
};

} // end namespace brille
#endif
