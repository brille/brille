/* This file is part of brille.

Copyright © 2020 Greg Tucker <greg.tucker@stfc.ac.uk>

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
#ifndef BRILLE_HALL_SYMBOL_HPP_
#define BRILLE_HALL_SYMBOL_HPP_
#include <cstring>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include "bravais.hpp"
#include "symmetry.hpp" // for class Symmetry
#include "pointsymmetry.hpp" // and PointSymmetry
namespace brille {

// template<class T> using Matrix = std::array<T,9>;
// template<class T> using Vector = std::array<T,3>;
// template<class T> using Matrices = std::vector<Matrix<T>>;
// template<class T> using Vectors = std::vector<Vector<T>>;

class SeitzSymbol {
  int N;
  std::string T;
  std::string A;
public:
  explicit SeitzSymbol(const int n=0, const std::string& t="", const std::string& a=""): N{n}, T{t}, A{a} {};
  SeitzSymbol(const int n, const std::string& t, const char a): N{n}, T{t} {A=a;};
  int set_order(const int n){ N = n; return N; }
  int get_order() const {return N;}
  std::string set_tran(const std::string& t){ T = t; return T; }
  std::string set_axis(const std::string& a){ A = a; return A; }
  std::string add_tran(const char& t){ T += t; return T; }
  std::string add_axis(const char& a){ A += a; return A; }
  std::string get_tran() const { return T; }
  std::string get_axis() const { return A; }
  std::string to_ascii() const;
  bool validate();
  Matrix<int> getr() const;
  Matrix<int> getr(const SeitzSymbol& pre) const;
  Vector<double> gett() const;
  Vector<double> gett(const SeitzSymbol& pre) const;
  char implicit_axis() const { return A.size() ? A[0] : 'z'; }
  char implicit_axis(const SeitzSymbol& pre) const {
    char a = '!';
    if (A.size()){
      a = A[0];
    } else {
      int tn = std::abs(this->N);
      int pn = std::abs(pre.get_order());
      if (2 == tn){
        if (2 == pn || 4 == pn) a = 'x';
        if (3 == pn || 6 == pn) a = '?';
      }
      if (3 == tn) a = '*';
    }
    return a;
  }
protected:
  Matrix<int> inner_getr(const char, const char) const;
  Vector<double> inner_gett(int, char, char) const;
};

class HallSymbol {
  Bravais L;
  bool centrosymmetric;
  std::vector<SeitzSymbol> symbols;
  Motion<double,double> V; // change of basis '4×4' matrix
public:
  explicit HallSymbol(const Bravais b, const bool c, const std::vector<SeitzSymbol>& ss):
    L{b}, centrosymmetric{c}, symbols{ss} {};
  explicit HallSymbol(const Bravais b=Bravais::_, const bool c=false): L{b}, centrosymmetric{c} {};
  HallSymbol(const std::string& strrep){ this->from_ascii(strrep); };
  Bravais getl() const { return L; }
  Bravais setl(const Bravais& b){ L = b; return L; }
  Bravais setl(const char b, const bool c){
    switch (b) {
      case 'P': L = Bravais::P; break;
      case 'A': L = Bravais::A; break;
      case 'B': L = Bravais::B; break;
      case 'C': L = Bravais::C; break;
      case 'I': L = Bravais::I; break;
      case 'R': L = Bravais::R; break;
      case 'F': L = Bravais::F; break;
      default: L = Bravais::_;
    }
    centrosymmetric = c;
    return L;
  }
  bool setcs(const bool& cs){ centrosymmetric = cs; return centrosymmetric; }
  std::vector<SeitzSymbol> setsymbols(const std::vector<SeitzSymbol>& s){ symbols = s; return symbols; }
  std::vector<SeitzSymbol> addsymbol(SeitzSymbol s){ symbols.push_back(s); return symbols; }
  bool from_ascii(const std::string&);
  std::string to_ascii() const;
  bool validate();
  Symmetry get_generators() const;
};

} // namespace brille
#endif
