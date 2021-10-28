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
/*! \file
    \author Greg Tucker
    \brief Classes to handle decoding Hall symbols
*/
#include <cstring>
// #include <string>
// #include <sstream>
// #include <iostream>
// #include <algorithm>
// #include "bravais.hpp"
#include "symmetry.hpp" // for class Symmetry
#include "pointsymmetry.hpp" // and PointSymmetry
namespace brille {

/*! \brief A class to represent a Seitz symbol

A Seitz matrix is a four-by-four matrix representing a symmetry operation and
is comprised of a three-by-three submatrix representation of a rotation or
rotoinversion matrix, and a three-element translation vector.

\f[
\begin{matrix}
R_{00} & R_{01} & R_{02} & v_0 \\
R_{10} & R_{11} & R_{12} & v_1 \\
R_{20} & R_{21} & R_{22} & v_2 \\
0      & 0      & 0      & 1
\end{matrix}
\f]

The matrix can be written in a shorthand notation as {R|v} and Seitz devised a
system to replace the three-by-three matrix R by a symbol denoting the
*order* of the matrix and the 'orientation' axis of the matrix.
\note In the case of a true rotation the orientation is the stationary axis of
      the matrix. In the case of a rotoinversion the orientation is the
      stationary axis of the rotation-part of the operation and is, therefore,
      *merely* inverted by the matrix.

In Seitz's notation the translation vector is written as '0' for the zero-vector
or as the (fractional) components of the vector directly.

This class does not handle true {R|v} Seitz symbols but instead can encode and
decode Hall's notation for a Seitz matrix, which uses the same principle but
different characters for the 'orientation' axis, and additionally encodes the
translation vector into the symbol.
*/
class SeitzSymbol {
  int N;         /**< Operation order */
  std::string T; /**< Translation part of symbol */
  std::string A; /**< Axis part of symbol */
public:
  explicit SeitzSymbol(const int n=0, const std::string& t="", const std::string& a=""): N{n}, T{t}, A{a} {};
  SeitzSymbol(const int n, const std::string& t, const char a): N{n}, T{t} {A=a;};
  //! Specify the rotation order of the symbol
  int set_order(const int n){ N = n; return N; }
  //! Get the rotation order of the symbol
  int get_order() const {return N;}
  //! Specify the encoded translation of the symbol
  std::string set_tran(const std::string& t){ T = t; return T; }
  //! Specify the encoded orientation axis of the symbol
  std::string set_axis(const std::string& a){ A = a; return A; }
  //! Append to the encoded translation of the symbol
  std::string add_tran(const char& t){ T += t; return T; }
  //! Append to the encoded orientation axis of the symbol
  std::string add_axis(const char& a){ A += a; return A; }
  //! Get the encoded translation of the symbol
  std::string get_tran() const { return T; }
  //! Get the encoded orientation axis of the symbol
  std::string get_axis() const { return A; }
  //! Construct the full encoded symbol in ASCII format
  std::string to_ascii() const;
  //! Check that the stored symbol parts are self consistent
  bool validate();
  /*! Decode the stored rotation order and orientation axis

  \warning Use only when no symbol preceeds this one in the Hall symbol.

  \return a three-by-three matrix representation of R
  */
  Matrix<int> getr() const;
  /*! \brief Decode the stored rotation order and orientation axis

  \param pre The preceeding symbol required to correctly infer missing symbol part(s)
  \return a three-by-three matrix representation of R
  */
  Matrix<int> getr(const SeitzSymbol& pre) const;
  /*! \brief Decode the stored translation vector

  \warning Use only when no symbol preceeds this one in the Hall symbol.

  \return a three-element vector representation of v
  */
  Vector<double> gett() const;
  /*! \brief Decode the stored translation vector

  \param pre The preceeding symbol required to correctly infer missing symbol part(s)
  \return a three-element vector representation of v
  */
  Vector<double> gett(const SeitzSymbol& pre) const;
  /*! \brief Determine the explicit or implicit orientation axis

  \warning Use only when no symbol preceeds this one in the Hall symbol.

  \return the stored orientation axis or 'z' which is always the first implicit axis
  */
  char implicit_axis() const { return A.size() ? A[0] : 'z'; }
  /*! \brief Determine the explicit or implicit orientation axis

  If a symbol does not contain an orientation axis then it uses an implicit
  axis which depends on the position of this symbol in the whole Hall symbol,
  this symbol's order, and possibly the preceeding symbol's order.

  If this is the first symbol the implicit axis is always (0,0,1), represented
  by character 'z'.
  If this symbol has order 3 the implicit axis is always (1,1,1), represented
  by character '*'.
  If this symbol has order 2 and the preceeding symbol has order 2 or 4 the
  implicit axis is (1,0,0), represented by character 'x'.
  If this symbol has order 2 and the preceeding symbol has order 3 or 6 the
  implixit axis is (1,-1,0), which does not have a character.

  \param pre The preceeding symbol required to correctly infer missing symbol part(s)
  \return the stored orientation axis or the inferred implict axis
  */
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

/*! \brief A class to interact with lattices in Hall notation

Hall proposed a compact notation to represent the Bravais lattice type and the
generators for the spacegroup of a lattice. In Hall's notation the Bravais
lattice type is represented by a single character and the generators are
expressed in a specialized notation for Seitz matrices which is allowed to
contain implicit information dependent on the preceeding Seitz matrix symbol.

Hall devised his notation to make the most use of the implicit rules and
thereby identified the most-compact set of symbols to describe all 530 unique
combinations of Bravais lattice, spacegroup symmetry, and origin choice.
However, as the generators of any group are generally not unique it is often
possible to construct more than one valid Hall symbol for the same unique
lattice.

This class contains methods to decode arbitrary Hall symbols and produce the
Bravais lattice type and the generators of the spacegroup.
*/
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
  //! Return the Bravais symbol
  Bravais getl() const { return L; }
  //! Specify the Bravais symbol
  Bravais setl(const Bravais& b){ L = b; return L; }
  /*! \brief Determine the Bravais symbol from a character

  \param b the Bravais character, one of the enumeration `brille::Bravais` keys
  \param c true for centrosymmetric, false for noncentrosymmetric spacegroup
  \return determined Bravais symbol
  */
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
  //! Specify the centrosymmetricity of the lattice
  bool setcs(const bool& cs){ centrosymmetric = cs; return centrosymmetric; }
  //! Specify the full list of partly-decoded Seitz symbol generators
  std::vector<SeitzSymbol> setsymbols(const std::vector<SeitzSymbol>& s){ symbols = s; return symbols; }
  //! Add a partly-decoded Seitz symbol to the list of generators
  std::vector<SeitzSymbol> addsymbol(SeitzSymbol s){ symbols.push_back(s); return symbols; }
  //! Decode a Hall symbol from an ASCII encoded string
  bool from_ascii(const std::string&);
  //! Encode this Hall symbol and return an ASCII string
  std::string to_ascii() const;
  //! Validate the internal consistency of this Hall symbol
  bool validate();
  //! Return a Symmetry object with the decoded generators in this Hall symbol
  Symmetry get_generators() const;
};

} // namespace brille
#endif
