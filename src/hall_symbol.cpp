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

#include "hall_symbol.hpp"
#include "bravais.hpp" // for bravais_letter
using namespace brille;

bool HallSymbol::from_ascii(const std::string& s){
  char c;
  bool hassubsup;
  this->L = Bravais::_;
  hassubsup = s.find('^') != std::string::npos;
  hassubsup |= s.find('_') != std::string::npos;
  // bool hasspace = s.find(' ') != std::string::npos;

  bool isneg = false, issup = false, issub = false;
  SeitzSymbol ssym;
  std::istringstream stream(s);
  while (stream.good()){
    c = char(stream.get());
    // if c is one of the Lattice symbols:
    if (std::strchr("PABCIRF",c)){
      this->setl(c, isneg);
      if (' '==stream.peek()) c = char(stream.get());
      isneg = false;
    }
    // the subscripts a, b, c, n, u, v, w, d, 1, 2, 3, 4, 5 are all translation
    // specifications. All but 2, 3, 4 are unambiguous so add them to the
    // Setiz symbol and reset the subscript flag.
    // if a subscript specification has been given or the temporary Seitz symbol
    // has non-zero order then a translation symbol should be added for sure
    //
    // non-subscripted 1, 2, 3, 4, 6 are all rotation order symbols
    // they *should* specify the start of a Seitz symbol, so if the temporary
    // Seitz symbol has non-zero order we should add it to the list and start a
    // new temporary Setiz symbol
    //
    // any of x, y, z, ', ", or * represents an axis specification
    // add the axis to the Seitz symbol and reset the superscript flag
    if (hassubsup){
      // superscript follows -- advance to the next character
      if ('^'==c) { issup = true; c = char(stream.get()); }
      // subscript follows -- advance to the next character
      if ('_'==c) { issub = true; c = char(stream.get()); }
      if (issub && std::strchr("abcnuvwd12345",c)) {
        ssym.add_tran(c);
        issub = false;
      } else if (std::strchr("12346",c)) {
        if (ssym.get_order()) this->addsymbol(ssym);
        ssym = SeitzSymbol((isneg?-1:1)*(c-'0'));
        isneg = false;
        issup = false;
      }
      if (issup && std::strchr("xyz'\"*",c)){
        ssym.add_axis(c);
        issup = false;
      }
    } else /* without ^ or _ we can only distinguish symbols by spaces */{
      if (ssym.get_order()){
        if (std::strchr("abcnuvwd12345",c)) ssym.add_tran(c);
        if (std::strchr("xyz'\"*",c)) ssym.add_axis(c);
      } else if (std::strchr("12346",c)) {
        ssym = SeitzSymbol((isneg?-1:1)*(c-'0'));
        isneg = false;
      }
    }
    // capture what should have been an overbar but in ASCII is a preceeding -
    if ('-'==c) isneg = true;
    // a space between parts of the Hall symbol signifies a completed (set of)
    // Seitz matrix(s)
    if (' '==c) {
      // if the temporary Seitz symbol has non-zero order store it away
      if (ssym.get_order()) this->addsymbol(ssym);
      // reset everything for a possible next symbol
      ssym = SeitzSymbol();
      issup = false;
      issub = false;
      isneg = false;
    }
    if ('('==c){
      std::string tmp;
      std::getline(stream, tmp, ')');
      this->V.from_ascii(tmp, true);
    }
  }
  // At the end of the string we might still have a valid temporary Seitz symbol
  if (ssym.get_order()) this->addsymbol(ssym);
  return this->validate();
}

bool HallSymbol::validate(){
  // verify that everything went ok
  bool ok = brille::bravais_is_known(this->L);
  for (auto ss: this->symbols) ok &= ss.validate();
  // further validation could be possible by finding the full spacegroup
  return ok;
}

Symmetry HallSymbol::get_generators() const {
  Symmetry gen;
  Matrix<int> r{{1,0,0, 0,1,0, 0,0,1}};
  Vectors<double> ts;
  double h{0.5}, t{1./3.}, d{2./3.};
  switch (L){
    case Bravais::P: ts = {{0,0,0}}; break;
    case Bravais::A: ts = {{0,0,0},{0,h,h}}; break;
    case Bravais::B: ts = {{0,0,0},{h,0,h}}; break;
    case Bravais::C: ts = {{0,0,0},{h,h,0}}; break;
    case Bravais::I: ts = {{0,0,0},{h,h,h}}; break;
    case Bravais::R: ts = {{0,0,0},{t,d,d},{d,t,t}}; break;
    case Bravais::F: ts = {{0,0,0},{0,h,h},{h,0,h},{h,h,0}}; break;
    default: throw std::runtime_error("Unknown lattice type");
  }
  for (auto & onet: ts) gen.add(r, onet);
  if (centrosymmetric) {
    r = {{-1,0,0, 0,-1,0, 0,0,-1}};
    for (auto & onet: ts) gen.add(r, onet);
  }
  // copy the symbols and ensure all axis specifiers are explicit:
  std::vector<SeitzSymbol> expsym;
  if (symbols.size()>0){
    expsym.push_back(SeitzSymbol(symbols[0].get_order(), symbols[0].get_tran(), symbols[0].implicit_axis()));
    for (size_t i=1; i<symbols.size(); ++i)
      expsym.push_back(SeitzSymbol(symbols[i].get_order(), symbols[i].get_tran(), symbols[i].implicit_axis(expsym[i-1])));
  }
  // with all axis symbols now (semi) explicit we can run this more easily
  if (expsym.size()>0){
    gen.add(expsym[0].getr(), expsym[0].gett());
    for (size_t i=1; i<expsym.size(); ++i)
      gen.add(expsym[i].getr(expsym[i-1]), expsym[i].gett(expsym[i-1]));
  }
  // perform the change of basis, if necessary:
  if (!this->V.has_identity_rotation())
    throw std::logic_error("HallSymbol can only handle origin shift change of basis, but can be extened if necessary.");
  if (!this->V.has_identity_translation()){
    Vector<double> invT = this->V.gett();
    for (auto & x: invT) x *= -1;
    Motion<double,double> invV(invT); // the rotation part is set to 𝟙
    Symmetry::Motions motions;
    // each motion is replaced by V M V⁻¹
    for (auto m: gen.getallm()) motions.push_back(this->V*(m*invV));
    gen = Symmetry(motions);
  }
  return gen;
}


std::string HallSymbol::to_ascii() const {
  std::string s;
  if (this->centrosymmetric) s += "-";
  s += brille::bravais_letter(this->L);
  for (auto z: this->symbols) s += " " + z.to_ascii();
  return s;
}

std::string SeitzSymbol::to_ascii() const {
  std::string s = std::to_string(this->N) + this->A + this->T;
  return s;
}

bool SeitzSymbol::validate() {
  // If we include E and -E, N can be ±1, ±2, ±3, ±4, or ±6
  if (N == 0 || N < -6 || N > 6 || N == -5 || N == 5) return false;
  // check the axis specification(s)
  int primes{0}, xyz{0};
  for (char c: A) switch (c) {
    case '*': primes += 3; if (xyz) return false; break;
    case '"': primes += 2; if (xyz) return false; break;
    case '\'': ++primes;   if (xyz) return false; break;
    case 'x':
    case 'y':
    case 'z': if (primes || xyz++) return false; break; // only one x, y, or z specification allowed
    default:
      return false; // non x, y, z, ', ", * not allowed
  }
  // replace valid combinations of ' and "
  if (primes==2) A = "\""; // ensure '' ≡ "
  if (primes==3) A = "*"; // and ''' ≡ *,
  if (primes>3) return false; // higher than 3 is a syntax error
  // if (xyz) -- everything should be fine

  // check the translation specification(s)
  std::string subscripts = "abcnuvwd12345";
  // a full check isn't possible at this point due to (possibly) implied
  // rotation axes. instead just ensure that there aren't multiples of any
  // of the allowed subscripts
  if (T.size()>0) for (char c: subscripts)
    if (std::count_if(T.begin(), T.end(), [c](const char& d){return d==c;}) > 1)
      return false;
  // getting this far *should* mean that the SeitzSymbol is valid
  return true;
}


Matrix<int> SeitzSymbol::getr() const {
  return this->inner_getr(this->implicit_axis(), '!');
}
Matrix<int> SeitzSymbol::getr(const SeitzSymbol& pre) const{
  return this->inner_getr(this->implicit_axis(pre), pre.implicit_axis());
}
Matrix<int> SeitzSymbol::inner_getr(const char ax, const char pax) const {
  Matrix<int> out{{0,0,0, 0,0,0, 0,0,0}};
  switch (std::abs(N)){
    case 1: out = { 1, 0, 0,  0, 1, 0,  0, 0, 1}; break;
    case 2: switch (ax) {
      case 'x': out = { 1, 0, 0,  0,-1, 0,  0, 0,-1}; break;
      case 'y': out = {-1, 0, 0,  0, 1, 0,  0, 0,-1}; break;
      case 'z': out = {-1, 0, 0,  0,-1, 0,  0, 0, 1}; break;
      // the special 2 preceeded by 3 or 6 case -> a-b (same as 2' preceeded by z)
      case '?': out = { 0,-1, 0, -1, 0, 0,  0, 0,-1}; break;
      case '\'': switch (pax) {
        case 'x': out = {-1, 0, 0,  0, 0,-1,  0,-1, 0}; break;
        case 'y': out = { 0, 0,-1,  0,-1, 0, -1, 0, 0}; break;
        case 'z': out = { 0,-1, 0, -1, 0, 0,  0, 0,-1}; break;
        default: throw std::runtime_error("Impossible pre-axis for ' in Seitz symbol");
      } break;
      case '"': switch (pax) {
        case 'x': out = {-1, 0, 0,  0, 0, 1,  0, 1, 0}; break;
        case 'y': out = { 0, 0, 1,  0,-1, 0,  1, 0, 0}; break;
        case 'z': out = { 0, 1, 0,  1, 0, 0,  0, 0,-1}; break;
        default: throw std::runtime_error("Impossible pre-axis for \" in Seitz symbol");
      } break;
    } break;
    case 3: switch (ax) {
      case 'x': out = { 1, 0, 0,  0, 0,-1,  0, 1,-1}; break;
      case 'y': out = {-1, 0, 1,  0, 1, 0, -1, 0, 0}; break;
      case 'z': out = { 0,-1, 0,  1,-1, 0,  0, 0, 1}; break;
      case '*': out = { 0, 0, 1,  1, 0, 0,  0, 1, 0}; break;
    } break;
    case 4: switch (ax) {
      case 'x': out = { 1, 0, 0,  0, 0,-1,  0, 1, 0}; break;
      case 'y': out = { 0, 0, 1,  0, 1, 0, -1, 0, 0}; break;
      case 'z': out = { 0,-1, 0,  1, 0, 0,  0, 0, 1}; break;
    } break;
    case 6: switch (ax) {
      case 'x': out = { 1, 0, 0,  0, 1,-1,  0, 1, 0}; break;
      case 'y': out = { 0, 0, 1,  0, 1, 0, -1, 0, 1}; break;
      case 'z': out = { 1,-1, 0,  1, 0, 0,  0, 0, 1}; break;
    } break;
    default: throw std::runtime_error("Impossible rotation order in Seitz symbol");
  }
  // space invert the rotation if necessary
  if (N<0) for (int i=0; i<9; ++i) out[i] *= -1;
  return out;
}

Vector<double> SeitzSymbol::gett() const {
  return this->inner_gett(std::abs(this->get_order()), this->implicit_axis(), '!');
}
Vector<double> SeitzSymbol::gett(const SeitzSymbol& pre) const {
  return this->inner_gett(std::abs(this->get_order()), this->implicit_axis(pre), pre.implicit_axis());
}

Vector<double> SeitzSymbol::inner_gett(const int an, const char ax, const char pax) const {
  Vector<double> a{{0,0,0}}, t{{0,0,0}};
  if (std::any_of(T.begin(), T.end(), [](const char c){return std::strchr("12345",c);})){
    switch (ax) {
      case 'x': a = {1,0,0}; break;
      case 'y': a = {0,1,0}; break;
      case 'z': a = {0,0,1}; break;
      case '*': a = {1,1,1}; break;
      case '\'': switch (pax) {
        case 'x': a = {0,1,-1}; break;
        case 'y': a = {1,0,-1}; break;
        case 'z': a = {1,-1,0}; break;
      } break;
      case '"': switch (pax) {
        case 'x': a = {0,1,1}; break;
        case 'y': a = {1,0,1}; break;
        case 'z': a = {1,1,0}; break;
      }
    }
  }
  double h{0.5}, q{0.25};
  for (auto & c: T) switch (c) {
    case 'a': t[0] += h; break;
    case 'b': t[1] += h; break;
    case 'c': t[2] += h; break;
    case 'n': t[0] += h; t[1] += h; t[2] += h; break;
    case 'u': t[0] += q; break;
    case 'v': t[1] += q; break;
    case 'w': t[2] += q; break;
    case 'd': t[0] += q; t[1] += q; t[2] += q; break;
    case '1': for (int i=0; i < 3; ++i) t[i] +=   a[i]/static_cast<double>(an); break;
    case '2': for (int i=0; i < 3; ++i) t[i] += 2*a[i]/static_cast<double>(an); break;
    case '3': for (int i=0; i < 3; ++i) t[i] += 3*a[i]/static_cast<double>(an); break;
    case '4': for (int i=0; i < 3; ++i) t[i] += 4*a[i]/static_cast<double>(an); break;
    case '5': for (int i=0; i < 3; ++i) t[i] += 5*a[i]/static_cast<double>(an); break;
  }
  // translations are bounded by [0,1)
  for (int i=0; i < 3; ++i) t[i] -= std::floor(t[i]);
  return t;
}
