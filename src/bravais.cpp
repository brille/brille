#include "bravais.hpp"


std::string bravais_string(const Bravais b){
  std::string repr;
  switch (b){
    case Bravais::_: repr = "centering error";        break;
    case Bravais::P: repr = "primitive";              break;
    case Bravais::A: repr = "A-face centred";         break;
    case Bravais::B: repr = "B-face centred";         break;
    case Bravais::C: repr = "C-face centred";         break;
    case Bravais::I: repr = "body centred";           break;
    case Bravais::F: repr = "face centred";           break;
    case Bravais::R: repr = "rhombohedrally centred"; break;
  }
  return repr;
}

char bravais_letter(const Bravais b){
  switch (b){
    case Bravais::_: return '!';
    case Bravais::P: return 'P';
    case Bravais::A: return 'A';
    case Bravais::B: return 'B';
    case Bravais::C: return 'C';
    case Bravais::I: return 'I';
    case Bravais::F: return 'F';
    case Bravais::R: return 'R';
  }
  return '\0';
}

bool bravais_is_known(const Bravais b){
  switch (b){
    case Bravais::P:
    case Bravais::A:
    case Bravais::B:
    case Bravais::C:
    case Bravais::I:
    case Bravais::F:
    case Bravais::R:
      return true;
    default:
      return false;
  }
  return false;
}
