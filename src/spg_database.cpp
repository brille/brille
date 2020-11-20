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

/* This file has evolved from spg_database.cpp distributed as part of spglib.
   Changes have been made to introduce C++ style classes as well as other
   modifications to suit the needs of brille.
   spglib was licensed under the following BSD 3-clause license:

 Copyright (C) 2010 Atsushi Togo
 All rights reserved.

 This file is part of spglib. https://github.com/atztogo/spglib

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions
 are met:

 * Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in
   the documentation and/or other materials provided with the
   distribution.

 * Neither the name of the phonopy project nor the names of its
   contributors may be used to endorse or promote products derived
   from this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 POSSIBILITY OF SUCH DAMAGE. */

#include <stdlib.h>
#include <cstring>
#include "spg_database.hpp"
#include "hall_symbol.hpp"
using namespace brille;

/* Modifications from spglib original:
' CENTERING_ERROR' --> 'Bravais::_'
'       PRIMITIVE' --> 'Bravais::P'
'          A_FACE' --> 'Bravais::A'
'          B_FACE' --> 'Bravais::B'  (not present)
'          C_FACE' --> 'Bravais::C'
'            BODY' --> 'Bravais::I'
'            FACE' --> 'Bravais::F'
'        R_CENTER' --> 'Bravais::R'
'            BASE' --> 'Bravais::X' (not present)
2019-05-28 G Tucker */
/* The ALL_SPACEGROUP Spacegroup array has been switched from fixed-length char
   arrays to variable length const std::strings (once initialized)            */

/* In Hall symbols (3rd column), the double-prime is an escaped double quote*/
static const Spacegroup ALL_SPACEGROUPS[] = {
  {  0, ""      , ""                , ""                               , ""                   , ""          , ""     , Bravais::_,  0,   0},
  {  1, "C1^1"  , "P 1"             , "P 1"                            , "P 1"                , "P1"        , ""     , Bravais::P,  1,   1},
  {  2, "Ci^1"  , "-P 1"            , "P -1"                           , "P -1"               , "P-1"       , ""     , Bravais::P,  2,   2},
  {  3, "C2^1"  , "P 2y"            , "P 2 = P 1 2 1"                  , "P 1 2 1"            , "P2"        , "b"    , Bravais::P,  3,   3},
  {  3, "C2^1"  , "P 2"             , "P 2 = P 1 1 2"                  , "P 1 1 2"            , "P2"        , "c"    , Bravais::P,  3,   4},
  {  3, "C2^1"  , "P 2x"            , "P 2 = P 2 1 1"                  , "P 2 1 1"            , "P2"        , "a"    , Bravais::P,  3,   5},
  {  4, "C2^2"  , "P 2yb"           , "P 2_1 = P 1 2_1 1"              , "P 1 2_1 1"          , "P2_1"      , "b"    , Bravais::P,  3,   6},
  {  4, "C2^2"  , "P 2c"            , "P 2_1 = P 1 1 2_1"              , "P 1 1 2_1"          , "P2_1"      , "c"    , Bravais::P,  3,   7},
  {  4, "C2^2"  , "P 2xa"           , "P 2_1 = P 2_1 1 1"              , "P 2_1 1 1"          , "P2_1"      , "a"    , Bravais::P,  3,   8},
  {  5, "C2^3"  , "C 2y"            , "C 2 = C 1 2 1"                  , "C 1 2 1"            , "C2"        , "b1"   , Bravais::C,  3,   9},
  {  5, "C2^3"  , "A 2y"            , "C 2 = A 1 2 1"                  , "A 1 2 1"            , "C2"        , "b2"   , Bravais::A,  3,  10},
  {  5, "C2^3"  , "I 2y"            , "C 2 = I 1 2 1"                  , "I 1 2 1"            , "C2"        , "b3"   , Bravais::I,  3,  11},
  {  5, "C2^3"  , "A 2"             , "C 2 = A 1 1 2"                  , "A 1 1 2"            , "C2"        , "c1"   , Bravais::A,  3,  12},
  {  5, "C2^3"  , "B 2"             , "C 2 = B 1 1 2 = B 2"            , "B 1 1 2"            , "C2"        , "c2"   , Bravais::P,  3,  13},
  {  5, "C2^3"  , "I 2"             , "C 2 = I 1 1 2"                  , "I 1 1 2"            , "C2"        , "c3"   , Bravais::I,  3,  14},
  {  5, "C2^3"  , "B 2x"            , "C 2 = B 2 1 1"                  , "B 2 1 1"            , "C2"        , "a1"   , Bravais::P,  3,  15},
  {  5, "C2^3"  , "C 2x"            , "C 2 = C 2 1 1"                  , "C 2 1 1"            , "C2"        , "a2"   , Bravais::C,  3,  16},
  {  5, "C2^3"  , "I 2x"            , "C 2 = I 2 1 1"                  , "I 2 1 1"            , "C2"        , "a3"   , Bravais::I,  3,  17},
  {  6, "Cs^1"  , "P -2y"           , "P m = P 1 m 1"                  , "P 1 m 1"            , "Pm"        , "b"    , Bravais::P,  4,  18},
  {  6, "Cs^1"  , "P -2"            , "P m = P 1 1 m"                  , "P 1 1 m"            , "Pm"        , "c"    , Bravais::P,  4,  19},
  {  6, "Cs^1"  , "P -2x"           , "P m = P m 1 1"                  , "P m 1 1"            , "Pm"        , "a"    , Bravais::P,  4,  20},
  {  7, "Cs^2"  , "P -2yc"          , "P c = P 1 c 1"                  , "P 1 c 1"            , "Pc"        , "b1"   , Bravais::P,  4,  21},
  {  7, "Cs^2"  , "P -2yac"         , "P c = P 1 n 1"                  , "P 1 n 1"            , "Pc"        , "b2"   , Bravais::P,  4,  22},
  {  7, "Cs^2"  , "P -2ya"          , "P c = P 1 a 1"                  , "P 1 a 1"            , "Pc"        , "b3"   , Bravais::P,  4,  23},
  {  7, "Cs^2"  , "P -2a"           , "P c = P 1 1 a"                  , "P 1 1 a"            , "Pc"        , "c1"   , Bravais::P,  4,  24},
  {  7, "Cs^2"  , "P -2ab"          , "P c = P 1 1 n"                  , "P 1 1 n"            , "Pc"        , "c2"   , Bravais::P,  4,  25},
  {  7, "Cs^2"  , "P -2b"           , "P c = P 1 1 b = P b"            , "P 1 1 b"            , "Pc"        , "c3"   , Bravais::P,  4,  26},
  {  7, "Cs^2"  , "P -2xb"          , "P c = P b 1 1"                  , "P b 1 1"            , "Pc"        , "a1"   , Bravais::P,  4,  27},
  {  7, "Cs^2"  , "P -2xbc"         , "P c = P n 1 1"                  , "P n 1 1"            , "Pc"        , "a2"   , Bravais::P,  4,  28},
  {  7, "Cs^2"  , "P -2xc"          , "P c = P c 1 1"                  , "P c 1 1"            , "Pc"        , "a3"   , Bravais::P,  4,  29},
  {  8, "Cs^3"  , "C -2y"           , "C m = C 1 m 1"                  , "C 1 m 1"            , "Cm"        , "b1"   , Bravais::C,  4,  30},
  {  8, "Cs^3"  , "A -2y"           , "C m = A 1 m 1"                  , "A 1 m 1"            , "Cm"        , "b2"   , Bravais::A,  4,  31},
  {  8, "Cs^3"  , "I -2y"           , "C m = I 1 m 1"                  , "I 1 m 1"            , "Cm"        , "b3"   , Bravais::I,  4,  32},
  {  8, "Cs^3"  , "A -2"            , "C m = A 1 1 m"                  , "A 1 1 m"            , "Cm"        , "c1"   , Bravais::A,  4,  33},
  {  8, "Cs^3"  , "B -2"            , "C m = B 1 1 m = B m"            , "B 1 1 m"            , "Cm"        , "c2"   , Bravais::P,  4,  34},
  {  8, "Cs^3"  , "I -2"            , "C m = I 1 1 m"                  , "I 1 1 m"            , "Cm"        , "c3"   , Bravais::I,  4,  35},
  {  8, "Cs^3"  , "B -2x"           , "C m = B m 1 1"                  , "B m 1 1"            , "Cm"        , "a1"   , Bravais::P,  4,  36},
  {  8, "Cs^3"  , "C -2x"           , "C m = C m 1 1"                  , "C m 1 1"            , "Cm"        , "a2"   , Bravais::C,  4,  37},
  {  8, "Cs^3"  , "I -2x"           , "C m = I m 1 1"                  , "I m 1 1"            , "Cm"        , "a3"   , Bravais::I,  4,  38},
  {  9, "Cs^4"  , "C -2yc"          , "C c = C 1 c 1"                  , "C 1 c 1"            , "Cc"        , "b1"   , Bravais::C,  4,  39},
  {  9, "Cs^4"  , "A -2yac"         , "C c = A 1 n 1"                  , "A 1 n 1"            , "Cc"        , "b2"   , Bravais::A,  4,  40},
  {  9, "Cs^4"  , "I -2ya"          , "C c = I 1 a 1"                  , "I 1 a 1"            , "Cc"        , "b3"   , Bravais::I,  4,  41},
  {  9, "Cs^4"  , "A -2ya"          , "C c = A 1 a 1"                  , "A 1 a 1"            , "Cc"        , "-b1"  , Bravais::A,  4,  42},
  {  9, "Cs^4"  , "C -2ybc"         , "C c = C 1 n 1"                  , "C 1 n 1"            , "Cc"        , "-b2"  , Bravais::C,  4,  43},
  {  9, "Cs^4"  , "I -2yc"          , "C c = I 1 c 1"                  , "I 1 c 1"            , "Cc"        , "-b3"  , Bravais::I,  4,  44},
  {  9, "Cs^4"  , "A -2a"           , "C c = A 1 1 a"                  , "A 1 1 a"            , "Cc"        , "c1"   , Bravais::A,  4,  45},
  {  9, "Cs^4"  , "B -2bc"          , "C c = B 1 1 n"                  , "B 1 1 n"            , "Cc"        , "c2"   , Bravais::P,  4,  46},
  {  9, "Cs^4"  , "I -2b"           , "C c = I 1 1 b"                  , "I 1 1 b"            , "Cc"        , "c3"   , Bravais::I,  4,  47},
  {  9, "Cs^4"  , "B -2b"           , "C c = B 1 1 b = B b"            , "B 1 1 b"            , "Cc"        , "-c1"  , Bravais::P,  4,  48},
  {  9, "Cs^4"  , "A -2ac"          , "C c = A 1 1 n"                  , "A 1 1 n"            , "Cc"        , "-c2"  , Bravais::A,  4,  49},
  {  9, "Cs^4"  , "I -2a"           , "C c = I 1 1 a"                  , "I 1 1 a"            , "Cc"        , "-c3"  , Bravais::I,  4,  50},
  {  9, "Cs^4"  , "B -2xb"          , "C c = B b 1 1"                  , "B b 1 1"            , "Cc"        , "a1"   , Bravais::P,  4,  51},
  {  9, "Cs^4"  , "C -2xbc"         , "C c = C n 1 1"                  , "C n 1 1"            , "Cc"        , "a2"   , Bravais::C,  4,  52},
  {  9, "Cs^4"  , "I -2xc"          , "C c = I c 1 1"                  , "I c 1 1"            , "Cc"        , "a3"   , Bravais::I,  4,  53},
  {  9, "Cs^4"  , "C -2xc"          , "C c = C c 1 1"                  , "C c 1 1"            , "Cc"        , "-a1"  , Bravais::C,  4,  54},
  {  9, "Cs^4"  , "B -2xbc"         , "C c = B n 1 1"                  , "B n 1 1"            , "Cc"        , "-a2"  , Bravais::P,  4,  55},
  {  9, "Cs^4"  , "I -2xb"          , "C c = I b 1 1"                  , "I b 1 1"            , "Cc"        , "-a3"  , Bravais::I,  4,  56},
  { 10, "C2h^1" , "-P 2y"           , "P 2/m = P 1 2/m 1"              , "P 1 2/m 1"          , "P2/m"      , "b"    , Bravais::P,  5,  57},
  { 10, "C2h^1" , "-P 2"            , "P 2/m = P 1 1 2/m"              , "P 1 1 2/m"          , "P2/m"      , "c"    , Bravais::P,  5,  58},
  { 10, "C2h^1" , "-P 2x"           , "P 2/m = P 2/m 1 1"              , "P 2/m 1 1"          , "P2/m"      , "a"    , Bravais::P,  5,  59},
  { 11, "C2h^2" , "-P 2yb"          , "P 2_1/m = P 1 2_1/m 1"          , "P 1 2_1/m 1"        , "P2_1/m"    , "b"    , Bravais::P,  5,  60},
  { 11, "C2h^2" , "-P 2c"           , "P 2_1/m = P 1 1 2_1/m"          , "P 1 1 2_1/m"        , "P2_1/m"    , "c"    , Bravais::P,  5,  61},
  { 11, "C2h^2" , "-P 2xa"          , "P 2_1/m = P 2_1/m 1 1"          , "P 2_1/m 1 1"        , "P2_1/m"    , "a"    , Bravais::P,  5,  62},
  { 12, "C2h^3" , "-C 2y"           , "C 2/m = C 1 2/m 1"              , "C 1 2/m 1"          , "C2/m"      , "b1"   , Bravais::C,  5,  63},
  { 12, "C2h^3" , "-A 2y"           , "C 2/m = A 1 2/m 1"              , "A 1 2/m 1"          , "C2/m"      , "b2"   , Bravais::A,  5,  64},
  { 12, "C2h^3" , "-I 2y"           , "C 2/m = I 1 2/m 1"              , "I 1 2/m 1"          , "C2/m"      , "b3"   , Bravais::I,  5,  65},
  { 12, "C2h^3" , "-A 2"            , "C 2/m = A 1 1 2/m"              , "A 1 1 2/m"          , "C2/m"      , "c1"   , Bravais::A,  5,  66},
  { 12, "C2h^3" , "-B 2"            , "C 2/m = B 1 1 2/m = B 2/m"      , "B 1 1 2/m"          , "C2/m"      , "c2"   , Bravais::P,  5,  67},
  { 12, "C2h^3" , "-I 2"            , "C 2/m = I 1 1 2/m"              , "I 1 1 2/m"          , "C2/m"      , "c3"   , Bravais::I,  5,  68},
  { 12, "C2h^3" , "-B 2x"           , "C 2/m = B 2/m 1 1"              , "B 2/m 1 1"          , "C2/m"      , "a1"   , Bravais::P,  5,  69},
  { 12, "C2h^3" , "-C 2x"           , "C 2/m = C 2/m 1 1"              , "C 2/m 1 1"          , "C2/m"      , "a2"   , Bravais::C,  5,  70},
  { 12, "C2h^3" , "-I 2x"           , "C 2/m = I 2/m 1 1"              , "I 2/m 1 1"          , "C2/m"      , "a3"   , Bravais::I,  5,  71},
  { 13, "C2h^4" , "-P 2yc"          , "P 2/c = P 1 2/c 1"              , "P 1 2/c 1"          , "P2/c"      , "b1"   , Bravais::P,  5,  72},
  { 13, "C2h^4" , "-P 2yac"         , "P 2/c = P 1 2/n 1"              , "P 1 2/n 1"          , "P2/c"      , "b2"   , Bravais::P,  5,  73},
  { 13, "C2h^4" , "-P 2ya"          , "P 2/c = P 1 2/a 1"              , "P 1 2/a 1"          , "P2/c"      , "b3"   , Bravais::P,  5,  74},
  { 13, "C2h^4" , "-P 2a"           , "P 2/c = P 1 1 2/a"              , "P 1 1 2/a"          , "P2/c"      , "c1"   , Bravais::P,  5,  75},
  { 13, "C2h^4" , "-P 2ab"          , "P 2/c = P 1 1 2/n"              , "P 1 1 2/n"          , "P2/c"      , "c2"   , Bravais::P,  5,  76},
  { 13, "C2h^4" , "-P 2b"           , "P 2/c = P 1 1 2/b = P 2/b"      , "P 1 1 2/b"          , "P2/c"      , "c3"   , Bravais::P,  5,  77},
  { 13, "C2h^4" , "-P 2xb"          , "P 2/c = P 2/b 1 1"              , "P 2/b 1 1"          , "P2/c"      , "a1"   , Bravais::P,  5,  78},
  { 13, "C2h^4" , "-P 2xbc"         , "P 2/c = P 2/n 1 1"              , "P 2/n 1 1"          , "P2/c"      , "a2"   , Bravais::P,  5,  79},
  { 13, "C2h^4" , "-P 2xc"          , "P 2/c = P 2/c 1 1"              , "P 2/c 1 1"          , "P2/c"      , "a3"   , Bravais::P,  5,  80},
  { 14, "C2h^5" , "-P 2ybc"         , "P 2_1/c = P 1 2_1/c 1"          , "P 1 2_1/c 1"        , "P2_1/c"    , "b1"   , Bravais::P,  5,  81},
  { 14, "C2h^5" , "-P 2yn"          , "P 2_1/c = P 1 2_1/n 1"          , "P 1 2_1/n 1"        , "P2_1/c"    , "b2"   , Bravais::P,  5,  82},
  { 14, "C2h^5" , "-P 2yab"         , "P 2_1/c = P 1 2_1/a 1"          , "P 1 2_1/a 1"        , "P2_1/c"    , "b3"   , Bravais::P,  5,  83},
  { 14, "C2h^5" , "-P 2ac"          , "P 2_1/c = P 1 1 2_1/a"          , "P 1 1 2_1/a"        , "P2_1/c"    , "c1"   , Bravais::P,  5,  84},
  { 14, "C2h^5" , "-P 2n"           , "P 2_1/c = P 1 1 2_1/n"          , "P 1 1 2_1/n"        , "P2_1/c"    , "c2"   , Bravais::P,  5,  85},
  { 14, "C2h^5" , "-P 2bc"          , "P 2_1/c = P 1 1 2_1/b = P 2_1/b", "P 1 1 2_1/b"        , "P2_1/c"    , "c3"   , Bravais::P,  5,  86},
  { 14, "C2h^5" , "-P 2xab"         , "P 2_1/c = P 2_1/b 1 1"          , "P 2_1/b 1 1"        , "P2_1/c"    , "a1"   , Bravais::P,  5,  87},
  { 14, "C2h^5" , "-P 2xn"          , "P 2_1/c = P 2_1/n 1 1"          , "P 2_1/n 1 1"        , "P2_1/c"    , "a2"   , Bravais::P,  5,  88},
  { 14, "C2h^5" , "-P 2xac"         , "P 2_1/c = P 2_1/c 1 1"          , "P 2_1/c 1 1"        , "P2_1/c"    , "a3"   , Bravais::P,  5,  89},
  { 15, "C2h^6" , "-C 2yc"          , "C 2/c = C 1 2/c 1"              , "C 1 2/c 1"          , "C2/c"      , "b1"   , Bravais::C,  5,  90},
  { 15, "C2h^6" , "-A 2yac"         , "C 2/c = A 1 2/n 1"              , "A 1 2/n 1"          , "C2/c"      , "b2"   , Bravais::A,  5,  91},
  { 15, "C2h^6" , "-I 2ya"          , "C 2/c = I 1 2/a 1"              , "I 1 2/a 1"          , "C2/c"      , "b3"   , Bravais::I,  5,  92},
  { 15, "C2h^6" , "-A 2ya"          , "C 2/c = A 1 2/a 1"              , "A 1 2/a 1"          , "C2/c"      , "-b1"  , Bravais::A,  5,  93},
  { 15, "C2h^6" , "-C 2ybc"         , "C 2/c = C 1 2/n 1"              , "C 1 2/n 1"          , "C2/c"      , "-b2"  , Bravais::C,  5,  94},
  { 15, "C2h^6" , "-I 2yc"          , "C 2/c = I 1 2/c 1"              , "I 1 2/c 1"          , "C2/c"      , "-b3"  , Bravais::I,  5,  95},
  { 15, "C2h^6" , "-A 2a"           , "C 2/c = A 1 1 2/a"              , "A 1 1 2/a"          , "C2/c"      , "c1"   , Bravais::A,  5,  96},
  { 15, "C2h^6" , "-B 2bc"          , "C 2/c = B 1 1 2/n"              , "B 1 1 2/n"          , "C2/c"      , "c2"   , Bravais::P,  5,  97},
  { 15, "C2h^6" , "-I 2b"           , "C 2/c = I 1 1 2/b"              , "I 1 1 2/b"          , "C2/c"      , "c3"   , Bravais::I,  5,  98},
  { 15, "C2h^6" , "-B 2b"           , "C 2/c = B 1 1 2/b = B 2/b"      , "B 1 1 2/b"          , "C2/c"      , "-c1"  , Bravais::P,  5,  99},
  { 15, "C2h^6" , "-A 2ac"          , "C 2/c = A 1 1 2/n"              , "A 1 1 2/n"          , "C2/c"      , "-c2"  , Bravais::A,  5, 100},
  { 15, "C2h^6" , "-I 2a"           , "C 2/c = I 1 1 2/a"              , "I 1 1 2/a"          , "C2/c"      , "-c3"  , Bravais::I,  5, 101},
  { 15, "C2h^6" , "-B 2xb"          , "C 2/c = B 2/b 1 1"              , "B 2/b 1 1"          , "C2/c"      , "a1"   , Bravais::P,  5, 102},
  { 15, "C2h^6" , "-C 2xbc"         , "C 2/c = C 2/n 1 1"              , "C 2/n 1 1"          , "C2/c"      , "a2"   , Bravais::C,  5, 103},
  { 15, "C2h^6" , "-I 2xc"          , "C 2/c = I 2/c 1 1"              , "I 2/c 1 1"          , "C2/c"      , "a3"   , Bravais::I,  5, 104},
  { 15, "C2h^6" , "-C 2xc"          , "C 2/c = C 2/c 1 1"              , "C 2/c 1 1"          , "C2/c"      , "-a1"  , Bravais::C,  5, 105},
  { 15, "C2h^6" , "-B 2xbc"         , "C 2/c = B 2/n 1 1"              , "B 2/n 1 1"          , "C2/c"      , "-a2"  , Bravais::P,  5, 106},
  { 15, "C2h^6" , "-I 2xb"          , "C 2/c = I 2/b 1 1"              , "I 2/b 1 1"          , "C2/c"      , "-a3"  , Bravais::I,  5, 107},
  { 16, "D2^1"  , "P 2 2"           , "P 2 2 2"                        , "P 2 2 2"            , "P222"      , ""     , Bravais::P,  6, 108},
  { 17, "D2^2"  , "P 2c 2"          , "P 2 2 2_1"                      , "P 2 2 2_1"          , "P222_1"    , ""     , Bravais::P,  6, 109},
  { 17, "D2^2"  , "P 2a 2a"         , "P 2_1 2 2"                      , "P 2_1 2 2"          , "P2_122"    , "cab"  , Bravais::P,  6, 110},
  { 17, "D2^2"  , "P 2 2b"          , "P 2 2_1 2"                      , "P 2 2_1 2"          , "P22_12"    , "bca"  , Bravais::P,  6, 111},
  { 18, "D2^3"  , "P 2 2ab"         , "P 2_1 2_1 2"                    , "P 2_1 2_1 2"        , "P2_12_12"  , ""     , Bravais::P,  6, 112},
  { 18, "D2^3"  , "P 2bc 2"         , "P 2 2_1 2_1"                    , "P 2 2_1 2_1"        , "P22_12_1"  , "cab"  , Bravais::P,  6, 113},
  { 18, "D2^3"  , "P 2ac 2ac"       , "P 2_1 2 2_1"                    , "P 2_1 2 2_1"        , "P2_122_1"  , "bca"  , Bravais::P,  6, 114},
  { 19, "D2^4"  , "P 2ac 2ab"       , "P 2_1 2_1 2_1"                  , "P 2_1 2_1 2_1"      , "P2_12_12_1", ""     , Bravais::P,  6, 115},
  { 20, "D2^5"  , "C 2c 2"          , "C 2 2 2_1"                      , "C 2 2 2_1"          , "C222_1"    , ""     , Bravais::C,  6, 116},
  { 20, "D2^5"  , "A 2a 2a"         , "A 2_1 2 2"                      , "A 2_1 2 2"          , "A2_122"    , "cab"  , Bravais::A,  6, 117},
  { 20, "D2^5"  , "B 2 2b"          , "B 2 2_1 2"                      , "B 2 2_1 2"          , "B22_12"    , "bca"  , Bravais::P,  6, 118},
  { 21, "D2^6"  , "C 2 2"           , "C 2 2 2"                        , "C 2 2 2"            , "C222"      , ""     , Bravais::C,  6, 119},
  { 21, "D2^6"  , "A 2 2"           , "A 2 2 2"                        , "A 2 2 2"            , "A222"      , "cab"  , Bravais::A,  6, 120},
  { 21, "D2^6"  , "B 2 2"           , "B 2 2 2"                        , "B 2 2 2"            , "B222"      , "bca"  , Bravais::P,  6, 121},
  { 22, "D2^7"  , "F 2 2"           , "F 2 2 2"                        , "F 2 2 2"            , "F222"      , ""     , Bravais::F,  6, 122},
  { 23, "D2^8"  , "I 2 2"           , "I 2 2 2"                        , "I 2 2 2"            , "I222"      , ""     , Bravais::I,  6, 123},
  { 24, "D2^9"  , "I 2b 2c"         , "I 2_1 2_1 2_1"                  , "I 2_1 2_1 2_1"      , "I2_12_12_1", ""     , Bravais::I,  6, 124},
  { 25, "C2v^1" , "P 2 -2"          , "P m m 2"                        , "P m m 2"            , "Pmm2"      , ""     , Bravais::P,  7, 125},
  { 25, "C2v^1" , "P -2 2"          , "P 2 m m"                        , "P 2 m m"            , "P2mm"      , "cab"  , Bravais::P,  7, 126},
  { 25, "C2v^1" , "P -2 -2"         , "P m 2 m"                        , "P m 2 m"            , "Pm2m"      , "bca"  , Bravais::P,  7, 127},
  { 26, "C2v^2" , "P 2c -2"         , "P m c 2_1"                      , "P m c 2_1"          , "Pmc2_1"    , ""     , Bravais::P,  7, 128},
  { 26, "C2v^2" , "P 2c -2c"        , "P c m 2_1"                      , "P c m 2_1"          , "Pcm2_1"    , "ba-c" , Bravais::P,  7, 129},
  { 26, "C2v^2" , "P -2a 2a"        , "P 2_1 m a"                      , "P 2_1 m a"          , "P2_1ma"    , "cab"  , Bravais::P,  7, 130},
  { 26, "C2v^2" , "P -2 2a"         , "P 2_1 a m"                      , "P 2_1 a m"          , "P2_1am"    , "-cba" , Bravais::P,  7, 131},
  { 26, "C2v^2" , "P -2 -2b"        , "P b 2_1 m"                      , "P b 2_1 m"          , "Pb2_1m"    , "bca"  , Bravais::P,  7, 132},
  { 26, "C2v^2" , "P -2b -2"        , "P m 2_1 b"                      , "P m 2_1 b"          , "Pm2_1b"    , "a-cb" , Bravais::P,  7, 133},
  { 27, "C2v^3" , "P 2 -2c"         , "P c c 2"                        , "P c c 2"            , "Pcc2"      , ""     , Bravais::P,  7, 134},
  { 27, "C2v^3" , "P -2a 2"         , "P 2 a a"                        , "P 2 a a"            , "P2aa"      , "cab"  , Bravais::P,  7, 135},
  { 27, "C2v^3" , "P -2b -2b"       , "P b 2 b"                        , "P b 2 b"            , "Pb2b"      , "bca"  , Bravais::P,  7, 136},
  { 28, "C2v^4" , "P 2 -2a"         , "P m a 2"                        , "P m a 2"            , "Pma2"      , ""     , Bravais::P,  7, 137},
  { 28, "C2v^4" , "P 2 -2b"         , "P b m 2"                        , "P b m 2"            , "Pbm2"      , "ba-c" , Bravais::P,  7, 138},
  { 28, "C2v^4" , "P -2b 2"         , "P 2 m b"                        , "P 2 m b"            , "P2mb"      , "cab"  , Bravais::P,  7, 139},
  { 28, "C2v^4" , "P -2c 2"         , "P 2 c m"                        , "P 2 c m"            , "P2cm"      , "-cba" , Bravais::P,  7, 140},
  { 28, "C2v^4" , "P -2c -2c"       , "P c 2 m"                        , "P c 2 m"            , "Pc2m"      , "bca"  , Bravais::P,  7, 141},
  { 28, "C2v^4" , "P -2a -2a"       , "P m 2 a"                        , "P m 2 a"            , "Pm2a"      , "a-cb" , Bravais::P,  7, 142},
  { 29, "C2v^5" , "P 2c -2ac"       , "P c a 2_1"                      , "P c a 2_1"          , "Pca2_1"    , ""     , Bravais::P,  7, 143},
  { 29, "C2v^5" , "P 2c -2b"        , "P b c 2_1"                      , "P b c 2_1"          , "Pbc2_1"    , "ba-c" , Bravais::P,  7, 144},
  { 29, "C2v^5" , "P -2b 2a"        , "P 2_1 a b"                      , "P 2_1 a b"          , "P2_1ab"    , "cab"  , Bravais::P,  7, 145},
  { 29, "C2v^5" , "P -2ac 2a"       , "P 2_1 c a"                      , "P 2_1 c a"          , "P2_1ca"    , "-cba" , Bravais::P,  7, 146},
  { 29, "C2v^5" , "P -2bc -2c"      , "P c 2_1 b"                      , "P c 2_1 b"          , "Pc2_1b"    , "bca"  , Bravais::P,  7, 147},
  { 29, "C2v^5" , "P -2a -2ab"      , "P b 2_1 a"                      , "P b 2_1 a"          , "Pb2_1a"    , "a-cb" , Bravais::P,  7, 148},
  { 30, "C2v^6" , "P 2 -2bc"        , "P n c 2"                        , "P n c 2"            , "Pnc2"      , ""     , Bravais::P,  7, 149},
  { 30, "C2v^6" , "P 2 -2ac"        , "P c n 2"                        , "P c n 2"            , "Pcn2"      , "ba-c" , Bravais::P,  7, 150},
  { 30, "C2v^6" , "P -2ac 2"        , "P 2 n a"                        , "P 2 n a"            , "P2na"      , "cab"  , Bravais::P,  7, 151},
  { 30, "C2v^6" , "P -2ab 2"        , "P 2 a n"                        , "P 2 a n"            , "P2an"      , "-cba" , Bravais::P,  7, 152},
  { 30, "C2v^6" , "P -2ab -2ab"     , "P b 2 n"                        , "P b 2 n"            , "Pb2n"      , "bca"  , Bravais::P,  7, 153},
  { 30, "C2v^6" , "P -2bc -2bc"     , "P n 2 b"                        , "P n 2 b"            , "Pn2b"      , "a-cb" , Bravais::P,  7, 154},
  { 31, "C2v^7" , "P 2ac -2"        , "P m n 2_1"                      , "P m n 2_1"          , "Pmn2_1"    , ""     , Bravais::P,  7, 155},
  { 31, "C2v^7" , "P 2bc -2bc"      , "P n m 2_1"                      , "P n m 2_1"          , "Pnm2_1"    , "ba-c" , Bravais::P,  7, 156},
  { 31, "C2v^7" , "P -2ab 2ab"      , "P 2_1 m n"                      , "P 2_1 m n"          , "P2_1mn"    , "cab"  , Bravais::P,  7, 157},
  { 31, "C2v^7" , "P -2 2ac"        , "P 2_1 n m"                      , "P 2_1 n m"          , "P2_1nm"    , "-cba" , Bravais::P,  7, 158},
  { 31, "C2v^7" , "P -2 -2bc"       , "P n 2_1 m"                      , "P n 2_1 m"          , "Pn2_1m"    , "bca"  , Bravais::P,  7, 159},
  { 31, "C2v^7" , "P -2ab -2"       , "P m 2_1 n"                      , "P m 2_1 n"          , "Pm2_1n"    , "a-cb" , Bravais::P,  7, 160},
  { 32, "C2v^8" , "P 2 -2ab"        , "P b a 2"                        , "P b a 2"            , "Pba2"      , ""     , Bravais::P,  7, 161},
  { 32, "C2v^8" , "P -2bc 2"        , "P 2 c b"                        , "P 2 c b"            , "P2cb"      , "cab"  , Bravais::P,  7, 162},
  { 32, "C2v^8" , "P -2ac -2ac"     , "P c 2 a"                        , "P c 2 a"            , "Pc2a"      , "bca"  , Bravais::P,  7, 163},
  { 33, "C2v^9" , "P 2c -2n"        , "P n a 2_1"                      , "P n a 2_1"          , "Pna2_1"    , ""     , Bravais::P,  7, 164},
  { 33, "C2v^9" , "P 2c -2ab"       , "P b n 2_1"                      , "P b n 2_1"          , "Pbn2_1"    , "ba-c" , Bravais::P,  7, 165},
  { 33, "C2v^9" , "P -2bc 2a"       , "P 2_1 n b"                      , "P 2_1 n b"          , "P2_1nb"    , "cab"  , Bravais::P,  7, 166},
  { 33, "C2v^9" , "P -2n 2a"        , "P 2_1 c n"                      , "P 2_1 c n"          , "P2_1cn"    , "-cba" , Bravais::P,  7, 167},
  { 33, "C2v^9" , "P -2n -2ac"      , "P c 2_1 n"                      , "P c 2_1 n"          , "Pc2_1n"    , "bca"  , Bravais::P,  7, 168},
  { 33, "C2v^9" , "P -2ac -2n"      , "P n 2_1 a"                      , "P n 2_1 a"          , "Pn2_1a"    , "a-cb" , Bravais::P,  7, 169},
  { 34, "C2v^10", "P 2 -2n"         , "P n n 2"                        , "P n n 2"            , "Pnn2"      , ""     , Bravais::P,  7, 170},
  { 34, "C2v^10", "P -2n 2"         , "P 2 n n"                        , "P 2 n n"            , "P2nn"      , "cab"  , Bravais::P,  7, 171},
  { 34, "C2v^10", "P -2n -2n"       , "P n 2 n"                        , "P n 2 n"            , "Pn2n"      , "bca"  , Bravais::P,  7, 172},
  { 35, "C2v^11", "C 2 -2"          , "C m m 2"                        , "C m m 2"            , "Cmm2"      , ""     , Bravais::C,  7, 173},
  { 35, "C2v^11", "A -2 2"          , "A 2 m m"                        , "A 2 m m"            , "A2mm"      , "cab"  , Bravais::A,  7, 174},
  { 35, "C2v^11", "B -2 -2"         , "B m 2 m"                        , "B m 2 m"            , "Bm2m"      , "bca"  , Bravais::P,  7, 175},
  { 36, "C2v^12", "C 2c -2"         , "C m c 2_1"                      , "C m c 2_1"          , "Cmc2_1"    , ""     , Bravais::C,  7, 176},
  { 36, "C2v^12", "C 2c -2c"        , "C c m 2_1"                      , "C c m 2_1"          , "Ccm2_1"    , "ba-c" , Bravais::C,  7, 177},
  { 36, "C2v^12", "A -2a 2a"        , "A 2_1 m a"                      , "A 2_1 m a"          , "A2_1ma"    , "cab"  , Bravais::A,  7, 178},
  { 36, "C2v^12", "A -2 2a"         , "A 2_1 a m"                      , "A 2_1 a m"          , "A2_1am"    , "-cba" , Bravais::A,  7, 179},
  { 36, "C2v^12", "B -2 -2b"        , "B b 2_1 m"                      , "B b 2_1 m"          , "Bb2_1m"    , "bca"  , Bravais::P,  7, 180},
  { 36, "C2v^12", "B -2b -2"        , "B m 2_1 b"                      , "B m 2_1 b"          , "Bm2_1b"    , "a-cb" , Bravais::P,  7, 181},
  { 37, "C2v^13", "C 2 -2c"         , "C c c 2"                        , "C c c 2"            , "Ccc2"      , ""     , Bravais::C,  7, 182},
  { 37, "C2v^13", "A -2a 2"         , "A 2 a a"                        , "A 2 a a"            , "A2aa"      , "cab"  , Bravais::A,  7, 183},
  { 37, "C2v^13", "B -2b -2b"       , "B b 2 b"                        , "B b 2 b"            , "Bb2b"      , "bca"  , Bravais::P,  7, 184},
  { 38, "C2v^14", "A 2 -2"          , "A m m 2"                        , "A m m 2"            , "Amm2"      , ""     , Bravais::A,  7, 185},
  { 38, "C2v^14", "B 2 -2"          , "B m m 2"                        , "B m m 2"            , "Bmm2"      , "ba-c" , Bravais::P,  7, 186},
  { 38, "C2v^14", "B -2 2"          , "B 2 m m"                        , "B 2 m m"            , "B2mm"      , "cab"  , Bravais::P,  7, 187},
  { 38, "C2v^14", "C -2 2"          , "C 2 m m"                        , "C 2 m m"            , "C2mm"      , "-cba" , Bravais::C,  7, 188},
  { 38, "C2v^14", "C -2 -2"         , "C m 2 m"                        , "C m 2 m"            , "Cm2m"      , "bca"  , Bravais::C,  7, 189},
  { 38, "C2v^14", "A -2 -2"         , "A m 2 m"                        , "A m 2 m"            , "Am2m"      , "a-cb" , Bravais::A,  7, 190},
  { 39, "C2v^15", "A 2 -2c"         , "A e m 2"                        , "A e m 2"            , "Aem2"      , ""     , Bravais::A,  7, 191},
  { 39, "C2v^15", "B 2 -2c"         , "B m e 2"                        , "B m e 2"            , "Bme2"      , "ba-c" , Bravais::P,  7, 192},
  { 39, "C2v^15", "B -2c 2"         , "B 2 e m"                        , "B 2 e m"            , "B2em"      , "cab"  , Bravais::P,  7, 193},
  { 39, "C2v^15", "C -2b 2"         , "C 2 m e"                        , "C 2 m e"            , "C2me"      , "-cba" , Bravais::C,  7, 194},
  { 39, "C2v^15", "C -2b -2b"       , "C m 2 e"                        , "C m 2 e"            , "Cm2e"      , "bca"  , Bravais::C,  7, 195},
  { 39, "C2v^15", "A -2c -2c"       , "A e 2 m"                        , "A e 2 m"            , "Ae2m"      , "a-cb" , Bravais::A,  7, 196},
  { 40, "C2v^16", "A 2 -2a"         , "A m a 2"                        , "A m a 2"            , "Ama2"      , ""     , Bravais::A,  7, 197},
  { 40, "C2v^16", "B 2 -2b"         , "B b m 2"                        , "B b m 2"            , "Bbm2"      , "ba-c" , Bravais::P,  7, 198},
  { 40, "C2v^16", "B -2b 2"         , "B 2 m b"                        , "B 2 m b"            , "B2mb"      , "cab"  , Bravais::P,  7, 199},
  { 40, "C2v^16", "C -2c 2"         , "C 2 c m"                        , "C 2 c m"            , "C2cm"      , "-cba" , Bravais::C,  7, 200},
  { 40, "C2v^16", "C -2c -2c"       , "C c 2 m"                        , "C c 2 m"            , "Cc2m"      , "bca"  , Bravais::C,  7, 201},
  { 40, "C2v^16", "A -2a -2a"       , "A m 2 a"                        , "A m 2 a"            , "Am2a"      , "a-cb" , Bravais::A,  7, 202},
  { 41, "C2v^17", "A 2 -2ac"        , "A e a 2"                        , "A e a 2"            , "Aea2"      , ""     , Bravais::A,  7, 203},
  { 41, "C2v^17", "B 2 -2bc"        , "B b e 2"                        , "B b e 2"            , "Bbe2"      , "ba-c" , Bravais::P,  7, 204},
  { 41, "C2v^17", "B -2bc 2"        , "B 2 e b"                        , "B 2 e b"            , "B2eb"      , "cab"  , Bravais::P,  7, 205},
  { 41, "C2v^17", "C -2bc 2"        , "C 2 c e"                        , "C 2 c e"            , "C2ce"      , "-cba" , Bravais::C,  7, 206},
  { 41, "C2v^17", "C -2bc -2bc"     , "C c 2 e"                        , "C c 2 e"            , "Cc2e"      , "bca"  , Bravais::C,  7, 207},
  { 41, "C2v^17", "A -2ac -2ac"     , "A e 2 a"                        , "A e 2 a"            , "Ae2a"      , "a-cb" , Bravais::A,  7, 208},
  { 42, "C2v^18", "F 2 -2"          , "F m m 2"                        , "F m m 2"            , "Fmm2"      , ""     , Bravais::F,  7, 209},
  { 42, "C2v^18", "F -2 2"          , "F 2 m m"                        , "F 2 m m"            , "F2mm"      , "cab"  , Bravais::F,  7, 210},
  { 42, "C2v^18", "F -2 -2"         , "F m 2 m"                        , "F m 2 m"            , "Fm2m"      , "bca"  , Bravais::F,  7, 211},
  { 43, "C2v^19", "F 2 -2d"         , "F d d 2"                        , "F d d 2"            , "Fdd2"      , ""     , Bravais::F,  7, 212},
  { 43, "C2v^19", "F -2d 2"         , "F 2 d d"                        , "F 2 d d"            , "F2dd"      , "cab"  , Bravais::F,  7, 213},
  { 43, "C2v^19", "F -2d -2d"       , "F d 2 d"                        , "F d 2 d"            , "Fd2d"      , "bca"  , Bravais::F,  7, 214},
  { 44, "C2v^20", "I 2 -2"          , "I m m 2"                        , "I m m 2"            , "Imm2"      , ""     , Bravais::I,  7, 215},
  { 44, "C2v^20", "I -2 2"          , "I 2 m m"                        , "I 2 m m"            , "I2mm"      , "cab"  , Bravais::I,  7, 216},
  { 44, "C2v^20", "I -2 -2"         , "I m 2 m"                        , "I m 2 m"            , "Im2m"      , "bca"  , Bravais::I,  7, 217},
  { 45, "C2v^21", "I 2 -2c"         , "I b a 2"                        , "I b a 2"            , "Iba2"      , ""     , Bravais::I,  7, 218},
  { 45, "C2v^21", "I -2a 2"         , "I 2 c b"                        , "I 2 c b"            , "I2cb"      , "cab"  , Bravais::I,  7, 219},
  { 45, "C2v^21", "I -2b -2b"       , "I c 2 a"                        , "I c 2 a"            , "Ic2a"      , "bca"  , Bravais::I,  7, 220},
  { 46, "C2v^22", "I 2 -2a"         , "I m a 2"                        , "I m a 2"            , "Ima2"      , ""     , Bravais::I,  7, 221},
  { 46, "C2v^22", "I 2 -2b"         , "I b m 2"                        , "I b m 2"            , "Ibm2"      , "ba-c" , Bravais::I,  7, 222},
  { 46, "C2v^22", "I -2b 2"         , "I 2 m b"                        , "I 2 m b"            , "I2mb"      , "cab"  , Bravais::I,  7, 223},
  { 46, "C2v^22", "I -2c 2"         , "I 2 c m"                        , "I 2 c m"            , "I2cm"      , "-cba" , Bravais::I,  7, 224},
  { 46, "C2v^22", "I -2c -2c"       , "I c 2 m"                        , "I c 2 m"            , "Ic2m"      , "bca"  , Bravais::I,  7, 225},
  { 46, "C2v^22", "I -2a -2a"       , "I m 2 a"                        , "I m 2 a"            , "Im2a"      , "a-cb" , Bravais::I,  7, 226},
  { 47, "D2h^1" , "-P 2 2"          , "P m m m"                        , "P 2/m 2/m 2/m"      , "Pmmm"      , ""     , Bravais::P,  8, 227},
  { 48, "D2h^2" , "P 2 2 -1n"       , "P n n n"                        , "P 2/n 2/n 2/n"      , "Pnnn"      , "1"    , Bravais::P,  8, 228},
  { 48, "D2h^2" , "-P 2ab 2bc"      , "P n n n"                        , "P 2/n 2/n 2/n"      , "Pnnn"      , "2"    , Bravais::P,  8, 229},
  { 49, "D2h^3" , "-P 2 2c"         , "P c c m"                        , "P 2/c 2/c 2/m"      , "Pccm"      , ""     , Bravais::P,  8, 230},
  { 49, "D2h^3" , "-P 2a 2"         , "P m a a"                        , "P 2/m 2/a 2/a"      , "Pmaa"      , "cab"  , Bravais::P,  8, 231},
  { 49, "D2h^3" , "-P 2b 2b"        , "P b m b"                        , "P 2/b 2/m 2/b"      , "Pbmb"      , "bca"  , Bravais::P,  8, 232},
  { 50, "D2h^4" , "P 2 2 -1ab"      , "P b a n"                        , "P 2/b 2/a 2/n"      , "Pban"      , "1"    , Bravais::P,  8, 233},
  { 50, "D2h^4" , "-P 2ab 2b"       , "P b a n"                        , "P 2/b 2/a 2/n"      , "Pban"      , "2"    , Bravais::P,  8, 234},
  { 50, "D2h^4" , "P 2 2 -1bc"      , "P n c b"                        , "P 2/n 2/c 2/b"      , "Pncb"      , "1cab" , Bravais::P,  8, 235},
  { 50, "D2h^4" , "-P 2b 2bc"       , "P n c b"                        , "P 2/n 2/c 2/b"      , "Pncb"      , "2cab" , Bravais::P,  8, 236},
  { 50, "D2h^4" , "P 2 2 -1ac"      , "P c n a"                        , "P 2/c 2/n 2/a"      , "Pcna"      , "1bca" , Bravais::P,  8, 237},
  { 50, "D2h^4" , "-P 2a 2c"        , "P c n a"                        , "P 2/c 2/n 2/a"      , "Pcna"      , "2bca" , Bravais::P,  8, 238},
  { 51, "D2h^5" , "-P 2a 2a"        , "P m m a"                        , "P 2_1/m 2/m 2/a"    , "Pmma"      , ""     , Bravais::P,  8, 239},
  { 51, "D2h^5" , "-P 2b 2"         , "P m m b"                        , "P 2/m 2_1/m 2/b"    , "Pmmb"      , "ba-c" , Bravais::P,  8, 240},
  { 51, "D2h^5" , "-P 2 2b"         , "P b m m"                        , "P 2/b 2_1/m 2/m"    , "Pbmm"      , "cab"  , Bravais::P,  8, 241},
  { 51, "D2h^5" , "-P 2c 2c"        , "P c m m"                        , "P 2/c 2/m 2_1/m"    , "Pcmm"      , "-cba" , Bravais::P,  8, 242},
  { 51, "D2h^5" , "-P 2c 2"         , "P m c m"                        , "P 2/m 2/c 2_1/m"    , "Pmcm"      , "bca"  , Bravais::P,  8, 243},
  { 51, "D2h^5" , "-P 2 2a"         , "P m a m"                        , "P 2_1/m 2/a 2/m"    , "Pmam"      , "a-cb" , Bravais::P,  8, 244},
  { 52, "D2h^6" , "-P 2a 2bc"       , "P n n a"                        , "P 2/n 2_1/n 2/a"    , "Pnna"      , ""     , Bravais::P,  8, 245},
  { 52, "D2h^6" , "-P 2b 2n"        , "P n n b"                        , "P 2_1/n 2/n 2/b"    , "Pnnb"      , "ba-c" , Bravais::P,  8, 246},
  { 52, "D2h^6" , "-P 2n 2b"        , "P b n n"                        , "P 2/b 2/n 2_1/n"    , "Pbnn"      , "cab"  , Bravais::P,  8, 247},
  { 52, "D2h^6" , "-P 2ab 2c"       , "P c n n"                        , "P 2/c 2_1/n 2/n"    , "Pcnn"      , "-cba" , Bravais::P,  8, 248},
  { 52, "D2h^6" , "-P 2ab 2n"       , "P n c n"                        , "P 2_1/n 2/c 2/n"    , "Pncn"      , "bca"  , Bravais::P,  8, 249},
  { 52, "D2h^6" , "-P 2n 2bc"       , "P n a n"                        , "P 2/n 2/a 2_1/n"    , "Pnan"      , "a-cb" , Bravais::P,  8, 250},
  { 53, "D2h^7" , "-P 2ac 2"        , "P m n a"                        , "P 2/m 2/n 2_1/a"    , "Pmna"      , ""     , Bravais::P,  8, 251},
  { 53, "D2h^7" , "-P 2bc 2bc"      , "P n m b"                        , "P 2/n 2/m 2_1/b"    , "Pnmb"      , "ba-c" , Bravais::P,  8, 252},
  { 53, "D2h^7" , "-P 2ab 2ab"      , "P b m n"                        , "P 2_1/b 2/m 2/n"    , "Pbmn"      , "cab"  , Bravais::P,  8, 253},
  { 53, "D2h^7" , "-P 2 2ac"        , "P c n m"                        , "P 2_1/c 2/n 2/m"    , "Pcnm"      , "-cba" , Bravais::P,  8, 254},
  { 53, "D2h^7" , "-P 2 2bc"        , "P n c m"                        , "P 2/n 2_1/c 2/m"    , "Pncm"      , "bca"  , Bravais::P,  8, 255},
  { 53, "D2h^7" , "-P 2ab 2"        , "P m a n"                        , "P 2/m 2_1/a 2/n"    , "Pman"      , "a-cb" , Bravais::P,  8, 256},
  { 54, "D2h^8" , "-P 2a 2ac"       , "P c c a"                        , "P 2_1/c 2/c 2/a"    , "Pcca"      , ""     , Bravais::P,  8, 257},
  { 54, "D2h^8" , "-P 2b 2c"        , "P c c b"                        , "P 2/c 2_1/c 2/b"    , "Pccb"      , "ba-c" , Bravais::P,  8, 258},
  { 54, "D2h^8" , "-P 2a 2b"        , "P b a a"                        , "P 2/b 2_1/a 2/a"    , "Pbaa"      , "cab"  , Bravais::P,  8, 259},
  { 54, "D2h^8" , "-P 2ac 2c"       , "P c a a"                        , "P 2/c 2/a 2_1/a"    , "Pcaa"      , "-cba" , Bravais::P,  8, 260},
  { 54, "D2h^8" , "-P 2bc 2b"       , "P b c b"                        , "P 2/b 2/c 2_1/b"    , "Pbcb"      , "bca"  , Bravais::P,  8, 261},
  { 54, "D2h^8" , "-P 2b 2ab"       , "P b a b"                        , "P 2_1/b 2/a 2/b"    , "Pbab"      , "a-cb" , Bravais::P,  8, 262},
  { 55, "D2h^9" , "-P 2 2ab"        , "P b a m"                        , "P 2_1/b 2_1/a 2/m"  , "Pbam"      , ""     , Bravais::P,  8, 263},
  { 55, "D2h^9" , "-P 2bc 2"        , "P m c b"                        , "P 2/m 2_1/c 2_1/b"  , "Pmcb"      , "cab"  , Bravais::P,  8, 264},
  { 55, "D2h^9" , "-P 2ac 2ac"      , "P c m a"                        , "P 2_1/c 2/m 2_1/a"  , "Pcma"      , "bca"  , Bravais::P,  8, 265},
  { 56, "D2h^10", "-P 2ab 2ac"      , "P c c n"                        , "P 2_1/c 2_1/c 2/n"  , "Pccn"      , ""     , Bravais::P,  8, 266},
  { 56, "D2h^10", "-P 2ac 2bc"      , "P n a a"                        , "P 2/n 2_1/a 2_1/a"  , "Pnaa"      , "cab"  , Bravais::P,  8, 267},
  { 56, "D2h^10", "-P 2bc 2ab"      , "P b n b"                        , "P 2_1/b 2/n 2_1/b"  , "Pbnb"      , "bca"  , Bravais::P,  8, 268},
  { 57, "D2h^11", "-P 2c 2b"        , "P b c m"                        , "P 2/b 2_1/c 2_1/m"  , "Pbcm"      , ""     , Bravais::P,  8, 269},
  { 57, "D2h^11", "-P 2c 2ac"       , "P c a m"                        , "P 2_1/c 2/a 2_1/m"  , "Pcam"      , "ba-c" , Bravais::P,  8, 270},
  { 57, "D2h^11", "-P 2ac 2a"       , "P m c a"                        , "P 2_1/m 2/c 2_1/a"  , "Pmca"      , "cab"  , Bravais::P,  8, 271},
  { 57, "D2h^11", "-P 2b 2a"        , "P m a b"                        , "P 2_1/m 2_1/a 2/b"  , "Pmab"      , "-cba" , Bravais::P,  8, 272},
  { 57, "D2h^11", "-P 2a 2ab"       , "P b m a"                        , "P 2_1/b 2_1/m 2/a"  , "Pbma"      , "bca"  , Bravais::P,  8, 273},
  { 57, "D2h^11", "-P 2bc 2c"       , "P c m b"                        , "P 2/c 2_1/m 2_1/b"  , "Pcmb"      , "a-cb" , Bravais::P,  8, 274},
  { 58, "D2h^12", "-P 2 2n"         , "P n n m"                        , "P 2_1/n 2_1/n 2/m"  , "Pnnm"      , ""     , Bravais::P,  8, 275},
  { 58, "D2h^12", "-P 2n 2"         , "P m n n"                        , "P 2/m 2_1/n 2_1/n"  , "Pmnn"      , "cab"  , Bravais::P,  8, 276},
  { 58, "D2h^12", "-P 2n 2n"        , "P n m n"                        , "P 2_1/n 2/m 2_1/n"  , "Pnmn"      , "bca"  , Bravais::P,  8, 277},
  { 59, "D2h^13", "P 2 2ab -1ab"    , "P m m n"                        , "P 2_1/m 2_1/m 2/n"  , "Pmmn"      , "1"    , Bravais::P,  8, 278},
  { 59, "D2h^13", "-P 2ab 2a"       , "P m m n"                        , "P 2_1/m 2_1/m 2/n"  , "Pmmn"      , "2"    , Bravais::P,  8, 279},
  { 59, "D2h^13", "P 2bc 2 -1bc"    , "P n m m"                        , "P 2/n 2_1/m 2_1/m"  , "Pnmm"      , "1cab" , Bravais::P,  8, 280},
  { 59, "D2h^13", "-P 2c 2bc"       , "P n m m"                        , "P 2/n 2_1/m 2_1/m"  , "Pnmm"      , "2cab" , Bravais::P,  8, 281},
  { 59, "D2h^13", "P 2ac 2ac -1ac"  , "P m n m"                        , "P 2_1/m 2/n 2_1/m"  , "Pmnm"      , "1bca" , Bravais::P,  8, 282},
  { 59, "D2h^13", "-P 2c 2a"        , "P m n m"                        , "P 2_1/m 2/n 2_1/m"  , "Pmnm"      , "2bca" , Bravais::P,  8, 283},
  { 60, "D2h^14", "-P 2n 2ab"       , "P b c n"                        , "P 2_1/b 2/c 2_1/n"  , "Pbcn"      , ""     , Bravais::P,  8, 284},
  { 60, "D2h^14", "-P 2n 2c"        , "P c a n"                        , "P 2/c 2_1/a 2_1/n"  , "Pcan"      , "ba-c" , Bravais::P,  8, 285},
  { 60, "D2h^14", "-P 2a 2n"        , "P n c a"                        , "P 2_1/n 2_1/c 2/a"  , "Pnca"      , "cab"  , Bravais::P,  8, 286},
  { 60, "D2h^14", "-P 2bc 2n"       , "P n a b"                        , "P 2_1/n 2/a 2_1/b"  , "Pnab"      , "-cba" , Bravais::P,  8, 287},
  { 60, "D2h^14", "-P 2ac 2b"       , "P b n a"                        , "P 2/b 2_1/n 2_1/a"  , "Pbna"      , "bca"  , Bravais::P,  8, 288},
  { 60, "D2h^14", "-P 2b 2ac"       , "P c n b"                        , "P 2_1/c 2_1/n 2/b"  , "Pcnb"      , "a-cb" , Bravais::P,  8, 289},
  { 61, "D2h^15", "-P 2ac 2ab"      , "P b c a"                        , "P 2_1/b 2_1/c 2_1/a", "Pbca"      , ""     , Bravais::P,  8, 290},
  { 61, "D2h^15", "-P 2bc 2ac"      , "P c a b"                        , "P 2_1/c 2_1/a 2_1/b", "Pcab"      , "ba-c" , Bravais::P,  8, 291},
  { 62, "D2h^16", "-P 2ac 2n"       , "P n m a"                        , "P 2_1/n 2_1/m 2_1/a", "Pnma"      , ""     , Bravais::P,  8, 292},
  { 62, "D2h^16", "-P 2bc 2a"       , "P m n b"                        , "P 2_1/m 2_1/n 2_1/b", "Pmnb"      , "ba-c" , Bravais::P,  8, 293},
  { 62, "D2h^16", "-P 2c 2ab"       , "P b n m"                        , "P 2_1/b 2_1/n 2_1/m", "Pbnm"      , "cab"  , Bravais::P,  8, 294},
  { 62, "D2h^16", "-P 2n 2ac"       , "P c m n"                        , "P 2_1/c 2_1/m 2_1/n", "Pcmn"      , "-cba" , Bravais::P,  8, 295},
  { 62, "D2h^16", "-P 2n 2a"        , "P m c n"                        , "P 2_1/m 2_1/c 2_1/n", "Pmcn"      , "bca"  , Bravais::P,  8, 296},
  { 62, "D2h^16", "-P 2c 2n"        , "P n a m"                        , "P 2_1/n 2_1/a 2_1/m", "Pnam"      , "a-cb" , Bravais::P,  8, 297},
  { 63, "D2h^17", "-C 2c 2"         , "C m c m"                        , "C 2/m 2/c 2_1/m"    , "Cmcm"      , ""     , Bravais::C,  8, 298},
  { 63, "D2h^17", "-C 2c 2c"        , "C c m m"                        , "C 2/c 2/m 2_1/m"    , "Ccmm"      , "ba-c" , Bravais::C,  8, 299},
  { 63, "D2h^17", "-A 2a 2a"        , "A m m a"                        , "A 2_1/m 2/m 2/a"    , "Amma"      , "cab"  , Bravais::A,  8, 300},
  { 63, "D2h^17", "-A 2 2a"         , "A m a m"                        , "A 2_1/m 2/a 2/m"    , "Amam"      , "-cba" , Bravais::A,  8, 301},
  { 63, "D2h^17", "-B 2 2b"         , "B b m m"                        , "B 2/b 2_1/m 2/m"    , "Bbmm"      , "bca"  , Bravais::P,  8, 302},
  { 63, "D2h^17", "-B 2b 2"         , "B m m b"                        , "B 2/m 2_1/m 2/b"    , "Bmmb"      , "a-cb" , Bravais::P,  8, 303},
  { 64, "D2h^18", "-C 2bc 2"        , "C m c e"                        , "C 2/m 2/c 2_1/e"    , "Cmce"      , ""     , Bravais::C,  8, 304},
  { 64, "D2h^18", "-C 2bc 2bc"      , "C c m e"                        , "C 2/c 2/m 2_1/e"    , "Ccme"      , "ba-c" , Bravais::C,  8, 305},
  { 64, "D2h^18", "-A 2ac 2ac"      , "A e m a"                        , "A 2_1/e 2/m 2/a"    , "Aema"      , "cab"  , Bravais::A,  8, 306},
  { 64, "D2h^18", "-A 2 2ac"        , "A e a m"                        , "A 2_1/e 2/a 2/m"    , "Aeam"      , "-cba" , Bravais::A,  8, 307},
  { 64, "D2h^18", "-B 2 2bc"        , "B b e m"                        , "B 2/b 2_1/e 2/m"    , "Bbem"      , "bca"  , Bravais::P,  8, 308},
  { 64, "D2h^18", "-B 2bc 2"        , "B m e b"                        , "B 2/m 2_1/e 2/b"    , "Bmeb"      , "a-cb" , Bravais::P,  8, 309},
  { 65, "D2h^19", "-C 2 2"          , "C m m m"                        , "C 2/m 2/m 2/m"      , "Cmmm"      , ""     , Bravais::C,  8, 310},
  { 65, "D2h^19", "-A 2 2"          , "A m m m"                        , "A 2/m 2/m 2/m"      , "Ammm"      , "cab"  , Bravais::A,  8, 311},
  { 65, "D2h^19", "-B 2 2"          , "B m m m"                        , "B 2/m 2/m 2/m"      , "Bmmm"      , "bca"  , Bravais::P,  8, 312},
  { 66, "D2h^20", "-C 2 2c"         , "C c c m"                        , "C 2/c 2/c 2/m"      , "Cccm"      , ""     , Bravais::C,  8, 313},
  { 66, "D2h^20", "-A 2a 2"         , "A m a a"                        , "A 2/m 2/a 2/a"      , "Amaa"      , "cab"  , Bravais::A,  8, 314},
  { 66, "D2h^20", "-B 2b 2b"        , "B b m b"                        , "B 2/b 2/m 2/b"      , "Bbmb"      , "bca"  , Bravais::P,  8, 315},
  { 67, "D2h^21", "-C 2b 2"         , "C m m e"                        , "C 2/m 2/m 2/e"      , "Cmme"      , ""     , Bravais::C,  8, 316},
  { 67, "D2h^21", "-C 2b 2b"        , "C m m e"                        , "C 2/m 2/m 2/e"      , "Cmme"      , "ba-c" , Bravais::C,  8, 317},
  { 67, "D2h^21", "-A 2c 2c"        , "A e m m"                        , "A 2/e 2/m 2/m"      , "Aemm"      , "cab"  , Bravais::A,  8, 318},
  { 67, "D2h^21", "-A 2 2c"         , "A e m m"                        , "A 2/e 2/m 2/m"      , "Aemm"      , "-cba" , Bravais::A,  8, 319},
  { 67, "D2h^21", "-B 2 2c"         , "B m e m"                        , "B 2/m 2/e 2/m"      , "Bmem"      , "bca"  , Bravais::P,  8, 320},
  { 67, "D2h^21", "-B 2c 2"         , "B m e m"                        , "B 2/m 2/e 2/m"      , "Bmem"      , "a-cb" , Bravais::P,  8, 321},
  { 68, "D2h^22", "C 2 2 -1bc"      , "C c c e"                        , "C 2/c 2/c 2/e"      , "Ccce"      , "1"    , Bravais::C,  8, 322},
  { 68, "D2h^22", "-C 2b 2bc"       , "C c c e"                        , "C 2/c 2/c 2/e"      , "Ccce"      , "2"    , Bravais::C,  8, 323},
  { 68, "D2h^22", "C 2 2 -1bc"      , "C c c e"                        , "C 2/c 2/c 2/e"      , "Ccce"      , "1ba-c", Bravais::C,  8, 324},
  { 68, "D2h^22", "-C 2b 2c"        , "C c c e"                        , "C 2/c 2/c 2/e"      , "Ccce"      , "2ba-c", Bravais::C,  8, 325},
  { 68, "D2h^22", "A 2 2 -1ac"      , "A e a a"                        , "A 2/e 2/a 2/a"      , "Aeaa"      , "1cab" , Bravais::A,  8, 326},
  { 68, "D2h^22", "-A 2a 2c"        , "A e a a"                        , "A 2/e 2/a 2/a"      , "Aeaa"      , "2cab" , Bravais::A,  8, 327},
  { 68, "D2h^22", "A 2 2 -1ac"      , "A e a a"                        , "A 2/e 2/a 2/a"      , "Aeaa"      , "1-cba", Bravais::A,  8, 328},
  { 68, "D2h^22", "-A 2ac 2c"       , "A e a a"                        , "A 2/e 2/a 2/a"      , "Aeaa"      , "2-cba", Bravais::A,  8, 329},
  { 68, "D2h^22", "B 2 2 -1bc"      , "B b e b"                        , "B 2/b 2/e 2/b"      , "Bbeb"      , "1bca" , Bravais::P,  8, 330},
  { 68, "D2h^22", "-B 2bc 2b"       , "B b c b"                        , "B 2/b 2/e 2/b"      , "Bbcb"      , "2bca" , Bravais::P,  8, 331},
  { 68, "D2h^22", "B 2 2 -1bc"      , "B b e b"                        , "B 2/b 2/e 2/b"      , "Bbeb"      , "1a-cb", Bravais::P,  8, 332},
  { 68, "D2h^22", "-B 2b 2bc"       , "B b e b"                        , "B 2/b 2/e 2/b"      , "Bbeb"      , "2a-cb", Bravais::P,  8, 333},
  { 69, "D2h^23", "-F 2 2"          , "F m m m"                        , "F 2/m 2/m 2/m"      , "Fmmm"      , ""     , Bravais::F,  8, 334},
  { 70, "D2h^24", "F 2 2 -1d"       , "F d d d"                        , "F 2/d 2/d 2/d"      , "Fddd"      , "1"    , Bravais::F,  8, 335},
  { 70, "D2h^24", "-F 2uv 2vw"      , "F d d d"                        , "F 2/d 2/d 2/d"      , "Fddd"      , "2"    , Bravais::F,  8, 336},
  { 71, "D2h^25", "-I 2 2"          , "I m m m"                        , "I 2/m 2/m 2/m"      , "Immm"      , ""     , Bravais::I,  8, 337},
  { 72, "D2h^26", "-I 2 2c"         , "I b a m"                        , "I 2/b 2/a 2/m"      , "Ibam"      , ""     , Bravais::I,  8, 338},
  { 72, "D2h^26", "-I 2a 2"         , "I m c b"                        , "I 2/m 2/c 2/b"      , "Imcb"      , "cab"  , Bravais::I,  8, 339},
  { 72, "D2h^26", "-I 2b 2b"        , "I c m a"                        , "I 2/c 2/m 2/a"      , "Icma"      , "bca"  , Bravais::I,  8, 340},
  { 73, "D2h^27", "-I 2b 2c"        , "I b c a"                        , "I 2/b 2/c 2/a"      , "Ibca"      , ""     , Bravais::I,  8, 341},
  { 73, "D2h^27", "-I 2a 2b"        , "I c a b"                        , "I 2/c 2/a 2/b"      , "Icab"      , "ba-c" , Bravais::I,  8, 342},
  { 74, "D2h^28", "-I 2b 2"         , "I m m a"                        , "I 2/m 2/m 2/a"      , "Imma"      , ""     , Bravais::I,  8, 343},
  { 74, "D2h^28", "-I 2a 2a"        , "I m m b"                        , "I 2/m 2/m 2/b"      , "Immb"      , "ba-c" , Bravais::I,  8, 344},
  { 74, "D2h^28", "-I 2c 2c"        , "I b m m"                        , "I 2/b 2/m 2/m"      , "Ibmm"      , "cab"  , Bravais::I,  8, 345},
  { 74, "D2h^28", "-I 2 2b"         , "I c m m"                        , "I 2/c 2/m 2/m"      , "Icmm"      , "-cba" , Bravais::I,  8, 346},
  { 74, "D2h^28", "-I 2 2a"         , "I m c m"                        , "I 2/m 2/c 2/m"      , "Imcm"      , "bca"  , Bravais::I,  8, 347},
  { 74, "D2h^28", "-I 2c 2"         , "I m a m"                        , "I 2/m 2/a 2/m"      , "Imam"      , "a-cb" , Bravais::I,  8, 348},
  { 75, "C4^1"  , "P 4"             , "P 4"                            , "P 4"                , "P4"        , ""     , Bravais::P,  9, 349},
  { 76, "C4^2"  , "P 4w"            , "P 4_1"                          , "P 4_1"              , "P4_1"      , ""     , Bravais::P,  9, 350},
  { 77, "C4^3"  , "P 4c"            , "P 4_2"                          , "P 4_2"              , "P4_2"      , ""     , Bravais::P,  9, 351},
  { 78, "C4^4"  , "P 4cw"           , "P 4_3"                          , "P 4_3"              , "P4_3"      , ""     , Bravais::P,  9, 352},
  { 79, "C4^5"  , "I 4"             , "I 4"                            , "I 4"                , "I4"        , ""     , Bravais::I,  9, 353},
  { 80, "C4^6"  , "I 4bw"           , "I 4_1"                          , "I 4_1"              , "I4_1"      , ""     , Bravais::I,  9, 354},
  { 81, "S4^1"  , "P -4"            , "P -4"                           , "P -4"               , "P-4"       , ""     , Bravais::P, 10, 355},
  { 82, "S4^2"  , "I -4"            , "I -4"                           , "I -4"               , "I-4"       , ""     , Bravais::I, 10, 356},
  { 83, "C4h^1" , "-P 4"            , "P 4/m"                          , "P 4/m"              , "P4/m"      , ""     , Bravais::P, 11, 357},
  { 84, "C4h^2" , "-P 4c"           , "P 4_2/m"                        , "P 4_2/m"            , "P4_2/m"    , ""     , Bravais::P, 11, 358},
  { 85, "C4h^3" , "P 4ab -1ab"      , "P 4/n"                          , "P 4/n"              , "P4/n"      , "1"    , Bravais::P, 11, 359},
  { 85, "C4h^3" , "-P 4a"           , "P 4/n"                          , "P 4/n"              , "P4/n"      , "2"    , Bravais::P, 11, 360},
  { 86, "C4h^4" , "P 4n -1n"        , "P 4_2/n"                        , "P 4_2/n"            , "P4_2/n"    , "1"    , Bravais::P, 11, 361},
  { 86, "C4h^4" , "-P 4bc"          , "P 4_2/n"                        , "P 4_2/n"            , "P4_2/n"    , "2"    , Bravais::P, 11, 362},
  { 87, "C4h^5" , "-I 4"            , "I 4/m"                          , "I 4/m"              , "I4/m"      , ""     , Bravais::I, 11, 363},
  { 88, "C4h^6" , "I 4bw -1bw"      , "I 4_1/a"                        , "I 4_1/a"            , "I4_1/a"    , "1"    , Bravais::I, 11, 364},
  { 88, "C4h^6" , "-I 4ad"          , "I 4_1/a"                        , "I 4_1/a"            , "I4_1/a"    , "2"    , Bravais::I, 11, 365},
  { 89, "D4^1"  , "P 4 2"           , "P 4 2 2"                        , "P 4 2 2"            , "P422"      , ""     , Bravais::P, 12, 366},
  { 90, "D4^2"  , "P 4ab 2ab"       , "P 4 2_1 2"                      , "P 4 2_1 2"          , "P42_12"    , ""     , Bravais::P, 12, 367},
  { 91, "D4^3"  , "P 4w 2c"         , "P 4_1 2 2"                      , "P 4_1 2 2"          , "P4_122"    , ""     , Bravais::P, 12, 368},
  { 92, "D4^4"  , "P 4abw 2nw"      , "P 4_1 2_1 2"                    , "P 4_1 2_1 2"        , "P4_12_12"  , ""     , Bravais::P, 12, 369},
  { 93, "D4^5"  , "P 4c 2"          , "P 4_2 2 2"                      , "P 4_2 2 2"          , "P4_222"    , ""     , Bravais::P, 12, 370},
  { 94, "D4^6"  , "P 4n 2n"         , "P 4_2 2_1 2"                    , "P 4_2 2_1 2"        , "P4_22_12"  , ""     , Bravais::P, 12, 371},
  { 95, "D4^7"  , "P 4cw 2c"        , "P 4_3 2 2"                      , "P 4_3 2 2"          , "P4_322"    , ""     , Bravais::P, 12, 372},
  { 96, "D4^8"  , "P 4nw 2abw"      , "P 4_3 2_1 2"                    , "P 4_3 2_1 2"        , "P4_32_12"  , ""     , Bravais::P, 12, 373},
  { 97, "D4^9"  , "I 4 2"           , "I 4 2 2"                        , "I 4 2 2"            , "I422"      , ""     , Bravais::I, 12, 374},
  { 98, "D4^10" , "I 4bw 2bw"       , "I 4_1 2 2"                      , "I 4_1 2 2"          , "I4_122"    , ""     , Bravais::I, 12, 375},
  { 99, "C4v^1" , "P 4 -2"          , "P 4 m m"                        , "P 4 m m"            , "P4mm"      , ""     , Bravais::P, 13, 376},
  {100, "C4v^2" , "P 4 -2ab"        , "P 4 b m"                        , "P 4 b m"            , "P4bm"      , ""     , Bravais::P, 13, 377},
  {101, "C4v^3" , "P 4c -2c"        , "P 4_2 c m"                      , "P 4_2 c m"          , "P4_2cm"    , ""     , Bravais::P, 13, 378},
  {102, "C4v^4" , "P 4n -2n"        , "P 4_2 n m"                      , "P 4_2 n m"          , "P4_2nm"    , ""     , Bravais::P, 13, 379},
  {103, "C4v^5" , "P 4 -2c"         , "P 4 c c"                        , "P 4 c c"            , "P4cc"      , ""     , Bravais::P, 13, 380},
  {104, "C4v^6" , "P 4 -2n"         , "P 4 n c"                        , "P 4 n c"            , "P4nc"      , ""     , Bravais::P, 13, 381},
  {105, "C4v^7" , "P 4c -2"         , "P 4_2 m c"                      , "P 4_2 m c"          , "P4_2mc"    , ""     , Bravais::P, 13, 382},
  {106, "C4v^8" , "P 4c -2ab"       , "P 4_2 b c"                      , "P 4_2 b c"          , "P4_2bc"    , ""     , Bravais::P, 13, 383},
  {107, "C4v^9" , "I 4 -2"          , "I 4 m m"                        , "I 4 m m"            , "I4mm"      , ""     , Bravais::I, 13, 384},
  {108, "C4v^10", "I 4 -2c"         , "I 4 c m"                        , "I 4 c m"            , "I4cm"      , ""     , Bravais::I, 13, 385},
  {109, "C4v^11", "I 4bw -2"        , "I 4_1 m d"                      , "I 4_1 m d"          , "I4_1md"    , ""     , Bravais::I, 13, 386},
  {110, "C4v^12", "I 4bw -2c"       , "I 4_1 c d"                      , "I 4_1 c d"          , "I4_1cd"    , ""     , Bravais::I, 13, 387},
  {111, "D2d^1" , "P -4 2"          , "P -4 2 m"                       , "P -4 2 m"           , "P-42m"     , ""     , Bravais::P, 14, 388},
  {112, "D2d^2" , "P -4 2c"         , "P -4 2 c"                       , "P -4 2 c"           , "P-42c"     , ""     , Bravais::P, 14, 389},
  {113, "D2d^3" , "P -4 2ab"        , "P -4 2_1 m"                     , "P -4 2_1 m"         , "P-42_1m"   , ""     , Bravais::P, 14, 390},
  {114, "D2d^4" , "P -4 2n"         , "P -4 2_1 c"                     , "P -4 2_1 c"         , "P-42_1c"   , ""     , Bravais::P, 14, 391},
  {115, "D2d^5" , "P -4 -2"         , "P -4 m 2"                       , "P -4 m 2"           , "P-4m2"     , ""     , Bravais::P, 14, 392},
  {116, "D2d^6" , "P -4 -2c"        , "P -4 c 2"                       , "P -4 c 2"           , "P-4c2"     , ""     , Bravais::P, 14, 393},
  {117, "D2d^7" , "P -4 -2ab"       , "P -4 b 2"                       , "P -4 b 2"           , "P-4b2"     , ""     , Bravais::P, 14, 394},
  {118, "D2d^8" , "P -4 -2n"        , "P -4 n 2"                       , "P -4 n 2"           , "P-4n2"     , ""     , Bravais::P, 14, 395},
  {119, "D2d^9" , "I -4 -2"         , "I -4 m 2"                       , "I -4 m 2"           , "I-4m2"     , ""     , Bravais::I, 14, 396},
  {120, "D2d^10", "I -4 -2c"        , "I -4 c 2"                       , "I -4 c 2"           , "I-4c2"     , ""     , Bravais::I, 14, 397},
  {121, "D2d^11", "I -4 2"          , "I -4 2 m"                       , "I -4 2 m"           , "I-42m"     , ""     , Bravais::I, 14, 398},
  {122, "D2d^12", "I -4 2bw"        , "I -4 2 d"                       , "I -4 2 d"           , "I-42d"     , ""     , Bravais::I, 14, 399},
  {123, "D4h^1" , "-P 4 2"          , "P 4/m m m"                      , "P 4/m 2/m 2/m"      , "P4/mmm"    , ""     , Bravais::P, 15, 400},
  {124, "D4h^2" , "-P 4 2c"         , "P 4/m c c"                      , "P 4/m 2/c 2/c"      , "P4/mcc"    , ""     , Bravais::P, 15, 401},
  {125, "D4h^3" , "P 4 2 -1ab"      , "P 4/n b m"                      , "P 4/n 2/b 2/m"      , "P4/nbm"    , "1"    , Bravais::P, 15, 402},
  {125, "D4h^3" , "-P 4a 2b"        , "P 4/n b m"                      , "P 4/n 2/b 2/m"      , "P4/nbm"    , "2"    , Bravais::P, 15, 403},
  {126, "D4h^4" , "P 4 2 -1n"       , "P 4/n n c"                      , "P 4/n 2/n 2/c"      , "P4/nnc"    , "1"    , Bravais::P, 15, 404},
  {126, "D4h^4" , "-P 4a 2bc"       , "P 4/n n c"                      , "P 4/n 2/n 2/c"      , "P4/nnc"    , "2"    , Bravais::P, 15, 405},
  {127, "D4h^5" , "-P 4 2ab"        , "P 4/m b m"                      , "P 4/m 2_1/b m"      , "P4/mbm"    , ""     , Bravais::P, 15, 406},
  {128, "D4h^6" , "-P 4 2n"         , "P 4/m n c"                      , "P 4/m 2_1/n c"      , "P4/mnc"    , ""     , Bravais::P, 15, 407},
  {129, "D4h^7" , "P 4ab 2ab -1ab"  , "P 4/n m m"                      , "P 4/n 2_1/m m"      , "P4/nmm"    , "1"    , Bravais::P, 15, 408},
  {129, "D4h^7" , "-P 4a 2a"        , "P 4/n m m"                      , "P 4/n 2_1/m m"      , "P4/nmm"    , "2"    , Bravais::P, 15, 409},
  {130, "D4h^8" , "P 4ab 2n -1ab"   , "P 4/n c c"                      , "P 4/n 2_1/c c"      , "P4/ncc"    , "1"    , Bravais::P, 15, 410},
  {130, "D4h^8" , "-P 4a 2ac"       , "P 4/n c c"                      , "P 4/n 2_1/c c"      , "P4/ncc"    , "2"    , Bravais::P, 15, 411},
  {131, "D4h^9" , "-P 4c 2"         , "P 4_2/m m c"                    , "P 4_2/m 2/m 2/c"    , "P4_2/mmc"  , ""     , Bravais::P, 15, 412},
  {132, "D4h^10", "-P 4c 2c"        , "P 4_2/m c m"                    , "P 4_2/m 2/c 2/m"    , "P4_2/mcm"  , ""     , Bravais::P, 15, 413},
  {133, "D4h^11", "P 4n 2c -1n"     , "P 4_2/n b c"                    , "P 4_2/n 2/b 2/c"    , "P4_2/nbc"  , "1"    , Bravais::P, 15, 414},
  {133, "D4h^11", "-P 4ac 2b"       , "P 4_2/n b c"                    , "P 4_2/n 2/b 2/c"    , "P4_2/nbc"  , "2"    , Bravais::P, 15, 415},
  {134, "D4h^12", "P 4n 2 -1n"      , "P 4_2/n n m"                    , "P 4_2/n 2/n 2/m"    , "P4_2/nnm"  , "1"    , Bravais::P, 15, 416},
  {134, "D4h^12", "-P 4ac 2bc"      , "P 4_2/n n m"                    , "P 4_2/n 2/n 2/m"    , "P4_2/nnm"  , "2"    , Bravais::P, 15, 417},
  {135, "D4h^13", "-P 4c 2ab"       , "P 4_2/m b c"                    , "P 4_2/m 2_1/b 2/c"  , "P4_2/mbc"  , ""     , Bravais::P, 15, 418},
  {136, "D4h^14", "-P 4n 2n"        , "P 4_2/m n m"                    , "P 4_2/m 2_1/n 2/m"  , "P4_2/mnm"  , ""     , Bravais::P, 15, 419},
  {137, "D4h^15", "P 4n 2n -1n"     , "P 4_2/n m c"                    , "P 4_2/n 2_1/m 2/c"  , "P4_2/nmc"  , "1"    , Bravais::P, 15, 420},
  {137, "D4h^15", "-P 4ac 2a"       , "P 4_2/n m c"                    , "P 4_2/n 2_1/m 2/c"  , "P4_2/nmc"  , "2"    , Bravais::P, 15, 421},
  {138, "D4h^16", "P 4n 2ab -1n"    , "P 4_2/n c m"                    , "P 4_2/n 2_1/c 2/m"  , "P4_2/ncm"  , "1"    , Bravais::P, 15, 422},
  {138, "D4h^16", "-P 4ac 2ac"      , "P 4_2/n c m"                    , "P 4_2/n 2_1/c 2/m"  , "P4_2/ncm"  , "2"    , Bravais::P, 15, 423},
  {139, "D4h^17", "-I 4 2"          , "I 4/m m m"                      , "I 4/m 2/m 2/m"      , "I4/mmm"    , ""     , Bravais::I, 15, 424},
  {140, "D4h^18", "-I 4 2c"         , "I 4/m c m"                      , "I 4/m 2/c 2/m"      , "I4/mcm"    , ""     , Bravais::I, 15, 425},
  {141, "D4h^19", "I 4bw 2bw -1bw"  , "I 4_1/a m d"                    , "I 4_1/a 2/m 2/d"    , "I4_1/amd"  , "1"    , Bravais::I, 15, 426},
  {141, "D4h^19", "-I 4bd 2"        , "I 4_1/a m d"                    , "I 4_1/a 2/m 2/d"    , "I4_1/amd"  , "2"    , Bravais::I, 15, 427},
  {142, "D4h^20", "I 4bw 2aw -1bw"  , "I 4_1/a c d"                    , "I 4_1/a 2/c 2/d"    , "I4_1/acd"  , "1"    , Bravais::I, 15, 428},
  {142, "D4h^20", "-I 4bd 2c"       , "I 4_1/a c d"                    , "I 4_1/a 2/c 2/d"    , "I4_1/acd"  , "2"    , Bravais::I, 15, 429},
  {143, "C3^1"  , "P 3"             , "P 3"                            , "P 3"                , "P3"        , ""     , Bravais::P, 16, 430},
  {144, "C3^2"  , "P 31"            , "P 3_1"                          , "P 3_1"              , "P3_1"      , ""     , Bravais::P, 16, 431},
  {145, "C3^3"  , "P 32"            , "P 3_2"                          , "P 3_2"              , "P3_2"      , ""     , Bravais::P, 16, 432},
  {146, "C3^4"  , "R 3"             , "R 3"                            , "R 3"                , "R3"        , "H"    , Bravais::R, 16, 433},
  {146, "C3^4"  , "P 3*"            , "R 3"                            , "R 3"                , "R3"        , "R"    , Bravais::P, 16, 434},
  {147, "C3i^1" , "-P 3"            , "P -3"                           , "P -3"               , "P-3"       , ""     , Bravais::P, 17, 435},
  {148, "C3i^2" , "-R 3"            , "R -3"                           , "R -3"               , "R-3"       , "H"    , Bravais::R, 17, 436},
  {148, "C3i^2" , "-P 3*"           , "R -3"                           , "R -3"               , "R-3"       , "R"    , Bravais::P, 17, 437},
  {149, "D3^1"  , "P 3 2"           , "P 3 1 2"                        , "P 3 1 2"            , "P312"      , ""     , Bravais::P, 18, 438},
  {150, "D3^2"  , "P 3 2\""         , "P 3 2 1"                        , "P 3 2 1"            , "P321"      , ""     , Bravais::P, 18, 439},
  {151, "D3^3"  , "P 31 2c (0 0 1)" , "P 3_1 1 2"                      , "P 3_1 1 2"          , "P3_112"    , ""     , Bravais::P, 18, 440},
  {152, "D3^4"  , "P 31 2\""        , "P 3_1 2 1"                      , "P 3_1 2 1"          , "P3_121"    , ""     , Bravais::P, 18, 441},
  {153, "D3^5"  , "P 32 2c (0 0 -1)", "P 3_2 1 2"                      , "P 3_2 1 2"          , "P3_212"    , ""     , Bravais::P, 18, 442},
  {154, "D3^6"  , "P 32 2\""        , "P 3_2 2 1"                      , "P 3_2 2 1"          , "P3_221"    , ""     , Bravais::P, 18, 443},
  {155, "D3^7"  , "R 3 2\""         , "R 3 2"                          , "R 3 2"              , "R32"       , "H"    , Bravais::R, 18, 444},
  {155, "D3^7"  , "P 3* 2"          , "R 3 2"                          , "R 3 2"              , "R32"       , "R"    , Bravais::P, 18, 445},
  {156, "C3v^1" , "P 3 -2\""        , "P 3 m 1"                        , "P 3 m 1"            , "P3m1"      , ""     , Bravais::P, 19, 446},
  {157, "C3v^2" , "P 3 -2"          , "P 3 1 m"                        , "P 3 1 m"            , "P31m"      , ""     , Bravais::P, 19, 447},
  {158, "C3v^3" , "P 3 -2\"c"       , "P 3 c 1"                        , "P 3 c 1"            , "P3c1"      , ""     , Bravais::P, 19, 448},
  {159, "C3v^4" , "P 3 -2c"         , "P 3 1 c"                        , "P 3 1 c"            , "P31c"      , ""     , Bravais::P, 19, 449},
  {160, "C3v^5" , "R 3 -2\""        , "R 3 m"                          , "R 3 m"              , "R3m"       , "H"    , Bravais::R, 19, 450},
  {160, "C3v^5" , "P 3* -2"         , "R 3 m"                          , "R 3 m"              , "R3m"       , "R"    , Bravais::P, 19, 451},
  {161, "C3v^6" , "R 3 -2\"c"       , "R 3 c"                          , "R 3 c"              , "R3c"       , "H"    , Bravais::R, 19, 452},
  {161, "C3v^6" , "P 3* -2n"        , "R 3 c"                          , "R 3 c"              , "R3c"       , "R"    , Bravais::P, 19, 453},
  {162, "D3d^1" , "-P 3 2"          , "P -3 1 m"                       , "P -3 1 2/m"         , "P-31m"     , ""     , Bravais::P, 20, 454},
  {163, "D3d^2" , "-P 3 2c"         , "P -3 1 c"                       , "P -3 1 2/c"         , "P-31c"     , ""     , Bravais::P, 20, 455},
  {164, "D3d^3" , "-P 3 2\""        , "P -3 m 1"                       , "P -3 2/m 1"         , "P-3m1"     , ""     , Bravais::P, 20, 456},
  {165, "D3d^4" , "-P 3 2\"c"       , "P -3 c 1"                       , "P -3 2/c 1"         , "P-3c1"     , ""     , Bravais::P, 20, 457},
  {166, "D3d^5" , "-R 3 2\""        , "R -3 m"                         , "R -3 2/m"           , "R-3m"      , "H"    , Bravais::R, 20, 458},
  {166, "D3d^5" , "-P 3* 2"         , "R -3 m"                         , "R -3 2/m"           , "R-3m"      , "R"    , Bravais::P, 20, 459},
  {167, "D3d^6" , "-R 3 2\"c"       , "R -3 c"                         , "R -3 2/c"           , "R-3c"      , "H"    , Bravais::R, 20, 460},
  {167, "D3d^6" , "-P 3* 2n"        , "R -3 c"                         , "R -3 2/c"           , "R-3c"      , "R"    , Bravais::P, 20, 461},
  {168, "C6^1"  , "P 6"             , "P 6"                            , "P 6"                , "P6"        , ""     , Bravais::P, 21, 462},
  {169, "C6^2"  , "P 61"            , "P 6_1"                          , "P 6_1"              , "P6_1"      , ""     , Bravais::P, 21, 463},
  {170, "C6^3"  , "P 65"            , "P 6_5"                          , "P 6_5"              , "P6_5"      , ""     , Bravais::P, 21, 464},
  {171, "C6^4"  , "P 62"            , "P 6_2"                          , "P 6_2"              , "P6_2"      , ""     , Bravais::P, 21, 465},
  {172, "C6^5"  , "P 64"            , "P 6_4"                          , "P 6_4"              , "P6_4"      , ""     , Bravais::P, 21, 466},
  {173, "C6^6"  , "P 6c"            , "P 6_3"                          , "P 6_3"              , "P6_3"      , ""     , Bravais::P, 21, 467},
  {174, "C3h^1" , "P -6"            , "P -6"                           , "P -6"               , "P-6"       , ""     , Bravais::P, 22, 468},
  {175, "C6h^1" , "-P 6"            , "P 6/m"                          , "P 6/m"              , "P6/m"      , ""     , Bravais::P, 23, 469},
  {176, "C6h^2" , "-P 6c"           , "P 6_3/m"                        , "P 6_3/m"            , "P6_3/m"    , ""     , Bravais::P, 23, 470},
  {177, "D6^1"  , "P 6 2"           , "P 6 2 2"                        , "P 6 2 2"            , "P622"      , ""     , Bravais::P, 24, 471},
  {178, "D6^2"  , "P 61 2 (0 0 -1)" , "P 6_1 2 2"                      , "P 6_1 2 2"          , "P6_122"    , ""     , Bravais::P, 24, 472},
  {179, "D6^3"  , "P 65 2 (0 0 1)"  , "P 6_5 2 2"                      , "P 6_5 2 2"          , "P6_522"    , ""     , Bravais::P, 24, 473},
  {180, "D6^4"  , "P 62 2c (0 0 1)" , "P 6_2 2 2"                      , "P 6_2 2 2"          , "P6_222"    , ""     , Bravais::P, 24, 474},
  {181, "D6^5"  , "P 64 2c (0 0 -1)", "P 6_4 2 2"                      , "P 6_4 2 2"          , "P6_422"    , ""     , Bravais::P, 24, 475},
  {182, "D6^6"  , "P 6c 2c"         , "P 6_3 2 2"                      , "P 6_3 2 2"          , "P6_322"    , ""     , Bravais::P, 24, 476},
  {183, "C6v^1" , "P 6 -2"          , "P 6 m m"                        , "P 6 m m"            , "P6mm"      , ""     , Bravais::P, 25, 477},
  {184, "C6v^2" , "P 6 -2c"         , "P 6 c c"                        , "P 6 c c"            , "P6cc"      , ""     , Bravais::P, 25, 478},
  {185, "C6v^3" , "P 6c -2"         , "P 6_3 c m"                      , "P 6_3 c m"          , "P6_3cm"    , ""     , Bravais::P, 25, 479},
  {186, "C6v^4" , "P 6c -2c"        , "P 6_3 m c"                      , "P 6_3 m c"          , "P6_3mc"    , ""     , Bravais::P, 25, 480},
  {187, "D3h^1" , "P -6 2"          , "P -6 m 2"                       , "P -6 m 2"           , "P-6m2"     , ""     , Bravais::P, 26, 481},
  {188, "D3h^2" , "P -6c 2"         , "P -6 c 2"                       , "P -6 c 2"           , "P-6c2"     , ""     , Bravais::P, 26, 482},
  {189, "D3h^3" , "P -6 -2"         , "P -6 2 m"                       , "P -6 2 m"           , "P-62m"     , ""     , Bravais::P, 26, 483},
  {190, "D3h^4" , "P -6c -2c"       , "P -6 2 c"                       , "P -6 2 c"           , "P-62c"     , ""     , Bravais::P, 26, 484},
  {191, "D6h^1" , "-P 6 2"          , "P 6/m m m"                      , "P 6/m 2/m 2/m"      , "P6/mmm"    , ""     , Bravais::P, 27, 485},
  {192, "D6h^2" , "-P 6 2c"         , "P 6/m c c"                      , "P 6/m 2/c 2/c"      , "P6/mcc"    , ""     , Bravais::P, 27, 486},
  {193, "D6h^3" , "-P 6c 2"         , "P 6_3/m c m"                    , "P 6_3/m 2/c 2/m"    , "P6_3/mcm"  , ""     , Bravais::P, 27, 487},
  {194, "D6h^4" , "-P 6c 2c"        , "P 6_3/m m c"                    , "P 6_3/m 2/m 2/c"    , "P6_3/mmc"  , ""     , Bravais::P, 27, 488},
  {195, "T^1"   , "P 2 2 3"         , "P 2 3"                          , "P 2 3"              , "P23"       , ""     , Bravais::P, 28, 489},
  {196, "T^2"   , "F 2 2 3"         , "F 2 3"                          , "F 2 3"              , "F23"       , ""     , Bravais::F, 28, 490},
  {197, "T^3"   , "I 2 2 3"         , "I 2 3"                          , "I 2 3"              , "I23"       , ""     , Bravais::I, 28, 491},
  {198, "T^4"   , "P 2ac 2ab 3"     , "P 2_1 3"                        , "P 2_1 3"            , "P2_13"     , ""     , Bravais::P, 28, 492},
  {199, "T^5"   , "I 2b 2c 3"       , "I 2_1 3"                        , "I 2_1 3"            , "I2_13"     , ""     , Bravais::I, 28, 493},
  {200, "Th^1"  , "-P 2 2 3"        , "P m -3"                         , "P 2/m -3"           , "Pm-3"      , ""     , Bravais::P, 29, 494},
  {201, "Th^2"  , "P 2 2 3 -1n"     , "P n -3"                         , "P 2/n -3"           , "Pn-3"      , "1"    , Bravais::P, 29, 495},
  {201, "Th^2"  , "-P 2ab 2bc 3"    , "P n -3"                         , "P 2/n -3"           , "Pn-3"      , "2"    , Bravais::P, 29, 496},
  {202, "Th^3"  , "-F 2 2 3"        , "F m -3"                         , "F 2/m -3"           , "Fm-3"      , ""     , Bravais::F, 29, 497},
  {203, "Th^4"  , "F 2 2 3 -1d"     , "F d -3"                         , "F 2/d -3"           , "Fd-3"      , "1"    , Bravais::F, 29, 498},
  {203, "Th^4"  , "-F 2uv 2vw 3"    , "F d -3"                         , "F 2/d -3"           , "Fd-3"      , "2"    , Bravais::F, 29, 499},
  {204, "Th^5"  , "-I 2 2 3"        , "I m -3"                         , "I 2/m -3"           , "Im-3"      , ""     , Bravais::I, 29, 500},
  {205, "Th^6"  , "-P 2ac 2ab 3"    , "P a -3"                         , "P 2_1/a -3"         , "Pa-3"      , ""     , Bravais::P, 29, 501},
  {206, "Th^7"  , "-I 2b 2c 3"      , "I a -3"                         , "I 2_1/a -3"         , "Ia-3"      , ""     , Bravais::I, 29, 502},
  {207, "O^1"   , "P 4 2 3"         , "P 4 3 2"                        , "P 4 3 2"            , "P432"      , ""     , Bravais::P, 30, 503},
  {208, "O^2"   , "P 4n 2 3"        , "P 4_2 3 2"                      , "P 4_2 3 2"          , "P4_232"    , ""     , Bravais::P, 30, 504},
  {209, "O^3"   , "F 4 2 3"         , "F 4 3 2"                        , "F 4 3 2"            , "F432"      , ""     , Bravais::F, 30, 505},
  {210, "O^4"   , "F 4d 2 3"        , "F 4_1 3 2"                      , "F 4_1 3 2"          , "F4_132"    , ""     , Bravais::F, 30, 506},
  {211, "O^5"   , "I 4 2 3"         , "I 4 3 2"                        , "I 4 3 2"            , "I432"      , ""     , Bravais::I, 30, 507},
  {212, "O^6"   , "P 4acd 2ab 3"    , "P 4_3 3 2"                      , "P 4_3 3 2"          , "P4_332"    , ""     , Bravais::P, 30, 508},
  {213, "O^7"   , "P 4bd 2ab 3"     , "P 4_1 3 2"                      , "P 4_1 3 2"          , "P4_132"    , ""     , Bravais::P, 30, 509},
  {214, "O^8"   , "I 4bd 2c 3"      , "I 4_1 3 2"                      , "I 4_1 3 2"          , "I4_132"    , ""     , Bravais::I, 30, 510},
  {215, "Td^1"  , "P -4 2 3"        , "P -4 3 m"                       , "P -4 3 m"           , "P-43m"     , ""     , Bravais::P, 31, 511},
  {216, "Td^2"  , "F -4 2 3"        , "F -4 3 m"                       , "F -4 3 m"           , "F-43m"     , ""     , Bravais::F, 31, 512},
  {217, "Td^3"  , "I -4 2 3"        , "I -4 3 m"                       , "I -4 3 m"           , "I-43m"     , ""     , Bravais::I, 31, 513},
  {218, "Td^4"  , "P -4n 2 3"       , "P -4 3 n"                       , "P -4 3 n"           , "P-43n"     , ""     , Bravais::P, 31, 514},
  {219, "Td^5"  , "F -4c 2 3"       , "F -4 3 c"                       , "F -4 3 c"           , "F-43c"     , ""     , Bravais::F, 31, 515},
  {220, "Td^6"  , "I -4bd 2c 3"     , "I -4 3 d"                       , "I -4 3 d"           , "I-43d"     , ""     , Bravais::I, 31, 516},
  {221, "Oh^1"  , "-P 4 2 3"        , "P m -3 m"                       , "P 4/m -3 2/m"       , "Pm-3m"     , ""     , Bravais::P, 32, 517},
  {222, "Oh^2"  , "P 4 2 3 -1n"     , "P n -3 n"                       , "P 4/n -3 2/n"       , "Pn-3n"     , "1"    , Bravais::P, 32, 518},
  {222, "Oh^2"  , "-P 4a 2bc 3"     , "P n -3 n"                       , "P 4/n -3 2/n"       , "Pn-3n"     , "2"    , Bravais::P, 32, 519},
  {223, "Oh^3"  , "-P 4n 2 3"       , "P m -3 n"                       , "P 4_2/m -3 2/n"     , "Pm-3n"     , ""     , Bravais::P, 32, 520},
  {224, "Oh^4"  , "P 4n 2 3 -1n"    , "P n -3 m"                       , "P 4_2/n -3 2/m"     , "Pn-3m"     , "1"    , Bravais::P, 32, 521},
  {224, "Oh^4"  , "-P 4bc 2bc 3"    , "P n -3 m"                       , "P 4_2/n -3 2/m"     , "Pn-3m"     , "2"    , Bravais::P, 32, 522},
  {225, "Oh^5"  , "-F 4 2 3"        , "F m -3 m"                       , "F 4/m -3 2/m"       , "Fm-3m"     , ""     , Bravais::F, 32, 523},
  {226, "Oh^6"  , "-F 4c 2 3"       , "F m -3 c"                       , "F 4/m -3 2/c"       , "Fm-3c"     , ""     , Bravais::F, 32, 524},
  {227, "Oh^7"  , "F 4d 2 3 -1d"    , "F d -3 m"                       , "F 4_1/d -3 2/m"     , "Fd-3m"     , "1"    , Bravais::F, 32, 525},
  {227, "Oh^7"  , "-F 4vw 2vw 3"    , "F d -3 m"                       , "F 4_1/d -3 2/m"     , "Fd-3m"     , "2"    , Bravais::F, 32, 526},
  {228, "Oh^8"  , "F 4d 2 3 -1cd"   , "F d -3 c"                       , "F 4_1/d -3 2/c"     , "Fd-3c"     , "1"    , Bravais::F, 32, 527},
  {228, "Oh^8"  , "-F 4cvw 2vw 3"   , "F d -3 c"                       , "F 4_1/d -3 2/c"     , "Fd-3c"     , "2"    , Bravais::F, 32, 528},
  {229, "Oh^9"  , "-I 4 2 3"        , "I m -3 m"                       , "I 4/m -3 2/m"       , "Im-3m"     , ""     , Bravais::I, 32, 529},
  {230, "Oh^10" , "-I 4bd 2c 3"     , "I a -3 d"                       , "I 4_1/a -3 2/d"     , "Ia-3d"     , ""     , Bravais::I, 32, 530},
};

bool brille::hall_number_ok(const int h) {return h>0 && h<531;}

int brille::international_number_to_hall_number(const int n, const std::string& c){
	if (n>=0 && n<230) for (int i=1; i<531; i++){
		Spacegroup spg(ALL_SPACEGROUPS[i]);
    if (n==spg.number && (0==c.size() || 0==c.compare(spg.choice))) return i;
	}
	return -1;
}

int brille::international_string_to_hall_number(const std::string& n, const std::string& c){
  for (int i=1; i<531; i++) {
	  Spacegroup spg(ALL_SPACEGROUPS[i]);
	  // now check for matching international table names
    if( (0==c.size() || 0==c.compare(spg.choice)) &&
        (  0==n.compare(spg.international)
        || 0==n.compare(spg.international_full)
        || 0==n.compare(spg.international_short))
      ) return i;
  }
  return 0; // no matching strings
}

int brille::hall_symbol_to_hall_number(const std::string& hsymbol){
  for (int i=1; i<531; i++){
    Spacegroup spg(ALL_SPACEGROUPS[i]);
    // only accept exact matches
    if (0==hsymbol.compare(spg.hall_symbol)) return i;
  }
  return 0; // no exact matches
}

// When we don't know if we're searching for a Hall Symbol or
// International Table Name/Symbol search for either, starting with Hall symbols
int brille::string_to_hall_number(const std::string& str, const std::string& ch){
  int tmp = brille::hall_symbol_to_hall_number(str);
  if (!tmp) tmp = brille::international_string_to_hall_number(str, ch);
  return tmp;
}

void Spacegroup::set_from_hall_number(const int hn){
  Spacegroup spg(ALL_SPACEGROUPS[brille::hall_number_ok(hn) ? hn : 0]);
  this->number              = spg.number;
  this->schoenflies         = spg.schoenflies;
  this->hall_symbol         = spg.hall_symbol;
  this->international       = spg.international;
  this->international_full  = spg.international_full;
  this->international_short = spg.international_short;
  this->choice              = spg.choice;
  this->bravais             = spg.bravais;
  this->pointgroup_number   = spg.pointgroup_number;
  this->hall_number         = spg.hall_number;
}
Symmetry Spacegroup::get_spacegroup_symmetry() const {
  HallSymbol hs(this->hall_symbol);
  Symmetry generators = hs.get_generators();
  Symmetry fullsym = generators.generate();
  return fullsym;
}
PointSymmetry Spacegroup::get_pointgroup_symmetry(const int time_reversal) const{
  Symmetry sym = this->get_spacegroup_symmetry();
  std::vector<std::array<int,9>> uniqrots = get_unique_rotations(sym.getallr(), time_reversal);
  return PointSymmetry(uniqrots);
}
