/* Copyright (C) 2010 Atsushi Togo
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

#ifndef __spg_database_H__
#define __spg_database_H__

#include<iostream>
#include "symmetry.h"

typedef struct {
  int number;
  int hall_number;
  int pointgroup_number;
  char schoenflies[7];
  char hall_symbol[17];
  char international[32];
  char international_long[20];
  char international_short[11];
  char choice[6];
  double bravais_lattice[3][3];
  double origin_shift[3];
} Spacegroup;

typedef enum {
  CENTERING_ERROR,
  PRIMITIVE,
  BODY,
  FACE,
  A_FACE,
  B_FACE,
  C_FACE,
  BASE,
  R_CENTER,
} Centering;

typedef struct {
  int number;
  char schoenflies[7];
  char hall_symbol[17];
  char international[32];
  char international_full[20];
  char international_short[11];
  char choice[6];
  Centering centering;
  int pointgroup_number;
} SpacegroupType;

int spgdb_remove_space(char * str, const int nchar);
int spgdb_equal_to_quote(char * str, const int len);

bool hall_number_ok(const int hall_number);
int spgdb_get_operation(int *rot, double *trans, const int idx);
void spgdb_get_operation_index(int indices[2], const int hall_number);
Symmetry spgdb_get_spacegroup_operations(const int hall_number);
SpacegroupType spgdb_get_spacegroup_type(const int hall_number);

int spgdb_international_number_to_hall_number(const int number);
int spgdb_international_to_hall_number(const char* itname);
size_t spgdb_get_symmetry_from_database(int *rotations, double *translations, const size_t space, const int hall_number);

#endif
