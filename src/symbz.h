/* Copyright (C) 2008 Atsushi Togo */
/* All rights reserved. */

/* This file is part of spglib. */

/* Redistribution and use in source and binary forms, with or without */
/* modification, are permitted provided that the following conditions */
/* are met: */

/* * Redistributions of source code must retain the above copyright */
/*   notice, this list of conditions and the following disclaimer. */

/* * Redistributions in binary form must reproduce the above copyright */
/*   notice, this list of conditions and the following disclaimer in */
/*   the documentation and/or other materials provided with the */
/*   distribution. */

/* * Neither the name of the phonopy project nor the names of its */
/*   contributors may be used to endorse or promote products derived */
/*   from this software without specific prior written permission. */

/* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS */
/* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT */
/* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS */
/* FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE */
/* COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, */
/* INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, */
/* BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; */
/* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER */
/* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT */
/* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN */
/* ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE */
/* POSSIBILITY OF SUCH DAMAGE. */

#ifndef __spglib_H__
#define __spglib_H__

#ifdef __cplusplus
extern "C" {
#endif
/* SPGCONST is used instead of 'const' so to avoid gcc warning. */
/* However there should be better way than this way.... */
#ifndef SPGCONST
#define SPGCONST
#endif

#include <stddef.h>

  typedef enum {
    SYMBZ_SUCCESS = 0,
    SYMBZ_SPACEGROUP_SEARCH_FAILED,
    SYMBZ_SYMMETRY_OPERATION_SEARCH_FAILED,
    SYMBZ_POINTGROUP_NOT_FOUND,
    SYMBZ_ARRAY_SIZE_SHORTAGE,
    SYMBZ_NONE,
  } SymBZError;

  typedef struct {
    int spacegroup_number;
    int hall_number;
    char international_symbol[11];
    char hall_symbol[17];
    char choice[6];
    double transformation_matrix[3][3];
    double origin_shift[3];
    int n_operations;
    int (*rotations)[3][3];
    double (*translations)[3];
    int n_atoms;
    int *wyckoffs;
    char (*site_symmetry_symbols)[7];
    int *equivalent_atoms;
    int *mapping_to_primitive;
    int n_std_atoms;
    double std_lattice[3][3];
    int *std_types;
    double (*std_positions)[3];
    double std_rotation_matrix[3][3];
    int *std_mapping_to_primitive;
    /* int pointgroup_number; */
    char pointgroup_symbol[6];
  } SpglibDataset;

  typedef struct {
    int number;
    char international_short[11];
    char international_full[20];
    char international[32];
    char schoenflies[7];
    char hall_symbol[17];
    char choice[6];
    char pointgroup_international[6];
    char pointgroup_schoenflies[4];
    int arithmetic_crystal_class_number;
    char arithmetic_crystal_class_symbol[7];
  } SpglibSpacegroupType;

  int symbz_get_major_version(void);
  int symbz_get_minor_version(void);
  int symbz_get_micro_version(void);

  SymBZError symbz_get_error_code(void);
  const char * symbz_get_error_message(SymBZError spglib_error);

  SpglibSpacegroupType symbz_get_spacegroup_type(const int hall_number);
  size_t symbz_get_symmetry_from_database(int *rotations,double *translations, const size_t size, const int hall_number);

  int symbz_get_pointgroup_rotations_hall_number(int *rotations, const int max_size, const int hall_number, const int is_time_reversal);
  int symbz_get_hall_number_from_international(const char* itname);

  int symbz_get_bz(const double *lengths, const double *angles, const int search_length, const int *max_sizes, double *verts, int *faces, int *fpv, int *counts);

  // size_t symbz_get_bz_step_grid_xyz(const double *lengths, const double *angles, const int search_length, const size_t *multiplicity, const size_t maxN, double *xyz);
  // size_t symbz_get_bz_step_grid_hkl(const double *lengths, const double *angles, const int search_length, const size_t *multiplicity, const size_t maxN, double *xyz);
  // size_t symbz_get_bz_inva_grid_xyz(const double *lengths, const double *angles, const int search_length, const double *stepsize, const int stepisrlu, const size_t maxN, double *xyz);
  // size_t symbz_get_bz_inva_grid_hkl(const double *lengths, const double *angles, const int search_length, const double *stepsize, const int stepisrlu, const size_t maxN, double *xyz);

#ifdef __cplusplus
}
#endif
#endif
