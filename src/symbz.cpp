#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>

#include "symbz.h"
#include "arithmetic.h"
#include "linear_algebra.h"
#include "pointgroup.h"
#include "spg_database.h"
#include "symmetry.h"
#include "bz.h"
#include "grid.h"
#include "bz_grid.h"


/*-------*/
/* error */
/*-------*/
static SymBZError symbz_error_code = SYMBZ_SUCCESS;

typedef struct {
  SymBZError error;
  const char *message;
} SymBZErrorMessage;

static SymBZErrorMessage spglib_error_message[] = {
  {SYMBZ_SUCCESS, "no error"},
  {SYMBZ_SPACEGROUP_SEARCH_FAILED, "spacegroup search failed"},
  {SYMBZ_SYMMETRY_OPERATION_SEARCH_FAILED, "symmetry operation search failed"},
  {SYMBZ_POINTGROUP_NOT_FOUND, "pointgroup not found"},
  {SYMBZ_ARRAY_SIZE_SHORTAGE, "array size shortage"},
  {SYMBZ_NONE, ""},
};


/*-------*/
/* error */
/*-------*/
SymBZError symbz_get_error_code(void){
  return symbz_error_code;
}
const char * symbz_get_error_message(SymBZError error){
  for (int i = 0; i < 100; i++) {
    if (SYMBZ_NONE == spglib_error_message[i].error) break;
    if (error == spglib_error_message[i].error) return spglib_error_message[i].message;
  }
  return nullptr;
}


/*---------*/
/* general */
/*---------*/
/* Return 0 if failed */
size_t symbz_get_symmetry_from_database(int *rotations, double *translations, const size_t space, const int hall_number){
  Symmetry sym;
  sym = spgdb_get_spacegroup_operations(hall_number);

  if ( sym.size()<= space ){
    for (size_t i=0; i<sym.size(); i++){
      sym.getrot(i, rotations + i*9);
      sym.gettran(i, translations + i*3);
    }
	symbz_error_code = SYMBZ_SUCCESS;
    return sym.size();
  }
  if (sym.size()>space)
    symbz_error_code = SYMBZ_ARRAY_SIZE_SHORTAGE;
  else
    symbz_error_code = SYMBZ_SPACEGROUP_SEARCH_FAILED;
  return 0;
}

/* Return spglibtype.number = 0 if failed */
SpglibSpacegroupType symbz_get_spacegroup_type(const int hall_number){
  SpglibSpacegroupType spglibtype;
  SpacegroupType spgtype;
  Pointgroup pointgroup;
  int arth_number;
  char arth_symbol[7];

  spglibtype.number = 0;
  strcpy(spglibtype.schoenflies, "");
  strcpy(spglibtype.hall_symbol, "");
  strcpy(spglibtype.choice, "");
  strcpy(spglibtype.international, "");
  strcpy(spglibtype.international_full, "");
  strcpy(spglibtype.international_short, "");
  strcpy(spglibtype.pointgroup_international, "");
  strcpy(spglibtype.pointgroup_schoenflies, "");
  spglibtype.arithmetic_crystal_class_number = 0;
  strcpy(spglibtype.arithmetic_crystal_class_symbol, "");

  if (0 < hall_number && hall_number < 531) {
    spgtype = spgdb_get_spacegroup_type(hall_number);
    spglibtype.number = spgtype.number;
    strcpy(spglibtype.schoenflies, spgtype.schoenflies);
    strcpy(spglibtype.hall_symbol, spgtype.hall_symbol);
    strcpy(spglibtype.choice, spgtype.choice);
    strcpy(spglibtype.international, spgtype.international);
    strcpy(spglibtype.international_full, spgtype.international_full);
    strcpy(spglibtype.international_short, spgtype.international_short);
    pointgroup = ptg_get_pointgroup(spgtype.pointgroup_number);
    strcpy(spglibtype.pointgroup_international, pointgroup.symbol);
    strcpy(spglibtype.pointgroup_schoenflies, pointgroup.schoenflies);
    arth_number = arth_get_symbol(arth_symbol, spgtype.number);
    spglibtype.arithmetic_crystal_class_number = arth_number;
    strcpy(spglibtype.arithmetic_crystal_class_symbol, arth_symbol);
    symbz_error_code = SYMBZ_SUCCESS;
  } else {
    symbz_error_code = SYMBZ_SPACEGROUP_SEARCH_FAILED;
  }

  return spglibtype;
}


int symbz_get_pointgroup_rotations_hall_number(int *rotations, const int max_size, const int hall_number, const int is_time_reversal){
  return get_pointgroup_rotations_hall_number(rotations,max_size,hall_number,is_time_reversal);
}

int symbz_get_hall_number_from_international(const char* itname){
	int hall_number;
	hall_number = spgdb_international_to_hall_number(itname);
	if (hall_number)
		symbz_error_code = SYMBZ_SUCCESS;
	else
		symbz_error_code = SYMBZ_SPACEGROUP_SEARCH_FAILED;
	return hall_number;
}


// Entirely new-to-SymBZ functions follow:

int symbz_get_bz(const double *lengths, const double *angles, const int search_length, const int *max_sizes, double *verts, int *faces, int *fpv, int *counts){
	Direct dlat(lengths,angles);
	BrillouinZone bz(dlat.star(),search_length);
	counts[0] = bz.vertices_count();
	counts[1] = bz.faces_count();

	if (max_sizes[0]<counts[0] || max_sizes[1]<counts[1]) return 1;

	bz.get_vertices_bare(max_sizes[0],verts);
	bz.get_faces_bare(max_sizes[1],faces);
	bz.get_faces_per_vertex_bare(max_sizes[0],fpv);
	return 0;
}

// size_t symbz_get_bz_step_grid_xyz(const double *lengths, const double *angles, const int search_length, const size_t *multiplicity, const size_t maxN, double *xyz){
// 	Direct d(lengths,angles);
// 	Reciprocal r = d.star();
// 	BrillouinZone bz(r,search_length);
// 	BrillouinZoneGrid3 bzg(r,multiplicity);
// 	size_t num = bzg.get_grid_xyz(maxN, xyz);
// 	return num;
// }
// size_t symbz_get_bz_step_grid_hkl(const double *lengths, const double *angles, const int search_length, const size_t *multiplicity, const size_t maxN, double *xyz){
// 	Direct d(lengths,angles);
// 	d.print();
// 	Reciprocal r = d.star();
// 	r.print();
// 	BrillouinZone bz(r,search_length);
// 	bz.print();
// 	printf("Generate BrillouinZoneGrid3 with multiplicity=[%lu,%lu,%lu] ...",multiplicity[0],multiplicity[1],multiplicity[2]);
// 	BrillouinZoneGrid3 bzg(r,multiplicity);
// 	printf(" done.\n Copy up to %lu (h,k,l) vectors which form the grid ...",maxN);
// 	size_t num = bzg.get_grid_hkl(maxN, xyz);
// 	printf(" done.\n Got back %lu grid points.\n",num);
// 	return num;
// }
// size_t symbz_get_bz_inva_grid_xyz(const double *lengths, const double *angles, const int search_length, const double *stepsize, const int stepisrlu, const size_t maxN, double *xyz){
// 	Direct d(lengths,angles);
// 	Reciprocal r = d.star();
// 	BrillouinZone bz(r,search_length);
// 	BrillouinZoneGrid3 bzg(r,stepsize,stepisrlu);
// 	size_t num = bzg.get_grid_xyz(maxN, xyz);
// 	return num;
// }
// size_t symbz_get_bz_inva_grid_hkl(const double *lengths, const double *angles, const int search_length, const double *stepsize, const int stepisrlu, const size_t maxN, double *xyz){
// 	Direct d(lengths,angles);
// 	Reciprocal r = d.star();
// 	BrillouinZone bz(r,search_length);
// 	BrillouinZoneGrid3 bzg(r,stepsize,stepisrlu);
// 	size_t num = bzg.get_grid_hkl(maxN, xyz);
// 	return num;
// }
