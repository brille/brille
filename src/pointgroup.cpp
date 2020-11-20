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

/* This file has evolved from pointgroup.cpp distributed as part of spglib.
   Changes have been made to introduce C++ style classes as well as other
   modifications to suit the needs of brille.
   spglib was licensed under the following BSD 3-clause license:

 Copyright (C) 2008 Atsushi Togo
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

#include <algorithm>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "pointgroup.hpp"
#include "symmetry.hpp"
#include "pointsymmetry.hpp"
#include "utilities.hpp"
#include "spg_database.hpp"

#define NUM_ROT_AXES 73
#define ZERO_PREC 1e-10

using namespace brille;

typedef struct {
  int table[10];
  std::string symbol;
  std::string schoenflies;
  Holohedry holohedry;
  Laue laue;
} PointgroupType;

static PointgroupType pointgroup_data[33] = {
	/*Table 6,̄4,̄3,̄2,̄1,1,2,3,4,6 frequencies*/
  {/* 0*/{0,0,0,0,0,0,0,0,0,0},""     ,""   ,Holohedry::_         ,Laue::_    ,},
  {/* 1*/{0,0,0,0,0,1,0,0,0,0},"1"    ,"C1" ,Holohedry::triclinic ,Laue::_1   ,},
  {/* 2*/{0,0,0,0,1,1,0,0,0,0},"-1"   ,"Ci" ,Holohedry::triclinic ,Laue::_1   ,},
  {/* 3*/{0,0,0,0,0,1,1,0,0,0},"2"    ,"C2" ,Holohedry::monoclinic,Laue::_2m  ,},
  {/* 4*/{0,0,0,1,0,1,0,0,0,0},"m"    ,"Cs" ,Holohedry::monoclinic,Laue::_2m  ,},
  {/* 5*/{0,0,0,1,1,1,1,0,0,0},"2/m"  ,"C2h",Holohedry::monoclinic,Laue::_2m  ,},
  {/* 6*/{0,0,0,0,0,1,3,0,0,0},"222"  ,"D2" ,Holohedry::orthogonal,Laue::_mmm ,},
  {/* 7*/{0,0,0,2,0,1,1,0,0,0},"mm2"  ,"C2v",Holohedry::orthogonal,Laue::_mmm ,},
  {/* 8*/{0,0,0,3,1,1,3,0,0,0},"mmm"  ,"D2h",Holohedry::orthogonal,Laue::_mmm ,},
  {/* 9*/{0,0,0,0,0,1,1,0,2,0},"4"    ,"C4" ,Holohedry::tetragonal,Laue::_4m  ,},
  {/*10*/{0,2,0,0,0,1,1,0,0,0},"-4"   ,"S4" ,Holohedry::tetragonal,Laue::_4m  ,},
  {/*11*/{0,2,0,1,1,1,1,0,2,0},"4/m"  ,"C4h",Holohedry::tetragonal,Laue::_4m  ,},
  {/*12*/{0,0,0,0,0,1,5,0,2,0},"422"  ,"D4" ,Holohedry::tetragonal,Laue::_4mmm,},
  {/*13*/{0,0,0,4,0,1,1,0,2,0},"4mm"  ,"C4v",Holohedry::tetragonal,Laue::_4mmm,},
  {/*14*/{0,2,0,2,0,1,3,0,0,0},"-42m" ,"D2d",Holohedry::tetragonal,Laue::_4mmm,},
  {/*15*/{0,2,0,5,1,1,5,0,2,0},"4/mmm","D4h",Holohedry::tetragonal,Laue::_4mmm,},
  {/*16*/{0,0,0,0,0,1,0,2,0,0},"3"    ,"C3" ,Holohedry::trigonal  ,Laue::_3   ,},
  {/*17*/{0,0,2,0,1,1,0,2,0,0},"-3"   ,"C3i",Holohedry::trigonal  ,Laue::_3   ,},
  {/*18*/{0,0,0,0,0,1,3,2,0,0},"32"   ,"D3 ",Holohedry::trigonal  ,Laue::_3m  ,},
  {/*19*/{0,0,0,3,0,1,0,2,0,0},"3m"   ,"C3v",Holohedry::trigonal  ,Laue::_3m  ,},
  {/*20*/{0,0,2,3,1,1,3,2,0,0},"-3m"  ,"D3d",Holohedry::trigonal  ,Laue::_3m  ,},
  {/*21*/{0,0,0,0,0,1,1,2,0,2},"6"    ,"C6" ,Holohedry::hexagonal ,Laue::_6m  ,},
  {/*22*/{2,0,0,1,0,1,0,2,0,0},"-6"   ,"C3h",Holohedry::hexagonal ,Laue::_6m  ,},
  {/*23*/{2,0,2,1,1,1,1,2,0,2},"6/m"  ,"C6h",Holohedry::hexagonal ,Laue::_6m  ,},
  {/*24*/{0,0,0,0,0,1,7,2,0,2},"622"  ,"D6" ,Holohedry::hexagonal ,Laue::_6mmm,},
  {/*25*/{0,0,0,6,0,1,1,2,0,2},"6mm"  ,"C6v",Holohedry::hexagonal ,Laue::_6mmm,},
  {/*26*/{2,0,0,4,0,1,3,2,0,0},"-6m2" ,"D3h",Holohedry::hexagonal ,Laue::_6mmm,},
  {/*27*/{2,0,2,7,1,1,7,2,0,2},"6/mmm","D6h",Holohedry::hexagonal ,Laue::_6mmm,},
  {/*28*/{0,0,0,0,0,1,3,8,0,0},"23"   ,"T"  ,Holohedry::cubic     ,Laue::_m3  ,},
  {/*29*/{0,0,8,3,1,1,3,8,0,0},"m-3"  ,"Th" ,Holohedry::cubic     ,Laue::_m3  ,},
  {/*30*/{0,0,0,0,0,1,9,8,6,0},"432"  ,"O"  ,Holohedry::cubic     ,Laue::_m3m ,},
  {/*31*/{0,6,0,6,0,1,3,8,0,0},"-43m" ,"Td" ,Holohedry::cubic     ,Laue::_m3m ,},
  {/*32*/{0,6,8,9,1,1,9,8,6,0},"m-3m" ,"Oh" ,Holohedry::cubic     ,Laue::_m3m ,}
};

std::string brille::my_to_string(const Holohedry& h){
  std::string str;
  switch(h){
    case Holohedry::triclinic:  str = "triclinic";  break;
    case Holohedry::monoclinic: str = "monoclinic"; break;
    case Holohedry::orthogonal: str = "orthogonal"; break;
    case Holohedry::tetragonal: str = "tetragonal"; break;
    case Holohedry::trigonal:   str = "trigonal";   break;
    case Holohedry::hexagonal:  str = "hexagonal";  break;
    case Holohedry::cubic:      str = "cubic";      break;
    default: str = "Holohedry Error";
  }
  return str;
}
std::string brille::my_to_string(const Laue& l){
  std::string str;
  switch(l){
    case Laue::_1:    str="1";    break;
    case Laue::_2m:   str="2m";   break;
    case Laue::_mmm:  str="mmm";  break;
    case Laue::_4m:   str="4m";   break;
    case Laue::_4mmm: str="4mmm"; break;
    case Laue::_3:    str="3";    break;
    case Laue::_3m:   str="3m";   break;
    case Laue::_6m:   str="6m";   break;
    case Laue::_6mmm: str="6mmm"; break;
    case Laue::_m3:   str="m3";   break;
    case Laue::_m3m:  str="m3m";  break;
    default: str = "Laue Error";
  }
  return str;
}

static int identity[9] = { 1, 0, 0,  0, 1, 0,  0, 0, 1};
static int inversion[9]= {-1, 0, 0,  0,-1, 0,  0, 0,-1};

/* Each rotation and rotoinversion symmetry opperation, W, has an associated
   vector, ⃗u, which remains unchanged by the symmetry operation, e.g.,
	 		⃗u = W ⃗u
	 Due to the limited number of unique (crystal) pointgroup operations there
	 are also a liminted number of associated vectors. Rather than solving the
	 eigenvalue problem (W-I) ⃗u = 0 when the ⃗u for a given W is required it is
	 faster(?) to check which of the possible rotation axes remains stationary.
*/
static int rot_axes[][3] = {
	/* 0*/ { 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1}, { 0, 1, 1}, { 1, 0, 1}, /* 4*/
	/* 5*/ { 1, 1, 0}, { 0, 1,-1}, {-1, 0, 1}, { 1,-1, 0}, { 1, 1, 1}, /* 9*/
	/*10*/ {-1, 1, 1}, { 1,-1, 1}, { 1, 1,-1}, { 0, 1, 2}, { 2, 0, 1}, /*14*/
	/*15*/ { 1, 2, 0}, { 0, 2, 1}, { 1, 0, 2}, { 2, 1, 0}, { 0,-1, 2}, /*19*/
	/*20*/ { 2, 0,-1}, {-1, 2, 0}, { 0,-2, 1}, { 1, 0,-2}, {-2, 1, 0}, /*24*/
	/*25*/ { 2, 1, 1}, { 1, 2, 1}, { 1, 1, 2}, { 2,-1,-1}, {-1, 2,-1}, /*29*/
	/*30*/ {-1,-1, 2}, { 2, 1,-1}, {-1, 2, 1}, { 1,-1, 2}, { 2,-1, 1}, /*34*/
	/*35*/ { 1, 2,-1}, {-1, 1, 2}, { 3, 1, 2}, { 2, 3, 1}, { 1, 2, 3}, /*39*/
	/*40*/ { 3, 2, 1}, { 1, 3, 2}, { 2, 1, 3}, { 3,-1, 2}, { 2, 3,-1}, /*44*/
	/*45*/ {-1, 2, 3}, { 3,-2, 1}, { 1, 3,-2}, {-2, 1, 3}, { 3,-1,-2}, /*49*/
	/*50*/ {-2, 3,-1}, {-1,-2, 3}, { 3,-2,-1}, {-1, 3,-2}, {-2,-1, 3}, /*54*/
	/*55*/ { 3, 1,-2}, {-2, 3, 1}, { 1,-2, 3}, { 3, 2,-1}, {-1, 3, 2}, /*59*/
	/*60*/ { 2,-1, 3}, { 1, 1, 3}, {-1, 1, 3}, { 1,-1, 3}, {-1,-1, 3}, /*64*/
	/*65*/ { 1, 3, 1}, {-1, 3, 1}, { 1, 3,-1}, {-1, 3,-1}, { 3, 1, 1}, /*69*/
	/*70*/ { 3, 1,-1}, { 3,-1, 1}, { 3,-1,-1},
};

static int get_pointgroup_number_by_rotations(const int *rotations,const int num_rotations);
static int get_pointgroup_number(const PointSymmetry & pointsym);
static int get_pointgroup_class_table(int table[10], const PointSymmetry & pointsym);
static int isometry_type(const int *rot);
// int rotation_order(const int *rot);
static int get_rotation_axis(const int *rot);
static int get_orthogonal_axis(int ortho_axes[], const int *proper_rot, const int rot_order);
static int laue2m(int axes[3], const PointSymmetry & pointsym);
static int laue_one_axis(int axes[3], const PointSymmetry & pointsym, const int rot_order);
static int lauennn(int axes[3], const PointSymmetry & pointsym, const int rot_order);
static int get_axes(int axes[3], const Laue laue, const PointSymmetry & pointsym);
static void get_proper_rotation(int *prop_rot, const int *rot);
static void set_transformation_matrix(int *tmat, const int axes[3]);


/* Retrun pointgroup.number = 0 if failed */
Pointgroup brille::ptg_get_transformation_matrix(int *transform_mat, const int *rotations, const int num_rotations)
{
  int i, pg_num;
  int axes[3];
  PointSymmetry pointsym;
  Pointgroup pointgroup;
  for (i=0; i<9; i++) transform_mat[i]=0;

  pg_num = get_pointgroup_number_by_rotations(rotations, num_rotations);
  pointgroup = Pointgroup(pg_num);
  if (pg_num > 0) {
    pointsym = brille::ptg_get_pointsymmetry(rotations, num_rotations);
    get_axes(axes, pointgroup.laue, pointsym);
    set_transformation_matrix(transform_mat, axes);
  }
  return pointgroup;
}

void Pointgroup::setup(void){
  if (number<0 || number>=33) throw std::runtime_error("Invalid pointgroup number");
  PointgroupType pgt = pointgroup_data[this->number];
  this->symbol = pgt.symbol;
  this->schoenflies = pgt.schoenflies;
  this->holohedry = pgt.holohedry;
  this->laue = pgt.laue;
}

PointSymmetry brille::ptg_get_pointsymmetry(const int *rotations, const int num_rotations)
{
  PointSymmetry pointsym(0);
  // copy the first rotation matrix
	pointsym.add(rotations);
  // so that our outer for loop can start from the second one
  for (int i = 1; i < num_rotations; i++) {
    bool isunique = true;
    for (int j = 0; j < num_rotations; j++)
      if ( brille::utils::equal_matrix(rotations+i*9, pointsym.data(j)) ){
       isunique = false;
       break;
	  }
    if (isunique) pointsym.add(rotations+i*9);
  }
  return pointsym;
}

static int get_pointgroup_number_by_rotations(const int *rotations, const int num_rotations)
{
	PointSymmetry psym = brille::ptg_get_pointsymmetry(rotations, num_rotations);
  return get_pointgroup_number(psym);
}

static int get_pointgroup_number(const PointSymmetry & pointsym)
{
  int table[10];
  PointgroupType pgtype;
  /* Get frequency of each isometry type for the point symmetry operations */
  if (! get_pointgroup_class_table(table, pointsym)) return 0;
	// and check it against the known frequency tables of the 32 crystallograpic pointgroups
  for (int i = 1; i < 33; i++) {
    int eq = 0;
    pgtype = pointgroup_data[i];
    for (int j = 0; j < 10; j++) if (pgtype.table[j] == table[j]) ++eq;
    if (10 == eq) return i; // all frequencies match, so this is pointgroup i.
  }
  return 0;
}

/*! Fill in an isometry class table such that each table[i] contains the number
of pointgroup operations with the i-encoded isometry type.

| index | isometry |
|-------|----------|
|   0   |    -6    |
|   1   |    -4    |
|   2   |    -3    |
|   3   |    -2    |
|   4   |    -1    |
|   5   |     1    |
|   6   |     2    |
|   7   |     3    |
|   8   |     4    |
|   9   |     6    |
*/
static int get_pointgroup_class_table(int table[10], const PointSymmetry & pointsym){
  for (size_t i = 0; i < 10; i++) { table[i] = 0; }
  for (size_t i = 0; i < pointsym.size(); i++) {
    int rot_type = isometry_type(pointsym.data(i));
    if (rot_type == -1) {
      return 0;
    } else {
      table[rot_type]++;
    }
  }
  return 1;
}

/*! \brief Return the encoded Isometry type of a 3×3 crystallographic pointgroup opperation.

	For a crystallographic pointgroup opperation W, find its isometry type
	(1, 2, 3, 4, 6, ̄1, ̄2, ̄3, ̄4, ̄6)
	and endode it in an integer.
	The isometry type depends on the determinant and trace of W and, as given in
	the International Tables section 1.2.2.4

| det(W)   | -1 ||||| +1 |||||
| tr(W)    | -2 | -1 |  0 |  1 | -3 | 3 | -1 | 0 | 1 | 2 |
|----------|----|----|----|----|----|----|----|---|---|---|
| isometry | -6 | -4 | -3 | -2 | -1 | 1 |  2 | 3 | 4 | 6 |
| oder     |  6 |  4 |  3 |  2 |  1 | 1 |  2 | 3 | 4 | 6 |
| encoded  |  0 |  1 |  2 |  3 |  4 | 5 |  6 | 7 | 8 | 9 |
	*/
static int isometry_type(const int *rot){
  int rot_type=-1; // error value
	int tr = brille::utils::trace(rot);
  if ( brille::utils::matrix_determinant(rot) == -1) {
    switch (tr) {
    case -2: /* -6 */ rot_type = 0; break;
    case -1: /* -4 */ rot_type = 1; break;
    case  0: /* -3 */ rot_type = 2; break;
    case  1: /* -2 */ rot_type = 3; break;
    case -3: /* -1 */ rot_type = 4; break;
    }
  } else {
    switch (tr) {
    case  3: /*  1 */ rot_type = 5; break;
    case -1: /*  2 */ rot_type = 6; break;
    case  0: /*  3 */ rot_type = 7; break;
    case  1: /*  4 */ rot_type = 8; break;
    case  2: /*  6 */ rot_type = 9; break;
    }
  }
  return rot_type;
}
int brille::isometry_value(const int *rot){
	int values[] = {-6,-4,-3,-2,-1,1,2,3,4,6};
	int isometry = isometry_type(rot);
	return (isometry>=0) ? values[isometry] : 0;
}
int brille::rotation_order(const int *rot){
	int orders[] = {6,4,3,2,1,1,2,3,4,6};
	int isometry = isometry_type(rot);
	return (isometry>=0) ? orders[isometry] : 0;
}

static int get_axes(int axes[3], const Laue laue, const PointSymmetry & pointsym)
{
  switch (laue) {
  case Laue::_1   : axes[0] = 0; axes[1] = 1; axes[2] = 2; break;
  case Laue::_2m  : laue2m(axes, pointsym); break;
  case Laue::_mmm : lauennn(axes, pointsym, 2); break;
  case Laue::_4m  : laue_one_axis(axes, pointsym, 4); break;
  case Laue::_4mmm: laue_one_axis(axes, pointsym, 4); break;
  case Laue::_3   : laue_one_axis(axes, pointsym, 3); break;
  case Laue::_3m  : laue_one_axis(axes, pointsym, 3); break;
  case Laue::_6m  : laue_one_axis(axes, pointsym, 3); break;
  case Laue::_6mmm: laue_one_axis(axes, pointsym, 3); break;
  case Laue::_m3  : lauennn(axes, pointsym, 2); break;
  case Laue::_m3m : lauennn(axes, pointsym, 4); break;
	default: return 0;
  }
  return 1;
}


static int laue2m(int axes[3], const PointSymmetry & pointsym)
{
  int num_ortho_axis, min_norm;
  int prop_rot[9], t_mat[9];
  int ortho_axes[NUM_ROT_AXES];

  for (size_t i = 0; i < pointsym.size(); i++) {
    get_proper_rotation(prop_rot, pointsym.data(i));
    // Find the first two-fold rotation
		// det(prop_rot)==1, so two-fold has tr==-1 for 2 and ̄2.
		if (brille::utils::trace(prop_rot)==-1){
			// It's rotation-invariant axis is axes[1]
	    axes[1] = get_rotation_axis(prop_rot);
	    break;
		}
  }
  // find all unique crystallographic point group rotation axes
	// which are perpendicular to the proper rotation (with rotation axis axes[1])
  num_ortho_axis = get_orthogonal_axis(ortho_axes, prop_rot, 2);
	// if there are no perpendicular axes something is wrong.
  if (!num_ortho_axis) return 0;

	int * norm_squared = nullptr;
  norm_squared = new int[num_ortho_axis]();
	const int max_norm=8;
	int below_count = 0;
	for (int i=0; i<num_ortho_axis; ++i){
		norm_squared[i] = brille::utils::vector_norm_squared(rot_axes[ortho_axes[i]]);
		if (norm_squared[i] < max_norm - ZERO_PREC) ++below_count;
	}
	if (below_count < 2) return 0;

	// check all orthogonal axes for the shortest one with squared length < 8
	/* the vectors [v: rot_axes] have |v|² ϵ {1, 2, 3, 5, 6, 11, 14} so we are
	   automatically rulling out 36 of the 73 vectors */
  // The shortest with squared norm less than 8 is the second axis.
  min_norm = max_norm;
	int min_i = -1;
  for (int i = 0; i < num_ortho_axis; i++) {
    if (norm_squared[i] < min_norm - ZERO_PREC) {
      min_norm = norm_squared[i];
      axes[0] = ortho_axes[i];
			min_i = i;
    }
  }
	norm_squared[min_i] = max_norm+1;
  /* The third axis */
	// check again all orthogonal axes for the second shortest
  min_norm = max_norm;
	for (int i = 0; i < num_ortho_axis; i++) {
    if (norm_squared[i] < min_norm - ZERO_PREC) {
      min_norm = norm_squared[i];
      axes[2] = ortho_axes[i];
    }
  }
	delete[] norm_squared;

  set_transformation_matrix(t_mat, axes);
  if (brille::utils::matrix_determinant(t_mat) < 0) {
    int tmpval = axes[0];
    axes[0] = axes[2];
    axes[2] = tmpval;
  }

  return 1;
}

static int laue_one_axis(int axes[3], const PointSymmetry & pointsym, const int rot_order)
{
  int num_ortho_axis, is_found=0;
  int axis_vec[3], prop_rot[9], t_mat[9];
  int ortho_axes[NUM_ROT_AXES];

  for (size_t i = 0; i < pointsym.size(); i++) {
    get_proper_rotation(prop_rot, pointsym.data(i));
    /* For order=4, look for a four-fold rotation, which has tr(W)==1 */
    if (rot_order == 4 && brille::utils::trace(prop_rot)==1) {
      axes[2] = get_rotation_axis(prop_rot);
      break;
    }
    /* For order=3, look for a three-fold rotation, which has tr(W)=0 */
    if (rot_order == 3 && brille::utils::trace(prop_rot) == 0) {
      axes[2] = get_rotation_axis(prop_rot);
      break;
    }
  }
	// find all unique crystallographic point group rotation axes
	// which are perpendicular to the proper rotation (with rotation axis axes[1])
	num_ortho_axis = get_orthogonal_axis(ortho_axes, prop_rot, rot_order);
	// if there are no perpendicular axes something is wrong.
	if (!num_ortho_axis) return 0;

	// look for a pair of vectors perpendicular to the rotation axis which have
	// axes[1] = R * axes[0]
  axes[1] = -1;
	bool axes_found=false;
  for (int i = 0; i < num_ortho_axis; i++) {
    axes[0] = ortho_axes[i];
    brille::utils::multiply_matrix_vector(axis_vec, prop_rot, rot_axes[axes[0]]);
    for (int j = 0; j < num_ortho_axis; j++) {
			if (brille::utils::equal_vector<int,3>(axis_vec, rot_axes[ortho_axes[j]]))
				is_found = 1;
			if (!is_found){
				for (int k=0; k<3; ++k)
					axis_vec[k] *= -1;
				if (brille::utils::equal_vector<int,3>(axis_vec, rot_axes[ortho_axes[j]]))
					is_found = -1;
			}
			if (is_found){
				axes[1] = ortho_axes[j] + (is_found < 0 ? NUM_ROT_AXES : 0);
				axes_found = true;
				break;
			}
    }
		if (axes_found){
	    set_transformation_matrix(t_mat, axes);
			/* to avoid F-center choice which has det=4 */
			if (abs(brille::utils::matrix_determinant(t_mat)) > 3) axes_found = false;
		}
		if (axes_found) break;
  }
	/* axes are not correctly found. */
	if (!axes_found) return 0;
  set_transformation_matrix(t_mat, axes);
  if (brille::utils::matrix_determinant(t_mat) < 0) {
    int tmpval = axes[0];
    axes[0] = axes[1];
    axes[1] = tmpval;
  }
  return 1;
}

static int lauennn(int axes[3], const PointSymmetry & pointsym, const int rot_order)
{
  int count=0, axis;
  int prop_rot[9];
  for (int i = 0; i < 3; i++) axes[i] = -1;
	// look for three two-fold or four-fold axes
  for (size_t i = 0; i < pointsym.size(); i++) {
    get_proper_rotation(prop_rot, pointsym.data(i));
		int tr = brille::utils::trace(prop_rot);
		// two-fold rotations have tr(W)==-1, four-fold have tr(W)==1
		if ((tr == -1 && rot_order == 2) || (tr==1 && rot_order==4)){
			axis = get_rotation_axis(prop_rot);
			if (axis!=axes[0] && axis!=axes[1] && axis!=axes[2]){
				if (count >= 3) return 0; // something has gone wrong if there are more than three unique n-fold axes
				axes[count++] = axis;
			}
		}
  }
	// sort the axes:
	if (axes[1] > axes[2] + ZERO_PREC){
		axis = axes[1]; axes[1]=axes[2]; axes[2] = axis;
	}
	if (axes[0] > axes[1] + ZERO_PREC){
		axis = axes[0]; axes[0]=axes[1]; axes[1] = axis;
	}
	if (axes[1] > axes[2] + ZERO_PREC){
		axis = axes[1]; axes[1]=axes[2]; axes[2] = axis;
	}
	set_transformation_matrix(prop_rot, axes);
	if (brille::utils::matrix_determinant(prop_rot) < 0){
		axis = axes[1]; axes[1]=axes[2]; axes[2] = axis;
	}
  return 1;
}

static int get_rotation_axis(const int *proper_rot)
{
  int axis = -1;
  // 1 and ̄1 have no associated rotation axis, so skip them.
  if (!brille::utils::equal_matrix(proper_rot, identity)){
		int vec[3];
	  // The rotation axis for any other isometry is the one which solves the
		// eigenvalue problem (R-I) ⃗u = 0. Since the unique rotation vectors of
		// crystallographic pointgroup operations is limited, it's easier to just
		// check which of the axes is stationary under R ( ⃗u = R ⃗u )
	  for (int i = 0; i < NUM_ROT_AXES; i++) {
	    brille::utils::multiply_matrix_vector(vec, proper_rot, rot_axes[i]);
			if (brille::utils::equal_vector<int,3>(vec, rot_axes[i])){
				axis = i;
				break;
			}
	  }
	}
	return axis;
}


static int orthogonal_to_axis(int ortho_axes[], const int axis){
	int num_ortho=0;
	for (int i=0; i<NUM_ROT_AXES; ++i){
		int i_dot_a = 0;
		for (int j=0; j<3; ++j) i_dot_a += rot_axes[i][j]*rot_axes[axis][j];
		if (0==i_dot_a) ortho_axes[num_ortho++] = i;
	}
	return num_ortho;
}
/*! \brief Determine which of the unique rotation axes are orthogonal to a given
proper rotation of specified order.

For the rotation matrix R, first calculate the matrix S = ∑ᵢ₌₁ⁿ Rⁱ (Rⁿ≡I) then
find the set of all possible rotation axes of crystallographic pointgroups which
are perpendicular to R since they have S ⃗v = 0 -- that is, ⃗v + R ⃗v + … + Rⁿ⁻¹ ⃗v
will be zero if ⃗v ⋅ ⃗u = 0, where ⃗u is the rotation axis of R.
*/
static int get_orthogonal_axis(int ortho_axes[], const int *proper_rot, const int rot_order)
{
  int i, num_ortho_axis=0;
  int vec[3];
  int sum_rot[9]={0,0,0, 0,0,0, 0,0,0}, rot[9]={1,0,0, 0,1,0, 0,0,1};
  for (i = 0; i < rot_order; i++) {
    brille::utils::multiply_matrix_matrix(rot, proper_rot, rot);
    brille::utils::add_matrix(sum_rot, rot, sum_rot);
  }
  for (i = 0; i < NUM_ROT_AXES; i++) {
    brille::utils::multiply_matrix_vector(vec, sum_rot, rot_axes[i]);
    if (vec[0] == 0 && vec[1] == 0 && vec[2] == 0) ortho_axes[num_ortho_axis++] = i;
  }
  return num_ortho_axis;
}

static void get_proper_rotation(int *prop_rot, const int *rot)
{
  if (brille::utils::matrix_determinant(rot) == -1) {
    brille::utils::multiply_matrix_matrix(prop_rot, inversion, rot);
  } else {
    brille::utils::copy_matrix(prop_rot, rot);
  }
}

static void set_transformation_matrix(int *tmat, const int *axes){
	// axes[i] > NUM_ROT_AXES indicates a rotoinversion (e.g, an improper rotation)
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      tmat[i+3*j] = (axes[j]<NUM_ROT_AXES ? 1 : -1) * rot_axes[axes[j]%NUM_ROT_AXES][i];
}


std::vector<std::array<int,9>> brille::get_unique_rotations(const std::vector<std::array<int,9>>& rotations, const int is_time_reversal)
{
	const size_t N = rotations.size();
	std::vector<std::array<int,9>> rot_tmp;
	rot_tmp.reserve(2*N);
	std::vector<int> unique_rot;
	unique_rot.reserve(2*N);

	// copy the rotations to rot_tmp
  std::copy(rotations.begin(), rotations.end(), std::back_inserter(rot_tmp));
	// so that we can add their inverses if time reversal symmetry is allowed
	if (is_time_reversal){
		rot_tmp.resize(2*N);
		for (size_t i=0; i<N; ++i)
		 brille::utils::multiply_matrix_matrix(rot_tmp[i+N].data(), inversion, rot_tmp[i].data());
	}
	// check for uniqueness of the rotations
  bool i_is_not_unique = false;
	unique_rot.push_back(0); // the first rotation is the first unique rotation
  for (size_t i=1; i<rot_tmp.size(); ++i) {
		// check all rotations against the already established-to-be-unique rotations
		for (int j: unique_rot){
			i_is_not_unique = brille::utils::equal_matrix(rot_tmp[j].data(), rot_tmp[i].data());
			if (i_is_not_unique) break;
		}
		// if i is unique, add it to the vector of unique rotation indices
		if (!i_is_not_unique) unique_rot.push_back(static_cast<int>(i));
  }

	// copy from our temporary vector just the unique arrays for output
	std::vector<std::array<int,9>> rot_unique;
	rot_unique.reserve(unique_rot.size());
	for (int i: unique_rot) rot_unique.push_back(rot_tmp[i]);
  // std::transform(unique_rot.begin(), unique_rot.end(), std::back_inserter(rot_unique), [](int i){return rot_tmp[i];});

  return rot_unique;
}

static int _internal_pointgroup_rotations(int *rotations, const int max_size, const std::vector<std::array<int,9>>& all_rots, const int is_time_reversal)
{
  std::vector<std::array<int,9>> unq_rots = brille::get_unique_rotations(all_rots,is_time_reversal);
	if (unq_rots.size() > static_cast<size_t>(max_size)){
    std::string msg = "Indicated maximum size " + std::to_string(max_size);
    msg += " is less than number of unique rotations (";
    msg += std::to_string(unq_rots.size()) + ")";
    throw std::out_of_range(msg);
  }	else {
		for (size_t i=0; i < unq_rots.size(); i++)
      std::copy(unq_rots[i].begin(), unq_rots[i].end(), rotations+i*9);
  }
	return static_cast<int>(unq_rots.size());
}

int brille::get_pointgroup_rotations_hall_number(int *rotations, const int max_size, const int hall_number, const int is_time_reversal){
	Symmetry sym = Spacegroup(hall_number).get_spacegroup_symmetry();
	int num_unique = _internal_pointgroup_rotations(rotations, max_size, sym.getallr(), is_time_reversal);
	return num_unique;
}

std::array<std::array<int,3>,3> brille::rotation_axis_and_perpendicular_vectors(const int* rot){
	int prop_rot[9];
	std::array<std::array<int,3>,3> out {{ {{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}} }};

	get_proper_rotation(prop_rot, rot);
	int axis = get_rotation_axis(prop_rot);
	if (axis < 0) return out;

	out[0] = {rot_axes[axis][0], rot_axes[axis][1], rot_axes[axis][2]};
	int ortho[NUM_ROT_AXES], norm[NUM_ROT_AXES];
	int northo = orthogonal_to_axis(ortho, axis);
	for (int i=0; i<northo; ++i) norm[i] = brille::utils::vector_norm_squared(rot_axes[ortho[i]]);

	int idx=0;
	for (int j=1; j<3; ++j){
		for (int i=0; i<northo; ++i) if (norm[i] < norm[idx]) idx = i;
		for (int k=0; k<3; ++k) out[j][k] = rot_axes[ortho[idx]][k];
		norm[idx] = 100; // none of rot_axes have squared norm above 14.
	}
	return out;
}
