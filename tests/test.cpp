#include "symbz.h"
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <random>

#include "arrayvector.h"
#include "lattice.h"
#include "latvec.h"
#include "bz.h"
#include "bz_grid.h"


static int test_symbz_get_spacegroup_type(void);
static int test_symbz_get_symmetry_from_database(void);
static void show_spacegroup_type(const SpglibSpacegroupType spgtype);

static int test_arrayvector();
static int test_bz_info();
static int test_lattice();
static int test_latvec();
static int test_bz_grid();
static int test_bz_moveinto();

int main(void)
{

  int (*funcs[])(void) = {
		test_arrayvector,
		test_lattice,
		test_latvec,
		test_bz_info,
		test_bz_grid,
		test_bz_moveinto,
    // test_symbz_get_spacegroup_type,
    // test_symbz_get_symmetry_from_database,
    NULL};

  int i, result;

  for (i = 0; i < 30; i++) {
    if (*funcs[i] == NULL) {
      break;
    } else {
      result = (funcs[i])();
      fflush(stdout);
    }
    if (result) {
      return 1;
    }
  }

  return 0;
}


static int test_symbz_get_spacegroup_type(void)
{
  SpglibSpacegroupType spgtype;
  spgtype = symbz_get_spacegroup_type(446);

  printf("*** symbz_get_spacegroup_type ***:\n");
  if (spgtype.number) {
    show_spacegroup_type(spgtype);
    return 0;
  } else {
    return 1;
  }
}

static int test_symbz_get_symmetry_from_database(void)
{
  int *rotations;
  double *translations;
  size_t size=192;
  rotations = new int[size*9]();
  translations = new double[size*3]();

  size = symbz_get_symmetry_from_database(rotations, translations, size, 460);

  if (size) {
    printf("*** symbz_get_symmetry_from_database ***:\n");
    for (size_t i = 0; i < size; i++) {
      printf("--- %d ---\n", i + 1);
      for (size_t j = 0; j < 3; j++){
        for (size_t k = 0; k<3; k++)
          printf("%2d ",rotations[i*9+j+k*3]);
        printf("\n");
	  }
	  for (size_t k=0; k<3; k++)
        printf("%6.4f ",translations[i*3+k]);
      printf("\n");
    }
  }
  delete[] rotations;
  delete[] translations;
  return (size ? 0 : 1);
}

static void show_spacegroup_type(const SpglibSpacegroupType spgtype)
{
  printf("Number:            %d\n", spgtype.number);
  printf("International:     %s\n", spgtype.international_short);
  printf("International:     %s\n", spgtype.international_full);
  printf("International:     %s\n", spgtype.international);
  printf("Schoenflies:       %s\n", spgtype.schoenflies);
  printf("Hall symbol:       %s\n", spgtype.hall_symbol);
  printf("Point group intl:  %s\n", spgtype.pointgroup_international);
  printf("Point group Schoe: %s\n", spgtype.pointgroup_schoenflies);
  printf("Arithmetic cc num. %d\n", spgtype.arithmetic_crystal_class_number);
  printf("Arithmetic cc sym. %s\n", spgtype.arithmetic_crystal_class_symbol);
}


int test_arrayvector(){
	printf("*** test_arrayvector ***\n");
	double tmp[9] = {1,2,3,4,5,6,7,8,9};
	ArrayVector<double> novalues(3,3);
	// novalues.print();
	ArrayVector<double> withvalue(3,3,tmp);
	// withvalue.print();

	novalues += 1.11111111;
	// novalues.print();

	novalues += withvalue;
	// novalues.print();

	ArrayVector<double> newav = novalues/withvalue;
	// newav.print();

	ArrayVector<double> setafterinit;
	setafterinit = newav*3.14159;
	// setafterinit.print();
	// printf("asignment works!\n");

	ArrayVector<double> newerav = newav + 1e-16;

	// if ( newav.isapprox( newerav ) )
		// printf("isapprox works!\n");

	ArrayVector<double> xtr = newav.extract(1u);
	// xtr.print();
	ArrayVector<double> nxt = xtr - newav;
	// nxt.print();

	return 0;
}


static int test_lattice(){
	printf("\n*** test_lattice ***\n");

  Direct d1(2*PI,2*PI,2*PI,PI/2.0,PI/2.0,PI/2.0);
  Reciprocal r1 = d1.star();
  Reciprocal r2(1.,1.,1.,PI/2.0,PI/2.0,PI/2.0);
  r1.print();
  r2.print();

    double eps = std::numeric_limits<double>::epsilon();
	if (! r1.issame(r2)){
		printf("delta (%g %g %g) (%g %g %g)\n",
		r1.get_a()-r2.get_a(),
		r1.get_b()-r2.get_b(),
		r1.get_c()-r2.get_c(),
		r1.get_alpha()-r2.get_alpha(),
		r1.get_beta()-r2.get_beta(),
		r1.get_gamma()-r2.get_gamma());
		printf("sum*eps (%g %g %g) (%g %g %g)\n",
		eps*(r1.get_a()+r2.get_a()),
		eps*(r1.get_b()+r2.get_b()),
		eps*(r1.get_c()+r2.get_c()),
		eps*(r1.get_alpha()+r2.get_alpha()),
		eps*(r1.get_beta()+r2.get_beta()),
		eps*(r1.get_gamma()+r2.get_gamma()));

        return 1;
    }

  return 0;
}

static int test_latvec(){
    printf("\n*** test_latvec ***\n");
    Direct dlat(1.,1.,1.);
    Reciprocal rlat = dlat.star();

    // printf("direct lattice: "); dlat.print();
    // printf("reciprocal lattice: "); rlat.print();

    double tmp[] = { 1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0 };

    LQVec<double> q(rlat, 3, tmp);
    // printf("reciprocal lattice vectors:\n");
    // q.print();

    // printf("with lengths (in inverse angtrom):\n");
		// norm(q).print();

    // printf("cross products:\n");
    // q.get(0).print();
    // printf("cross\n");
    // q.get(1).print();
    // printf("=\n");
    // q.cross(0,1).print();

		// printf("accessing by index\n");
		LQVec<double> b = q[0];
		// b.print();

    return 0;
}
int test_bz_info(){
	printf("*** test_bz_info ***\n");

	// double len1[3] = {3.95,3.95,12.9};
	double len1[3] = {2.0*PI,2.0*PI,2.0*PI};
	double ang1[3] = {PI/2.0, PI/2.0, PI/2.0};
	Direct dlat(len1,ang1);
	Reciprocal rlat = dlat.star();

	int max_search_index = 1;
	BrillouinZone bz(rlat,max_search_index);
	size_t fc = bz.faces_count();
	size_t vc = bz.vertices_count();

	double *verts = safealloc<double>(3*vc);
	int *faces = safealloc<int>(3*fc);
	int *faces_per_vertex = safealloc<int>(3*vc);

	bz.get_faces_bare(3*fc,faces);
	bz.get_vertices_bare(3*vc,verts);
	bz.get_faces_per_vertex_bare(3*vc,faces_per_vertex);


	// printf("Brillouin Zone with %u vertices, %u faces\n",vc,fc);
	// if (vc){
	// 	printf("Vertices at:\n");
	// 	for(size_t i=0; i<vc; i++)
	// 		// printf("[% 7.4f % 7.4f % 7.4f]\n",verts[i*3],verts[i*3+1],verts[i*3+2]);
	// 		printf(" % 7.4f % 7.4f % 7.4f\n",verts[i*3],verts[i*3+1],verts[i*3+2]);
	// }
	// if (fc){
	// 	printf("Brillouin Zone constucted from:\n");
	// 	for(size_t i=0; i<fc; i++)
	// 		// printf("(% 2d % 2d % 2d)\n",faces[i*3],faces[i*3+1],faces[i*3+2]);
	// 		printf(" % 2d % 2d % 2d\n",faces[i*3],faces[i*3+1],faces[i*3+2]);
	// }

	delete[] verts;
	delete[] faces;
	delete[] faces_per_vertex;

	return (fc&&vc) ? 0: 1;
}

static int test_bz_moveinto(){
	printf("\n*** test_bz_moveinto ***\n");
	Direct d(3.0,3.0,9.0,PI/2,PI/2,2*PI/3);
	Reciprocal rlat = d.star();
	BrillouinZone bz(rlat);

	// printf("Face vectors\n"); bz.get_faces().print();

	std::default_random_engine generator;
	std::uniform_real_distribution<double> qrand(0.0,1.0); // [0,1)
	std::uniform_int_distribution<int> taurand(-10,10); // [0,10]

	double q[30], Q[30];
	int tau[30];
	for (int i=0; i<30; i++){
		q[i] = qrand(generator);
		tau[i] = taurand(generator);
		Q[i] = q[i]+tau[i];
	}
  LQVec<double> Qv(rlat,10,Q);
	LQVec<double> qv(rlat,10);
	LQVec<int> tv(rlat,10);

	if (!bz.moveinto(&Qv,&qv,&tv)) return 1;
	// printf("Q =\n"); Qv.print();
	// printf("q =\n"); qv.print();
	// printf("tau =\n"); tv.print();
	// so the following is untested so far
	int numwrong = 0;
	for (int i=0; i<10; ++i) for (int j=0; j<3; ++j) {
	  if ( abs(qv.getvalue(i,j)+tv.getvalue(i,j) - q[i*3+j]-tau[i*3+j]) > 1e-15) {
      // printf("% 8.5f%+3d != % 8.5f%+3d (%g)\n",qv.getvalue(i,j),tv.getvalue(i,j),q[i*3+j],tau[i*3+j],qv.getvalue(i,j)+tv.getvalue(i,j)-q[i*3+j]-tau[i*3+j]);
      ++numwrong;
    }
	}
  if (numwrong) printf("%d of %u q and/or tau elements are wrong!\n",numwrong,3*qv.size());
	return numwrong;
}


static int test_bz_grid(){
	printf("\n*** test_bz_grid ***\n");
	size_t mg3N[3] = {3,4,2};
	MapGrid3<int> mg3;
	mg3.resize(mg3N);
	mg3.set_map();


	Direct dlat(9.,9.,9.,PI/2.,PI/2.,2*PI/3.);
	Reciprocal rlat= dlat.star();
	BrillouinZone bz(rlat);

	size_t N[3] = {10,10,0};

	BrillouinZoneGrid3 bzg(bz,N);
	bzg.print_map();

	ArrayVector<double> x = bzg.get_grid_x();
	ArrayVector<double> y = bzg.get_grid_y();
	ArrayVector<double> z = bzg.get_grid_z();

	// printf("Brillouin Zone Grid x values:\n"); x.print();
	// printf("Brillouin Zone Grid y values:\n"); y.print();
	// printf("Brillouin Zone Grid z values:\n"); z.print();

	ArrayVector<double> xyz = bzg.get_grid_xyz();
	// printf("Brillouin Zone Grid x,y,z values (in inverse Angstrom):\n");
	// xyz.print();

	ArrayVector<double> hkl = bzg.get_grid_hkl();
	// printf("Brillouin Zone Grid h,k,l values (in rlu):\n");
	// hkl.print();

	ArrayVector<double> map_xyz = bzg.get_mapped_xyz();


	return 0;
}
