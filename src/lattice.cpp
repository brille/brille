#include <math.h>
#include "linear_algebra.h"
#include "lattice.h"
#include <string>

Lattice::Lattice(const double *lengths, const double *angles){
	for (int i=0;i<3;i++){
		this->len[i] = lengths[i];
		this->ang[i] = angles[i];
	}
	volume=calculatevolume();
}
Lattice::Lattice(double la, double lb, double lc, double al, double bl, double cl, double vol){
	this->len[0] = la;
	this->len[1] = lb;
	this->len[2] = lc;
	this->ang[0] = al;
	this->ang[1] = bl;
	this->ang[2] = cl;
	this->volume = vol;
}
Lattice::Lattice(double la, double lb, double lc, double al, double bl, double cl){
	this->len[0] = la;
	this->len[1] = lb;
	this->len[2] = lc;
	this->ang[0] = al;
	this->ang[1] = bl;
	this->ang[2] = cl;
	this->volume = calculatevolume();
}

double Lattice::calculatevolume(){
	double tmp=1;
	double *a = this->ang;
	double *l = this->len;
	tmp *= sin(( a[0] +a[1] +a[2])/2.0);
	tmp *= sin((-a[0] +a[1] +a[2])/2.0);
	tmp *= sin(( a[0] -a[1] +a[2])/2.0);
	tmp *= sin(( a[0] +a[1] -a[2])/2.0);
	tmp = sqrt(tmp);
	tmp *= 2*l[0]*l[1]*l[2];
	this->volume = tmp;
	return tmp;
}
Lattice Lattice::inner_star() const {
	const double *a = this->ang;
	const double *l = this->len;

	double sas, sbs, scs;
	sas = 2*PI*l[1]*l[2]*sin(a[0])/this->volume;
	sbs = 2*PI*l[2]*l[0]*sin(a[1])/this->volume;
	scs = 2*PI*l[0]*l[1]*sin(a[2])/this->volume;

	double cosa, cosb, cosc, sina, sinb, sinc;
	cosa = cos(a[0]);
	cosb = cos(a[1]);
	cosc = cos(a[2]);
	sina = sin(a[0]);
	sinb = sin(a[1]);
	sinc = sin(a[2]);
	double saa, sbb, scc;
	saa = acos( (cosb*cosc-cosa)/(sinb*sinc) );
	sbb = acos( (cosc*cosa-cosb)/(sinc*sina) );
	scc = acos( (cosa*cosb-cosc)/(sina*sinb) );

	return Lattice(sas, sbs, scs, saa, sbb, scc, 8*PICUBED/this->volume);
}
void Lattice::get_metric_tensor(double * mt) const {
	const double *a = this->ang;
	const double *l = this->len;

	double cosa, cosb, cosc;
	cosa = cos(a[0]);
	cosb = cos(a[1]);
	cosc = cos(a[2]);

	mt[0] = l[0]*l[0];
	mt[1] = l[1]*l[0]*cosc;
	mt[2] = l[2]*l[0]*cosb;

	mt[3] = l[0]*l[1]*cosc;
	mt[4] = l[1]*l[1];
	mt[5] = l[2]*l[1]*cosa;

	mt[6] = l[0]*l[2]*cosb;
	mt[7] = l[1]*l[2]*cosa;
	mt[8] = l[2]*l[2];
}
void Lattice::get_covariant_metric_tensor(double *mt) const {
	this->get_metric_tensor(mt);
}
void Lattice::get_contravariant_metric_tensor(double *mt) const {
	double *tmp = new double[9]();
	this->get_metric_tensor(tmp);
	matrix_inverse(mt,tmp);
	delete[] tmp;
}

bool Lattice::issame(const Lattice lat) const{
	const double *a1 = this->ang;
	const double *l1 = this->len;
	const double *a2 = lat.ang;
	const double *l2 = lat.len;
	double eps = 2*std::numeric_limits<double>::epsilon();
	for (int i=0; i<3; i++){
		if ( my_abs(l1[i]-l2[i]) > eps*my_abs(l1[i]+l2[i]) ) return false;
		if ( my_abs(a1[i]-a2[i]) > eps*my_abs(a1[i]+a2[i]) ) return false;
	}
	return true;
}

void Lattice::print(){
	printf("(%g %g %g) " ,this->len[0], this->len[1], this->len[2]);
	printf("(%g %g %g)\n",this->ang[0]/PI*180, this->ang[1]/PI*180, this->ang[2]/PI*180);
}
std::string lattice2string(const Lattice& l, const std::string lenunit="", const std::string angunit="°"){
	std::string repr;
	repr = "(" + std::to_string(l.get_a()) + " "
						 + std::to_string(l.get_b()) + " "
						 + std::to_string(l.get_c()) + ")" + lenunit
			 +" (" + std::to_string(l.get_alpha()/PI*180) + " "
						 + std::to_string(l.get_beta()/PI*180)  + " "
						 + std::to_string(l.get_gamma()/PI*180) + ")" + angunit;
	return repr;
}
std::string Lattice::string_repr(){ return lattice2string(*this); }
std::string Direct::string_repr(){return lattice2string(*this,"Å");}
std::string Reciprocal::string_repr(){return lattice2string(*this,"Å⁻¹");}
// 	std::string repr;
// 	repr = "(" + std::to_string(this->len[0]) + " "
// 	           + std::to_string(this->len[1]) + " "
// 						 + std::to_string(this->len[2]) + ")"
// 			 +" (" + std::to_string(this->ang[0]/PI*180) + " "
// 			       + std::to_string(this->ang[1]/PI*180) + " "
// 						 + std::to_string(this->ang[2]/PI*180) + ")°";
// 	return repr;
// }

Reciprocal Direct::star() const {
	return Reciprocal(this->inner_star());
}
Direct Reciprocal::star() const {
	return Direct(this->inner_star());
}



void Direct::get_xyz_transform(double *toxyz) const {
	// there are infinite possibilities for your choice of axes.
	// the original spglib used x along a and y in the (a,b) plane
	// here we're going with x along astar and y in the (astar, bstar) plane -- this is the "B-matrix"
	Reciprocal r = this->star();
	double *B = new double[9]();
	r.get_B_matrix(B);
	// we want toxyz to be 2*pi*transpose(inverse(B));
	// use it as a buffer
	matrix_inverse(toxyz,B);
	matrix_transpose(toxyz); // in-place transpose using swap
	for (int i=0; i<9; i++) toxyz[i] *= 2*PI;
	delete[] B;
}

void Reciprocal::get_B_matrix(double *B) const {
	//Calculate the B-matrix as in Acta Cryst. (1967). 22, 457
	// http://dx.doi.org/10.1107/S0365110X67000970
	Direct d = this->star();

	double asb, asg;
	asb = sin(this->get_beta());
	if (asb<0) asb*=-1.0;
	asg = sin(this->get_gamma());
	if (asg<0) asg*=-1.0;

	// Be careful about indexing B. Make sure each vector goes in as a column, not a row!
	// if you mix this up, you'll only notice for non-orthogonal space groups

	// a-star along x -- first column = [0,3,6] in C ([0,1,2] in MATLAB)
	B[0] = this->get_a();
	B[3] = 0.0;
	B[6] = 0.0;

	// b-star in the x-y plane -- second column = [1,4,7] in C ([3,4,5] in MATLAB)
	B[1] = this->get_b()*cos(this->get_gamma());
	B[4] = this->get_b()*asg;
	B[7] = 0.0;

	// and c-star -- third column = [2,5,8] in C ([6,7,8] in MATLAB)
	B[2] = this->get_c()*cos(this->get_beta());
	B[5] = -1.0*(this->get_c())*asb*cos(d.get_alpha());
	B[8] = 2*PI/d.get_c();
}
void Reciprocal::get_xyz_transform(double *toxyz) const {
	// there are infinite possibilities for your choice of axes.
	// the original spglib used x along a and y in the (a,b) plane
	// here we're going with x along astar and y in the (astar, bstar) plane -- this is the "B-matrix"
	this->get_B_matrix(toxyz);
}

// We have to define these separately from Lattice since Lattice.star() doesn't exist.
bool Direct::isstar(const Reciprocal latt) const{
	// two lattices are the star of each other if the star of one is the same as the other
	// we need to check both ways in case rounding errors in the inversion cause a problem
	return ( this->issame( latt.star() ) || latt.issame( this->star() ) );
}
bool Reciprocal::isstar(const Direct latt) const{
	// two lattices are the star of each other if the star of one is the same as the other
	// we need to check both ways in case rounding errors in the inversion cause a problem
	return ( this->issame( latt.star() ) || latt.issame( this->star() ) );
}

//bool Direct::issame(const Reciprocal) const {printf("call to Direct::issame(Reciprocal)\n"); return false;}
bool Direct::isstar(const Direct) const {return false;}
//bool Reciprocal::issame(const Direct) const {printf("call to Reciprocal::issame(Direct)\n"); return false;}
bool Reciprocal::isstar(const Reciprocal) const {return false;}

void Direct::print(){
	printf("(%g %g %g)A " ,this->len[0], this->len[1], this->len[2]);
	printf("(%g %g %g)\n",this->ang[0]/PI*180, this->ang[1]/PI*180, this->ang[2]/PI*180);
}
void Reciprocal::print(){
	printf("(%g %g %g)/A " ,this->len[0], this->len[1], this->len[2]);
	printf("(%g %g %g)\n",this->ang[0]/PI*180, this->ang[1]/PI*180, this->ang[2]/PI*180);
}
