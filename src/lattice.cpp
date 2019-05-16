#include <cmath>
#include "linear_algebra.h"
#include "lattice.h"
#include <string>

Lattice::Lattice(const double* latmat, const int h){
	double l[3]={0,0,0}, a[3]={0,0,0}, n[9];
	// compute the dot product of each column with itself
	for (int i=0; i<3; ++i)	for (int j=0; j<3; ++j)	l[i] += latmat[i*3+j]*latmat[i*3+j];
	// the lattice vector lengths are the square root of this
	for (int i=0; i<3; ++i) l[i] = std::sqrt(l[i]);
	// normalize the column vectors, leaving only angle information
	for (int i=0; i<3; ++i) for (int j=0; j<3; ++j) n[i*3+j] = latmat[i*3+j]/l[i];
	// take the dot product between cyclically permuted columns: 0=1⋅2, 1=2⋅0, 2=0⋅1
	for (int i=0; i<3; ++i)	for (int j=0; j<3; ++j)	a[i] += n[3*((i+1)%3)+j]*n[3*((i+2)%3)+j];
	// the lattice angles are the arccosines of these dot products of normalized lattice vectors
	for (int i=0; i<3; ++i) a[i] = std::acos(a[i]);

	this->set_len_pointer(l);
	this->set_ang_pointer(a);
	this->volume=this->calculatevolume();
	this->check_hall_number(h);
}

Lattice::Lattice(const double* lengths, const double* angles, const int h){
	this->set_len_pointer(lengths);
	this->set_ang_pointer(angles);
  this->volume=this->calculatevolume();
	this->check_hall_number(h);
}
Lattice::Lattice(const double la, const double lb, const double lc, const double al, const double bl, const double cl, const int h){
	this->set_len_scalars(la,lb,lc);
	this->set_ang_scalars(al,bl,cl);
	this->volume = this->calculatevolume();
  this->check_hall_number(h);
}
Lattice::Lattice(const double *lengths, const double *angles, const std::string itname){
	this->set_len_pointer(lengths);
	this->set_ang_pointer(angles);
  this->volume=this->calculatevolume();
	this->check_IT_name(itname);
}
Lattice::Lattice(const double la, const double lb, const double lc, const double al, const double bl, const double cl, const std::string itname){
	this->set_len_scalars(la,lb,lc);
	this->set_ang_scalars(al,bl,cl);
	this->volume = this->calculatevolume();
  this->check_IT_name(itname);
}
void Lattice::set_len_pointer(const double *lvec){
	for (int i=0;i<3;i++) this->len[i] = lvec[i];
}
void Lattice::set_ang_pointer(const double *avec){
	for (int i=0;i<3;i++)	this->ang[i] = avec[i];
}
void Lattice::set_len_scalars(const double a, const double b, const double c){
	this->len[0] = a;
	this->len[1] = b;
	this->len[2] = c;
}
void Lattice::set_ang_scalars(const double a, const double b, const double g){
	this->ang[0] = a;
	this->ang[1] = b;
	this->ang[2] = g;
}
void Lattice::check_hall_number(const int h){
	this->hall = hall_number_ok(h) ? h : 0;
}
void Lattice::check_IT_name(const std::string itname){
	this->hall = spgdb_international_to_hall_number(itname);
}
double Lattice::unitvolume() const{
	// The volume of a parallelpiped with unit length sides and our body angles
	double c[3];
	for (int i=0;i<3;++i) c[i]=std::cos(this->ang[i]);
	double sos=0, prd=2;
	for (int i=0;i<3;++i){
		sos+=c[i]*c[i];
		prd*=c[i];
	}
	return std::sqrt( 1 - sos + prd );
}
double Lattice::calculatevolume(){
	// we could replace this by l[0]*l[1]*l[2]*this->unitvolume()
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

	return Lattice(sas, sbs, scs, saa, sbb, scc, this->hall);
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

Reciprocal Direct::star() const {
	return Reciprocal(this->inner_star());
}
Direct Reciprocal::star() const {
	return Direct(this->inner_star());
}

void Direct::get_lattice_matrix(double *latmat) const{
	// define the lattice basis vectors using the same convention as spglib
	double c0=std::cos(this->ang[0]), c1=std::cos(this->ang[1]), c2=std::cos(this->ang[2]);
	double s2=std::sin(this->ang[2]);

	double xhat[3]={1,0,0}, yhat[3]={c2,s2,0};
	double zhat[3]={c1, (c0-c2*c1)/s2, this->unitvolume()/s2};

	for (int i=0;i<3;++i){
		latmat[i*3+0] = this->len[0]*xhat[i];
		latmat[i*3+1] = this->len[1]*yhat[i];
		latmat[i*3+2] = this->len[2]*zhat[i];
	}

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

void Reciprocal::get_lattice_matrix(double *latmat) const{
	Direct d = this->star();
	double m[9];
	d.get_lattice_matrix(m);
	matrix_inverse(latmat,m);
	matrix_transpose(latmat);
	for (int i=0; i<9; ++i) latmat[i]*=2*PI;
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


Direct Direct::primitive(void) const{
	double plm[9], lm[9];
	SpacegroupType spgt = spgdb_get_spacegroup_type(this->hall);
	PrimitiveTransform P(spgt.centering);
	if (P.does_anything()){
		this->get_lattice_matrix(lm);
		std::array<double,9> Parray = P.get_to_primitive();
		multiply_matrix_matrix<double,double,double,3>(plm,lm,Parray.data());
		return Direct(plm);
	}
	return *this;
}
Reciprocal Reciprocal::primitive(void) const{
	return this->star().primitive().star();
}
