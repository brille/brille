/* Copyright 2019 Greg Tucker
//
// This file is part of brille.
//
// brille is free software: you can redistribute it and/or modify it under the
// terms of the GNU Affero General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// brille is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with brille. If not, see <https://www.gnu.org/licenses/>.            */

#include "lattice.hpp"
#include "hall_symbol.hpp"

Lattice::Lattice(const double* latmat, const int h){
  double l[3]={0,0,0}, a[3]={0,0,0};
  latmat_to_lenang(latmat,3,1,l,a); // (3,1) -> assume row-ordered matrix
  this->set_len_pointer(l,1);
  this->set_ang_pointer(a,1, AngleUnit::radian);
  this->volume=this->calculatevolume();
  this->check_hall_number(h);
}
Lattice::Lattice(const double* lengths, const double* angles, const int h, const AngleUnit au){
  this->set_len_pointer(lengths,1);
  this->set_ang_pointer(angles,1, au);
  this->volume=this->calculatevolume();
  this->check_hall_number(h);
}
Lattice::Lattice(const double la, const double lb, const double lc, const double al, const double bl, const double cl, const int h):
  len{{la,lb,lc}}, ang{{al,bl,cl}} {
  // this->set_len_scalars(la,lb,lc);
  // this->set_ang_scalars(al,bl,cl, AngleUnit::not_provided);
  this->check_ang(AngleUnit::not_provided);
  this->volume = this->calculatevolume();
  this->check_hall_number(h);
}
Lattice::Lattice(const double* latmat, const std::string& itname, const std::string& choice){
  double l[3]={0,0,0}, a[3]={0,0,0};
  latmat_to_lenang(latmat,3,1,l,a);
  this->set_len_pointer(l,1);
  this->set_ang_pointer(a,1, AngleUnit::radian);
  this->volume=this->calculatevolume();
  this->check_IT_name(itname, choice);
}
Lattice::Lattice(const double *lengths, const double *angles, const std::string& itname, const std::string& choice, const AngleUnit au){
  this->set_len_pointer(lengths,1);
  this->set_ang_pointer(angles,1, au);
  this->volume=this->calculatevolume();
  this->check_IT_name(itname, choice);
}
Lattice::Lattice(const double la, const double lb, const double lc, const double al, const double bl, const double cl, const std::string& itname, const std::string& choice):
len({{la,lb,lc}}), ang({{al,bl,cl}}) {
  // this->set_len_scalars(la,lb,lc);
  // this->set_ang_scalars(al,bl,cl, AngleUnit::not_provided);
  this->check_ang(AngleUnit::not_provided);
  this->volume = this->calculatevolume();
  this->check_IT_name(itname, choice);
}
// void Lattice::set_len_scalars(const double a, const double b, const double c){
//   this->len[0] = a;
//   this->len[1] = b;
//   this->len[2] = c;
// }
// void Lattice::set_ang_scalars(const double a, const double b, const double g, const AngleUnit angle_unit){
//   this->ang[0] = a;
//   this->ang[1] = b;
//   this->ang[2] = g;
//   this->check_ang(angle_unit);
// }
void Lattice::check_ang(const AngleUnit angle_unit){
  AngleUnit au = angle_unit;
  if (au == AngleUnit::not_provided){
    double minang = (std::numeric_limits<double>::max)();
    double maxang = (std::numeric_limits<double>::lowest)();
    for (int i=0; i<3; ++i){
      if (this->ang[i] < minang) minang = this->ang[i];
      if (this->ang[i] > maxang) maxang = this->ang[i];
    }
    if (minang < 0.) throw std::runtime_error("Unexpected negative inter-facial cell angle");
    au = (maxang < 2.0*PIOVERTWO) ? AngleUnit::radian : AngleUnit::degree;
  }
  if (au != AngleUnit::radian){
    double conversion = PIOVERTWO*((AngleUnit::degree == au) ? 1.0/90.0 : 2.0);
    for (int i=0;i<3;i++) this->ang[i] *= conversion;
  }
}
void Lattice::check_hall_number(const int h){
  this->spg = Spacegroup(h); // if h is invalid the next three lines might fail`
  this->ptg = this->spg.get_pointgroup();
  this->spgsym = this->spg.get_spacegroup_symmetry();
  this->ptgsym = this->spg.get_pointgroup_symmetry();
}
void Lattice::check_IT_name(const std::string& itname, const std::string& choice){
  int hall_number = string_to_hall_number(itname, choice);
  if (hall_number > 0 && hall_number < 531){
    this->check_hall_number(hall_number);
  } else {
    // maybe itname is actually a non-standard Hall symbol?
    // possibly leave Spacegroup and Pointgroup (partially) unset
    HallSymbol hs(itname);
    if (hs.validate()){
      Symmetry generators = hs.get_generators();
      this->spgsym = generators.generate();
      this->spg.set_hall_symbol(hs.to_ascii()); // use a standardized form
      this->spg.set_bravais_type(hs.getl());
    } else {
      // last-ditch effort: maybe x,y,z notation Seitz matrices were passed?
      std::istringstream stream(itname);
      std::vector<Motion<int,double>> motions;
      for (std::string m; std::getline(stream, m, ';'); ) motions.push_back(Motion<int,double>(m));
      Symmetry mgens(motions);
      Motion<int,double> mone; // initalised to {𝟙|0}
      if (!mgens.has(mone)) mgens.add(mone); // make sure {𝟙|0}≡E is present
      this->spgsym = mgens.generate();
    }
    this->ptgsym = PointSymmetry(get_unique_rotations(this->spgsym.getallr(),0));
  }
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
  double *a = this->ang.data();
  double *l = this->len.data();
  tmp *= sin(( a[0] +a[1] +a[2])/2.0);
  tmp *= sin((-a[0] +a[1] +a[2])/2.0);
  tmp *= sin(( a[0] -a[1] +a[2])/2.0);
  tmp *= sin(( a[0] +a[1] -a[2])/2.0);
  tmp = sqrt(tmp);
  tmp *= 2*l[0]*l[1]*l[2];
  this->volume = tmp;
  if (std::isnan(tmp)){
    std::string msg = "Invalid lattice unit cell ";
    if (std::isnan(this->unitvolume())){
      msg += "angles [";
      for (int i=0; i<3; ++i) msg += " " + std::to_string(a[0]/PI*180.);
      msg += " ]";
    } else {
      msg += "lengths [";
      for (int i=0; i<3; ++i) msg += " " + std::to_string(l[i]);
      msg += " ]";
    }
    throw std::domain_error(msg);
  }
  return tmp;
}
Lattice Lattice::inner_star() const {
  std::array<double,3> lstar, astar, cosang, sinang;
  for (size_t i=0; i<3u; ++i){
    cosang[i] = std::cos(this->ang[i]);
    sinang[i] = std::sin(this->ang[i]);
  }
  for (size_t i=0; i<3u; ++i){
    size_t j = (i+1)%3;
    size_t k = (i+2)%3;
    lstar[i] = 2*PI*sinang[i]*this->len[j]*this->len[k]/this->volume;
    astar[i] = std::acos((cosang[j]*cosang[k]-cosang[i])/(sinang[j]*sinang[k]));
  }
  return Lattice(lstar, astar, this->spg, this->spgsym, this->ptg, this->ptgsym, this->basis);
}
void Lattice::get_metric_tensor(double * mt) const {
  const double *a = this->ang.data();
  const double *l = this->len.data();

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
  double *tmp = nullptr;
  tmp = new double[9]();
  this->get_metric_tensor(tmp);
  matrix_inverse(mt,tmp);
  delete[] tmp;
}

bool Lattice::issame(const Lattice& lat) const{
  return approx_vector(3, this->ang.data(), lat.ang.data()) && approx_vector(3, this->len.data(), lat.len.data());
}

bool Lattice::isapprox(const Lattice& lat) const {
  return this->ispermutation(lat)==0 ? false : true;
}
int Lattice::ispermutation(const Lattice& lat) const {
  double a[3], l[3];
  int i, j, ap[3]{0,2,1}; // for anti-permutations
  for (j=0; j<3; ++j){
    for (i=0; i<3; ++i){
      a[i] = this->ang[(i+j)%3];
      l[i] = this->len[(i+j)%3];
    }
    if (approx_vector(3, a, lat.ang.data()) && approx_vector(3, l, lat.len.data()))
      return j+1;
  }
  for (j=0; j<3; ++j){
    for (i=0; i<3; ++i){
      a[i] = this->ang[ap[(i+j)%3]];
      l[i] = this->len[ap[(i+j)%3]];
    }
    if (approx_vector(3, a, lat.ang.data()) && approx_vector(3, l, lat.len.data()))
      return -j-1;
  }
  return 0;
}


void Lattice::print(){
  printf("(%g %g %g) " ,this->len[0], this->len[1], this->len[2]);
  printf("(%g %g %g)\n",this->ang[0]/PI*180, this->ang[1]/PI*180, this->ang[2]/PI*180);
}
std::string lattice2string(const Lattice& l, const std::string& lenunit="", const std::string& angunit="°"){
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
  this->get_lattice_matrix(latmat, 3u, 1u);
}
template<class I> void Direct::get_lattice_matrix(double *latmat, std::vector<I>& strides) const {
  this->get_lattice_matrix(latmat, static_cast<size_t>(strides[0]/sizeof(double)), static_cast<size_t>(strides[1]/sizeof(double)));
}
void Direct::get_lattice_matrix(double *latmat, const size_t c, const size_t r) const{
  // define the lattice basis vectors using the same convention as spglib
  double c0=std::cos(this->ang[0]), c1=std::cos(this->ang[1]), c2=std::cos(this->ang[2]);
  double s2=std::sin(this->ang[2]);

  double xhat[3]={1,0,0}, yhat[3]={c2,s2,0};
  double zhat[3]={c1, (c0-c2*c1)/s2, this->unitvolume()/s2};

  for (int i=0;i<3;++i){
    // latmat[i*c+0*r] = this->len[0]*xhat[i];
    // latmat[i*c+1*r] = this->len[1]*yhat[i];
    // latmat[i*c+2*r] = this->len[2]*zhat[i];
    latmat[0*c+i*r] = this->len[0]*xhat[i]; // ̂x direction, first row
    latmat[1*c+i*r] = this->len[1]*yhat[i]; // ̂y direction, second row
    latmat[2*c+i*r] = this->len[2]*zhat[i]; // ̂z direction, third row
  }
}
void Direct::get_xyz_transform(double* fromxyz) const{
  this->get_xyz_transform(fromxyz, 3u, 1u);
}
template<class I> void Direct::get_xyz_transform(double* x, std::vector<I>& s) const {
  this->get_xyz_transform(x, static_cast<size_t>(s[0]/sizeof(double)), static_cast<size_t>(s[1]/sizeof(double)));
}
void Direct::get_xyz_transform(double *toxyz, const size_t c, const size_t r) const {
  // there are infinite possibilities for your choice of axes.
  // the original spglib used x along a and y in the (a,b) plane
  // here we're going with x along astar and y in the (astar, bstar) plane -- this is the "B-matrix"
  double B[9], t[9];
  this->star().get_B_matrix(B, 3u, 1u);
  // we want toxyz to be 2*pi*transpose(inverse(B));
  // use t as a buffer
  matrix_inverse(t,B);
  matrix_transpose(t); // in-place transpose using swap
  for (int i=0; i<3; ++i) for (int j=0; j<3; ++j) toxyz[i*c+j*r] = 2*PI*t[3*i+j];
}
void Direct::get_inverse_xyz_transform(double* fromxyz) const{
  this->get_inverse_xyz_transform(fromxyz, 3u, 1u);
}
template<class I> void Direct::get_inverse_xyz_transform(double* x, std::vector<I>& s) const {
  this->get_inverse_xyz_transform(x, static_cast<size_t>(s[0]/sizeof(double)), static_cast<size_t>(s[1]/sizeof(double)));
}
void Direct::get_inverse_xyz_transform(double *fromxyz, const size_t c, const size_t r) const{
  double toxyz[9], tmp[9];
  this->get_xyz_transform(toxyz,3u,1u);
  if(!matrix_inverse(tmp, toxyz))
    throw std::runtime_error("xyz_transform matrix has zero determinant");
  for (size_t i=0; i<3u; ++i) for (size_t j=0; j<3u; ++j) fromxyz[i*c+j*r] = tmp[i*3+j];
}
void Reciprocal::get_lattice_matrix(double *latmat) const{
  this->get_lattice_matrix(latmat, 3u, 1u);
}
template<class I> void Reciprocal::get_lattice_matrix(double* x, std::vector<I>& strides) const {
  this->get_lattice_matrix(x, static_cast<size_t>(strides[0]/sizeof(double)), static_cast<size_t>(strides[1]/sizeof(double)));
}
void Reciprocal::get_lattice_matrix(double *latmat, const size_t c, const size_t r) const{
  Direct d = this->star();
  double m[9], tmp[9];
  d.get_lattice_matrix(m, 3u, 1u);
  matrix_inverse(tmp,m);
  matrix_transpose(tmp);
  for (int i=0; i<9; ++i) tmp[i]*=2*PI;
  for (size_t i=0; i<3u; ++i) for (size_t j=0; j<3u; ++j) latmat[c*i+r*j] = tmp[3*i+j];
}
void Reciprocal::get_B_matrix(double* B) const {
  this->get_B_matrix(B, 3u, 1u);
}
template<class I> void Reciprocal::get_B_matrix(double* B, std::vector<I>& strides) const {
  this->get_B_matrix(B, static_cast<size_t>(strides[0]/sizeof(double)), static_cast<size_t>(strides[1]/sizeof(double)));
}
void Reciprocal::get_B_matrix(double *B, const size_t c, const size_t r) const {
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

  // a-star along x -- first column (constant row index = 0)
  B[0*c+0*r] = this->get_a();
  B[1*c+0*r] = 0.0;
  B[2*c+0*r] = 0.0;
  // b-star in the x-y plane -- second column
  B[0*c+1*r] = this->get_b()*cos(this->get_gamma());
  B[1*c+1*r] = this->get_b()*asg;
  B[2*c+1*r] = 0.0;
  // and c-star -- third column
  B[0*c+2*r] = this->get_c()*cos(this->get_beta());
  B[1*c+2*r] = -1.0*(this->get_c())*asb*cos(d.get_alpha());
  B[2*c+2*r] = 2*PI/d.get_c();
}
void Reciprocal::get_xyz_transform(double *toxyz) const { this->get_xyz_transform(toxyz, 3u, 1u); }
template<class I> void Reciprocal::get_xyz_transform(double* toxyz, std::vector<I>& strides) const {
  this->get_xyz_transform(toxyz, static_cast<size_t>(strides[0]/sizeof(double)), static_cast<size_t>(strides[1]/sizeof(double)));
}
void Reciprocal::get_xyz_transform(double *toxyz, const size_t c, const size_t r) const {
  // there are infinite possibilities for your choice of axes.
  // the original spglib used x along a and y in the (a,b) plane
  // here we're going with x along astar and y in the (astar, bstar) plane -- this is the "B-matrix"
  this->get_B_matrix(toxyz, c, r);
}
void Reciprocal::get_inverse_xyz_transform(double *fromxyz) const{
  this->get_inverse_xyz_transform(fromxyz, 3u, 1u);
}
template<class I> void Reciprocal::get_inverse_xyz_transform(double* x, std::vector<I>& strides) const {
  this->get_inverse_xyz_transform(x, static_cast<size_t>(strides[0]/sizeof(double)), static_cast<size_t>(strides[1]/sizeof(double)));
}
void Reciprocal::get_inverse_xyz_transform(double *fromxyz, const size_t c, const size_t r) const {
  double toxyz[9], tmp[9];
  this->get_xyz_transform(toxyz, 3u, 1u);
  if(!matrix_inverse(tmp, toxyz))
    throw std::runtime_error("xyz_transform matrix has zero determinant");
  for (size_t i=0; i<3u; ++i) for(size_t j=0; j<3u; ++j) fromxyz[i*c+j*r] = tmp[i*3+j];
}

// We have to define these separately from Lattice since Lattice.star() doesn't exist.
bool Direct::isstar(const Reciprocal& latt) const{
  // two lattices are the star of each other if the star of one is the same as the other
  // we need to check both ways in case rounding errors in the inversion cause a problem
  return ( this->issame( latt.star() ) || latt.issame( this->star() ) );
}
bool Reciprocal::isstar(const Direct& latt) const{
  // two lattices are the star of each other if the star of one is the same as the other
  // we need to check both ways in case rounding errors in the inversion cause a problem
  return ( this->issame( latt.star() ) || latt.issame( this->star() ) );
}

//bool Direct::issame(const Reciprocal) const {printf("call to Direct::issame(Reciprocal)\n"); return false;}
bool Direct::isstar(const Direct&) const {return false;}
//bool Reciprocal::issame(const Direct) const {printf("call to Reciprocal::issame(Direct)\n"); return false;}
bool Reciprocal::isstar(const Reciprocal&) const {return false;}

void Direct::print(){
  printf("(%g %g %g)A " ,this->len[0], this->len[1], this->len[2]);
  printf("(%g %g %g)\n",this->ang[0]/PI*180, this->ang[1]/PI*180, this->ang[2]/PI*180);
}
void Reciprocal::print(){
  printf("(%g %g %g)/A " ,this->len[0], this->len[1], this->len[2]);
  printf("(%g %g %g)\n",this->ang[0]/PI*180, this->ang[1]/PI*180, this->ang[2]/PI*180);
}


Direct Direct::primitive(void) const{
  PrimitiveTransform P(this->spg.get_bravais_type());
  if (P.does_anything()){
    double plm[9], lm[9];
    this->get_lattice_matrix(lm); // now returns *row* vectors!
    // The transformation matrix P gives us the primitive basis column-vector
    // matrix Aₚ from the standard basis column-vector matrix Aₛ by
    // Aₚ = Aₛ P. But here we have Aₛᵀ and want Aₚᵀ:
    // Aₚᵀ = (Aₛ P)ᵀ = Pᵀ Aₛᵀ
    // So we need the transpose of P from the PrimitiveTransform object:
    multiply_matrix_matrix<double,double,double,3>(plm,P.get_Pt().data(),lm);
    return Direct(plm);
  }
  return *this;
}
Reciprocal Reciprocal::primitive(void) const{
  // We could try to use the fact that for column vectors Bₚ = Bₛ (P⁻¹)ᵀ
  // but it's probably safer to continue using the reciprocal of the primitive
  // of the reciprocal of this lattice.
  return this->star().primitive().star();
}
