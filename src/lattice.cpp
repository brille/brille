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
#include "lattice.hpp"
#include "hall_symbol.hpp"
#include "approx.hpp"

using namespace brille;

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
Lattice::Lattice(const double la, const double lb, const double lc, const double al, const double bl, const double cl, const int h){
  double l[3]{la,lb,lc}, a[3]{al,bl,cl};
  this->set_len_pointer(l,1);
  this->set_ang_pointer(a,1, AngleUnit::not_provided);
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
Lattice::Lattice(const double la, const double lb, const double lc, const double al, const double bl, const double cl, const std::string& itname, const std::string& choice){
  double l[3]{la,lb,lc}, a[3]{al,bl,cl};
  this->set_len_pointer(l,1);
  this->set_ang_pointer(a,1, AngleUnit::not_provided);
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
    au = (maxang < brille::math::pi) ? AngleUnit::radian : AngleUnit::degree;
  }
  if (au != AngleUnit::radian){
    double conversion = brille::math::pi/((AngleUnit::degree == au) ? 180.0 : 1.0);
    for (int i=0;i<3;i++) this->ang[i] *= conversion;
  }
}
void Lattice::check_hall_number(const int h){
  Spacegroup spg(h); // if h is invalid the next three lines might fail
  this->bravais = spg.get_bravais_type();
  this->spgsym = spg.get_spacegroup_symmetry();
  this->ptgsym = spg.get_pointgroup_symmetry();
  this->check_parameter_symmetry();
}
void Lattice::check_IT_name(const std::string& itname, const std::string& choice){
  int hall_number = brille::string_to_hall_number(itname, choice);
  if (hall_number > 0 && hall_number < 531){
    this->check_hall_number(hall_number);
  } else {
    // maybe itname is actually a non-standard Hall symbol?
    HallSymbol hs(itname);
    Symmetry gens;
    if (hs.validate()){
      gens = hs.get_generators();
      bravais = hs.getl();
    } else {
      // last-ditch effort: maybe x,y,z notation Seitz matrices were passed?
      gens.from_ascii(itname);
      bravais = Bravais::P; // without any extra information this must be primitive
    }
    this->set_spacegroup_symmetry(gens);
  }
}
double Lattice::unitvolume() const{
  // The volume of a parallelpiped with unit length sides and our body angles
  // Part of equation 2.6.6 in https://www.iucr.org/__data/assets/pdf_file/0019/15823/22.pdf
  double sos=0, prd=2;
  for (int i=0;i<3;++i){
    sos+=cosines[i]*cosines[i];
    prd*=cosines[i];
  }
  return std::sqrt( 1 - sos + prd );
}
double Lattice::calculatevolume(){
  // from equation 2.6.6 in https://www.iucr.org/__data/assets/pdf_file/0019/15823/22.pdf
  volume = len[0]*len[1]*len[2] * this->unitvolume(); // is this off by 2pi?
  if (std::isnan(volume)){
    std::string msg = "Invalid lattice unit cell ";
    if (std::isnan(this->unitvolume())){
      msg += "angles [";
      for (int i=0; i<3; ++i) msg += " " + std::to_string(ang[0]/brille::math::pi*180.);
      msg += " ]";
    } else {
      msg += "lengths [";
      for (int i=0; i<3; ++i) msg += " " + std::to_string(len[i]);
      msg += " ]";
    }
    throw std::domain_error(msg);
  }
  return volume;
}
Lattice Lattice::inner_star() const {
  std::array<double,3> lstar{}, astar{}, cstar{}, sstar{};
//  for (size_t i=0; i<3u; ++i){
//    cosang[i] = std::cos(this->ang[i]);
//    sinang[i] = std::sin(this->ang[i]);
//  }
  for (size_t i=0; i<3u; ++i){
    size_t j = (i+1)%3;
    size_t k = (i+2)%3;
    // the grouping of len[j]*len[k] is done to avoid rounding errors when
    // e.g., a* and b* *should be*  the same but order-of-operations causes
    // rounding errors
    lstar[i] = brille::math::two_pi * sines[i] * (len[j] * len[k]) / this->volume;
    cstar[i] = (cosines[j]*cosines[k]-cosines[i])/(sines[j]*sines[k]);
    astar[i] = std::acos(cstar[i]);
    sstar[i] = std::sin(astar[i]);
  }
  return Lattice(lstar, astar, cstar, sstar, this->bravais, this->spgsym, this->ptgsym, this->basis);
}
void Lattice::get_metric_tensor(double * mt) const {
  mt[0] = len[0]*len[0];
  mt[1] = len[1]*len[0]*cosines[2];
  mt[2] = len[2]*len[0]*cosines[1];

  mt[3] = len[0]*len[1]*cosines[2];
  mt[4] = len[1]*len[1];
  mt[5] = len[2]*len[1]*cosines[0];

  mt[6] = len[0]*len[2]*cosines[1];
  mt[7] = len[1]*len[2]*cosines[0];
  mt[8] = len[2]*len[2];
}
void Lattice::get_covariant_metric_tensor(double *mt) const {
  this->get_metric_tensor(mt);
}
void Lattice::get_contravariant_metric_tensor(double *mt) const {
  double *tmp = nullptr;
  tmp = new double[9]();
  this->get_metric_tensor(tmp);
  brille::utils::matrix_inverse(mt,tmp);
  delete[] tmp;
}

std::vector<double> Lattice::get_metric_tensor() const {
  std::vector<double> m(9);
  this->get_metric_tensor(m.data());
  return m;
}
std::vector<double> Lattice::get_covariant_metric_tensor() const {
  std::vector<double> m(9);
  this->get_covariant_metric_tensor(m.data());
  return m;
}
std::vector<double> Lattice::get_contravariant_metric_tensor() const {
  std::vector<double> m(9);
  this->get_contravariant_metric_tensor(m.data());
  return m;
}

bool Lattice::issame(const Lattice& lat) const{
  return brille::approx::vector(3, this->ang.data(), lat.ang.data()) && brille::approx::vector(3, this->len.data(), lat.len.data());
}

bool Lattice::isapprox(const Lattice& lat) const {
  return this->ispermutation(lat) != 0;
}
int Lattice::ispermutation(const Lattice& lat) const {
  double a[3], l[3];
  int i, j, ap[3]{0,2,1}; // for anti-permutations
  for (j=0; j<3; ++j){
    for (i=0; i<3; ++i){
      a[i] = this->ang[(i+j)%3];
      l[i] = this->len[(i+j)%3];
    }
    if (brille::approx::vector(3, a, lat.ang.data()) && brille::approx::vector(3, l, lat.len.data()))
      return j+1;
  }
  for (j=0; j<3; ++j){
    for (i=0; i<3; ++i){
      a[i] = this->ang[ap[(i+j)%3]];
      l[i] = this->len[ap[(i+j)%3]];
    }
    if (brille::approx::vector(3, a, lat.ang.data()) && brille::approx::vector(3, l, lat.len.data()))
      return -j-1;
  }
  return 0;
}


void Lattice::print(){
  printf("(%g %g %g) " ,this->len[0], this->len[1], this->len[2]);
  printf("(%g %g %g)\n",this->ang[0]/brille::math::pi*180, this->ang[1]/brille::math::pi*180, this->ang[2]/brille::math::pi*180);
}
std::string lattice2string(const Lattice& l, const std::string& lenunit="", const std::string& angunit="°"){
  std::string repr;
  repr = "(" + std::to_string(l.get_a()) + " "
             + std::to_string(l.get_b()) + " "
             + std::to_string(l.get_c()) + ")" + lenunit
       +" (" + std::to_string(l.get_alpha()/brille::math::pi*180) + " "
             + std::to_string(l.get_beta()/brille::math::pi*180)  + " "
             + std::to_string(l.get_gamma()/brille::math::pi*180) + ")" + angunit;
  return repr;
}
std::string Lattice::string_repr(){ return lattice2string(*this); }

Bravais Lattice::set_bravais_type(const Bravais b) {
    this->bravais = b;
    return this->get_bravais_type();
}

void Lattice::check_parameter_symmetry() {
  // the point-group symmetry applies to the basis vectors
  // (the space-group translations do not affect the vectors)
  // furthermore, only those operations with an order above 1 matter
  auto to_check = ptgsym.higher(1);
  std::array<int,3> a{1,0,0}, b{0,1,0}, c{0,0,1};
  // look for an operation that maps a to b, b to c, or a to c:
  bool a_to_b{false}, a_to_c{false}, b_to_a{false}, b_to_c{false}, c_to_a{false}, c_to_b{false};
  for (ind_t i=0; i < to_check.size(); ++i){
    auto op = to_check.get(i); // these are the operations for the Direct lattice, can we get away with ignoring that?
    std::array<int,3> t{0,0,0};
    utils::mul_mat_vec(t.data(), 3u, op.data(), a.data());
    if (t[0] == b[0] && t[1] == b[1] && t[2] == b[2]) a_to_b = true;
    if (t[0] == c[0] && t[1] == c[1] && t[2] == c[2]) a_to_c = true;
    utils::mul_mat_vec(t.data(), 3u, op.data(), b.data());
    if (t[0] == a[0] && t[1] == a[1] && t[2] == a[2]) b_to_a = true;
    if (t[0] == c[0] && t[1] == c[1] && t[2] == c[2]) b_to_c = true;
    utils::mul_mat_vec(t.data(), 3u, op.data(), c.data());
    if (t[0] == a[0] && t[1] == a[1] && t[2] == a[2]) c_to_a = true;
    if (t[0] == b[0] && t[1] == b[1] && t[2] == b[2]) c_to_b = true;
  }
  if (a_to_b ^ b_to_a || a_to_c ^ c_to_a || b_to_c ^ c_to_b){
    throw std::runtime_error("Only half of the symmetry found?!");
  }
  if ((a_to_b && b_to_c && !a_to_c) || (a_to_b && !b_to_c && a_to_c) || (!a_to_b && b_to_c && a_to_c)) {
    throw std::runtime_error("Two axis pair equivalencies requires all three pairs!");
  }
  if ((a_to_b && b_to_c && a_to_c) && (len[0] != len[1] || len[1] != len[2])) {
    auto abc = (len[0] + len[1] + len[2]) / 3;
    info_update("The provided symmetry requires equal basis vector lengths but not all of", len, " are the same. All will be replaced by ", abc);
    info_update("Check whether the provided angles are appropriate, ", ang, " radian");
    len[0] = abc;
    len[1] = abc;
    len[2] = abc;
  }
  if ((a_to_b && !b_to_c && !a_to_c) && (len[0] != len[1])){
    auto ab = (len[0] + len[1] ) / 2;
    info_update("The provided symmetry requires equal a & b basis vector lengths but ", len, " are not the same. All will be replaced by ", ab);
    info_update("Check whether the provided angles are appropriate, ", ang, " radian");
    info_update("with cos(ang)=", cosines, " and sin(ang)=", sines);
    len[0] = ab;
    len[1] = ab;
  }
  if ((!a_to_b && b_to_c && !a_to_c) && (len[1] != len[2])){
    auto bc = (len[1] + len[2] ) / 2;
    info_update("The provided symmetry requires equal b & c basis vector lengths but ", len, " are not the same. All will be replaced by ", bc);
    info_update("Check whether the provided angles are appropriate, ", ang, " radian");
    info_update("with cos(ang)=", cosines, " and sin(ang)=", sines);
    len[1] = bc;
    len[2] = bc;
  }
  if ((!a_to_b && !b_to_c && a_to_c) && (len[0] != len[2])){
    auto ac = (len[0] + len[2] ) / 2;
    info_update("The provided symmetry requires equal a & c basis vector lengths but ", len, " are not the same. All will be replaced by ", ac);
    info_update("Check whether the provided angles are appropriate, ", ang, " radian");
    info_update("with cos(ang)=", cosines, " and sin(ang)=", sines);
    len[0] = ac;
    len[2] = ac;
  }
}


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
//  double c0=std::cos(this->ang[0]), c1=std::cos(this->ang[1]), c2=std::cos(this->ang[2]);
//  double s2=std::sin(this->ang[2]);

  double x[3]{1,0,0};
  double y[3]{cosines[2],sines[2],0};
  double z[3]{cosines[1], (cosines[0]-cosines[2]*cosines[1])/sines[2], this->unitvolume()/sines[2]};

  for (int i=0;i<3;++i){
    // latmat[i*c+0*r] = this->len[0]*x[i];
    // latmat[i*c+1*r] = this->len[1]*y[i];
    // latmat[i*c+2*r] = this->len[2]*z[i];
    latmat[0*c+i*r] = len[0]* x[i]; // ̂x direction, first row
    latmat[1*c+i*r] = len[1]* y[i]; // ̂y direction, second row
    latmat[2*c+i*r] = len[2]* z[i]; // ̂z direction, third row
  }
}
void Direct::get_xyz_transform(double* toxyz) const{
  this->get_xyz_transform(toxyz, 3u, 1u);
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
  brille::utils::matrix_inverse(t,B);
  brille::utils::matrix_transpose(t); // in-place transpose using swap
  for (int i=0; i<3; ++i) for (int j=0; j<3; ++j) toxyz[i*c+j*r] = brille::math::two_pi*t[3*i+j];
}
std::vector<double> Direct::get_xyz_transform() const {
  std::vector<double> mat(9);
  this->get_xyz_transform(mat.data());
  return mat;
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
  if(!brille::utils::matrix_inverse(tmp, toxyz))
    throw std::runtime_error("xyz_transform matrix has zero determinant");
  for (size_t i=0; i<3u; ++i) for (size_t j=0; j<3u; ++j) fromxyz[i*c+j*r] = tmp[i*3+j];
}
std::vector<double> Direct::get_inverse_xyz_transform() const {
  std::vector<double> mat(9);
  this->get_inverse_xyz_transform(mat.data());
  return mat;
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
  brille::utils::matrix_inverse(tmp,m);
  brille::utils::matrix_transpose(tmp);
  for (int i=0; i<9; ++i) tmp[i]*=brille::math::two_pi;
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
  asb = sines[1];
  if (asb<0) asb*=-1.0;
  asg = sines[2];
  if (asg<0) asg*=-1.0;
  // Be careful about indexing B. Make sure each vector goes in as a column, not a row!
  // if you mix this up, you'll only notice for non-orthogonal space groups
  auto d_cosines = d.get_cosines();

  // a-star along x -- first column (constant row index = 0)
  B[0*c+0*r] = this->get_a();
  B[1*c+0*r] = 0.0;
  B[2*c+0*r] = 0.0;
  // b-star in the x-y plane -- second column
  B[0*c+1*r] = this->get_b()*cosines[2];
  B[1*c+1*r] = this->get_b()*asg;
  B[2*c+1*r] = 0.0;
  // and c-star -- third column
  B[0*c+2*r] = this->get_c()*cosines[1];
  B[1*c+2*r] = -1.0*(this->get_c())*asb*d_cosines[0]; //cos(d.get_alpha());
  B[2*c+2*r] = brille::math::two_pi/d.get_c();
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
std::vector<double> Reciprocal::get_xyz_transform() const {
  std::vector<double> mat(9);
  this->get_xyz_transform(mat.data());
  return mat;
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
  if(!brille::utils::matrix_inverse(tmp, toxyz))
    throw std::runtime_error("xyz_transform matrix has zero determinant");
  for (size_t i=0; i<3u; ++i) for(size_t j=0; j<3u; ++j) fromxyz[i*c+j*r] = tmp[i*3+j];
}
std::vector<double> Reciprocal::get_inverse_xyz_transform() const {
  std::vector<double> mat(9);
  this->get_inverse_xyz_transform(mat.data());
  return mat;
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
  printf("(%g %g %g)\n",this->ang[0]/brille::math::pi*180, this->ang[1]/brille::math::pi*180, this->ang[2]/brille::math::pi*180);
}
void Reciprocal::print(){
  printf("(%g %g %g)/A " ,this->len[0], this->len[1], this->len[2]);
  printf("(%g %g %g)\n",this->ang[0]/brille::math::pi*180, this->ang[1]/brille::math::pi*180, this->ang[2]/brille::math::pi*180);
}


Direct Direct::primitive() const{
  PrimitiveTransform P(this->bravais);
  if (P.does_anything()){
    double plm[9], lm[9];
    this->get_lattice_matrix(lm); // now returns *row* vectors!
    // The transformation matrix P gives us the primitive basis column-vector
    // matrix Aₚ from the standard basis column-vector matrix Aₛ by
    // Aₚ = Aₛ P. But here we have Aₛᵀ and want Aₚᵀ:
    // Aₚᵀ = (Aₛ P)ᵀ = Pᵀ Aₛᵀ
    // So we need the transpose of P from the PrimitiveTransform object:
    brille::utils::multiply_matrix_matrix<double,PrimitiveTraits::sixP,double,3>(plm,P.get_6Pt().data(),lm);
    for (double & i : plm) i /= 6.0; // correct for having used 6*P
    return Direct(plm);
  }
  return *this;
}
Reciprocal Reciprocal::primitive() const{
  // We could try to use the fact that for column vectors Bₚ = Bₛ (P⁻¹)ᵀ
  // but it's probably safer to continue using the reciprocal of the primitive
  // of the reciprocal of this lattice.
  return this->star().primitive().star();
}
