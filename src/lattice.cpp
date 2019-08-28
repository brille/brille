#include "lattice.h"

Lattice::Lattice(const double* latmat, const int h){
  double l[3]={0,0,0}, a[3]={0,0,0};
  latmat_to_lenang(latmat,3,1,l,a); // (3,1) -> assume row-ordered matrix
  this->set_len_pointer(l,1);
  this->set_ang_pointer(a,1);
  this->volume=this->calculatevolume();
  this->check_hall_number(h);
}
Lattice::Lattice(const double* lengths, const double* angles, const int h){
  this->set_len_pointer(lengths,1);
  this->set_ang_pointer(angles,1);
  this->volume=this->calculatevolume();
  this->check_hall_number(h);
}
Lattice::Lattice(const double la, const double lb, const double lc, const double al, const double bl, const double cl, const int h){
  this->set_len_scalars(la,lb,lc);
  this->set_ang_scalars(al,bl,cl);
  this->volume = this->calculatevolume();
  this->check_hall_number(h);
}
Lattice::Lattice(const double* latmat, const std::string itname){
  double l[3]={0,0,0}, a[3]={0,0,0};
  latmat_to_lenang(latmat,3,1,l,a);
  this->set_len_pointer(l,1);
  this->set_ang_pointer(a,1);
  this->volume=this->calculatevolume();
  this->check_IT_name(itname);
}
Lattice::Lattice(const double *lengths, const double *angles, const std::string itname){
  this->set_len_pointer(lengths,1);
  this->set_ang_pointer(angles,1);
  this->volume=this->calculatevolume();
  this->check_IT_name(itname);
}
Lattice::Lattice(const double la, const double lb, const double lc, const double al, const double bl, const double cl, const std::string itname){
  this->set_len_scalars(la,lb,lc);
  this->set_ang_scalars(al,bl,cl);
  this->volume = this->calculatevolume();
  this->check_IT_name(itname);
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
  this->hall = international_string_to_hall_number(itname);
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
    if ( std::abs(l1[i]-l2[i]) > eps*std::abs(l1[i]+l2[i]) ) return false;
    if ( std::abs(a1[i]-a2[i]) > eps*std::abs(a1[i]+a2[i]) ) return false;
  }
  return true;
}

bool Lattice::isapprox(const Lattice lat, const double tol) const {
  return this->ispermutation(lat, tol)==0 ? false : true;
}
int Lattice::ispermutation(const Lattice lat, const double tol) const {
  const double *a1 = this->ang;
  const double *l1 = this->len;
  const double *a2 = lat.ang;
  const double *l2 = lat.len;
  double dif, sum;
  bool perm_equiv = true;
  // check if the a,b,c axes are permuted
  for (int j=0; j<3; ++j){
    perm_equiv = true;
    for (int i=0; i<3; i++){
      dif = std::abs(l1[i]-l2[(i+j)%3]);
      sum = std::abs(l1[i]+l2[(i+j)%3]);
      if ( dif > tol*sum || dif > tol ) perm_equiv = false;
      dif = std::abs(a1[i]-a2[(i+j)%3]);
      sum = std::abs(a1[i]+a2[(i+j)%3]);
      if ( dif > tol*sum || dif > tol ) perm_equiv = false;
    }
    if (perm_equiv) return j+1;
  }
  // antipermutations are possible too
  int ap[3]{0,2,1};
  for (int j=0; j<3; ++j){
    perm_equiv = true;
    for (int i=0; i<3; i++){
      dif = std::abs(l1[i]-l2[ap[(i+j)%3]]);
      sum = std::abs(l1[i]+l2[ap[(i+j)%3]]);
      if ( dif > tol*sum || dif > tol ) perm_equiv = false;
      dif = std::abs(a1[i]-a2[ap[(i+j)%3]]);
      sum = std::abs(a1[i]+a2[ap[(i+j)%3]]);
      if ( dif > tol*sum || dif > tol ) perm_equiv = false;
    }
    if (perm_equiv) return -j-1;
  }
  return 0;
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
  PrimitiveTransform P(this->hall);
  if (P.does_anything()){
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
