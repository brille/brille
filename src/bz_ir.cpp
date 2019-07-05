#include "bz_ir.h"


void IrreducibleBrillouinZone::set_vertices(ArrayVector<double> newverts){
  if (newverts.numel()!=3u) throw "IrreducibleBrillouinZone objects only take 3-D vectors for vertices";
  this->vertices = newverts;
}
void IrreducibleBrillouinZone::set_faces(ArrayVector<int> newfaces){
  if (newfaces.numel()!=3u) throw "IrreducibleBrillouinZone objects only take 3-D vectors for faces";
  this->faces = newfaces;
}
void IrreducibleBrillouinZone::set_faces_per_vertex(ArrayVector<int> newfpv){
  if (newfpv.numel()!=3u) throw "IrreducibleBrillouinZone vertices are the intersections of three faces each";
  this->faces_per_vertex = newfpv;
}

void IrreducibleBrillouinZone::wedge_search(const int time_reversal){
  // First we need the pointgroup (or at least its psuedo-rotation matrices)
  PointSymmetry psym = get_pointgroup_symmetry(this->lattice.hall,time_reversal);
  LQVec<double> normals(this->lattice, 2*psym.size());
  LQvec<double> uxy(this->lattice,4);
  // then we should find the highest (or lowest?) order rotation/rotoinversion

  // or just go through them in order for now
  std::array<std::array<int,3>,3> axis_plane_vecs;
  size_t nidx = 0;
  double a_over_b;
  for (size_t i=0; i<psym.size(); ++i){
    // get the rotation axis and two perpendicular vectors
    axis_plane_vecs = rotation_axis_and_perpendicular_vectors(psym.get(idx));
    // skip ̄1 and 1 'rotations' as they don't have a rotation axis:
    if (0==vector_norm_squared(axis_plane_vecs[0].data())) continue;
    // if the second or third vector is ⃗0 we can't do anything either
    if (0==vector_norm_squared(axis_plane_vecs[1].data())) continue;
    if (0==vector_norm_squared(axis_plane_vecs[2].data())) continue;
    // turn the hkl indices into fullfledged vectors, for dot and cross products
    for (size_t k=0; k<3; ++k)
      for (size_t j=0; j<3; ++j)
        uxy.insert(static_cast<double>(axis_plane_vecs[k][j]), k, j);
    // plus make a copy of the first perpendicular vector
    for (size_t j=0; j<3; ++j) uxy.insert(uxy.getvalue(1,j),3,j);
    if (nidx){
      // For the first rotation, stick with the first perpendicular vector,
      // otherwise we need to make sure we pick a vector which has a positive
      // dot product with all previous normal vectors.
      /*  Since we can always describe a vector as a linear comibnation of
          a set of basis vectors, every vector in the perpendicular plane can
          be expressed as ⃗v = a ⃗x + b ⃗y (as long as ⃗x and ⃗y span the plane).
          We want to pick (a,b) such that ⃗v ⋅ ̂nᵢ ≥ 0 for i=0,…,nidx-1.
          This is equivalent to (a ⃗x⋅n̂ᵢ + b ⃗y⋅n̂ᵢ) ≥ 0 → a/b ≥ -( ⃗x⋅n̂ᵢ )/( ⃗y⋅n̂ᵢ ).
          If we arbitrarily choose b≡1 then we can find the minimum a value
          required to satisfy each constraint, the maximum value of which will
          satisfy all constraints simuntaneously.
      */
      a_over_b = std::numeric_limits<double>::min();

    }
  }
  // its rotation axis and two perpendicular vectors
  // for the first we can blindly take the first perpendicular vector
  for (int i=0; i<3; ++i) x.insert(axis_plane_vecs[1][i],0,i);

}

void IrreducibleBrillouinZone::vertex_search(const int extent){
  // LQVec<int> tau(this->lattice);
  // int ntau = make_all_indices(&tau,extent);
  LQVec<int> tau(this->lattice, make_relative_neighbour_indices(extent) );
  int ntau = (int)tau.size();
  // the number of unique combinations of 3-taus is the number that we need to check
  size_t ntocheck=0;
  // there is probably a better way to do this, but brute force never hurt anyone
  for (int i=0; i<(ntau-2); i++) for (int j=i+1; j<(ntau-1); j++) for (int k=j+1; k<ntau; k++) ntocheck++;

  LQVec<double> all_vertices(this->lattice,ntocheck);
  ArrayVector<int> all_ijk(3,ntocheck);

  // LQVec<double> tauhat(this->lattice), halftau(this->lattice);
  // ArrayVector<double> lentau = tau.norm();
  ArrayVector<double> lentau = norm(tau);
  LQVec<double> tauhat = tau/lentau;
  double two = 2;
  LQVec<double> halftau(tau/two);

  ArrayVector<double> tauhat_xyz;
  tauhat_xyz = tauhat.get_xyz();

  int count=0;
  for (int i=0; i<(ntau-2); i++){
    for (int j=i+1; j<(ntau-1); j++){
      for (int k=j+1; k< ntau   ; k++){
        if ( three_plane_intersection(&tauhat, &halftau, &tauhat_xyz, i,j,k, &all_vertices, count) ){
          all_ijk.insert(i, count, 0); //insert value i at position (count,0)
          all_ijk.insert(j, count, 1);
          all_ijk.insert(k, count, 2);
          count++;
        }
      }
    }
  }
  // there are count intersections of three planes (strictly count<=ntocheck, but probably count < ntocheck/2)

  // next we need to check for the kernel of intersection points which are closer to the origin than any (non-intersection-defining) planes
  LQVec<double> in_verts(this->lattice,count);
  ArrayVector<int> in_ijk(3,count);
  int in_cnt = 0;
  for (int i=0; i<count; i++){
    // this between_origin_and_plane expects all vectors in an orthonormal frame
    if ( between_origin_and_plane( &halftau, &all_vertices, &all_ijk, i, &in_verts, in_cnt, 1e-10 ) ){
      in_ijk.set(in_cnt++, all_ijk.datapointer(i));
    }
  }
  // if ( in_cnt > 0) in_verts.print(0,in_cnt-1);

  // it's possible that multiple three-plane intersections have given the same
  // intersection point -- e.g., for a cubic system the intersection points of
  // (100),(010),(001); (110),(010),(001); (100),(110),(001); (100),(010),(011)
  // (101),(010),(001); ... are all the same point, (111).
  // The true vertex of the Brillouin Zone is the intersection of the three
  // planes with smallest norm

  // First, find the first vertex which is unique of equivalent vertices
  bool *vertisunique = new bool[in_cnt]();
  for (int i=0; i<in_cnt; i++) vertisunique[i] = true;
  for (int i=0; i<in_cnt-1; i++){
    if (vertisunique[i]){
      for (int j=i+1;j<in_cnt; j++){
        if (vertisunique[j] && in_verts.isapprox(i,j)){
          // printf("vert %d == %d\n",i,j);
           vertisunique[j] = false;
         }
      }
    }
  }
  // count up the unique vertices and keep track of their indices
  int unqcnt=0, *unqidx = new int[in_cnt]();
  for (int i=0; i<in_cnt; i++) if (vertisunique[i]) unqidx[unqcnt++]=i;

  //printf("%d unique vertices\n",unqcnt);

  // create a mapping which holds the indices of all equivalent vertices
  // so unqidxmap[1,:] are all of the indices which equal the first unique vertex
  int *unqidxmap = new int[in_cnt*unqcnt]();
  int *numidxmap = new int[unqcnt]();
  for (int i=0; i<unqcnt; i++){
    numidxmap[i]=1;
    unqidxmap[i*in_cnt] = unqidx[i];
    for (int j=0; j<in_cnt; j++){
      if (unqidx[i]!=j && !vertisunique[j] && in_verts.isapprox(unqidx[i],j))
        unqidxmap[ i*in_cnt + numidxmap[i]++ ] = j;
    }
  }
  delete[] vertisunique;
  // and determine the "length" of the planes which define each vertex, maintaining
  // the equivalence relationship already established, but allocating only the
  // memory actually needed by first finding the maximum number intersections
  // which gave the same vertex
  int maxequiv = 0;
  for (int i=0; i<unqcnt; i++) if ( numidxmap[i]>maxequiv) maxequiv=numidxmap[i];
  double *unqlenmap = new double[maxequiv*unqcnt](); // no need to allocate in_cnt*unqcnt memory when maxequiv is known
  for (int i=0; i<unqcnt; i++){
    for (int j=0; j<numidxmap[i]; j++){
      unqlenmap[i*maxequiv + j] = 0;
      for (int k=0; k<3; k++){
        unqlenmap[i*maxequiv + j] += halftau.norm( in_ijk.getvalue(unqidxmap[i*in_cnt+j],k));
      }
    }
  }
  // use the "length" information to select which equivalent vertex we should keep
  int *minequividx = new int[unqcnt]();
  double *minequivlen = new double[unqcnt]();
  for (int i=0; i<unqcnt; i++){
    minequivlen[i] = std::numeric_limits<double>::max(); // better than 1./0.
    for (int j=0; j<numidxmap[i]; j++){
      if ( unqlenmap[i*maxequiv +j] < minequivlen[i]){
        minequividx[i] = unqidxmap[i*in_cnt+j];
        minequivlen[i] = unqlenmap[i*maxequiv +j];
      }
    }
  }
  delete[] unqidx;
  delete[] numidxmap;
  delete[] unqidxmap;
  delete[] unqlenmap;
  delete[] minequivlen;

  if (unqcnt == 0)   throw std::runtime_error("No unique vertices found?!");

  LQVec<double> unq_vrt(this->lattice, unqcnt);
  ArrayVector<int> unq_ijk(3,unqcnt);
  for (int i=0; i<unqcnt; i++){
      unq_vrt.set( i, in_verts.datapointer(minequividx[i]) );
      unq_ijk.set( i, in_ijk.datapointer(minequividx[i]) );
  }
  delete[] minequividx;

  // store the reciprocal space positions of the vertices of the first Brillouin Zone
  this->set_vertices(unq_vrt); // does this work with LQVec smlst_vrt?

  // determine which of the taus actually contribute to at least one vertex
  int ncontrib=0, *contrib = new int[ntau]();
  for (int i=0; i<ntau; i++)
  for (int j=0; j<unqcnt; j++)
  if ( unq_ijk.getvalue(j,0) == i || unq_ijk.getvalue(j,1) == i || unq_ijk.getvalue(j,2) == i ){
    contrib[ncontrib++]=i;
    break;
  }
  ArrayVector<int> faces(3,ncontrib);
  for (int i=0; i<ncontrib; i++) faces.set(i, tau.datapointer(contrib[i]));
  this->set_faces(faces);

  // Each vertex is the intersection of three faces -- smlst_ijk contains the indexes into tau
  // but since tau contains planes which do not contribute to the first Brillouin Zone
  // we still have work to do. Replace the indices in smlst_ijk with their equivalent
  // indices into this->faces (using contrib as the map)
  ArrayVector<int> fpv(3,unqcnt);
  for (int i=0; i<unqcnt; i++)
  for (int j=0; j<3; j++)
  for (int k=0; k<ncontrib; k++){
    if ( unq_ijk.getvalue(i,j) == contrib[k] ){
      fpv.insert(k,i,j);
      break;
    }
  }
  this->set_faces_per_vertex(fpv);

  delete[] contrib;
}

LQVec<double> IrreducibleBrillouinZone::get_vertices(void) const {
  LQVec<double> lqverts(this->lattice,this->vertices);
  if (this->isprimitive())
    lqverts = transform_from_primitive(this->outerlattice,lqverts);
  return lqverts;
}
LQVec<double> IrreducibleBrillouinZone::get_primitive_vertices(void) const {
  return LQVec<double>(this->lattice, this->vertices);
}
LQVec<int>    IrreducibleBrillouinZone::get_faces   (void) const {
  LQVec<int> lqfaces(this->lattice,this->faces);
  if (this->isprimitive())
    lqfaces = transform_from_primitive(this->outerlattice,lqfaces);
  return lqfaces;
}
LQVec<int> IrreducibleBrillouinZone::get_primitive_faces(void) const {
  return LQVec<int>(this->lattice, this->faces   );
}
ArrayVector<int> IrreducibleBrillouinZone::get_faces_per_vertex(void) const {
  ArrayVector<int> out = this->faces_per_vertex; // make sure we return a copy, not the internal object
  return out;
}

void IrreducibleBrillouinZone::print() const {
  printf("IrreducibleBrillouinZone with %u vertices and %u faces\n",this->vertices_count(),this->faces_count());
}

template<typename T> ArrayVector<bool> IrreducibleBrillouinZone::isinside(const LQVec<T> &p){
  bool already_same = this->lattice.issame(p.get_lattice());
  bool transform_needed = this->outerlattice.issame(p.get_lattice());
  LQVec<T> pprim(this->lattice);
  if (!(already_same || transform_needed))
    throw std::runtime_error("Q points provided to IrreducibleBrillouinZone::isinside must be in the standard or primitive lattice used to define the IrreducibleBrillouinZone object");
  if (transform_needed)
    pprim = transform_to_primitive(this->outerlattice,p);
  const LQVec<T> & psl = transform_needed ? pprim : p;
  // this IrreducibleBrillouinZone object has already been instantiated, meaning that it
  // knows *which* reciprocal lattice points define it!
  // we just need to check whether the point(s) in p are closer to the origin than the planes which make-up the Brillouin Zone
  ArrayVector<bool> out(1u,psl.size());
  LQVec<double> fv(this->lattice, (this->faces)/2.0); // this->faces is a) ArrayVector<integer> and b) reciprocal lattice points, and we want their halves
  bool tmp = true;
  ArrayVector<double> fv2 = dot(fv,fv);
  ArrayVector<double> fvp(1u,fv.size());
  T vvj, vpj, tol = std::numeric_limits<T>::epsilon(); // zero for integer-type T
  for (size_t i=0; i<psl.size(); i++){
    fvp = dot( fv, psl.get(i) );
    tmp = true;
    for (int j=0; j<fv.size(); j++){
      vvj = fv2.getvalue(j);
      vpj = fvp.getvalue(j);
      // if fⱼ⋅(pᵢ-fⱼ) is larger than ϵ×fⱼ⋅(pᵢ+fⱼ) then pⱼ is outside of the IrreducibleBrillouinZone
      if ( (vpj-vvj) > (vpj+vvj)*tol && (vpj-vvj) > tol ){ // both in case vpj+vvj < 1
        // std::cout << std::to_string(vpj-vvj) << " > " << std::to_string((vpj+vvj)*tol);
        // std::cout << " and " << std::to_string(tol) << " ==> ";
        // std::cout <<  p.to_string(i) << " is outside" << std::endl;
        tmp = false;
        break;
      }
    }
    out.insert(tmp,i);
  }
  return out;
}

bool IrreducibleBrillouinZone::moveinto(const LQVec<double>& Q, LQVec<double>& q, LQVec<int>& tau){
  bool already_same = this->lattice.issame(Q.get_lattice());
  LQVec<double> Qprim(this->lattice), qprim(this->lattice);
  LQVec<int> tauprim(this->lattice);
  PrimitiveTransform PT(this->outerlattice.get_hall());
  bool transform_needed = ( PT.does_anything() && this->outerlattice.issame(Q.get_lattice()) );
  if (!(already_same || transform_needed))
    throw std::runtime_error("Q points provided to IrreducibleBrillouinZone::isinside must be in the standard or primitive lattice used to define the IrreducibleBrillouinZone object");

  if (transform_needed)  Qprim = transform_to_primitive(this->outerlattice,Q);
  const LQVec<double> & Qsl = transform_needed ? Qprim : Q;
  LQVec<double> & qsl = transform_needed ? qprim : q;
  LQVec<int> & tausl = transform_needed? tauprim : tau;

  // Determine which points in Q are already inside the first BZ
  ArrayVector<bool> allinside = this->isinside(Qsl);
  // ensure that qsl and tausl can hold each qi and taui
  qsl.resize(Qsl.size());
  tausl.resize(Qsl.size());

  LQVec<int> facehkl(this->lattice,this->faces);
  ArrayVector<double> facelen = norm(facehkl);
  LQVec<double> facenrm = facehkl/facelen;
  LQVec<double> qi;
  LQVec<int> taui;
  ArrayVector<double> q_dot_facenrm;
  ArrayVector<int> Nhkl;
  size_t maxat = 0;
  int maxnm = 0;
  size_t count =0;
  for (size_t i=0; i<Qsl.size(); i++){
    count = 0;
    qi = Qsl.get(i);
    taui = 0*tausl.get(i);
    while (!allinside.getvalue(i) && count++ < 50*facelen.size()){
      // std::cout << "Moving q = " << qi.to_string() << std::endl;
      q_dot_facenrm = dot( qi , facenrm );
      Nhkl = (q_dot_facenrm/facelen).round();
      // std::cout << "Nhkl = " << Nhkl.to_string() << std::endl;
      if ( Nhkl.areallzero() ) {allinside.insert(true,i); break;} // qi is *on* the Brilluoin Zone surface (or inside) so break.
      maxnm = 0;
      maxat = 0;
      for (size_t j=0; j<Nhkl.size(); ++j){
        if (Nhkl.getvalue(j)>=maxnm && (maxnm==0 || q_dot_facenrm.getvalue(j)>q_dot_facenrm.getvalue(maxat)) ){
          maxnm = Nhkl.getvalue(j);
          maxat = j;
        }
      }
      // std::cout << "Of which, the maximum is vector " << std::to_string(maxat);
      // std::cout << " with value " << facehkl.to_string(maxat) << " " << std::to_string(maxnm);
      // std::cout << " of which will be removed." << std::endl;

      qi -= facehkl[maxat] * (double)(maxnm); // ensure we subtract LQVec<double>
      taui += facehkl[maxat] * maxnm; // but add LQVec<int>

      allinside.insert(this->isinside(qi).getvalue(0), i);
    }
    qsl.set(i, &qi);
    tausl.set(i, &taui);
  }
  if (!allinside.arealltrue()){
    std::string msg;
    for (size_t i=0; i<Qsl.size(); ++i)
      if (!allinside.getvalue(i))
        msg += "Q=" + Qsl.to_string(i) + " is outside of the IrreducibleBrillouinZone "
            + " : tau = " + tausl.to_string(i) + " , q = " + qsl.to_string(i) + "\n";
    throw std::runtime_error(msg);
  }
  if (transform_needed){ // then we need to transform back q and tau
    q   = transform_from_primitive(this->outerlattice,qsl);
    tau = transform_from_primitive(this->outerlattice,tausl);
  }
  return allinside.arealltrue(); // return false if any points are still outside of the first Brilluoin Zone
}

bool three_plane_intersection(const LQVec<double> *n,                // plane normals
                              const LQVec<double> *p,                // a point on each plane
                              const ArrayVector<double> *xyz,        // the plane normals in a inverse Angstrom orthonormal coordinate system
                              const int i, const int j, const int k, // indices into n, p, xyz
                              LQVec<double> *iat,                    // output storage array of intersection points
                              const int idx)                         // the index where the intersection is inserted if found
                              {
  // we need to check whether the matrix formed by the orthonormal-frame components of the three planes is nearly-singular
  double *M = new double[9];
  xyz->get(i, M);
  xyz->get(j, M+3);
  xyz->get(k, M+6);
  double detM;
  detM = matrix_determinant(M);
  delete[] M;
  if ( my_abs(detM) > 1e-10 ){ // this 1e-10 provides a cutoff for far-from-origin intersections
    LQVec<double> ni,nj,nk, pi,pj,pk, cij,cjk,cki, tmp;
    ni=n->get(i);  nj=n->get(j);  nk=n->get(k);
    pi=p->get(i);  pj=p->get(j);  pk=p->get(k);
    // cij=ni.cross(&nj);  cjk=nj.cross(&nk);  cki=nk.cross(&ni);
    cij=cross(ni,nj);  cjk=cross(nj,nk);  cki=cross(nk,ni);

    // tmp = cjk*(pi.dot(&ni)) + cki*(pj.dot(&nj)) + cij*(pk.dot(&nk));
    tmp = cjk*dot(pi,ni) + cki*dot(pj,nj) + cij*dot(pk,nk);
    tmp /= detM;

    iat->set(idx, tmp.datapointer() );
    return true;
  }
  return false;
}

bool between_origin_and_plane(const LQVec<double> *p,
                              const LQVec<double> *v,
                              const ArrayVector<int> *ijk,
                              const int idx,
                              LQVec<double> *inv,
                              const int store_at,
                              const double tol){
  // p and v should be the points defining each plane and the vertices of the intersections of three planes
  ArrayVector<double> v_p = dot(v->get(idx),*p), p_p = dot(*p,*p);
  // we want to skip over the planes which gave us this intersection point
  size_t i, skip1, skip2, skip3;
  skip1 = (size_t) ijk->getvalue(idx,0);
  skip2 = (size_t) ijk->getvalue(idx,1);
  skip3 = (size_t) ijk->getvalue(idx,2);

  double vpi, ppi, eps = std::numeric_limits<double>::epsilon();
  for (i=0; i < p->size(); i++){
    if ( !(i==skip1||i==skip2||i==skip3) ){
      vpi = v_p.getvalue(i);
      ppi = p_p.getvalue(i);
      if ((vpi-ppi) > (vpi+ppi)*eps && (vpi-ppi) > eps && (vpi-ppi) > tol)
        return false;
    }
  }
  // none of p are closer to the origin than v(i)
  inv->set(store_at, v->datapointer(idx));
  return true;
}
