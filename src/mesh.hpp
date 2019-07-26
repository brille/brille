
template<class T> void Mesh3<T>::check_elements(void){
  size_t total_elements = 1u;
  // scalar + eigenvector + vector + matrix*matrix elements
  size_t known_elements = static_cast<size_t>(this->elements[0])
                        + static_cast<size_t>(this->elements[1])
                        + static_cast<size_t>(this->elements[2])
                        + static_cast<size_t>(this->elements[3])*static_cast<size_t>(this->elements[3]);
  // no matter what, shape[0] should be the number of gridded points
  if (shape.size()>2){
    // if the number of dimensions of the shape array is greater than two,
    // the second element is the number of modes per point                    */
    this->branches = shape.getvalue(1u);
    for (size_t i=2u; i<this->shape.size(); ++i) total_elements *= shape.getvalue(i);
  } else {
    // shape is [n_points, n_elements] or [n_points,], so there is only one mode
    this->branches = 1u;
    total_elements = shape.size() > 1 ? shape.getvalue(1u) : 1u;
  }
  if (0 == known_elements)
    this->elements[0] = total_elements;
  if (known_elements && known_elements != total_elements){
    std::string msg ="Inconsistent element counts: "
                    + std::to_string(known_elements) + " = "
                    + std::to_string(this->elements[0]) + "+"
                    + std::to_string(this->elements[1]) + "+"
                    + std::to_string(this->elements[2]) + "+"
                    + std::to_string(this->elements[3]) + "² ≠ "
                    + std::to_string(total_elements);
    throw std::runtime_error(msg);
  }
}
template<class T> int Mesh3<T>::replace_data(const ArrayVector<T>& newdata,
                                             const ArrayVector<size_t>& newshape,
                                             const std::array<unsigned,4>& newelements)
{
  this->data = newdata;
  this->shape = newshape;
  this->elements = newelements;
  this->check_elements();
  return 0;
}
template<class T> int Mesh3<T>::replace_data(const ArrayVector<T>& newdata,
                                             const std::array<unsigned,4>& newelements)
{
  ArrayVector<size_t> shape(1,2);
  shape.insert(newdata.size(),0);
  shape.insert(newdata.numel(),1);
  return this->replace_data(newdata, shape, newelements);
}


// Perform sanity checks before attempting to interpolate
template<typename T> template<typename R> unsigned int Mesh3<T>::check_before_interpolating(const ArrayVector<R>& x) const{
  unsigned int mask = 0u;
  if (this->data.size()==0)
    throw std::runtime_error("The mesh must be filled before interpolating!");
  if (x.numel()!=3u)
    throw std::runtime_error("Mesh3 requires x values which are three-vectors.");
  return mask;
};
//! Perform linear interpolating at the specified points in the mesh's orthonormal frame
template<typename T> template<typename R> ArrayVector<T> Mesh3<T>::interpolate_at(const ArrayVector<R>& x) const{
  this->check_before_interpolating(x);
  ArrayVector<T> out(this->data.numel(), x.size());

  Delaunay::Point point;
  Delaunay::Cell_handle cell;
  // Delaunay::Vertex vertex;
  Delaunay::Locate_type type;
  int v0, v1;
  size_t corners[4];
  ArrayVector<double> verts(3u, 4u);
  std::vector<double> weights;
  for (size_t i=0; i<x.size(); ++i){
    point = Delaunay::Point(x.getvalue(i,0), x.getvalue(i,1), x.getvalue(i,2));
    cell = this->mesh.locate(point, type, v0, v1);
    switch (type){
      case Delaunay::VERTEX: // v0 is the index into cell of the exactly-matching vertex
        {
          verts.resize(1u);
          auto vertex = cell->vertex(v0);
          corners[0] = vertex->info();
          point = vertex->point();
          for (size_t j=0; j<3; ++j) verts.insert(point[j], 0, j);
        }
        break;
      case Delaunay::EDGE: // v0 and v1 index into cell of the vertices forming the line the point is on
        {
          verts.resize(2u);
          auto vertex = cell->vertex(v0);
          corners[0] = vertex->info();
          point = vertex->point();
          for (size_t j=0; j<3; ++j) verts.insert(point[j], 0, j);
          vertex = cell->vertex(v1);
          corners[1] = vertex->info();
          point = vertex->point();
          for (size_t j=0; j<3; ++j) verts.insert(point[j], 1, j);
        }
        break;
      case Delaunay::FACET: // v0 is the vertex opposite the triangular face in which the point lies (v1 is meaningless)
        {
          verts.resize(3u);
          v1=0; // abuse v1 as a counter
          for (int k=0; k<4; ++k) if (k != v0) {
            auto vertex = cell->vertex(k);
            corners[v1] = vertex->info();
            for (size_t j=0; j<3; ++j) verts.insert(point[j], v1, j);
            ++v1;
          }
        }
        break;
      case Delaunay::CELL: // the point is within the volume of the cell (v0 and v1 are meaningless)
        {
          verts.resize(4u);
          for (int k=0; k<4; ++k){
            auto vertex = cell->vertex(k);
            corners[k] = vertex->info();
            for (size_t j=0; j<3; ++j) verts.insert(point[j], k, j);
          }
        }
        break;
      case Delaunay::OUTSIDE_CONVEX_HULL:
      default:
        throw std::runtime_error("Mesh3<T>::interpolate_at: This should not be possible");
    }
    weights = tetrahedron_weights(x.extract(i), verts);
    new_unsafe_interpolate_to(this->data, this->elements, this->branches, verts.size(), corners, weights, out, i);
  }
  return out;
}


template<class T>
template<class R, class S>
ArrayVector<S> Mesh3<T>::debye_waller_sum(const LQVec<R>& Q, const R t_K) const{
  return this->debye_waller_sum(Q.get_xyz(), t_K);
}

template<class T>
template<class R, class S>
ArrayVector<S> Mesh3<T>::debye_waller_sum(const ArrayVector<R>& Q, const R t_K) const{
  const S hbar = 6.582119569E-13; // meV⋅s
  const S kB   = 8.617333252E-2; // meV⋅K⁻¹
  if (Q.numel() != 3)
    throw std::runtime_error("Debye-Waller factor requires 3-vector Q.");
  if (this->elements[0] != 1u)
    throw std::runtime_error("Debye-Waller factor requires one scalar (energy) per mode.");
  size_t nIons = this->elements[1] / 3u;
  if (0 == nIons || this->elements[2]*3u != nIons)
    throw std::runtime_error("Debye-Waller factor requires 3-vector eigenvector(s).");
  size_t nQ = Q.size();
  ArrayVector<S> WdQ(nIons,nQ); // Wᵈ(Q) has nIons entries per Q point

  S coth_en, Q_dot_e_2;
  size_t span = 1u + nIons*3u + this->elements[2] + this->elements[3]*this->elements[3];
  size_t nq = this->shape.getvalue(0u);

  const S beta = kB*t_K; // meV
  const S pref{hbar*hbar/static_cast<S>(2*nq)}; // meV²⋅s²

  S qj_sum;
  // for each input Q point
  for (size_t Qidx=0; Qidx<nQ; ++Qidx){
    // and each ion
    for (size_t d=0; d<nIons; ++d){
      qj_sum = S(0);
      // sum over all reduced q in the first Brillouin zone
      for (size_t q=0; q<nq; ++q){
        // and over all 3*nIon branches at each q
        for (size_t j=0; j<this->branches; ++j){
          // for each branch energy, find <2nₛ+1>/ħωₛ ≡ coth(2ħωₛβ)/ħωₛ
          coth_en = coth_over_en(this->data.getvalue(q,j*span), beta);
          // and find |Q⋅ϵₛ|². Note: vector_product(x,y) *is* |x⋅y|²
          Q_dot_e_2 = vector_product(3u, Q.datapointer(Qidx), this->data.datapointer(q,j*span+1u+3u*d));
          // adding |Q⋅ϵₛ|²coth(2ħωₛβ)/ħωₛ to the sum over s for [Qidx, d]
          qj_sum += Q_dot_e_2 * coth_en;
        }
      }
      // with the sum over s complete, normalize by ħ²/2 divided by the number
      // of points in the Brillouin zone and store the result at W[Qidx, d];
      WdQ.insert(qj_sum*pref, Qidx, d);
    }
  }
  return WdQ;
}

template<class T>
template<class R, template<class> class A, class S>
ArrayVector<S> Mesh3<T>::debye_waller(const A<R>& Q, const std::vector<R>& M, const R t_K) const{
  size_t nIons = this->elements[1] / 3u;
  if (0 == nIons || this->elements[1]*3u != nIons)
    throw std::runtime_error("Debye-Waller factor requires 3-vector eigenvector(s).");
  if (M.size() != nIons)
    throw std::runtime_error("Debye-Waller factor requires an equal number of ions and masses.");
  ArrayVector<S> WdQ = this->debye_waller_sum(Q, t_K);
  ArrayVector<S> factor(1u, Q.size());
  S d_sum;
  for (size_t Qidx=0; Qidx<Q.size(); ++Qidx){
    d_sum = S(0);
    for (size_t d=0; d<nIons; ++d){
      d_sum += std::exp(WdQ.getvalue(Qidx, d)/M[d]);
    }
    factor.insert(d_sum*d_sum, Qidx);
  }
  return factor;
}
