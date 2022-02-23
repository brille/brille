template<class T> template<class R> bool LQPolyhedron<T>::operator!=(const LQPolyhedron<R>& that) const {
  bool vertices_permuted{false};
  if (_vertices != that._vertices){
    vertices_permuted = _vertices.is_permutation(that._vertices);
    if (!vertices_permuted) return true;
  }
  const auto & that_faces{that._faces};
  if (_faces.size() != that_faces.size()) return true;

  auto add_to_no = [](const size_t & x, const std::vector<ind_t>& v){return x + v.size();};
  auto this_no = std::accumulate(_faces.begin(), _faces.end(), 0u, add_to_no);
  auto that_no = std::accumulate(that_faces.begin(), that_faces.end(), 0u, add_to_no);
  if (that_no != this_no) return true;

  std::vector<ind_t> permutation(_vertices.size(0));
  if (vertices_permuted){
    permutation = _vertices.permutation_vector(that._vertices);
  } else {
    std::iota(permutation.begin(), permutation.end(), 0u);
  }
  std::vector<std::vector<ind_t>> faces;
  faces.reserve(that_faces.size());
  auto permute = [permutation](const auto & i){return permutation[i];};
  for (const auto & that_face: that_faces){
    std::vector<ind_t> face;
    face.reserve(that_face.size());
    std::transform(that_face.begin(), that_face.end(), std::back_inserter(face), permute);
    faces.push_back(face);
  }
  for (const auto & face: _faces){
    auto at = faces.end();
    for (auto itr = faces.begin(); itr != at; ++itr){
      auto itr_size = (*itr).size();
      if (itr_size != face.size()) continue;
      for (size_t roll=0; roll < itr_size; ++roll){
        std::vector<ind_t> one(itr_size);
        for (size_t i=0; i < itr_size; ++i) one[i] = (*itr)[(i + roll) % itr_size];
        if (std::equal(face.begin(), face.end(), one.begin())){
          at = itr;
          break;
        }
      }
      if (at == itr) break;
    }
    if (at != faces.end()){
      faces.erase(at);  // make sure we don't match the same face twice
    } else {
      return true;
    }
  }
  return false;
}


template<class T> [[nodiscard]] LQVec<T> LQPolyhedron<T>::face_points() const {
  LQVec<T> p(_vertices.get_lattice(), face_count(), 3u);
  size_t idx{0};
  for (const auto & face: _faces){
    p.set(idx++, _vertices.extract(face).sum(0) / static_cast<T>(face.size()));
  }
  return p;
}

template<class T>[[nodiscard]] LQVec<T> LQPolyhedron<T>::face_normals() const {
  LQVec<T> n(_vertices.get_lattice(), face_count(), 3u);
  size_t idx{0};
  for (const auto & face: _faces) n.set(idx++, three_point_normal(_vertices, face));
  return n;
}

template<class T> [[nodiscard]] Array2<ind_t> LQPolyhedron<T>::edges() const {
  auto add_to_no = [](const size_t & x, const std::vector<ind_t>& v){return x + v.size();};
  auto no = std::accumulate(_faces.begin(), _faces.end(), 0u, add_to_no);
  std::vector<bool> unseen(no*no, true);
  Array2<ind_t> edges(no>>1, 2u);
  ind_t found{0}, a, b;
  for (const auto & face: _faces) for (size_t i=0; i<face.size(); ++i) {
    a = face[i];
    b = face[(i+1) % face.size()];
    if (unseen[a * no + b] && unseen[a + b * no]) {
      unseen[a * no + b] = unseen[a + b * no] = false;
      edges[{found, 0}] = a;
      edges[{found, 1}] = b;
      ++found;
    }
  }
  if (found != no >> 1) {
    std::string msg = "Found " + std::to_string(found) + " edge index pairs ";
    msg += "but expected to find " + std::to_string(no >> 1);
    if (found < no >> 1) msg += ". Is the LQPolyhedron open?";
    throw std::runtime_error(msg);
  }
  return edges;
}

template<class T> [[nodiscard]] LQVec<T> LQPolyhedron<T>::half_edges() const {
  auto indexes = this->edges();
  LQVec<T> half(_vertices.get_lattice(), indexes.size(0), 3u);
  for (ind_t i=0; i<indexes.size(0); ++i){
    half.set(i, (_vertices.view(indexes[{i, 0}]) + _vertices.view(indexes[{i, 1}])) / T(2));
  }
  return half;
}

template<class T> [[nodiscard]] T LQPolyhedron<T>::volume() const {
  /* per, e.g., http://wwwf.imperial.ac.uk/~rn/centroid.pdf
  For a polyhedron with N triangular faces, each with ordered vertices
  (aᵢ, bᵢ, cᵢ), one can define nᵢ = (bᵢ-aᵢ)×(cᵢ-aᵢ) for each face and then
  find that the volume of the polyhedron is V = 1/6 ∑ᵢ₌₁ᴺ aᵢ⋅ nᵢ
  */
  T volume{0}, sub{0};
  for (const auto & face: _faces){
    const auto a{_vertices.view(face[0])};
    for (ind_t i=1; i < face.size() - 1; ++i){
      volume += dot(cross(_vertices.view(face[i]) - a, _vertices.view(face[i+1]) - a), a).sum();
    }
  }
  return volume / T(6);
}

template<class T> [[nodiscard]] LQVec<T> LQPolyhedron<T>::centroid() const {
  LQVec<T> cen(_vertices.get_lattice(), {1u, 3u}, T(0));
  for (const auto & face: _faces) {
    const auto a{_vertices.view(face[0])};
    for (ind_t i=1; i < face.size() - 1; ++i) {
      const auto b{_vertices.view(face[i])};
      const auto c{_vertices.view(face[i + 1])};
      cen += cross(b-a, c-a) * ((a+b)*(a+b) + (b+c)*(b+c) + (c+a)*(c+a));
    }
  }
  return cen / (T(48) * this->volume());
}

template<class T>[[nodiscard]] T LQPolyhedron<T>::circumsphere_radius() const {
  auto c2v = _vertices - this->centroid();
  return norm(c2v).max(0).sum();
}

template<class T>[[nodiscard]] LQVec<T> LQPolyhedron<T>::rand_rejection(const ind_t n, const unsigned int seed) const {
  auto tics = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed > 0 ? seed : static_cast<unsigned int>(tics));
  std::uniform_real_distribution<T> distribution(T(0), T(1));
  auto min = _vertices.min(0);
  auto delta = _vertices.max(0) - min;
  LQVec<T> points(_vertices.get_lattice(), n, 3u);
  for (ind_t i=0; i < n; ){
    points.set(i, min + delta * distribution(generator));
    if (this->contains(points.view(i))[0]) ++i;
  }
  return points;
}

template<class T> template<class R> LQPolyhedron<T> LQPolyhedron<T>::operator+(const LQPolyhedron<R>& that) const {
  // combine vertices:
  auto vertices = cat(0, _vertices, that._vertices);
  // combine faces
  faces_t faces(_faces);
  faces.reserve(_faces.size() + that._faces.size());
  auto add_vertex_offset = [offset{_vertices.size(0)}](const ind_t i){return i + offset;};
  for (auto face: that._faces) {
    face_t one;
    one.reserve(face.size());
    std::transform(face.begin(), face.end(), std::back_inserter(one), add_vertex_offset);
    faces.push_back(one);
  }

  // check for duplicate vertices
  std::tie(vertices, faces) = remove_duplicate_points_and_update_face_indexing(vertices, faces);

  // check for overlapping faces
  /* such a check is hard, so we'll take a shortcut which is hopefully good enough:
   *    if any two faces contain three of the same vertices, they must be co-planar
   * we are, of course, discounting the possibility that two faces are coplanar and overlap without sharing points
   * but that seems *way* too difficult to handle at the moment. */
  for (ind_t i=0; i < faces.size()-1;  ++i) {
    const auto & a{faces[i]};
    for (ind_t j = i + 1; j < faces.size(); ++j) {
      const auto & b{faces[j]};
      auto c = std::count_if(a.begin(), a.end(), [&b](const auto & v){return std::count(b.begin(), b.end(), v) > 0;});
//      if (c == a.size()){
//        // they share all of their points, they're either duplicates or face each other
//        // in either case, we need to check for equal vertex ordering with rolling the starting index of one
//        // e.g., [0, 1, 2] == [0, 1, 2] || [1, 2, 3] || [2, 3, 0] || [3, 0, 1] for 4-point duplicate check
//        // or    [0, 2, 1] == [0, 1, 2] || [1, 2, 3] || [2, 3, 0] || [3, 0, 1] for 4-point opposites check
//      }
      if (c > 2) {
        std::string msg = "The faces " + std::to_string(i) + " and " + std::to_string(j);
        msg += " contain " + std::to_string(c) + " common vertex indexes, which indicates that they are coplanar!";
        throw std::runtime_error(msg);
        // If this gets thrown we should really think about handling this
      }
    }
  }
  return {vertices, faces};
}

// return a modified copy of this LQPolyhedron
template<class T> template<class R> LQPolyhedron<T> LQPolyhedron<T>::rotate(const std::array<R,9>& rot) const {
  T B[9], invB[9], invB_R[9];
  _vertices.get_lattice().get_B_matrix(B);
  utils::matrix_inverse(invB, B);
  utils::multiply_matrix_matrix(invB_R, invB, rot.data());
  auto xyz = _vertices.get_xyz(); // this is B * (hkl)
  for (ind_t i=0; i < xyz.size(0); ++i) {
    utils::mul_mat_vec_inplace(3u, invB_R, xyz.ptr(i)); // rotate and convert back to (hkl)
  }
  return {xyz, _faces};
}

template<class T> LQPolyhedron<T> LQPolyhedron<T>::apply(const PointSymmetry& ps, const ind_t index) const {
  // a point symmetry describes how the real space lattice transforms, but we have vertices in the reciprocal space
  // So we need to get the transposed rotation matrix
  auto Rt = transpose(ps.get(index));
  LQVec<T> rotated(_vertices.get_lattice(), _vertices.size(0), 3u);
  for (ind_t i=0; i < _vertices.size(0); ++i){
    utils::multiply_matrix_vector(rotated.ptr(i), Rt.data(), _vertices.ptr(i));
  }
  return {rotated, _faces};
}


// geometric properties in relation to another point or polyhedron
template<class T> template<class R> [[nodiscard]] std::vector<bool> LQPolyhedron<T>::contains(const LQVec<R>& x) const {
  std::vector<bool> out;
  out.reserve(x.size(0));
  auto n = this->normals();
  auto p = this->points();
  for (ind_t i=0; i<x.size(0); ++i){
    // FIXME, consider increasing the tolerance here!
    out.push_back(dot(n, x.view(i) - p).all(cmp::le, 0.));
  }
  return out;
}

template<class T> template<class R> [[nodiscard]] bool LQPolyhedron<T>::intersects(const LQPolyhedron<R>& that) const {
  return !approx::scalar(this->intersection(that).get_volume(), 0.);
}

template<class T> template<class R> [[nodiscard]] LQPolyhedron<T> LQPolyhedron<T>::intersection(const LQPolyhedron<R>& that) const {
  auto v = that.vertices();
  std::vector<std::vector<ind_t>> abc(3u);
  for (auto & i: abc) i.reserve(that.faces_count());
  for (const auto & face: that.faces()){
    for (ind_t i=0; i<3u; ++i) abc[i].push_back(face[i]);
  }
  return this->bisect(v.extract(abc[0]), v.extract(abc[1]), v.extract(abc[2]));
}

template<class T> template<class R> [[nodiscard]] LQPolyhedron<T> LQPolyhedron<T>::divide(const LQVec<R>& a, const LQVec<R>& b, const LQVec<R>& c) const {
  // this all seems unnecessary
  auto z = this->centroid();
  LQPolyhedron<T> centred(_vertices - z, _faces);
  auto divided = centred.bisect(a - z, b - z, c - z);
  return divided.translate(z);
}

template<class T> template<class R> [[nodiscard]] size_t LQPolyhedron<T>::face_index(const LQVec<R>& a, const LQVec<R>& b, const LQVec<R>& c) const {
  auto match{_faces.size()};
  auto ahkl = a.get_hkl();
  auto bhkl = b.get_hkl();
  auto chkl = c.get_hkl();
  auto vhkl = _vertices.get_hkl();
  for (size_t i=0; i< _faces.size(); ++i) {
    const auto & face{_faces[i]};
    std::vector<double> o3d;
    o3d.reserve(face.size());
    for (const auto & vertex: face){
      o3d.push_back(orient3d(ahkl.ptr(0), bhkl.ptr(0), chkl.ptr(0), vhkl.ptr(vertex)));
    }
    auto all_zero = std::all_of(o3d.begin(), o3d.end(), [](const auto & x){return approx::scalar(x, 0.);});
    if (all_zero) {
      all_zero = dot(three_point_normal(_vertices, face[0], face[1], face[2]), three_point_normal(a, b, c)).all(cmp::gt, 0.);
    }
    if (all_zero) match = i;
  }
  return match;
}

template<class T> template<class R> [[nodiscard]] bool LQPolyhedron<T>::none_beyond(const LQVec<R>& a, const LQVec<R>& b, const LQVec<R>& c) const {
  auto match = face_index(a, b, c);
  face_t v(_vertices.size(0));
  std::iota(v.begin(), v.end(), 0);
  if (match < _faces.size()){
    auto on_face = [f=_faces[match]](const auto & x){return std::find(f.begin(), f.end(), x) != f.end();};
    auto itr = std::remove_if(v.begin(), v.end(), on_face);
    v.erase(itr, v.end());
  }
  auto ahkl = a.get_hkl();
  auto bhkl = b.get_hkl();
  auto chkl = c.get_hkl();
  auto vhkl = _vertices.get_hkl();
  auto is_negative = [&](const auto & i){
    return orientd3(ahkl.ptr(0), bhkl.ptr(0), chkl.ptr(0), vhkl.ptr(i)) < 0;
  };
  return std::all_of(v.begin(), v.end(), is_negative);
}

template<class T> template<class R> LQPolyhedron<T> LQPolyhedron<T>::cut(const LQVec<R>& a, const LQVec<R>& b, const LQVec<R>& c) const {
  assert(a.size(0) == b.size(0) && a.size(0) == c.size(0) && a.size(0) == 1);
  if (this->none_beyond(a, b, c)) return *this;
  auto keep = point_inside_plane(a, b, c, _vertices);
  if (std::find(keep.begin(), keep.end(), false) == keep.end()) return *this;
  auto cut = keep_to_cut_list(keep, _faces);
  face_t new_face;
  ind_t face_count{0}, vertex_count{0};
  auto v = this->vertices(); // make a copy,  I hope
  auto f = this->faces(); // ditto
  for (size_t j=0; j < cut.size() - 1; ++j) for (size_t k=j+1; k < cut.size(); ++k){
    auto at = edge_plane_intersection(v, f[cut[j]], f[cut[k]], a, b, c);
    if (at.size(0) == 1){
      auto index = v.first(cmp::eq, at); // == v.size(0) if not found
      if (index < v.size(0)) {
        if (index >= v.size(0) - vertex_count){
          info_update("We have matched a vertex added by the cut; sort-out averaging the repeatedly found vertex");
        }
        keep[index] = true; // insurance against removing a needed vertex
      } else {
        v.append(0, at);
        ++vertex_count;
      }
      utils::add_if_missing(f[cut[j]], index);
      utils::add_if_missing(f[cut[k]], index);
      new_face.push_back(index);
      ++face_count;
    }
  }
  // add the new face to the list of faces, if it is not already present
  if (!new_face.empty() && !utils::unordered_list_in_lists(new_face, f)){
    f.push_back(sort_convex_polygon_face(a, b, c, new_face, v));
  }
  // extend to cover the newly found vertices
  for (ind_t z=0; z<vertex_count; ++z) keep.push_back(true);

  std::tie(v, f) = remove_points_and_update_face_indexing(keep, v, f);

  // remove any faces without three vertices
  auto itr = std::remove_if(f.begin(), f.end(), [](const auto & face){return unique(face).size() < 3;});
  f.erase(itr, f.end());

  // remove dangling faces
  bool again{false};
  do {
    face_t adjacent_face_count(v.size(0), 0u);
    for (ind_t j = 0; j < v.size(0); ++j) {
      for (const auto &face: f) {
        if (std::find(face.begin(), face.end(), j) != face.end()) ++adjacent_face_count[j];
      }
    }
    auto check = [c = adjacent_face_count](const auto &face) { return !is_not_dangling(c, face); };
    std::vector<bool> not_ok;
    std::transform(f.begin(), f.end(), std::back_inserter(not_ok), check);
    again = std::find(not_ok.begin(), not_ok.end(), true) != not_ok.end();
    if (again) {
      f.erase(std::remove_if(f.begin(), f.end(), check), f.end());
    }
  } while (again);

  // remove any vertices not on a face
  std::vector<bool> kv(v.size(0), false);
  for (const auto & face: f) for (const auto & x: face) kv[x] = true;
  if (std::find(kv.begin(), kv.end(), false) != kv.end()){
    face_t kv_map;
    ind_t kv_count{0};
    for (auto && j : kv) kv_map.push_back(j ? kv_count++ : v.size(0));
    for (auto & face: f){
      // a direct transform is OK since no face vertices are beyond the valid mapping
      std::transform(face.begin(), face.end(), face.begin(), [kv_map](const auto & i){return kv_map[i];});
    }
    v = v.extract(kv);
  }

  if (v.size(0) < 4 || f.size() < 4) return LQPolyhedron<T>();
  return {v, f};
}

template<class T> template<class R> LQPolyhedron<T> LQPolyhedron<T>::bisect(const LQVec<R>& p_a, const LQVec<R>& p_b, const LQVec<R>&  p_c) const {
  // this is, by far, the worst possible implementation of this algorithm
  auto plane_count = p_a.size(0);
  assert(p_a.ndim()==2 && p_b.ndim()==2 && p_c.ndim() == 2 && p_a.size(1)==3 && p_b.size(1)==(3) && p_c.size(1) == 3);
  assert(p_a.size(0) == p_b.size(0) && p_b.size(0) == p_c.size(0));
  LQPolyhedron<T> out(*this);
  std::vector<ind_t> vertex_map;

  auto pv = out.vertices();
  auto pn = out.normals();
  auto faces = out.faces();

  for (ind_t plane=0; plane < plane_count; ++plane){
    out = out.cut(p_a.view(plane), p_b.view(plane), p_c.view(plane));
    if (approx::scalar(out.volume(), 0.)) break;
  }
  return out;
}
