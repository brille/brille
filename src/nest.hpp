template<typename T>
void Nest<T>::construct(const Polyhedron& poly, const size_t max_branchings, const double max_volume){
  SimpleTet root_tet(poly);
  double exponent;
  exponent = std::log(root_tet.maximum_volume()/max_volume)/std::log(static_cast<double>(max_branchings));
  // we want more than one tetrahedra at the root node
  if (root_tet.number_of_tetrahedra()<2 && max_branchings > 0){
    root_tet = SimpleTet(poly, max_volume*std::pow(static_cast<double>(max_branchings)-1, exponent));
    exponent = std::log(root_tet.maximum_volume()/max_volume)/std::log(static_cast<double>(max_branchings));
  }
  // copy-over the vertices of the root node
  vertices_ = root_tet.get_vertices();
  // we need to make a guess about how many vertices we'll need:
  size_t number_density = static_cast<size_t>(poly.get_volume()/max_volume);
  size_t nVerts = vertices_.size();
  vertices_.resize(number_density + nVerts);
  // copy over the per-tetrahedron vertex indices to the root's branches
  const ArrayVector<size_t>& tvi{root_tet.get_vertices_per_tetrahedron()};
  for (size_t i=0; i<tvi.size(); ++i){
    std::array<size_t,4> single;
    for (size_t j=0; j<4u; ++j) single[j] = tvi.getvalue(i,j); // no need to adjust indices at this stage
    root_.branches().push_back(NestNode(single));
  }

  if (max_branchings > 0)
  for (auto b: root_.branches()) if (b.volume(vertices_) > max_volume)
    b = this->subdivide(b, 1u, max_branchings, max_volume, exponent, nVerts);

  // ensure that we only keep actual vertices:
  if (vertices_.size() > nVerts) vertices_.resize(nVerts);
}

template<typename T>
NestNode Nest<T>::subdivide(
  NestNode& node, const size_t nBr, const size_t maxBr,
  const double max_volume, const double exp, size_t& nVerts
){
  if (!node.is_terminal()) return node; // we can only branch un-branched nodes
  Polyhedron poly(vertices_.extract(node.boundary().vertices()));
  double mult = (maxBr > nBr) ? std::pow(static_cast<double>(maxBr-nBr),exp) : 1.0;
  SimpleTet node_tet(poly, max_volume*mult);
  // add any new vertices to the object's array, keep a mapping for all:
  std::vector<size_t> map;
  const ArrayVector<double>& ntv{node_tet.get_vertices()};
  for (size_t i=0; i<ntv.size(); ++i){
    const ArrayVector<double> vi{ntv.extract(i)};
    std::vector<size_t> idx = find(norm(vertices_.first(nVerts)-vi).is_approx("==",0.));
    if (idx.size()>1)
      throw std::runtime_error("Multiple matching vertices?!");
    if (idx.size()==1){
      map.push_back(idx[0]);
    } else {
      // new vertex. make sure we have room to store it:
      if (vertices_.size() < nVerts+1) vertices_.resize(nVerts*2);
      vertices_.set(nVerts, vi);
      map.push_back(nVerts++);
    }
  }
  // copy over the per-tetrahedron vertex indices, applying the mapping as we go
  const ArrayVector<size_t>& tvi{node_tet.get_vertices_per_tetrahedron()};
  for (size_t i=0; i<tvi.size(); ++i){
    std::array<size_t,4> single;
    for (size_t j=0; j<4u; ++j) single[j] = map[tvi.getvalue(i,j)];
    node.branches().push_back(NestNode(single));
  }
  // If there are still levels to go:
  if (nBr < maxBr)
  for (auto b: node.branches()) if (b.volume(vertices_) > max_volume)
    b = this->subdivide(b, nBr+1u, maxBr, max_volume, exp, nVerts);
  return node;
}
