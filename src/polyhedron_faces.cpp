#include <set>
#include <map>
#include "polyhedron_faces.hpp"

using namespace brille;
using namespace brille::polyhedron;

bool Faces::operator!=(const Faces& that) const {
  const auto & that_faces{that._faces};
  if (_faces.size() != that_faces.size()) return true;
  auto add_to_no = [](const size_t & x, const std::vector<ind_t>& v){return x + v.size();};
  auto this_no = std::accumulate(_faces.begin(), _faces.end(), 0u, add_to_no);
  auto that_no = std::accumulate(that_faces.begin(), that_faces.end(), 0u, add_to_no);
  if (that_no != this_no) return true;

  std::vector<std::vector<ind_t>> faces;
  faces.reserve(that_faces.size());
  for (const auto & that_face: that_faces){
    std::vector<ind_t> face;
    face.reserve(that_face.size());
    std::copy(that_face.begin(), that_face.end(), std::back_inserter(face));
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

[[nodiscard]] std::vector<std::vector<std::pair<ind_t, ind_t>>> Faces::edges_per_face() const {
  std::vector<std::vector<std::pair<ind_t, ind_t>>> epf;
  for (const auto & face: _faces) {
    std::vector<std::pair<ind_t, ind_t>> ep;
    for (size_t i=0; i<face.size(); ++i) {
      ep.emplace_back(face[i], face[(i+1) % face.size()]);
    }
    epf.push_back(ep);
  }
  return epf;
}

[[nodiscard]] Array2<ind_t> Faces::edges() const {
  // The number of faces plus the number of vertices *is* the number of edges plus two:
  std::map<ind_t, ind_t> known;
  ind_t n_vert{0};
  for (const auto & face: _faces) for (const auto & v: face) if(known.count(v) == 0) known[v] = n_vert++;
  assert(known.size() == n_vert);
  size_t edge_count = _faces.size() + known.size() - 2u;
  Array2<ind_t> edges(static_cast<ind_t>(edge_count), 2u);
//  auto add_to_no = [](const size_t & x, const std::vector<ind_t>& v){return x + v.size();};
//  auto no = std::accumulate(_faces.begin(), _faces.end(), 0u, add_to_no);
//  std::vector<bool> unseen(no*no, true);
//  Array2<ind_t> edges(no>>1, 2u);
  std::vector<bool> unseen(n_vert * n_vert, true);
  ind_t found{0};
  for (const auto & face: _faces) for (size_t i=0; i<face.size(); ++i) {
      auto a = known[face[i]];
      auto b = known[face[(i+1) % face.size()]];
      if (unseen[a * n_vert + b] && unseen[a + b * n_vert]) {
        unseen[a * n_vert + b] = unseen[a + b * n_vert] = false;
        edges[{found, 0}] = face[i];
        edges[{found, 1}] = face[(i+1) % face.size()];
        ++found;
      }
    }
  if (found != edge_count) {
    std::string msg = "Found " + std::to_string(found) + " edge index pairs ";
    msg += "but expected to find " + std::to_string(edge_count);
    if (found < edge_count) msg += ". Is the LQPolyhedron open?";
    throw std::runtime_error(msg);
  }
  return edges;
}

[[nodiscard]] Array2<ind_t> Faces::planes() const {
  auto fs = static_cast<ind_t>(_faces.size());
  Array2<ind_t> planes(fs, 3u);
  for (ind_t i=0; i<fs; ++i){
    for (ind_t j=0; j<3; ++j) {
      planes[{i, j}] = _faces[i][j];
    }
  }
  return planes;
}


Faces Faces::combine(const Faces& that, const ind_t offset) const {
  faces_t faces(_faces);
  faces.reserve(size() + that.size());
  auto add_offset = [offset](const auto i){return i + offset;};
  for (const auto & face: that._faces) {
    face_t one;
    one.reserve(face.size());
    std::transform(face.begin(), face.end(), std::back_inserter(one), add_offset);
    faces.push_back(one);
  }
  return Faces(faces);
}