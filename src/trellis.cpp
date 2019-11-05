#include "trellis.h"

bool CubeNode::indices_weights(const ArrayVector<double>& vertices,
                               const ArrayVector<double>& x,
                               std::vector<size_t>& indices,
                               std::vector<double>& weights) const
{
  // The CubeNode object contains the indices into `vertices` necessary to find
  // the 8 corners of the cube. Those indices should be ordered
  // (000) (100) (110) (010) (101) (001) (011) (111)
  // so that vertex_indices[i] and vertex_indices[7-i] are connected by a body diagonal
  ArrayVector<double> node_verts = vertices.extract(vertex_indices);
  double node_volume = abs(node_verts.extract(0)-node_verts.extract(7)).prod(1).getvalue(0,0);
  ArrayVector<double> w = abs(x - node_verts).prod(1)/node_volume; // the normalised volume of each sub-parallelpiped
  // If any normalised weights are greater than 1+eps() the point isn't in this node
  if (w.any_approx(">",1.)) return false;
  ArrayVector<bool> needed = w.is_approx(">", 0.);
  indices.clear();
  weights.clear();
  for (size_t i=0; i<8u; ++i) if (needed.getvalue(i)) {
    // the weight corresponds to the vertex opposite the one used to find the partial volume
    indices.push_back(vertex_indices[7-i]);
    weights.push_back(w.getvalue(i));
  }
  return true;
}
bool PolyNode::indices_weights(const ArrayVector<double>& vertices,
                               const ArrayVector<double>& x,
                               std::vector<size_t>& indices,
                               std::vector<double>& weights) const
{
  indices.clear();
  weights.clear();
  std::array<double,4> w{0,0,0,0};
  for (size_t i=0; i<vi_t.size(); ++i)
  if (this->tetrahedra_contains(i, vertices, x, w)){
    for (size_t j=0; j<4u; ++j) if (!approx_scalar(w[j],0.)){
      indices.push_back(vi_t[i][j]);
      weights.push_back(w[j]);
    }
    return true;
  }
  return false;
}
bool PolyNode::tetrahedra_contains(const size_t t,
                                   const ArrayVector<double>& v,
                                   const ArrayVector<double>& x,
                                   std::array<double,4>& w) const
{
  double vol6 = orient3d( v.data(vi_t[t][0u]), v.data(vi_t[t][1u]), v.data(vi_t[t][2u]), v.data(vi_t[t][3u]) );
  w[0] = orient3d( x.data(),            v.data(vi_t[t][1u]), v.data(vi_t[t][2u]), v.data(vi_t[t][3u]) )/vol6;
  w[1] = orient3d( v.data(vi_t[t][0u]), x.data(),            v.data(vi_t[t][2u]), v.data(vi_t[t][3u]) )/vol6;
  w[2] = orient3d( v.data(vi_t[t][0u]), v.data(vi_t[t][1u]), x.data(),            v.data(vi_t[t][3u]) )/vol6;
  w[3] = orient3d( v.data(vi_t[t][0u]), v.data(vi_t[t][1u]), v.data(vi_t[t][2u]), x.data()            )/vol6;
  if (std::any_of(w.begin(), w.end(), [](double z){return z < 0. && !approx_scalar(z, 0.);}))
    return false;
  return true;
}
