#ifndef _TRIANGULATION_H_
#define _TRIANGULATION_H_
/*----------------------------- HEADER INCLUDES ------------------------------*/
// combined:
#include <vector>
#include <array>
#include <cassert>
// Tetgen header
#include "tetgen/tetgen.h"
// debugging output
#include "debug.h"

// template<class T, template<class> class L, typename=typename std::enable_if<std::is_base_of<ArrayVector<T>,L<T>>::value>::type>
class TetrahedralTriangulation{
  // L<T> vertex_positions; // (nVertices, 3), ArrayVector<T>, LQVec<T>, or LDVec<T>
  size_t nVertices;
  size_t nTetrahedra;
  ArrayVector<double> vertex_positions; // (nVertices, 3)
  ArrayVector<size_t> vertices_per_tetrahedron; // (nTetrahedra, 4)
  std::vector<std::vector<size_t>> tetrahedra_per_vertex; // (nVertices,)(1+,)
  std::vector<std::vector<size_t>> neighbours_per_tetrahedron; // (nTetrahedra,)(1+,)
public:
  size_t number_of_vertices(void) const {return nVertices;};
  size_t number_of_tetrahedra(void) const {return nTetrahedra;};
  const ArrayVector<double>& get_vertex_positions(void) const {return vertex_positions;};
  const ArrayVector<size_t>& get_vertices_per_tetrahedron(void) const {return vertices_per_tetrahedron;};

  TetrahedralTriangulation(void): nVertices(0), nTetrahedra(0), vertex_positions({3u,0u}), vertices_per_tetrahedron({4u,0u}){};
  TetrahedralTriangulation(const tetgenio& tgio): vertex_positions({3u,0u}), vertices_per_tetrahedron({4u,0u}){
    nVertices = static_cast<size_t>(tgio.numberofpoints);
    nTetrahedra = static_cast<size_t>(tgio.numberoftetrahedra);
    // copy-over all vertex positions:
    vertex_positions.resize(nVertices);
    for (size_t i=0; i<nVertices; ++i)
    for (size_t j=0; j<3u; ++j)
    vertex_positions.insert(tgio.pointlist[3*i+j], i,j);
    // copy-over all tetrahedron vertex indices
    vertices_per_tetrahedron.resize(nTetrahedra);
    for (size_t i=0; i<nTetrahedra; ++i)
    for (size_t j=0; j<4u; ++j)
    vertices_per_tetrahedron.insert(static_cast<size_t>(tgio.tetrahedronlist[i*tgio.numberofcorners+j]),i,j);
    // Construct the tetrahedra per vertex vector of vectors
    tetrahedra_per_vertex.resize(nVertices);
    for (size_t i=0; i<nVertices; ++i)
    for (size_t j=0; j<nTetrahedra; ++j)
    for (size_t k=0; k<4u; ++k)
    if (vertices_per_tetrahedron.getvalue(j,k)==i)
      tetrahedra_per_vertex[i].push_back(j);
    // Construct the neighbours per tetrahedron vector of vectors
    neighbours_per_tetrahedron.resize(nTetrahedra);
    for (size_t i=0; i<nTetrahedra; ++i)
    for (size_t j=0; j<4u; ++j)
    if (tgio.neighborlist[i*4u+j] >= 0)
      neighbours_per_tetrahedron[i].push_back(static_cast<size_t>(tgio.neighborlist[i*4u+j]));
  }
  // emulate the locate functionality of CGAL -- for which we need an enumeration
  enum Locate_Type{ VERTEX, EDGE, FACET, CELL, OUTSIDE_CONVEX_HULL };
  size_t locate(const ArrayVector<double>& x, Locate_Type& type, size_t& v0, size_t& v1) const{
    if (x.numel() != 3u || x.size() != 1u){
      std::string msg = "locate requires a single 3-element vector.";
      throw std::runtime_error(msg);
    }
    // find the closest vertex to our test-point:
    ArrayVector<double> d = norm(vertex_positions - x); // or does this need to be (vertex_position-x).norm()?
    double min = std::numeric_limits<double>::max();
    size_t idx = nVertices + 1;
    for (size_t i=0; i<d.size(); ++i) if (d.getvalue(i)<min){
      min = d.getvalue(i);
      idx = i;
    }
    // if the minimum distance is zero we're done:
    if (approx_scalar(min, 0.0)){
      type = VERTEX;
      return tetrahedra_per_vertex[idx][0];
    }
    // Check each of the tetrahedra which contain this vertex to see if the
    // point is inside/on the surface/on an edge of the tetrahedron
    double o[5];
    size_t zero_count;
    for (auto tet_idx: tetrahedra_per_vertex[idx]){
      // ~6× the volume of this tetrahedron -- really here in case the tetrahedron vertices are not ordered correctly
      o[4] = orient3d(vertex_positions.datapointer(vertices_per_tetrahedron.getvalue(tet_idx,0u)),
                      vertex_positions.datapointer(vertices_per_tetrahedron.getvalue(tet_idx,1u)),
                      vertex_positions.datapointer(vertices_per_tetrahedron.getvalue(tet_idx,2u)),
                      vertex_positions.datapointer(vertices_per_tetrahedron.getvalue(tet_idx,3u)));
      // replace the first vertex by x
      o[0] = orient3d(x.datapointer(),
                      vertex_positions.datapointer(vertices_per_tetrahedron.getvalue(tet_idx,1u)),
                      vertex_positions.datapointer(vertices_per_tetrahedron.getvalue(tet_idx,2u)),
                      vertex_positions.datapointer(vertices_per_tetrahedron.getvalue(tet_idx,3u)));
      // replace the second vertex by x
      o[1] = orient3d(vertex_positions.datapointer(vertices_per_tetrahedron.getvalue(tet_idx,0u)),
                      x.datapointer(),
                      vertex_positions.datapointer(vertices_per_tetrahedron.getvalue(tet_idx,2u)),
                      vertex_positions.datapointer(vertices_per_tetrahedron.getvalue(tet_idx,3u)));
      // replace the third vertex by x
      o[2] = orient3d(vertex_positions.datapointer(vertices_per_tetrahedron.getvalue(tet_idx,0u)),
                      vertex_positions.datapointer(vertices_per_tetrahedron.getvalue(tet_idx,1u)),
                      x.datapointer(),
                      vertex_positions.datapointer(vertices_per_tetrahedron.getvalue(tet_idx,3u)));
      // replace the fourth vertex by x
      o[3] = orient3d(vertex_positions.datapointer(vertices_per_tetrahedron.getvalue(tet_idx,0u)),
                      vertex_positions.datapointer(vertices_per_tetrahedron.getvalue(tet_idx,1u)),
                      vertex_positions.datapointer(vertices_per_tetrahedron.getvalue(tet_idx,2u)),
                      x.datapointer());
      if ((o[0]/o[4])<0 || (o[1]/o[4])<0 || (o[2]/o[4])<0 || (o[3]/o[4])<0)
        continue; // the point is outside of this tetrahedron
      zero_count = 0;
      for (size_t i=1; i<5u; ++i) if (o[i]==0) ++zero_count; // use approx_scalar instead?
      switch (zero_count){
        case 0:
          type = CELL;
          break;
        case 1:
          type = FACET;
          v0 = (o[0]==0) ? 0 : (o[1]==0) ? 1: (o[2]==0) ? 2 : 3;
          break;
        case 2:
          type = EDGE;
          v0 = (o[0]!=0) ? 0 : (o[0]!=0) ? 1 : 2;
          v1 = (o[v1+1]!=0) ? v1+1 : (o[v1+2]!=0) ? v1+2 : v1+3;
          break;
        case 3:
          type = VERTEX;
          v0 = (o[0]!=0) ? 0 : (o[1]!=0) ? 1 : (o[2]!=0) ? 2 : 3;
          break;
        default:
          throw std::logic_error("Intentionally unreachable switch statement.");
      }
      return tet_idx;
    }
    type = OUTSIDE_CONVEX_HULL;
    return 0;
  };
  // Make a new locator which slots into the interpolation routine more easily
  std::vector<size_t> locate_for_interpolation(const ArrayVector<double>& x) const {
    if (x.numel()!=3u && x.size()!=1u){
      std::string msg = "locate requires a single 3-element vector.";
      throw std::runtime_error(msg);
    }
    size_t v0, v1;
    Locate_Type type;
    size_t tet_idx = this->locate(x, type, v0, v1);
    std::vector<size_t> vert_idx;
    switch (type){
      case VERTEX:
        vert_idx.push_back(vertices_per_tetrahedron.getvalue(tet_idx,v0));
        break;
      case EDGE:
        vert_idx.push_back(vertices_per_tetrahedron.getvalue(tet_idx,v0));
        vert_idx.push_back(vertices_per_tetrahedron.getvalue(tet_idx,v1));
        break;
      case FACET:
        for (size_t i=0; i<4u; ++i) if (i!=v0)
        vert_idx.push_back(vertices_per_tetrahedron.getvalue(tet_idx, i));
        break;
      case CELL:
        for (size_t i=0; i<4u; ++i)
        vert_idx.push_back(vertices_per_tetrahedron.getvalue(tet_idx, i));
        break;
      case OUTSIDE_CONVEX_HULL:
        throw std::logic_error("Interpolation attempted for a point outside the convex hull");
        break;
      default:
        throw std::logic_error("Intentionally unreachable switch statement..");
    }
    return vert_idx;
  };
  /* If the following function is more useful than the preceeding, it could be
     advantageous to replicate the above code in this function's for loop.    */
  std::vector<std::vector<size_t>> locate_all_for_interpolation(const ArrayVector<double>& x) const {
    if (x.numel()!=3u){
      std::string msg = "locate_all requires 3-element vector(s)";
      throw std::runtime_error(msg);
    }
    std::vector<std::vector<size_t>> all_idx(x.size());
    for (size_t i=0; i<x.size(); ++i)
      all_idx[i] = this->locate_for_interpolation(x.extract(i));
    return all_idx;
  }
};

template <typename T>
TetrahedralTriangulation triangulate(const ArrayVector<T>& verts,
                                     const std::vector<std::vector<int>>& vpf,
                                     const double max_cell_size=-1.0,
                                     const double min_dihedral=-1.0,
                                     const double max_dihedral=-1.0,
                                     const double radius_edge_ratio=-1.0,
                                     const int max_mesh_points=-1
) {
  assert(verts.numel() == 3); // otherwise we can't make a 3-D triangulation
  // create the tetgenbehavior object which contains all options/switches for tetrahedralize
  status_update("Creating `tetgenbehavior` object");
  tetgenbehavior tgb;
  tgb.plc = 1; // we will always tetrahedralize a piecewise linear complex
  tgb.quality = 1; // we will (almost) always improve the tetrahedral mesh
  tgb.neighout = 1; // we *need* the neighbour information to be stored into tgo.
  if (max_cell_size > 0){
    // volume constraint with a specified size
    tgb.fixedvolume = 1;
    tgb.maxvolume = max_cell_size;
  } else{
    //volume constraint without a specified size?
    tgb.varvolume = 1;
  }
  if (max_mesh_points>0) tgb.steinerleft = max_mesh_points;
  if (radius_edge_ratio>0) tgb.minratio = radius_edge_ratio;
  if (min_dihedral>0) tgb.mindihedral = min_dihedral;
  if (max_dihedral>0) tgb.optmaxdihedral = max_dihedral;
  #ifndef VERBOSE_MESHING
  tgb.quiet = 1;
  #endif
  #ifdef VERBOSE_MESHING
  tgb.verbose = 10000;
  #endif

  // make the input and output tetgenio objects and fill the input with our polyhedron
  status_update("Creating input and output `tetgenio` objects");
  tetgenio tgi, tgo;
  // we have to handle initializing points/facets, but tetgenio has a destructor
  // which handles deleting all non-NULL fields.
  status_update("Initialize and fill the input object's pointlist parameter");
  tgi.numberofpoints = static_cast<int>(verts.size());
  tgi.pointlist = new double[3*tgi.numberofpoints];
  tgi.pointmarkerlist = new int[tgi.numberofpoints];
  //tgi.point2tetlist = new int[tgi.numberofpoints];
  int idx=0;
  for (size_t i=0; i<verts.size(); ++i){
    tgi.pointmarkerlist[i] = static_cast<int>(i);
    for (size_t j=0; j<verts.numel(); ++j)
      tgi.pointlist[idx++] = verts.getvalue(i,j);
  }
  status_update("Initialize and fill the input object's facetlist parameter");
  tgi.numberoffacets = static_cast<int>(vpf.size());
  tgi.facetlist = new tetgenio::facet[tgi.numberoffacets];
  tgi.facetmarkerlist = new int[tgi.numberoffacets];
  for (size_t i=0; i<vpf.size(); ++i){
    tgi.facetmarkerlist[i] = static_cast<int>(i);
    tgi.facetlist[i].numberofpolygons = 1;
    tgi.facetlist[i].polygonlist = new tetgenio::polygon[1];
    tgi.facetlist[i].polygonlist[0].numberofvertices = static_cast<int>(vpf[i].size());
    tgi.facetlist[i].polygonlist[0].vertexlist = new int[tgi.facetlist[i].polygonlist[0].numberofvertices];
    for (size_t j=0; j<vpf[i].size(); ++j)
      tgi.facetlist[i].polygonlist[0].vertexlist[j] = vpf[i][j];
  }
  // The input is now filled with the piecewise linear complex information.
  // so we can call tetrahedralize:
  status_update("Calling tetgen::tetrahedralize");
  try {
      tetrahedralize(&tgb, &tgi, &tgo);
  } catch (const std::logic_error& e) {
    std::string msg = "tetgen threw a logic error with message\n" + std::string(e.what());
    throw std::runtime_error(msg);
  } catch (const std::runtime_error& e) {
    std::string msg = "tetgen threw a runtime_error with message\n" + std::string(e.what());
    throw std::runtime_error(msg);
  } catch (...) {
    std::string msg = "tetgen threw an undetermined error";
    throw std::runtime_error(msg);
  }
  status_update("Constructing TetrahedralTriangulation object");
  return TetrahedralTriangulation(tgo);
}

#endif // _TRIANGULATION_H_
