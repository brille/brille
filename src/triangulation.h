#ifndef _TRIANGULATION_H_
#define _TRIANGULATION_H_
/*----------------------------- HEADER INCLUDES ------------------------------*/
// combined:
#include <vector>
#include <array>
#include <cassert>
// CGAL headers
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
// Tetgen header
#include "tetgen/tetgen.h"
/*--------------------------------- TYPEDEFS ---------------------------------*/
// for CGAL Delaunay (re)Triangulation w/ Fast_location and an index per point
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_3<size_t, K>      Vertex_Base;
typedef CGAL::Delaunay_triangulation_cell_base_3<K>                  Cell_Base;
typedef CGAL::Triangulation_data_structure_3<Vertex_Base, Cell_Base> Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location>  Delaunay;

/*! Check whether two vectors contain consecutive equivalent elements in any
position, even [last, first].
E.g.,  [4 5 3 1] and [2 4 5] both have consecutive element values 4 and 5;
[9 3 2 0] and [3 4 8 0 9 10] both have consecutive element values 0 and 9.
*/
template<typename T, typename=typename std::enable_if<std::is_integral<T>::value>::type>
bool equal_consecutive_elements(const std::vector<T>& a, const std::vector<T>& b){
  for (size_t i=0; i<a.size(); ++i)
  for (size_t j=0; j<b.size(); ++j)
  if (a[i]==b[j] && a[(i+1)%a.size()]==b[(j+1)%b.size()]) return true;
  return false;
}

template<typename T, typename=typename std::enable_if<std::is_integral<T>::value>::type>
std::vector<std::vector<T>> check_halfedge_order(const std::vector<std::vector<T>>& a){
  std::vector<std::vector<T>> b(a.size());
  if (a.size()<1) return b;
  b[0] = a[0];
  bool flip;
  std::vector<T> ai, aj;
  // check each vector past the first
  for (size_t i=1; i<a.size(); ++i){
    // against all vectors before it
    for (size_t j=0; j<i; ++j){
      // whether any two consecutive elements in each vector are equivalent
      flip = equal_consecutive_elements(a[i], a[j]);
      if (flip) break;
    }
    if (flip) // equivalent consecutive elements, so we need to flip the vector order
    for (size_t ii=a[i].size(); ii>0; --ii) b[i].push_back(a[i][ii-1]);
    else
    b[i] = a[i];
  }
  // check that the halfedge order is now ok:
  for (size_t i=1; i<b.size(); ++i) for (size_t j=0; j<i; ++j) if (equal_consecutive_elements(b[i], b[j])){
    std::string msg = "The half-edge orders [";
    for (auto k: a[i]) msg += " " + std::to_string(k);
    msg += " ] and [";
    for (auto k: a[j]) msg += " " + std::to_string(k);
    msg += "\nwhich have been reordered to [";
    for (auto k: b[i]) msg += " " + std::to_string(k);
    msg += " ] and [";
    for (auto k: b[j]) msg += " " + std::to_string(k);
    msg += " ] are incompatible.";
    throw std::runtime_error(msg);
  }
  return b;
}

//a functor that returns a std::pair<Point,unsigned>.
//the unsigned integer is incremented at each call to
//operator()
struct Auto_count : public CGAL::cpp98::unary_function<const Delaunay::Point&,std::pair<Delaunay::Point,unsigned> >{
  mutable unsigned i;
  Auto_count() : i(0){}
  std::pair<Delaunay::Point,unsigned> operator()(const Delaunay::Point& p) const {
    return std::make_pair(p,i++);
  }
};

#define VERBOSE_MESHING

inline static void status_update(const std::string& status){
  #ifdef VERBOSE_MESHING
  std::cout << status << std::endl;
  #endif
}

template <typename T>
Delaunay triangulate(const ArrayVector<T>& verts,
                     const std::vector<std::vector<int>>& vpf,
                     const double max_cell_size=-1.0,
                     const double min_dihedral=-1.0,
                     const double max_dihedral=-1.0,
                     const double radius_edge_ratio=-1.0,
                     const int max_mesh_points=-1
                   )
{
  assert(verts.numel() == 3); // otherwise we can't make a 3-D triangulation
  // create the tetgenbehavior object which contains all options/switches for tetrahedralize
  status_update("Creating `tetgenbehavior` object");
  tetgenbehavior tgb;
  tgb.plc = 1; // we will always tetrahedralize a piecewise linear complex
  tgb.quality = 1; // we will (almost) always improve the tetrahedral mesh
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

  // get just the vertex locations from the generated mesh
  status_update("Pulling out the mesh vertex points" );
  std::vector<Delaunay::Point> points;
  for (int i=0; i < tgo.numberofpoints; ++i)
    points.push_back(Delaunay::Point(tgo.pointlist[3*i], tgo.pointlist[3*i+1], tgo.pointlist[3*i+2]));
  status_update("In order to make a new fast-access CGAL triangulation with an index per vertex.");
  // following 5.3.3 https://doc.cgal.org/latest/Triangulation_3/index.html#Chapter_3D_Triangulations
  // to add an unsigned integer for each point in the triangulation.
  Delaunay dtri(boost::make_transform_iterator(points.begin(),Auto_count()),
                boost::make_transform_iterator(points.end(),  Auto_count()));
  status_update("Triangulation complete.");
  return dtri;
}

#endif // _TRIANGULATION_H_
