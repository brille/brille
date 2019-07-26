#ifndef _TRIANGULATION_H_
#define _TRIANGULATION_H_
/*----------------------------- HEADER INCLUDES ------------------------------*/
// combined:
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <vector>
#include <array>
#include <cassert>
// polyhedron creation
#include <CGAL/Polyhedron_incremental_builder_3.h>
// mesh generation
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>
/*--------------------------------- TYPEDEFS ---------------------------------*/
// combined:
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
// Mesh Domain
typedef CGAL::Mesh_polyhedron_3<K>::type Polyhedron; // This inherits from Polyhedron_3, I think.
typedef Polyhedron::HalfedgeDS HalfedgeDS; // needed for the polyhedron incremental builder
typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> Mesh_domain;
#ifdef CGAL_CONCURRENT_MESH_3
  typedef CGAL::Parallel_tag Concurrency_tag;
#else
  typedef CGAL::Sequential_tag Concurrency_tag;
#endif
// Mesh Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default,Concurrency_tag>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr,Mesh_domain::Corner_index,Mesh_domain::Curve_index> C3t3;
// Mesh Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
// Delaunay (re)Triangulation w/ Fast_location and an index per point
typedef CGAL::Triangulation_vertex_base_with_info_3<size_t, K>      Vertex_Base;
typedef CGAL::Delaunay_triangulation_cell_base_3<K>                  Cell_Base;
typedef CGAL::Triangulation_data_structure_3<Vertex_Base, Cell_Base> Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location>  Delaunay;

template <class HDS>
class Build_IrBZ_Polyhedron : public CGAL::Modifier_base<HDS> {
  typedef typename HDS::Vertex::Point Point;
  std::vector<Point> verts;
  std::vector<std::vector<int>> facets;
public:
  Build_IrBZ_Polyhedron(const ArrayVector<double>& v, const std::vector<std::vector<int>>& vpf): facets(vpf) {
    for (size_t i=0; i<v.size(); ++i)
      verts.push_back(Point(v.getvalue(i,0), v.getvalue(i,1), v.getvalue(i,2)));
  }
  void operator()(HDS& hds){
    // initialize the incremental builder
    CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);
    // initalize the polyhedron with number of vertices and number of facets
    B.begin_surface(verts.size(), facets.size());
    // add all vertices
    for (auto v: verts) B.add_vertex(v);
    // then add the facets
    for (auto facet: facets){
      B.begin_facet();
      for (auto i: facet) B.add_vertex_to_facet(i);
      B.end_facet();
    }
    // and finalize the polyhedron
    B.end_surface();
  }
};

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

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

#define VERBOSE_MESHING

template <typename T>
Delaunay triangulate(const ArrayVector<T>& verts,
                     const std::vector<std::vector<int>>& vpf,
                     const double in_max_cell_size,
                     const double in_min_facet_angle,
                     const double in_cell_radius_edge_ratio)
{
  // default values taken from CGAL/Mesh_criteria_3.h which has each as a
  // Tr::Geom_traits::FT([value]) where the CGAL documentation indicates that
  // undefined named parameters take the value `ignored`
  double max_cell_size = in_max_cell_size > 0 ? in_max_cell_size: 0.; // std::cbrt(3*in_max_vol/4/PI) : 0.; // cbrt == cube-root
  double min_facet_angle = in_min_facet_angle > 0 ? in_min_facet_angle : 0.;
  double max_cell_radius_edge_ratio = in_cell_radius_edge_ratio > 0 ? in_cell_radius_edge_ratio : 0.;
  assert(verts.numel() == 3); // otherwise we can't make a 3-D triangulation
  // create the polyhedron builder
  Build_IrBZ_Polyhedron<HalfedgeDS> irbz_polyhedron(verts, check_halfedge_order(vpf));
  // and the polyhedron itself
  Polyhedron polyhedron;
  polyhedron.delegate(irbz_polyhedron);
  // create the mesh domain
  Mesh_domain domain(polyhedron);
  // automatically find and add the sharp polyhedron features
  domain.detect_features();
  // setup the meshing criteria:
  /* The parameters are named parameters and can be passed in any order provided
  their name is given (see example below). The name of each parameter is the one
  that is written in the description of the function
  (e.g. parameters::facet_size).

  The description of each parameter is as follows:
  edge_size: a scalar field (resp. a constant) providing a space varying
             (resp. a uniform) upper bound for the lengths of curve edges.
             This parameter has to be set to a positive value when 1-dimensional
             features protection is used.
  facet_angle: a lower bound for the angles (in degrees) of the surface mesh
               facets.
  facet_size: a scalar field (resp. a constant) describing a space varying
              (resp. a uniform) upper-bound or for the radii of the surface
              Delaunay balls.
  facet_distance: a scalar field (resp. a constant) describing a space varying
                  (resp. a uniform) upper bound for the distance between the
                  facet circumcenter and the center of its surface Delaunay
                  ball.
  facet_topology: the set of topological constraints which have to be verified
                  by each surface facet. The default value is
                  CGAL::FACET_VERTICES_ON_SURFACE. See Mesh_facet_topology
                  manual page to get all possible values.
  cell_radius_edge_ratio: an upper bound for the radius-edge ratio of the mesh
                          tetrahedra.
  cell_size: a scalar field (resp. a constant) describing a space varying
             (resp. a uniform) upper-bound for the circumradii of the mesh
             tetrahedra.
  Note that each size or distance parameter can be specified using two ways:
  either as a scalar field or as a numerical value when the field is uniform.

  Each parameter has a special default value ignored which means that the
  corresponding criterion will be ignored. Numerical sizing or distance values,
  as well as scalar fields should be given in the unit used for coordinates of
  points in the mesh domain class of the mesh generation process.             */
  Mesh_criteria criteria(facet_angle=min_facet_angle,
                         cell_radius_edge_ratio=max_cell_radius_edge_ratio,
                         cell_size=max_cell_size
                        );
  // Mesh_criteria criteria(edge_size = 0.025, facet_angle = 25, facet_size = 0.05,
  //                        facet_distance = 0.005, cell_radius_edge_ratio =3,
  //                        cell_size = max_vol);
  // Perform the mesh generation
  #ifdef VERBOSE_MESHING
    std::cout << "Calling CGAL::make_mesh_3" << std::endl;
  #endif
  C3t3 c3t3;
  try {
      c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);
  } catch (const std::overflow_error& e) {
    std::string msg = "make_mesh_3 threw an overflow_error with message\n" + std::string(e.what());
    throw std::runtime_error(msg);
  } catch (const std::runtime_error& e) {
    std::string msg = "make_mesh_3 threw a runtime_error with message\n" + std::string(e.what());
    throw std::runtime_error(msg);
  } catch (const std::exception& e) {
    std::string msg = "make_mesh_3 threw an exception with message\n" +std::string(e.what());
    throw std::runtime_error(msg);
  } catch (...) {
    std::string msg = "make_mesh_3 threw an undetermined error";
  }

  // get just the triangulation from the generated mesh
  #ifdef VERBOSE_MESHING
    std::cout << "Getting it's triangulation" << std::endl;
  #endif
  Tr tri = c3t3.triangulation();
  // and pull out the vertex points
  #ifdef VERBOSE_MESHING
    std::cout << "Pulling out the mesh vertex points" << std::endl;
  #endif
  std::vector<Delaunay::Point> points;
  for (auto vit = tri.finite_vertices_begin(); vit != tri.finite_vertices_end(); ++vit){
    auto point = vit->point();
    points.push_back(Delaunay::Point(point[0], point[1], point[2]));
  }
  #ifdef VERBOSE_MESHING
    std::cout << "In order to make a new fast-access triangulation with an index per vertex." << std::endl;
  #endif
  // following 5.3.3 https://doc.cgal.org/latest/Triangulation_3/index.html#Chapter_3D_Triangulations
  // to add an unsigned integer for each point in the triangulation.
  Delaunay dtri(boost::make_transform_iterator(points.begin(),Auto_count()),
                boost::make_transform_iterator(points.end(),  Auto_count()));
  #ifdef VERBOSE_MESHING
    std::cout << "Triangulation complete." << std::endl;
  #endif
  return dtri;
}

#endif // _TRIANGULATION_H_
