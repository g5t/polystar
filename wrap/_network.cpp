#include <pybind11/pybind11.h>
#include "_network.hpp"
void wrap_polygon_network(pybind11::module &m){
  using namespace polystar;
  using namespace polystar::polygon;
  namespace py = pybind11;

  py::class_<Network<Wire,double,Array2>> bare(m, "Network");
  define_polygon_network_inits<double>(bare);
  define_polygon_network<double>(bare);

  py::class_<Network<Wire,int,Array2>> coord(m, "CoordinateNetwork");
  define_polygon_network_inits<int>(coord);
  define_polygon_network<int>(coord);

//  py::class_<Poly<double,LVec>> lvec(m, "LPolyhedron");
//  define_polygon<double>(lvec);
//  define_polygon_lvec<double>(lvec);
}
