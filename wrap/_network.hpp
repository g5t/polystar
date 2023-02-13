#ifndef _POLYSTAR_POLYGON_NETWORK_HPP_
#define _POLYSTAR_POLYGON_NETWORK_HPP_

#include <pybind11/pybind11.h>
#include "_array.hpp"
#include "_c_to_python.hpp"
#include "polygon.hpp"
#include "utilities.hpp"

template<class T, class A>
void define_polygon_network_inits(py::class_<A> & cls){
  namespace py = pybind11;
  using namespace polystar;
  using namespace polystar::polygon;
  using namespace pybind11::literals;

  cls.def(py::init([](const Poly<T,Array2> p){
    return p.triangulate();
  }), "polygon"_a);
}

template<class T, class A>
void define_polygon_network(py::class_<A> & cls){
  namespace py = pybind11;
  using namespace polystar;
  using namespace polystar::polygon;
  using namespace pybind11::literals;

  cls.def("simplify",&A::simplify);
  cls.def("wires",&A::wires);
  cls.def("polygons",&A::polygons);
  cls.def("path",[](const A& net, const py::array_t<T>& pyfrom, const py::array_t<T>& pyto){
    auto from = py2a2(pyfrom);
    auto to = py2a2(pyto);
    return a2py(net.path(from, to));
  });
}

#endif