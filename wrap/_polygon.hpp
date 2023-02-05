#ifndef _POLYSTAR_POLYGON_HPP_
#define _POLYSTAR_POLYGON_HPP_
#include <pybind11/pybind11.h>
#include "_array.hpp"
#include "_c_to_python.hpp"
#include "polygon.hpp"
#include "utilities.hpp"

template<class T, class A>
void define_polygon_inits(py::class_<A> & cls){
  namespace py = pybind11;
  using namespace polystar;
  using namespace polystar::polygon;
  using namespace pybind11::literals;
  cls.def(py::init([](const py::array_t<T> &pyv) {
    return A(py2a2(pyv));
  }), "vertices"_a);

  cls.def(py::init([](const py::array_t<T> &pyv, const std::vector<int> &border) {
    auto v = py2a2(pyv);
    return A(v, Wires(border));
  }), "vertices"_a, "border"_a);

  cls.def(py::init([](const py::array_t<T> &pyv, const std::vector<int> & border, const std::vector<std::vector<int>> &wires) {
    auto v = py2a2(pyv);
    return A(v, Wires(border, wires));
  }), "vertices"_a, "border"_a, "wires"_a);
}

template<class T, class A>
void define_polygon(py::class_<A> & cls){
  namespace py = pybind11;
  using namespace polystar;
  using namespace polystar::polygon;
  using namespace pybind11::literals;
  cls.def_property_readonly("vertices", [](const A& p){return a2py(p.vertices());});
  cls.def_property_readonly("border", [](const A& p){return p.wires().border();});
  cls.def_property_readonly("wires", [](const A& p){return p.wires().wires();});
//  cls.def_property_readonly("border", [](const A& p){return a2py(p.wires().border());});
//  cls.def_property_readonly("wires", [](const A& p){return a2py(p.wires().wires());});
  cls.def_property_readonly("area", &A::area);
  cls.def_property_readonly("mirror",&A::mirror);
  cls.def_property_readonly("centroid",&A::centroid);
//  cls.def("intersection",[](const A& p, const A& o){return p.intersection(o);});
}

#endif