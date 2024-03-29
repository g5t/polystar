#ifndef _POLYSTAR_POLYHEDRON_HPP_
#define _POLYSTAR_POLYHEDRON_HPP_
#include <pybind11/pybind11.h>
#include "_array.hpp"
#include "_c_to_python.hpp"
#include "polyhedron_flex.hpp"
#include "utilities.hpp"

template<class T, class A>
void define_polyhedron_inits(py::class_<A> & cls){
  namespace py = pybind11;
  using namespace polystar;
  using namespace polystar::polyhedron;
  using namespace pybind11::literals;
  cls.def(py::init([](const py::array_t<T> &pyv) {
    return A(py2a2(pyv));
  }), "vertices"_a);

  cls.def(py::init([](const py::array_t<T> &pyv, const std::vector<std::vector<int>> &faces) {
    return A(py2a2(pyv), Faces(faces));
  }), "vertices"_a, "faces"_a);
}

template<class T, class A>
void define_polyhedron_lvec(py::class_<A> & cls){
  namespace py = pybind11;
  using namespace polystar;
  using namespace polystar::polyhedron;
  using namespace pybind11::literals;
  cls.def("rotate", [](const A &p, const std::array<T,9> &s) {
    return p.apply(s);
  });
}

template<class T, class A>
void define_polyhedron(py::class_<A> & cls){
  namespace py = pybind11;
  using namespace polystar;
  using namespace polystar::polyhedron;
  using namespace pybind11::literals;
  cls.def_property_readonly("vertices", [](const A& p){return a2py(p.vertices());});
  cls.def_property_readonly("points", [](const A& p){return a2py(p.face_points());});
  cls.def_property_readonly("normals", [](const A& p){return a2py(p.face_normals());});
  cls.def_property_readonly("faces", [](const A& p){return p.faces().faces();});
  cls.def_property_readonly("volume", &A::volume);
  cls.def_property_readonly("mirror",&A::mirror);
  cls.def_property_readonly("centre",&A::centre);
  cls.def("intersection",[](const A& p, const A& o){return p.intersection(o);});
}

#endif