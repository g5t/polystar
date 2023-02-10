#ifndef _POLYSTAR_BITMAP_HPP_
#define _POLYSTAR_BITMAP_HPP_
#include <pybind11/pybind11.h>
#include "_array.hpp"
#include "_c_to_python.hpp"
#include "bitmap.hpp"
#include "bitmap_path.hpp"
#include "utilities.hpp"

template<class T, class A>
void define_bitmap_inits(py::class_<A> & cls){
  namespace py = pybind11;
  using namespace polystar;
  using namespace polystar::bitmap;
  using namespace pybind11::literals;
  cls.def(py::init([](const py::array_t<T> &pyv) {
    return A(py2a2(pyv));
  }), "values"_a);
}

template<class T, class A>
void define_bitmap(py::class_<A> & cls){
  namespace py = pybind11;
  using namespace polystar;
  using namespace polystar::bitmap;
  using namespace pybind11::literals;
  cls.def_property_readonly("values", [](const A& p){return a2py(p.values());});
  cls.def("dilate", [](const A& p, int n){return p.dilate(n);});
  cls.def("erode", [](const A& p, int n){return p.erode(n);});
  cls.def("dilate_edges", [](const A& p, int n){return p.dilate_edges(n);});
  cls.def("erode_edges", [](const A& p, int n){return p.erode_edges(n);});
  cls.def("extract_image_polygons", [](const A& p, T level){return p.extract_polygons(coordinates::system::y_down_x_right, level);});
  cls.def("extract_array_polygons", [](const A& p, T level){return p.extract_polygons(coordinates::system::y_up_x_right, level);});
  cls.def("write_image", [](const A&p, std::string filename){return write(p.image(coordinates::system::y_down_x_right), filename);});
  cls.def("write_array", [](const A&p, std::string filename){return write(p.image(coordinates::system::y_up_x_right), filename);});
}

#endif