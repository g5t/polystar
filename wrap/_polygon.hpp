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
  cls.def_property_readonly("border", [](const A& p){return Wire::base_t(p.wires().border());});
  cls.def_property_readonly("wires", [](const A& p){
    auto wires = p.wires().wires();
    // strip off the Wire type information
    std::vector<Wire::base_t> ow;
    ow.reserve(wires.size());
    for (const auto & wire: wires) ow.emplace_back(wire);
    return ow;
  });
//  cls.def_property_readonly("border", [](const A& p){return a2py(p.wires().border());});
//  cls.def_property_readonly("wires", [](const A& p){return a2py(p.wires().wires());});
  cls.def_property_readonly("area", &A::area);
  cls.def_property_readonly("mirror",&A::mirror);
  cls.def_property_readonly("inverse",&A::inverse);
  cls.def_property_readonly("is_convex",&A::is_convex);
  cls.def_property_readonly("centroid",[](const A& p){return a2py(p.centroid());});
  cls.def("contains",[](const A& p, const py::array_t<T>& pyv){
    auto v = py2a2(pyv);
    return p.contains(v);
  });
  cls.def("translate",[](const A& p, const py::array_t<T>& pyv){
    auto v = py2a2(pyv);
    return p.translate(v);
  });
  cls.def("intersects",[](const A & p, const A & o){return p.intersects(o);});

  cls.def("intersection",[](const A& p, const A& o){return polygon_intersection(p, o);});
  cls.def("convex_intersection",[](const A& p, const A& o){return polygon_convex_intersection(p, o);});

  cls.def("simplify",&A::simplify);
  cls.def("triangulate",&A::triangulate);
  cls.def("write_svg", [](const A& p, const std::string filename, const std::string fill, const std::string stroke){
    p.write_svg(filename, std::make_optional(fill), std::make_optional(stroke));
  }, "filename"_a, py::kw_only(), "fill"_a="none", "stroke"_a="black");
  cls.def("to_svg", [](const A& p, const std::string fill, const std::string stroke){
    return p.to_svg(std::make_optional(fill), std::make_optional(stroke));
  }, py::kw_only(), "fill"_a="none", "stroke"_a="black");

  cls.def("transform", [](const A & p, py::array_t<T>& pyv){
    auto info = pyv.request();
    if (info.ndim != 2 || info.size != 4 || info.shape[0] != 2 || info.shape[1] != 2) {
      throw std::runtime_error("polystar::Polygon::transform takes a 2x2 array");
    }
    std::array<T, 4> matrix;
    matrix[0] = ((T *) info.ptr)[0];
    matrix[1] = ((T *) info.ptr)[info.strides[1]/info.itemsize];
    matrix[2] = ((T *) info.ptr)[info.strides[0]/info.itemsize];
    matrix[3] = ((T *) info.ptr)[(info.strides[0] + info.strides[1])/info.itemsize];
    return p.transform(matrix);
  }, "matrix"_a);
  cls.def("skew", &A::skew);

// overloaded operations:
  cls.def("__eq__", [](const A& p, const A& o){return p == o;});
  cls.def("__neq__", [](const A& p, const A& o){return p != o;});
//  cls.def("__add__", [](const A& p, const A& o){return p + o;});
//  cls.def("__add__", [](const A& p, const std::vector<A> & v){return p.combine_all(v);});

}
//
//template<class T, class A>
//void define_svg_animation(pybind11::module &m, py::class_<A> & cls){
//  namespace py = pybind11;
//  using namespace polystar;
//  using namespace polystar::polygon;
//  using namespace pybind11::literals;
//
//  m.def("animated_svg",[](const std::vector<A>& polygons,
//      const std::string & fill, const std::string & stroke,
//      const double fps){
//    std::vector<SVG::SVG> svgs;
//    svgs.reserve(polygons.size());
//    for (const auto & polygon: polygons) svgs.push_back(polygon.to_svg(std::make_optional(fill), std::make_optional(stroke)));
//    auto animated = SVG::frame_animate(svgs, fps);
//    return animated;
//  }, "polygons"_a, py::kw_only(), "fill"_a="none", "stroke"_a="black");
//}

#endif