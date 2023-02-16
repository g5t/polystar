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
  cls.def("polygons",[](const A& net){
    std::vector<Poly<T,Array2>> p;
    p.reserve(net.size());
    for (const auto & w: net.wires()) {
      p.emplace_back(net.vertices(), Wires(w));
    }
    return p;
  });
  cls.def("path",[](const A& net, const py::array_t<T>& pyfrom, const py::array_t<T>& pyto){
    auto from = py2a2(pyfrom);
    auto to = py2a2(pyto);
    return a2py(net.path(from, to));
  });
  cls.def("path_to_each",[](const A& net, const py::array_t<T>& pyfrom, const py::array_t<T>& pyto){
    auto from = py2a2(pyfrom);
    auto to = py2a2(pyto);
    auto costs_paths = net.path_to_each(from, to);
    std::vector<double> costs;
    std::vector<py::array_t<T>> paths;
    costs.reserve(costs_paths.size());
    paths.reserve(costs_paths.size());
    for (const auto & [cost, path]: costs_paths){
      costs.push_back(cost);
      paths.push_back(a2py(path));
    }
    return std::make_tuple(costs, paths);
  }, "from"_a, "to"_a);

  cls.def("write_svg", [](const A& p, const std::string filename, const std::string fill, const std::string stroke){
    p.write_svg(filename, std::make_optional(fill), std::make_optional(stroke));
  }, "filename"_a, py::kw_only(), "fill"_a="none", "stroke"_a="black");
  cls.def("to_svg", [](const A& p, const std::string fill, const std::string stroke){
    return p.to_svg(std::make_optional(fill), std::make_optional(stroke));
  }, py::kw_only(), "fill"_a="none", "stroke"_a="black");

}

#endif