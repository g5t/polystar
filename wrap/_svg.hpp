#pragma once
#include <pybind11/pybind11.h>
#include "_array.hpp"
#include "_c_to_python.hpp"
#include "svg.hpp"
#include "utilities.hpp"

template<class A>
void define_svg_animation(pybind11::module &m){
  namespace py = pybind11;
  using namespace polystar;
  using namespace pybind11::literals;

  m.def("animated_svg",[](const std::vector<A>& polygons,
                          const std::string & fill, const std::string & stroke,
                          const double fps){
    std::vector<SVG::SVG> svgs;
    svgs.reserve(polygons.size());
    for (const auto & polygon: polygons) {
      svgs.push_back(polygon.to_svg(std::make_optional(fill), std::make_optional(stroke)));
      svgs.back().autoscale();
    }
    auto animated = SVG::frame_animate(svgs, fps);
    animated.autoscale();
    return animated;
  }, "polygons"_a, py::kw_only(), "fill"_a="none", "stroke"_a="black", "fps"_a=60.);
}
