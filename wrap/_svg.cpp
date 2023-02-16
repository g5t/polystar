#include <pybind11/pybind11.h>
//#include <pybind11/stl.h>
#include "_svg.hpp"

void wrap_svg(pybind11::module &m){
  namespace py = pybind11;
  using namespace pybind11::literals;

  py::class_<SVG::SVG> svg(m, "SVG");
//
  svg.def("__str__",[](SVG::SVG & a){return std::string(a);});
  svg.def("write",[](SVG::SVG & a, const std::string& filename){
    auto of = std::ofstream(filename);
    of << std::string(a);
  },"filename"_a);

//  m.def("animate1", [](std::vector<SVG::SVG> & a, double fps){
//    return SVG::frame_animate(a, fps);
//  }, "frames"_a, "frames_per_second"_a=60);
//  m.def("animate", &SVG::frame_animate, "frames"_a, "frames_per_second"_a=60);
}
