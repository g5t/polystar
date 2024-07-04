#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "polygon_wire.hpp"

void wrap_wire(pybind11::module &m){
  namespace py = pybind11;
  using namespace pybind11::literals;

  py::class_<polystar::polygon::Wire> wire(m, "Wire");
  wire.def("__len__", [](polystar::polygon::Wire& w){return w.size();});
  wire.def("__getitem__", [](polystar::polygon::Wire& w, size_t i){return w.edge(i);});
  wire.def("edge", &polystar::polygon::Wire::edge, "index"_a);


//
//  svg.def("__str__",[](SVG::SVG & a){return std::string(a);});
//  svg.def("write",[](SVG::SVG & a, const std::string& filename){
//    auto of = std::ofstream(filename);
//    of << std::string(a);
//  },"filename"_a);

//  m.def("animate1", [](std::vector<SVG::SVG> & a, double fps){
//    return SVG::frame_animate(a, fps);
//  }, "frames"_a, "frames_per_second"_a=60);
//  m.def("animate", &SVG::frame_animate, "frames"_a, "frames_per_second"_a=60);
}
