/* This file is part of polystar.

Copyright © 2019,2020 Greg Tucker <greg.tucker@stfc.ac.uk>

polystar is free software: you can redistribute it and/or modify it under the
terms of the GNU Affero General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option)
any later version.

polystar is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with polystar. If not, see <https://www.gnu.org/licenses/>.            */
#include <pybind11/pybind11.h>
#include "_polygon.hpp"
#include "_svg.hpp"
void wrap_polygon(pybind11::module &m){
  using namespace polystar;
  using namespace polystar::polygon;
  namespace py = pybind11;

  py::class_<Poly<double,Array2>> bare(m, "Polygon");
  define_polygon_inits<double>(bare);
  define_polygon<double>(bare);

  py::class_<Poly<int,Array2>> coord(m, "CoordinatePolygon");
  define_polygon_inits<int>(coord);
  define_polygon<int>(coord);

  define_svg_animation<Poly<double,Array2>>(m);
  define_svg_animation<Poly<int,Array2>>(m);

//  py::class_<Poly<double,LVec>> lvec(m, "LPolyhedron");
//  define_polygon<double>(lvec);
//  define_polygon_lvec<double>(lvec);
}
