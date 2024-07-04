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
#include "version.hpp"

#include <pybind11/pybind11.h>

void wrap_debug(pybind11::module &);
void wrap_polyhedron(pybind11::module &);
void wrap_approx(pybind11::module &);
void wrap_polygon(pybind11::module &);
void wrap_bitmap(pybind11::module &);
void wrap_polygon_network(pybind11::module &);
void wrap_svg(pybind11::module &);
void wrap_wire(pybind11::module &);

void wrap_version(pybind11::module & m){
  using namespace polystar::version;
  m.attr("__version__") = version_number;
  m.attr("version") = meta_version;
  m.attr("git_revision") = git_revision;
  m.attr("git_branch") = git_branch;
  m.attr("build_datetime") = build_datetime;
  m.attr("build_hostname") = build_hostname;
}

PYBIND11_MODULE(_polystar, m){
  m.doc() = R"pbdoc(
    pybind11 module :py:mod:`polystar._polystar`
    ----------------------------------------
    This module provides the interface to the C++ library.

    All of the symbols defined within :py:mod:`polystar._polystar` are imported by
    :py:mod:`polystar` to make using them easier.
    If in doubt, the interfaced classes can be accessed via their submodule
    syntax.

    .. currentmodule:: polystar._polystar

    .. autosummary::
      :toctree: _generate

  )pbdoc";
  wrap_version(m);
  wrap_svg(m);
  wrap_polyhedron(m);
  wrap_debug(m);
  wrap_approx(m);
  wrap_polygon(m);
  wrap_bitmap(m);
  wrap_polygon_network(m);
  wrap_wire(m);
}
