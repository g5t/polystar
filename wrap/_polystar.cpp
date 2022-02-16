/* This file is part of brille.

Copyright © 2022 Greg Tucker <greg.tucker@ess.eu>

brille is free software: you can redistribute it and/or modify it under the
terms of the GNU Affero General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option)
any later version.

brille is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with brille. If not, see <https://www.gnu.org/licenses/>.            */
#include "version.hpp"

#include <pybind11/pybind11.h>

void wrap_debug(pybind11::module &);
void wrap_polyhedron(pybind11::module &);

void wrap_version(pybind11::module & m){
  using namespace brille::version;
  m.attr("__version__") = version_number;
  std::string v = version_number;
  if (!std::string(git_revision).empty()){
    v += "+" + std::string(git_branch);
    v += "." + std::string(git_revision).substr(0,7);
  }
  m.attr("version") = v;
  m.attr("git_revision") = git_revision;
  m.attr("build_datetime") = build_datetime;
  m.attr("build_hostname") = build_hostname;
}

PYBIND11_MODULE(_polystar, m){
  m.doc() = R"pbdoc(
    pybind11 module :py:mod:`polystar._polystar`
    ----------------------------------------
    This module provides the interface to the C++ library.

    All of the symbols defined within :py:mod:`polystar._polystar`
		are imported by :py:mod:`polystar` to make using them easier.
    If in doubt, the interfaced classes can be accessed via their submodule
    syntax.

    .. code-block:: python

			from polystar._polystar import Polyhedron
      from polystar.plotting import plot

			poly = Polyhedron(...)

      plot(poly)

    .. currentmodule::polystar._polystar

    .. autosummary::
      :toctree: _generate

  )pbdoc";
  wrap_version(m);
  wrap_polyhedron(m);
  wrap_debug(m);
}
