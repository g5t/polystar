# Copyright © 2023 Greg Tucker <greggory.tucker@gmail.com>
#
# This file is part of polystar.
#
# polystar is free software: you can redistribute it and/or modify it under the
# terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# polystar is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.
#
# See the GNU Affero General Public License for more details.
# You should have received a copy of the GNU Affero General Public License
# along with polystar. If not, see <https://www.gnu.org/licenses/>.

"""Python module :py:mod:`polystar`
=================================

This module provides access to the C++ polystar library which can be used to
interact with polyhedron and polygon routines.

.. currentmodule:: polystar

.. autosummary::
    :toctree: _generate
"""
from .bound import (
    __version__,
    version,
    ApproxConfig,
    Polyhedron,
    Polygon,
    CoordinatePolygon,
    BitmapI,
    BitmapF,
    BitmapD,
    Network,
    CoordinateNetwork,
    SVG,
    animated_svg,
    __polygon_types__,
    __polyhedron_types__,
    __bitmap_types__,
    __network_types__,
)

__all__ = [
    "Polyhedron",
    "Polygon",
    "CoordinatePolygon",
    "BitmapI",
    "BitmapF",
    "BitmapD",
    "Network",
    "CoordinateNetwork",
    "SVG",
    "animated_svg",
]
