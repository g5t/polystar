# The compiled module is installed at the base of the package directory structure such that its properties
# can be included here. If the module is loaded, however, without being installed first these import statements will
# raise a ModuleNotFoundError, which prevents automated pre-installation testing from taking place.
try:
    from ._polystar import __version__, version, ApproxConfig
    from ._polystar import Polyhedron
    from ._polystar import Polygon, CoordinatePolygon
    from ._polystar import BitmapI, BitmapF, BitmapD
    from ._polystar import Network, CoordinateNetwork
    from ._polystar import SVG, animated_svg

    # Store types for use in, e.g., the plotting routines
    __polygon_types__ = (Polygon, CoordinatePolygon)
    __polyhedron_types__ = (Polyhedron, )
    __bitmap_types__ = (BitmapI, BitmapF, BitmapD)
    __network_types__ = (Network, CoordinateNetwork)

except ModuleNotFoundError:
    pass
