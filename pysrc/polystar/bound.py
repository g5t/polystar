from importlib.util import find_spec
bound_spec = find_spec('polystar._polystar')
if bound_spec is not None:
    # The compiled module is installed at the base of the package directory structure such that its properties
    # can be included here. If the module is loaded, however, without being installed first these import statements will
    # raise a ModuleNotFoundError, which prevents automated pre-installation testing from taking place.
    try:
        from ._polystar import (
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
        )

        # Store types for use in, e.g., the plotting routines
        __polygon_types__ = (
            Polygon,
            CoordinatePolygon,
        )
        __polyhedron_types__ = (
            Polyhedron,
        )
        __bitmap_types__ = (
            BitmapI,
            BitmapF,
            BitmapD,
        )
        __network_types__ = (
            Network,
            CoordinateNetwork,
        )

        __all__ = [
            "__version__",
            "version",
            "ApproxConfig",
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
            "__polygon_types__",
            "__polyhedron_types__",
            "__bitmap_types__",
            "__network_types__",
        ]

    except ImportError as import_error:
        import platform

        if platform.system == "Windows" and "DLL" in import_error.msg:
            msg = "You may be missing the latest Visual C++ redistributable package"
            msg += " install it from Microsoft @ https://support.microsoft.com/en-us/help/2977003"
            msg += " before trying to import polystar again"
            print(msg)
        raise import_error
else:
    __version__ = '0.0.0'
    __all__ = ['__version__']
