def network_plot(network, ax=None, set_limits=False, show=False, **kwargs):
    import matplotlib.pyplot as pplot
    if ax is None:
        fig, ax = pplot.subplots()
    for polygon in network.polygons():
        polygon_plot(polygon, ax=ax, **kwargs)
    if set_limits:
        from numpy import min, max, vstack
        vertices = vstack([p.vertices for p in network.polygons()])
        ll = min(vertices, axis=0)
        ul = max(vertices, axis=0)
        ax.set_xlim(ll[0], ul[0])
        ax.set_ylim(ll[1], ul[1])
        ax.set_aspect(1.0)
    if show:
        pplot.show()
    return ax


def polygon_plot(polygon, ax=None, set_limits=False, show=False, **kwargs):
    import matplotlib.pyplot as pplot
    if ax is None:
        fig, ax = pplot.subplots()

    ax.add_patch(polygon_patch(polygon, **kwargs))
    if set_limits:
        from numpy import min, max
        ll = min(polygon.vertices, axis=0)
        ul = max(polygon.vertices, axis=0)
        ax.set_xlim(ll[0], ul[0])
        ax.set_ylim(ll[1], ul[1])
        ax.set_aspect(1.0)
    if show:
        pplot.show()
    return ax


def polygon_patch(polygon, **kwargs):
    from polystar.bound import __polygon_types__
    if not isinstance(polygon, __polygon_types__):
        print("Only plotting of polygon types supported")

    import matplotlib.path as mpath
    import matplotlib.patches as mpatches

    import numpy
    codes_vertices = [wire_codes_vertices(polygon.vertices, polygon.border)]
    codes_vertices.extend([wire_codes_vertices(polygon.vertices, wire) for wire in polygon.wires])
    codes = numpy.hstack([c for c, v in codes_vertices])
    verts = numpy.vstack([v for c, v in codes_vertices])

    path = mpath.Path(verts, codes)
    patch = mpatches.PathPatch(path, **kwargs)
    return patch


def wire_codes_vertices(all_vertices, wire):
    import numpy
    import matplotlib.path as mpath
    codes = numpy.ones(len(wire)+1, dtype=mpath.Path.code_type) * mpath.Path.LINETO
    codes[0] = mpath.Path.MOVETO
    vertices = all_vertices[wire]
    vertices = numpy.vstack((vertices, all_vertices[wire[0]]))
    return codes, vertices




def make_colours(n, color=None):
    """Construct a list of colors for use in displaying Polygon objects

    Parameters
    ----------
    n : int
        The number of colors required
    color : Union[List[Union[str, tuple[numbe, number, number]], Union[str, tuple[number, number, number]]
        If color is not provided, the list of all colors known to matplotlib will be used

    Returns
    -------
    List[colors]
        A length-n list of colors. If the starting value of `color` is less than length-n, it will be tiled to length-n.

    Examples
    --------
    >>> make_colours(7, ['red', 'blue', 'green'])
    ['red', 'blue', 'green', 'red', 'blue', 'green', 'red']

    >>> make_colours(4, 'black')
    ['black', 'black', 'black', 'black']
    """
    if color is None:
        from matplotlib.colors import get_named_colors_mapping
        color = get_named_colors_mapping()

    from collections.abc import Iterable
    if isinstance(color, Iterable):
        color = list(color)
    if isinstance(color, str) or (isinstance(color, (list, tuple)) and len(color) == 3):
        color = [color]

    from numpy import ndarray, array, tile
    if not isinstance(color, ndarray):
        color = array(color)
    if color.shape[0] < n:
        color = tile(color, 1+n//color.shape[0])
    return color[0:n]
