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
