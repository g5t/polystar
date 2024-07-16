import numpy as np


def gauss2d(nx, ny, dtype='float'):
    x, y = np.meshgrid(np.linspace(-2, 2, nx), np.linspace(-2, 2, ny), indexing='ij')
    x2 = x * x
    y2 = y * y
    z = np.exp(-x2 - y2)
    if 'int' in dtype:
        z = 255 * (z - np.min(z)) / (np.max(z) - np.min(z))
    return z.astype(dtype)


def peaks2d(nx, ny, dtype='float'):
    x, y = np.meshgrid(np.linspace(-2, 2, nx), np.linspace(-2, 2, ny), indexing='ij')
    x2 = x * x
    x3 = x2 * x
    y2 = y * y
    y5 = y2 * y * y * y
    z = 3 * (1 - x) ** 2 * np.exp(-x2 - (y + 1) ** 2) - 10 * (x / 5 - x3 - y5) * np.exp(-x2 - y2) - 1 / 3 * np.exp(
        -(x + 1) ** 2 - y ** 2)
    if 'int' in dtype:
        z = 255 * (z - np.min(z)) / (np.max(z) - np.min(z))
    return z.astype(dtype)


def square(nx, ny, dtype='float'):
    x, y = np.meshgrid(np.linspace(-2, 2, nx), np.linspace(-2, 2, ny), indexing='ij')
    z = 0 * x
    z[np.logical_and(np.abs(x) < 1, np.abs(y) < 1)] = 1
    return z.astype(dtype)


def circle(nx, ny, dtype='float'):
    x, y = np.meshgrid(np.linspace(-2, 2, nx), np.linspace(-2, 2, ny), indexing='ij')
    z = 0 * x
    z[(x * x + y * y) < 1] = 1
    return z.astype(dtype)


def diamond(nx, ny, dtype='float'):
    x, y = np.meshgrid(np.linspace(-2, 2, nx), np.linspace(-2, 2, ny), indexing='ij')
    z = 0 * x
    z[np.logical_and(np.abs(x + y) < 1, np.abs(x - y) < 1)] = 1
    return z.astype(dtype)


def multi(nx, ny, dtype='float'):
    x, y = np.meshgrid(np.linspace(-2, 2, nx), np.linspace(-2, 2, ny), indexing='ij')
    z = 0 * x
    # Diamond at (1,1)
    z[np.logical_and(np.abs((x - 1) + (y - 1)) < 0.2, np.abs((x - 1) - (y - 1)) < 0.2)] = 1
    # Square at (-0.5, 0.8)
    z[np.logical_and(np.abs(x + 0.5) < 0.3, np.abs(y - 0.8) < 0.3)] = 1
    # Circle at (0, -1.1)
    z[(x * x + (y + 1.1) ** 2) < 0.25] = 1
    return z.astype(dtype)


def board():
    from polystar import BitmapI, CoordinatePolygon
    bm = BitmapI(multi(300, 300, 'int'))
    polygons = bm.extract_image_polygons(1)
    holes = [p.inverse for p in polygons]
    board = CoordinatePolygon(np.array([[0, 0], [0, 300], [300, 300], [300, 0]], dtype='int'))
    return board + holes


def array():
    from polystar import BitmapI, CoordinatePolygon
    bm = BitmapI(multi(300, 300, 'int'))
    polygons = bm.extract_array_polygons(1)
    holes = [p.inverse for p in polygons]
    board = CoordinatePolygon(np.array([[0, 0], [0, 300], [300, 300], [300, 0]], dtype='int'))
    return board + holes


def polylineplot(pl, *args, **kwargs):
    import matplotlib.pyplot as pp
    pp.plot(pl[:, 0], pl[:, 1], *args, **kwargs)


def polyplot(p, indexes=False):
    import matplotlib.pyplot as pp
    b = p.border
    b.append(b[0])
    polylineplot(p.vertices[b])
    for w in p.wires:
        w.append(w[0])
        polylineplot(p.vertices[w])
    if indexes:
        for i, x in enumerate(p.vertices):
            pp.text(*x, f'{i}')


# # This 'worked' before switching to Triangle.
# # Now the board *SHOULD NOT* be simplified before triangulating, as we lose hole-information
#def from_file(filename, dilate=0):
#  import _polystar as p
#  import numpy as np
#  loaded = np.load(filename)
#  bitmap = p.BitmapI(loaded)
#  if dilate > 0:
#    bitmap = bitmap.dilate(dilate)
#
#  holes = [poly.inverse for poly in bitmap.extract_image_polygons(1)]
#  lx, ly = loaded.shape
#  border = np.array([[0,0],[0,lx],[ly,lx],[ly,0]])
#  board = p.CoordinatePolygon(border) + holes
#  simple = board.simplify()
#  triangulated = simple.triangulate()
#  return board, simple, triangulated


def from_file(filename, dilate=0):
    import polystar as p
    import numpy as np
    loaded = np.load(filename)
    bitmap = p.BitmapI(loaded)
    if dilate > 0:
        bitmap = bitmap.dilate(dilate)

    holes = [poly.inverse for poly in bitmap.extract_image_polygons(1)]
    lx, ly = loaded.shape
    border = np.array([[0, 0], [0, lx], [ly, lx], [ly, 0]])
    board = p.CoordinatePolygon(border) + holes
    triangulated = board.triangulate()
    return board, triangulated


def from_numpy_stack_file(filename):
    import polystar as p

    def combine(hs, nx, ny):
        border = p.CoordinatePolygon([[0, 0], [0, nx], [ny, ny], [nx, 0]])
        return border + hs

    loaded = np.load(filename)  #
    nx, ny, time_series = loaded.shape

    bitmaps = [p.BitmapI(o) for o in np.einsum('ijk->kij', loaded)]

    holes = [[poly.inverse for poly in b.extract_image_polygons(1)] for b in bitmaps]

    combined = [combine(h, nx, ny) for h in holes]

    networks = [c.triangulate() for c in combined]

    return {'loaded': loaded, 'bitmap': bitmaps, 'hole': holes, 'combined': combined, 'network': networks}


def clamp(x, minimum=0, maximum=1):
    if x < minimum:
        return minimum
    if x > maximum:
        return maximum
    return x


def local_map_path(map: np.ndarray, x, y, dilate: int = 0, filenames: dict = None):
    import polystar
    nx, ny = map.shape
    # bitmap = polystar.BitmapI(map.copy(order='C')) # copy in case that's the problem?
    bitmap = polystar.BitmapI(map)  # copy in case that's the problem?
    if dilate:
        bitmap = bitmap.dilate(dilate)
    holes = [poly.inverse for poly in bitmap.extract_image_polygons(1)]
    combined = polystar.CoordinatePolygon([[0, 0], [0, ny], [nx, ny], [nx, 0]]) + holes
    network = combined.triangulate()
    # clamp the goal points to be inside the array
    path = network.path([[nx // 2, ny // 2]],
                        [[clamp(x, minimum=1, maximum=nx - 1), clamp(y, minimum=1, maximum=ny - 1)]])
    if np.prod(path.shape) < 4:
        print(f"!!! {path}")
        if 'array' in filenames:
            np.save(filenames['array'].name, map)
        if 'bitmap' in filenames:
            bitmap.write_image(filenames['bitmap'].name)
        if 'combined' in filenames:
            combined.write_svg(filenames['combined'].name, fill='tan')
        if 'network' in filenames:
            network.write_svg(filenames['network'].name, fill='tan')
        if 'path' in filenames:
            np.save(filenames['path'].name, path)

    return path
