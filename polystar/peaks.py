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
    z = 3 * (1 - x)**2 * np.exp(-x2-(y+1)**2) - 10*(x/5-x3-y5)*np.exp(-x2-y2) - 1/3*np.exp(-(x+1)**2-y**2)
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
    z[np.logical_and(np.abs(x+y) < 1, np.abs(x-y) < 1)] = 1
    return z.astype(dtype)


def multi(nx, ny, dtype='float'):
    x, y = np.meshgrid(np.linspace(-2, 2, nx), np.linspace(-2, 2, ny), indexing='ij')
    z = 0 * x
    # Diamond at (1,1)
    z[np.logical_and(np.abs((x-1)+(y-1)) < 0.2, np.abs((x-1)-(y-1)) < 0.2)] = 1
    # Square at (-0.5, 0.8)
    z[np.logical_and(np.abs(x+0.5) < 0.3, np.abs(y-0.8) < 0.3)] = 1
    # Circle at (0, -1.1)
    z[(x * x + (y + 1.1)**2) < 0.25] = 1
    return z.astype(dtype)


def board():
    from _polystar import BitmapI, CoordinatePolygon
    bm = BitmapI(multi(300, 300, 'int'))
    polygons = bm.extract_image_polygons(1)
    holes = [p.inverse for p in polygons]
    board = CoordinatePolygon(np.array([[0, 0], [0,300], [300, 300], [300, 0]], dtype='int'))
    return board + holes

def polylineplot(pl, *args, **kwargs):
    import matplotlib.pyplot as pp
    pp.plot(pl[:,0], pl[:,1], *args, **kwargs)

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


