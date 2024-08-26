from unittest import TestCase


class PolygonTestCase(TestCase):

    def test_polygon_creation(self):
        from numpy import array, allclose
        from polystar import Polygon
        square = Polygon([[0, 0], [1, 0], [1, 1], [0, 1]])

        self.assertEqual(square.area, 1.0)

        centroid = square.centroid
        self.assertTrue(allclose(centroid, array([[0.5, 0.5]])))

        self.assertEqual(square.contains(centroid), [True])

    def test_polygon_triangulation(self):
        from polystar import Polygon
        square = Polygon([[0, 0], [1, 0], [1, 1], [0, 1]])
        net = square.triangulate()

        self.assertEqual(len(net.polygons()), 2)
        for polygon in net.polygons():
            self.assertEqual(len(polygon.wires), 0, "The triangles of a square must not have holes")

        self.assertEqual(len(net.wires()), 2)
        for wire in net.wires():
            self.assertEqual(len(wire), 3, "The wire of a triangle has three entries")

    def memory_layout_checks(self, vertices, border, area):
        from polystar import Polygon
        poly = Polygon(vertices)
        poly_area = poly.area
        self.assertAlmostEqual(poly_area, area)
        self.assertAlmostEqual(Polygon(vertices, border).area, area)

    def test_vertices_memory_layout_auto_conversion(self):
        from numpy import array, ascontiguousarray, asfortranarray, allclose
        v_list = [[0, 0.1], [0, 0.001], [1, 0.001], [1, 0.1]]
        border = [0, 1, 2, 3]
        area = (0.1 - 0.001) * (1 - 0)
        # pybind11 converts nested lists as expected
        self.memory_layout_checks(v_list, border, area)

    def test_vertices_memory_layout_numpy_array_conversion(self):
        from numpy import array, ascontiguousarray, asfortranarray, allclose
        v_list = [[0, 0.1], [0, 0.001], [1, 0.001], [1, 0.1]]
        border = [0, 1, 2, 3]
        area = (0.1 - 0.001) * (1 - 0)
        # numpy.array provides C contiguous arrays
        v_array = array(v_list)
        self.assertTrue(v_array.flags['C_CONTIGUOUS'])
        self.memory_layout_checks(v_array, border, area)

    def test_vertices_memory_layout_numpy_ascontiguousarray(self):
        from numpy import array, ascontiguousarray, asfortranarray, allclose
        v_list = [[0, 0.1], [0, 0.001], [1, 0.001], [1, 0.1]]
        border = [0, 1, 2, 3]
        area = (0.1 - 0.001) * (1 - 0)
        # we can also be explicit
        v_c = ascontiguousarray(v_list)
        v_array = array(v_list)
        self.assertTrue(v_c.flags['C_CONTIGUOUS'])
        self.assertTrue(allclose(v_c, v_array))
        self.memory_layout_checks(v_c, border, area)

    def test_vertices_memory_layout_numpy_ascfortranarray(self):
        from numpy import array, ascontiguousarray, asfortranarray, allclose
        v_list = [[0, 0.1], [0, 0.001], [1, 0.001], [1, 0.1]]
        border = [0, 1, 2, 3]
        area = (0.1 - 0.001) * (1 - 0)
        # providing column ordered (FORTRAN style) data should not make any difference
        v_f = asfortranarray(v_list)
        v_array = array(v_list)
        self.assertTrue(v_f.flags['F_CONTIGUOUS'])
        self.assertTrue(allclose(v_f, v_array))
        self.memory_layout_checks(v_f, border, area)

    def test_vertices_memory_layout_numpy_strided(self):
        from numpy import array, ascontiguousarray, asfortranarray, allclose
        v_list = [[0, 0.1], [0, 0.001], [1, 0.001], [1, 0.1]]
        border = [0, 1, 2, 3]
        area = (0.1 - 0.001) * (1 - 0)
        # we should also be able to provide strided arrays
        v_extended = array([[v[0], -1, v[1]] for v in v_list])
        v_stride = v_extended[:, ::2]
        v_array = array(v_list)
        self.assertFalse(v_stride.flags['C_CONTIGUOUS'])
        self.assertFalse(v_stride.flags['F_CONTIGUOUS'])
        self.assertTrue(allclose(v_stride, v_array))
        self.memory_layout_checks(v_stride, border, area)

    def test_found_intersection_error(self):
        from polystar import Polygon
        base = Polygon([[0.   , 0.001], [0.   , 0.013], [0.01 , 0.013], [0.01 , 0.001]])
        reg0 = Polygon([[ 0.02 ,  0.001], [-0.1  ,  0.013], [-0.15 ,  0.013], [-0.03 ,  0.001]])
        reg1 = Polygon([[ 1.2e-01,  1.0e-03], [ 2.77555756e-17,  1.3e-02], [-5.0e-02,  1.3e-02], [ 7.0e-02,  1.0e-03]])

        exp0 = Polygon([[0.01, 0.002], [0.  , 0.003], [0.  , 0.001], [0.01, 0.001]])
        exp1 = Polygon([[0, 0.013], [0, 0.008], [0.01, 0.007], [0.01, 0.012]])

        pb0 = base.intersection(reg0)
        pb1 = base.intersection(reg1)
        p0b = reg0.intersection(base)
        p1b = reg1.intersection(base)

        for p in (pb0, pb1, p0b, p1b):
            self.assertEqual(len(p), 1)

        self.assertAlmostEqual(pb0[0].area, p0b[0].area)
        self.assertAlmostEqual(pb1[0].area, p1b[0].area)
        self.assertAlmostEqual(exp0.area, pb0[0].area)
        self.assertAlmostEqual(exp1.area, pb1[0].area)

        self.assertEqual(pb0[0], p0b[0])
        self.assertEqual(pb1[0], p1b[0])
        self.assertEqual(exp0, pb0[0])
        self.assertEqual(exp1, pb1[0])

    def test_self_intersection(self):
        from polystar import Polygon
        poly = Polygon([(1, 7), (1, 4), (2, 3), (2, 6)], [0, 1, 2, 3])

        self.assertTrue(poly.is_convex)
        self.assertTrue(poly.intersects(poly))

        inter = poly.intersection(poly)
        self.assertEqual(len(inter), 1)
        self.assertEqual(inter[0], poly)

        cinter = poly.convex_intersection(poly)
        self.assertEqual(cinter, poly)


if __name__ == '__main__':
    import unittest
    unittest.main()
