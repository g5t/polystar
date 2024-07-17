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


if __name__ == '__main__':
    import unittest
    unittest.main()
