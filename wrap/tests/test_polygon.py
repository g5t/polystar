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
