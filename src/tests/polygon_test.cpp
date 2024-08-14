#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "array_.hpp"

#include "polygon.hpp"

using namespace polystar;
using namespace polystar::polygon;

TEST_CASE("Polygon area", "[polygon]"){
  double h{2}, w{3};

  std::vector<std::array<double,2>> va_vertices{{0, 0}, {0, h}, {w, h}, {w, 0}};
  std::vector<ind_t> border{{0, 3, 2, 1}};

  auto vertices = bArray<double>::from_std(va_vertices);
  auto poly = Poly(vertices, static_cast<Wire>(border));

  REQUIRE_THAT(poly.area(), Catch::Matchers::WithinRel(h*w, 1e-12));

  auto hull = Poly(vertices); // w/o border information, the Convex Hull is found
  REQUIRE_THAT(hull.area(), Catch::Matchers::WithinRel(h*w, 1e-12));
  REQUIRE(poly == hull);

  Wire ordered;
  ordered.resize(vertices.size(0));
  std::iota(ordered.begin(), ordered.end(), 0);
  // force clockwise ordering of the vertices
  auto inv_poly = Poly(vertices, ordered);
  REQUIRE_THAT(inv_poly.area(), Catch::Matchers::WithinRel(-h*w, 1e-12));

  std::vector<std::array<double, 2>> va_tri_vertices{{10, 20}, {20, 20}, {15, 30}};
  auto tri_vertices = bArray<double>::from_std(va_tri_vertices);
  auto triangle = Poly(tri_vertices);
  REQUIRE_THAT(triangle.area(), Catch::Matchers::WithinRel(50, 1e-12));
}


TEST_CASE("Non-convex Polygon area", "[polygon]"){
  std::vector<std::array<double, 2>> va_vertices{
      {1, 1}, {2, 1}, {3, 2}, {4, 1},
      {5, 1}, {5, 2}, {4, 3}, {5, 4},
      {5, 5}, {4, 5}, {3, 4}, {2, 5},
      {1, 5}, {1, 4}, {2, 3}, {1, 2}
  };
  auto vertices = bArray<double>::from_std(va_vertices);

  Wire ordered;
  ordered.resize(vertices.size(0));
  std::iota(ordered.begin(), ordered.end(), 0);

  auto poly = Poly(vertices, ordered);
  auto hull = Poly(vertices);

  SECTION("Area"){
    REQUIRE_THAT(poly.area(), Catch::Matchers::WithinRel(12.0, 1e-12));
    REQUIRE_THAT(hull.area(), Catch::Matchers::WithinRel(16.0, 1e-12));
  }

  SECTION("Centroid"){
    auto centroid = poly.centroid();
    REQUIRE_THAT(centroid.val(0, 0), Catch::Matchers::WithinRel(3., 1e-12));
    REQUIRE_THAT(centroid.val(0, 1), Catch::Matchers::WithinRel(3., 1e-12));
  }

  SECTION("Intersection"){
    auto overlap = polygon_intersection(poly, hull);
    REQUIRE(overlap.size() == 1u);
    REQUIRE(overlap[0] == poly);

    auto o1 = poly.intersection(hull);
    auto o2 = hull.intersection(poly);
    REQUIRE(o1.size() == 1u);
    REQUIRE(o2.size() == 1u);
    REQUIRE(o1[0] == o2[0]);
    REQUIRE(overlap[0] == o1[0]);
  }
}