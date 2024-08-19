#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "array_.hpp"

#include "polygon.hpp"

using namespace polystar;
using namespace polystar::polygon;
using namespace Catch::Matchers;

TEST_CASE("Polygon area", "[polygon]"){
  double h{2}, w{3};

  std::vector<std::array<double,2>> va_vertices{{0, 0}, {0, h}, {w, h}, {w, 0}};
  std::vector<ind_t> border{{0, 3, 2, 1}};

  auto vertices = bArray<double>::from_std(va_vertices);
  auto poly = Poly(vertices, static_cast<Wire>(border));

  REQUIRE_THAT(poly.area(), WithinRel(h*w, 1e-12));

  auto hull = Poly(vertices); // w/o border information, the Convex Hull is found
  REQUIRE_THAT(hull.area(), WithinRel(h*w, 1e-12));
  REQUIRE(poly == hull);

  Wire ordered;
  ordered.resize(vertices.size(0));
  std::iota(ordered.begin(), ordered.end(), 0);
  // force clockwise ordering of the vertices
  auto inv_poly = Poly(vertices, ordered);
  REQUIRE_THAT(inv_poly.area(), WithinRel(-h*w, 1e-12));

  std::vector<std::array<double, 2>> va_tri_vertices{{10, 20}, {20, 20}, {15, 30}};
  auto tri_vertices = bArray<double>::from_std(va_tri_vertices);
  auto triangle = Poly(tri_vertices);
  REQUIRE_THAT(triangle.area(), WithinRel(50, 1e-12));
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
    REQUIRE_THAT(poly.area(), WithinRel(12.0, 1e-12));
    REQUIRE_THAT(hull.area(), WithinRel(16.0, 1e-12));
  }

  SECTION("Centroid"){
    auto centroid = poly.centroid();
    REQUIRE_THAT(centroid.val(0, 0), WithinRel(3., 1e-12));
    REQUIRE_THAT(centroid.val(0, 1), WithinRel(3., 1e-12));
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


TEST_CASE("Edge intersection", "[polygon][edge]"){
  std::vector<std::array<double, 2>> va_vertices {
      {1, 1}, {2, 3}, {1, 3}, {4, 1}, {4, 3}, {2, 2}, {5, 0}, {1, 0}, {-2, -1}, {7, 5},
      {-1, 0.25}, {1, 0.25}, {0.4, 0.4}, {0.4, -0.4}
  };
  auto vertices = bArray<double>::from_std(va_vertices);
  SECTION("Mid-edge intersection"){
    auto edge_1 = std::make_pair<ind_t, ind_t>(0, 1);
    auto edge_2 = std::make_pair<ind_t, ind_t>(2, 3);
    REQUIRE(intersect2d(vertices, edge_1, vertices, edge_2));
    auto [flag, at] = intersection2d(vertices, edge_1, vertices, edge_2);
    REQUIRE(flag == 1);
    REQUIRE_THAT(at.val(0, 0), WithinRel(1.75, 1e-12));
    REQUIRE_THAT(at.val(0, 1), WithinRel(2.5, 1e-12));
  }
  SECTION("Error-prone edge intersection"){
    auto edge_1 = std::make_pair<ind_t, ind_t>(10, 11);
    auto edge_2 = std::make_pair<ind_t, ind_t>(12, 13);
    REQUIRE(intersect2d(vertices, edge_1, vertices, edge_2));
    auto [flag, at] = intersection2d(vertices, edge_1, vertices, edge_2);
    REQUIRE(flag == 1);
    REQUIRE_THAT(at.val(0, 0), WithinRel(0.4, 1e-12));
    REQUIRE_THAT(at.val(0, 1), WithinRel(0.25, 1e-12));
  }
  SECTION("Edge intersection at vertex"){
    auto edge_1 = std::make_pair<ind_t, ind_t>(1, 0);
    auto edge_2 = std::make_pair<ind_t, ind_t>(4, 0);
    auto edge_1_reverse = std::make_pair(edge_1.second, edge_1.first);
    auto edge_2_reverse = std::make_pair(edge_2.second, edge_2.first);
    // An opinionated, and perhaps incorrect, choice was to exclude the first vertex of both edges
    // from the possible intersection point.
    REQUIRE(!intersect2d(vertices, edge_1_reverse, vertices, edge_2_reverse));
    // Flipping one edge still doesn't find the intersection since the other first-vertex is excluded
    REQUIRE(!intersect2d(vertices, edge_1_reverse, vertices, edge_2));
    REQUIRE(!intersect2d(vertices, edge_1, vertices, edge_2_reverse));
    // With the common vertex as the second in each edge, they are found to intersect
    REQUIRE(intersect2d(vertices, edge_1, vertices, edge_2));
    auto [flag, at] = intersection2d(vertices, edge_1, vertices, edge_2);
    REQUIRE(flag == 1);
    REQUIRE_THAT(at.val(0, 0), WithinRel(1., 1e-12));
    REQUIRE_THAT(at.val(0, 1), WithinRel(1., 1e-12));
    REQUIRE(at.row_is(cmp::eq, vertices.view(edge_2.second)).all());
  }
  SECTION("Edge no-intersection"){
    // The infinite lines through segments
    // {1, 1} -- {2, 3} and {2, 2} -- {5, 0}
    // intersect, inside the first segment but outside the second.
    auto edge_1 = std::make_pair<ind_t, ind_t>(0, 1);
    auto edge_2 = std::make_pair<ind_t, ind_t>(5, 6);
    REQUIRE(!intersect2d(vertices, edge_1, vertices, edge_2));
    REQUIRE(!intersect2d(vertices, edge_2, vertices, edge_1));
  }
  SECTION("Parallel edges"){
    auto edge_1 = std::make_pair<ind_t, ind_t>(0, 1);
    auto edge_2 = std::make_pair<ind_t, ind_t>(5, 7);
    REQUIRE(!intersect2d(vertices, edge_1, vertices, edge_2));
  }
  SECTION("Colinear edges"){
    auto edge_1 = std::make_pair<ind_t, ind_t>(0, 4);
    auto edge_2 = std::make_pair<ind_t, ind_t>(8, 9);
    REQUIRE(intersect2d(vertices, edge_1, vertices, edge_2));
    REQUIRE(intersect2d(vertices, edge_2, vertices, edge_1));
    auto [flag, at] = intersection2d(vertices, edge_1, vertices, edge_2);
    REQUIRE(flag == 2);
    // at contains both vertices of edge_1
    REQUIRE(at.row_is(cmp::eq, vertices.view(edge_1.first)).sum() == 1);
    REQUIRE(at.row_is(cmp::eq, vertices.view(edge_1.second)).sum() == 1);
  }

}

auto make_poly = [](const std::vector<std::array<double, 2>>& va_vertices){
  auto vertices = bArray<double>::from_std(va_vertices);
  Wire border;
  border.resize(vertices.size(0));
  std::iota(border.begin(), border.end(), 0);
  return Poly(vertices, border);
};

TEST_CASE("Polygon intersections", "[polygon]"){
  std::vector<std::array<double, 2>> va_vertices_1, va_vertices_2, va_vertices_r;
  double area_r{0};
  SECTION("2x2 Square"){
    va_vertices_1 = {{0, 0}, {2, 0}, {2, 2}, {0, 2}};

    SECTION("Diagonal offset, corner intersection"){
      va_vertices_2 = {{1, 1}, {3, 1}, {3, 3}, {1, 3}};
      va_vertices_r = {{1, 1}, {2, 1}, {2, 2}, {1, 2}};
      area_r = 1.0;
    }
    SECTION("Horizontal offset, edge intersection"){
      va_vertices_2 = {{1, 0}, {3, 0}, {3, 2}, {1, 2}};
      va_vertices_r = {{1, 0}, {2, 0}, {2, 2}, {1, 2}};
      area_r = 2.0;
    }
    SECTION("Intersection at vertex is not a valid Polygon"){
      va_vertices_2 = {{2, 2}, {4, 2}, {4, 4}, {2, 4}};
    }
    SECTION("Intersection along one edge is not a valid Polygon"){
      va_vertices_2 = {{2, 1}, {4, 1}, {4, 3}, {2, 3}};
    }
    SECTION("No intersection"){
      va_vertices_2 = {{3, 3}, {5, 3}, {5, 5}, {3, 5}};
    }
    SECTION("Pinch point"){
      va_vertices_1 = {{0, 0}, {2, 0}, {2, 2}, {4, 2}, {4, 4}, {2, 4}, {2, 2}, {0, 2}};
      va_vertices_2 = {{1, 1}, {3, 1}, {3, 3}, {1, 3}};
      va_vertices_r = {{1, 1,}, {2, 1}, {2, 2}, {3, 2}, {3, 3}, {2, 3}, {2, 2}, {1, 2}};
      area_r = 2.0;
    }
  }
  auto poly_1 = make_poly(va_vertices_1);
  auto poly_2 = make_poly(va_vertices_2);
  auto result = polygon_intersection(poly_1, poly_2);

  if (area_r > 0){
    auto poly_r = make_poly(va_vertices_r);
    REQUIRE(result.size() == 1u);
    REQUIRE_THAT(result[0].area(), WithinRel(poly_r.area(), 1e-12));
    REQUIRE_THAT(result[0].area(), WithinRel(area_r, 1e-12));
    REQUIRE(poly_r == result[0].without_extraneous_vertices());
  } else {
    REQUIRE(result.empty());
  }

}


TEST_CASE("Multiple polygon intersections", "[polygon]"){
  auto perform_tests = [](auto va1, auto va2, auto var, auto ar){
    auto p1 = make_poly(va1);
    auto p2 = make_poly(va2);
    auto rs = polygon_intersection(p1, p2);
    if (ar.empty()){
      REQUIRE(rs.empty());
    } else {
      std::vector<Poly<double, polystar::Array2>> polys_r;
      polys_r.reserve(var.size());
      for (const auto& va_vertices : var){
        polys_r.push_back(make_poly(va_vertices));
      }
      // std::cout << rs.size() << " results: \n";
      // for (const auto & r: rs) std::cout << r.without_extraneous_vertices() << "\n\n";
      REQUIRE(rs.size() == var.size());
      for (size_t i=0; i<rs.size(); ++i){
        REQUIRE_THAT(rs[i].area(), WithinRel(polys_r[i].area(), 1e-12));
        REQUIRE_THAT(rs[i].area(), WithinRel(ar[i], 1e-12));
        REQUIRE(polys_r[i] == rs[i].without_extraneous_vertices());
      }
    }
  };

  std::vector<std::array<double, 2>> va_vertices_1, va_vertices_2;
  std::vector<std::vector<std::array<double, 2>>> va_vertices_r;
  std::vector<double> area_r;

  SECTION("Forked"){
    va_vertices_1 = {{0, 0}, {3, 0}, {3, 4}, {0, 4}};
    va_vertices_2 = {{2, 1}, {5, 1}, {5, 3}, {2, 3}, {4, 2}};
    va_vertices_r = {
        {{2, 1}, {3, 1}, {3, 1.5}},
        {{2, 3}, {3, 2.5}, {3, 3}},
    };
    area_r = {0.25, 0.25};
    perform_tests(va_vertices_1, va_vertices_2, va_vertices_r, area_r);
  }
  SECTION("Inverted big G and rectangle"){
    va_vertices_1 = {{5, 5}, {8, 5}, {8, 3}, {7, 1}, {3, 1}, {2, 3}, {2, 8},
                     {3, 9}, {7, 9}, {8, 8}, {8, 7}, {7, 7}, {6, 8}, {4, 8},
                     {3, 7}, {3, 3}, {4, 2}, {6, 2}, {7, 3}, {7, 4}, {5, 4}};
    va_vertices_2 = {{4, 3}, {6, 3}, {6, 10}, {4, 10}};
    va_vertices_r = {{{4, 9}, {6, 9}, {6, 10}, {4, 10}},
                     {{6, 8}, {4, 8}, {4, 3}, {6, 3}, {6, 4}, {5 ,4}, {5, 5}, {6, 5}}};
    area_r = {2.0, 9.0};
    perform_tests(va_vertices_1, va_vertices_2, va_vertices_r, area_r);
  }
  SECTION("Big G and rectangle"){
    va_vertices_1 = {{5, 5}, {5, 4}, {7, 4}, {7, 3}, {6, 2}, {4, 2}, {3, 3},
                     {3, 7}, {4, 8}, {6, 8}, {7, 7}, {8, 7}, {8, 8}, {7, 9},
                     {3, 9}, {2, 8}, {2, 3}, {3, 1}, {7, 1}, {8, 3}, {8, 5}};
    va_vertices_2 = {{4, 3}, {6, 3}, {6, 10}, {4, 10}};
    va_vertices_r = {{{4, 8}, {6, 8}, {6, 9}, {4, 9}},
                     {{6, 5}, {5, 5}, {5, 4}, {6, 4}}};
    area_r = {2.0, 1.0};
    perform_tests(va_vertices_1, va_vertices_2, va_vertices_r, area_r);
  }

}