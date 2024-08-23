#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "array_.hpp"

#include "polygon.hpp"

using namespace polystar;
using namespace polystar::polygon;
using namespace Catch::Matchers;


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
  SECTION("Collinear edges"){
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


TEST_CASE("Collinear edge intersection reciprocity", "[polygon][edge][reciprocity]"){
  auto run = [](const auto & var, const auto & ea, const auto & eb){
    auto v = bArray<double>::from_std(var);
    REQUIRE(intersect2d(v, ea, v, eb));
    REQUIRE(intersect2d(v, eb, v, ea));
    auto [v1, a1] = intersection2d(v, ea, v, eb);
    auto [v2, a2] = intersection2d(v, eb, v, ea);
    REQUIRE(v1 == v2);
    REQUIRE(v1 == 2);
    REQUIRE(a1.row_is(cmp::eq, v.view(eb.first)).sum() == 1);
    REQUIRE(a1.row_is(cmp::eq, v.view(eb.second)).sum() == 1);
    REQUIRE(a2.row_is(cmp::eq, v.view(eb.first)).sum() == 1);
    REQUIRE(a2.row_is(cmp::eq, v.view(eb.second)).sum() == 1);
  };
  std::vector<std::array<double ,2>> va;
  std::pair<std::pair<ind_t, ind_t>, std::pair<ind_t, ind_t>> edges;
  SECTION("Polygon intersection reciprocity error"){
    // this is the case found when testing polygon intersection reciprocity
    va = {{-0.03, 0.001}, {0, 0.001}, {0.01, 0.001}, {0.02, 0.001}};
    edges = {{0, 3}, {1, 2}};
    run(va, edges.first, edges.second);
    edges.second = {2, 1};
    run(va, edges.first, edges.second);

  }
}