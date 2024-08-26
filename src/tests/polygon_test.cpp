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


auto make_poly = [](const std::vector<std::array<double, 2>>& va_vertices){
  auto vertices = bArray<double>::from_std(va_vertices);
  Wire border;
  border.resize(vertices.size(0));
  std::iota(border.begin(), border.end(), 0);
  return Poly(vertices, border);
};

TEST_CASE("Polygon intersections", "[polygon]"){
  auto do_tests = [](const auto & name, const auto & va1, const auto & va2, const auto & var, const auto var_area){
    auto poly_1 = make_poly(va1);
    auto poly_2 = make_poly(va2);
    auto result = polygon_intersection(poly_1, poly_2);

    if (var_area > 0){
      auto poly_r = make_poly(var);
      REQUIRE(result.size() == 1u);
      auto r = result[0].without_extraneous_vertices();
      if (std::abs(r.area() - poly_r.area()) > 1e-12){
        std::cout << "Test " << name << "area mismatch \nresult:\n" << r << "\ndiffers from expected:\n" << poly_r << "\n";
      }
      REQUIRE_THAT(r.area(), WithinRel(poly_r.area(), 1e-12));
      REQUIRE_THAT(r.area(), WithinRel(var_area, 1e-12));
      REQUIRE(poly_r == r.without_extraneous_vertices());
    } else {
      REQUIRE(result.empty());
    }
  };

  auto do_convex_tests = [](const auto & va1, const auto & va2, const auto & var, const auto var_area) {
    auto poly_1 = make_poly(va1);
    auto poly_2 = make_poly(va2);
    // all polygons are convex, so we can use the other algorithm
    REQUIRE(poly_1.is_convex());
    REQUIRE(poly_2.is_convex());
    auto result = polygon_convex_intersection(poly_1, poly_2);

    REQUIRE_THAT(result.area(), WithinRel(var_area, 1e-12));
    if (var_area > 0){
      auto poly_r = make_poly(var);
      REQUIRE_THAT(result.area(), WithinRel(poly_r.area(), 1e-12));
      REQUIRE(poly_r == result.without_extraneous_vertices());
    }
  };

  std::vector<std::array<double, 2>> va_vertices_1, va_vertices_2, va_vertices_r;
  double area_r{0};
  SECTION("2x2 Square"){
    va_vertices_1 = {{0, 0}, {2, 0}, {2, 2}, {0, 2}};

    SECTION("Diagonal offset, corner intersection"){
      va_vertices_2 = {{1, 1}, {3, 1}, {3, 3}, {1, 3}};
      va_vertices_r = {{1, 1}, {2, 1}, {2, 2}, {1, 2}};
      area_r = 1.0;
      do_tests("Diagonal offset, corner intersection", va_vertices_1, va_vertices_2, va_vertices_r, area_r);
      do_convex_tests(va_vertices_1, va_vertices_2, va_vertices_r, area_r);
    }
    SECTION("Horizontal offset, edge intersection"){
      va_vertices_2 = {{1, 0}, {3, 0}, {3, 2}, {1, 2}};
      va_vertices_r = {{1, 0}, {2, 0}, {2, 2}, {1, 2}};
      area_r = 2.0;
      do_tests("Horizontal offset, edge intersection", va_vertices_1, va_vertices_2, va_vertices_r, area_r);
      do_convex_tests(va_vertices_1, va_vertices_2, va_vertices_r, area_r);
    }
    SECTION("Intersection at vertex is not a valid Polygon"){
      va_vertices_2 = {{2, 2}, {4, 2}, {4, 4}, {2, 4}};
      do_tests("Intersection at vertex is not a valid Polygon", va_vertices_1, va_vertices_2, va_vertices_r, area_r);
      do_convex_tests(va_vertices_1, va_vertices_2, va_vertices_r, area_r);
    }
    SECTION("Intersection along one edge is not a valid Polygon"){
      va_vertices_2 = {{2, 1}, {4, 1}, {4, 3}, {2, 3}};
      do_tests("Intersection along one edge is not a valid Polygon", va_vertices_1, va_vertices_2, va_vertices_r, area_r);
      do_convex_tests(va_vertices_1, va_vertices_2, va_vertices_r, area_r);
    }
    SECTION("No intersection"){
      va_vertices_2 = {{3, 3}, {5, 3}, {5, 5}, {3, 5}};
      do_tests("No intersection", va_vertices_1, va_vertices_2, va_vertices_r, area_r);
      do_convex_tests(va_vertices_1, va_vertices_2, va_vertices_r, area_r);
    }
    SECTION("Pinch point"){
      va_vertices_1 = {{0, 0}, {2, 0}, {2, 2}, {4, 2}, {4, 4}, {2, 4}, {2, 2}, {0, 2}};
      va_vertices_2 = {{1, 1}, {3, 1}, {3, 3}, {1, 3}};
      va_vertices_r = {{1, 1,}, {2, 1}, {2, 2}, {3, 2}, {3, 3}, {2, 3}, {2, 2}, {1, 2}};
      area_r = 2.0;
      do_tests("Pinch point", va_vertices_1, va_vertices_2, va_vertices_r, area_r);
      // no convex tests since the pinch point makes it non-convex
    }
  }
  SECTION("1x5 Rectangle"){
    va_vertices_1 = {{3, 0}, {4, 0}, {4, 5}, {3, 5}};
    SECTION("Parallelogram top edge-vertex intersection"){
      va_vertices_2 = {{5, 0}, {9, 0}, {4, 5}, {0, 5}};
      va_vertices_r = {{3, 5}, {3, 2}, {4, 1}, {4, 5}};
      area_r = 3.5;
      do_tests("Parallelogram top edge-vertex intersection", va_vertices_1, va_vertices_2, va_vertices_r, area_r);
      do_convex_tests(va_vertices_1, va_vertices_2, va_vertices_r, area_r);
    }
    SECTION("Parallelogram bottom edge-vertex intersection"){
      va_vertices_2 = {{1, 0}, {5, 0}, {0, 5}, {-4, 5}};
      va_vertices_r = {{3, 2}, {3, 0}, {4, 0}, {4, 1}};
      area_r = 1.5;
      do_tests("Parallelogram bottom edge-vertex intersection", va_vertices_1, va_vertices_2, va_vertices_r, area_r);
      do_convex_tests(va_vertices_1, va_vertices_2, va_vertices_r, area_r);
    }
    SECTION("Right parallelogram bottom edge-vertex intersection"){
      va_vertices_2 = {{0, 0}, {4, 0}, {9, 5}, {5, 5}};
      va_vertices_r = {{3, 3}, {3, 0}, {4, 0}, {4, 4}};
      area_r = 3.5;
      do_tests("Right parallelogram bottom edge-vertex intersection", va_vertices_1, va_vertices_2, va_vertices_r, area_r);
      do_convex_tests(va_vertices_1, va_vertices_2, va_vertices_r, area_r);
    }
  }
  SECTION("1x3 Parallelogram"){
    va_vertices_1 = {{1, 7}, {1, 4}, {2, 3}, {2, 6}};
    SECTION("Self intersection"){
      area_r = 3.0;
      do_tests("1x3 Parallelogram self intersection", va_vertices_1, va_vertices_1, va_vertices_1, area_r);
      do_convex_tests(va_vertices_1, va_vertices_1, va_vertices_1, area_r);
    }
    SECTION("Similar lozenge"){
      va_vertices_2 = {{1.5, 7.}, {1., 7.}, {1., 4}, {1.5, 3.}, {2., 3.}, {2., 6.}};
      area_r = 3.0;
      do_tests("1x3 Parallelogram with larger polygon",
               va_vertices_1, va_vertices_2, va_vertices_1, area_r);
      do_convex_tests(va_vertices_1, va_vertices_2, va_vertices_1, area_r);
    }
    SECTION("Similar parallelogram"){
      va_vertices_2 = {{1.00000095, 7.}, {1., 7.}, {1., 3.99999905}, {1.99999905, 3.}, {2., 3.}, {2., 6.00000095}};
      area_r = 3.0;
      do_tests("1x3 Parallelogram with slightly larger similar parallelogram",
               va_vertices_1, va_vertices_2, va_vertices_1, area_r);
      do_convex_tests(va_vertices_1, va_vertices_2, va_vertices_1, area_r);
    }
  }

}

TEST_CASE("Polygon intersection reciprocity", "[polygon][reciprocity]"){

  auto do_test = [](const auto &v1, const auto &vA, const auto &vR){
    auto poly_1 = make_poly(v1);
    auto poly_A = make_poly(vA);

    REQUIRE(poly_1.area() > 0);
    REQUIRE(poly_A.area() > 0);
//    std::cout << name << " Poly1:\n" << poly_1 << "\n" << name << " PolyA:\n" << poly_A << "\n";

    auto int_1A = polygon_intersection(poly_1, poly_A);
    auto int_A1 = polygon_intersection(poly_A, poly_1);

    REQUIRE(int_1A.size() == int_A1.size());
    for (size_t i=0; i < int_1A.size(); ++i){
      auto p1A = int_1A[i].without_extraneous_vertices();
      auto pA1 = int_A1[i].without_extraneous_vertices();
//      std::cout << name << " " << i << " 1A:\n" << p1A << "\n" << name << " " << i << " A1:\n" << pA1 << "\n";
      REQUIRE_THAT(p1A.area(), WithinRel(p1A.area(), 1e-14));
      REQUIRE(p1A == pA1);
    }
    if (!vR.empty()){
      REQUIRE(vR.size() == int_1A.size());
      for (size_t i = 0; i < vR.size(); ++i){
        auto poly_R = make_poly(vR[i]);
        REQUIRE_THAT(int_1A[i].area(), WithinRel(poly_R.area(), 1e-12));
        REQUIRE(int_1A[i].without_extraneous_vertices() == poly_R);
      }
    }
  };
  auto do_convex_tests = [](const auto & v1, const auto & vA, const auto & vR) {
    auto poly_1 = make_poly(v1);
    auto poly_A = make_poly(vA);
    // all polygons are convex, so we can use the other algorithm
    REQUIRE(poly_1.is_convex());
    REQUIRE(poly_A.is_convex());
    auto r_1A = polygon_convex_intersection(poly_1, poly_A);
    auto r_A1 = polygon_convex_intersection(poly_A, poly_1);
    REQUIRE_THAT(r_1A.area(), WithinRel(r_A1.area(), 1e-12));
    REQUIRE(r_1A == r_A1);

    auto poly_r = make_poly(vR);
    REQUIRE_THAT(r_1A.area(), WithinRel(poly_r.area(), 1e-12));
    REQUIRE(poly_r == r_1A);
  };


  std::vector<std::array<double, 2>> va1, vaA;
  std::vector<std::vector<std::array<double, 2>>> vaR;
  SECTION("Rectangle and parallelogram") {
    va1 = {{0,    0.013},
           {0,    0.001},
           {0.01, 0.001},
           {0.01, 0.013}};
    SECTION("On-border points, lower") {
      vaA = {{0.02,  0.001},
             {-0.1,  0.013},
             {-0.15, 0.013},
             {-0.03, 0.001}};
      vaR = {
          {{0.01, 0.002}, {0, 0.003}, {0, 0.001}, {0.01, 0.001}}
      };
      do_test(va1, vaA, vaR);
      do_convex_tests(va1, vaA, vaR[0]);
    }
    SECTION("On-border points, upper") {
      vaA = {{-0.05,          0.013},
             {0.012,          0.001},
             {0.07,           0.001},
             {2.77555756E-17, 0.013},};
      vaR = {
          {{0, 0.013}, {0, 0.003322580645161289}, {0.01, 0.0013870967741935471}, {0.01, 0.011285714285714291}}
      };
      do_test(va1, vaA, vaR);
      do_convex_tests(va1, vaA, vaR[0]);
    }
  }
}

//TEST_CASE("Self-intersecting Polygon intersection reciprocity", "[polygon]"){
//
//  auto do_test = [](const auto & name, const auto &v1, const auto &vA, const auto &vR){
//    auto poly_1 = make_poly(v1);
//    auto poly_A = make_poly(vA);
//
//    std::cout << name << " Poly1:\n" << poly_1 << "\n" << name << " PolyA:\n" << poly_A << "\n";
//
//    auto int_1A = polygon_intersection(poly_1, poly_A);
//    auto int_A1 = polygon_intersection(poly_A, poly_1);
//
//    REQUIRE(int_1A.size() == int_A1.size());
//    for (size_t i=0; i < int_1A.size(); ++i){
//      REQUIRE_THAT(int_1A[i].area(), WithinRel(int_A1[i].area(), 1e-14));
//      auto p1A = int_1A[i].without_extraneous_vertices();
//      auto pA1 = int_A1[i].without_extraneous_vertices();
//      REQUIRE(p1A == pA1);
//      std::cout << name << " " << i << " 1A:\n" << p1A << "\n" << name << " " << i << " A1:\n" << pA1 << "\n";
//    }
//
//    if (!vR.empty()){
//      REQUIRE(vR.size() == int_1A.size());
//      for (size_t i = 0; i < vR.size(); ++i){
//        auto poly_R = make_poly(vR[i]);
//        REQUIRE_THAT(int_1A[i].area(), WithinRel(poly_R.area(), 1e-12));
//        REQUIRE(int_1A[i].without_extraneous_vertices() == poly_R);
//      }
//
//    }
//
//  };
//
//  std::vector<std::array<double, 2>> va1, vaA;
//  std::vector<std::vector<std::array<double, 2>>> vaR;
//  SECTION("Hour glass, lower"){
//  va1 = {{0, 0.001}, {0, 0.013}, {0.01, 0.013}, {0.01, 0.001}};
//  vaA = {{0.02, 0.001}, {-0.15, 0.013}, {-0.1, 0.013}, {-0.03, 0.001}};
//  vaR = {{{0.01, 0.002}, {0, 0.003}, {0, 0.001}, {0.01, 0.001}}}; // wrong
//  do_test("lower", va1, vaA, vaR); // error thrown because areas of intersections don't match
//  }
//  SECTION("Hour glass, upper"){
//  va1 = {{0, 0.001}, {0, 0.013}, {0.01, 0.013}, {0.01, 0.001}};
//  vaA = {{0.012, 0.001}, {2.77555756E-17, 0.013}, {-0.05, 0.013}, {0.07, 0.001}};
//  vaR = {{{0, 0.013}, {0, 0.008}, {0.01, 0.007}, {0.01, 0.012}}}; // wrong
//  do_test("upper", va1, vaA, vaR);
//  }
//}


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