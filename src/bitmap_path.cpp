#include <optional>
#include "bitmap.hpp"

using namespace polystar;
using namespace polystar::bitmap;
using namespace polystar::bitmap::coordinates;
using namespace polystar::bitmap::path;

std::array<std::pair<Point<int>,Point<int>>, 4> cs_neighbours(const coordinates::system cs){
  switch (cs){
    case coordinates::system::y_up_x_right:
      return {{ {{1, -1}, {1, 0}}, {{0, -1}, {1, -1}}, {{-1, 0}, {-1, -1}}, {{1, 1}, {0, 1}} }};
    case coordinates::system::y_down_x_right:
      return {{ {{1, -1}, {1, 0}}, {{0, -1}, {-1, -1}}, {{-1, 0}, {-1, -1}}, {{-1, 1}, {0, 1}} }};
    default:
      throw std::runtime_error("cs_neighbours:: Not implemented");
  }
}
std::array<std::array<Point<int>, 4>, 4> cs_corners(const coordinates::system cs){
  switch (cs){
    case coordinates::system::y_up_x_right:
      return {{
                {{{0, -1}, {0,  0}, {1, -1}, {1,  0}}},
                {{{0,  0}, {1,  0}, {0, -1}, {1, -1}}},
                {{{1,  0}, {1, -1}, {0,  0}, {0, -1}}},
                {{{1, -1}, {0, -1}, {1,  0}, {0,  0}}}
      }};
    case coordinates::system::y_down_x_right:
      return {{
                {{{-1,  0}, {-1, -1}, { 0,  0}, { 0, -1}}},
                {{{ 0,  0}, {-1,  0}, { 0, -1}, {-1, -1}}},
                {{{ 0, -1}, { 0,  0}, {-1, -1}, {-1, 0}}},
                {{{-1, -1}, { 0, -1}, {-1,  0}, { 0,  0}}}
              }};
    default:
      throw std::runtime_error("cs_corners:: Not Implemented");
  }
}
std::array<size_t, 4> cs_left_direction(const coordinates::system cs){
  switch (cs) {
    case coordinates::system::y_up_x_right:
      return {{1, 2, 3, 0}};
    case coordinates::system::y_down_x_right:
      return {{3, 0, 1, 2}};
    default:
      throw std::runtime_error("cs_left_direction:: Not implemented");
  }
}
std::array<size_t, 4> cs_right_direction(const coordinates::system cs){
  switch (cs) {
    case coordinates::system::y_up_x_right:
      return {{3, 0, 1, 2}};
    case coordinates::system::y_down_x_right:
      return {{1, 2, 3, 0}};
    default:
      throw std::runtime_error("cs_left_direction:: Not implemented");
  }
}

template<class T>
Array2<T> flip_range(Array2<T> image, int x0, int y0, int xz){
  if (xz < x0) std::swap(x0, xz);
  for (int i=x0; i<xz; ++i) image.val(y0, i) = !(image.val(y0, i));
  return image;
}

template<class T, class Ptr>
Array2<T> flip_inside(Array2<T> image, Ptr first, Ptr last){
  auto yz = (last - 1)->coord()[0];
  auto xa = first->coord()[1];
  for (; first != last; ++first){
    auto c = first->coord();
    if (c[0] != yz) {
      image = flip_range(image, c[1], std::min(c[0], yz), xa);
      yz = c[0];
    }
  }
  return image;
}


std::tuple<int, std::vector<Point<int>>>
polystar::bitmap::trace_path(Array2<bool> input, const coordinates::system cs, const Point<int> & start) {
auto in_bounds = [&](const auto & p){
  auto c = p.coord();
  return c[0] >= 0 && c[1] >= 0 && c[0] < static_cast<int>(input.size(0)) && c[1] < static_cast<int>(input.size(1));
};
auto is_white = [&](const auto & p){return !in_bounds(p) || !input[as<ind_t>(p).coord()];};
/* We picked the starting point as the index where a transition from white to black has just occurred */
const std::array<Point<int>, 4> directions {{{1, 0}, {0, -1}, {-1, 0}, {0, 1}}};

auto abcd_direction = cs_corners(cs);
auto left = cs_left_direction(cs);
auto right = cs_right_direction(cs);

std::optional<size_t> dir{std::nullopt};
if (is_white(start + Point<int>(0, -1))) {
  // white to left of black -- go along vertical axis first
  dir = 0;
} else if (is_white(start + Point<int>(-1, 0))) {
  // white above black -- go along horizontal axis first
  dir = 3;
}
if (!dir.has_value()) throw std::runtime_error("Unsupported starting edge detection direction!");

std::vector<Point<int>> vertices;
vertices.push_back(start);

auto next_direction = [&](const auto & p, const size_t direction){
  auto abcd = abcd_direction[direction];
  auto a = is_white(p + abcd[0]);
  auto b = is_white(p + abcd[1]);
  auto c = is_white(p + abcd[2]);
  auto d = is_white(p + abcd[3]);
  if (a && b && c && d){
    return std::optional<size_t>(std::nullopt);
  }
  if ((a && !b && !c && !d) || (!a && b && c && d)) {
    return std::make_optional(left[direction]);
  }
  if ((!a && b && !c && !d) || (a && !b && c && d)) {
    return std::make_optional(right[direction]);
  }
  if ((a && !b && c && !d) || (!a && b && !c && d)) {
    return std::make_optional(direction);
  }
  // TODO Implement more complicated choices here?
  return std::make_optional(direction);
};

// calculate path area along the way:
int area{0};

//auto visited = Array2<bool>(input.size(0), input.size(1), false);
do {
  // take a step
  vertices.push_back(vertices.back() + directions[dir.value()]);
  // add to the area
  area += directions[dir.value()].coord()[0] * vertices.back().coord()[1];
  // check for a newly-closed loop and fill in the enclosed area with its opposite color
  if (vertices.size() > 4) { // we need 5 vertices to enclose *1* pixel
  auto ptr = std::find(vertices.begin(), vertices.end() - 1, vertices.back());
  if (std::distance(ptr, vertices.end()) > 4) {
      input = flip_inside(input, ptr, vertices.end());
    }
  }
  dir = next_direction(vertices.back(), dir.value());
} while (dir.has_value() && start != vertices.back());

return std::make_tuple(area, vertices);
}

//
//bool is_cyclic(int a, int b, int c){
//  return a <= c ? (a <= b && b < c) : (a <= b || b<c);
//}
//
//int local_mod(int a, int n){
//  return a >= n ? a % n : a >= 0 ? a : n - 1 - (-1 - a) % n;
//}
//
//int floor_division(int a, int n){
//  return a >= 0 ? a / n : -1 - (-1 - a) / n;
//}
//
//std::vector<int> path_point_longest(const std::vector<Point<int>> & inp){
//  std::vector<int> nc;
//  nc.reserve(inp.size());
//  auto cj = inp.front().coord();
//  size_t j{0};
//  for (size_t i=inp.size(); i-->0; ){
//    auto ci = inp[i].coord();
//    if (ci[0] != cj[0] && ci[1] != cj[1]) {
//      j = i + 1; // i+1 != N since inp[N-1] == inp[0]
//      cj = inp[j].coord();
//    }
//    nc.push_back(j);
//  }
//  std::array<int,4> count{{0, 0, 0, 0}};
//  std::vector<int> pivots;
//  pivots.reserve(inp.size());
//  coordinates::Point<int> encode{1,3}, cont[2]{{0, 0}, {0, 0}};
//
//  int no = static_cast<int>(inp.size());
//
//  for (size_t i=inp.size(); i-->0; ){
//    count[0] = count[1] = count[2] = count[3] = 0;
//
//    count[static_cast<size_t>((3 + encode * (inp[static_cast<size_t>(local_mod(i+1, no))] - inp[i])) / 2)]++;
//
//    cont[0] = coordinates::Point<int>(0, 0);
//    cont[1] = coordinates::Point<int>(0, 0);
//
//    auto j = nc[i];
//    auto j1 = i;
//    while (true) {
//      count[static_cast<size_t>((3 + encode * (inp[j] - inp[j1]).clamp(-1, 1)) / 2)]++;
//      if (count[0] && count[1] && count[2] && count[3]) {
//        pivots.push_back(j1);
//        break;
//      }
//      auto current = inp[j] - inp[i];
//      if (cross(cont[0], current) < 0 || cross(cont[1], current) > 0) {
//        auto deltaj = (inp[j] - inp[j1]).clamp(-1, 1);
//        current = inp[j1] - inp[i];
//        auto a = cross(cont[0], current);
//        auto b = cross(cont[0], deltaj);
//        auto c = cross(cont[1], current);
//        auto d = cross(cont[1], deltaj);
//        j = b < 0 ? floor_division(a, -b) : d > 0 ? floor_division(-c, d) : -100000000;
//        pivots.push_back(local_mod(j1 + j, no));
//        break;
//      }
//      auto ac = abs(current).coord();
//      if (ac[0] > 1 || ac[1] > 1){
//        auto cc = current.coord();
//        auto oy = cc[0] + ((cc[1] <= 0 && (cc[1] < 0 || cc[0] < 0)) ? 1 : -1);
//        auto ox = cc[1] + ((cc[0] >= 0 && (cc[0] > 0 || cc[1] < 0)) ? 1 : -1);
//        Point<int> of0{oy, ox};
//        if (cross(cont[0], of0) >= 0) {
//          cont[0] = of0;
//        }
//        oy = cc[0] + ((cc[1] >= 0 && (cc[1] > 0 || cc[0] < 0)) ? 1 : -1);
//        ox = cc[1] + ((cc[0] <= 0 && (cc[0] < 0 || cc[1] < 0)) ? 1 : -1);
//        Point<int> of1{oy, ox};
//        if (cross(cont[1], of1) <= 0) {
//          cont[1] = of1;
//        }
//      }
//      j1 = j;
//      j = nc[j1];
//      if (!is_cyclic(j, i, j1)) {
//        pivots.push_back(no);
//        break;
//      }
//    }
//  }
//
//  j = pivots.back();
//  std::vector<int> longest(inp.size(), 0);
//  longest.back() = j;
//  for (size_t i=inp.size()-1; i-->0; ){
//    if (is_cyclic(static_cast<int>(i+1), pivots[i], j)) j = pivots[j];
//    longest[i] = j;
//  }
//  for (size_t i=inp.size()-1; is_cyclic(local_mod(static_cast<int>(i+1), no), j, longest[i]); --i){
//    longest[i] = j;
//  }
//
//  return longest;
//}

std::vector<Point<int>> polystar::bitmap::fix_path(const std::vector<Point<int>> & p) {
  // This doesn't produce optimal edges in the potrace sense, but is probably sufficient for now:
  auto stop_segment = [](const auto & s){
    auto l = s.back() - s.front();
    auto pc = as<double>(l).coord();
    auto perp = Point(pc[1], -pc[0]);
    auto pn = perp / std::sqrt(perp * perp);
    for (size_t i=1; i<s.size()-1; ++i){
      if (std::abs((s[i] - s[0]) * pn) >= 1) return true;
    }
    return false;
  };
  std::vector<Point<int>> out;
  out.reserve(p.size());
  out.push_back(p[0]);
  std::vector<Point<int>> segment;
  segment.push_back(p[0]);
  for (size_t i=1; i<p.size(); ++i){
    segment.push_back(p[i]);
    if (stop_segment(segment)){
      out.push_back(p[i-1]);
      segment.clear();
      segment.push_back(p[i-1]);
      segment.push_back(p[i]);
    }
  }
  return out;
}