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
  /* For the corner of a pixel (top left in y_down_x_right)
   * Return the movement required to reach the pixel data
   *    {behind and to the left, behind and to the right, ahead and to the left, ahead and to the right}
   * for each facing direction
   *    {positive y, negative x, negative y, positive x}
   *
   * Note, the coordinates are all (y, x) :/
   * */
  switch (cs){
    case coordinates::system::y_up_x_right:
      return {{
                {{{-1, -1}, {-1,  0}, { 0, -1}, { 0,  0}}},
                {{{-1,  0}, { 0,  0}, {-1, -1}, { 0, -1}}},
                {{{ 0,  0}, { 0, -1}, {-1,  0}, {-1, -1}}},
                {{{ 0, -1}, {-1, -1}, { 0,  0}, {-1,  0}}}
              }};
    case coordinates::system::y_down_x_right:
      return {{// back left, back right, ahead left, ahead right
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

// Along with the definitions of the coordinate systems below, this is (surprisingly) independent
// of *which* way positive y faces
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

// The next direction chosen depends on whether the coordinate system is right- or left-handed
// These rules appear to only work for left-handed coordinate systems, strangely.
// A (temporary?) fix for right-handed coordinate systems is to flip the left and right indexing used below
auto next_direction = [&](const auto & p, const size_t direction){
  auto abcd = abcd_direction[direction];
  auto a = is_white(p + abcd[0]); // back and to the left
  auto b = is_white(p + abcd[1]); // back and to the right
  auto c = is_white(p + abcd[2]); // forward and to the left
  auto d = is_white(p + abcd[3]); // forward and to the right
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
  vertices.push_back(vertices.back() + directions[*dir]);
  // add to the area
  area += directions[*dir].coord()[0] * vertices.back().coord()[1];
  // check for a newly-closed loop and fill in the enclosed area with its opposite color
  if (vertices.size() > 4) { // we need 5 vertices to enclose *1* pixel
  auto ptr = std::find(vertices.begin(), vertices.end() - 1, vertices.back());
  if (std::distance(ptr, vertices.end()) > 4) {
    input = flip_inside(input, ptr, vertices.end());
    }
  }
  dir = next_direction(vertices.back(), *dir);
} while (dir.has_value() && start != vertices.back());

return std::make_tuple(area, vertices);
}

std::vector<Point<int>> polystar::bitmap::fix_path(const std::vector<Point<int>> & p) {
  // This doesn't produce optimal edges in the potrace sense, but is probably sufficient for now:
  auto stop_segment = [](const auto & s){
    auto l = s.back() - s.front();
    auto pc = as<double>(l).coord();
    auto perp = Point(pc[1], -pc[0]);
    auto pn = perp / std::sqrt(perp * perp);
    // force horizontal or vertical lines once we exceed some length?
    if (s.size() > 20) for (size_t i=19; i<s.size()-1; ++i) if (std::abs((s[i] - s[0]) * pn) > 0.5) return true;
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