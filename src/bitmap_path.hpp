#ifndef POLYSTAR_BITMAP_PATH_HPP
#define POLYSTAR_BITMAP_PATH_HPP

#include <vector>
#include "bitmap_coord.hpp"

namespace polystar::bitmap::path {
  template<class T>
  class PathPartialSums {
  public:
    using value_t = std::vector<T>;
  private:
    value_t x_;
    value_t y_;
    value_t xx_;
    value_t xy_;
    value_t yy_;
  public:
    explicit PathPartialSums(const std::vector<coordinates::Point<T>> & inp) {
      auto c0 = inp.front().coord();
      reserve(inp.size());
      for (const auto &p: inp) {
        auto c = p.coord();
        auto cx = c[1] - c0[1];
        auto cy = c[0] - c0[0];
        extend_back(cx, cy);
      }
    }
    void reserve(size_t n){
      x_.reserve(n);
      y_.reserve(n);
      xx_.reserve(n);
      xy_.reserve(n);
      yy_.reserve(n);
    }
    void extend_back(T x, T y){
      T xb{0}, yb{0}, xxb{0}, xyb{0}, yyb{0};
      if (!x_.empty()){
        xb = x_.back();
        yb = y_.back();
        xxb = xx_.back();
        xyb = xy_.back();
        yyb = yy_.back();
      }
      x_.push_back(xb + x);
      y_.push_back(yb + y);
      xx_.push_back(xxb + x * x);
      xy_.push_back(xyb + x * y);
      yy_.push_back(yyb + y * y);
    }
    template<class I>
    double penalty(const std::vector<coordinates::Point<I>> & p, size_t a, size_t b) const {
      auto beyond = b >= p.size();
      if (beyond) {
        b -= p.size();
      }
      auto jx = x_[b+1] - x_[a];
      auto jy = y_[b+1] - y_[a];
      auto jxx = xx_[b+1] - xx_[a];
      auto jxy = xy_[b+1] - xy_[a];
      auto jyy = yy_[b+1] - yy_[a];
      auto jk = b + 1 - a;
      if (beyond) {
        jx += x_[p.size()];
        jy += y_[p.size()];
        jxx += xx_[p.size()];
        jxy += xy_[p.size()];
        jyy += yy_[p.size()];
        jk += p.size();
      }
      auto v = (coordinates::as<double>(p[a] + p[b]) * 0.5 - coordinates::as<double>(p[0])).coord();
      auto e = (coordinates::as<double>(p[b] - p[a])).coord();
      auto cxx = ((jxx - 2 * jx * v[1]) / jk + v[1] * v[1]);
      auto cxy = ((jxy - jx * v[0] - jy * v[1]) / jk + v[0] * v[1]);
      auto cyy = ((jyy - 2 * jy * v[0]) / jk + v[0] * v[0]);
      return std::sqrt(e[1] * e[1] * cxx - 2 * e[1] * e[0] * cxy + e[0] * e[0] * cyy);
    }

  };

//  /* Trace a constant-value region of a bitmap to extract its border as a path.
//   * With inspiration from the description of the potrace algorithm.
//   * https://potrace.sourceforge.net/potrace.pdf
//   * */
//
//  /*
//   *   0 1 2 3 4
//   *    0 1 2 3
//   * 0 +-+-+-+-+
//   *  0|w|w|w|w|
//   * 1 +-+-+-+-+
//   *  1|w|b|b|b|
//   * 2 +-+-+-+-+
//   *  2|b|b|b|w|
//   * 3 +-+-+-+-+
//   *
//   * */
//
//  template<class T> std::tuple<int, std::vector<coordinates::Point<int>>>
//  trace_path(Bitmap<T> input, const coordinates::system cs,
//             const coordinates::Point<int> & start,
//             const T black, const T white, const T tol, const int dig) {
//    using namespace approx_float;
//    using namespace coordinates;
//    auto is_white = [&](const auto & p){return !image.in_bounds(p) || scalar(input.value(p), white, tol, tol, dig); };
//    // from our starting point find the orientation which gives black on left and white on right
//    /* We picked the starting point as the index where a transition from white to black has just occurred */
//    const std::array<Point<int>, 4> directions {{1, 0}, {0, -1}, {-1, 0}, {0, 1}};
////    const std::array<Point<int>, 4> reverse {{-1, 0}, {0, 1}, {1, 0}, {0, -1}};
//    const std::array<std::array<size_t, 3>, 4> next_directions {{1,2,3}, {0,2,3}, {0,1,3}, {0,1,2}};
//    const std::array<std::pair<Point<int>,Point<int>>, 4> neighbors {
//      {{1, -1}, {1, 0}}, {{0, -1}, {1, -1}}, {{-1, 0}, {-1, -1}}, {{1, 1}, {1, 0}}
//    };
//    size_t dir = 4;
//    if (is_white(start + Point<int>(0, -1))) {
//      // white to left of black -- go along vertical axis first
//      dir = 0;
//    } else if (is_white(start + Point<int>(-1, 0))) {
//      // white above black -- go along horizontal axis first
//      dir = 3;
//    }
//    if (dir >= 4) throw std::runtime_error("Unsupported starting edge detection direction!");
//
//    std::vector<Point<int>> vertices;
//    vertices.push_back(start);
//
//    auto is_possible = [&](const size_t i){
//      return is_white(vertices.back() + neighbors[i].first) ^ is_white(vertices.back() + neighbors[i].second);
//    };
//
//    auto next_possible = [&](const auto & p, const size_t d){
//      std::vector<size_t> np;
//      np.reserve(3);
//      for (const auto & i: next_directions[d]) {
////      for (size_t i=0; i<4; ++i) if (directions[i] != reverse[d]) {
//        // this direction is only possible if *one* neighbour is white
//        if (is_white(p + neighbors[i].first) ^ is_white(p + neighbors[i].second)) np.push_back(i);
//      }
//      return np;
//    };
//
//    auto flip_range = [&](int x0, int y0, int xz){
//      if (xz > xa) std::swap(xa, xz);
//      for (int i=x0; i<xz; ++i) input.value(y0, x0) = is_white(Point<int>(y0, x0)) ? black : white;
//    };
//
//    auto flip_inside = [&](auto first, const auto last){
//      auto yz = *(last - 1).coord()[0];
//      auto xa = *first.coord()[1];
//      for (; first != last; ++first){
//        auto c = *first.coord();
//        if (c[0] != yz) {
//          flip_range(c[1], std::min(c[0], yz), xa);
//          yz = c[0];
//        }
//      }
//    };
//
//    // calculate path area along the way:
//    int area{0};
//
//    do {
//      // take a step
//      vertices.push_back(vertices.back() + directions[dir]);
//      // add to the area
//      area += directions[dir].coord()[0] * vertices.back().coord()[1];
//      // check for a newly-closed loop and fill in the enclosed area with its opposite color
//      if (vertices.size() > 4) { // we need 5 vertices to enclose *1* pixel
//        auto ptr = std::find(vertices.begin(), vertices.end() - 1, vertices.back());
//        if (std::distance(ptr, vertices.end()) > 4) {
//          /* Since the Bitmap class uses an Array2 for storage, its imperative that the data it contains is not needed
//           * for any other operation, as overwriting the 'copy' of Bitmap here will overwrite the original's Array2
//           * container. (At some point the Bitmap.map_ should have been `.decouple()`ed.) */
//          flip_inside(ptr, vertices.end());
//        }
//      }
//      // decide whether we turn left, go straight, or turn right based on which direction give left black right white
//      auto options = next_possible(vertices.back(), dir);
////      std::vector<size_t> options;
////      std::transform(next_directions[d].begin(), next_directions[d].end(), std::back_inserter(options), is_possible);
//      if (options.empty()) break;
//      // always take the first option even if there is more tha one?
//      dir = options[0];
//    } while (start != vertices.back());
//
//    return std::make_tuple(area, vertices);
//  }


}

#endif