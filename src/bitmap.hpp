#ifndef POLYSTAR_BITMAP_HPP
#define POLYSTAR_BITMAP_HPP

#include <vector>
#include <deque>
#include <functional>
#include <optional>

#include "array_.hpp"
#include "polygon.hpp"
#include "bitmap_color.hpp"
#include "bitmap_coord.hpp"
#include "bitmap_path.hpp"
#include "bitmap_show.hpp"

namespace polystar::bitmap {
  template<class T>
  class Bitmap {
  public:
    using map_t = Array2<T>;
  private:
    map_t map_;

  public:
    explicit Bitmap(map_t m): map_(m) {}
    explicit Bitmap(std::vector<std::vector<T>> m): map_(m) {}

    Bitmap(ind_t rows, ind_t cols, T value=T(0)){
      map_ = Array<T>(rows, cols, value);
    }

    ind_t rows() const {return map_.size(0u);}
    ind_t cols() const {return map_.size(1u);}
    const T & value(ind_t row, ind_t col) const {return map_.val(row, col);}
    T & value(ind_t row, ind_t col) {return map_.val(row, col);}
    template<class I> T value(const std::array<I, 2>& c) const {
      return map_.val(c);
    }
    template<class I> T value(const coordinates::Point<I> & p) const {
      return map_.val(coordinates::as<ind_t>(p).coord());
    }
    template<class I> bool in_bounds(const coordinates::Point<I> & p) const {
      auto ci = p.coord();
      return ci[0] >= 0 && ci[0] < static_cast<I>(map_.size(0)) && ci[1] >= 0 && ci[1] < static_cast<I>(map_.size(1));
    }

    std::vector<std::vector<Color>> image(coordinates::system cs = coordinates::system::y_up_x_right) const {
      return image(colormaps::cividis, cs);
    }
    std::vector<std::vector<Color>> image(const std::map<T, Color>& map, coordinates::system cs = coordinates::system::y_up_x_right) const {
      std::vector<std::vector<Color>> image;
      image.reserve(rows());
      switch(cs){
        case coordinates::system::y_up_x_right: {
          for (ind_t i=rows(); i-->0;){
            std::vector<Color> row;
            row.reserve(cols());
            for (ind_t j=0; j<cols(); ++j) row.push_back(map[map_.val(i,j)]);
            image.push_back(row);
          }
        }
          break;
        case coordinates::system::y_down_x_right: {
          for (ind_t i=0; i<rows(); ++i){
            std::vector<Color> row;
            row.reserve(cols());
            for (ind_t j=0; j<cols(); ++j) row.push_back(map[map_.val(i,j)]);
            image.push_back(row);
          }
        }
          break;
        default:
          throw std::runtime_error("right-to-left coordinate systems not implemented");
      }
      return image;
    }
    template<size_t N>
    std::vector<std::vector<Color>> image(const std::array<Color, N>& map, coordinates::system cs) const {
      T min_val{map_.val(0,0)}, max_val{map_.val(0,0)};
      for (const auto & val: map_.valItr()) {
        if (val < min_val) min_val = val;
        if (val > max_val) max_val = val;
      }
      auto range = max_val - min_val;
      std::vector<std::vector<Color>> image;
      const auto R{rows()};
      const auto C{cols()};
      image.reserve(R);
      auto which = [&](const T v){
        auto w = static_cast<size_t>((static_cast<T>(N)-1)*(v - min_val)/range);
        return std::min(N-1, w);
      };
      auto right = [&](const size_t i){
        std::vector<Color> row; row.reserve(C);
        for (ind_t j=0; j<C; ++j) row.push_back(map[which(map_.val(i, j))]);
        return row;
      };
      auto left = [&](const size_t i){
        std::vector<Color> row; row.reserve(C);
        for (ind_t j=C; j-->0;) row.push_back(map[which(map_.val(i, j))]);
        return row;
      };
      auto down = [&](auto & row){for (ind_t i=0; i<R; ++i) image.push_back(row(i));};
      auto up = [&](auto & row){for (ind_t i=R; i-->0; ) image.push_back(row(i));};
      switch(cs){
        case coordinates::system::y_up_x_right: up(right); break;
        case coordinates::system::y_down_x_right: down(right); break;
        case coordinates::system::y_up_x_left: up(left); break;
        case coordinates::system::y_down_x_left: down(left); break;
        default: throw std::runtime_error("coordinate systems not implemented");
      }
      return image;
    }


    Array2<T> values() const {return map_;}
    Bitmap<T> dilate(const ind_t n) const;
    Bitmap<T> erode(const ind_t n) const;
    Bitmap<T> dilate_edges(const ind_t n) const;
    Bitmap<T> erode_edges(const ind_t n) const;

    std::vector<polygon::Poly<int, Array2>> extract_polygons(coordinates::system, T level, T tol=T(0), int dig=1) const;
//    polygon::Poly<int, Array2> extract_polygon(Array2<bool>&, coordinates::Point<ind_t>, coordinates::system, T, T, int) const;
  };

  template<class T> using filter_queue_t = std::deque<std::pair<T, ind_t>>;
  /* filters from https://hal.science/hal-00692897/document
 * */
//  template<class T>
//  std::optional<T> filter1d(std::function<bool(T,T)> comparer,
//                            ind_t read, ind_t write, T value, ind_t left, ind_t right, ind_t N,
//                            filter_queue_t<T> & fifo){
  template<class T, class Comparer>
  std::optional<T> filter1d(Comparer comparer,
                            ind_t read, ind_t write, T value, ind_t left, ind_t right, ind_t N,
                            filter_queue_t<T> & fifo){
    // remove un-needed values (comparer is <= for dilation; >= for erosion)
    while (!fifo.empty() && comparer(fifo.back().first, value)) fifo.pop_back();
    // delete out-of-window value, if present
    if (!fifo.empty() && (left + fifo.front().second) < write) fifo.pop_front();
    // add the current value to the queue
    fifo.emplace_back(value, read);
    // return the front value if we've reached the end of the window or line
    return std::min(N-1, write + right) == read ? std::make_optional(fifo.front().first) : std::nullopt;
  }

  /*              comparer         pad
   * For dilation    <=     numeric_limits<T>::min
   *     erosion     >=     numeric_limits<T>::max
   * */
//  template<class T>
//  Bitmap<T> filter2d(std::function<bool(T,T)> comparer, T pad, const Bitmap<T>& image,
//                     ind_t left, ind_t right, ind_t top, ind_t bottom){
  template<class T, class Comparer>
  Bitmap<T> filter2d(Comparer comparer, T pad, const Bitmap<T>& image,
                     ind_t left, ind_t right, ind_t top, ind_t bottom){
    const auto M{image.rows()};
    const auto N{image.cols()};
    Bitmap<T> out(M, N);
    std::vector<filter_queue_t<T>> queues(image.cols());
    ind_t line_read{0}, line_write{0};
    while (line_write < M){
      bool wrote_to_line{false};
      filter_queue_t<T> fifo;
      ind_t column_read{0}, column_write{0};
      while (column_write < N){
        // horizontal filter on the lines_written line
        auto horizontal = std::make_optional(pad);
        if (line_read < M){
          auto value = column_read < N ? image.value(line_read, column_read++) : pad;
          horizontal = filter1d(comparer, std::min(column_read, N-1), column_write+1, value, left, right, N, fifo);
        }
        // vertical filter on the columns_written column
        if (horizontal.has_value()){
          auto vertical = filter1d(comparer, std::min(line_read, M-1), line_write, *horizontal, top, bottom, M, queues[column_write]);
          if (vertical.has_value()) {
            out.value(line_write, column_write) = *vertical;
            wrote_to_line = true;
          }
          column_write++;
        }
      }
      line_read++;
      if (wrote_to_line) line_write++;
    }
    return out;
  }

  template<class T> bool dilate_compare(const T & a, const T & b) {return a <= b;}
  template<class T> bool erode_compare(const T & a, const T & b) {return a >= b;}

  template<class T> Bitmap<T> Bitmap<T>::dilate(ind_t n) const {
    return filter2d(dilate_compare<T>, (std::numeric_limits<T>::min)(), *this, n, n, n, n);
  }
  template<class T> Bitmap<T> Bitmap<T>::erode(ind_t n) const {
    return filter2d(erode_compare<T>, (std::numeric_limits<T>::max)(), *this, n, n, n, n);
  }
  template<class T> Bitmap<T> Bitmap<T>::dilate_edges(ind_t n) const {
    return Bitmap<T>(dilate(n).values() - map_);
  }
  template<class T> Bitmap<T> Bitmap<T>::erode_edges(ind_t n) const {
    return Bitmap<T>(map_ - erode(n).values());
  }



//  template<class T> std::vector<polygon::Poly<int, Array2>>
//    Bitmap<T>::extract_polygons(const coordinates::system cs, const T level, const T tol, const int dig) const {
//    // go searching for part of a polygon by rastering the image:
//    Array2<bool> visited(map_.size(0), map_.size(1), false);
//    std::vector<polygon::Poly<int, Array2>> out;
//    for (const auto & sub: map_.subItr()) if (!visited[sub]) {
//      if (approx_float::scalar(map_[sub], level, tol, tol, dig)) {
//        out.push_back(extract_polygon(visited, sub, cs, level, tol, dig));
//      }
//    }
//    return out;
//  }
//  template<class T> polygon::Poly<int, Array2> Bitmap<T>::extract_polygon(
//    Array2<bool>& visited, coordinates::Point<ind_t> sub0, const coordinates::system cs, const T level, const T tol, const int dig
//    ) const {
//    using namespace coordinates;
//    // try to walk the border of this polygon
//    auto inbound = [s0=map_.size(0),s1=map_.size(1)](auto & x){
//      auto c = x.coord();
//      if (c[0] < 0) return false;
//      if (static_cast<ind_t>(c[0]) >= s0) return false;
//      if (c[1] < 0) return false;
//      if (static_cast<ind_t>(c[1]) >= s1) return false;
//      return true;
//    };
//    auto on_level = [&](auto & x){
//      return approx_float::scalar(map_[as<ind_t>(x).coord()], level, tol, tol, dig);
//    };
//
//    std::vector<std::array<int, 2>> vertices;
//    auto directions = ccw_directions<int>(cs);
//    size_t dir{0};
//    auto first_pos = as<int>(sub0);
//    auto last_pos = first_pos;
//    auto next_pos = last_pos.coord()[0] > 0 ? last_pos + directions[dir] : last_pos;
//    bool stepped{true};
//    size_t step{0}; // protect against the first iteration not going anywhere
////    std::cout << "Walk edge from " << first_pos << "\n";
//    do {
//      bool turn{false};
//      if (inbound(next_pos)) {
//        visited[as<ind_t>(next_pos).coord()] = true;
//        if (on_level(next_pos)) {
//          last_pos = next_pos;
//          next_pos += directions[dir];
//          stepped = true;
//          step++;
//        } else {
//          turn = true;
//        }
//      } else {
//        turn = true;
//      }
//      if (turn) {
//        dir = (dir + 1) % directions.size();
//        next_pos = last_pos + directions[dir];
//        if (stepped){
//          stepped = false;
//          vertices.push_back(as<int>(last_pos).coord());
//        }
//      }
//    } while (step < 1 || last_pos != first_pos);
//
//    auto poly_vertices = Array2<int>::from_std(vertices);
//    std::vector<ind_t> border(static_cast<size_t>(poly_vertices.size(0)));
//    std::iota(border.begin(), border.end(), 0u);
//    return polygon::Poly<int, Array2>(poly_vertices, polygon::Wires(border));
//  }

  std::tuple<int, std::vector<coordinates::Point<int>>>
  trace_path(Array2<bool> input, const coordinates::system cs, const coordinates::Point<int> & start);

  std::vector<coordinates::Point<int>> fix_path(const std::vector<coordinates::Point<int>> &);

  template<class T> std::vector<polygon::Poly<int, Array2>>
  Bitmap<T>::extract_polygons(const coordinates::system cs, const T level, const T tol, const int dig) const {
    // make a new bitmap which is 1 for 'level' and 0 otherwise
    Array2<bool> bitmap(map_.size(0), map_.size(1), false);
    using approx_float::scalar;
    for (const auto & sub: map_.subItr()) bitmap[sub] = scalar(map_[sub], level, tol, tol, dig);

    int minimum_area = 2;
    using namespace coordinates;

    std::vector<polygon::Poly<int, Array2>> out;
    for (const auto & sub: bitmap.subItr()) if (bitmap[sub]) {
        auto [area, points] = trace_path(bitmap, cs, as<int>(Point(sub)));
        if (std::abs(area) >= minimum_area){
          // convert the pixel path to a polygon path
//          auto path = points;
          auto path = fix_path(points);
          std::vector<std::array<ind_t,2>> vertices;
          vertices.reserve(path.size());
          std::transform(path.begin(), path.end(), std::back_inserter(vertices),
                         [](const auto & p){return as<ind_t>(p).coord();});
          std::vector<ind_t> border(path.size());
          std::iota(border.begin(), border.end(), 0u);
          out.emplace_back(Array2<int>::from_std(vertices), polygon::Wires(border));
        }
      }
    return out;
  }

} // end namespace polystar::bitmap

#endif