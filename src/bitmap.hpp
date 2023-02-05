#ifndef POLYSTAR_BITMAP_HPP
#define POLYSTAR_BITMAP_HPP

#include <vector>
#include <deque>
#include <functional>
#include <optional>

#include "array_.hpp"
#include "polygon.hpp"

namespace polystar::bitmap {
  namespace coordinates {
    enum class system {y_up_x_right, y_down_x_right, y_up_x_left, y_down_x_left, none};

    template<class T>
    class Point {
    public:
      using coord_t = std::array<T, 2>;
    protected:
      coord_t coord_;
    public:
      Point(T a, T b): coord_({a, b}) {}
      Point(coord_t c): coord_(std::move(c)) {}
      Point<T> operator+(const Point<T>& o) const { return {coord_[0] + o.coord_[0], coord_[1] + o.coord_[1]}; }
      Point<T> operator-(const Point<T>& o) const { return {coord_[0] - o.coord_[0], coord_[1] - o.coord_[1]}; }
      Point<T>& operator+=(const Point<T>& o) {
        coord_[0] += o.coord_[0];
        coord_[1] += o.coord_[1];
        return *this;
      }
      Point<T>& operator-=(const Point<T>& o) {
        coord_[0] -= o.coord_[0];
        coord_[1] -= o.coord_[1];
        return *this;
      }
      bool operator!=(const Point<T>& o) const {return (coord_[0] != o.coord_[0]) || (coord_[1] != o.coord_[1]);}
      coord_t coord() const {return coord_;}
    };

    template<class R, class T> Point<R> as(Point<T> a){
      auto c = a.coord();
      return {static_cast<R>(c[0]), static_cast<R>(c[1])};
    }

    template<class T> std::vector<Point<T>> ccw_directions(system value){
      std::vector<Point<T>> yuxr{{0, -1}, {-1, -1}, {-1, 0}, {-1, 1}, {0, 1}, {1, 1}, {1, 0}, {1, -1}},
      ydxr{{0, -1}, {1, -1}, {1, 0}, {1, 1}, {0, 1}, {-1, 1}, {-1, 0}, {-1, -1}},
      yuxl{{0, 1}, {-1, 1}, {-1, 0}, {-1, -1}, {0, -1}, {1, -1}, {1, 0}, {1, 1}},
      ydxl{{0, 1}, {1, 1}, {1, 0}, {1, -1}, {0, -1}, {-1, -1}, {-1, 0}, {-1, 1}};
      switch (value){
        case system::y_up_x_right: return yuxr;
        case system::y_down_x_right: return ydxr;
        case system::y_up_x_left: return yuxl;
        case system::y_down_x_left: return ydxl;
        default: throw std::runtime_error("Unknown enumerated coordinate system value");
      }
    }
  }

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

    Array2<T> values() const {return map_;}
    Bitmap<T> dilate(const ind_t n) const;
    Bitmap<T> erode(const ind_t n) const;
    Bitmap<T> dilate_edges(const ind_t n) const;
    Bitmap<T> erode_edges(const ind_t n) const;

    std::vector<polygon::Poly<ind_t, Array2>> extract_polygons(coordinates::system, T level, T tol=T(0), int dig=1) const;
    polygon::Poly<ind_t, Array2> extract_polygon(Array2<bool>&, coordinates::Point<ind_t>, coordinates::system, T, T, int) const;
  };

  template<class T> using filter_queue_t = std::deque<std::pair<T, ind_t>>;
  /* filters from https://hal.science/hal-00692897/document
 * */
  template<class T>
  std::optional<T> filter1d(std::function<bool(T,T)> comparer,
                            ind_t read, ind_t write, T value, ind_t left, ind_t right, ind_t N,
                            filter_queue_t<T> & fifo){
    // remove un-needed values (comparer is <= for dilation; >= for erosion)
    while (!fifo.empty() && comparer(fifo.back().first, value)) fifo.pop_back();
    // delete out-of-window value, if present
    if (!fifo.empty() && write > left + fifo.front().second) fifo.pop_front();
    // add the current value to the queue
    fifo.emplace_back(value, read);
    // return the front value if we've reached the end of the window or line
    return std::min(N, write + right) == read ? std::make_optional(fifo.front().first) : std::nullopt;
  }

  /*              comparer         pad
   * For dilation    <=     numeric_limits<T>::min
   *     erosion     >=     numeric_limits<T>::max
   * */
  template<class T>
  Bitmap<T> filter2d(std::function<bool(T,T)> comparer, T pad, const Bitmap<T>& image,
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
          horizontal = filter1d(comparer, std::min(column_read, N), column_write, value, left, right, N, fifo);
        }
        // vertical filter on the columns_written column
        if (horizontal.has_value()){
          auto vertical = filter1d(comparer, std::min(line_read, M), line_write, horizontal.value(), top, bottom, M, queues[column_write]);
          if (vertical.has_value()) {
            out.value(line_write, column_write) = vertical.value();
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
    auto e = n >> 1;
    return filter2d(dilate_compare<T>, (std::numeric_limits<T>::min)(), *this, e, e, e, e);
  }
  template<class T> Bitmap<T> Bitmap<T>::erode(ind_t n) const {
    auto e = n >> 1;
    return filter2d(erode_compare<T>, (std::numeric_limits<T>::max)(), *this, e, e, e, e);
  }
  template<class T> Bitmap<T> Bitmap<T>::dilate_edges(ind_t n) const {
    return Bitmap<T>(dilate(n).values() - map_);
  }
  template<class T> Bitmap<T> Bitmap<T>::erode_edges(ind_t n) const {
    return Bitmap<T>(map_ - erode(n).values());
  }






  template<class T> std::vector<polygon::Poly<ind_t, Array2>> Bitmap<T>::extract_polygons(const coordinates::system cs, const T level, const T tol, const int dig) const {
    // go searching for part of a polygon by rastering the image:
    Array2<bool> visited(map_.size(0), map_.size(1), false);
    std::vector<polygon::Poly<ind_t, Array2>> out;
    for (const auto & sub: map_.subItr()) if (!visited[sub]) {
      if (approx_float::scalar(map_[sub], level, tol, tol, dig)) out.push_back(extract_polygon(visited, sub, cs, level, tol, dig));
    }
    return out;
  }
  template<class T> polygon::Poly<ind_t, Array2> Bitmap<T>::extract_polygon(
    Array2<bool>& visited, coordinates::Point<ind_t> sub0, const coordinates::system cs, const T level, const T tol, const int dig
    ) const {
    using namespace coordinates;
    // try to walk the border of this polygon
    auto inbound = [s0=map_.size(0),s1=map_.size(1)](auto & x){
      auto c = x.coord();
      if (c[0] < 0) return false;
      if (static_cast<ind_t>(c[0]) >= s0) return false;
      if (c[1] < 0) return false;
      if (static_cast<ind_t>(c[1]) >= s1) return false;
      return true;
    };
    auto on_level = [&](auto & x){
      return approx_float::scalar(map_[as<ind_t>(x).coord()], level, tol, tol, dig);
    };

    std::vector<std::array<ind_t, 2>> vertices;
    auto directions = ccw_directions<int>(cs);
    size_t dir{0};
    auto last_dir = directions[0];
    auto first_pos = as<int>(sub0);
    auto last_pos = first_pos;
    auto next_pos = last_pos.coord()[0] > 0 ? last_pos + directions[dir] : last_pos;
    bool stepped{true};
    size_t step{0}; // protect against the first iteration not going anywhere
    do {
      bool turn{false};
      if (inbound(next_pos)) {
        visited[as<ind_t>(next_pos).coord()] = true;
        if (on_level(next_pos)) {
          last_pos = next_pos;
          next_pos += directions[dir];
          stepped = true;
          step++;
        } else {
          turn = true;
        }
      } else {
        turn = true;
      }
      if (turn) {
        dir = (dir + 1) % directions.size();
        next_pos += directions[dir];
        if (stepped){
          stepped = false;
          vertices.push_back(as<ind_t>(last_dir).coord());
        }
      }
    } while (step < 1 && last_pos != first_pos);

    auto poly_vertices = Array2<ind_t>::from_std(vertices);
    std::vector<ind_t> border(static_cast<size_t>(poly_vertices.size(0)));
    std::iota(border.begin(), border.end(), 0u);
    return polygon::Poly<ind_t, Array2>(poly_vertices, polygon::Wires(border));
  }

} // end namespace polystar::bitmap

#endif