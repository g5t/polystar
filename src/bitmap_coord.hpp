//
// Created by gst on 08/02/23.
//

#ifndef POLYSTAR_NEW_BITMAP_COORD_HPP
#define POLYSTAR_NEW_BITMAP_COORD_HPP

#include <array>
#include <vector>
namespace polystar::bitmap {
  namespace coordinates {
    enum class system {
      y_up_x_right, y_down_x_right, y_up_x_left, y_down_x_left, none
    };

    template<class T> class Point;

    template<class R, class T>
    Point<R> as(Point<T> a) {
      auto c = a.coord();
      return {static_cast<R>(c[0]), static_cast<R>(c[1])};
    }

    template<class T>
    class Point {
    public:
      using coord_t = std::array<T, 2>;
    protected:
      coord_t coord_;
    public:
      Point(T a, T b) : coord_({a, b}) {}

      Point(coord_t c) : coord_(std::move(c)) {}

      Point<T> operator+(const Point<T> &o) const { return {coord_[0] + o.coord_[0], coord_[1] + o.coord_[1]}; }

      Point<T> operator-(const Point<T> &o) const { return {coord_[0] - o.coord_[0], coord_[1] - o.coord_[1]}; }

      Point<T> &operator+=(const Point<T> &o) {
        coord_[0] += o.coord_[0];
        coord_[1] += o.coord_[1];
        return *this;
      }

      Point<T> &operator-=(const Point<T> &o) {
        coord_[0] -= o.coord_[0];
        coord_[1] -= o.coord_[1];
        return *this;
      }

      Point<T> clamp(const T minimum, const T maximum) const {
        auto cl = [&](const auto x){return x < minimum ? minimum : x > maximum ? maximum : x;};
        return {cl(coord_[0]), cl(coord_[1])};
      }

      template<class R, class S=std::common_type_t<T,R>>
      S operator*(const Point<R> & o) const {
        auto c = as<S>(o).coord();
        return static_cast<S>(coord_[0]) * c[0] + static_cast<S>(coord_[1]) * c[1];
      }

      template<class R, class S=std::common_type_t<T,R>>
      Point<S> operator*(const R r) const {
        return {static_cast<S>(coord_[0]) * static_cast<S>(r),
                static_cast<S>(coord_[1]) * static_cast<S>(r)};
      }
      template<class R, class S=std::common_type_t<T,R>>
      Point<S> operator/(const R r) const {
        return {static_cast<S>(coord_[0]) / static_cast<S>(r),
                static_cast<S>(coord_[1]) / static_cast<S>(r)};
      }

      bool operator!=(const Point<T> &o) const { return (coord_[0] != o.coord_[0]) || (coord_[1] != o.coord_[1]); }
      bool operator==(const Point<T> &o) const { return (coord_[0] == o.coord_[0]) && (coord_[1] == o.coord_[1]); }

      coord_t coord() const { return coord_; }

      friend std::ostream &operator<<(std::ostream &os, const Point<T> &p) {
        auto c = p.coord();
        os << "(" << c[0] << "," << c[1] << ")";
        return os;
      }
    };

    template<class T> T cross(const Point<T> & a, const Point<T> & b) {
      auto ca = a.coord(); // {y, x}
      auto cb = b.coord();
      return ca[1] * cb[0] - ca[0] * cb[1];
    }
    template<class T> Point<T> abs(const Point<T>& a) {
      if constexpr (std::is_unsigned_v<T>) {return a;}
      auto c = a.coord();
      return {std::abs(c[0]), std::abs(c[1])};
    }

    template<class T>
    std::vector <Point<T>> ccw_directions(system value) {
      std::vector <Point<T>> yuxr{{0,  -1},
                                  {-1, -1},
                                  {-1, 0},
                                  {-1, 1},
                                  {0,  1},
                                  {1,  1},
                                  {1,  0},
                                  {1,  -1}},
        ydxr{{0,  -1},
             {1,  -1},
             {1,  0},
             {1,  1},
             {0,  1},
             {-1, 1},
             {-1, 0},
             {-1, -1}},
        yuxl{{0,  1},
             {-1, 1},
             {-1, 0},
             {-1, -1},
             {0,  -1},
             {1,  -1},
             {1,  0},
             {1,  1}},
        ydxl{{0,  1},
             {1,  1},
             {1,  0},
             {1,  -1},
             {0,  -1},
             {-1, -1},
             {-1, 0},
             {-1, 1}};
      switch (value) {
        case system::y_up_x_right:
          return yuxr;
        case system::y_down_x_right:
          return ydxr;
        case system::y_up_x_left:
          return yuxl;
        case system::y_down_x_left:
          return ydxl;
        default:
          throw std::runtime_error("Unknown enumerated coordinate system value");
      }
    }
  }
}
#endif //POLYSTAR_NEW_BITMAP_COORD_HPP
