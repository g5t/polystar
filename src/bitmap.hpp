#ifndef POLYSTAR_BITMAP_HPP
#define POLYSTAR_BITMAP_HPP

#include <vector>
#include <deque>
#include <functional>
#include <optional>
#include <map>
#include <fstream>

#include "array_.hpp"
#include "polygon.hpp"

namespace polystar::bitmap {
  class Color {
  public:
    using t = uint8_t;
  private:
    t r_=0;
    t g_=0;
    t b_=0;
  public:
    Color(t r, t g, t b): r_(r), g_(g), b_(b) {}
    [[nodiscard]] t r() const {return r_;}
    [[nodiscard]] t g() const {return g_;}
    [[nodiscard]] t b() const {return b_;}
    [[nodiscard]] const char * R() const {return (const char *)&r_;}
    [[nodiscard]] const char * B() const {return (const char *)&b_;}
    [[nodiscard]] const char * G() const {return (const char *)&g_;}
    [[nodiscard]] uint64_t gray() const {
      // this could be color-model dependent
      return static_cast<uint64_t>(r_) + static_cast<uint64_t>(g_) + static_cast<uint64_t>(b_);
    }
    void write(std::ofstream& of) const {
      of.write((char*)&b_, sizeof(uint8_t));
      of.write((char*)&g_, sizeof(uint8_t));
      of.write((char*)&r_, sizeof(uint8_t));
    }
    bool operator<(Color o) const {return gray() < o.gray();}
    bool operator>(Color o) const {return gray() > o.gray();}
    bool operator==(Color o) const {return r_==o.r() && g_==o.g() && b_==o.b();}
    bool operator!=(Color o) const {return r_!=o.r() || g_!=o.g() || b_!=o.b();}
  };

  void write(const std::vector<std::vector<Color>> & image, const std::string & filename);

  template<class T>
  void write(const std::map<T, Color>& map, const std::vector<std::vector<T>>& data, const std::string & filename){
    std::vector<std::vector<Color>> image;
    image.reserve(data.size());
    auto width = data[0].size();
    for (const auto & row: data) if (row.size() > width) width = row.size();
    for (const auto & row: data){
      std::vector<Color> values(width, Color(0,0,0));
      for (size_t i=0; i<row.size(); ++i) values[i] = map.at(row[i]);
      image.push_back(values);
    }
    write(image, filename);
  }
  namespace colormaps {
      // similar to original from https://doi.org/10.1371/journal.pone.0199239 but modified for uint8 values
      static const
      std::array<Color,256> cividis {{{  0,  35,  78}, {  0,  35,  80}, {  0,  36,  81}, {  0,  37,  83},
                                     {  0,  38,  85}, {  0,  38,  86}, {  0,  39,  88}, {  0,  40,  90},
                                     {  0,  40,  92}, {  0,  41,  93}, {  0,  42,  95}, {  0,  42,  97},
                                     {  0,  43,  99}, {  0,  44, 101}, {  0,  44, 102}, {  0,  45, 104},
                                     {  0,  46, 106}, {  0,  46, 108}, {  0,  47, 110}, {  0,  48, 111},
                                     {  0,  48, 113}, {  0,  49, 113}, {  0,  50, 113}, {  1,  50, 113},
                                     {  5,  51, 113}, {  8,  52, 113}, { 12,  52, 113}, { 15,  53, 113},
                                     { 18,  54, 112}, { 20,  54, 112}, { 22,  55, 112}, { 24,  56, 112},
                                     { 26,  56, 112}, { 28,  57, 111}, { 30,  58, 111}, { 32,  59, 111},
                                     { 33,  59, 111}, { 35,  60, 111}, { 37,  61, 110}, { 38,  61, 110},
                                     { 39,  62, 110}, { 41,  63, 110}, { 42,  63, 110}, { 44,  64, 110},
                                     { 45,  65, 110}, { 46,  66, 109}, { 47,  66, 109}, { 49,  67, 109},
                                     { 50,  68, 109}, { 51,  68, 109}, { 52,  69, 109}, { 53,  70, 109},
                                     { 55,  70, 109}, { 56,  71, 109}, { 57,  72, 109}, { 58,  73, 108},
                                     { 59,  73, 108}, { 60,  74, 108}, { 61,  75, 108}, { 62,  75, 108},
                                     { 63,  76, 108}, { 64,  77, 108}, { 65,  77, 108}, { 66,  78, 108},
                                     { 68,  79, 108}, { 69,  79, 108}, { 70,  80, 108}, { 71,  81, 108},
                                     { 72,  82, 108}, { 73,  82, 108}, { 73,  83, 108}, { 74,  84, 108},
                                     { 75,  84, 108}, { 76,  85, 109}, { 77,  86, 109}, { 78,  86, 109},
                                     { 79,  87, 109}, { 80,  88, 109}, { 81,  89, 109}, { 82,  89, 109},
                                     { 83,  90, 109}, { 84,  91, 109}, { 85,  91, 109}, { 86,  92, 109},
                                     { 87,  93, 110}, { 88,  93, 110}, { 89,  94, 110}, { 89,  95, 110},
                                     { 90,  96, 110}, { 91,  96, 110}, { 92,  97, 110}, { 93,  98, 111},
                                     { 94,  98, 111}, { 95,  99, 111}, { 96, 100, 111}, { 97, 100, 111},
                                     { 97, 101, 112}, { 98, 102, 112}, { 99, 103, 112}, {100, 103, 112},
                                     {101, 104, 112}, {102, 105, 113}, {103, 105, 113}, {104, 106, 113},
                                     {105, 107, 113}, {105, 108, 114}, {106, 108, 114}, {107, 109, 114},
                                     {108, 110, 114}, {109, 110, 115}, {110, 111, 115}, {111, 112, 115},
                                     {111, 113, 115}, {112, 113, 116}, {113, 114, 116}, {114, 115, 116},
                                     {115, 115, 117}, {116, 116, 117}, {116, 117, 117}, {117, 118, 118},
                                     {118, 118, 118}, {119, 119, 119}, {120, 120, 119}, {121, 121, 119},
                                     {122, 121, 120}, {122, 122, 120}, {123, 123, 120}, {124, 124, 120},
                                     {125, 124, 121}, {126, 125, 121}, {127, 126, 121}, {128, 126, 121},
                                     {129, 127, 121}, {130, 128, 121}, {131, 129, 121}, {132, 129, 121},
                                     {133, 130, 121}, {134, 131, 121}, {134, 132, 121}, {135, 132, 121},
                                     {136, 133, 121}, {137, 134, 121}, {138, 135, 121}, {139, 135, 121},
                                     {140, 136, 121}, {141, 137, 121}, {142, 138, 121}, {143, 138, 121},
                                     {144, 139, 120}, {145, 140, 120}, {146, 141, 120}, {147, 141, 120},
                                     {148, 142, 120}, {149, 143, 120}, {150, 144, 120}, {151, 145, 120},
                                     {152, 145, 119}, {153, 146, 119}, {154, 147, 119}, {155, 148, 119},
                                     {156, 148, 119}, {157, 149, 119}, {158, 150, 118}, {159, 151, 118},
                                     {160, 152, 118}, {161, 152, 118}, {162, 153, 117}, {163, 154, 117},
                                     {164, 155, 117}, {165, 156, 117}, {166, 156, 116}, {167, 157, 116},
                                     {168, 158, 116}, {169, 159, 116}, {170, 159, 115}, {171, 160, 115},
                                     {172, 161, 115}, {173, 162, 114}, {174, 163, 114}, {175, 164, 114},
                                     {176, 164, 113}, {177, 165, 113}, {178, 166, 113}, {179, 167, 112},
                                     {180, 168, 112}, {181, 168, 112}, {182, 169, 111}, {183, 170, 111},
                                     {184, 171, 110}, {185, 172, 110}, {186, 173, 109}, {188, 173, 109},
                                     {189, 174, 109}, {190, 175, 108}, {191, 176, 108}, {192, 177, 107},
                                     {193, 178, 107}, {194, 178, 106}, {195, 179, 106}, {196, 180, 105},
                                     {197, 181, 105}, {198, 182, 104}, {199, 183, 103}, {200, 184, 103},
                                     {201, 184, 102}, {202, 185, 102}, {203, 186, 101}, {204, 187, 101},
                                     {205, 188, 100}, {207, 189,  99}, {208, 190,  99}, {209, 191,  98},
                                     {210, 191,  97}, {211, 192,  97}, {212, 193,  96}, {213, 194,  95},
                                     {214, 195,  94}, {215, 196,  94}, {216, 197,  93}, {217, 198,  92},
                                     {218, 199,  91}, {220, 199,  90}, {221, 200,  90}, {222, 201,  89},
                                     {223, 202,  88}, {224, 203,  87}, {225, 204,  86}, {226, 205,  85},
                                     {227, 206,  84}, {228, 207,  83}, {230, 208,  82}, {231, 209,  81},
                                     {232, 209,  80}, {233, 210,  79}, {234, 211,  78}, {235, 212,  77},
                                     {236, 213,  76}, {237, 214,  74}, {239, 215,  73}, {240, 216,  72},
                                     {241, 217,  71}, {242, 218,  69}, {243, 219,  68}, {244, 220,  66},
                                     {246, 221,  65}, {247, 222,  63}, {248, 223,  62}, {249, 224,  60},
                                     {250, 225,  58}, {252, 226,  57}, {253, 226,  55}, {254, 227,  53},
                                     {255, 229,  52}, {255, 230,  53}, {255, 231,  54}, {255, 233,  56}}};
    }

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