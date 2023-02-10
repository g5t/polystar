#ifndef POLYSTAR_POLYGON_WIRE_HPP
#define POLYSTAR_POLYGON_WIRE_HPP
#include <vector>
#include <array>
#include <optional>
#include <numeric>
#include <limits>

#include "comparisons.hpp"
#include "array_.hpp"
#include "geometry.hpp"
#include "approx_float.hpp"

namespace polystar::polygon {
  using ind_t = polystar::ind_t;

  class Wire : public std::vector<ind_t> {
  public:
    using edge_t = std::pair<ind_t, ind_t>;
    using base_t = std::vector<ind_t>;
  private:
    inline ind_t vi(size_t i) const { return operator[](i); }

    inline ind_t vj(size_t i) const { return operator[]((i + 1) % size()); }

  public:
    inline edge_t edge(size_t i) const { return std::make_pair(vi(i), vj(i)); }

    // apparently the empty constructor is needed for explicitness?
    explicit Wire() {}

    // type conversion (require that T is integer?)
    template<class T>
    explicit Wire(const std::vector <T> &indexes) {
      reserve(indexes.size());
      std::transform(indexes.begin(), indexes.end(),
                     std::back_insert_iterator < std::vector < ind_t >> (*this),
                     [](const T &x) { return static_cast<ind_t>(x); });
    }

    template<class T, template<class> class A, class U = typename A<T>::shape_t>
    explicit Wire(const A<T> &vertices, const T tol = T(0), const int dig = 1) {
      // Convex Hull construction
      auto unique = vertices.is_unique(tol, dig);
      auto v = vertices.extract(unique);
      auto chi = find_convex_hull_edges(v, tol, dig);
      if (!chi.empty()) {
        // the wire should index vertices but indexes the unique vertices instead at the moment:
        // the Convex Hull algorithm only returned strictly-necessary-point indexes, so no extra work required
        std::vector <ind_t> index_map(v.size(0), vertices.size(0));
        for (ind_t i = 0; i < v.size(0); ++i) index_map[i] = vertices.first(cmp::eq, v.view(i));
        reserve(chi.size());
        for (auto &index: chi) push_back(index_map[index]);
      }
    }

    std::back_insert_iterator <base_t> appender() { return std::back_insert_iterator<base_t>(*this); }

    template<class T, template<class> class A>
    T area(const A<T> &x) const {
      T a{0};
      for (size_t i = 0; i < size(); ++i) {
        a += cross2d(x.view(vi(i)), x.view(vj(i))).val(0, 0);
      }
      return a / T(2);
    }

    template<class T, template<class> class A>
    A<T> centroid(const A<T> &x) const {
      auto c = 0 * x.view(0);
      auto twice_area = T(0);
      for (size_t i = 0; i < size(); ++i) {
        const auto &a{x.view(vi(i))};
        const auto &b{x.view(vj(i))};
        auto ta = cross2d(a, b).val(0);
        c += (a + b) * ta;
        twice_area += ta;
      }
      return c / (T(3) * twice_area);
    }

    Wire mirror() const {
      auto out = Wire();
      out.reserve(size());
      std::copy(crend(), crbegin(), std::back_insert_iterator<base_t>(out));
      return out;
    }

    template<class T, template<class> class A>
    bool contains(const A<T> &point, const A<T> &x) const {
      auto segment = cat(0, point, (std::numeric_limits<T>::max)() + T(0) * point);
      edge_t se{0, 1};
      size_t crossings{0};
      for (size_t i = 0; i < size(); ++i) if (intersect2d(segment, se, x, edge(i))) ++crossings;
      return (crossings % 2u) == 1u;
    }

    template<class T, template<class> class A>
    bool intersects(const A<T> &other, const edge_t &oe, const A<T> &x) const {
      for (size_t i = 0; i < size(); ++i) if (intersect2d(other, oe, x, edge(i))) return true;
      return false;
    }

    template<class T, template<class> class A>
    A<T> intersections(const A<T> &other, const edge_t &oe, const A<T> &x) const {
      size_t count{0};
      auto out = A<T>();
      for (size_t i = 0; i < size(); ++i) {
        auto n_x = intersection2d(other, oe, x, edge(i));
        if (n_x.first > 0) out = (count++) ? cat(0, out, n_x.second) : n_x.second;
      }
      if (count > 1) {
        // multiple edges intersected with the other edge:
        // if an intersection point *is* a vertex then it would be the same for two neighbouring edges
        // check for repeated vertices (or for vertices in x)?
        out = out.extract(out.unique());
      }
      return out;
    }

    template<class T>
    Wire replace(const std::vector <T> &replacement) const {
      Wire p;
      p.reserve(this->size());
      std::transform(this->cbegin(), this->cend(), std::back_insert_iterator < std::vector < ind_t >> (p),
                     [&replacement](const T &index) { return replacement[index]; });
      return p;
    }

    template<class T, template<class> class A>
    [[nodiscard]] T circumscribed_radius(const A<T> &x) const {
      auto c = centroid(x);
      T r{0};
      for (const auto &i: *this) if (auto t = norm(c - x.view(i)); t > r) r = t;
      return r;
    }

    Wire offset(const ind_t value) const {
      Wire o;
      o.reserve(size());
      std::transform(cbegin(), cend(), o.appender(), [value](const ind_t &i) { return i + value; });
      return o;
    }

    Wire combine(const Wire &that, const ind_t value) const {
      Wire o;
      o.reserve(size() + that.size());
      std::copy(cbegin(), cend(), o.appender());
      auto t = that.offset(value);
      std::copy(t.cbegin(), t.cend(), o.appender());
      return o;
    }

    template<class T, class R, template<class> class A,
      template<class> class B>
    std::tuple <A<T>, Wire> cut(const A<T> &x, const B<R> &a, const B<R> &b) const {
//      auto is_b = is_beyond(x, a, b);
//      A<T> xo(0u, 2u);
//      if (!std::any_of(is_b.begin(), is_b.end(), true)) return std::make_tuple(x, Wire(*this));
//      if (!std::all_of(is_b.begin(), is_b.end(), true)) return std::make_tuple(A<T>(), Wire());
//      // otherwise, go-around the edges and find intersections
//      auto ov = cat(0, a, b);
//      edge_t oe{0,1};
//      for (size_t i=0; i<size(); ++i){
//        auto [n, p] = intersection2d(ov, oe, x, edge(i));
//        if (n == 1) {
//
//        }
//      }
      return std::make_tuple(x, Wire(*this));
    }

    template<class T, class R, template<class> class A,
      template<class> class B>
    bool any_beyond(const A<T> &x, const B<R> &a, const B<R> &b) const {
      assert(a.size(0) == b.size(0) && a.size(0) == 1);
      auto v = b - a;
      for (const auto &i: *this) if (cross2d(v, x.view(i)).val(0, 0) < 0) return true;
      return false;
    }

    template<class T, class R, template<class> class A,
      template<class> class B>
    std::vector<bool> is_beyond(const A<T> &x, const B<R> &a, const B<R> &b) const {
      assert(a.size(0) == b.size(0) && a.size(0) == 1);
      auto v = b - a;
      std::vector<bool> out;
      out.reserve(size());
      std::transform(cbegin(), cend(), std::back_inserter(out),
                     [&x, &v](const ind_t i) { return cross2d(v, x.view(i)).val(0, 0) < 0; });
      return out;
    }

    std::string python_string() const {
      std::stringstream s;
      s << "[";
      size_t count{0};
      for (const auto &i: *this) s << i << (++count < size() ? "," : "");
      s << "]";
      return s.str();
    }

    friend std::ostream &operator<<(std::ostream &os, Wire w) {
      os << w.python_string();
      return os;
    }

  };

} // polystar::polygon
#endif