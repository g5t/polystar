#ifndef POLYSTAR_POLYGON_WIRE_HPP
#define POLYSTAR_POLYGON_WIRE_HPP
#include <vector>
#include <array>
#include <optional>
#include <numeric>
#include <limits>
#include <list>

#include "comparisons.hpp"
#include "array_.hpp"
#include "geometry.hpp"
#include "approx_float.hpp"
#include "polygon_network.hpp"
#include "svg.hpp"

namespace polystar::polygon {
  using ind_t = polystar::ind_t;

  class Wire : public std::vector<ind_t> {
  public:
    using edge_t = std::pair<ind_t, ind_t>;
    using base_t = std::vector<ind_t>;
  private:
    /// \brief Return the ith entry, respecting cyclic boundary conditions
    [[nodiscard]] inline ind_t vi(size_t i) const { return operator[](i % size()); }
    /// \brief Return (i+1)th entry, respecting cyclic boundary conditions
    [[nodiscard]] inline ind_t vj(size_t i) const { return operator[]((i + 1) % size()); }

  public:
    /// \brief Return the ith edge pair, which is the ith and (i+1)th entries in that order
    [[nodiscard]] inline edge_t edge(size_t i) const { return std::make_pair(vi(i), vj(i)); }

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

    template<size_t N>
    explicit Wire(const std::array<ind_t, N> w) {
      reserve(N);
      for (const auto & x: w) push_back(x);
    }

    std::back_insert_iterator <base_t> appender() { return std::back_insert_iterator<base_t>(*this); }

    /// \brief Check whether the wire would give an out-of-bounds indexing error if accessing the provided Array2
    template<class T, template<class> class A>
    bool indexing_error(const A<T>& x) const {
      if (empty()) return false;
      auto max_ptr = std::max_element(begin(), end());
      return *max_ptr >= x.size(0u);
    }

    template<class T, template<class> class A>
    bool is_convex(const A<T>& x) const {
//      auto a{x.view(vj(0)) - x.view(vi(0))};
//      auto b{x.view(vj(1)) - x.view(vi(1))};
//      T first = cross2d(a, b).val(0,0);
//      a = b;
//      for (size_t i=2; i<size(); ++i){
//        auto second = cross2d(a, b).val(0, 0);
//        b = x.view(vj(i)) - a.view(vi(i));
//        if (first *  second < 0) return false;
//        std::swap(a, b);
//      }
//      return true;
      auto wa = wire_winding_number_over_pi(x);
      debug_update_if(abs(wa) != 2, "Winding number ", wa, "pi is not +/-2pi, so this wire is not convex");
      return abs(wa) == 2;
    }

    bool operator==(const Wire & that) const {
      return size() == that.size()
          && std::is_permutation(begin(), end(), that.begin())
          && is_positive_permutation(base_t(*this), base_t(that));
    }
    bool operator!=(const Wire & that) const {
      return size() != that.size()
          || !std::is_permutation(begin(), end(), that.begin())
          || !is_positive_permutation(base_t(*this), base_t(that));
    }

    bool shares_edge(const Wire& other, bool backwards_allowed = false) const {
      // in two dimensions any triangulated (or other-shape-ed) polygon will have internal edges
      // with the property that two shapes share opposite-vertex-ordered edges
      for (size_t i=0; i<size(); ++i) for (size_t j=0; j<other.size(); ++j) {
        auto ei = edge(i);
        auto ej = other.edge(j);
        if (ei.first == ej.second && ei.second == ej.first) return true;
        if (backwards_allowed && ei.first == ej.first && ei.second == ej.second) return true;
      }
      return false;
    }

    template<class T, template<class> class A>
      bool uses(const A<T> & point, const A<T> & x, const T tol = T(0), const int dig = 1) const {
      for (const auto & i: *this) if (x.view(i).count(polystar::cmp::eq, point, tol, tol, dig) == 1) return true;
      return false;
    }

    template<class T, template<class> class A>
    Network<Wire,T,A> triangulate(const A<T>& x) const {
      using group_t = std::array<ind_t,3>;
      std::vector<group_t> groups;
      groups.reserve(size());
      for (auto ptr = begin(); ptr != end(); ++ptr){
        auto prev = ptr == begin() ? end() - 1 : ptr - 1;
        auto next = ptr == end() - 1 ? begin() : ptr + 1;
        groups.push_back(std::array<ind_t,3>({*prev, *ptr, *next}));
      }
      std::vector<group_t> convex;
      std::vector<group_t> concave;
      for (const auto & g: groups){
        if (cross2d(x.view(g[1]) - x.view(g[0]), x.view(g[2]) - x.view(g[1])).val(0,0) > 0){
          convex.push_back(g);
        } else {
          concave.push_back(g);
        }
      }
      std::vector<group_t> ears;

      auto update_concave_convex = [&](const auto & g, const auto & ng){
        // if g[1] was concave, it might not be any longer
        if (auto cp = std::find(concave.begin(), concave.end(), g); cp != concave.end()){
          if (cross2d(x.view(ng[1]) - x.view(ng[0]), x.view(ng[2]) - x.view(ng[1])).val(0,0) > 0) {
            convex.push_back(ng);
            concave.erase(cp);
          } else {
            // update the group membership
            *cp = ng;
          }
        } else {
          // the point remains convex but needs to have its indexing updated
          auto xp = std::find(convex.begin(), convex.end(), g);
          if (xp == convex.end()) throw std::runtime_error("Point was not convex but should have been");
          *xp = ng;
        }
      };
      auto update_ears = [&](const auto & ip, const auto & g, const auto & ng){
        // this is called *after* update_convex_concave -- so we *must* check for ng in convex!
        // if g[i] is convex, it might have changed ear-ness:
        bool not_tri{ng[0] == ng[2]};
        if (!not_tri && std::find(convex.begin(), convex.end(), ng) != convex.end()) {
          auto tri = Wire(ng);
//          std::cout << "Inner (concave) points" << ip;
          auto cinc = tri.contains(ip, x);
          // Complex polygons can visit the same point multiple times, and it may be concave once and convex another
          for (ind_t i=0; i<ip.size(0); ++i) if (cinc[i] && tri.uses(ip.view(i), x)) cinc[i] = false;
//          for (const auto & q: cinc) std::cout << (q ? "1" : "0");
//          std::cout << "\n";
          auto is_ear = std::count(cinc.begin(), cinc.end(), true) == 0;
          auto ptr = std::find(ears.begin(), ears.end(), g);
//          std::cout << "Old group (" << g[0] << " " << g[1]  << " " << g[2] << ") new group (";
//          std::cout << ng[0] << " " << ng[1] << " " << ng[2] << ")";
          if (!is_ear){
            if (ptr != ears.end()) {
              ears.erase(ptr);
//              std::cout << " was an ear but no longer is\n";
//            } else {
//              std::cout << " remains a non-ear\n";
            }
          } else {
            if (ptr == ears.end()) {
              ears.push_back(ng);
//              std::cout << " has become an ear\n";
            } else {
              *ptr = ng;
//              std::cout << " remains an ear, and has its indexing updated\n";
            }
          }
        }
        // remove any non-triangular groups from ears
        if (not_tri) if (auto ptr = std::find(ears.begin(), ears.end(), g); ptr != ears.end()) ears.erase(ptr);
      };
      std::vector<ind_t> middle;
      middle.reserve(concave.size());
      for (const auto & g: concave) middle.push_back(g[1]);
      auto inner_points = x.extract(middle);
      for (const auto & g: groups) update_ears(inner_points, g, g);

      auto update_groups = [&](const auto & g, const auto & ng){
        auto ptr = std::find(groups.begin(), groups.end(), g);
        if (ptr == groups.end()) throw std::runtime_error("Group not found");
//        std::cout << "Update group (" << g[0] << " " << g[1]  << " " << g[2] << ") to (";
//        std::cout << ng[0] << " " << ng[1] << " " << ng[2] << ")\n";
        *ptr = ng;
      };

      auto ear_measure = [&](const auto & ear){
        const auto & v0{x.view(ear[0])};
        const auto & v1{x.view(ear[1])};
        const auto & v2{x.view(ear[2])};
        auto a = norm(v0 - v1).sum();
        auto b = norm(v1 - v2).sum();
        auto c = norm(v2 - v0).sum();
        auto area2 = (a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c);
        auto area = std::sqrt(area2);
        // true circumscribed diameter, but we don't need this
        // auto d = std::min(std::max(std::max(a, b), c), a * b * c / 2 / area);
        return area * area2 / (a * b * c);
      };


      // now we only need to work with the known ears (and check the neighbors after removing each)
      std::vector<Wire> triangles;
      while (!ears.empty()){
//        std::cout << "\nRemaining groups: ";
//        for (const auto & e: groups) std::cout << "(" << e[0] << " " << e[1] << " " << e[2] << ") ";
//        std::cout << "\nRemaining convex: ";
//        for (const auto & e: convex) std::cout << "(" << e[0] << " " << e[1] << " " << e[2] << ") ";
//        std::cout << "\nRemaining concave: ";
//        for (const auto & e: concave) std::cout << "(" << e[0] << " " << e[1] << " " << e[2] << ") ";
//        std::cout << "\nRemaining ear indexes: ";
//        for (const auto & e: ears) std::cout << "(" << e[0] << " " << e[1] << " " << e[2] << ") ";
//        std::cout << "\n";
        // pick the most compact ear to be the first triangle, erase it from the ears at this point:
        std::partial_sort(ears.begin(), ears.begin()+1, ears.end(), [&](const auto & e1, const auto & e2){
          return ear_measure(e1) > ear_measure(e2);
        });
        group_t i = ears.front();
        ears.erase(ears.begin());
        // find its group
        auto ptr = std::find(groups.begin(), groups.end(), i);
        if (ptr == groups.end()) {
          std::cout << "Remaining groups do not have (" << i[0] << " " << i[1] << " " << i[2]  << ")? \n";
          for (const auto & g: groups) std::cout << "(" << g[0] << " " << g[1] << " " << g[2]  << ")\n";
          throw std::runtime_error("triangle indexing error");
        }
        std::cout << "Selected is " << ptr->operator[](0) << " " << ptr->operator[](1) << " " << ptr->operator[](2) << " with area ratio " << ear_measure(*ptr) << "\n";
        // construct the triangle
        auto tri = Wire(*ptr);
        triangles.push_back(tri);
        // connect the neighbouring groups
        auto p{*(ptr == groups.begin() ? groups.end() - 1 : ptr - 1)};
        auto n{*(ptr == groups.end() - 1 ? groups.begin() : ptr + 1)};
        group_t np{{p[0], p[1], n[1]}};
        group_t nn{{p[1], n[1], n[2]}};
//        std::cout << "Update {prev|next} from {";
//        std::cout << p[0] << " " << p[1] << " " << p[2] << "|" << n[0] << " " << n[1] << " " << n[2] << "} to {";
//        std::cout << np[0] << " " << np[1] << " " << np[2] << "|" << nn[0] << " " << nn[1] << " " << nn[2] << "}\n";

        // check for change of concave-ness, updating the concave and convex lists
        update_concave_convex(p, np);
        update_concave_convex(n, nn);
        // update the inner points subset of vertices
        middle.clear();
        middle.reserve(concave.size());
        for (const auto & g: concave) middle.push_back(g[1]);
        inner_points = x.extract(middle);
        // (possibly) update the ears
        update_ears(inner_points, p, np);
        update_ears(inner_points, n, nn);
        // update the groups:
        update_groups(p, np);
        update_groups(n, nn);
        // the update calls above may have invalidated pointers, so find the group again
        ptr = std::find(groups.begin(), groups.end(), i);
        groups.erase(ptr);

      }
      return Network<Wire,T,A>(triangles, x);
    }

    template<class T, template<class> class A>
    T area(const A<T> &x) const {
      if (empty()) return 0;
      auto max_ptr = std::max_element(begin(), end());
      if (*max_ptr >= x.size(0u)) {
        debug_update("Wire points\n", x.to_string(),"and wire", *this);
        throw std::runtime_error("Wire::area attempting to access out-of-bounds element of Array2");
      }
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
      std::copy(crbegin(), crend(), std::back_insert_iterator<base_t>(out));
      return out;
    }

    Wire inverse() const {
      return mirror();
    }

    template<class T, template<class> class A>
      A<double> find_contained_point(const A<T>& x) const {
      // Return *a* point within the closed border of the polygon wire
      for (size_t i=0; i<size()-1; ++i){
        auto pi = x.view(vi(i));
        // check the vertices not next to this one
        for (size_t j=0; j<size(); ++j) if (i > j ? i - j > 1 : j - i > 1) {
          // make sure there isn't a crossing we should worry about
          if (std::count(begin(), end(), vi(j)) == 1){
            auto mid = 0.5 * (pi + x.view(vi(j)));
            if (contains(mid, x)[0]) return mid;
          }
        }
      }
      // No found contained point. This seems very unlikely, but I don't know what to do so punt
      return 1.0*x;
    }

    template<class T, class S, template<class> class A>
    std::vector<bool> contains(const A<S> & point, const A<T> & x) const{
      auto wn = winding_number(point, x);
      std::vector<bool> out;
      out.reserve(wn.size());
      std::transform(wn.begin(), wn.end(), std::back_inserter(out), [](const auto & x){return x!=0;});
      return out;
    }

    template<class T, class R, template<class> class A, class S=std::common_type_t<T,R>>
    std::vector<size_t> winding_number(const A<R> & point, const A<T> & x) const {
      std::vector<size_t> out;
      out.reserve(point.size(0));
      std::array<ind_t, 2> y_coord {0, 1};
      for (ind_t i=0; i<point.size(0); ++i) {
        size_t wn{0};
        auto point_i = point.view(i);
        auto point_y = static_cast<S>(point_i.val(y_coord));
        for (size_t j=0; j<size(); ++j){
          auto e = edge(j);
          auto first = x.view(e.first);
          auto second = x.view(e.second);
          auto first_y = static_cast<S>(first.val(y_coord));
          auto second_y = static_cast<S>(second.val(y_coord));
          if (first_y <= point_y) {
            if(second_y > point_y) {
              // upwards crossing with infinite horizontal line from point
              if (cross2d(second - first, point_i - first).sum() > 0) ++wn;
            }
          } else {
            if (second_y <= point_y) {
              if (cross2d(second - first, point_i - first).sum() < 0) --wn;
            }
          }
        }
        out.push_back(wn);
      }
      return out;
    }

    template<class T, template<class> class A>
    T wire_winding_number_over_pi(const A<T> & vertices) const {
      auto edge_dir = [&](ind_t i){
        using polystar::math::pi;
        auto e = edge(i);
        auto d = vertices.view(e.second) - vertices.view(e.first);
        return std::atan2(d.val(0, 1),  d.val(0, 0)) / pi;
      };
      T no{0};
      std::optional<T> orientation;
      auto dir = edge_dir(size()-1); // starting direction is given by _last_ edge
      for (ind_t i=0; i<size(); ++i){
        auto old_dir = dir;
        dir = edge_dir(i);
        auto angle = dir - old_dir;
        debug_update("i=",i," dir=",dir," old_dir=",old_dir," angle=",angle);
        if (angle <= -1){
          angle += 2;
        } else if (angle > 1){
          angle -= 2;
        }
        if (orientation.has_value()){
          if (orientation.value() * angle <= 0.) return 0;
        } else {
          if (angle == 0.) return 0;
          orientation = angle > 0 ? 1. : -1.;
        }
        no += angle;
        debug_update("->angle=",angle," no=",no);
      }
      return no;
    }



    template<class T, class R, template<class> class A>
    std::vector<bool> is_on(const A<R> & point, const A<T> & x) const {
      std::vector<bool> out;
      out.reserve(point.size(0));
      for (ind_t i=0; i<point.size(0); ++i){
        auto pi = point.view(i);
        bool on_edge{false};
        for (size_t j=0; j<size(); ++j){
          auto e = edge(j);
          auto a = x.view(e.first);
          auto v0 = pi - a;
          auto v1 = x.view(e.second) - a;
          auto v0n = norm(v0).val(0, 0);
          auto v1n = norm(v1).val(0, 0);
          if (v0n == 0 || (cross2d(v0, v1).sum() == 0 && dot(v0, v1).val(0, 0) > 0  && v0n <= v1n)){
            on_edge=true;
            break;
          }
        }
        out.push_back(on_edge);
      }
      return out;
    }

    template<class T, template<class> class A>
    std::vector<size_t> crossing_number(const A<T> &point, const A<T> &x, end_type inclusive = end_type::second, bool verbose=false) const {
      std::vector<size_t> out;
      out.reserve(point.size(0));
      auto x_max = 2 * x.max(0);
      edge_t segment_edge{0, 1};
      for (ind_t i=0; i<point.size(0); ++i){
        auto segment = cat(0, point.view(i), x_max);
        size_t crossings{0};
        for (size_t j=0; j<size(); ++j) if (intersect2d(segment, segment_edge, x, edge(j), inclusive, verbose)) ++crossings;
        out.push_back(crossings);
      }
      return out;
    }

    template<class T, template<class> class A>
    bool intersects(const A<T> &other, const edge_t &oe, const A<T> &x, end_type inclusive = end_type::both) const {
      for (size_t i = 0; i < size(); ++i) if (intersect2d(other, oe, x, edge(i), inclusive)) return true;
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
      for (const auto &i: *this) if (auto t = norm(c - x.view(i)).val(0,0); t > r) r = t;
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

//    template<class T, class R, template<class> class A,
//      template<class> class B>
//    std::tuple <A<T>, Wire> cut(const A<T> &x, const B<R> &a, const B<R> &b) const {
////      auto is_b = is_beyond(x, a, b);
////      A<T> xo(0u, 2u);
////      if (!std::any_of(is_b.begin(), is_b.end(), true)) return std::make_tuple(x, Wire(*this));
////      if (!std::all_of(is_b.begin(), is_b.end(), true)) return std::make_tuple(A<T>(), Wire());
////      // otherwise, go-around the edges and find intersections
////      auto ov = cat(0, a, b);
////      edge_t oe{0,1};
////      for (size_t i=0; i<size(); ++i){
////        auto [n, p] = intersection2d(ov, oe, x, edge(i));
////        if (n == 1) {
////
////        }
////      }
//      return std::make_tuple(x, Wire(*this));
//    }

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

    template<class T, template<class> class A>
    void add_to_svg(SVG::Path & path, const A<T>& x) const {
      // don't attempt to write a non-closed wire path
      if (size() < 3) return;
      auto first_and_last = vi(0);
      if (first_and_last >= x.size(0)) {
        std::cout << "Wire to path conversion failure due to out-of-bounds vertex index " << first_and_last
                  << " for length " << x.size(0) << " vertex vector\n";
        return;
      }
      path.move_to(x.val(first_and_last, 0), x.val(first_and_last, 1));
      size_t skip=0;
      for (const auto & i: *this) if (skip++) {
        if (i < x.size(0)){
          path.line_to(x.val(i, 0), x.val(i, 1));
        } else {
          std::cout << "Wire to path conversion failure due to out-of-bounds vertex index " << i
                    << " for length " << x.size(0) << " vertex vector\n";
        }
      }
      // close the wire
      path.line_to(x.val(first_and_last, 0), x.val(first_and_last, 1));
    }

  };

  Wire wire_merge(const Wire & a, const Wire & b);

} // polystar::polygon
#endif