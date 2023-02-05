#ifndef POLYSTAR_POLYGON_HPP
#define POLYSTAR_POLYGON_HPP
#include <vector>
#include <array>
#include <optional>
#include <numeric>
#include <limits>

#include "comparisons.hpp"
#include "array_.hpp"
#include "geometry.hpp"
#include "approx_float.hpp"

namespace polystar::polygon{
  using ind_t = polystar::ind_t;

  class Wire : public std::vector<ind_t> {
  public:
    using edge_t = std::pair<ind_t, ind_t>;
    using base_t = std::vector<ind_t>;
  private:
    inline ind_t vi(size_t i) const {return operator[](i);}
    inline ind_t vj(size_t i) const {return operator[]((i+1) % size());}
  public:
    inline edge_t edge(size_t i) const {return std::make_pair(vi(i), vj(i));}

    // apparently the empty constructor is needed for explicitness?
    explicit Wire() {}

    // type conversion (require that T is integer?)
    template<class T>
    explicit Wire(const std::vector<T>& indexes) {
      reserve(indexes.size());
      std::transform(indexes.begin(), indexes.end(),
                     std::back_insert_iterator<std::vector<ind_t>>(*this),
                     [](const T & x){return static_cast<ind_t>(x);});
    }

    template<class T, template<class> class A, class U = typename A<T>::shape_t>
      explicit Wire(const A<T>& vertices, const T tol=T(0), const int dig=1){
      // Convex Hull construction
      auto unique = vertices.is_unique(tol, dig);
      auto v = vertices.extract(unique);
      auto chi= find_convex_hull_edges(v, tol, dig);
      if (!chi.empty()){
        // the wire should index vertices but indexes the unique vertices instead at the moment:
        // the Convex Hull algorithm only returned strictly-necessary-point indexes, so no extra work required
        std::vector<ind_t> index_map(v.size(0), vertices.size(0));
        for (ind_t i=0; i<v.size(0); ++i) index_map[i] = vertices.first(cmp::eq, v.view(i));
        reserve(chi.size());
        for (auto & index: chi) push_back(index_map[index]);
      }
    }

    std::back_insert_iterator<base_t> appender() {return std::back_insert_iterator<base_t>(*this);}

    template<class T, template<class> class A>
    T area(const A<T>& x) const {
      T a{0};
      for (size_t i=0; i<size(); ++i) a += cross2d(x.view(vi(i)), x.view(vj(i))).val(0, 0);
      return a / T(2);
    }
    template<class T, template<class> class A>
    A<T> centroid(const A<T>& x) const {
      auto c = 0 * x.view(0);
      auto twice_area = T(0);
      for (size_t i=0; i<size(); ++i) {
        const auto & a{x.view(vi(i))};
        const auto & b{x.view(vj(i))};
        auto ta = cross2d(a, b).val(0);
        c += (a+b) * ta;
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
    bool contains(const A<T>& point, const A<T>& x) const {
      auto segment = cat(0, point, (std::numeric_limits<T>::max)() + T(0) * point);
      edge_t se{0,1};
      size_t crossings{0};
      for (size_t i=0; i<size(); ++i) if (intersect2d(segment, se, x, edge(i))) ++crossings;
      return (crossings % 2u) == 1u;
    }

    template<class T, template<class> class A>
    bool intersects(const A<T>& other, const edge_t& oe, const A<T>& x) const {
      for (size_t i=0; i<size(); ++i) if (intersect2d(other, oe, x, edge(i))) return true;
      return false;
    }

    template<class T, template<class> class A>
    A<T> intersections(const A<T> & other, const edge_t & oe, const A<T> & x) const {
      size_t count{0};
      auto out = A<T>();
      for (size_t i=0; i<size(); ++i) {
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
    Wire replace(const std::vector<T>& replacement) const {
      Wire p;
      p.reserve(this->size());
      std::transform(this->cbegin(), this->cend(), std::back_insert_iterator<std::vector<ind_t>>(p),
                     [&replacement](const T &index){return replacement[index];});
      return p;
    }
    template<class T, template<class> class A>
    [[nodiscard]] T circumscribed_radius(const A<T>& x) const {
      auto c = centroid(x);
      T r{0};
      for (const auto& i: *this) if (auto t = norm(c - x.view(i)); t > r) r = t;
      return r;
    }

    Wire offset(const ind_t value) const {
      Wire o;
      o.reserve(size());
      std::transform(cbegin(), cend(), o.appender(), [value](const ind_t & i){return i + value;});
      return o;
    }
    Wire combine(const Wire & that, const ind_t value) const {
      Wire o;
      o.reserve(size() + that.size());
      std::copy(cbegin(), cend(), o.appender());
      auto t = that.offset(value);
      std::copy(t.cbegin(), t.cend(), o.appender());
      return o;
    }
    template<class T, class R, template<class> class A, template<class> class B>
    std::tuple<A<T>, Wire> cut(const A<T> & x, const B<R>& a, const B<R>& b) const {
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
    template<class T, class R, template<class> class A, template<class> class B>
    bool any_beyond(const A<T> & x, const B<R>& a, const B<R>& b) const {
      assert(a.size(0) == b.size(0) && a.size(0) == 1);
      auto v = b - a;
      for (const auto & i: *this) if (cross2d(v, x.view(i)).val(0, 0) < 0) return true;
      return false;
    }
    template<class T, class R, template<class> class A, template<class> class B>
    std::vector<bool> is_beyond(const A<T> & x, const B<R>& a, const B<R>& b) const {
      assert(a.size(0) == b.size(0) && a.size(0) == 1);
      auto v = b - a;
      std::vector<bool> out;
      out.reserve(size());
      std::transform(cbegin(), cend(), std::back_inserter(out), [&x, &v](const ind_t i){return cross2d(v, x.view(i)).val(0,0) < 0;});
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
    friend std::ostream & operator<<(std::ostream& os, Wire w) {
      os << w.python_string();
      return os;
    }

  };

  class Wires{
  public:
    using wire_t = Wire;
    using proto_t = typename std::vector<wire_t>;
  protected:
    wire_t border_;
    std::optional<proto_t> wires_ = std::nullopt;
  public:
    explicit Wires(): border_(), wires_() {}
    explicit Wires(wire_t b): border_(std::move(b)) {}
    template<class T>
    explicit Wires(const std::vector<T>& b): border_(wire_t(b)) {}
    explicit Wires(wire_t b, proto_t w): border_(std::move(b)), wires_(std::move(w)) {}
    template<class T>
    explicit Wires(const std::vector<T>& b, const std::vector<std::vector<T>>& w): border_(wire_t(b)) {
      proto_t w_;
      w_.reserve(w.size());
      for (const auto & x: w) w_.push_back(wire_t(x));
      wires_ = w_;
    }
    template<class T, template<class> class A, class U = typename A<T>::shape_t>
      explicit Wires(const A<T>& vertices, const T tol=T(0), const int dig=1): border_(wire_t(vertices, tol, dig)){}

    template<class T, template<class> class A, class U = typename A<T>::shape_t>
      Wires(wire_t border, proto_t wires, const A<T>& vertices){
      border_=std::move(border);
      auto v = T(1) * vertices;
      proto_t w;
      std::tie(v, border_, w) = remove_duplicate_points_and_update_wire_indexing(v, border_, wires);
      std::tie(border_, w) = polygon_border_wires(border_, w, v);
      if (polygon_edge_vertex_perge(v, border_, w)){
        std::vector<ind_t> map(v.size(0), vertices.size(0));
        for (ind_t i=0; i<v.size(0); ++i){
          map[i] = vertices.first(polystar::cmp::eq, v.view(i));
        }
        for (auto & idx: border_) idx = map[idx];
        for (auto & wire: w) for (auto & idx: wire) idx = map[idx];
        wires_ = w;
      }
    }

    Wires(const Wires & that) = default;
    Wires(Wires && that) noexcept: border_(std::move(that.border_)), wires_(std::move(that.wires_)) {}
    Wires& operator=(const Wires& that)= default;

    bool operator!=(const Wires&) const;
    bool operator==(const Wires& that) const {return !this->operator!=(that);}

    template<class I> Wires permute(const std::vector<I>& permutation) const {
      auto border = border_.replace(permutation);
      if (wires_.has_value()){
        proto_t wires;
        wires.reserve(wires_.value().size());
        for (const auto & wire: wires_.value()) wires.push_back(wire.replace(permutation));
        return Wires(border, wires);
      }
      return Wires(border);
    }

    // direct property accessors
    [[nodiscard]] size_t border_size() const {return border_.size();}
    [[nodiscard]] wire_t border() const {return border_;}
    [[nodiscard]] size_t wire_count() const {return wires_.has_value() ? wires_.value().size() : 0;}
    [[nodiscard]] proto_t wires() const {return wires_.has_value() ? wires_.value() : proto_t();}
    [[nodiscard]] wire_t wire(const ind_t& i) const {
      assert(wires_.has_value());
      assert(i < wire_count());
      return wires_.value()[i];
    }
    // calculated property accessors
    [[nodiscard]] Wires mirror() const {
      auto b = border_.mirror();
      proto_t w;
      w.reserve(wires_.has_value() ? wires_.value().size() : 0u);
      if (wires_.has_value()) for (const auto & w_: wires_.value()) w.push_back(w_.mirror());
      return wires_.has_value() ? Wires(b, w) : Wires(b);
    }
    [[nodiscard]] wire_t::base_t indexes() const {
      wire_t::base_t all;
      auto not_in_all = [&all](const auto &i){return std::find(all.begin(), all.end(), i) == all.end();};
      for (const auto & x: border_) if (not_in_all(x)) all.push_back(x);
      if (wires_.has_value()) for (const auto & wire: wires_.value()) for (const auto & x: wire) if (not_in_all(x)) all.push_back(x);
      return all;
    }
    [[nodiscard]] polystar::Array2<ind_t> border_edges() const {
      auto n = border_.size();
      polystar::Array2<ind_t> be(static_cast<ind_t>(n+1), 2u);
      for (ind_t i=0; i<n; ++i) {
        auto edge = border_.edge(i);
        be.val(i, 0) = edge.first;
        be.val(i, 1) = edge.second;
      }
      return be;
    }
    template<class T, template<class> class A, class U = typename A<T>::shape_t>
    [[nodiscard]] std::enable_if_t<polystar::isArray<T,A>, std::tuple<A<T>, A<T>>>
    borders(const A<T>& x) const {
      auto n = static_cast<ind_t>(border_.size());
      auto a = 0 * x.view(0); a.resize(n);
      auto b = 0 * x.view(0); b.resize(n);
      for (size_t i=0; i < border_.size(); ++i){
        auto e = border_.edge(i);
        a.set(i, x.view(e.first));
        b.set(i, x.view(e.second));
      }
      return std::make_tuple(a, b);
    }
    template<class T, template<class> class A, class U = typename A<T>::shape_t>
    [[nodiscard]] T area(const A<T>& x) const {
      auto result = border_.area(x);
      // result should be positive, and any wires should be holes which reduce the area
      if (wires_.has_value()) for (const auto & w: wires_.value()) result += w.area(x);
      return result;
    }
    template<class T, template<class> class A, class U = typename A<T>::shape_t>
      [[nodiscard]] A<T> centroid(const A<T>& x) const {
      auto result = border_.centroid(x);
      if (wires_.has_value()) for (const auto & w: wires_.value()) result += w.centroid(x);
      return result;
    }
    template<class T, template<class> class A, class U = typename A<T>::shape_t>
      [[nodiscard]] T circumscribed_radius(const A<T>& x) const {
      // only the *border* contributes to the circumscribed radius; holes inside would only shift the centroid
      return border_.circumscribed_radius(x);
    }

    Wires combine(const Wires & that, const ind_t offset) const {
      auto b = border_.combine(that.border(), offset);
      auto w = wires_.has_value() ? wires_.value() : proto_t();
      w.reserve(w.size() + that.wire_count());
      if (that.wire_count()) for (const auto & tw: that.wires()) w.push_back(tw.offset(offset));
      return w.empty() ? Wires(b) : Wires(b, w);
    }

    template<class T, class R, template<class> class A, template<class> class B>
    std::tuple<A<T>, Wires> cut(const A<T>& x, const B<R>& a, const B<R>& b, const R tol=R(0), const int dig=1) const {
      assert(a.ndim() ==2 && b.ndim() == 2);
      assert(a.size(1) == 2 && b.size(1) == 2);
      assert(a.size(0) == b.size(0));
      auto ov = (T(1) * x).decouple();
      Wires ow(*this);
      for (ind_t i=0; i<a.size(0); ++i) std::tie(ov, ow) = ow.one_cut(ov, a.view(i), b.view(i), tol, dig);
      return std::make_tuple(ov, ow);
    }
    template<class T, class R, template<class> class A, template<class> class B>
    std::tuple<A<T>, Wires> one_cut(const A<T>& x, const B<R>& a, const B<R>& b, const R tol=R(0), const int dig=1) const {
      assert(a.ndim() ==2 && b.ndim() == 2);
      assert(a.size(1) == 2 && b.size(1) == 2);
      assert(a.size(0) == b.size(0) && a.size(0) == 1);
      auto ov = (T(1) * x).decouple();
      if (none_beyond(x, a, b)) return std::make_tuple(x, Wires());
      // otherwise find the vertices beyond a->b and remove them

      return std::make_tuple(x, Wires(*this));
    }

    template<class T, class R, template<class> class A, template<class> class B>
    bool none_beyond(const A<T>& x, const B<R>& a, const B<R>& b) const {
      if (border_.any_beyond(x, a, b)) return false;
      if (wires_.has_value()) for (const auto & w: wires_.value()) if (w.any_beyond(x, a, b)) return false;
      return true;
    }

    std::string python_string() const {
      std::stringstream s;
      s << border_.python_string();
      s << ",[";
      size_t count{0};
      if (wires_.has_value()) for (const auto & w: wires_.value()) s << w.python_string() << (++count < wires_.value().size() ? "," : "");
      s << "]";
      return s.str();
    }
    friend std::ostream & operator<<(std::ostream & os, Wires w) {
      os << w.python_string();
      return os;
    }
    template<class T, template<class> class A>
    bool contains(const A<T>& point, const A<T>& x) const {
      auto segment = cat(0, point, (std::numeric_limits<T>::max)() + T(0) * point);
      wire_t::edge_t se{0,1};
      size_t crossings{0};
      for (size_t i=0; i<border_.size(); ++i) if (intersect2d(segment, se, x, border_.edge(i))) ++crossings;
      if (wires_.has_value()) for (const auto & w: wires_.value())
        for (size_t i=0; i<w.size(); ++i) if (intersect2d(segment, se, x, w.edge(i))) ++crossings;
      return (crossings % 2u) == 1u;
    }

#ifdef USE_HIGHFIVE
    template<class H> std::enable_if_t<std::is_base_of_v<HighFive::Object, H>, bool>
      to_hdf(H& obj, const std::string & entry) const {
      auto group = overwrite_group(obj, entry);
      bool ok{true};
      ok &= list_to_hdf<ind_t>(border_, group, "border");
      std::vector<std::vector<ind_t>> copied_wires; // not great, but vector<Wire> can´t be figured-out by the compiler
      if (wires_.has_value()) for (const auto & w: wires_.value()) copied_wires.push_back(w);
      ok &= lists_to_hdf(copied_wires, group, "wires");
      return ok;
    }
    [[nodiscard]] bool to_hdf(const std::string & filename, const std::string & dataset, unsigned perm=HighFive::File::OpenOrCreate) const {
      HighFive::File file(filename, perm);
      return to_hdf(file, dataset);
    }
    template<class H> static std::enable_if_t<std::is_base_of_v<HighFive::Object, H>, Wires>
      from_hdf(H& obj, const std::string& entry){
      auto group = obj.getGroup(entry);
      auto b = list_from_hdf<ind_t>(group, "border");
      auto w = lists_from_hdf<ind_t>(group, "wires");
      return w.empty() ? Wires(b) : Wires(b, w);
    }
    static Wires from_hdf(const std::string& filename, const std::string& dataset){
      HighFive::File file(filename, HighFive::File::ReadOnly);
      return Wires::from_hdf(file, dataset);
    }
#endif
  };

  template<class T, template<class> class A>
  class Poly {
  public:
    using vertex_t = A<T>;
    using wires_t = Wires;
  protected:
    vertex_t vertices_;
    wires_t wires_;
  public:
    explicit Poly(): vertices_(), wires_() {}
    // Convex Hull constructors
    explicit Poly(const vertex_t & v, T tol=T(0), int dig=1): vertices_(v), wires_(v, tol, dig) {
      finish_convex_hull(tol, dig);
    }
    explicit Poly(vertex_t && v, T tol=T(0), int dig=1): vertices_(v), wires_(v, tol, dig){
      finish_convex_hull(tol, dig);
    }
    // Simple polygon constructors
    Poly(const vertex_t & v, wires_t::wire_t w): vertices_(v), wires_(w) {}
    Poly(vertex_t && v, wires_t::wire_t && w): vertices_(std::move(v)), wires_(std::move(w)) {}
    // Complex polygon constructors
    Poly(const vertex_t & v, wires_t::wire_t b, const wires_t::proto_t & w): vertices_(v), wires_(b, w) {}
    Poly(vertex_t && v, wires_t::wire_t && b, wires_t::proto_t && w): vertices_(std::move(v)), wires_(std::move(b), std::move(w)) {}
    Poly(const vertex_t & v, const wires_t & w): vertices_(v), wires_(w) {}
    Poly(vertex_t && v, wires_t && w): vertices_(std::move(v)), wires_(std::move(w)) {}
    // Copy constructor
    Poly(const Poly<T,A> & that): vertices_(that.vertices_), wires_(that.wires_) {}
    Poly(Poly<T,A> & that) noexcept: vertices_(that.vertices_), wires_(std::move(that.wires_)) {}
    // Copy assignment
    Poly<T,A> & operator=(const Poly<T,A> & that){
      vertices_ = that.vertices_;
      wires_ = that.wires_;
      return *this;
    }
    Poly<T,A> & operator=(Poly<T,A> && that){
      vertices_ = that.vertices_;
      wires_ = std::move(that.wires_);
      return *this;
    }
    // direct property accessor
    [[nodiscard]] ind_t vertex_count() const {return vertices_.size(0);}
    [[nodiscard]] size_t face_count() const {return 1u;}
    [[nodiscard]] vertex_t vertices() const {return vertices_;}
    [[nodiscard]] wires_t wires() const {return wires_;}
    // calculated property accessors
    [[nodiscard]] T area() const {return wires_.area(vertices_);}
    [[nodiscard]] vertex_t centroid() const {return wires_.centroid(vertices_);}
    [[nodiscard]] T circumscribed_radius() const {return wires_.circumscribed_radius(vertices_);}

    // methods
    [[nodiscard]] Poly<T,A> convex_hull() const {return Poly(vertices_);}
    [[nodiscard]] bool is_not_approx(const Poly<T,A> & that, const T tol=T(0), const int dig=1) const {
      bool permuted{false};
      if (vertices_ != that.vertices){
        permuted = vertices_.is_permutation(that.vertices_, tol, tol, dig);
        if (!permuted) return true;
      }
      if (permuted){
        auto permutation = vertices_.permutation_vector(that.vertices_, tol, tol, dig);
        auto permuted_wires = wires_.permute(permutation);
        return permuted_wires != that.wires_;
      }
      return wires_ != that.wires_;
    }
    [[nodiscard]] bool is_approx(const Poly<T,A>& that, const T tol=T(0), const int dig=1) const {
      return !is_not_approx(that, tol, dig);
    }
    [[nodiscard]] bool operator!=(const Poly<T,A>& that) const {return is_not_approx(that);}
    [[nodiscard]] bool operator==(const Poly<T,A>& that) const {return !is_not_approx(that);}

    Poly<T,A> operator+(const Poly<T,A>& that) const {return combine(that);}
    Poly<T,A> combine(const Poly<T,A>& that, const T tol=T(0), const int dig=1) const {
      // combine vertices
      auto v = cat(0, vertices_, that.vertices_);
      // combined wires
      auto w = wires_.combine(that.wires_, vertices_.size(0)).wires();
      // check for duplicate vertices
      std::tie(v, w) = remove_duplicate_points_and_update_wire_indexing(v, w, tol, dig);
      // look for overlapping wires now that vertex indexing has been modified
      std::vector<bool> needed(w.size(), true);
      for (ind_t i=0; i < w.size()-1; ++i) if (needed[i]) {
        const auto & a{w[i]};
        for (ind_t j=i+1; j < w.size(); ++j) if (needed[j]) {
          const auto & b{w[j]};
          auto c = std::count_if(a.begin(), a.end(), [&b](const auto & z){
            return std::count(b.begin(), b.end(), z) > 0;
          });
          if (c && std::is_permutation(a.begin(), a.end(), b.begin())) {
            needed[i] = is_positive_permutation(a, b);
            needed[j] = false;
          }
        }
      }
      if (std::find(needed.begin(), needed.end(), false) != needed.end()) {
        for (ind_t i=0; i<w.size(); ++i) if (!needed[i]) w[i].clear();
        w.erase(std::remove_if(w.begin(), w.end(), [](const auto & x){return x.empty();}), w.end());
      }
      return {v, w};
    }

    Poly<T,A> mirror() const {return {T(-1) * vertices_, wires_.mirror()};}
    Poly<T,A> centre() const {return {vertices_ - centroid(), wires_};}
    Poly<T,A> translate(const A<T>& v) const {return {vertices_ + v, wires_};}

    template<class R, template<class> class B, class... P>
      [[nodiscard]] std::vector<bool> contains(const B<R>& x, P... p) const {
      return wires_.contains(vertices_, x, p...);
    }
    template<class R, class... P>
      [[nodiscard]] std::vector<bool> contains(const std::vector<std::array<R,2>> & x, P... p) const {
      return contains(from_std_like(vertices_, x), p...);
    }
    template<class R, template<class> class B>
      [[nodiscard]] std::enable_if_t<isArray<R, B>, bool>
      intersects(const Poly<R, B> & that, const R tol=R(0), const int dig=1) const {
      auto overlap = intersection(that, tol, dig);
      if (!approx_float::scalar(overlap.area() / (area() + overlap.area()), R(0), tol, tol, dig)) {
        return true;
      }
      return false;
    }
    template<class R, template<class> class B>
      [[nodiscard]] std::enable_if_t<isArray<R,B>, Poly<T,A>> intersection(const Poly<R,B>& that, const R tol=R(0), const int dig=1) const {
      auto [a, b] = that.edges();
      return cut(a, b, tol, dig);
    }
    template<class R, template<class> class B>
    [[nodiscard]] std::enable_if_t<isArray<R,B>, Poly<T,A>> cut(const B<R>& a, const B<R>& b, const R tol=R(0), const int dig=1) const {
      auto [v, w] = wires_.cut(vertices_, a, b, tol, dig);
      return {v, w};
    }
//    template<class R, template<class> class B>
//    [[nodiscard]] std::enable_if_t<isArray<R,B>, size_t> edge_index(const B<R>& a, const B<R>& b) const {
//      return wires_.edge_index(vertices_, a, b);
//    }
//    template<class R, template<class> class B>
//    [[nodiscard]] std::enable_if_t<isArray<R,B>, bool> has_edge(const B<R>& a, const B<R>& b) const {
//      return wires_.has_edge(vertices_, a, b);
//    }
    template<class R, template<class> class B>
    [[nodiscard]] std::enable_if_t<isArray<R,B>, bool> none_beyond(const B<R>& a, const B<R>& b) const {
      return wires_.none_beyond(vertices_, a, b);
    }

    [[nodiscard]] std::string python_string() const {
      return "np.array(" + get_xyz(vertices_).to_string()+"), " + wires_.python_string();
    }
    friend std::ostream & operator<<(std::ostream & os, const Poly<T,A>& p){
      os << p.python_string();
      return os;
    }

#ifdef USE_HIGHFIVE
    template<class H> std::enable_if_t<std::is_base_of_v<HighFive::Object, H>, bool>
      to_hdf(H& obj, const std::string & entry) const {
      auto group = overwrite_group(obj, entry);
      bool ok{true};
      ok &= vertices_.to_hdf(group, "vertices");
      ok &= wires_.to_hdf(group, "wires");
      return ok;
    }
    [[nodiscard]] bool to_hdf(const std::string & filename, const std::string & dataset, unsigned perm=HighFive::File::OpenOrCreate) const {
      HighFive::File file(filename, perm);
      return to_hdf(file, dataset);
    }
    template<class H> static std::enable_if_t<std::is_base_of_v<HighFive::Object, H>, Poly<T,A>>
      from_hdf(H& obj, const std::string& entry){
      auto group = obj.getGroup(entry);
      auto v = vertex_t::from_hdf(group, "vertices");
      auto w = wires_t::from_hdf(group, "wires");
      return Poly<T,A>(v, w);
    }
    static Poly<T,A> from_hdf(const std::string& filename, const std::string& dataset){
      HighFive::File file(filename, HighFive::File::ReadOnly);
      return Poly<T,A>::from_hdf(file, dataset);
    }
#endif

  private:
    void finish_convex_hull(const T, const int){
      // check for unused vertices and remove them
      auto present = wires_.indexes(); // unordered list of wires_ indexes present/used
      std::sort(present.begin(), present.end());
      std::vector<ind_t> indexes(vertices_.size(0));
      std::iota(indexes.begin(), indexes.end(), 0u);
      if (std::includes(present.begin(), present.end(), indexes.begin(), indexes.end())) return;

      // make map from old to new indexes, and an extraction list
      std::vector<bool> keep(indexes.size(), false);
      std::fill(indexes.begin(), indexes.end(), indexes.size());
      ind_t kept{0};
      for (const auto & x: present){
        indexes[x] = kept++;
        keep[x] = true;
      }
      // update wire vectors
      wires_ = wires_.permute(indexes); // I think this is the same as replacing all values in Wires with indexes[values]
      vertices_ = vertices_.extract(keep);
    }
  };

  template<class T, template<class> class A>
    std::enable_if_t<isArray<T,A>, Poly<T,A>>
    bounding_box(const A<T>& points){
      auto min = get_xyz(points).min(0);
      auto max = get_xyz(points).max(0);
      std::vector<std::array<T,2>> v{
        {min[{0, 0}], min[{0, 1}]}, // 00
        {max[{0, 0}], min[{0, 1}]}, // 10
        {max[{0, 0}], max[{0, 1}]}, // 11
        {min[{0, 0}], max[{0, 1}]}, // 01
      };
      auto wires = typename Poly<T,A>::wires_t({{0,1,2,3}});
      auto vert = from_xyz_like(points, bArray<T>::from_std(v));
      return {vert, wires};
    }
}

#endif