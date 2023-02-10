#ifndef POLYSTAR_POLYGON_WIRES_HPP
#define POLYSTAR_POLYGON_WIRES_HPP
#include <vector>
#include <array>
#include <optional>
#include <numeric>
#include <limits>

#include "comparisons.hpp"
#include "array_.hpp"
#include "geometry.hpp"
#include "approx_float.hpp"
#include "polygon_wire.hpp"

namespace polystar::polygon {
  using ind_t = polystar::ind_t;

  class Wires {
  public:
    using wire_t = Wire;
    using proto_t = typename std::vector<wire_t>;
  protected:
    wire_t border_;
    std::optional <proto_t> wires_ = std::nullopt;
  public:
    explicit Wires() : border_(), wires_() {}

    explicit Wires(wire_t b) : border_(std::move(b)) {}

    template<class T>
    explicit Wires(const std::vector <T> &b): border_(wire_t(b)) {}

    explicit Wires(wire_t b, proto_t w) : border_(std::move(b)), wires_(std::move(w)) {}

    template<class T>
    explicit Wires(const std::vector <T> &b, const std::vector <std::vector<T>> &w): border_(wire_t(b)) {
      proto_t w_;
      w_.reserve(w.size());
      for (const auto &x: w) w_.push_back(wire_t(x));
      wires_ = w_;
    }

    template<class T, template<class> class A, class U = typename A<T>::shape_t>
    explicit Wires(const A<T> &vertices, const T tol = T(0), const int dig = 1): border_(wire_t(vertices, tol, dig)) {}

    template<class T, template<class> class A, class U = typename A<T>::shape_t>
    Wires(wire_t border, proto_t wires, const A<T> &vertices) {
      border_ = std::move(border);
      auto v = T(1) * vertices;
      proto_t w;
      std::tie(v, border_, w) = remove_duplicate_points_and_update_wire_indexing(v, border_, wires);
      std::tie(border_, w) = polygon_border_wires(border_, w, v);
      if (polygon_edge_vertex_perge(v, border_, w)) {
        std::vector <ind_t> map(v.size(0), vertices.size(0));
        for (ind_t i = 0; i < v.size(0); ++i) {
          map[i] = vertices.first(polystar::cmp::eq, v.view(i));
        }
        for (auto &idx: border_) idx = map[idx];
        for (auto &wire: w) for (auto &idx: wire) idx = map[idx];
        wires_ = w;
      }
    }

    Wires(const Wires &that) = default;

    Wires(Wires &&that)

    noexcept: border_(std::move(that.border_)), wires_(std::move(that
    .wires_)) {}

    Wires &operator=(const Wires &that) = default;

    bool operator!=(const Wires &) const;

    bool operator==(const Wires &that) const { return !this->operator!=(that); }

    template<class I>
    Wires permute(const std::vector <I> &permutation) const {
      auto border = border_.replace(permutation);
      if (wires_.has_value()) {
        proto_t wires;
        wires.reserve(wires_.value().size());
        for (const auto &wire: wires_.value()) wires.push_back(wire.replace(permutation));
        return Wires(border, wires);
      }
      return Wires(border);
    }

    // direct property accessors
    [[nodiscard]] size_t border_size() const { return border_.size(); }

    [[nodiscard]] wire_t border() const { return border_; }

    [[nodiscard]] size_t wire_count() const { return wires_.has_value() ? wires_.value().size() : 0; }

    [[nodiscard]] proto_t wires() const { return wires_.has_value() ? wires_.value() : proto_t(); }

    [[nodiscard]] wire_t wire(const ind_t &i) const {
      assert(wires_.has_value());
      assert(i < wire_count());
      return wires_.value()[i];
    }

    // calculated property accessors
    [[nodiscard]] Wires mirror() const {
      auto b = border_.mirror();
      proto_t w;
      w.reserve(wires_.has_value() ? wires_.value().size() : 0u);
      if (wires_.has_value()) for (const auto &w_: wires_.value()) w.push_back(w_.mirror());
      return wires_.has_value() ? Wires(b, w) : Wires(b);
    }

    [[nodiscard]] wire_t::base_t indexes() const {
      wire_t::base_t all;
      auto not_in_all = [&all](const auto &i) { return std::find(all.begin(), all.end(), i) == all.end(); };
      for (const auto &x: border_) if (not_in_all(x)) all.push_back(x);
      if (wires_.has_value())
        for (const auto &wire: wires_.value())
          for (const auto &x: wire)
            if (not_in_all(x))
              all.push_back(x);
      return all;
    }

    [[nodiscard]] polystar::Array2<ind_t> border_edges() const {
      auto n = border_.size();
      polystar::Array2<ind_t> be(static_cast<ind_t>(n + 1), 2u);
      for (ind_t i = 0; i < n; ++i) {
        auto edge = border_.edge(i);
        be.val(i, 0) = edge.first;
        be.val(i, 1) = edge.second;
      }
      return be;
    }

    template<class T, template<class> class A, class U = typename A<T>::shape_t>
    [[nodiscard]] std::enable_if_t <polystar::isArray<T, A>, std::tuple<A<T>, A<T>>>
    borders(const A<T> &x) const {
      auto n = static_cast<ind_t>(border_.size());
      auto a = 0 * x.view(0);
      a.resize(n);
      auto b = 0 * x.view(0);
      b.resize(n);
      for (size_t i = 0; i < border_.size(); ++i) {
        auto e = border_.edge(i);
        a.set(i, x.view(e.first));
        b.set(i, x.view(e.second));
      }
      return std::make_tuple(a, b);
    }

    template<class T, template<class> class A, class U = typename A<T>::shape_t>
    [[nodiscard]] T area(const A<T> &x) const {
      auto result = border_.area(x);
      // result should be positive, and any wires should be holes which reduce the area
      if (wires_.has_value()) for (const auto &w: wires_.value()) result += w.area(x);
      return result;
    }

    template<class T, template<class> class A, class U = typename A<T>::shape_t>
    [[nodiscard]] A<T> centroid(const A<T> &x) const {
      auto result = border_.centroid(x);
      if (wires_.has_value()) for (const auto &w: wires_.value()) result += w.centroid(x);
      return result;
    }

    template<class T, template<class> class A, class U = typename A<T>::shape_t>
    [[nodiscard]] T circumscribed_radius(const A<T> &x) const {
      // only the *border* contributes to the circumscribed radius; holes inside would only shift the centroid
      return border_.circumscribed_radius(x);
    }

    Wires combine(const Wires &that, const ind_t offset) const {
      auto b = border_.combine(that.border(), offset);
      auto w = wires_.has_value() ? wires_.value() : proto_t();
      w.reserve(w.size() + that.wire_count());
      if (that.wire_count()) for (const auto &tw: that.wires()) w.push_back(tw.offset(offset));
      return w.empty() ? Wires(b) : Wires(b, w);
    }

    template<class T, class R, template<class> class A,
      template<class> class B>
    std::tuple <A<T>, Wires>
    cut(const A<T> &x, const B<R> &a, const B<R> &b, const R tol = R(0), const int dig = 1) const {
      assert(a.ndim() == 2 && b.ndim() == 2);
      assert(a.size(1) == 2 && b.size(1) == 2);
      assert(a.size(0) == b.size(0));
      auto ov = (T(1) * x).decouple();
      Wires ow(*this);
      for (ind_t i = 0; i < a.size(0); ++i) std::tie(ov, ow) = ow.one_cut(ov, a.view(i), b.view(i), tol, dig);
      return std::make_tuple(ov, ow);
    }

    template<class T, class R, template<class> class A,
      template<class> class B>
    std::tuple <A<T>, Wires>
    one_cut(const A<T> &x, const B<R> &a, const B<R> &b, const R tol = R(0), const int dig = 1) const {
      assert(a.ndim() == 2 && b.ndim() == 2);
      assert(a.size(1) == 2 && b.size(1) == 2);
      assert(a.size(0) == b.size(0) && a.size(0) == 1);
      auto ov = (T(1) * x).decouple();
      if (none_beyond(x, a, b)) return std::make_tuple(x, Wires());
      // otherwise find the vertices beyond a->b and remove them

      return std::make_tuple(x, Wires(*this));
    }

    template<class T, class R, template<class> class A,
      template<class> class B>
    bool none_beyond(const A<T> &x, const B<R> &a, const B<R> &b) const {
      if (border_.any_beyond(x, a, b)) return false;
      if (wires_.has_value()) for (const auto &w: wires_.value()) if (w.any_beyond(x, a, b)) return false;
      return true;
    }

    std::string python_string() const {
      std::stringstream s;
      s << border_.python_string();
      s << ",[";
      size_t count{0};
      if (wires_.has_value())
        for (const auto &w: wires_.value())
          s << w.python_string() << (++count < wires_.value().size() ? "," : "");
      s << "]";
      return s.str();
    }

    friend std::ostream &operator<<(std::ostream &os, Wires w) {
      os << w.python_string();
      return os;
    }

    template<class T, template<class> class A>
    bool contains(const A<T> &point, const A<T> &x) const {
      auto segment = cat(0, point, (std::numeric_limits<T>::max)() + T(0) * point);
      wire_t::edge_t se{0, 1};
      size_t crossings{0};
      for (size_t i = 0; i < border_.size(); ++i) if (intersect2d(segment, se, x, border_.edge(i))) ++crossings;
      if (wires_.has_value())
        for (const auto &w: wires_.value())
          for (size_t i = 0; i < w.size(); ++i) if (intersect2d(segment, se, x, w.edge(i))) ++crossings;
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

} // namespace polystar::polygon
#endif