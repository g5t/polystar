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
#include "svg.hpp"

namespace polystar::polygon {
  using ind_t = polystar::ind_t;
  class Wire;

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

    Wires(Wires &&that) noexcept: border_(std::move(that.border_)), wires_(std::move(that.wires_)) {}

    Wires &operator=(const Wires &that) = default;

    bool operator!=(const Wires & that) const {
      if (border_ != that.border()) return true;
      auto ow = that.wires();
      if (wires_.has_value()) {
        if (ow.size() != wires_->size()) return true;
        /* wires_.value().size() is unavailable for some macOS variants
         * We've already checked that wires_ has a value, so we can safely
         * use pointer dereferencing to access that value.
         * https://stackoverflow.com/questions/44217316/how-do-i-use-stdoptional-in-c
         * */
        size_t count{0};
        for (const auto &w: *wires_){
          if (std::find_if(ow.begin(), ow.end(), [&](const auto &x) { return x == w; }) != ow.end()) count++;
        }
        return ow.size() != count;
      }
      // if that.wires is empty then they're the same
      return !ow.empty();
    }

    bool operator==(const Wires &that) const { return !this->operator!=(that); }

    template<class I>
    Wires permute(const std::vector <I> &permutation) const {
      auto border = border_.replace(permutation);
      if (wires_.has_value()) {
        proto_t wires;
        wires.reserve(wires_->size());
        for (const auto &wire: *wires_) wires.push_back(wire.replace(permutation));
        return Wires(border, wires);
      }
      return Wires(border);
    }

    // Remove any wires_ by cleverly appending them to border_
    template<class T, template<class> class A>
      Wires simplify(const A<T>& vertices) const {
      // if there are no internal wires, there's nothing to do:
      if (!wires_.has_value()) return Wires(border_);
      // HUGE ASSUMPTION: All wires must be fully inside the border, arranged to have 'negative' area (they're holes)
      //                  and do not overlap with one another.
      /* 1. Pick an internal wire which has not had an edge failure
            If there are no such wires, the polygon can not be simplified.
         2. Find an edge between one of its vertices and a border vertex which is inside the polygon
            AND WHICH IS NOT ALREADY IN THE WIRE TWICE!
            (The edge must be *inside* the border and *outside* all wires!)
            If there are no such edges, record the failure and go back to 1.
         3. Insert the wire in the border (the chosen border vertex and wire vertex are visited twice each)
         4. If there are any remaining wires, clear the failed-edge list and go back to 1.
      */
      auto b = border_;
      auto ws = *wires_;
//      size_t failures{0u};
      decltype(ws)::difference_type failures{0u};

      using p_t = typename wire_t::edge_t;
      using pp_t = std::pair<p_t, p_t>;

      // Count how many times each vertex appears in the border -- update this going along
      std::vector<size_t> border_counts(vertices.size(0), 0u); // the border and wires index vertices_, so this should be safe
      for (const auto & x: b) border_counts[x]++;
      size_t max_counts{2};
      // allow for the maximum connections at one point to be increased if we otherwise fail
      while (static_cast<size_t>(failures) < ws.size() && max_counts < 10u){
        std::cout << "border counts: ";
        for (const auto & bc: border_counts) std::cout << bc << " ";
        std::cout << "\n";
        auto w = ws[failures];
        // find the possible joining edges
        std::vector<pp_t> choices;
        choices.reserve(w.size() * b.size());
        for (ind_t i=0; i < w.size(); ++i) if (border_counts[w[i]] < max_counts) {
          for (ind_t j=0; j < b.size(); ++j) if (border_counts[b[j]] < max_counts) {
            p_t se{w[i], b[j]};
            if (!intersects(vertices, se, vertices, end_type::neither)) choices.emplace_back(std::make_pair(i, j), se);
          }
        }
        if (choices.empty()) {
          ++failures;
        } else {
          // pick *which* possible edge we want, maybe the shortest one is best?
          std::vector<std::pair<double, pp_t>> lengths;
          lengths.reserve(choices.size());
          std::transform(choices.begin(), choices.end(), std::back_inserter(lengths), [&](const auto & ch){
            return std::make_pair(norm(vertices.view(ch.second.first) - vertices.view(ch.second.second)).sum(), ch);
          });
          std::cout << "possible edges:\n";
          for (const auto & [l, e]: lengths)
            std::cout << "[" << e.first.first << "][" << e.first.second << "]("
                      << e.second.first << "--" << e.second.second << ") " << l
                      << " with " << border_counts[e.second.first]  << " and " << border_counts[e.second.second]
                      << " pre existing connections \n";

          auto ptr = std::min_element(lengths.begin(), lengths.end(),
                                      [](const auto & a, const auto & b){return a.first < b.first;});
          auto best = ptr->second.first; // (i, j) from above, index into w and b respectively
          std::cout << "best edge indexes: (" << best.first << "--" << best.second << ")\n";
          // Add the edge border[bi] -- wire[failures][wi]
          // The new border is b[:bi] -- w[wi:] -- w[:wi] -- b[bi:]
          // Where x[:n] means [0,n] {upper bound inclusive} and x[n:] means [n,x.size()) {upper bound exclusive}
          typename wire_t::base_t new_b;
          new_b.reserve(b.size() + w.size());
          for (size_t i=0; i<std::min(static_cast<size_t>(best.second+1), b.size()); ++i) new_b.push_back(b[i]);
          for (size_t i=best.first; i<w.size(); ++i) new_b.push_back(w[i]);
          for (size_t i=0; i<std::min(static_cast<size_t>(best.first+1), w.size()); ++i) new_b.push_back(w[i]);
          for (size_t i=best.second; i<b.size(); ++i) new_b.push_back(b[i]);
          // record that the wires are now part of the border
          for (const auto & wi: w) border_counts[wi]+=1;
          // and that both the connected wire point and connected border point now have an extra connection:
          border_counts[w[best.first]]+=1;
          border_counts[b[best.second]]+=1;

          // the border is now new_b
          b = Wire(new_b);
          // remove ws[failures]
          ws.erase(ws.begin()+failures);
          // reset the failures counter
          failures = 0u;
        }
        if (static_cast<size_t>(failures) >= ws.size()) { // failures is always positive, so hopefully this is ok
          // relax the repeated border point restriction (this makes triangulation harder)
          max_counts++;
          // start again from the front of the wires list
          failures = 0;
        }
      }
      if (!ws.empty()) {
        std::cout << "Simplification of the complex polygon is not possible\n";
        return Wires(b, ws);
      }
      return Wires(b);
    }

    template<class T, template<class> class A>
      Network<Wire,T,A> triangulate(const A<T>& v) const {
      auto s = simplify(v);
      auto tri = s.border().triangulate(v);
      tri.off_limits(border_);
      return tri;
    }


    template<class T, template<class> class A>
      bool intersects(const A<T>& other, const typename wire_t::edge_t edge, const A<T> & ours, end_type inclusive = end_type::both) const {
      if (border_.intersects(other, edge, ours, inclusive)) return true;
      if (wires_.has_value()) for (const auto & w: *wires_) if (w.intersects(other, edge, ours, inclusive)) return true;
      // no intersection ... but the edge could be *inside* the border (... and not inside a hole) -- push this to another function
      return false;
    }

    // direct property accessors
    [[nodiscard]] size_t border_size() const { return border_.size(); }

    [[nodiscard]] wire_t border() const { return border_; }

    [[nodiscard]] size_t wire_count() const { return wires_.has_value() ? wires_->size() : 0; }

    [[nodiscard]] proto_t wires() const { return wires_.has_value() ? *wires_ : proto_t(); }

    [[nodiscard]] wire_t wire(const ind_t &i) const {
      assert(wires_.has_value());
      assert(i < wire_count());
      return (*wires_)[i];
    }

    // calculated property accessors
    [[nodiscard]] Wires mirror() const {
      auto b = border_.mirror();
      proto_t w;
      w.reserve(wires_.has_value() ? wires_->size() : 0u);
      if (wires_.has_value()) for (const auto &w_: *wires_) w.push_back(w_.mirror());
      return wires_.has_value() ? Wires(b, w) : Wires(b);
    }

    [[nodiscard]] Wires inverse() const{
      return mirror();
    }

    [[nodiscard]] wire_t::base_t indexes() const {
      wire_t::base_t all;
      auto not_in_all = [&all](const auto &i) { return std::find(all.begin(), all.end(), i) == all.end(); };
      for (const auto &x: border_) if (not_in_all(x)) all.push_back(x);
      if (wires_.has_value())
        for (const auto &wire: *wires_)
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
    [[nodiscard]] bool indexing_error(const A<T> &x) const {
      if (border_.indexing_error(x)) return true;
      if (wires_.has_value()) for (const auto & w: *wires_) if (w.indexing_error(x)) return true;
      return false;
    }

    template<class T, template<class> class A, class U = typename A<T>::shape_t>
    [[nodiscard]] T area(const A<T> &x) const {
      auto result = border_.area(x);
      // result should be positive, and any wires should be holes which reduce the area
      if (wires_.has_value()) for (const auto &w: *wires_) result += w.area(x);
      return result;
    }

    template<class T, template<class> class A, class U = typename A<T>::shape_t>
    [[nodiscard]] A<T> centroid(const A<T> &x) const {
      auto result = border_.centroid(x);
      if (wires_.has_value()) for (const auto &w: *wires_) result += w.centroid(x);
      return result;
    }

    template<class T, template<class> class A, class U = typename A<T>::shape_t>
    [[nodiscard]] T circumscribed_radius(const A<T> &x) const {
      // only the *border* contributes to the circumscribed radius; holes inside would only shift the centroid
      return border_.circumscribed_radius(x);
    }

    [[nodiscard]] Wires combine(const Wires &that, const ind_t offset) const {
      auto b = border_.combine(that.border(), offset);
      auto w = wires_.has_value() ? *wires_ : proto_t();
      w.reserve(w.size() + that.wire_count());
      if (that.wire_count()) for (const auto &tw: that.wires()) w.push_back(tw.offset(offset));
      return w.empty() ? Wires(b) : Wires(b, w);
    }

    [[nodiscard]] Wires insert_wire(const wire_t & w) const {
      auto ws = wires_.has_value() ? *wires_ : proto_t();
      ws.push_back(w);
      return Wires(border_, ws);
    }

//    template<class T, class R, template<class> class A,
//      template<class> class B>
//    std::tuple <A<T>, Wires>
//    cut(const A<T> &x, const B<R> &a, const B<R> &b, const R tol = R(0), const int dig = 1) const {
//      assert(a.ndim() == 2 && b.ndim() == 2);
//      assert(a.size(1) == 2 && b.size(1) == 2);
//      assert(a.size(0) == b.size(0));
//      auto ov = (T(1) * x).decouple();
//      Wires ow(*this);
//      for (ind_t i = 0; i < a.size(0); ++i) std::tie(ov, ow) = ow.one_cut(ov, a.view(i), b.view(i), tol, dig);
//      return std::make_tuple(ov, ow);
//    }
//
//    template<class T, class R, template<class> class A,
//      template<class> class B>
//    std::tuple <A<T>, Wires>
//    one_cut(const A<T> &x, const B<R> &a, const B<R> &b, const R tol = R(0), const int dig = 1) const {
//      assert(a.ndim() == 2 && b.ndim() == 2);
//      assert(a.size(1) == 2 && b.size(1) == 2);
//      assert(a.size(0) == b.size(0) && a.size(0) == 1);
//      auto ov = (T(1) * x).decouple();
//      if (none_beyond(x, a, b)) return std::make_tuple(x, Wires());
//      // otherwise find the vertices beyond a->b and remove them
//
//      return std::make_tuple(x, Wires(*this));
//    }

    template<class T, class R, template<class> class A,
      template<class> class B>
    bool none_beyond(const A<T> &x, const B<R> &a, const B<R> &b) const {
      if (border_.any_beyond(x, a, b)) return false;
      if (wires_.has_value()) for (const auto &w: *wires_) if (w.any_beyond(x, a, b)) return false;
      return true;
    }

    [[nodiscard]] std::string python_string() const {
      std::stringstream s;
      s << border_.python_string();
      s << ",[";
      size_t count{0};
      if (wires_.has_value())
        for (const auto &w: *wires_)
          s << w.python_string() << (++count < wires_->size() ? "," : "");
      s << "]";
      return s.str();
    }

    friend std::ostream &operator<<(std::ostream &os, const Wires& w) {
      os << w.python_string();
      return os;
    }

    template<class T, template<class> class A, class R, template<class> class B>
    std::vector<bool> border_contains(const B<R> &point, const A<T> &x, end_type inclusive = end_type::second) const {
      std::vector<bool> out;
      out.reserve(static_cast<size_t>(point.size(0)));
      wire_t::edge_t se{0, 1};
      for (ind_t p=0; p<point.size(0); ++p) {
        auto segment = cat(0, point.view(p), T(2) * x.max(0));
        size_t crossings{0};
        for (size_t i = 0; i < border_.size(); ++i)
          if (intersect2d(segment, se, x, border_.edge(i), inclusive)) ++crossings;
        out.push_back((crossings % 2u) == 1u);
      }
      return out;
    }

    template<class T, template<class> class A, class R, template<class> class B>
    std::vector<bool> contains(const B<R> &point, const A<T> &x, end_type inclusive = end_type::second) const {
      std::vector<bool> out;
      out.reserve(static_cast<size_t>(point.size(0)));
      wire_t::edge_t se{0, 1};
      for (ind_t p=0; p<point.size(0); ++p) {
        auto segment = cat(0, point.view(p), T(2) * x.max(0));
        size_t crossings{0};
        for (size_t i = 0; i < border_.size(); ++i)
          if (intersect2d(segment, se, x, border_.edge(i), inclusive))
            ++crossings;
        if (wires_.has_value())
          for (const auto &w: *wires_)
            for (size_t i = 0; i < w.size(); ++i) if (intersect2d(segment, se, x, w.edge(i), inclusive)) ++crossings;
        out.push_back((crossings % 2u) == 1u);
      }
      return out;
    }

    template<class T, template<class> class A>
    bool overlaps(const A<T>& other, const typename wire_t::edge_t edge, const A<T>& ours, end_type inclusive = end_type::both) const {
      // whether the edge crosses through this polygon *or* either end is *inside* the polygon without being in a hole
      if (intersects(other, edge, ours)) return true;
      const auto p0{other.view(edge.first)};
      const auto p1{other.view(edge.second)};
      if (!border_.contains(p0, ours, inclusive) && !border_.contains(p1, ours, inclusive)) return false;
      if (wires_.has_value()) for (const auto & w: *wires_)
        if (w.contains(p0, ours, inclusive) && w.contains(p1, ours, inclusive)) return false;
      return true;
    }

//    template<class P, class T, template<class> class A>
//    P intersection(const P& that, const A<T> & ours, const R tol=R(0), const int dig=1) const {
//
//    }

    template<class T, template<class> class A>
    void add_to_svg(SVG::Path & path, const A<T>& x) const {
      border_.add_to_svg(path, x);
      if (wires_.has_value()) for (const auto & w: *wires_) w.add_to_svg(path, x);
    }


#ifdef USE_HIGHFIVE
    template<class H> std::enable_if_t<std::is_base_of_v<HighFive::Object, H>, bool>
      to_hdf(H& obj, const std::string & entry) const {
      auto group = overwrite_group(obj, entry);
      bool ok{true};
      ok &= list_to_hdf<ind_t>(border_, group, "border");
      std::vector<std::vector<ind_t>> copied_wires; // not great, but vector<Wire> can´t be figured-out by the compiler
      if (wires_.has_value()) for (const auto & w: *wires_) copied_wires.push_back(w);
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
  std::enable_if_t<isArray<T,A>, std::tuple<A<T>, Wires>>
  remove_duplicate_points_and_update_wire_indexing(const A<T>& points, const Wires & ws, const T tol=T(0), const int dig=1){
    auto are_unique = points.is_unique(tol, dig);
    if(std::find(are_unique.begin(), are_unique.end(), false) != are_unique.end()){
      // find unique *existing* indexes
      auto index = points.unique_idx(tol, dig);
      // and update them to point into reduced, only unique, vertices
      ind_t cnt{0};
      for (ind_t j=0; j < points.size(0); ++j){
        index[j] = are_unique[j] ? cnt++ : index[index[j]];
      }
      auto b = ws.border();
      auto hs = ws.wires();
      for (auto & x: b) x = index[x];
      for (auto & h: hs) for (auto & x: h) x= index[x];
      return std::make_tuple(points.extract(are_unique), Wires(b, hs));
    }
    return std::make_tuple(points, ws);
  }

  Wires remove_extraneous_wires(const Wires & ws);


} // namespace polystar::polygon
#endif