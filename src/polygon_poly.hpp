#ifndef POLYSTAR_POLYGON_POLY_HPP
#define POLYSTAR_POLYGON_POLY_HPP
#include <vector>
#include <array>
#include <optional>
#include <numeric>
#include <limits>

#include "comparisons.hpp"
#include "array_.hpp"
#include "geometry.hpp"
#include "approx_float.hpp"
#include "polygon_wires.hpp"
#include "polygon_network.hpp"

namespace polystar::polygon{
  using ind_t = polystar::ind_t;

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
    [[nodiscard]] Poly<T,A> simplify() const {return Poly(vertices_, wires_.simplify(vertices));}
    Network<T,A> triangulate() const {
      return wires_.triangulate(vertices_);
    }
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