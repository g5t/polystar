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
#include "polygon_clipping.hpp"

#include "svg.hpp"
#include "triangle.hpp"

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
    [[nodiscard]] Poly<T,A> simplify() const {return Poly(vertices_, wires_.simplify(vertices_));}
//    Network<Wire,T,A> triangulate() const {
//      return wires_.triangulate(vertices_);
//    }
    Network<Wire,T,A> triangulate() const { return polystar::triangle::triangulate<Wire>(vertices_, wires_); }
    [[nodiscard]] bool is_not_approx(const Poly<T,A> & that, const T tol=T(0), const int dig=1) const {
      bool permuted{false};
      if (vertices_ != that.vertices()){
        permuted = vertices_.is_permutation(that.vertices(), tol, tol, dig);
        if (!permuted) return true;
      }
      if (permuted){
        auto permutation = vertices_.permutation_vector(that.vertices(), tol, tol, dig);
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

    Poly<T,A> operator+(const Poly<T,A>& that) const {return intersection(that);}
    Poly<T,A> combine(const Poly<T,A>& that, const T tol=T(0), const int dig=1) const {
      // combine vertices
      auto v = cat(0, vertices_, that.vertices());
      // combined wires
      Wires w = wires_.combine(that.wires(), vertices_.size(0));
      // check for duplicate vertices
      std::tie(v, w) = remove_duplicate_points_and_update_wire_indexing(v, w, tol, dig);
      // look for overlapping wires now that vertex indexing has been modified
      return Poly<T,A>(v, remove_extraneous_wires(w));
    }
    Poly<T,A> combine_all(const std::vector<Poly<T,A>> & others, const T tol=T(0), const int dig=1) const {
      if (others.empty()) return *this;
      auto out = intersection(others.front(), tol, dig);
      for (auto ptr = others.begin()+1; ptr != others.end(); ++ptr) out = out.intersection(*ptr, tol, dig);
      return out;
    }

    Poly<T,A> mirror() const {return {T(-1) * vertices_, wires_.mirror()};}
    Poly<T,A> inverse() const {return {vertices_, wires_.inverse()};}
    Poly<T,A> centre() const {return {vertices_ - centroid(), wires_};}
    Poly<T,A> translate(const A<T>& v) const {return {vertices_ + v, wires_};}

    Poly<T,A> transform(const std::array<T,4>& matrix) const {return {matrix * vertices_, wires_};}
    Poly<T,A> skew(T factor, int source, int sink) const {
      if (source == sink) {
        throw std::logic_error("polystar::polygon::Poly::skew: source can not equal sink");
      }
      if (source < 0 || source > 1 || sink < 0 || sink > 1) {
        throw std::logic_error("polystar::polygon::Poly::skew: source and sink must be in bounds");
      }
      std::array<T, 4> matrix{{1, 0, 0, 1}};
      matrix[source * 2 + sink] = factor;
      return transform(matrix);
    }

    template<class R, template<class> class B>
      [[nodiscard]] std::vector<bool> border_contains(const B<R>& x) const {
      return wires_.border_contains(x, vertices_);
    }
    template<class R, template<class> class B>
      [[nodiscard]] std::vector<bool> contains(const B<R>& x) const {
      return wires_.contains(x, vertices_);
    }
    template<class R>
      [[nodiscard]] std::vector<bool> contains(const std::vector<std::array<R,2>> & x) const {
      return contains(from_std_like(vertices_, x));
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

    [[nodiscard]] Poly<T,A> intersection(const Poly<T,A>& that, T tol = T(0), int dig = 0) const {
      // simple-case checks first: ignore the possibility of negative polygon inclusion?
      auto that_in_this = border_contains(that.vertices());
      if (std::all_of(that_in_this.begin(), that_in_this.end(), [](const auto x){return x;})) {
        if (that.area() > T(0) && that.wires().wire_count() == 0) return *this;
        return insert_hole(that.simplify(), tol, dig);
      }
      auto this_in_that = that.border_contains(vertices_);
      if (std::all_of(this_in_that.begin(), this_in_that.end(), [](const auto x){return x;})) {
        if (area() > T(0) && wires_.wire_count() == 0) return that;
        return that.insert_hole(simplify(), tol, dig);
      }
      if (wires_.wire_count() || that.wires().wire_count()){
        return intersection(triangulate(), that.triangulate());
      }
      return intersection(vertices_, wires_.border(), that.vertices(), that.wires().border());
    }
//    template<class R, template<class> class B>
//    [[nodiscard]] std::enable_if_t<isArray<R,B>, Poly<T,A>> cut(const B<R>& a, const B<R>& b, const R tol=R(0), const int dig=1) const {
//      auto [v, w] = wires_.cut(vertices_, a, b, tol, dig);
//      return {v, w};
//    }
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

    void add_to_svg(SVG::Path & path) const {
      wires_.add_to_svg(path, vertices_);
    }
    SVG::SVG to_svg(const std::optional<std::string_view> fill, const std::optional<std::string_view> stroke) const {
      SVG::SVG svg;
      svg.style("path").set_attr("fill", std::string(fill.value_or("tan")));
      svg.style("path").set_attr("stroke",std::string(stroke.value_or("black")));
      auto path = svg.add_child<SVG::Path>();
      add_to_svg(*path);
      svg.autoscale();
      return svg;
    }
    SVG::SVG to_svg(const std::string_view fill) const {
      return to_svg(std::make_optional(fill), std::nullopt);
    }
//    SVG::SVG to_svg(const std::string_view fill, const std::string_view stroke) const {
//      return to_svg(std::make_optional(fill), std::make_optional(stroke));
//    }
    SVG::SVG to_svg() const {
      return to_svg(std::nullopt, std::nullopt);
    }
    template<class... Args>
    void write_svg(const std::string & filename, Args... args) const {
      auto svg = to_svg(args...);
      auto of = std::ofstream(filename);
      of << std::string(svg);
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

    Poly<T,A> insert_hole(const Poly<T,A> & hole, const T tol, const int dig) const {
      if (hole.area() > 0) return *this;
      auto ov = vertices_.decouple();
      // combine the vertices from both polygons into a new vertex list,
      auto hv = hole.vertices();
      auto hb = hole.wires().border();
      typename wires_t::wire_t hw;
      hw.reserve(hb.size());
      for (const auto b: hb){
        // find the index of hv.view(b) in ov
        auto oi = ov.first(polystar::cmp::eq, hv.view(b), tol, tol, dig); // does this work like I think it does?
        // if hv.view(b) is not in ov, oi == oi.size(0) -- append it in that case
        if (oi >= ov.size(0)) ov = cat(0, ov, hv.view(b));
        // push the vertex to hw
        hw.push_back(oi);
      }
      return {ov, wires_.insert_wire(hw)};

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

  namespace utils {
    template<class T, template<class> class A>
    bool
    could_not_overlap(const A<T> &va, const Wire &wa, const A<T> &vb, const Wire &wb) {
      auto a_radius = wa.circumscribed_radius(va);
      auto b_radius = wb.circumscribed_radius(vb);
      auto a_centre = wa.centroid(va);
      auto b_centre = wb.centroid(vb);
      auto d = a_centre.distance(b_centre);
      return d > a_radius + b_radius;
    }


    template<class T, template<class> class A>
    bool
    contains(const A<T> &v0, const Wire &w0, const A<T> &v1, const Wire &w1) {
      auto one_in_zero = w0.contains(v1.extract(w1), v0);
      // if all are contained, then the whole is contained
      return std::all_of(one_in_zero.begin(), one_in_zero.end(), [](const auto x) { return x; });
    }
  }


  /// \brief Find all intersection points for two polygon wires, then return their intersection polygons
  template<class T, template<class> class A>
  std::vector<Poly<T, A>>
  intersection(const A<T> & va, const Wire & wa, const A<T> & vb, const Wire & wb){
    // If the two polygons are too far apart, return nothing
    if (utils::could_not_overlap(va, wa, vb, wb)) return {};
    // If all of A is inside B, return A
    if (utils::contains(vb, wb, va, wa)) return {Poly<T, A>(va, wa)};
    // If all of B is inside A, return b
    if (utils::contains(va, wa, vb, wb)) return {Poly<T, A>(vb, wb)};
    // Otherwise, find the intersection points:
    // Combine the two sets of vertices into a single list
    auto v = cat(0, va.extract(wa), vb.extract(wb));
    // Make the new wires for A and B, since we may have reduced or re-ordered their vertices:
    Wire new_wa, new_wb;
    new_wa.resize(wa.size());
    new_wb.resize(wb.size());
    std::iota(new_wa.begin(), new_wa.end(), 0);
    std::iota(new_wb.begin(), new_wb.end(), wa.size());
    // Construct the dual doubly-linked lists of vertex indices
    auto lists = clip::VertexLists(new_wa, new_wb);
    auto no_problems = clip::weiler_atherton(v, lists);
    if (no_problems){
      std::cerr << "Warning: " << no_problems << " problems found in the Weiler-Atherton algorithm" << std::endl;
    }
    auto wires = lists.intersection_wires();
    std::vector<Poly<T, A>> result;
    result.reserve(wires.size());
    for (const auto & w: wires) {
      auto p = Poly<T, A>(v, w);
      if (p.area()) result.push_back(p);
    }
    return result;
  }

  /// \brief Find the intersection of two polygons, triangulating if either has any holes
  template<class T, template<class> class A>
  std::vector<Poly<T, A>>
  intersection(const Poly<T, A> & a, const Poly<T, A> & b) {
    if (a.wires().wire_count() || b.wires().wire_count()){
      return intersection(a.triangulate(), b.triangulate());
    }
    return intersection(a.vertices(), a.wires().border(), b.vertices(), b.wires().border());
  }

  /// \brief Find the intersection of two triangulated polygon networks by considering all pairs of triangles
  template<class T, template<class> class A, class W=Wire>
  std::vector<Poly<T, A>>
  intersection(const Network<W, T, A> & a, const Network<W, T, A> & b) {
    std::vector<Poly<T, A>> result;
    result.reserve(std::max(a.size(), b.size()));
    for (const auto & wa: a.wires()) {
      for (const auto & wb: b.wires()) {
        auto ps = intersection(a.vertices(), wa, b.vertices(), wb);
        for (const auto & p: ps) if (p.area()) result.push_back(p);
      }
    }
    // TODO go through the resulting polygons and combine any that are touching?
    return result;
  }


  /// \brief Find all intersection points for two polygon wires, then return their union polygons
  template<class T, template<class> class A>
  std::vector<Poly<T, A>>
  union_polygons(const A<T> & va, const Wire & wa, const A<T> & vb, const Wire & wb){
    // If the two polygons are too far apart, return both
    if (utils::could_not_overlap(va, wa, vb, wb)) return {{va, wa}, {vb, wb}};
    // If all of A is inside B, return B
    auto a_in_b = wb.contains(va.extract(wa), vb);
    if (utils::contains(vb, wb, va, wa)) return {{vb, wb}};
    // If all of B is inside A, return A
    if (utils::contains(va, wa, vb, wb)) return {{va, wa}};
    // Otherwise, find the intersection points:
    // Combine the two sets of vertices into a single list
    auto v = cat(0, va.extract(wa), vb.extract(wb));
    // Make the new wires for A and B, since we may have reduced or re-ordered their vertices:
    Wire new_wa, new_wb;
    new_wa.resize(wa.size());
    new_wb.resize(wb.size());
    std::iota(new_wa.begin(), new_wa.end(), 0);
    std::iota(new_wb.begin(), new_wb.end(), wa.size());
    // Construct the dual doubly-linked lists of vertex indices
    auto lists = clip::VertexLists(new_wa, new_wb);
    auto no_problems = clip::weiler_atherton(v, lists);
    if (no_problems){
      std::cerr << "Warning: " << no_problems << " problems found in the Weiler-Atherton algorithm" << std::endl;
    }
    auto wires = lists.union_wires();
    // the union wires returned vector is {border, hole0, hole1, ...}
    std::vector<Poly<T, A>> result;
    result.reserve(1);
    auto border = wires.front();
    auto holes = std::vector<Wire>(wires.begin() + 1, wires.end());
    result.push_back(Poly<T, A>(v, border, holes));
    return result;
  }

  /// \brief Find the intersection of two polygons, triangulating if either has any holes
  template<class T, template<class> class A>
  std::vector<Poly<T, A>>
  union_polygons(const Poly<T, A> & a, const Poly<T, A> & b) {
    if (a.wires().wire_count() || b.wires().wire_count()){
      return union_polygons(a.triangulate(), b.triangulate());
    }
    return union_polygons(a.vertices(), a.wires().border(), b.vertices(), b.wires().border());
  }

  /// \brief Find the intersection of two triangulated polygon networks by considering all pairs of triangles
  template<class T, template<class> class A, class W=Wire>
  std::vector<Poly<T, A>>
  union_polygons(const Network<W, T, A> & a, const Network<W, T, A> & b) {
    std::vector<Poly<T, A>> result;
    result.reserve(std::max(a.size(), b.size()));
    for (const auto & wa: a.wires()) {
      for (const auto & wb: b.wires()) {
        auto ps = union_polygons(a.vertices(), wa, b.vertices(), wb);
        for (const auto & p: ps) if (p.area()) result.push_back(p);
      }
    }
    // TODO go through the resulting polygons and combine any that are touching?
    return result;
  }


}

#endif