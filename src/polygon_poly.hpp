#ifndef POLYSTAR_POLYGON_POLY_HPP
#define POLYSTAR_POLYGON_POLY_HPP

#ifndef POLYSTAR_POLYGON_HPP
#pragma message("polygon_poly.hpp included without first including polygon.hpp"\
                " Network defines missing"\
)
#endif

#include <utility>
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
    Poly(const vertex_t & v, wires_t::wire_t w): vertices_(v), wires_(std::move(w)) {
      check_wire_indexing();
    }
    Poly(vertex_t && v, wires_t::wire_t && w): vertices_(std::move(v)), wires_(std::move(w)) {
      check_wire_indexing();
    }
    // Complex polygon constructors
    Poly(const vertex_t & v, wires_t::wire_t b, const wires_t::proto_t & w): vertices_(v), wires_(std::move(b), w) {
      check_wire_indexing();
    }
    Poly(vertex_t && v, wires_t::wire_t && b, wires_t::proto_t && w): vertices_(std::move(v)), wires_(std::move(b), std::move(w)) {
      check_wire_indexing();
    }
    Poly(const vertex_t & v, const wires_t & w): vertices_(v), wires_(w) {
      check_wire_indexing();
    }
    Poly(vertex_t && v, wires_t && w): vertices_(std::move(v)), wires_(std::move(w)) {
      check_wire_indexing();
    }
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
    [[nodiscard]] bool is_convex() const {
      debug_update_if(wires_.wire_count(), "Non-convex due to internal wires");
      if (wires_.wire_count()) return false;
      return wires_.border().is_convex(vertices_);
    }
    [[nodiscard]] Poly<T,A> convex_hull() const {return Poly(vertices_);}
    /// \brief Remove any internal wires by adding them to the border via new zero-width gaps.
    [[nodiscard]] Poly<T,A> simplify() const {return Poly(vertices_, wires_.simplify(vertices_));}
    /// \brief Ensure only indexed vertices are kept in the returned Poly object
    [[nodiscard]] Poly<T,A> without_extraneous_vertices() const {
      // find the unique vertex indices in the border and all holes
      auto b = wires_.border();
      for (const auto & w: wires_.wires()){
        for (const auto & i: w) {
          if (std::find(b.begin(), b.end(), i) == b.end()) {
            b.push_back(i);
          }
        }
      }
      // extract only those vertices
      auto v = vertices_.extract(b);
      // find the map from old vertex indexing to new
      std::map<ind_t, ind_t> map;
      for (const auto & x: b){
        auto m = v.row_is(cmp::eq, vertices_.view(x)).first();
        map[x] = m;
      }
      // use the map to translate
      auto translate = [&map](const Wire & from, Wire & to){
        to.reserve(from.size());
        for (const auto & i: from) to.push_back(map[i]);
      };
      // create the new border and holes
      Wire nb;
      translate(wires_.border(), nb);
      std::vector<Wire> nh;
      nh.reserve(wires_.wire_count());
      for (const auto & w: wires_.wires()){
        Wire nw;
        translate(w, nw);
        nh.push_back(nw);
      }
      return Poly(v, nb, nh);
    }
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

    //Poly<T,A> operator+(const Poly<T,A>& that) const {return intersection(that);}
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
//    Poly<T,A> combine_all(const std::vector<Poly<T,A>> & others, const T tol=T(0), const int dig=1) const {
//      if (others.empty()) return *this;
//      auto out = intersection(others.front(), tol, dig);
//      for (auto ptr = others.begin()+1; ptr != others.end(); ++ptr) out = out.intersection(*ptr, tol, dig);
//      return out;
//    }

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
      auto a = area();
      if (approx_float::scalar(a, R(0), tol, tol, dig)) return false;
      auto is_not_empty = [tol, dig, a](auto p){
        auto pa = p.area();
        return !approx_float::scalar(pa / (pa + a), R(0), tol, tol, dig);
      };
      auto overlap = intersection(that, tol, dig);
      return std::any_of(overlap.begin(), overlap.end(), is_not_empty);
    }

    [[nodiscard]] std::vector<Poly<T,A>> intersection(const Poly<T,A>& that, T tol = T(0), int dig = 0) const {
      if (is_convex() && that.is_convex()){
        return {wires_convex_intersection(vertices_, wires_.border(), that.vertices(), that.wires().border())};
      }
      // simple-case checks first: ignore the possibility of negative polygon inclusion?
      auto that_in_this = border_contains(that.vertices());
      if (std::all_of(that_in_this.begin(), that_in_this.end(), [](const auto x){return x;})) {
        if (that.area() > T(0) && that.wires().wire_count() == 0) return {*this};
        return {insert_hole(that.simplify(), tol, dig)};
      }
      auto this_in_that = that.border_contains(vertices_);
      if (std::all_of(this_in_that.begin(), this_in_that.end(), [](const auto x){return x;})) {
        if (area() > T(0) && wires_.wire_count() == 0) return {that};
        return {that.insert_hole(simplify(), tol, dig)};
      }
      if (wires_.wire_count() || that.wires().wire_count()){
        return networks_intersection(triangulate(), that.triangulate());
      }
      return wires_intersection(vertices_, wires_.border(), that.vertices(), that.wires().border());
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
    void check_wire_indexing() const {
      // Ensure that only indexable vertices are included in the wires -- protect against user error.
      if (wires_.indexing_error(vertices_)) {
        debug_update("Maximum vertex index = ", vertices_.size(0u), " but wires are ", wires_.border());
        if (wires_.wire_count()) {
          for (const auto & w: wires_.wires()) debug_update(w);
        }
        throw std::runtime_error("Check input wire information for indexing error");
      }
    }
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
      auto r = wa.circumscribed_radius(va) + wb.circumscribed_radius(vb);
      auto v = wa.centroid(va) - wb.centroid(vb);
      auto d2 = dot(v, v).val(0, 0);
      return d2 > r * r;
    }


    template<class T, template<class> class A>
    bool
    contains(const A<T> &v0, const Wire &w0, const A<T> &v1, const Wire &w1) {
      auto pts = v1.extract(w1);
      auto one_in_zero = w0.contains(pts, v0);
      // if all are contained, then the whole is contained
      return std::all_of(one_in_zero.begin(), one_in_zero.end(), [](const auto x) { return x; });
    }

    template<class T, template<class> class A> bool
    contains_or_on_border(const A<T> &v0, const Wire &w0, const A<T> &v1, const Wire &w1){
      auto pts = v1.extract(w1);
      auto in = w0.contains(pts, v0);
      auto on = w0.is_on(pts, v0);
      auto all_either = [](const auto & x, const auto & y){
        return std::transform_reduce(x.begin(), x.end(), y.begin(), true,
                                     /*reducer*/ [](const auto a, const auto b){return a && b;},
                                     /*transformer*/ [](const auto a, const auto b){return a || b;}
                                     );
      };
      auto any_in = std::any_of(in.begin(), in.end(), [](const auto x){return x;});
      auto all_on = std::all_of(on.begin(), on.end(), [](const auto x){return x;});
      // rounding errors can cause a point to be both in and on?!
      debug_update_if(all_on, "The wires described by\n", v0.to_string(), w0, "\nand\n", v1.to_string(), w1,
                      "\nare coincident?\n", "any_in=", any_in, " all_either=", all_either(in, on));
      return any_in && !all_on && all_either(in, on);
    }

    template<class T, template<class> class A> bool
    second_in_first(const A<T> &v0, const Wire &w0, const A<T> &v1, const Wire &w1){
      auto pts = v1.extract(w1);
      auto in = w0.contains(pts, v0);
      auto on = w0.is_on(pts, v0);
      auto all_either = [](const auto & x, const auto & y){
        return std::transform_reduce(x.begin(), x.end(), y.begin(), true,
            /*reducer*/ [](const auto a, const auto b){return a && b;},
            /*transformer*/ [](const auto a, const auto b){return a || b;}
        );
      };
      auto any_in = std::any_of(in.begin(), in.end(), [](const auto x){return x;});
      return any_in && all_either(in, on);
    }

    template<class T, template<class> class A> bool
    equivalent(const A<T> & v0, const Wire & w0, const A<T> & v1, const Wire & w1){
      auto all = [](const auto & x){return std::all_of(x.begin(), x.end(), [](const auto y){return y;});};
      auto p1 = v1.extract(w1);
      auto on = w0.is_on(p1, v0);
      if (all(on)){
        auto p2 = v0.extract(w0);
        on = w1.is_on(p2, v1);
        return all(on) && p2.row_is(cmp::eq, p1).all();
      }
      return false;
    }
  }


  /// \brief Find all intersection points for two polygon wires, then return their intersection polygons
  template<class T, template<class> class A>
  std::vector<Poly<T, A>>
  wires_intersection(const A<T> & va, const Wire & wa, const A<T> & vb, const Wire & wb){
    // If the two polygons are too far apart, return nothing
    if (utils::could_not_overlap(va, wa, vb, wb)) return {};
    // If A and B are the same polygon, the following would fail so return either one now
    if (utils::equivalent(vb, wb, va, wa)){
      return {Poly<T, A>(va, wa)};
    }
    // If any of A is inside B and all of A is inside _or_ on the border of B, return A
    if (utils::contains_or_on_border(vb, wb, va, wa)) {
      return {Poly<T, A>(va, wa)};
    }
    // If any of B is inside A and all of B is inside _or_ on the border of A, return b
    if (utils::contains_or_on_border(va, wa, vb, wb)) {
      return {Poly<T, A>(vb, wb)};
    }
    // Otherwise, find the intersection points:
    // Combine the two sets of vertices into a single list, but keep only unique vertices
    auto v = cat(0, va.extract(wa), vb.extract(wb)).unique();
    // Make the new wires for A and B, since we may have reduced or re-ordered their vertices:
    Wire new_wa, new_wb;
    new_wa.reserve(wa.size());
    new_wb.reserve(wb.size());
    for (const auto & idx: wa){
      new_wa.push_back(v.row_is(cmp::eq, va.view(idx)).first());
    }
    for (const auto & idx: wb){
      new_wb.push_back(v.row_is(cmp::eq, vb.view(idx)).first());
    }
    // Construct the dual doubly-linked lists of vertex indices
    auto lists = clip::VertexLists(new_wa, new_wb);
    v = clip::weiler_atherton(v, lists); // updates v with intersection points
    debug_update("Result of Weiler-Atherton, vertices\n", v.to_string(), "lists:\n", lists);
    auto wires = lists.intersection_wires();
    std::vector<Poly<T, A>> result;
    result.reserve(wires.size());
    for (const auto & w: wires) {
      auto p = Poly<T, A>(v, w).without_extraneous_vertices();
      if (p.area()) result.push_back(p);
    }
    if (wires.empty()){
      // if there was no intersection found, we need to make sure one is not inside the other
      // now allowing for the possibility of all points being on the boundary
      if (utils::second_in_first(vb, wb, va, wa)) return {Poly<T, A>(va, wa)};
      if (utils::second_in_first(va, wa, vb, wb)) return {Poly<T, A>(vb, wb)};
    }
    return result;
  }

  /// \brief Find and return the intersection polygon for two convex polygon wires
  template<class T, template<class> class A>
  Poly<T, A>
  wires_convex_intersection(const A<T> & va, const Wire & wa, const A<T> & vb, const Wire & wb){
    auto [vertices, wire] = clip::orourke_chien_olsen_naddor(va, wa, vb, wb);
    return Poly(vertices, wire);
  }

  template<class T, template<class> class A>
  Poly<T, A>
  polygon_convex_intersection(const Poly<T, A> & a, const Poly<T, A> & b){
    if (!a.is_convex() || !b.is_convex()){
      throw std::logic_error("Provided polygons are not both convex!");
    }
    return wires_convex_intersection(a.vertices(), a.wires().border(), b.vertices(), b.wires().border());
  }

  /// \brief Find the intersection of two polygons, triangulating if either has any holes
  template<class T, template<class> class A>
  std::vector<Poly<T, A>>
  polygon_intersection(const Poly<T, A> & a, const Poly<T, A> & b) {
    if (a.wires().wire_count() || b.wires().wire_count()){
      return networks_intersection(a.triangulate(), b.triangulate());
    }
    return wires_intersection(a.vertices(), a.wires().border(), b.vertices(), b.wires().border());
  }

  /// \brief Find the intersection of two triangulated polygon networks by considering all pairs of triangles
  template<class T, template<class> class A, class W=Wire>
  std::vector<Poly<T, A>>
  networks_intersection(const Network<W, T, A> & a, const Network<W, T, A> & b) {
    std::vector<Poly<T, A>> result;
    result.reserve(std::max(a.size(), b.size()));
    for (const auto & wa: a.wires()) {
      for (const auto & wb: b.wires()) {
        auto ps = wires_intersection(a.vertices(), wa, b.vertices(), wb);
        for (const auto & p: ps) if (p.area()) result.push_back(p);
      }
    }
    // TODO go through the resulting polygons and combine any that are touching?
    return result;
  }


  /// \brief Find all intersection points for two polygon wires, then return their union polygons
  template<class T, template<class> class A>
  std::vector<Poly<T, A>>
  wires_union(const A<T> & va, const Wire & wa, const A<T> & vb, const Wire & wb){
    // If the two polygons are too far apart, return both
    if (utils::could_not_overlap(va, wa, vb, wb)) return {{va, wa}, {vb, wb}};
    // If all of A is inside B, return B
    auto a_in_b = wb.contains(va.extract(wa), vb);
    if (utils::contains(vb, wb, va, wa)) return {{vb, wb}};
    // If all of B is inside A, return A
    if (utils::contains(va, wa, vb, wb)) return {{va, wa}};
    // Otherwise, find the intersection points:
    // Combine the two sets of vertices into a single list
    auto v = cat(0, va.extract(wa), vb.extract(wb)).unique();
    // Make the new wires for A and B, since we may have reduced or re-ordered their vertices:
    Wire new_wa, new_wb;
    new_wa.reserve(wa.size());
    new_wb.reserve(wb.size());
    for (const auto & idx: wa){
      new_wa.push_back(v.row_is(cmp::eq, va.view(idx)).first());
    }
    for (const auto & idx: wb){
      new_wb.push_back(v.row_is(cmp::eq, vb.view(idx)).first());
    }
    // Construct the dual doubly-linked lists of vertex indices
    auto lists = clip::VertexLists(new_wa, new_wb);
    v = clip::weiler_atherton(v, lists);
//    std::cout << "Result of Weiler-Atherton, vertices:\n" << v.to_string() << "lists:\n" << lists << "\n";
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
  polygons_union(const Poly<T, A> & a, const Poly<T, A> & b) {
    if (a.wires().wire_count() || b.wires().wire_count()){
      return networks_union(a.triangulate(), b.triangulate());
    }
    return wires_union(a.vertices(), a.wires().border(), b.vertices(), b.wires().border());
  }

  /// \brief Find the intersection of two triangulated polygon networks by considering all pairs of triangles
  template<class T, template<class> class A, class W=Wire>
  std::vector<Poly<T, A>>
  networks_union(const Network<W, T, A> & a, const Network<W, T, A> & b) {
    std::vector<Poly<T, A>> result;
    result.reserve(std::max(a.size(), b.size()));
    for (const auto & wa: a.wires()) {
      for (const auto & wb: b.wires()) {
        auto ps = wires_union(a.vertices(), wa, b.vertices(), wb);
        for (const auto & p: ps) if (p.area()) result.push_back(p);
      }
    }
    // TODO go through the resulting polygons and combine any that are touching?
    return result;
  }


}

#endif