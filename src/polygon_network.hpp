#ifndef POLYSTAR_POLYGON_NETWORK_HPP
#define POLYSTAR_POLYGON_NETWORK_HPP
#include <memory>
#include <map>
#include <vector>
#include <algorithm>
#include <optional>
#include "graph.hpp"
#include "svg.hpp"
#include "geometry.hpp"

namespace polystar::polygon {
  template<class W, class T, template<class> class A>
  class Network {
  public:
    using wire_t = std::shared_ptr<W>;
    using ref_t = std::weak_ptr<W>;
    using list_t = std::vector<ref_t>;
    using map_t = std::map<wire_t, list_t>;
    using vertex_t = A<T>;
    using edge_t = std::pair<ind_t, ind_t>;
    using off_t = std::optional<W>;
  private:
    map_t map_;
    vertex_t vertices_;
    std::vector<edge_t> external_;
    off_t off_limits_=std::nullopt;
  public:
    Network(const std::vector<W> & wires, const Array2<T> & v): vertices_(v) {
      for (const auto & w: wires){
        auto sw = std::make_shared<W>(w);
        list_t wl;
        wl.reserve(wires.size()-1);
        for (auto & [ow, ol]: map_){
          if (ow.get()->shares_edge(w)){
            ol.emplace_back(sw);
            wl.emplace_back(ow);
          }
        }
        map_.emplace(sw, wl);
      }
      find_external_edges();
    }
//    void off_limits(const W & of) {
//      // only accept setting the off limits positions if it does not remove all known positions
//      for (const auto & i: indexes()) if (std::find(of.begin(), of.end(), i) == of.end()) {
//        off_limits_=of;
//        break;
//      }
//    }
    void off_limits(const W &) {
    }

    [[nodiscard]] size_t size() const { return map_.size(); }


    [[nodiscard]] std::vector<W> wires() const {
      std::vector<W> w;
      w.reserve(size());
      for (const auto & [x, l]: map_) w.emplace_back(*(x.get()));
      return w;
    }
    vertex_t vertices() const {return vertices_;}

    size_t erase_wire(wire_t w){
      auto old_size = map_.size();
      for (auto i = map_.begin(), last = map_.end(); i != last; ) {
        auto const & [wire, list] = *i;
        if (wire == w) {
          i = map_.erase(i);
        } else {
          ++i;
        }
      }
      return old_size - map_.size();
    }

    list_t neighbours(const wire_t & w) const {
      for (const auto &[p, l]: map_) if (p == w) return l;
      return list_t();
    }

    std::optional<wire_t> containing_wire(const vertex_t & x) const;

    Network<W,T,A> simplify();

    vertex_t path(const vertex_t & from, const vertex_t & to) const {
      auto last = containing_wire(to);
      auto first = containing_wire(from);
//      std::cout << "A polygon " << (first.has_value() ? "does" : "does not") << " contain the first point\n";
//      std::cout << "A polygon " << (last.has_value() ? "does" : "does not") << " contain the last point\n";
      if (!first.has_value() || !last.has_value()) return from;

      auto no_v = vertices_.size(0);
      auto graph = graph::Graph<double>(no_v+2); // no more than 2 + len(vertices_) nodes (probably exactly this)
      auto present = indexes();
      // remove the off-limits positions (if there are any)
      if (off_limits_.has_value() && off_limits_->size()) {
        const auto & ol{*off_limits_};
        present.erase(std::remove_if(present.begin(), present.end(), [&](const auto &i) {
          return std::find(ol.begin(), ol.end(), i) != ol.end();}), present.end());
      }
      for (const auto & i: present) {
        for (const auto & j: present) {
          if (i != j && connected(i, j)) {
//            std::cout << i << "--" << j << ": " << distance(i, j) << "\n";
            graph.add_bidirectional(i, j, distance(i, j));
          }
        }
      }
      // add the from and to connections
      const auto & fv{*first};
      const auto & lv{*last};
      for (const auto & i: present) {
        if (std::find(fv->begin(), fv->end(), i) != fv->end()) {
          // we can only come *from* the source
          graph.add_directional(no_v, i, norm(from - vertices_.view(i)).sum());
        }
        if (std::find(lv->begin(), lv->end(), i) != lv->end()) {
          // we only go *towards* the sink
          graph.add_directional(i, no_v+1, norm(to - vertices_.view(i)).sum());
        }
      }
      auto max_distance = 100* norm(vertices_.min(0) - vertices_.max(0)).sum();
      auto path = graph.shortestPath(no_v, no_v+1, max_distance);
//      std::cout << "Path (in vertices_) ";
//      for (const auto & x: path) std::cout << x << " ";
//      std::cout << "\n";

      // path should start at `from` and end at `to` (which are not in vertices_)
      if (path.front() != no_v) throw std::runtime_error("path doesn't start at the start?");
      if (path.back() != no_v+1) throw std::runtime_error("path does not end at the goal?");
      path.erase(path.begin());
      path.erase(path.end()-1);
      auto out = from.decouple();
      for (size_t i=0; i<path.size(); ++i){
        // see if skipping to the end is possible:
        auto pg = cat(0, out.view(out.size(0) - 1), to);
        if (!contains(pg, {0, 1})){
          // can't skip the whole way -- check if we can skip part way
//          i = skip_ahead_to(out, i, path);
          out = cat(0, out, vertices_.view(path[i]));
        }
      }
      // add the final point
      out = cat(0, out, to);
      return out;
    }

    /** \brief Find all paths and their costs from one vertex to many vertices
     * */
      std::vector<std::pair<double, vertex_t>> path_to_each(const vertex_t & from, const vertex_t & to) const {
        auto first = containing_wire(from);
        std::vector<std::optional<wire_t>> destinations;
        destinations.reserve(to.size(0));
        for (ind_t i=0; i<to.size(0); ++i) destinations.push_back(containing_wire(to.view(i)));
        if (!first.has_value()) {
          std::vector<std::pair<double, vertex_t>> no_paths;
          no_paths.reserve(to.size(0));
          for (ind_t i=0; i<to.size(0); ++i) no_paths.emplace_back((std::numeric_limits<double>::max)(), vertex_t());
          return no_paths;
        }
        auto no_v = vertices_.size(0);
        auto graph = graph::Graph<double>(no_v+1+destinations.size());
        auto present = indexes();
        // remove the off-limits positions (if there are any)
        if (off_limits_.has_value() && off_limits_->size()) {
          const auto & ol{*off_limits_};
          present.erase(std::remove_if(present.begin(), present.end(), [&](const auto &i) {
            return std::find(ol.begin(), ol.end(), i) != ol.end();}), present.end());
        }
        for (const auto & i: present) {
            for (const auto & j: present) {
                if (i != j && connected(i, j)) {
                    graph.add_bidirectional(i, j, distance(i, j));
                }
            }
        }
        // add the from and to connections
        std::vector<size_t> extra_indexes(1+destinations.size(), no_v);
        const auto & fv{*first};
        for (const auto & i: present) {
            if (std::find(fv->begin(), fv->end(), i) != fv->end()) {
                // we can only come *from* the source
                graph.add_directional(no_v, i, norm(from - vertices_.view(i)).sum());
            }
        }
        size_t destination_count{0};
        for (const auto & last: destinations) {
          extra_indexes[1+destination_count] = no_v + 1 + destination_count;
          if (last.has_value()) {
            const auto &lv{*last};
            for (const auto &i: present) {
              if (std::find(lv->begin(), lv->end(), i) != lv->end()) {
                // we only go *towards* the sink(s)
                graph.add_directional(i, extra_indexes[1+destination_count],
                                      norm(to.view(destination_count) - vertices_.view(i)).sum());
              }
            }
          }
          ++destination_count;
        }
        auto max_distance = 100* norm(vertices_.min(0) - vertices_.max(0)).sum();
        // find *all* paths from the source -- the value from each sink can be followed-backward to the source.
        auto [costs, prev] = graph.find_paths(no_v, max_distance);

        std::vector<std::pair<double, vertex_t>> costs_paths;
        costs_paths.reserve(destinations.size());
        destination_count = 0;
        for (const auto & last: destinations){
          if (last.has_value()) {
            costs_paths.push_back(filter_back_pointer_path(prev, from, to.view(destination_count), extra_indexes[0], extra_indexes[1+destination_count]));
          } else {
            costs_paths.emplace_back((std::numeric_limits<double>::max)(), vertex_t());
          }
          ++destination_count;
        }
        return costs_paths;
      }

    void add_to_svg(SVG::SVG & svg) const {
      for (const auto & [w, l]: map_){
        auto path = svg.add_child<SVG::Path>();
        w->add_to_svg(*path, vertices_);
      }
    }
    SVG::SVG to_svg(const std::optional<std::string_view> fill, const std::optional<std::string_view> stroke) const {
      SVG::SVG svg;
      svg.style("path").set_attr("fill", std::string(fill.value_or("tan")));
      svg.style("path").set_attr("stroke",std::string(stroke.value_or("black")));
      add_to_svg(svg);
      svg.autoscale();
      return svg;
    }
    SVG::SVG to_svg(const std::string_view fill) const {
      return to_svg(std::make_optional(fill), std::nullopt);
    }
    SVG::SVG to_svg() const {
      return to_svg(std::nullopt, std::nullopt);
    }
    template<class... Args>
    void write_svg(const std::string & filename, Args... args) const {
      auto svg = to_svg(args...);
      auto of = std::ofstream(filename);
      of << std::string(svg);
    }

  private:
    template<class I> std::pair<double, vertex_t> filter_back_pointer_path(const std::vector<I> & back_pointers, const vertex_t & first, const vertex_t & last, size_t from, size_t to) const {
        size_t next{to};
        auto reverse = last.decouple();
        while (next != from){
          // update next (which we know we can reach)
          next = back_pointers[next];
          // check if first is reachable from here
          auto pts = cat(0, reverse.view(reverse.size(0)-1), first);
          if (!contains(pts, {0,1})){
            // we can't skip to the end -- so check how far we *can* skip
            size_t beyond{back_pointers[next]};
            while (beyond != from && can_skip_to(reverse, beyond)) {
              next = beyond;
              beyond = back_pointers[next];
            }
            // we *can* move to next
            reverse = cat(0, reverse, vertices_.view(next));
          }
        }
        // reverse now contains the reverse path from `to` to one before `from` -- flip it around and calculate the cost
        auto forward = first.decouple();
        double cost{0};
        for (ind_t i=reverse.size(0); i-->0;){
          cost += norm(reverse.view(i) - forward.view(forward.size(0)-1)).sum();
          forward = cat(0, forward, reverse.view(i));
        }
        return std::make_pair(cost, forward);
    }
    bool can_skip_to(const A<T>& far, const size_t to) const {
      auto pts = cat(0, far.view(far.size(0)-1), vertices_.view(to));
      return contains(pts, {0, 1});
    }
    template<class I>
    size_t skip_ahead_to(A<T>& far, const size_t i, const std::vector<I>& path) const {
      if (i==path.size()-1) return i;
      std::vector<std::pair<size_t, bool>> skip;
      skip.reserve(path.size()-i);
      for (size_t j=i+1; j<path.size(); ++j){
        auto pts = cat(0, far.view(far.size(0)-1), vertices_.view(path[j]));
        skip.emplace_back(j, contains(pts, {0, 1}));
      }
      auto ptr = std::find_if(skip.begin(), skip.end(), [](const auto & s){return s.second;});
      return ptr != skip.end() ? ptr->first : i;
    }

    edge_t common_edge(const wire_t & asp, const wire_t & bsp) const ;
    template<class R>
    bool contains(A<R>& points, const std::pair<ind_t, ind_t> & edge) const {
      // check that the edge does not cross any *EXTERNAL* edges of the network
//      std::cout << "Does the edge " << edge.first << "--" << edge.second << " in\n" << points;
//      std::cout << " intersect with any of ";
//      for (const auto & ee: external_) std::cout << ee.first << "--" << ee.second << ", ";
//      std::cout << "in\n" << vertices_;
      bool contained{true};
      for (const auto & ee: external_) {
        if (intersect2d(points, edge, vertices_, ee, end_type::neither)) {
          contained = false;
        }
      }
//      std::cout << " " << (contained ? "no" : "yes") << "\n";
      return contained;
    }
    void find_external_edges() {
      external_.clear();
      for (const auto & [p, l]: map_){
        std::vector<edge_t> p_edges;
        auto pg = p.get();
        p_edges.reserve(pg->size());
        for (size_t i=0; i<pg->size(); ++i) p_edges.push_back(pg->edge(i));
        for (const auto & o: l){
          auto ce = common_edge(p, o.lock());
          auto ptr = std::find_if(p_edges.begin(), p_edges.end(), [&](const auto & x){
            return (ce.first == x.first && ce.second == x.second) || (ce.second == x.first && ce.first == x.second);});
          if (ptr != p_edges.end()) p_edges.erase(ptr);
        }
        if (!p_edges.empty()) std::copy(p_edges.begin(), p_edges.end(), std::back_inserter(external_));
      }
    }

    bool connected(ind_t i, ind_t j) const {
      for (const auto & [wire, refs]: map_){
        if(std::find(wire->begin(), wire->end(), i) != wire->end()
           && std::find(wire->begin(), wire->end(), j) != wire->end()) return true;
      }
      return false;
    }

    double distance(ind_t i, ind_t j) const {
      return norm(vertices_.view(i) - vertices_.view(j)).sum();
    }

    std::vector<ind_t> indexes() const {
      std::vector<ind_t> out;
      out.reserve(vertices_.size(0));
      for (const auto & [wire, refs]: map_){
        std::copy_if(wire->begin(), wire->end(), std::back_inserter(out), [&](const auto & x){
          return std::find(out.begin(), out.end(), x) == out.end();
        });
      }
      return out;
    }
  };
}

#endif