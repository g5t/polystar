#ifndef POLYSTAR_POLYGON_NETWORK_HPP
#define POLYSTAR_POLYGON_NETWORK_HPP
#include <memory>
#include <map>
#include <vector>
#include <algorithm>
//#include "polygon.hpp"
//#include "polygon_wires.hpp"
//#include "polygon_wire.hpp"

namespace polystar::polygon {
  // pre-declare interdependent classes
  class Wire;
  class Wires;
  template<class T, template<class> class A> class Poly;
  Wire wire_merge(const Wire &, const Wire &);

    template<class wire_t>
  class NetworkCost {
  public:
      using cost_t = double;
      using path_t = std::vector<ind_t>;
      using elem_t = std::pair<cost_t, path_t>;
      using map_t = std::map<wire_t, elem_t>;
  private:
      map_t map_;
  public:
      bool known(wire_t p) const {
          return std::find_if(map_.begin(), map_.end(), [p](const auto & x){
              const auto & [w, z] = x;
              return p == w;
          }) != map_.end();
      }
      elem_t cost(wire_t p) const {
          auto ptr = std::find_if(map_.begin(), map_.end(), [p](const auto & x){
              const auto & [w, z] = x;
              return p == w;
          });
          if (ptr != map_.end()) {
              const auto & [w, z] = *ptr;
              return z;
          }
          return std::make_pair<cost_t, path_t>(-1, path_t());
      }
      const elem_t & operator[](const wire_t & w) const {return map_[w];}
      elem_t & operator[](const wire_t & w) {return map_[w];}
  };

  template<class T, template<class> class A>
  class Network {
  public:
    using wire_t = std::shared_ptr<Wire>;
    using ref_t = std::weak_ptr<Wire>;
    using list_t = std::vector<ref_t>;
    using map_t = std::map<wire_t, list_t>;
    using vertex_t = A<T>;
    using edge_t = std::pair<ind_t, ind_t>;
  private:
    map_t map_;
    vertex_t vertices_;
    std::vector<edge_t> external_;
  public:
    Network(const std::vector<Wire> & wires, const Array2<T> & v): vertices_(v) {
      for (const auto & w: wires){
        auto sw = std::make_shared<Wire>(w);
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

    [[nodiscard]] size_t size() const { return map_.size(); }

    [[nodiscard]] std::vector<Wire> wires() const {
      std::vector<Wire> w;
      w.reserve(size());
      for (const auto & [x, l]: map_) w.emplace_back(*(x.get()));
      return w;
    }

    std::vector<Poly<T,A>> polygons() const {
      std::vector<Poly<T,A>> p;
      p.reserve(size());
      for (const auto & w: wires()) {
        auto ws = Wires(w);
        p.emplace_back(vertices_, ws);
      }
      return p;
    }

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

    Network<T,A> simplify();

    bool search_path(NetworkCost<wire_t> & costs, list_t border, wire_t goal, const A<T>& in_first) const {
        while (!border.empty()){
            // find the cheapest border polygon
            std::partial_sort(border.begin(), border.begin()+1, border.end(),
                              [&](const auto & a, const auto & b) -> bool {
                auto a_cost = costs[a.lock()];
                auto b_cost = costs[b.lock()];
                return (a_cost.first < b_cost.first);
            });
            const auto cbp{border.front().lock()};
            // shortcut once we've reached the goal polygon with a cost
            if (goal == cbp) return true;
            // remove the cheapest from the border
            border.erase(border.begin());

            // get the connected polygons from the cheapest
            auto next = neighbours(cbp);
            // by definition, we can reach all connected polygons, but we do not want to backtrack
            next.erase(std::remove_if(next.begin(), next.end(),
                                      [&](const auto & np){return costs.known(np.lock());}), next.end());

            // grab the cost information up to the current location
            auto [cc, cp] = costs[cbp];
            // now decide which of next should be visited
            for (auto & n: next){
                // find the connecting edge between where the cheapest border and this next polygon
                auto edge = common_edge(cbp, n.lock());
                auto e0 = vertices_.view(edge.first);
                auto e1 = vertices_.view(edge.second);
                // figure out the range of costs to move to the extremes of the edge.
                auto d0 =  cp.empty() ? norm(e0 - in_first).sum() : norm(e0 - vertices_.view(cp.back())).sum();
                auto d1 =  cp.empty() ? norm(e1 - in_first).sum() : norm(e1 - vertices_.view(cp.back())).sum();
                // the cost information then is ...
                std::vector<ind_t> path;
                path.reserve(cp.size()+1);
                std::copy(cp.begin(), cp.end(), std::back_inserter(path));
                path.push_back(d0 < d1 ? edge.first : edge.second);
                auto cost_info = std::make_pair(std::min(d0, d1), path);
                // store the cost to get to n
                if (auto nl=n.lock()) if (!costs.known(nl)) /* always true due to the filter above */ {
                  costs[nl] = cost_info;
                  // and add the weak_ptr to the border
                  if (std::find_if(border.begin(), border.end(), [&](const auto & b){
                    if (auto bl = b.lock()) return bl == nl;
                    return false;
                  }) == border.end()) border.push_back(n);
                }
            }
        }
        return false;
    }

    vertex_t path(const vertex_t & from, const vertex_t & to) const {
      auto first = containing_wire(to);
      auto last = containing_wire(from);
      std::cout << "A polygon " << (first.has_value() ? "does" : "does not") << " contain the final point\n";
      std::cout << "A polygon " << (last.has_value() ? "does" : "does not") << " contain the first point\n";
      if (!first.has_value() || !last.has_value()) return from;

      NetworkCost<wire_t> costs;
      list_t border;
      border.emplace_back(first.value()); // make a weak_ptr from the shared_ptr
      costs[first.value()] = std::make_pair<double, std::vector<ind_t>>(0, {});
      auto success = search_path(costs, border, last.value(), to);

      if (!success) std::cout << "No connection between points?\n";

      auto path = costs[last.value()].second;
      auto out = from.decouple();

      std::cout << "found vertex index path ";
      for (const auto & p: path) std::cout << p << ", ";
      std::cout << "\n";

      // cost_info is pair<pair<double, double>, pair<vector<ind_t>, vector<ind_t>>>
      // the vectors contain the edge-pair vertex indices in reverse order for the path that should be taken.

      edge_t edge{0,1};
      for (size_t i=path.size(); i-->0;) {
        // try to skip all the way to the end from here
        auto pg = cat(0, out.view(out.size(0) - 1), to);
        // only add points in the event to the end wasn't possible
        if (!contains(pg, edge)) {
          std::cout << "We didn't skip to the end\n";
          // try to skip over as many edges as we can
          i = skip_to(out, i, path);
          std::cout << "next point " << i << " is " << path[i] << "\n";
          out = cat(0, out, vertices_.view(path[i]));
        }
      }
      // whether we managed to skip any points, add the final point:
      out = cat(0, out, to);
      return out;
    }

  private:
    size_t skip_to(A<T>& far, const size_t i, const std::vector<ind_t>& path) const {
      std::cout << "skip_to with i=" << i << "\n";
      if (i==0) return i;
      std::vector<std::pair<size_t, bool>> skip;
      skip.reserve(i);
      edge_t e{0, 1};
      for (size_t j=i; j-->0;){
        auto pts = cat(0, far.view(far.size(0)-1), vertices_.view(path[j]));
        skip.emplace_back(j, contains(pts, e));
      }
      std::cout << "skip: ";
      for (const auto & p: skip) std::cout << p.first << ":" << (p.second ? "T" : "F") << " ";
      std::cout << "\n";
      // we can iterate forwards through skip because we constructed it backwards
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
  };
}

#endif