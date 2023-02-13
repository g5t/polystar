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
      using cost_t = std::pair<double, double>;
      using path_t = std::pair<std::vector<ind_t>, std::vector<ind_t>>;
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
          return std::make_pair(
            std::make_pair<double, double>(-1, -1),
            std::make_pair(std::vector<ind_t>(), std::vector<ind_t>())
            );
      }
      const elem_t & operator[](const wire_t & w) const {return map_[w];}
      elem_t & operator[](const wire_t & w) {return map_[w];}
  };

  template<class T, template<class> class A>
  class Network {
  public:
    using wire_t = std::shared_ptr<Wire>;
    using list_t = std::vector<std::weak_ptr<Wire>>;
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

    bool simplify();

    bool search_path(NetworkCost<wire_t> & costs, std::vector<wire_t>&  border, wire_t goal, const A<T>& in_first) const {
        while (!border.empty()){
            // find the cheapest border polygon
            std::partial_sort(border.begin(), border.begin()+1, border.end(),
                              [&](const auto & a, const auto & b) -> bool {
                auto a_cost = costs[a];
                auto b_cost = costs[b];
                if (a_cost.first.second <= b_cost.first.first) return true;
                if (a_cost.first.first >= b_cost.first.second) return false;
                return a_cost.first.first < b_cost.first.first;
            });
            const auto cbp{border.front()};
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
                double closest{-1}, farthest{-1};
                // figure out the range of costs to move to the extremes of the edge.
                if (cp.first.empty() || cp.second.empty()){
                    // we're in the first polygon -- calculate distances from in_first to the edge points
                    auto d0 = norm(e0 - in_first).sum();
                    auto d1 = norm(e1 - in_first).sum();
                    closest = std::min(d0, d1);
                    farthest = std::max(d0, d1);
                } else {
                  auto c0 = vertices_.view(cp.first.back());
                  auto c1 = vertices_.view(cp.second.back());
                  auto d00 = norm(e0 - c0).sum();
                  auto d01 = norm(e0 - c1).sum();
                  auto d10 = norm(e1 - c0).sum();
                  auto d11 = norm(e1 - c1).sum();
                  closest = std::min(std::min(d00, d01), std::min(d10, d11));
                  farthest = std::max(std::max(d00, d01), std::max(d10, d11));
                }
                // the cost information then is ...
                std::vector<ind_t> short_path, long_path;
                short_path.reserve(cp.first.size()+1);
                long_path.reserve(cp.second.size()+1);
                std::copy(cp.first.begin(), cp.first.end(), std::back_inserter(short_path));
                std::copy(cp.second.begin(), cp.second.end(), std::back_inserter(long_path));
                short_path.push_back(edge.first);
                long_path.push_back(edge.second);
                auto cost_range = std::make_pair(cc.first + closest, cc.second + farthest);
                auto cost_info = std::make_pair(cost_range, std::make_pair(short_path, long_path));
                // store the cost to get to n
                if (auto nl=n.lock()) if (!costs.known(nl)) /* always true due to the filter above */ {
                  costs[nl] = cost_info;
                  // and add the polygon to the border
                  if (std::find(border.begin(), border.end(), nl) == border.end()) border.push_back(nl);
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
      std::vector<wire_t> border{first.value()};
      std::pair<double,double> distances{0, 0};
      std::pair<std::vector<ind_t>, std::vector<ind_t>> empty_paths{{},{}};
      costs[first.value()] = std::make_pair(distances, empty_paths);
      auto success = search_path(costs, border, last.value(), to);

      if (!success) std::cout << "No connection between points?\n";

      auto paths = costs[last.value()].second;
      auto out = from.decouple();

      auto skip_to = [&](const size_t i, edge_t e){
        std::vector<std::pair<size_t, bool>> skip;
        skip.reserve(i);
        for (size_t j=i; j-->0;){
          auto pts = cat(0, out.view(out.size(0)-1), vertices_.view(paths.first[j]), vertices_.view(paths.second[j]));
          skip.emplace_back(j, contains(pts, e));
        }
        auto ptr = std::find_if(skip.rbegin(), skip.rend(), [](const auto & s){return s.second;});
        return ptr != skip.rend() ? ptr->first : i;
      };

      // cost_info is pair<pair<double, double>, pair<vector<ind_t>, vector<ind_t>>>
      // the vectors contain the edge-pair vertex indices in reverse order for the path that should be taken.

      std::pair<ind_t, ind_t> edge_left{0,1}, edge_right{0,2};
      for (size_t i=paths.first.size(); i-->0;) {
        // try to skip all the way to the end from here
        auto pg = cat(0, out.view(out.size(0) - 1), to);
        if (!contains(pg, edge_left)) {
          // only add points in the event to the end wasn't possible
          if (i == 0) {
            // pick the shorter path to the next edge?
            out = cat(0, out, vertices_.view(paths.first[i]));
          } else {
            // try to skip over as many edges as we can
            bool is_left{false};
            if (auto j = skip_to(i, edge_left); j < i){
              i = j;
              is_left = true;
            } else {
              i = skip_to(i, edge_right);
            }
            out = cat(0, out, vertices_.view(is_left ? paths.first[i] : paths.second[i]));
          }
        }
      }
      // whether we managed to skip any points, add the final point:
      out = cat(0, out, to);
      return out;
    }

  private:
    edge_t common_edge(const wire_t & asp, const wire_t & bsp) const ;
    template<class R>
    bool contains(A<R>& points, const std::pair<ind_t, ind_t> & edge) const {
      // check that the edge does not cross any *EXTERNAL* edges of the network
      for (const auto & ee: external_) if (intersect2d(points, edge, vertices_, ee, end_type::neither)) return false;
      return true;
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