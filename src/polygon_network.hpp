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

    template<class W>
  class NetworkCost {
  public:
      using wire_t = std::shared_ptr<W>;
      using cost_t = std::pair<double, double>;
      using path_t = std::pair<std::vector<ind_t>, std::vector<ind_t>>;
      using elem_t = std::pair<cost_t, path_t>;
      using map_t = std::map<wire_t, elem_t>;
  private:
      map_t map_;
  public:
      bool known(wire_t p) const {
          return std::find(map_.begin(), map_.end(), [p](const auto & x){
              const auto & [w, z] = x;
              return p == w;
          }) != map_.end();
      }
      elem_t cost(wire_t p) const {
          auto ptr = std::find(map_.begin(), map_.end(), [p](const auto & x){
              const auto & [w, z] = x;
              return p == w;
          });
          if (ptr != map_.end()) {
              const auto & [w, z] = *ptr;
              return z;
          }
          return std::make_pair(std::make_pair<double, double>(-1, -1), std::make_pair(std::vector<ind_t>(), std::vector<ind_t>()));
      }
      const elem_t & operator[](const wire_t & w) const {return map_[w];}
      elem_t & operator[](const wire_t & w) {return map_[w];}
  };

  template<class T, template<class> class A>
  class Network {
  public:
    using wire_t = std::shared_ptr<Wire>;
    using list_t = std::vector<wire_t>;
    using map_t = std::map<wire_t, list_t>;
    using vertex_t = A<T>;
  private:
    map_t map_;
    vertex_t vertices_;
  public:
    Network(const std::vector<Wire> & wires, const Array2<T> & v): vertices_(v) {
      for (const auto & w: wires){
        auto sw = std::make_shared<Wire>(w);
        list_t wl;
        wl.reserve(wires.size()-1);
        for (auto & [ow, ol]: map_){
          if (ow.get()->shares_edge(w)){
            ol.push_back(sw);
            wl.push_back(ow);
          }
        }
        map_.emplace(sw, wl);
      }
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


    std::optional<wire_t> containing_wire(const vertex_t & x);

    bool simplify();

    bool search_path(NetworkCost<wire_t, T, A> & costs, double minimum, double maximum, std::vector<wire_t>&  border, wire_t goal, A<T> short_path, A<T> long_path, const A<T>& in_first, const A<T>& in_last){
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
            auto cbp = border.front();
            // remove the cheapest from the border
            border.erase(border.begin());

            // get the connected polygons from the cheapest
            auto next = map_[cbp];
            // by definition we can reach all connected polygons, but we do not want to backtrack
            next.erase(std::remove_if(next.begin(), next.end(),
                                      [&](const auto & np){return costs.known(np);}), next.end());

            // grab the cost information up to the current location
            auto cbp_cost = costs[cbp];
            // now decide which of next should be visited
            for (auto & n: next){
                // find the connecting edge between where the cheapest border and this next polygon
                auto edge = common_edge(cbp, n);
                auto e0 = vertices_.view(edge.first);
                auto e1 = vertices_.view(edge.second);
                // figure out the range of costs to move to the extremes of the edge.
                if (cbp_cost.second.first.empty() || cbp_cost.second.second.empty()){
                    // we're in the first polygon -- calculate from in_first to the edge points
                    
                }

            }

        }
    }

    vertex_t path(const vertex_t & from, const vertex_t & to) const {
      /* 1. Find the ending poylygon
       * 2. Perform a search until the start is found
       * 3. Remove extraneous points on the path
       *    (if the path is A-B-C-..., and A-C intersects the edge containing B, B is extraneous
       */
      auto first = containing_wire(to);
      auto last = containing_wire(from);
      if (!first.has_value() || !last.has_value()) return from;

      NetworkCost<wire_t, T, A> costs;
      std::vector<wire_t> border{first};
      auto success = search_path(costs, 0, 0, border, last, from.decouple(), from.decouple(), to, from);


      // add all the intervening points
      // add the final point:
      out = cat(0, out, to);
      return out;
    }
  private:
      std::pair<ind_t, ind_t> common_edge(const wire_t & asp, const wire_t & bsp){
          // error checking to ensure that they're in each others' connected lists? No, since this is private
          auto a = asp.get();
          auto b = bsp.get();
          auto as = static_cast<ind_t>(a->size());
          auto bs = static_cast<ind_t>(b->size());
          for (ind_t i=0; i<as; ++i){
              for (ind_t j=0; j<bs; ++j){
                  if (a->operator[](i)  == b->operator[]((j+1)%bs) && b->operator[](j) == a->operator[]((i+1)%as)) {
                      return std::make_pair(a->operator[](i), b->operator[](j));
                  }
              }
          }
          return std::make_pair(0, 0);
      }
  };
}

#endif