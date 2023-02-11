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

//  template<class C, class P>
//  size_t erase_if(C c, P pred){
//    auto old_size = c.size();
//    for (auto i = c.begin(), last = c.end(); i != last; ) {
//      if (pred(*i)) {
//        i = c.erase(i);
//      } else {
//        ++i;
//      }
//    }
//    return old_size - c.size();
//  }

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

    size_t size() const { return map_.size(); }

    std::vector<Wire> wires() const {
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

    bool simplify();
  };
}

#endif