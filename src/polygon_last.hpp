#ifndef POLYSTAR_POLYGON_NETWORK_TPP
#define POLYSTAR_POLYGON_NETWORK_TPP
#include "polygon_network.hpp"
#include "polygon_wire.hpp"

/** \file Implements Network methods that depend on Wire
 *
 * The two classes have interdependent methods, Wire includes Network
 * so Network _can't_ include Wire.
 * This solution works, even if it is somewhat convoluted.
 */

//static
//std::tuple<std::array<polystar::ind_t,3>, std::array<polystar::ind_t,3>>
//merged_corners(const polystar::polygon::Wire & a, const polystar::polygon::Wire & b) {
//    using namespace polystar;
//    using c_t = std::array<ind_t, 3>;
//    auto as = static_cast<ind_t>(a.size());
//    auto bs = static_cast<ind_t>(b.size());
//    for (ind_t i=0; i < as; ++i){
//      auto a0 = a[i];
//      auto a1 = a[(i + 1) % as];
//      auto a2 = a[(i + 2) % as];
//      auto a3 = a[(i + 3) % as];
//      for (ind_t j=0; j < bs; ++j){
//        auto b0 = b[j];
//        auto b1 = b[(j + 1) % bs];
//        auto b2 = b[(j + 2) % bs];
//        auto b3 = b[(j + 3) % bs];
//        if (a1 == b2 && a2 == b1){
//            c_t va{{a0, a1, b3}};
//            c_t vb{{b0, b1, a3}};
//            return std::make_tuple(va, vb);
//        }
//      }
//    }
//    throw std::runtime_error("merged_corners called for non-neighbors?");
//    return std::make_tuple(c_t(), c_t());
//}

template<class W>
std::tuple<std::array<polystar::ind_t,3>, std::array<polystar::ind_t,3>>
merged_corners(const std::shared_ptr<W> asp, const std::shared_ptr<W> bsp) {
  auto a = asp.get();
  auto b = bsp.get();
  using namespace polystar;
  using c_t = std::array<ind_t, 3>;
  auto as = static_cast<ind_t>(a->size());
  auto bs = static_cast<ind_t>(b->size());
  for (ind_t i=0; i < as; ++i){
    auto a1 = a->operator[]((i + 1) % as);
    auto a2 = a->operator[]((i + 2) % as);
    for (ind_t j=0; j < bs; ++j){
      auto b1 = b->operator[]((j + 1) % bs);
      auto b2 = b->operator[]((j + 2) % bs);
      if (a1 == b2 && a2 == b1){
        c_t va{{a->operator[](i), a1, b->operator[]((j + 3) % bs)}};
        c_t vb{{b->operator[](j), b1, a->operator[]((i + 3) % as)}};
        return std::make_tuple(va, vb);
      }
    }
  }
  throw std::runtime_error("merged_corners called for non-neighbors?");
  return std::make_tuple(c_t(), c_t());
}

template<class W>
std::tuple<std::array<polystar::ind_t,3>, std::array<polystar::ind_t,3>>
merged_corners(const W & a, const std::shared_ptr<W> bsp) {
  auto b = bsp.get();
  using namespace polystar;
  using c_t = std::array<ind_t, 3>;
  auto as = static_cast<ind_t>(a.size());
  auto bs = static_cast<ind_t>(b->size());
  for (ind_t i=0; i < as; ++i){
    auto a1 = a.operator[]((i + 1) % as);
    auto a2 = a.operator[]((i + 2) % as);
    for (ind_t j=0; j < bs; ++j){
      auto b1 = b->operator[]((j + 1) % bs);
      auto b2 = b->operator[]((j + 2) % bs);
      if (a1 == b2 && a2 == b1){
        c_t va{{a.operator[](i), a1, b->operator[]((j + 3) % bs)}};
        c_t vb{{b->operator[](j), b1, a.operator[]((i + 3) % as)}};
        return std::make_tuple(va, vb);
      }
    }
  }
  throw std::runtime_error("merged_corners called for non-neighbors?");
  return std::make_tuple(c_t(), c_t());
}


template<class W, class T, template<class> class A>
bool would_be_convex(const std::shared_ptr<W> a, const std::shared_ptr<W> b, const A<T>& x){
  // We already know the neighbours share an edge, otherwise they would not be neighbours.
  // So find the edge, then check for convex-ness of the new corners
  auto [va, vb] = merged_corners(a, b);
  // we've thus far required that every polygon in the network is convex -- we exploit that now
  auto is_convex = [&](const auto & g){
      return cross2d(x.view(g[1]) - x.view(g[0]), x.view(g[2]) - x.view(g[1])).val(0,0) > 0;
  };
  return is_convex(va) && is_convex(vb);
}

template<class W, class T, template<class> class A>
bool would_be_convex(const W & a, const std::shared_ptr<W> b, const A<T>& x){
  // We already know the neighbours share an edge, otherwise they would not be neighbours.
  // So find the edge, then check for convex-ness of the new corners
  auto [va, vb] = merged_corners(a, b);
  // we've thus far required that every polygon in the network is convex -- we exploit that now
  auto is_convex = [&](const auto & g){
    return cross2d(x.view(g[1]) - x.view(g[0]), x.view(g[2]) - x.view(g[1])).val(0,0) > 0;
  };
  return is_convex(va) && is_convex(vb);
}

template<class W, class T, template<class> class A>
polystar::polygon::Network<W,T,A> polystar::polygon::Network<W,T,A>::simplify() {
  std::vector<wire_t> merged_away;
  auto is_in = [](const auto & v, const auto & x){
    return std::find(v.begin(), v.end(), x) != v.end();
  };
  std::vector<Wire> simple_wires;
  for (const auto & [w, l]: map_) if (!is_in(merged_away, w)) {
    auto simple_wire = *(w.get()); // copy
    // try combining this wire with its connected wires
    size_t count{0};
    for (const auto & t: l)
      if (auto tp=t.lock())
        if (!is_in(merged_away, tp))
          if (would_be_convex(simple_wire, tp, vertices_)){
            simple_wire = polystar::polygon::wire_merge(simple_wire, *(tp.get()));
            // Keep pointers to the now-merged wires, so we don't attempt to merge them again
            if (count++ == 0) merged_away.push_back(w);
            merged_away.push_back(tp);
          }
    simple_wires.push_back(simple_wire);
  }
  // use the constructor to sort-out connectivity information
  return polystar::polygon::Network<W,T,A>(simple_wires, vertices_);
}

template<class W, class T, template<class> class A>
std::optional<typename polystar::polygon::Network<W,T,A>::wire_t>
  polystar::polygon::Network<W,T,A>::containing_wire(const typename polystar::polygon::Network<W,T,A>::vertex_t & x) const {
    for (const auto & [w, l]: map_) {
      auto x_in = w->contains(x, vertices_);
      if (std::find(x_in.begin(), x_in.end(), false) == x_in.end()) return std::make_optional(w);
    }
    return std::nullopt;
}

template<class W, class T, template<class> class A>
  typename polystar::polygon::Network<W,T,A>::edge_t
  polystar::polygon::Network<W,T,A>::common_edge(const wire_t & asp, const wire_t & bsp) const {
  // error checking to ensure that they're in each other's connected lists? No, since this is private
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

#endif
