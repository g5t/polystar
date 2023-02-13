#ifndef POLYSTAR_POLYGON_NETWORK_TPP
#define POLYSTAR_POLYGON_NETWORK_TPP

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

template<class T, template<class> class A>
bool polystar::polygon::Network<T,A>::simplify() {
  std::vector<wire_t> merged_away;
  auto is_in = [](const auto & v, const auto & x){
    return std::find(v.begin(), v.end(), x) != v.end();
  };
  auto weak_is_in = [](const auto & v, const auto & x){
    auto xl = x.lock();
    if (!xl) return false;
    return std::find_if(v.begin(), v.end(), [&](const auto & y){return y.lock() == xl;}) != v.end();
  };
  auto remove_and_replace = [&](const auto & o, const auto & remove, const auto & replace){
    if (auto ol = o.lock()) {
      auto l = map_[ol];
      l.erase(std::remove_if(l.begin(), l.end(), [&](const auto &x) { return x.lock() == remove; }), l.end());
      l.emplace_back(replace);
    }
  };
  for (auto ptr = map_.begin(); ptr != map_.end(); ++ptr){
    auto & [w, l] = *ptr;
    // skip this one if it has already been merged into another
    if (is_in(merged_away, w)) continue;
    // try combining this wire with its connected wires
    for (const auto & t: l) if (auto tp=t.lock()) if (!is_in(merged_away, tp)) if (would_be_convex(w, tp, vertices_)){
      auto merged = polystar::polygon::wire_merge(*(w.get()), *(tp.get()));
      // We can combine the wires -- which requires combining their connected lists too
      auto wl = map_[w];
      auto tl = map_[tp];
      std::copy_if(tl.begin(), tl.end(), std::back_inserter(wl), [&](const auto & x){return !weak_is_in(wl, x);});
      // Keep pointers to the now-merged wires, so we don't attempt to merge them again
      merged_away.push_back(w);
      merged_away.push_back(tp);
      // go through the connections to this wire *and* that wire, remove pointers to the pre-merged wires
      auto new_w = std::make_shared<Wire>(merged);
      std::for_each(wl.begin(), wl.end(), [&](const auto & x){remove_and_replace(x, w, new_w);});
      std::for_each(tl.begin(), tl.end(), [&](const auto & x){remove_and_replace(x, tp, new_w);});
      // remove the reference to the other wire from this list
      wl.erase(std::remove_if(wl.begin(), wl.end(), [&](const auto & x){return x.lock()==tp;}), wl.end());
      // remove the other wire and its reference list from the map
      erase_wire(tp);
      // add the new wire with the merged result
      map_[new_w] = wl;
    }
  }
  // remove the merged polygons
  for (auto ptr = map_.begin(), last = map_.end(); ptr != last; ){
    const auto & [p, l] = *ptr;
    if (is_in(merged_away, p)) {
      ptr = map_.erase(ptr);
    } else {
      ++ptr;
    }
  }
  return !merged_away.empty();
}

template<class T, template<class> class A>
std::optional<typename polystar::polygon::Network<T,A>::wire_t>
  polystar::polygon::Network<T,A>::containing_wire(const typename polystar::polygon::Network<T,A>::vertex_t & x) const {
    for (const auto & [w, l]: map_) {
      auto x_in = w->contains(x, vertices_, end_type::second);
      if (std::find(x_in.begin(), x_in.end(), false) == x_in.end()) return std::make_optional(w);
    }
    return std::nullopt;
}

template<class T, template<class> class A>
  typename polystar::polygon::Network<T,A>::edge_t
  polystar::polygon::Network<T,A>::common_edge(const wire_t & asp, const wire_t & bsp) const {
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
