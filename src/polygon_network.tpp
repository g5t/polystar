#ifndef POLYSTAR_POLYGON_NETWORK_TPP
#define POLYSTAR_POLYGON_NETWORK_TPP

static
std::tuple<std::array<polystar::ind_t,3>, std::array<polystar::ind_t,3>>
merged_corners(const polystar::polygon::Wire & a, const polystar::polygon::Wire & b) {
    using namespace polystar;
    using c_t = std::array<ind_t, 3>;
    auto as = static_cast<ind_t>(a.size());
    auto bs = static_cast<ind_t>(b.size());
    for (ind_t a1=0; a1<as; ++a1){
        auto a0 = (a1 + as - 1) % as;
        auto a2 = (a1 + 1) % as;
        auto a3 = (a1 + 2) % as;
        for (ind_t b1=0; b1<bs; ++b1){
            auto b0 = (b1 + bs - 1) % bs;
            auto b2 = (b1 + 1) % bs;
            auto b3 = (b1 + 2) % bs;
            if (a1 == b2 && a2 == b1){
                c_t va{{a0, a1, b3}};
                c_t vb{{b0, b1, a3}};
                return std::make_tuple(va, vb);
            }
        }
    }
    return std::make_tuple(c_t(), c_t());
}

template<class T, template<class> class A>
bool would_be_convex(const polystar::polygon::Wire & a, const polystar::polygon::Wire & b, const A<T>& x){
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
  list_t merged_away;
  auto is_in = [](const auto & v, const auto & x){
    return std::find(v.begin(), v.end(), x) != v.end();
  };
  auto remove_and_replace = [&](const auto & o, const auto & remove, const auto & replace){
    auto l = map_[o];
    l.erase(std::remove_if(l.begin(), l.end(), [&](const auto & x){return x==remove;}), l.end());
    l.push_back(replace);
  };
  for (auto ptr = map_.begin(); ptr != map_.end(); ++ptr){
    auto & [w, l] = *ptr;
    // skip this one if it has already been merged into another
    if (is_in(merged_away, w)) continue;
    // try combining this wire with its connected wires
    for (const auto & t: l) if (!is_in(merged_away, t)) if (would_be_convex(*(w.get()), *(t.get()), vertices_)){
      auto merged = polystar::polygon::wire_merge(*(w.get()), *(t.get()));
      // We can combine the wires -- which requires combining their connected lists too
      auto tl = map_[t];
      std::copy_if(tl.begin(), tl.end(), std::back_inserter(l), [&](const auto & x){return !is_in(l, x);});
      // Keep pointers to the now-merged wires, so we don't attempt to merge them again
      merged_away.push_back(w);
      merged_away.push_back(t);
      // go through the connections to this wire *and* that wire, remove pointers to the pre-merged wires
      auto new_w = std::make_shared<Wire>(merged);
      std::for_each(l.begin(), l.end(), [&](const auto & x){remove_and_replace(x, w, new_w);});
      std::for_each(tl.begin(), tl.end(), [&](const auto & x){remove_and_replace(x, t, new_w);});
      // remove the reference to the other wire from this list
      l.erase(std::remove_if(l.begin(), l.end(), [&](const auto & x){return x==t;}), l.end());
      // remove the other wire and its reference list from the map
      erase_wire(t);
      // add the new wire with the merged result
      map_[new_w] = l;
      // remove the current pre-merged wire
      map_.erase(ptr);
    }
  }
//      for (auto & [w, l] : map_){
//        // skip this one if it has already been merged into another
//        if (is_in(merged_away, w)) continue;
//        // try combining this wire with its connected wires
//        for (const auto & t: l){
//          if (is_in(merged_away, t)) continue;
//          auto test = merge(*(w.get()), *(t.get()));
//          if (test.is_convex(vertices_)){
//            // We can combine the wires -- which requires combining their connected lists too
//            auto tl = map_[t];
//            std::copy_if(tl.begin(), tl.end(), std::back_inserter(l), [&](const auto & x){return !is_in(l, x);});
//            // Keep pointers to the now-merged wires, so we don't attempt to merge them again
//            merged_away.push_back(w);
//            merged_away.push_back(t);
//            // go through the connections to this wire *and* that wire, remove pointers to the pre-merged wires
//            auto new_w = std::make_shared<Wire>(test);
//            std::for_each(l.begin(), l.end(), [&](const auto & x){remove_and_replace(x, w, new_w);});
//            std::for_each(tl.begin(), tl.end(), [&](const auto & x){remove_and_replace(x, t, new_w);});
//            // remove the reference to the other wire from this list
//            l.erase(std::remove_if(l.begin(), l.end(), [&](const auto & x){return x==t;}), l.end());
//            // remove the other wire and its reference list from the map
//            erase_wire(t);
//            // replace this wire with the merged result
//            std::swap(w, new_w);
//            // and
//          }
//        }
//      }
  return !merged_away.empty();
}

template<class T, template<class> class A>
std::optional<typename polystar::polygon::Network<T,A>::wire_t> polystar::polygon::Network<T,A>::containing_wire(const typename polystar::polygon::Network<T,A>::vertex_t & x){
    for (const auto & [w, l]: map_) {
        auto x_in = w.get()->contains(x, vertices_, false);
        if (std::count(x_in.begin(), x_in.end(), true) == static_cast<size_t>(x.size(0))) return std::make_optional(w);
    }
    return std::nullopt;
}

#endif
