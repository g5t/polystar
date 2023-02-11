#ifndef POLYSTAR_POLYGON_NETWORK_TPP
#define POLYSTAR_POLYGON_NETWORK_TPP

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
    for (const auto & t: l){
      if (is_in(merged_away, t)) continue;
      auto test = polystar::polygon::wire_merge(*(w.get()), *(t.get()));
      if (test.is_convex(vertices_)){
        // We can combine the wires -- which requires combining their connected lists too
        auto tl = map_[t];
        std::copy_if(tl.begin(), tl.end(), std::back_inserter(l), [&](const auto & x){return !is_in(l, x);});
        // Keep pointers to the now-merged wires, so we don't attempt to merge them again
        merged_away.push_back(w);
        merged_away.push_back(t);
        // go through the connections to this wire *and* that wire, remove pointers to the pre-merged wires
        auto new_w = std::make_shared<Wire>(test);
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

#endif
