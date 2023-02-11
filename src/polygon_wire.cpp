#include "polygon_wire.hpp"

using namespace polystar;
using namespace polystar::polygon;

Wire polystar::polygon::wire_merge(const Wire & a, const Wire & b) {
  typename Wire::edge_t ei{0,0}, ej{0,0};
  size_t i{0}, j{0};
  bool found{false};
  for (i=0; i<a.size(); ++i){
    ei = a.edge(i);
    for (j=0; j<b.size(); ++j) {
      ej = b.edge(j);
      found = ei.first == ej.second && ei.second == ej.first;
      if (found) break;
    }
    if (found) break;
  }
  // join together the edges if a.size()==6, b.size()==10, i==3 and j==5, then we want
  // a[0], a[1], a[2],                                                 a[3], a[4], a[5]
  //             b[6], b[7], b[8], b[9], b[0], b[1], b[2], b[3], b[4], b[5]
  // ----------------------------------------------------------------------------------
  // a[0], a[1], a[2], b[7], b[8], b[9], b[0], b[1], b[2], b[3], b[4], a[3], a[4], a[5]
  Wire out;
  out.reserve(a.size() + b.size() - 2);
  for (size_t q=0; q<=i; ++q) out.push_back(a[q]);
  for (size_t q=2; q<b.size(); ++q) out.push_back(b[(q+j)%b.size()]); // start at 2 to skip the edge
  for (size_t q=i+1; q<a.size(); ++q) out.push_back(a[q]); // we already included point i above
  return out;
}