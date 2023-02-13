#include "polygon_wires.hpp"

using namespace polystar;
using namespace polystar::polygon;

Wires polystar::polygon::remove_extraneous_wires(const Wires & input){
  auto ws = input.wires();
  if (ws.empty()) return Wires(input.border());

  std::vector<bool> needed(ws.size(), true);
  for (size_t i=0; i<ws.size()-1; ++i) if (needed[i]) {
    const auto & a{ws[i]};
    for (size_t j=i+1; j<ws.size(); ++j) if (needed[j]) {
      const auto & b{ws[j]};
      auto count = std::count_if(a.begin(), a.end(), [&b](const auto & val){
        return std::count(b.begin(), b.end(), val) > 0;
      });
      if (count && std::is_permutation(a.begin(), a.end(), b.begin())) {
        needed[i] = is_positive_permutation(a, b);
        needed[j] = false;
      }
    }
  }
  if (std::find(needed.begin(), needed.end(), false) != needed.end()) {
    for (size_t i=0; i<ws.size(); ++i) if (!needed[i]) ws[i].clear();
    ws.erase(std::remove_if(ws.begin(), ws.end(), [](const auto & x){return x.empty();}), ws.end());
  }
  return Wires(input.border(), ws);
}