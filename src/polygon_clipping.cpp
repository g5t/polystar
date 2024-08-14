#include "polygon_clipping.hpp"

using namespace polystar::polygon;
using namespace polystar::polygon::clip;


std::vector<Wire> VertexLists::intersection_wires() const {
  std::vector<Wire> combined;
  if (A.is_empty() && B.is_empty()) return combined;
  if (A.is_empty()) return {B.wire(On::B)};
  if (B.is_empty()) return {A.wire(On::A)};
  // Use the Weiler-Atherton algorithm to combine the two wires
  // https://en.wikipedia.org/wiki/Weiler%E2%80%93Atherton_clipping_algorithm
  // https://liorsinai.github.io/mathematics/2023/09/30/polygon-clipping.html

  // Walk around wire B, looking for an intersection with wire A that is pointing into A
  auto v = B.first();
  do {
    v->visited(true);
    if (v->type() == Type::entry){
      // We have found an intersection point, now walk around wire B until we find the exit point
      Wire wire;
      auto start = v;
      do {
        v->visited(true);
        wire.push_back(v->value());
        v = v->next(On::B);
      } while (v != start && v->type() != Type::exit);
      // If we found the exit point, keep walking around A until we find the entry point again
      if (v != start) {
        v = v->next(On::A);
        while (v != start){
          v->visited(true);
          wire.push_back(v->value());
          v = v->next(On::A);
        }
      }
      combined.push_back(wire);
    }
    // Fast-forward over the already-visited vertices
    while (v->visited() && v != B.first()) v = v->next(On::B);
  } while (v != B.first());

  return combined;
}

std::vector<Wire> VertexLists::union_wires() const {
  std::vector<Wire> combined;
  if (A.is_empty() && B.is_empty()) return combined;
  if (A.is_empty()) return {B.wire(On::B)};
  if (B.is_empty()) return {A.wire(On::A)};
  // Use the Weiler-Atherton algorithm to combine the two wires
  // https://en.wikipedia.org/wiki/Weiler%E2%80%93Atherton_clipping_algorithm

  // Walk around wire B, looking for an intersection with wire A that is pointing into A
  auto v = B.first();
  Wire wire;
  On on = On::B;
  do {
    v->visited(true);
    wire.push_back(v->value());
    if (v->type() == Type::entry) on = On::A;
    else if (v->type() == Type::exit) on = On::B;
    v = v->next(on);
  } while (!v->visited() && v != B.first());
  combined.push_back(wire);

  // Walk around wire B again, looking for unvisited intersection vertices (which should be holes)
  v = B.first();
  do {
    if (v->type() == Type::entry && !v->visited()){
      Wire hole;
      auto start = v;
      do {
        v->visited(true);
        hole.push_back(v->value());
        v = v->prev(On::B);
      } while (v != start && v->type() != Type::exit);
      if (v != start) {
        v = v->prev(On::A);
        while (v != start){
          v->visited(true);
          hole.push_back(v->value());
          v = v->prev(On::A);
        }
      }
      combined.push_back(hole);
    }
    v = v->next(On::B);
  } while (v != B.first());

  return combined;
}

std::string polystar::polygon::clip::to_string(Type type){
  switch (type){
    case Type::unknown: return "?";
    case Type::entry: return "v";
    case Type::exit: return "^";
    case Type::original: return "o";
    case Type::edge: return "-";
    default: return "!";
  }
}

std::string polystar::polygon::clip::to_string(On on){
  switch (on){
    case On::neither: return "_";
    case On::A: return "A";
    case On::B: return "B";
    case On::both: return "+";
    default: return "!";
  }
}