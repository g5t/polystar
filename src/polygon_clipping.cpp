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
  // Modified for anticlockwise vertex ordering (start on A instead of B)

  // Walk around wire A, looking for an intersection with wire B that is pointing into A
  auto v = A.first();
  auto first_start = v;
  On on = On::A;
  Wire wire;
  bool recording{false};
  bool started{false};
  auto start = v;
  do {
    // With a latch on whether this is the 'first' go-round, if we return to the start, this polygon is done
    if (started && v == start) {
      // stop recording
      recording = false;
      // go back to A, no matter what
      on = On::A;
      if (!wire.empty()) {
//        std::cout << "Resulting wire is: [ ";
//        for (const auto & w: wire) std::cout << w << ", ";
//        std::cout << "]\n";
        combined.push_back(wire);
        wire = Wire();
      } else {
        // we're stuck in a loop -- punt by returning to the start?
        started = false;
        start = first_start;
        v = first_start;
      }
    }
    //
    if (!recording && !v->visited() && v->vertex_type() == Type::entry) {
//      std::cout << "Enter wire at " << v << "\n";
      // We have found an intersection point, now walk around wire A until we find the exit point
      start = v;
      recording = true;
      started=true;
    }
    if (recording) {
      // only switch between wires if we're following a valid intersection boundary
      if (!v->visited() && v->vertex_type() == Type::entry) on = On::A;
      if (!v->visited() && v->vertex_type() == Type::exit) on = On::B;
      // record the vertex index, protecting against repeated indexes
      // (from input-polygon vertices also being intersection points)
      if (wire.empty() || wire.back() != v->value()) {
//        std::cout << "Recording " << v << "\n";
        wire.push_back(v->value());
      }
      v->visited(true);
//    } else {
//      std::cout << "Passed " << v << "\n";
    }

    v = v->next(on);
  } while (recording || v != first_start);

  return combined;
}


std::vector<Wire> VertexLists::union_wires() const {
  std::vector<Wire> combined;
  if (A.is_empty() && B.is_empty()) return combined;
  if (A.is_empty()) return {B.wire(On::B)};
  if (B.is_empty()) return {A.wire(On::A)};
  // Use the Weiler-Atherton algorithm to combine the two wires
  // https://en.wikipedia.org/wiki/Weiler%E2%80%93Atherton_clipping_algorithm

  // FIXME As with intersections, there can be a closed loop that does not revisit the start vertex

  // Walk around wire B, looking for an intersection with wire A that is pointing into A
  auto v = B.first();
  Wire wire;
  On on = On::B;
  do {
    v->visited(true);
    wire.push_back(v->value());
    if (v->vertex_type() == Type::entry) on = On::A;
    else if (v->vertex_type() == Type::exit) on = On::B;
    v = v->next(on);
  } while (!v->visited() && v != B.first());
  combined.push_back(wire);

  // Walk around wire B again, looking for unvisited intersection vertices (which should be holes)
  v = B.first();
  do {
    if (v->vertex_type() == Type::entry && !v->visited()){
      Wire hole;
      auto start = v;
      do {
        v->visited(true);
        hole.push_back(v->value());
        v = v->prev(On::B);
      } while (v != start && v->vertex_type() != Type::exit);
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