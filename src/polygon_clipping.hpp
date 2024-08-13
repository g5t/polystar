#pragma once
#include "polygon_wire.hpp"

namespace polystar::polygon::clip
{
  using ind_t = polystar::ind_t;

  enum class Type {unknown, entry, exit, original, edge};
  enum class On {neither, A, B, both};

  class Vertex{
  public:
    using ptr = std::shared_ptr<Vertex>;
  private:
    ind_t _value{0};
    Type _type = Type::unknown;
    bool _visited = false;
  protected:
    ptr prev_A, next_A, prev_B, next_B;

  public:
    Vertex() = default;
    explicit Vertex(ind_t i, Type type = Type::original) : _value(i), _type(type) {}

    /// \brief Construct a 'normal' vertex
    Vertex(ind_t i, On on, const ptr& prev, const ptr& next): _value(i), _type(Type::original) {
      if (on == On::A) {
        prev_A = prev;
        next_A = next;
      } else if (on == On::B) {
        prev_B = prev;
        next_B = next;
      }
    }

    /// \param Construct a common (intersection) vertex
    Vertex(ind_t i, Type type, Vertex * prev_A, Vertex * next_A, Vertex * prev_B, Vertex * next_B)
      : _value(i), _type(type), prev_A(prev_A), next_A(next_A), prev_B(prev_B), next_B(next_B) {}

    [[nodiscard]] ind_t value() const { return _value; }
    [[nodiscard]] Type type() const { return _type; }
    [[nodiscard]] bool visited() const { return _visited; }
    void visited(bool v) { _visited = v; }

    [[nodiscard]] bool is_A() const { return next_A != nullptr && prev_A != nullptr; }
    [[nodiscard]] bool is_B() const { return next_B != nullptr && prev_B != nullptr; }
    [[nodiscard]] bool is_Both() const { return is_A() && is_B(); }

    void prev(On on, const ptr& v) {
      if (on == On::A || on == On::both) prev_A = v;
      if (on == On::B || on == On::both) prev_B = v;
    }
    void next(On on, const ptr& v) {
      if (on == On::A || on == On::both) next_A = v;
      if (on == On::B || on == On::both) next_B = v;
    }
    [[nodiscard]] ptr next(On on, Type type = Type::unknown) const {
      // step to the next ptr on A or B
      auto n = next_on(on);
      // if we take any type, or this one happens to be what we want, return it
      if (Type::unknown == type || n->type() == type) return n;
      // keep a reference to know where we started
      auto stop = n;
      do {
        // get the next ptr on A or B
        n = n->next_on(on);
        // return it if it's the right type
        if (n->type() == type) return n;
        // continue until we get a nullptr or back to where we started
      } while (n != nullptr && n != stop);
      return nullptr;
    }
    [[nodiscard]] ptr prev(On on, Type type = Type::unknown) const {
      // step to the prev ptr on A or B
      auto n = prev_on(on);
      // if we take any type, or this one happens to be what we want, return it
      if (Type::unknown == type || n->type() == type) return n;
      // keep a reference to know where we started
      auto stop = n;
      do {
        // get the next ptr on A or B
        n = n->prev_on(on);
        // return it if it's the right type
        if (n->type() == type) return n;
        // continue until we get a nullptr or back to where we started
      } while (n != nullptr && n != stop);
      return nullptr;
    }

  private:
    [[nodiscard]] ptr next_on(On on) const {
      if (on == On::A) return next_A;
      if (on == On::B) return next_B;
      return nullptr;
    }
    [[nodiscard]] ptr prev_on(On on) const {
      if (on == On::A) return prev_A;
      if (on == On::B) return prev_B;
      return nullptr;
    }

  };

  class VertexList{
    Vertex::ptr head;

  public:
    VertexList(const polystar::polygon::Wire & p, On on) {
      if (p.empty()) return;
      head = std::make_shared<Vertex>(p[0]);
      head->next(on, head);
      auto prev = head;
      for (ind_t i = 1; i < p.size(); ++i) {
        auto v = std::make_shared<Vertex>(p[i]);
        // insert the new vertex into the list
        v->prev(on, prev);
        v->next(on, prev->next(on));
        prev->next(on, v);
        // move to the next vertex
        prev = v;
      }
      head->prev(on, prev);
    }

    [[nodiscard]] polystar::polygon::Wire wire(On on) const {
      Wire res;
      auto v = head;
      if (v == nullptr) return res;
      do {
        res.push_back(v->value());
        v = v->next(on);
      } while (v != head);
      return res;
    }

    [[nodiscard]] bool is_empty() const { return head == nullptr || (head == head->next(On::A) && head == head->next(On::B));}
    [[nodiscard]] Vertex::ptr first() const { return head; }
  };

  class VertexLists{
    VertexList A, B;

  public:
    VertexLists(const polystar::polygon::Wire & a, const polystar::polygon::Wire & b): A(a, On::A), B(b, On::B) {}

    [[nodiscard]] Vertex::ptr first(On on) const {
      if (on == On::A) return A.first();
      if (on == On::B) return B.first();
      return nullptr;
    }

    [[nodiscard]] std::vector<Wire> intersection_wires() const;
    [[nodiscard]] std::vector<Wire> union_wires() const;
  };

  template<class T, template<class> class A>
  int weiler_atherton(const A<T> & v, VertexLists & lists){
    int no_problems{0};
    // Walk through both lists of vertices, looking for intersections
    auto first_a = lists.first(clip::On::A);
    auto first_b = lists.first(clip::On::B);
    auto ptr_a = first_a;
    do {
      // get the next original wire vertex on A
      ptr_a = ptr_a->next(clip::On::A, clip::Type::original);
      // Find Edge A
      auto edge_a = std::make_pair(ptr_a->value(), ptr_a->next(clip::On::A, clip::Type::original)->value());
      auto ptr_b = first_b;
      do {
        // get the next original wire vertex on B
        ptr_b = ptr_b->next(clip::On::B, clip::Type::original);
        // Find Edge B
        auto edge_b = std::make_pair(ptr_b->value(), ptr_b->next(clip::On::B, clip::Type::original)->value());
        // Find the intersection point of edge A and edge B
        auto [valid, at] = intersection2d(v, edge_a, v, edge_b, end_type::second); // include only one of the two ends
        //
        if (valid){
          //  1. Add it to the list of all vertices
          auto index = v.size();
          auto match = v.row_is(cmp::eq, at);
          if (match.any()){
            // the intersection point is already in the list of vertices
            index = match.first();
          } else {
            // the intersection point is not in the list of vertices
            index = v.size(0);
            v = v.append(0, at);
          }

          //  2. Find whether it points into or out of A.
          auto a_0 = v.view(edge_a.first);
          auto r = v.view(edge_b.second) - a_0; // we want to see if edge_b.second is to the right of edge_a
          auto s = v.view(edge_a.second) - a_0;
          auto cross = cross2d(r, s).val(0,0);
          auto type = cross > 0 ? clip::Type::entry : cross < 0 ? clip::Type::exit : clip::Type::edge;
          if (clip::Type::edge == type){
            ++no_problems;
            // use the other point on the edge of B to decide if this points in our out
            r = v.view(edge_b.first) - a_0;
            cross = cross2d(r, s).val(0,0);
            type = cross > 0 ? clip::Type::exit : cross < 0 ? clip::Type::entry : clip::Type::edge;
          }
          if (clip::Type::edge == type) ++no_problems;

          //  3. Insert it into the doubly-linked lists of vertices on both A and B
          auto ptr = std::make_shared<clip::Vertex>(index, type);
          ptr->prev(clip::On::A, ptr_a);
          ptr->next(clip::On::A, ptr_a->next(clip::On::A));
          ptr->prev(clip::On::B, ptr_b);
          ptr->next(clip::On::B, ptr_b->next(clip::On::B));
          ptr_a->next(clip::On::A, ptr);
          ptr_b->next(clip::On::B, ptr);
        }
      } while (ptr_b != first_b);
    } while (ptr_a != first_a);
    return no_problems;
  }
}