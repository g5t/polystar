#pragma once
#include "polygon_wire.hpp"

namespace polystar::polygon::clip
{
  using ind_t = polystar::ind_t;

  enum class Type {unknown, entry, exit, original, edge};
  std::string to_string(Type type);
  enum class On {neither, A, B, both};
  std::string to_string(On on);

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
    [[nodiscard]] Type vertex_type() const { return _type; }
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
      if (Type::unknown == type || n->vertex_type() == type) return n;
      // keep a reference to know where we started
      auto stop = n;
      do {
        // get the next ptr on A or B
        n = n->next_on(on);
        // return it if it's the right type
        if (n->vertex_type() == type) return n;
        // continue until we get a nullptr or back to where we started
      } while (n != nullptr && n != stop);
      return nullptr;
    }
    [[nodiscard]] ptr prev(On on, Type type = Type::unknown) const {
      // step to the prev ptr on A or B
      auto n = prev_on(on);
      // if we take any type, or this one happens to be what we want, return it
      if (Type::unknown == type || n->vertex_type() == type) return n;
      // keep a reference to know where we started
      auto stop = n;
      do {
        // get the next ptr on A or B
        n = n->prev_on(on);
        // return it if it's the right type
        if (n->vertex_type() == type) return n;
        // continue until we get a nullptr or back to where we started
      } while (n != nullptr && n != stop);
      return nullptr;
    }

//    friend std::ostream & operator<<(std::ostream & os, const Vertex & v);
    friend std::ostream & operator<<(std::ostream & os, const Vertex::ptr & ptr){
      auto is_A = ptr->is_A();
      auto is_B = ptr->is_B();
      auto on = is_A && is_B ? On::both : is_A ? On::A : is_B ? On::B : On::neither;
      os << ptr->value() << ":" << to_string(ptr->vertex_type()) << to_string(on);
      return os;
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
    friend std::ostream & operator<<(std::ostream & os, const VertexList & vertex_list){
      auto h = vertex_list.first();
      os << "VertexList:\non A [ ";
      auto p = h;
      do {
        if (p->is_A()) os << p << ", ";
        p = p->next(On::A);
      } while (p && p != h);
      os << "]\non B [ ";
      p = h;
      do {
        if (p->is_B()) os << p << ", ";
        p = p->next(On::B);
      } while (p && p != h);
      os << "]\n";
      return os;
    }
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

    [[nodiscard]] VertexList a() const {return A;}
    [[nodiscard]] VertexList b() const {return B;}

    [[nodiscard]] std::vector<Wire> intersection_wires() const;
    [[nodiscard]] std::vector<Wire> union_wires() const;
    friend std::ostream & operator<<(std::ostream & os, const VertexLists & vertex_lists){
      os << "VertexLists\n";
      os << "A: " << vertex_lists.a();
      os << "B: " << vertex_lists.b();
      return os;
    }
  };


  template<class T, template<class> class A>
  void insert(On on, const Vertex::ptr & from, std::pair<ind_t, ind_t> edge, const Vertex::ptr & point, A<T> & v){
    auto first = v.view(edge.first);
    auto vec = v.view(edge.second) - first;
    vec /= norm(vec);
    auto dist = dot(v.view(point->value()) - first, vec).sum();
    auto p_dist = 0 * dist;
    auto p = from;
    do {
      p = p->next(on);
      p_dist = dot(v.view(p->value()) - first, vec).sum();
    } while (p_dist < dist && p->value() != edge.second);
    // p is now either the end of the edge (which must be farther from the start than insertion point is)
    // or is the first point on the edge farther than the one to be inserted
    // So insert the provided point between p's predecessor and p:
    auto predecessor = p->prev(on);
    // the inserted point links back to the predecessor, and forward to p
    point->prev(on, predecessor);
    point->next(on, p);
    // p points back to the inserted point
    p->prev(on, point);
    // the predecessor points forward to the inserted point
    predecessor->next(on, point);
  }

  template<class T, template<class> class A>
  A<T> weiler_atherton(A<T> & v, VertexLists & lists){
    // Walk through both lists of vertices, looking for intersections
    auto first_a = lists.first(clip::On::A);
    auto first_b = lists.first(clip::On::B);
    auto ptr_a = first_a;
    do {
      // Find Edge A
      auto edge_a = std::make_pair(ptr_a->value(), ptr_a->next(clip::On::A, clip::Type::original)->value());
      auto ptr_b = first_b;
      do {
        // Find Edge B
        auto edge_b = std::make_pair(ptr_b->value(), ptr_b->next(clip::On::B, clip::Type::original)->value());
        // Find the intersection point of edge A and edge B
        auto [valid, at] = intersection2d(v, edge_a, v, edge_b);
        //
        if (valid > 1){
          // the intersection is of two parallel lines, and two intersection points were returned.
          // they must be one of
          // (a.1st, a.2nd), (a.1st, b.1st), (a.1st, b.2nd), (a.2nd, b.1st), (a.2nd, b.2nd), (b.1st, b.2nd)
          // or the reverse or any pair.
          // We do _not_ want a.1st, so check for any of (a.2nd, b.1st, or b.2nd)
          /* if a.1st is one of the two vertices, we can unambiguously take the other one
           * otherwise, pick the vertex which is closest to a1.
           * */
          auto a1 = v.view(edge_a.first);
          if (at.row_is(cmp::eq, a1).any()){
            at = at.extract(at.row_is(cmp::neq, a1));
          } else {
            auto diff = norm(at - a1); // a (2, 1) Array2
            at = at.extract(diff.row_is(cmp::eq, diff.min(0u)));
          }
        }
        if (valid){
          //  1. Add it to the list of all vertices
          auto index = v.size(0);
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
            // use the other point on the edge of B to decide if this points in our out
            r = v.view(edge_b.first) - a_0;
            cross = cross2d(r, s).val(0,0);
            type = cross > 0 ? clip::Type::exit : cross < 0 ? clip::Type::entry : clip::Type::edge;
          }
          if (clip::Type::edge == type) {
            // use the following vertex of B to decide
            auto follow = ptr_b->next(clip::On::B, clip::Type::original)->next(clip::On::B, clip::Type::original)->value();
            r = v.view(follow) - a_0;
            cross = cross2d(r, s).val(0,0);
            type = cross > 0 ? clip::Type::entry : cross < 0 ? clip::Type::exit : clip::Type::edge;
          }
          // if the type is still 'edge' we can safely(?) skip including this 'intersection' point in the output

          //  3. Insert it into the doubly-linked lists of vertices on both A and B
          auto ptr = std::make_shared<clip::Vertex>(index, type);
          insert(clip::On::A, ptr_a, edge_a, ptr, v);
          insert(clip::On::B, ptr_b, edge_b, ptr, v);
        }
        // move to the next original wire vertex on B
        ptr_b = ptr_b->next(clip::On::B, clip::Type::original);
      } while (ptr_b != first_b);
      // move to the next original wire vertex on A
      ptr_a = ptr_a->next(clip::On::A, clip::Type::original);
    } while (ptr_a != first_a);
    return v;
  }
}

