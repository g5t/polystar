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
    void vertex_type(Type type) {_type = type;}
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
    auto insert_point = [&v, &lists](const auto & point, const auto & eA, const auto & eB, auto & pA, auto & pB){
      debug_update("which is", point.to_string(0));
      //  1. Add it to the list of all vertices
      auto index = v.size(0);
      auto match = v.row_is(cmp::eq, point);
      std::shared_ptr<clip::Vertex> match_ptr{nullptr};
      if (match.any()){
        // the intersection point is already in the list of vertices
        index = match.first();
        // see if this point was already added to lists as entry or exit type, and shortcut out if it has
        // we can walk either wire
        auto w = lists.first(On::A);
        auto w_first = w;
        do {
          if (w->value() == index){
            auto t = w->vertex_type();
            if (clip::Type::exit == t || clip::Type::entry == t) return;
            if (clip::Type::edge == t) match_ptr = w;
          }
          w = w->next(On::A);
        } while (w != w_first && w != nullptr);
      } else {
        // the intersection point is not in the list of vertices
        index = v.size(0);
        v = v.append(0, point);
      }
      //  2. Find whether it points into or out of A.
      auto a_0 = v.view(eA.first);
      auto r = v.view(eB.second) - a_0; // we want to see if edge_b.second is to the right of edge_a
      auto s = v.view(eA.second) - a_0;
      auto cross = cross2d(r, s).val(0,0);
      auto type = cross > 0 ? clip::Type::entry : cross < 0 ? clip::Type::exit : clip::Type::edge;
      if (match_ptr){
        if (clip::Type::edge != type) match_ptr->vertex_type(type);
        debug_update("Updated ", match_ptr);
//        std::cout << "Updated " << match_ptr << "\n";
        return;
      }

      //  3. Insert it into the doubly-linked lists of vertices on both A and B
      auto ptr = std::make_shared<clip::Vertex>(index, type);
      insert(clip::On::A, pA, eA, ptr, v);
      insert(clip::On::B, pB, eB, ptr, v);
      debug_update("Inserted ", ptr);
//      std::cout << "Inserted " << ptr << "\n";
    };

    // Walk through both lists of vertices, looking for intersections
    auto first_a = lists.first(clip::On::A);
    auto first_b = lists.first(clip::On::B);
    auto ptr_a = first_a;
    debug_update("vertices ", v.to_string());
//    std::cout << "vertices " << v.to_string() << "\n";
    do {
      // Find Edge A
      auto edge_a = std::make_pair(ptr_a->value(), ptr_a->next(clip::On::A, clip::Type::original)->value());
      auto ptr_b = first_b;
      do {
        // Find Edge B
        auto edge_b = std::make_pair(ptr_b->value(), ptr_b->next(clip::On::B, clip::Type::original)->value());
        // Find the intersection point of edge A and edge B
        debug_update("Look for intersection of edge ", edge_a, " and ", edge_b);
//        std::cout << "Look for intersection of edge (";
//        std::cout << edge_a.first << "," << edge_a.second << ") and (";
//        std::cout << edge_b.first << "," << edge_b.second << ");";
        auto [valid, at] = intersection2d(v, edge_a, v, edge_b);
        debug_update("Found ", valid, " intersection", (valid != 1 ? "s" : ""));
//        std::cout << " found " << valid << " intersection" << (valid != 1 ? "s": "") << "\n";

        if (valid) {
          for (ind_t i=0; i<at.size(0); ++i){
            insert_point(at.view(i), edge_a, edge_b, ptr_a, ptr_b);
          }
        }
        // move to the next original wire vertex on B
        ptr_b = ptr_b->next(clip::On::B, clip::Type::original);
      } while (ptr_b != first_b);
      // move to the next original wire vertex on A
      ptr_a = ptr_a->next(clip::On::A, clip::Type::original);
    } while (ptr_a != first_a);

    // try to handle any edge vertices
    int excess{0}, edge_vertices{0}, total_vertices{0};
    auto p = lists.first(On::A);
    auto p_first = p;
    do{
      switch (p->vertex_type()){
        case clip::Type::entry: ++excess; break;
        case clip::Type::exit: --excess; break;
        case clip::Type::edge: ++edge_vertices; break;
        default: ++total_vertices;
      }
      p = p->next(On::A);
    } while (p != p_first && p != nullptr);

    while (excess && edge_vertices > 0){
      p = p_first;
      do{
        if (p->vertex_type() == clip::Type::edge) {
          p->vertex_type(excess > 0 ? clip::Type::exit : clip::Type::entry);
          --edge_vertices;
          excess > 0 ? --excess : ++excess;
        }
        p = p->next(On::A);
      } while (p != p_first && p != nullptr && excess && edge_vertices > 0);
    }

    return v;
  }

  enum class Inside {neither, A, B};

  /** \brief Implement the edge chasing algorithm from Computer Graphics and Image Processing 19, 384 (1982)
   *
   * An algorithm is outlined in 'A New Linear Algorithm for Intersecting Convex Polygons' which can find the
   * intersection wire of two convex polygons (with N and M points, respectively) in no worse than 2(N+M) steps.
   * It relies on finding pairwise edge intersections between the two polygons, progressing through edge pairs via
   * a heuristic that favors the edge that points _towards_ the other polygon.
   */
  template<class T, template<class> class A>
  std::pair<A<T>, Wire>
  orourke_chien_olsen_naddor(const A<T> & va, const Wire & wa, const A<T> & vb, const Wire & wb){
    auto N = wa.size();
    auto M = wb.size();
    A<T> found(0, 2); // no more than (N-1)*(M-1) edge intersections ... but there is no constructor which allocates this
    Wire wire;
    size_t maximum = 2 * (N+M);
    wire.reserve(maximum);

    auto half_plane = [](const auto & v, const auto & e, const auto & x){
      auto p_minus = v.view(e.first);
      auto p = v.view(e.second) - p_minus;
      auto y = x - p_minus;
      return cross2d(p, y).val(0, 0) >= 0;
    };
    auto add_found = [&found](const auto & x){
      auto match = found.row_is(cmp::eq, x.view(0));
      if (match.any()){
        return match.first();
      }
      found.append(0, x.view(0));
//      std::cout << found.to_string() << "\n";
      return found.size(0u) - 1;
    };
    auto add_wire = [&wire](const auto index){
      if (wire.empty() || wire.back() != index) wire.push_back(index);
    };
    auto add = [&add_wire, &add_found](const auto & x){
      add_wire(add_found(x));
    };

    size_t p{0}, q{0};
    Inside inside{Inside::neither};
    do {
      auto ea = wa.edge(p);
      auto eb = wb.edge(q);
//      std::cout << "p=" << p << " q=" << q << " ";
//      std::cout << "Edges (" << ea.first << ", " << ea.second << ") and (" << eb.first << ", " << eb.second << ")\n";
      auto [intersections, at] = intersection2d(va, ea, vb, eb, end_type::both, end_type::both);
      if (intersections == 1){
        // not collinear and intersecting:
        auto index = add_found(at.view(0));
//        std::cout << intersections << " Intersection of (" << ea.first << ", " << ea.second << ") and (";
//        std::cout << eb.first << ", " << eb.second << ") at " << at.to_string(0) << " has index " << index << "\n";
        if (!wire.empty() && 0 != wire.back() && 0 == index){
          // we're done -- return the vertices found up to this point
          return {found, wire};
        }
        // append this vertex, if we haven't repeated a point:
        add_wire(index);
        // do half-plane check
        inside = half_plane(vb, eb, va.view(ea.second)) ? Inside::A : Inside::B; // edge_a ends in edge_b's half plane?
//        std::cout << "Such that " << (Inside::A == inside ? "A" : "B") << " is inside now\n";
      }
      // if B x A >= 0, advance(B if A in half plane of B else A), else advance(A if B in half-plane of A else B)
      auto q_vec = vb.view(eb.second) - vb.view(eb.first);
      auto p_vec = va.view(ea.second) - va.view(ea.first);
      auto c = cross2d(q_vec, p_vec).val(0, 0);
      auto a_in_b = half_plane(vb, eb, va.view(ea.second));
      auto b_in_a = half_plane(va, ea, vb.view(eb.second));
//      std::cout << "Edge cross product " << c << " and half plane inclusions (a in b) = " << a_in_b << " and (b in a) = " << b_in_a << "\n";
      auto advance_b = c < 0 ? !b_in_a : a_in_b;
      if (advance_b){
//        std::cout << "Advance b";
        if (Inside::B == inside) {
//          std::cout << " inserting " << eb.second;
          add(vb.view(eb.second));
        }
        ++q;
//        std::cout << " such that q is " << q << "\n";
      } else { // advance_a
//        std::cout << "Advance a";
        if (Inside::A == inside) {
//          std::cout << " inserting " << ea.second;
          add(va.view(ea.second));
        }
        ++p;
//        std::cout << " such that p is " << p << "\n";
      }
    } while (p + q < maximum);
    if (wire.empty()){
      // no intersections -- return the _inner_ polygon (if there is one)
      if (wa.contains(vb.view(wb.front()), va).front()) return {vb, wb};
      if (wb.contains(va.view(wa.front()), vb).front()) return {va, wa};
    }
    return {found, wire};
  }
}

