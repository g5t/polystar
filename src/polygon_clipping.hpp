#pragma once
#include "polygon_poly.hpp"
#include "polygon_wire.hpp"

namespace polystar::polygon::clip
{
  using ind_t = polystar::ind_t;

  enum class Type {unknown, entry, exit, original};
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
    [[nodiscard]] ptr next(On on) const {
      if (on == On::A) return next_A;
      if (on == On::B) return next_B;
      return nullptr;
    }
    [[nodiscard]] ptr prev(On on) const {
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

    [[nodiscard]] std::vector<Wire> intersection_wires() const;
    [[nodiscard]] std::vector<Wire> union_wires() const;
  };

  // TODO Implement the line segment intersection finder to setup the VertexLists object for two Poly objects
}