#ifndef POLYSTAR_GRAPH_PQ_HPP
#define POLYSTAR_GRAPH_PQ_HPP

#include <queue>
#include <vector>
#include <list>

namespace polystar::graph {
  template<class cost_t>
  class Graph {
  public:
    using next_t = std::pair<cost_t, size_t>;
  private:
    size_t no_;
    std::vector<std::list<next_t>> adj_;
  public:
    Graph(size_t no);
    void add_bidirectional(size_t u, size_t v, cost_t weight);
    void add_directional(size_t u, size_t v, cost_t weight);
    std::vector<size_t> shortestPath(size_t s, size_t f, cost_t m);
  };

  template<class cost_t>
  Graph<cost_t>::Graph(size_t no): no_(no) {
    adj_.clear();
    adj_.resize(no_);
  }
  template<class cost_t>
  void Graph<cost_t>::add_bidirectional(size_t u, size_t v, cost_t w) {
    adj_[u].emplace_back(w, v);
    adj_[v].emplace_back(w, u);
  }
  template<class cost_t>
  void Graph<cost_t>::add_directional(size_t u, size_t v, cost_t w) {
    adj_[u].emplace_back(w, v);
  }
  template<class cost_t>
  std::vector<size_t> Graph<cost_t>::shortestPath(size_t src, size_t snk, cost_t max_cost) {
    std::priority_queue<next_t, std::vector<next_t>, std::greater<next_t>> pq;
    std::vector<cost_t> dist(no_, max_cost);
    // Insert source itself in priority queue and initialize its distance as 0.
    pq.push(std::make_pair(0., src));
    dist[src] = 0;
    // keep track of which node was the previous node before each
    std::vector<size_t> prev(no_, no_+1);

    while (!pq.empty()) {
      auto u = pq.top().second;
      pq.pop();
      // loop over all adjacent vertices of a vertex
      for (auto i = adj_[u].begin(); i != adj_[u].end(); ++i) {
        auto weight = (*i).first;
        auto v = (*i).second;
        // If there is shorter path to v through u.
        if (dist[v] > dist[u] + weight) {
          dist[v] = dist[u] + weight;
          pq.push(std::make_pair(dist[v], v));
          prev[v] = u;
        }
      }
    }
    // follow the path backwards
    std::vector<size_t> rpath;
    rpath.reserve(no_);
    rpath.push_back(snk);
    while (rpath.back() != src) rpath.push_back(prev[rpath.back()]);
    // and reverse the path to the direction we want
    std::vector<size_t> path;
    path.reserve(rpath.size());
    std::copy(rpath.rbegin(), rpath.rend(), std::back_inserter(path));
    return path;
  }

}

#endif