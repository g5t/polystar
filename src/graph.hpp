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
    std::pair<std::vector<cost_t>, std::vector<size_t>> find_paths(size_t, cost_t);
    std::pair<cost_t, std::vector<size_t>> reverse_path(const std::vector<cost_t>&, const std::vector<size_t>&, size_t, size_t);
    std::pair<cost_t, std::vector<size_t>> forward_path(const std::vector<cost_t>&, const std::vector<size_t>&, size_t, size_t);
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
  std::pair<std::vector<cost_t>, std::vector<size_t>> Graph<cost_t>::find_paths(size_t src, cost_t max_cost){
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
      return std::make_pair(dist, prev);
  }
  template<class cost_t>
  std::pair<cost_t, std::vector<size_t>> Graph<cost_t>::reverse_path(const std::vector<cost_t>& dist, const std::vector<size_t>& prev, size_t last, size_t first){
      std::vector<size_t> reversed;
      reversed.reserve(no_);
      reversed.push_back(last);
      while (reversed.back() != first) reversed.push_back(prev[reversed.back()]);
      return std::make_pair(dist[last], reversed);
  }
  template<class cost_t>
  std::pair<cost_t, std::vector<size_t>> Graph<cost_t>::forward_path(const std::vector<cost_t>& dist, const std::vector<size_t>& prev, size_t last, size_t first){
      auto [cost, reversed] = reverse_path(dist, prev, last, first);
      std::vector<size_t> path;
      path.reserve(reversed.size());
      std::copy(reversed.rbegin(), reversed.rend(), std::back_inserter(path));
      return std::make_pair(cost, path);
  }
  template<class cost_t>
  std::vector<size_t> Graph<cost_t>::shortestPath(size_t src, size_t snk, cost_t max_cost) {
    auto [dist, prev] = find_paths(src, max_cost);
    auto [cost, path] = forward_path(dist, prev, snk, src);
    return path;
  }
}

#endif