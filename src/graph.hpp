#ifndef POLYSTAR_GRAPH_HPP
#define POLYSTAR_GRAPH_HPP
#include <memory>
#include <vector>
#include <algorithm>
#include <exception>
namespace polystar::graph {
    template<class cost_t> class Node;
    template<class cost_t>
    class Node {
    public:
        using node_t = std::shared_ptr<Node<cost_t>>;
        using prev_t = std::weak_ptr<Node<cost_t>>;
        using tree_t = std::vector<std::pair<cost_t, prev_t>>;
    private:
        node_t node_;
        tree_t tree_;
        prev_t prev_;
        cost_t cost_=cost_t(0);
        bool visited_=false;
    public:
        void insert(node_t n, cost_t c) {
            auto ptr = std::find_if(tree_.begin(), tree_.end(), [&](const auto & wt){
                if (auto t = wt.second.lock()) return t == n;
                return false;
            });
            if (ptr != tree_.end()) {
                ptr->first = c;
                ptr->second = n;
            } else {
                tree_.emplace_back(c, n);
            }
        }
        node_t node() const {return node_;}
        cost_t cost() const {return cost_;}
        cost_t set_cost(cost_t c) {cost_ = c; return cost_;}
        bool visited() const {return visited_;}
        void visit() {visited_=true;}
        tree_t tree() const {return tree_;}
        prev_t points_at() const {return prev_;}
        void point_at(node_t n) {prev_ = n;}

    };

    template<class cost_t>
    class Graph {
    public:
        using node_t = typename Node<cost_t>::node_t;
        using tree_t = typename Node<cost_t>::tree_t;
        using nodes_t = std::vector<node_t>;
    private:
        nodes_t nodes_;
    public:
        template<class T, class Connected, class Cost>
        Graph(const std::vector<T> & items, Connected connected, Cost cost){
            // create a node shared pointer for every input item
            nodes_.clear();
            nodes_.reserve(items.size());
            for (const auto & x: items) nodes_.push_back(std::make_shared<node_t>());
            // connected(item[i], item[j]) returns true if there is an edge connecting the two items
            // cost(item[i], item[j]) returns the pairwise cost associated with the edge
            for (size_t i=0; i<items.size()-1; ++i){
                for (size_t j=i+1; j<items.size(); ++j){
                    if (connected(items[i], items[j])){
                        auto c = cost(items[i], items[j]);
                        nodes_[i].insert(nodes_[j], c);
                        nodes_[j].insert(nodes_[i], c);
                    }
                }
            }
        }
        nodes_t nodes() const {return nodes_;}

        // using this will require keeping an external bidirectional mapping from item to shared_ptr<Node<T>>
        nodes_t dijkstra(node_t from, node_t to) const {
            from->visit();
            nodes_t border{from};
            while (!border.empty()){
                // find the lowest-cost border node
                std::partial_sort(border.begin(), border.begin()+1, border.end(),
                                  [](const auto & a, const auto &b){return a->cost() < b->cost();});
                auto cheapest = border.front();
                // it can not get cheaper, so remove it from the border
                border.erase(border.begin());
                // and add it to the visited nodes
                cheapest->visit();

                // get the list of edge connected nodes
                auto next = cheapest->tree();
                // for every (cost, weak_ptr) that still points at a valid shared_ptr:
                for (const auto & cw: next) if (auto n=cw.second.lock()) {
                    // the cost to get to the next node from the current cheapest
                    auto cost = cheapest->cost() + cw.first;
                    // if we have no idea of its cost, or think it should cost more
                    if (!n->visited() || n->cost > cost) {
                        // update its cost
                        n->set_cost(cost);
                        // point it back towards the current cheapest node
                        n->point_at(cheapest);
                        // and add it to the border stack
                        if (std::find(border.begin(), border.end(), n) == border.end()) border.push_back(n);
                    }
                }
            }
            // the solution exists. Now find the path by walking the points_at weak_ptr references backwards:
            nodes_t rev_path{to};
            while (rev_path.back() != from){
                rev_path.push_back(rev_path.back()->points_at().lock());
            }
            nodes_t path;
            path.reserve(rev_path.size());
            std::copy(rev_path.rbegin(), rev_path.rend(), std::back_inserter(path));
            return path;
        }

//        tree_t & operator[](node_t node) {
//            auto ptr = std::find_if(nodes_.begin(), nodes_.end(), [&](const auto & n){return n.node() == node;});
//            if (ptr != nodes_.end()) return *ptr;
//            throw std::runtime_error("Node not in graph");
//        }
//        const tree_t & operator[](node_t node) const {
//            auto ptr = std::find_if(nodes_.cbegin(), nodes_.cend(), [&](const auto & n){return n.node() == node;});
//            if (ptr != nodes_.cend()) return *ptr;
//            throw std::runtime_error("Node not in graph");
//        }

//        // Allow for insertion of externally defined nodes -- but their connections must already be known
//        // ! This doesn't include the reverse of each edge ! Do not allow for now.
//        void insert(node_t node, tree_t tree){
//            auto ptr = std::find_if(nodes_.begin(), nodes_.end(), [&](const auto & n){return n.node() == node;});
//            if (ptr != nodes_.end()){
//                *ptr = tree;
//            } else {
//                nodes_.push_back(tree);
//            }
//        }
//        void remove(node_t node){
//            auto ptr = std::find_if(nodes_.begin(), nodes_.end(), [&](const auto & n){return n.node() == node;});
//            if (ptr != nodes_.end()) nodes_.erase(ptr);
//        }

    };
}

#endif