#include "blimit.hpp"

#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <sstream>
#include <string>
#include <queue>
#include <algorithm>
#include <set>

typedef std::map<int, std::set<int>> graph;
void add_edge(graph& current_graph, int node1, int node2);
void add_weight(std::map<std::pair<int, int>, int>& weights, int node1, int node2, int weight);
void construct_graph_from_file(std::ifstream& infile,
                               graph& current_graph,
                               std::set<int>& nodes,
                               std::map<std::pair<int, int>, int>& weights);
class NodeComparator
{
    std::map<std::pair<int, int>, int> weights;
    int parent_node;

public:

    NodeComparator(const int& t_parent_node, const std::map<std::pair<int, int>, int>& t_weights) {
        weights = t_weights;
        parent_node = t_parent_node;
    }

    bool operator() (const int& lhs, const int& rhs) {
        return weights[std::make_pair(parent_node, lhs)] > weights[std::make_pair(parent_node, rhs)];
    }
};




int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "usage: " << argv[0] << " thread-count inputfile b-limit" << std::endl;
        return 1;
    }

    int thread_count = std::stoi(argv[1]);
    int b_limit = std::stoi(argv[3]);
    std::string input_filename{argv[2]};
    std::ifstream infile(input_filename);
    graph current_graph;
    std::set<int> nodes;
    std::map<std::pair<int, int>, int> weights;
    construct_graph_from_file(infile, current_graph, nodes, weights);

    for (int b_method = 0; b_method < b_limit + 1; b_method++) {
        std::vector<int> Q(nodes.begin(), nodes.end());
        std::vector<int> R;
        std::map<int, std::vector<int>> T;
        std::map<int, std::priority_queue<int, std::vector<int>, NodeComparator>> S;
        for (int node : nodes) {
            std::priority_queue<int, std::vector<int>, NodeComparator> q(NodeComparator(node, weights));
            S[node] = q;
            T[node] = std::vector<int>();
        }

        while (!Q.empty()) {
            for (int current_node : Q) {
                int b = bvalue(b_method, current_node);
                while (T[current_node].size() < b) {
                    // Find the best candidate
                    int current_best_candidate = -1;
                    for (int candidate_node : current_graph[current_node]) {
//                        if (!(contains(T[current_node], candidate_edge.to))) {
                        // Get candidate_node's least attractive suitor



//                        }
//                            if ((S[candidate_node.to].empty()) ||
//                                    (weights[std::make_pair(candidate_node, current_node)] < S[candidate_node].top().weight)) {
//                                break;
//                            }
//                            if (first_iter) {
//                                current_best_candidate_edge = candidate_edge;
//                                first_iter = false;
//                            }
//                            if (candidate_edge.weight > current_best_candidate_edge.weight) {
//                                current_best_candidate_edge = candidate_edge;
//                            }
//                        }
//                    }
//                    if (first_iter) {
//                        break;
//                    } else { // current_node will adorate current_best_candidate_edge
//                        int candidate_node = current_best_candidate_edge.to;
//                        if (S[candidate_node].size() > bvalue(b_method, candidate_node)) {
//                            // Annuling
//                            edge annulled_edge = S[candidate_node].top();
//                            std::vector<int>::iterator position = std::find(T[annulled_edge.to].begin(),
//                                                                            T[annulled_edge.to].end(), candidate_node);
//                            if (position != T[annulled_edge.to].end())
//                                T[annulled_edge.to].erase(position);
//                            R.push_back(annulled_edge.to);
//                        }
//
//                        T[current_node].push_back(current_best_candidate_edge.to);
//
//                        // Find corresponding edge in the opposite direction
//                        std::vector<edge> edges_of_candidate = current_graph[candidate_node];
//                        for (edge e : edges_of_candidate) {
//                            if (e.to == current_node) {
//                                S[candidate_node].push(e);
//                                break;
//                            }
//                        }

//                    }
//                }
//            }
//            Q = R;
//            R.clear();
//        }
//        std::vector<int> already_printed;
//        for (std::pair<const int, std::priority_queue<edge>> &e : S) {
//            if (!(contains(already_printed, e.first))) {
//                already_printed.push_back(e.first);
//                int counter = 0;
//                while (!e.second.empty()) {
//                    counter += e.second.top().weight;
//                    e.second.pop();
//                }
//                if (counter != 0)
//                    std::cout << counter << std::endl;
//            }
//        }

                        // Print out a sum of weights of b-matchings

                    }
                }
            }
        }
    }
}

void construct_graph_from_file(const std::ifstream& infile,
                               graph& current_graph,
                               std::set<int>& nodes,
                               std::map<std::pair<int, int>, int>& weights) {
    int node1, node2, weight;
    std::string line;
    while (getline(infile, line)) {
        if (line[0] == '#') { continue; }
        std::istringstream iss(line);
        if (!(iss >> node1 >> node2 >> weight)) {
            std::cout << "ERROR" << std::endl;
            break;
        }
        nodes.insert(node1);
        nodes.insert(node2);
        add_weight(weights, node1, node2, weight);
        add_weight(weights, node2, node1, weight);
        add_edge(current_graph, node1, node2);
        add_edge(current_graph, node2, node1);

    }
}

void add_weight(std::map<std::pair<int, int>, int>& weights,
                const int node1,
                const int node2,
                const int weight) {
    weights.insert(std::make_pair(std::make_pair(node1, node2), weight));
}


void add_edge(graph& current_graph, const int node1, const int node2) {
    graph::iterator it = current_graph.find(node1);
    if (it == current_graph.end()) { // node1 doesnt exist in graph
        std::set<int> neighbors;
        neighbors.insert(node2);
        current_graph[node1] = neighbors;
    } else { // node1 exists
        current_graph[node1].insert(node2);
    }
}

//bool queue_order(const std::map<std::pair<int, int>, int>& weights,
//                 const std::pair<int, int> pair1,
//                 const std::pair<int, int> pair2) {
//    return weights[pair1] > weights[pair2];
//}
template<class C, class T>
auto contains(const C& v, const T& x)
-> decltype(end(v), true)
{
    return end(v) != std::find(begin(v), end(v), x);
}