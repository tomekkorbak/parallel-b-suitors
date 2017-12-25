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


struct edge {
    int to;
    int weight;
    bool operator < (const edge& rhs) const {
        return (weight >= rhs.weight);
    }
};

typedef std::map<int, std::vector<edge>> graph;

void add_edge(int node1, int node2, int weight, graph &graph);
void construct_graph_from_file(std::ifstream &infile, graph &graph, std::set<int> &nodes);

template<class C, class T>
auto contains(const C& v, const T& x)
-> decltype(end(v), true)
{
    return end(v) != std::find(begin(v), end(v), x);
}

int main(int argc, char* argv[]) {
    bool debug = true;


    if (argc != 4) {
        std::cerr << "usage: "<<argv[0]<<" thread-count inputfile b-limit"<< std::endl;
        return 1;
    }

    int thread_count = std::stoi(argv[1]);
    int b_limit = std::stoi(argv[3]);
    std::string input_filename{argv[2]};
    std::ifstream infile(input_filename);
    graph current_graph;
    std::set<int> nodes;
    construct_graph_from_file(infile, current_graph, nodes);

    for (int b_method = 0; b_method < b_limit + 1; b_method++) {
        if (debug) {
            std::cerr << "Start processing for b_method " << b_method << std::endl;
        }
        // this is just to show the blimit with which the program is linked
//        std::cerr << "bvalue node 44: " << bvalue(b_method, 44) << std::endl;

        // TODO: implement b-adorators here
        std::map<int, std::vector<int>> T;
        std::vector<int> R;
        std::map<int, std::priority_queue<edge>> S;

        for(int node : nodes) {
            S[node] = std::priority_queue<edge>();
            T[node] = std::vector<int>();
        }
        std::vector<int> Q(nodes.begin(), nodes.end());
        while (!Q.empty()) {
            for (int current_node : Q) {
                if (debug) {
                    std::cerr << " Looking for mates for " << current_node << std::endl;
                }
                int b = bvalue(b_method, current_node);
                while (T[current_node].size() < b) {
                    // Find the best candidate
                    bool found_candidate = false;
                    int current_best_weight = 0;
                    edge current_best_candidate_edge;
                    for (edge candidate_edge : current_graph[current_node]) {
                        if (!(contains(T[current_node], candidate_edge.to))) {
                            if ((!S[candidate_edge.to].empty() && (candidate_edge.weight > S[candidate_edge.to].top().weight))
                                || (!S[candidate_edge.to].empty() && S[candidate_edge.to].size() < bvalue(b_method, candidate_edge.to))
                                   || S[candidate_edge.to].empty()) {
                                found_candidate = true;
                                if (candidate_edge.weight > current_best_weight) {
                                    current_best_candidate_edge = candidate_edge;
                                    current_best_weight = candidate_edge.weight;
                                    if (debug) {
                                        std::cerr << "   " << candidate_edge.to << " with weight " << candidate_edge.weight
                                                  << " is current best candidate for " << current_node << ", whose b = "
                                                  << b << std::endl;
                                    }
                                }
                                else {
                                    if (debug) {
                                        std::cerr << "   " << candidate_edge.to << " with weight " << candidate_edge.weight
                                                  << " is not a suitable mate for " << current_node << ", whose b = "
                                                  << b << std::endl;
                                    }
                                }
                            }
                        }
                    }
                    if (!(found_candidate)) {
                        break;
                    }
                    else { // current_node will adorate current_best_candidate_edge
                        int candidate_node = current_best_candidate_edge.to;
                        if (debug) {
                            std::cerr << "  " << candidate_node << " is being married with " << current_node << std::endl;
                        }

                        if (S[candidate_node].size() == bvalue(b_method, candidate_node)) {
                            // Annuling
                            edge annulled_edge = S[candidate_node].top();
                            S[candidate_node].pop();
                            if (debug) {
                                std::cerr << "   Annuling marriage between " << candidate_node << " and "
                                          << annulled_edge.to << " with weight " << annulled_edge.weight << ", beacause "
                                          << current_node << " is a better party with weight "
                                          << current_best_candidate_edge.weight << std::endl;
                            }
                            std::vector<int>::iterator position = std::find(T[annulled_edge.to].begin(),
                                                                            T[annulled_edge.to].end(),
                                                                            candidate_node);
                            if (position != T[annulled_edge.to].end())
                                T[annulled_edge.to].erase(position);
                            R.push_back(annulled_edge.to);
                        }
                        T[current_node].push_back(candidate_node);

                        // Find corresponding edge in the opposite direction
                        std::vector<edge> edges_of_candidate = current_graph[candidate_node];
                        for (edge e : edges_of_candidate) {
                            if (e.to == current_node) {
                                S[candidate_node].push(e);
                                break;
                            }
                        }
                    }
                }
            }
            if (debug) {
                for (auto a : T) {
                    std::cerr << " T[" << a.first << "] = ";
                    for (auto b : a.second) {
                        std::cerr << b << " ";
                    }
                    std::cerr << std::endl;
                }
                std::map<int, std::priority_queue<edge>> temp = S;
                for (auto a : temp) {
                    std::cerr << " S[" << a.first << "] = ";
                    while (!a.second.empty()) {
                        std::cerr << a.second.top().to << " ";
                        a.second.pop();
                    }
                    std::cerr << std::endl;
                }

            }
            Q = R;
            R.clear();
        }

        if (debug) {
            std::cerr << "Done processing for b_method " << b_method << std::endl;
        }

        std::set<int> already_printed;
        int counter = 0;
        for (std::pair<const int, std::priority_queue<edge>> & e : S) {
            if (!(contains(already_printed, e.first))) {
                already_printed.insert(e.first);

                while(!e.second.empty()) {
                    edge top = e.second.top();
                    if (!(contains(already_printed, top.to))) {
                        counter += top.weight;
                    }

                    e.second.pop();
                }
            }
        }
        std::cout << counter << std::endl;


    }
}

void construct_graph_from_file(std::ifstream& infile, graph& current_graph, std::set<int>& nodes) {
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
        add_edge(node1, node2, weight, current_graph);
        add_edge(node2, node1, weight, current_graph);

    }
}

void add_edge(int node1, int node2, int weight, graph& graph) {
    edge current_edge{node2, weight};
    graph::iterator it = graph.find(node1);
    if (it == graph.end()) { // node1 doesnt exist in graph
            std::vector<edge> adjList;
            adjList.push_back(current_edge);
            graph[node1] = adjList;
        } else { // node1 exists
            graph[node1].push_back(current_edge);
        }
}

