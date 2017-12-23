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
    bool operator<(const edge& rhs) const
    {
        return (weight > rhs.weight) || !(weight==rhs.weight);
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
                int b = bvalue(b_method, current_node);
                while (T[current_node].size() < b) {
                    // Find the best candidate
                    edge current_best_candidate_edge;
                    current_best_candidate_edge.weight = 0;
                    bool first_iter = true;
                    for (edge candidate_edge : current_graph[current_node]) {
                        if (!(contains(T[current_node], candidate_edge.to))) {
                            if (!S[candidate_edge.to].empty() && (candidate_edge.weight < S[candidate_edge.to].top().weight)) {
                                break;
                            }
                            if (first_iter) {
                                current_best_candidate_edge = candidate_edge;
                                first_iter = false;
                            }
                            if (candidate_edge.weight > current_best_candidate_edge.weight) {
                                current_best_candidate_edge = candidate_edge;
                            }
                        }
                    }
                    if (first_iter) {
                        break;
                    }
                    else { // current_node will adorate current_best_candidate_edge
                        int candidate_node = current_best_candidate_edge.to;
                        if (S[candidate_node].size() > bvalue(b_method, candidate_node)) {
                            // Annuling
                            edge annulled_edge = S[candidate_node].top();
                            std::vector<int>::iterator position = std::find(T[annulled_edge.to].begin(), T[annulled_edge.to].end(), candidate_node);
                            if (position != T[annulled_edge.to].end())
                                T[annulled_edge.to].erase(position);
                            R.push_back(annulled_edge.to);
                        }

                        T[current_node].push_back(current_best_candidate_edge.to);

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
            Q = R;
            R.clear();
        }
        std::vector<int> already_printed;
        for (std::pair<const int, std::priority_queue<edge>> & e : S) {
            if (!(contains(already_printed, e.first))) {
                already_printed.push_back(e.first);
                int counter = 0;
                while(!e.second.empty()) {
                    counter += e.second.top().weight;
                    e.second.pop();
                }
                if (counter != 0)
                std::cout << counter << std::endl;
            }
        }

        // Print out a sum of weights of b-matchings

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

