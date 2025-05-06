#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <tuple>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <random>  // Required for std::mt19937 and std::shuffle

using namespace std;

// Random integer in [low, high]
int rand_int(int low, int high) {
    return low + rand() % (high - low + 1);
}

int main() {
    srand(static_cast<unsigned>(time(nullptr)));

    const int num_nodes = 10; // Must be even to keep edges integer with avg degree 3
    const int num_edges = (3 * num_nodes) / 2;
    const int num_deletions = 2;
    const int num_insertions = 2;
    const int max_weight = 10;

    set<tuple<int, int>> edge_set;
    vector<tuple<int, int, int>> edges;

    // Generate unique undirected edges
    while (edges.size() < static_cast<size_t>(num_edges)) {
        int u = rand_int(0, num_nodes - 1);
        int v = rand_int(0, num_nodes - 1);
        if (u == v) continue;
        if (u > v) std::swap(u, v);

        if (edge_set.emplace(u, v).second) {
            int w = rand_int(1, max_weight);
            edges.emplace_back(u, v, w);
        }
    }

    // Write graph.txt
    ofstream graph_out("graph.txt");
    for (const auto& [u, v, w] : edges) {
        graph_out << u << " " << v << " " << w << "\n";
    }
    graph_out.close();
    cout << "Created graph.txt with " << edges.size() << " edges and average degree 3\n";

    // Prepare deletions
    // Prepare deletions
    vector<tuple<int, int, int>> deletions;
    std::mt19937 rng(static_cast<unsigned>(time(nullptr)));
    std::shuffle(edges.begin(), edges.end(), rng);

    for (int i = 0; i < num_deletions && i < static_cast<int>(edges.size()); ++i) {
        deletions.push_back(edges[i]);
        edge_set.erase({min(get<0>(edges[i]), get<1>(edges[i])), max(get<0>(edges[i]), get<1>(edges[i]))});
    }

    // Prepare insertions
    vector<tuple<int, int, int>> insertions;
    while (insertions.size() < static_cast<size_t>(num_insertions)) {
        int u = rand_int(0, num_nodes - 1);
        int v = rand_int(0, num_nodes - 1);
        if (u == v) continue;
        if (u > v) std::swap(u, v);
        if (edge_set.count({u, v}) == 0) {
            int w = rand_int(1, max_weight);
            insertions.emplace_back(u, v, w);
            edge_set.emplace(u, v);
        }
    }

    // Write changes.txt
    ofstream changes_out("changes.txt");
    changes_out << deletions.size() << "\n";
    for (const auto& [u, v, w] : deletions) {
        changes_out << u << " " << v << " " << w << "\n";
    }
    changes_out << insertions.size() << "\n";
    for (const auto& [u, v, w] : insertions) {
        changes_out << u << " " << v << " " << w << "\n";
    }
    changes_out.close();
    cout << "Created changes.txt with " << deletions.size() << " deletions and " << insertions.size() << " insertions\n";

    return 0;
}
