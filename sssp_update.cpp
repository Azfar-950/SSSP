#include <mpi.h>
#include <omp.h>
#include <metis.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <set>
#include <climits>
#include <string>
#include <sstream>

#define INF INT_MAX

struct Edge {
    int u, v, weight;
    char op; // 'I' for insert, 'D' for delete
};

struct Vertex {
    int parent;
    int dist;
    bool affected;
    bool affected_del;
    std::vector<std::pair<int, int>> adj; // (neighbor, weight)
};

// Read graph in METIS format
void read_graph(const std::string& filename, int& n, int& m, std::vector<Vertex>& graph) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error opening graph file" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    file >> n >> m;
    graph.resize(n);
    std::string line;
    std::getline(file, line); // Skip newline
    for (int i = 0; i < n; ++i) {
        std::getline(file, line);
        std::istringstream iss(line);
        int v, w;
        while (iss >> v >> w) {
            graph[i].adj.push_back({v - 1, w}); // 0-based indexing
        }
    }
    file.close();
}

// Read changes file
void read_changes(const std::string& filename, std::vector<Edge>& changes) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error opening changes file" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        char op;
        int u, v, w;
        iss >> op >> u >> v >> w;
        changes.push_back({u - 1, v - 1, w, op}); // 0-based indexing
    }
    file.close();
}

// Partition graph using METIS
void partition_graph(int n, std::vector<Vertex>& graph, int nparts, std::vector<int>& part) {
    idx_t nvtxs = n;
    idx_t ncon = 1;
    std::vector<idx_t> xadj(n + 1), adjncy, adjwgt;
    xadj[0] = 0;
    for (int i = 0; i < n; ++i) {
        for (const auto& [v, w] : graph[i].adj) {
            adjncy.push_back(v);
            adjwgt.push_back(w);
        }
        xadj[i + 1] = adjncy.size();
    }
    idx_t objval;
    part.resize(n);
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
    int ret = METIS_PartGraphKway(&nvtxs, &ncon, xadj.data(), adjncy.data(), nullptr, nullptr,
                                  adjwgt.data(), &nparts, nullptr, nullptr, options, &objval, part.data());
    if (ret != METIS_OK) {
        std::cerr << "METIS partitioning failed" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

// Identify ghost nodes
void identify_ghost_nodes(const std::vector<Vertex>& graph, const std::vector<int>& part, int rank, int nparts,
                          std::set<int>& ghost_nodes, std::vector<std::vector<int>>& local_vertices) {
    local_vertices.resize(nparts);
    for (int i = 0; i < graph.size(); ++i) {
        local_vertices[part[i]].push_back(i);
        if (part[i] == rank) {
            for (const auto& [v, _] : graph[i].adj) {
                if (part[v] != rank) {
                    ghost_nodes.insert(v);
                }
            }
        }
    }
}

// Process edge changes (Step 1)
void process_changes(std::vector<Vertex>& graph, const std::vector<Edge>& changes, int rank,
                     const std::vector<int>& part, std::vector<bool>& affected, std::vector<bool>& affected_del) {
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < changes.size(); ++i) {
        const auto& e = changes[i];
        if (part[e.u] == rank || part[e.v] == rank) {
            int x, y;
            if (graph[e.u].dist > graph[e.v].dist) {
                x = e.v;
                y = e.u;
            } else {
                x = e.u;
                y = e.v;
            }
            if (e.op == 'D') {
                bool is_tree_edge = false;
                for (const auto& [v, w] : graph[x].adj) {
                    if (v == y && w == e.weight && graph[x].parent == y || graph[y].parent == x) {
                        is_tree_edge = true;
                        break;
                    }
                }
                if (is_tree_edge) {
                    #pragma omp critical
                    {
                        graph[y].dist = INF;
                        graph[y].parent = -1;
                        affected_del[y] = true;
                        affected[y] = true;
                    }
                }
            } else if (e.op == 'I') {
                if (graph[y].dist > graph[x].dist + e.weight) {
                    #pragma omp critical
                    {
                        graph[y].dist = graph[x].dist + e.weight;
                        graph[y].parent = x;
                        affected[y] = true;
                    }
                    // Update graph
                    #pragma omp critical
                    {
                        graph[x].adj.push_back({y, e.weight});
                        graph[y].adj.push_back({x, e.weight});
                    }
                }
            }
        }
    }
}

// Update affected subgraphs (Step 2, Asynchronous)
void update_sssp(std::vector<Vertex>& graph, std::vector<bool>& affected, std::vector<bool>& affected_del,
                 int rank, const std::vector<int>& part, const std::set<int>& ghost_nodes,
                 const std::vector<std::vector<int>>& local_vertices) {
    int level_async = 2; // Level of asynchrony
    bool change = true;
    while (change) {
        change = false;
        // Process deletions
        #pragma omp parallel
        {
            std::queue<int> q;
            int level;
            bool local_change = false;
            #pragma omp for schedule(dynamic) nowait
            for (int v : local_vertices[rank]) {
                if (affected_del[v]) {
                    q.push(v);
                    level = 0;
                    while (!q.empty() && level <= level_async) {
                        int x = q.front();
                        q.pop();
                        for (const auto& [c, _] : graph[x].adj) {
                            if (graph[c].parent == x) {
                                graph[c].dist = INF;
                                graph[c].parent = -1;
                                affected[c] = true;
                                affected_del[c] = true;
                                local_change = true;
                                if (level < level_async) {
                                    q.push(c);
                                }
                            }
                        }
                        level++;
                    }
                }
            }
            if (local_change) {
                #pragma omp critical
                change = true;
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, &change, 1, MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD);
        if (!change) break;
    }

    change = true;
    while (change) {
        change = false;
        #pragma omp parallel
        {
            std::queue<int> q;
            int level;
            bool local_change = false;
            #pragma omp for schedule(dynamic) nowait
            for (int v : local_vertices[rank]) {
                if (affected[v]) {
                    affected[v] = false;
                    q.push(v);
                    level = 0;
                    while (!q.empty() && level <= level_async) {
                        int x = q.front();
                        q.pop();
                        for (const auto& [n, w] : graph[x].adj) {
                            if (graph[x].dist > graph[n].dist + w) {
                                graph[x].dist = graph[n].dist + w;
                                graph[x].parent = n;
                                affected[x] = true;
                                local_change = true;
                                if (level < level_async) {
                                    q.push(x);
                                }
                            }
                            if (graph[n].dist > graph[x].dist + w) {
                                graph[n].dist = graph[x].dist + w;
                                graph[n].parent = x;
                                affected[n] = true;
                                local_change = true;
                                if (level < level_async && part[n] == rank) {
                                    q.push(n);
                                }
                            }
                        }
                        level++;
                    }
                }
            }
            if (local_change) {
                #pragma omp critical
                change = true;
            }
        }
        // Synchronize ghost nodes
        for (int v : ghost_nodes) {
            int owner = part[v];
            if (rank != owner) {
                int data[2] = {graph[v].dist, graph[v].parent};
                MPI_Bcast(data, 2, MPI_INT, owner, MPI_COMM_WORLD);
                graph[v].dist = data[0];
                graph[v].parent = data[1];
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, &change, 1, MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD);
    }
}

// Initialize SSSP tree (sequential Dijkstra for simplicity)
void init_sssp(std::vector<Vertex>& graph, int source) {
    int n = graph.size();
    std::vector<bool> visited(n, false);
    for (auto& v : graph) {
        v.dist = INF;
        v.parent = -1;
        v.affected = false;
        v.affected_del = false;
    }
    graph[source].dist = 0;
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<>> pq;
    pq.push({0, source});
    while (!pq.empty()) {
        int u = pq.top().second;
        pq.pop();
        if (visited[u]) continue;
        visited[u] = true;
        for (const auto& [v, w] : graph[u].adj) {
            if (!visited[v] && graph[u].dist + w < graph[v].dist) {
                graph[v].dist = graph[u].dist + w;
                graph[v].parent = u;
                pq.push({graph[v].dist, v});
            }
        }
    }
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 2) {
        if (rank == 0) std::cerr << "Exactly 2 processes required" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Read graph
    int n, m;
    std::vector<Vertex> graph;
    if (rank == 0) {
        read_graph("graph.metis", n, m, graph);
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != 0) graph.resize(n);
    // Broadcast graph structure (simplified for small graphs)
    for (int i = 0; i < n; ++i) {
        int adj_size = rank == 0 ? graph[i].adj.size() : 0;
        MPI_Bcast(&adj_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (rank != 0) graph[i].adj.resize(adj_size);
        for (int j = 0; j < adj_size; ++j) {
            int v = rank == 0 ? graph[i].adj[j].first : 0;
            int w = rank == 0 ? graph[i].adj[j].second : 0;
            MPI_Bcast(&v, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&w, 1, MPI_INT, 0, MPI_COMM_WORLD);
            graph[i].adj[j] = {v, w};
        }
    }

    // Partition graph
    std::vector<int> part;
    if (rank == 0) {
        partition_graph(n, graph, size, part);
    }
    part.resize(n);
    MPI_Bcast(part.data(), n, MPI_INT, 0, MPI_COMM_WORLD);

    // Identify ghost nodes
    std::set<int> ghost_nodes;
    std::vector<std::vector<int>> local_vertices;
    identify_ghost_nodes(graph, part, rank, size, ghost_nodes, local_vertices);

    // Read changes
    std::vector<Edge> changes;
    if (rank == 0) {
        read_changes("changes.txt", changes);
    }
    int num_changes = rank == 0 ? changes.size() : 0;
    MPI_Bcast(&num_changes, 1, MPI_INT, 0, MPI_COMM_WORLD);
    changes.resize(num_changes);
    for (int i = 0; i < num_changes; ++i) {
        int data[3] = {rank == 0 ? changes[i].u : 0, rank == 0 ? changes[i].v : 0, rank == 0 ? changes[i].weight : 0};
        char op = rank == 0 ? changes[i].op : ' ';
        MPI_Bcast(data, 3, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&op, 1, MPI_CHAR, 0, MPI_COMM_WORLD);
        changes[i] = {data[0], data[1], data[2], op};
    }

    // Initialize SSSP tree
    if (rank == 0) {
        init_sssp(graph, 0); // Source vertex 0
    }
    for (int i = 0; i < n; ++i) {
        int data[2] = {rank == 0 ? graph[i].dist : 0, rank == 0 ? graph[i].parent : 0};
        MPI_Bcast(data, 2, MPI_INT, 0, MPI_COMM_WORLD);
        graph[i].dist = data[0];
        graph[i].parent = data[1];
    }

    // Process changes
    std::vector<bool> affected(n, false), affected_del(n, false);
    double start_time = MPI_Wtime();
    process_changes(graph, changes, rank, part, affected, affected_del);

    // Update SSSP
    update_sssp(graph, affected, affected_del, rank, part, ghost_nodes, local_vertices);
    double end_time = MPI_Wtime();

    // Gather results
    std::vector<int> distances(n), parents(n);
    for (int i = 0; i < n; ++i) {
        distances[i] = graph[i].dist;
        parents[i] = graph[i].parent;
    }
    MPI_Reduce(distances.data(), distances.data(), n, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(parents.data(), parents.data(), n, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

    // Output results
    if (rank == 0) {
        std::ofstream out("sssp_output.txt");
        out << "Vertex\tDistance\tParent\n";
        for (int i = 0; i < n; ++i) {
            out << i << "\t" << (distances[i] == INF ? "INF" : std::to_string(distances[i]))
                << "\t" << parents[i] << "\n";
        }
        out.close();
        std::cout << "Execution time: " << (end_time - start_time) << " seconds\n";
    }

    MPI_Finalize();
    return 0;
}