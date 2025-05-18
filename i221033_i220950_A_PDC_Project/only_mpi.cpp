#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <fstream>
#include <sstream>
#include <limits>
#include <string>
#include <algorithm>
#include <chrono>

using namespace std;

// Constants
const int MPI_TAG_GHOST_UPDATE = 100;
const int MAX_ITERATIONS = 1000; // Safety limit
const int INF = numeric_limits<int>::max();

// Edge structure
struct Edge
{
    int to;
    int weight;
    Edge(int t, int w) : to(t), weight(w) {}
};

// Graph structure
class Graph
{
public:
    map<int, vector<Edge>> adj_list;
    void add_edge(int from, int to, int weight)
    {
        remove_edge(from, to); // Prevent duplicates
        adj_list[from].emplace_back(to, weight);
        adj_list[to].emplace_back(from, weight); // Undirected graph
    }
    void remove_edge(int from, int to)
    {
        auto &edges = adj_list[from];
        edges.erase(remove_if(edges.begin(), edges.end(),
                              [to](const Edge &e)
                              { return e.to == to; }),
                    edges.end());
        auto &reverse_edges = adj_list[to];
        reverse_edges.erase(remove_if(reverse_edges.begin(), reverse_edges.end(),
                                      [from](const Edge &e)
                                      { return e.to == from; }),
                            reverse_edges.end());
    }
    const vector<Edge> &get_edges(int vertex) const
    {
        static vector<Edge> empty;
        auto it = adj_list.find(vertex);
        return it != adj_list.end() ? it->second : empty;
    }
};

// Change operation structure
struct ChangeOperation
{
    int u, v;
    int weight;
    bool is_deletion;
};

// Read subgraph from file
Graph read_subgraph(const string &filename, int rank)
{
    Graph graph;
    ifstream file(filename);
    if (!file.is_open())
    {
        if (rank == 0)
            cerr << "Error opening subgraph file: " << filename << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    string line;
    while (getline(file, line))
    {
        istringstream iss(line);
        int u, v, w;
        if (iss >> u >> v >> w)
        {
            graph.add_edge(u, v, w);
        }
    }
    file.close();
    return graph;
}

// Read nodes from subgraph_<rank>_nodes.txt
void read_nodes(const string &filename, int rank, set<int> &local_vertices, set<int> &ghost_vertices, int &max_vertex)
{
    ifstream file(filename);
    if (!file.is_open())
    {
        if (rank == 0)
            cerr << "Error opening nodes file: " << filename << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    string line;
    bool reading_local = false, reading_ghost = false;
    max_vertex = -1;
    while (getline(file, line))
    {
        if (line == "# Local nodes")
        {
            reading_local = true;
            reading_ghost = false;
            continue;
        }
        else if (line == "# Ghost nodes")
        {
            reading_local = false;
            reading_ghost = true;
            continue;
        }
        istringstream iss(line);
        int node;
        if (iss >> node)
        {
            max_vertex = max(max_vertex, node);
            if (reading_local)
            {
                local_vertices.insert(node);
            }
            else if (reading_ghost)
            {
                ghost_vertices.insert(node);
            }
        }
    }
    file.close();
}

// Read changes from file
vector<ChangeOperation> read_changes(const string &filename, int rank)
{
    vector<ChangeOperation> changes;
    ifstream file(filename);
    if (!file.is_open())
    {
        if (rank == 0)
            cerr << "Error opening changes file: " << filename << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    string line;
    while (getline(file, line))
    {
        istringstream iss(line);
        string op;
        int u, v, w = 0;
        if (iss >> op >> u >> v)
        {
            ChangeOperation change;
            change.u = u;
            change.v = v;
            change.is_deletion = (op == "DEL");
            if (!change.is_deletion)
            {
                iss >> w;
                change.weight = w;
            }
            else
            {
                change.weight = 0;
            }
            changes.push_back(change);
        }
    }
    file.close();
    if (rank == 0)
    {
        cout << "[Rank " << rank << "] Read " << changes.size() << " changes." << endl;
    }
    return changes;
}

// Process pending GHOST_UPDATE messages
void process_ghost_updates(int rank, int size, const map<int, int> &global_to_local_map,
                           const map<int, int> &global_to_owner_map,
                           map<int, int> &distance,
                           map<int, int> &parent,
                           map<int, int> &ghost_distances,
                           map<int, bool> &affected,
                           bool &ghost_change,
                           const set<int> &local_vertices,
                           const set<int> &ghost_vertices,
                           const Graph &local_graph)
{
    MPI_Status status;
    ghost_change = false;
    set<int> affected_neighbors;

    // Estimate number of messages (conservative: one update per ghost vertex per other rank)
    int max_messages = ghost_vertices.size() * (size - 1);
    for (int i = 0; i < max_messages; ++i)
    {
        int flag;
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_TAG_GHOST_UPDATE, MPI_COMM_WORLD, &flag, &status);
        if (!flag)
            break;
        int global_id, new_dist;
        MPI_Recv(&global_id, 1, MPI_INT, status.MPI_SOURCE, MPI_TAG_GHOST_UPDATE, MPI_COMM_WORLD, &status);
        cout << "[Rank " << rank << "] Received GHOST_UPDATE for vertex " << global_id << " from rank " << status.MPI_SOURCE << endl;
        MPI_Recv(&new_dist, 1, MPI_INT, status.MPI_SOURCE, MPI_TAG_GHOST_UPDATE, MPI_COMM_WORLD, &status);
        if (ghost_distances.count(global_id) && global_to_owner_map.at(global_id) != rank)
        {
            if (new_dist < ghost_distances[global_id])
            {
                ghost_distances[global_id] = new_dist;
                affected[global_id] = true;
                ghost_change = true;
                cout << "[Rank " << rank << "] Updated ghost distance for vertex " << global_id << " to " << new_dist << endl;
                // Mark local vertices that have this ghost vertex as a neighbor
                for (int u : local_vertices)
                {
                    for (const auto &edge : local_graph.get_edges(u))
                    {
                        if (edge.to == global_id)
                        {
                            affected_neighbors.insert(u);
                        }
                    }
                }
            }
        }
    }
    // Mark affected local vertices
    for (int u : affected_neighbors)
    {
        affected[u] = true;
    }
}
void distribute_graph(int rank, int size, const string &subgraph_file, const string &nodes_file,
                      const vector<ChangeOperation> &changes,
                      map<int, int> &global_to_local_map,
                      map<int, int> &global_to_owner_map,
                      set<int> &local_vertices,
                      set<int> &ghost_vertices,
                      Graph &local_graph,
                      int &max_vertex,
                      int &num_global_vertices)
{
    // Read local and ghost vertices
    read_nodes(nodes_file, rank, local_vertices, ghost_vertices, max_vertex);

    // Compute global max_vertex across all ranks
    int global_max_vertex;
    MPI_Allreduce(&max_vertex, &global_max_vertex, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    max_vertex = global_max_vertex;
    num_global_vertices = max_vertex + 1;

    // Read local subgraph
    local_graph = read_subgraph(subgraph_file, rank);

    // Apply changes to local graph
    for (const auto &change : changes)
    {
        if (local_vertices.count(change.u) || local_vertices.count(change.v))
        {
            if (change.is_deletion)
            {
                local_graph.remove_edge(change.u, change.v);
                cout << "[Rank " << rank << "] Removed edge " << change.u << "->" << change.v << endl;
            }
            else
            {
                local_graph.add_edge(change.u, change.v, change.weight);
                cout << "[Rank " << rank << "] Added edge " << change.u << "->" << change.v << " (weight " << change.weight << ")" << endl;
            }
        }
    }

    // Recompute ghost vertices after applying changes
    ghost_vertices.clear();
    for (int u : local_vertices)
    {
        for (const auto &edge : local_graph.get_edges(u))
        {
            int v = edge.to;
            if (local_vertices.count(v) == 0)
            {
                ghost_vertices.insert(v);
            }
        }
    }

    // Populate global_to_owner_map and global_to_local_map
    for (int r = 0; r < size; ++r)
    {
        set<int> rank_local_vertices, rank_ghost_vertices;
        int rank_max_vertex;
        read_nodes("subgraph_" + to_string(r) + "_nodes.txt", rank, rank_local_vertices, rank_ghost_vertices, rank_max_vertex);
        for (int v : rank_local_vertices)
        {
            global_to_owner_map[v] = r;
            if (r == rank)
            {
                global_to_local_map[v] = v;
            }
        }
    }

    // Debug local graph
    cout << "[Rank " << rank << "] Local graph edges after changes:" << endl;
    for (int u : local_vertices)
    {
        for (const auto &edge : local_graph.get_edges(u))
        {
            cout << "[Rank " << rank << "] Edge " << u << "->" << edge.to << " (weight " << edge.weight << ")" << endl;
        }
    }
    cout << "[Rank " << rank << "] Graph distribution completed: " << local_vertices.size() << " local vertices, "
         << ghost_vertices.size() << " ghost vertices, " << num_global_vertices << " global vertices." << endl;
}
// Initialize SSSP
void initialize_sssp(int rank, int source, const set<int> &local_vertices,
                     const set<int> &ghost_vertices,
                     map<int, int> &distance, map<int, int> &parent,
                     map<int, int> &ghost_distances, map<int, bool> &affected)
{
    distance.clear();
    parent.clear();
    ghost_distances.clear();
    affected.clear();
    for (int u : local_vertices)
    {
        distance[u] = (u == source) ? 0 : INF;
        parent[u] = -1;
        affected[u] = (u == source);
    }
    for (int u : ghost_vertices)
    {
        ghost_distances[u] = INF;
        affected[u] = false;
    }
    if (ghost_vertices.count(source))
    {
        ghost_distances[source] = 0;
        affected[source] = true;
    }
    cout << "[Rank " << rank << "] SSSP initialization complete. Local vertices: ";
    for (int u : local_vertices)
        cout << u << " ";
    cout << "\n[Rank " << rank << "] Initial ghost distances: ";
    for (int u : ghost_vertices)
        cout << u << ":" << (ghost_distances[u] == INF ? "INF" : to_string(ghost_distances[u])) << " ";
    cout << endl;
}

// Apply changes (Step 1)
void apply_changes_step1(int rank, int size, const vector<ChangeOperation> &changes,
                         const map<int, int> &global_to_owner_map,
                         const set<int> &local_vertices,
                         Graph &local_graph, map<int, bool> &affected)
{
    for (const auto &change : changes)
    {
        if (local_vertices.count(change.u) || local_vertices.count(change.v))
        {
            if (change.is_deletion)
            {
                local_graph.remove_edge(change.u, change.v);
                affected[change.u] = true;
                affected[change.v] = true;
                cout << "[Rank " << rank << "] Removed edge " << change.u << "->" << change.v << endl;
            }
            else
            {
                local_graph.add_edge(change.u, change.v, change.weight);
                affected[change.u] = true;
                affected[change.v] = true;
                cout << "[Rank " << rank << "] Added edge " << change.u << "->" << change.v << " (weight " << change.weight << ")" << endl;
            }
        }
    }
    // Rebuild local graph edges for affected vertices
    cout << "[Rank " << rank << "] Local graph edges after applying changes:" << endl;
    for (const auto &[u, edges] : local_graph.adj_list)
    {
        for (const auto &edge : edges)
        {
            cout << "[Rank " << rank << "] Edge " << u << "->" << edge.to << " (weight " << edge.weight << ")" << endl;
        }
    }
    cout << "[Rank " << rank << "] Applied relevant changes." << endl;
}

// Relax edges distributed
void relax_edges_distributed(int rank, int size, const Graph &local_graph,
                             const map<int, int> &global_to_local_map,
                             const map<int, int> &global_to_owner_map,
                             map<int, int> &distance,
                             map<int, int> &parent,
                             map<int, bool> &affected,
                             map<int, int> &ghost_distances,
                             bool &local_change,
                             const set<int> &local_vertices,
                             const set<int> &ghost_vertices)
{
    bool ghost_change = false;
    process_ghost_updates(rank, size, global_to_local_map, global_to_owner_map,
                          distance, parent, ghost_distances, affected, ghost_change,
                          local_vertices, ghost_vertices, local_graph);

    local_change = false;
    vector<int> updated_vertices;
    map<int, bool> new_affected;
    set<int> sent_updates;

    for (auto it = local_vertices.begin(); it != local_vertices.end(); ++it)
    {
        int v_global = *it;
        bool vertex_changed = false;
        for (const auto &[u_global, edges] : local_graph.adj_list)
        {
            for (const auto &edge : edges)
            {
                if (edge.to == v_global)
                {
                    int weight = edge.weight;
                    int dist_u = global_to_local_map.count(u_global) ? distance[u_global] : ghost_distances.at(u_global);
                    int dist_v = distance[v_global];
                    if (dist_u != INF && dist_v > dist_u + weight)
                    {
#pragma omp critical
                        {
                            int new_dist = dist_u + weight;
                            distance[v_global] = new_dist;
                            parent[v_global] = u_global;
                            new_affected[v_global] = true;
                            updated_vertices.push_back(v_global);
                            vertex_changed = true;
                            local_change = true;
                            cout << "[Rank " << rank << "] Relax " << u_global << "->" << v_global
                                 << ": Dist[" << v_global << "] updated to " << new_dist << endl;
                        }
                    }
                }
            }
        }
        for (const auto &edge : local_graph.get_edges(v_global))
        {
            int w_global = edge.to;
            int weight = edge.weight;
            int dist_v = distance[v_global];
            int dist_w = global_to_local_map.count(w_global) ? distance[w_global] : ghost_distances.at(w_global);
            if (dist_v != INF && dist_w > dist_v + weight)
            {
#pragma omp critical
                {
                    int new_dist = dist_v + weight;
                    if (global_to_local_map.count(w_global))
                    {
                        distance[w_global] = new_dist;
                        parent[w_global] = v_global;
                        new_affected[w_global] = true;
                        updated_vertices.push_back(w_global);
                        vertex_changed = true;
                        local_change = true;
                        cout << "[Rank " << rank << "] Relax " << v_global << "->" << w_global
                             << ": Dist[" << w_global << "] updated to " << new_dist << endl;
                    }
                    else
                    {
                        if (ghost_distances[w_global] > new_dist)
                        {
                            ghost_distances[w_global] = new_dist;
                            new_affected[w_global] = true;
                            vertex_changed = true;
                            local_change = true;
                            cout << "[Rank " << rank << "] Updated ghost distance for vertex " << w_global
                                 << " to " << new_dist << endl;
                        }
                    }
                }
            }
        }
        if (vertex_changed)
        {
            new_affected[v_global] = true;
        }
    }

    vector<MPI_Request> send_requests;
    for (int u_global : updated_vertices)
    {
        if (sent_updates.find(u_global) == sent_updates.end())
        {
            sent_updates.insert(u_global);
            for (int r = 0; r < size; ++r)
            {
                if (r != rank)
                {
                    send_requests.emplace_back();
                    MPI_Isend(&u_global, 1, MPI_INT, r, MPI_TAG_GHOST_UPDATE, MPI_COMM_WORLD, &send_requests.back());
                    send_requests.emplace_back();
                    MPI_Isend(&distance[u_global], 1, MPI_INT, r, MPI_TAG_GHOST_UPDATE, MPI_COMM_WORLD, &send_requests.back());
                    cout << "[Rank " << rank << "] Sending GHOST_UPDATE for vertex " << u_global << " to rank " << r << endl;
                }
            }
        }
    }
    if (!send_requests.empty())
    {
        MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUSES_IGNORE);
    }

    affected = new_affected;

    cout << "[Rank " << rank << "] Distances: ";
    for (const auto &[v, d] : distance)
    {
        cout << v << ":" << (d == INF ? "INF" : to_string(d)) << " ";
    }
    cout << endl;
    cout << "[Rank " << rank << "] Ghost distances: ";
    for (const auto &[v, d] : ghost_distances)
    {
        cout << v << ":" << (d == INF ? "INF" : to_string(d)) << " ";
    }
    cout << endl;
    cout << "[Rank " << rank << "] Affected vertices: ";
    for (const auto &[v, a] : affected)
    {
        if (a)
            cout << v << " ";
    }
    cout << endl;
    cout << "[Rank " << rank << "] local_change: " << (local_change ? "YES" : "NO") << ", ghost_change: " << (ghost_change ? "YES" : "NO") << endl;
}

// Function to save the modified graph to a new graph.txt
void save_modified_graph(const map<int, vector<Edge>>& local_graph, 
    const set<int>& local_vertices, 
    int rank, int size, 
    const string& filename) {
// Step 1: Gather local edges
vector<string> local_edges;
for (const auto& vertex : local_vertices) {
for (const auto& edge : local_graph.at(vertex)) {
local_edges.push_back(to_string(vertex) + " " + 
            to_string(edge.to) + " " + 
            to_string(edge.weight));
}
}

// Step 2: Compute number of edges per rank
int local_edge_count = local_edges.size();
vector<int> edge_counts(size);
MPI_Allgather(&local_edge_count, 1, MPI_INT, 
edge_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

// Step 3: Compute displacements for gathering edges
vector<int> displs(size);
int total_edges = 0;
for (int i = 0; i < size; ++i) {
displs[i] = total_edges;
total_edges += edge_counts[i];
}

// Step 4: Serialize local edges into a single string
string local_edge_str;
for (const auto& edge : local_edges) {
local_edge_str += edge + "\n";
}

// Step 5: Compute send buffer size
vector<int> str_lengths(size);
int local_str_length = local_edge_str.size();
MPI_Allgather(&local_str_length, 1, MPI_INT, 
str_lengths.data(), 1, MPI_INT, MPI_COMM_WORLD);

// Step 6: Compute displacements for string data
vector<int> str_displs(size);
int total_str_length = 0;
for (int i = 0; i < size; ++i) {
str_displs[i] = total_str_length;
total_str_length += str_lengths[i];
}

// Step 7: Gather all edges at Rank 0
vector<char> all_edges_str(total_str_length);
MPI_Gatherv(local_edge_str.c_str(), local_str_length, MPI_CHAR,
all_edges_str.data(), str_lengths.data(), str_displs.data(),
MPI_CHAR, 0, MPI_COMM_WORLD);

// Step 8: Write to file on Rank 0
if (rank == 0) {
ofstream out_file(filename);
if (!out_file.is_open()) {
cerr << "[Rank 0] Error: Could not open " << filename << " for writing" << endl;
MPI_Abort(MPI_COMM_WORLD, 1);
}

string all_edges(all_edges_str.begin(), all_edges_str.end());
istringstream iss(all_edges);
string line;
set<string> unique_edges; // To avoid duplicates
while (getline(iss, line)) {
if (!line.empty()) {
unique_edges.insert(line);
}
}

for (const auto& edge : unique_edges) {
out_file << edge << endl;
}

out_file.close();
cout << "[Rank 0] Modified graph saved to " << filename << endl;
}
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Get number of OpenMP threads
    int num_threads = omp_get_max_threads();

    if (argc != 3)
    {
        if (rank == 0)
        {
            cerr << "Usage: " << argv[0] << " <changes_file> <source_vertex>" << endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    string changes_file = argv[1];
    int source = stoi(argv[2]);

    string subgraph_file = "subgraph_" + to_string(rank) + ".txt";
    string nodes_file = "subgraph_" + to_string(rank) + "_nodes.txt";
    vector<ChangeOperation> changes = read_changes(changes_file, rank);

    // Start timing
    auto start_time = chrono::high_resolution_clock::now();

    map<int, int> global_to_local_map, global_to_owner_map;
    set<int> local_vertices, ghost_vertices;
    Graph local_graph;
    int max_vertex, num_global_vertices;
    distribute_graph(rank, size, subgraph_file, nodes_file, changes, global_to_local_map, global_to_owner_map,
                     local_vertices, ghost_vertices, local_graph, max_vertex, num_global_vertices);

    // Check if source vertex exists
    if (global_to_owner_map.find(source) == global_to_owner_map.end()) {
        if (rank == 0) {
            cerr << "Error: Source vertex " << source << " does not exist in the graph." << endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    map<int, int> distance, parent, ghost_distances;
    map<int, bool> affected;
    initialize_sssp(rank, source, local_vertices, ghost_vertices, distance, parent, ghost_distances, affected);

    apply_changes_step1(rank, size, changes, global_to_owner_map, local_vertices, local_graph, affected);

    bool global_change = true;
    int iteration = 0;
    while (global_change && iteration < MAX_ITERATIONS)
    {
        bool local_change = false;
        relax_edges_distributed(rank, size, local_graph, global_to_local_map,
                                global_to_owner_map, distance, parent,
                                affected, ghost_distances, local_change,
                                local_vertices, ghost_vertices);

        int local_flag = local_change ? 1 : 0;
        int total_change = 0;
        MPI_Allreduce(&local_flag, &total_change, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
        global_change = (total_change != 0);
        iteration++;
        cout << "[Rank " << rank << "] Iteration " << iteration << " completed. Global change: "
             << (global_change ? "YES" : "NO") << ", total_change: " << total_change << endl;
    }

    // End timing
    auto end_time = chrono::high_resolution_clock::now();
    double execution_time = chrono::duration<double>(end_time - start_time).count();

    // Log performance data on rank 0
    if (rank == 0)
    {
        ofstream out_file("performance_log.txt", ios::app);
        if (!out_file.is_open())
        {
            cerr << "[Rank 0] Error: Could not open performance_log.txt for writing" << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        out_file << "MPI Processes: " << size
                 << ", OpenMP Threads: " << num_threads
                 << ", Execution Time (s): " << execution_time << endl;
        out_file.close();
        cout << "[Rank 0] Performance data appended to performance_log.txt" << endl;
    }

    cout << "[Rank " << rank << "] Final ghost distances: ";
    for (const auto &[v, d] : ghost_distances)
    {
        cout << v << ":" << (d == INF ? "INF" : to_string(d)) << " ";
    }
    cout << endl;

    vector<int> all_distances(num_global_vertices, INF);
    vector<int> all_parents(num_global_vertices, -1);
    for (int u : local_vertices)
    {
        all_distances[u] = distance[u];
        all_parents[u] = parent[u];
    }
    vector<int> global_distances(num_global_vertices);
    vector<int> global_parents(num_global_vertices);
    MPI_Reduce(all_distances.data(), global_distances.data(), num_global_vertices, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(all_parents.data(), global_parents.data(), num_global_vertices, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0)
    {
        cout << "*** Final SSSP Distances ***" << endl;
        for (int u = 0; u < num_global_vertices; ++u)
        {
            if (global_distances[u] != INF)
            {
                cout << "Vertex " << u << ": Distance = " << global_distances[u]
                     << ", Parent = " << global_parents[u] << endl;
            }
        }
    }

    MPI_Finalize();
    return 0;
}