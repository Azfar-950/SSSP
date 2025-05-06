#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <limits>
#include <ctime>
#include <unistd.h>
#include <stdarg.h>
#include <mpi.h>
#include <omp.h>
#include <algorithm>
#include <cstring>

#define INF std::numeric_limits<int>::max()
#define ROOT_RANK 0
#define MPI_TAG_GHOST_DIST_REQUEST 1001
#define MPI_TAG_GHOST_DIST_REPLY   1002
#define MPI_TAG_GHOST_UPDATE       1003
#define MPI_TAG_DEL_PROPAGATE      1004

#define DEBUG_LEVEL 3

void print_debug_impl(int level, int rank, const char *func, int line, const char *fmt, ...) {
    if (level <= DEBUG_LEVEL) {
        char hostname[MPI_MAX_PROCESSOR_NAME];
        int len;
        MPI_Get_processor_name(hostname, &len);
        std::cout << "[Rank " << rank << " on " << hostname << " | " << func << ":" << line << "] ";
        va_list args;
        va_start(args, fmt);
        vprintf(fmt, args);
        va_end(args);
        std::cout << std::endl;
    }
}
#define PRINT_DEBUG(level, ...) print_debug_impl(level, world_rank, __func__, __LINE__, __VA_ARGS__)

int world_rank = -1;

using idx_t = int;

struct Edge {
    idx_t target_global;
    int weight;
};

struct LocalVertex {
    idx_t global_id;
    int partition_id;
    int distance;
    idx_t parent_global;
    std::vector<Edge> edges;
    bool affected;
    bool affected_del;
    bool just_updated;
    bool is_ghost;
};

struct Change {
    idx_t u;
    idx_t v;
    int weight;
    std::string type;

    void serialize(std::vector<char>& buffer) const {
        size_t type_len = type.size();
        buffer.insert(buffer.end(), reinterpret_cast<const char*>(&u), reinterpret_cast<const char*>(&u) + sizeof(u));
        buffer.insert(buffer.end(), reinterpret_cast<const char*>(&v), reinterpret_cast<const char*>(&v) + sizeof(v));
        buffer.insert(buffer.end(), reinterpret_cast<const char*>(&weight), reinterpret_cast<const char*>(&weight) + sizeof(weight));
        buffer.insert(buffer.end(), reinterpret_cast<const char*>(&type_len), reinterpret_cast<const char*>(&type_len) + sizeof(type_len));
        buffer.insert(buffer.end(), type.begin(), type.end());
    }

    void deserialize(const char* buffer, size_t& offset) {
        memcpy(&u, buffer + offset, sizeof(u));
        offset += sizeof(u);
        memcpy(&v, buffer + offset, sizeof(v));
        offset += sizeof(v);
        memcpy(&weight, buffer + offset, sizeof(weight));
        offset += sizeof(weight);
        size_t type_len;
        memcpy(&type_len, buffer + offset, sizeof(type_len));
        offset += sizeof(type_len);
        type.assign(buffer + offset, buffer + offset + type_len);
        offset += type_len;
    }
};

struct LocalGraphPartition {
    idx_t num_local_vertices;
    idx_t num_ghost_vertices;
    idx_t total_vertices;
    std::vector<LocalVertex> vertices;
    std::vector<idx_t> local_to_global_id;
    std::map<idx_t, idx_t> global_to_local;
    idx_t num_global_vertices;
    std::vector<idx_t> partition_assignment;
};

bool read_changes(const std::string& filename, std::vector<Change>& changes) {
    std::ifstream f(filename);
    if (!f.is_open()) {
        std::cerr << "Error opening changes file: " << filename << std::endl;
        return false;
    }
    PRINT_DEBUG(2, "Reading changes from %s", filename.c_str());

    std::string line, type;
    idx_t u, v;
    int w;
    int line_num = 0;
    while (std::getline(f, line)) {
        ++line_num;
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        if (line.empty()) {
            PRINT_DEBUG(2, "Skipping empty line %d", line_num);
            continue;
        }
        std::istringstream iss(line);
        if (!(iss >> type)) {
            PRINT_DEBUG(2, "Skipping invalid line %d: %s", line_num, line.c_str());
            continue;
        }
        Change ch;
        ch.type = type;
        if (type == "INS" && iss >> u >> v >> w) {
            ch.u = u;
            ch.v = v;
            ch.weight = w;
            changes.push_back(ch);
        } else if (type == "DEL" && iss >> u >> v) {
            ch.u = u;
            ch.v = v;
            ch.weight = 0;
            changes.push_back(ch);
        } else {
            PRINT_DEBUG(2, "Malformed change at line %d: %s", line_num, line.c_str());
            continue;
        }
    }

    f.close();
    PRINT_DEBUG(2, "Read %zu changes.", changes.size());
    return true;
}

void distribute_graph(int rank, int size, LocalGraphPartition& local_partition) {
    PRINT_DEBUG(2, "Starting graph distribution process for rank %d.", rank);

    std::string nodes_filename = "subgraph_" + std::to_string(rank) + "_nodes.txt";
    std::string edges_filename = "subgraph_" + std::to_string(rank) + ".txt";

    std::ifstream nodes_file(nodes_filename);
    if (!nodes_file.is_open()) {
        std::cerr << "[Rank " << rank << "] Error opening nodes file: " << nodes_filename << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    std::vector<idx_t> local_nodes, ghost_nodes;
    std::string line;
    bool reading_local = false, reading_ghost = false;
    idx_t max_vertex = -1;
    while (std::getline(nodes_file, line)) {
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        if (line == "# Local nodes") {
            reading_local = true;
            reading_ghost = false;
        } else if (line == "# Ghost nodes") {
            reading_local = false;
            reading_ghost = true;
        } else if (!line.empty() && (reading_local || reading_ghost)) {
            idx_t node;
            std::istringstream iss(line);
            if (iss >> node) {
                if (reading_local) local_nodes.push_back(node);
                else ghost_nodes.push_back(node);
                max_vertex = std::max(max_vertex, node);
            }
        }
    }
    nodes_file.close();

    PRINT_DEBUG(2, "Local nodes for rank %d:", rank);
    for (idx_t node : local_nodes) {
        PRINT_DEBUG(2, "Node %ld", static_cast<long>(node));
    }
    PRINT_DEBUG(2, "Ghost nodes for rank %d:", rank);
    for (idx_t node : ghost_nodes) {
        PRINT_DEBUG(2, "Node %ld", static_cast<long>(node));
    }

    local_partition.num_local_vertices = local_nodes.size();
    local_partition.num_ghost_vertices = ghost_nodes.size();
    local_partition.total_vertices = local_partition.num_local_vertices + local_partition.num_ghost_vertices;

    std::ifstream edges_file(edges_filename);
    if (!edges_file.is_open()) {
        std::cerr << "[Rank " << rank << "] Error opening edges file: " << edges_filename << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    std::vector<std::tuple<idx_t, idx_t, int>> edges;
    idx_t u, v;
    int w;
    while (edges_file >> u >> v >> w) {
        if (w < 0) {
            std::cerr << "[Rank " << rank << "] Error: Negative weight found in edge (" << u << ", " << v << "): " << w << std::endl;
        }
        edges.emplace_back(u, v, w);
        max_vertex = std::max(max_vertex, std::max(u, v));
    }
    edges_file.close();

    idx_t global_max_vertex;
    MPI_Allreduce(&max_vertex, &global_max_vertex, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    local_partition.num_global_vertices = global_max_vertex + 1;
    PRINT_DEBUG(2, "Global max_vertex: %ld, num_global_vertices: %ld", static_cast<long>(global_max_vertex), static_cast<long>(local_partition.num_global_vertices));

    local_partition.vertices.resize(local_partition.total_vertices);
    local_partition.local_to_global_id.resize(local_partition.total_vertices);
    local_partition.global_to_local.clear();

    idx_t idx = 0;
    for (idx_t node : local_nodes) {
        local_partition.local_to_global_id[idx] = node;
        local_partition.vertices[idx].global_id = node;
        local_partition.vertices[idx].partition_id = rank;
        local_partition.vertices[idx].is_ghost = false;
        local_partition.global_to_local[node] = idx;
        idx++;
    }
    for (idx_t node : ghost_nodes) {
        local_partition.local_to_global_id[idx] = node;
        local_partition.vertices[idx].global_id = node;
        local_partition.vertices[idx].is_ghost = true;
        local_partition.global_to_local[node] = idx;
        idx++;
    }

    for (idx_t i = 0; i < local_partition.num_local_vertices; ++i) {
        idx_t u_global = local_partition.local_to_global_id[i];
        PRINT_DEBUG(2, "Loading edges for local vertex %ld", static_cast<long>(u_global));
        for (const auto& [u, v, w] : edges) {
            if (u == u_global) {
                Edge e;
                e.target_global = v;
                e.weight = w;
                local_partition.vertices[i].edges.push_back(e);
                PRINT_DEBUG(2, "Edge: %ld -> %ld (weight=%d)", static_cast<long>(u_global), static_cast<long>(v), w);
            } else if (v == u_global) {
                Edge e;
                e.target_global = u;
                e.weight = w;
                local_partition.vertices[i].edges.push_back(e);
                PRINT_DEBUG(2, "Edge: %ld -> %ld (weight=%d)", static_cast<long>(u_global), static_cast<long>(u), w);
            }
        }
        PRINT_DEBUG(2, "Vertex %ld has %zu edges", static_cast<long>(u_global), local_partition.vertices[i].edges.size());
    }

    local_partition.partition_assignment.resize(local_partition.num_global_vertices, -1);
    std::vector<idx_t> send_buffer(local_partition.num_global_vertices, -1);
    for (idx_t node : local_nodes) {
        if (node >= 0 && node < local_partition.num_global_vertices) {
            send_buffer[node] = rank;
        }
    }

    PRINT_DEBUG(2, "Send buffer before MPI_Allreduce:");
    for (idx_t i = 0; i < local_partition.num_global_vertices; ++i) {
        PRINT_DEBUG(2, "Vertex %ld -> Rank %ld", static_cast<long>(i), static_cast<long>(send_buffer[i]));
    }

    std::vector<idx_t> recv_buffer(local_partition.num_global_vertices);
    MPI_Allreduce(send_buffer.data(), recv_buffer.data(), local_partition.num_global_vertices, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    local_partition.partition_assignment = recv_buffer;

    PRINT_DEBUG(2, "Partition assignment before verification:");
    for (idx_t i = 0; i < local_partition.num_global_vertices; ++i) {
        PRINT_DEBUG(2, "Vertex %ld -> Rank %ld", static_cast<long>(i), static_cast<long>(local_partition.partition_assignment[i]));
    }

    for (idx_t i = 0; i < local_partition.num_global_vertices; ++i) {
        if (local_partition.partition_assignment[i] == -1) {
            std::cerr << "[Rank " << rank << "] Error: Vertex " << i << " has no owner in partition assignment" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        if (local_partition.partition_assignment[i] >= size) {
            std::cerr << "[Rank " << rank << "] Error: Invalid partition ID " << local_partition.partition_assignment[i] << " for vertex " << i << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    for (idx_t i = local_partition.num_local_vertices; i < local_partition.total_vertices; ++i) {
        idx_t global_id = local_partition.local_to_global_id[i];
        local_partition.vertices[i].partition_id = local_partition.partition_assignment[global_id];
    }

    PRINT_DEBUG(1, "Graph distribution completed: %ld local vertices, %ld ghost vertices, %ld total vertices, %ld global vertices.",
                static_cast<long>(local_partition.num_local_vertices), static_cast<long>(local_partition.num_ghost_vertices),
                static_cast<long>(local_partition.total_vertices), static_cast<long>(local_partition.num_global_vertices));
}

void cleanup_local_partition(LocalGraphPartition& lp) {
    lp.vertices.clear();
    lp.local_to_global_id.clear();
    lp.global_to_local.clear();
    lp.partition_assignment.clear();
    PRINT_DEBUG(2, "Cleaned up local graph partition data.");
}

void initialize_sssp(LocalGraphPartition& lp, idx_t source_global_id) {
    PRINT_DEBUG(2, "Initializing SSSP distances and parents.");
    #pragma omp parallel for
    for (idx_t i = 0; i < lp.total_vertices; ++i) {
        lp.vertices[i].distance = (lp.vertices[i].global_id == source_global_id) ? 0 : INF;
        lp.vertices[i].parent_global = -1;
        lp.vertices[i].affected = (lp.vertices[i].global_id == source_global_id);
        lp.vertices[i].affected_del = false;
        lp.vertices[i].just_updated = false;
    }
    PRINT_DEBUG(2, "SSSP initialization complete.");
}

idx_t find_local_index(const LocalGraphPartition& lp, idx_t global_id) {
    auto it = lp.global_to_local.find(global_id);
    return (it != lp.global_to_local.end()) ? it->second : -1;
}

int find_owner_rank(const LocalGraphPartition& lp, idx_t global_id) {
    if (global_id < 0 || global_id >= lp.num_global_vertices) {
        PRINT_DEBUG(1, "Error: Global ID %ld out of bounds (0 to %ld)", static_cast<long>(global_id), static_cast<long>(lp.num_global_vertices - 1));
        return -1;
    }
    return lp.partition_assignment[global_id];
}

void propagate_deletions_distributed(LocalGraphPartition& lp, int rank, int size) {
    PRINT_DEBUG(1, "Starting deletion propagation phase.");

    std::set<idx_t> local_affected_del;
    for (idx_t i = 0; i < lp.total_vertices; ++i) {
        if (lp.vertices[i].affected_del && !lp.vertices[i].is_ghost) {
            local_affected_del.insert(lp.vertices[i].global_id);
        }
    }

    while (!local_affected_del.empty()) {
        std::set<idx_t> next_level;
        for (const auto& v_global : local_affected_del) {
            idx_t v_local = find_local_index(lp, v_global);
            if (v_local == -1) continue;

            for (idx_t j = 0; j < lp.total_vertices; ++j) {
                if (lp.vertices[j].parent_global == v_global) {
                    idx_t c_global = lp.vertices[j].global_id;
                    int owner = find_owner_rank(lp, c_global);
                    if (owner == rank) {
                        lp.vertices[j].distance = INF;
                        lp.vertices[j].parent_global = -1;
                        lp.vertices[j].affected_del = true;
                        lp.vertices[j].affected = true;
                        next_level.insert(c_global);
                    } else {
                        MPI_Send(&c_global, 1, MPI_INT, owner, MPI_TAG_DEL_PROPAGATE, MPI_COMM_WORLD);
                    }
                }
            }
        }
        local_affected_del = std::move(next_level);

        MPI_Status status;
        int flag;
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_TAG_DEL_PROPAGATE, MPI_COMM_WORLD, &flag, &status);
        while (flag) {
            idx_t c_global;
            MPI_Recv(&c_global, 1, MPI_INT, status.MPI_SOURCE, MPI_TAG_DEL_PROPAGATE, MPI_COMM_WORLD, &status);
            idx_t c_local = find_local_index(lp, c_global);
            if (c_local != -1 && !lp.vertices[c_local].is_ghost) {
                lp.vertices[c_local].distance = INF;
                lp.vertices[c_local].parent_global = -1;
                lp.vertices[c_local].affected_del = true;
                lp.vertices[c_local].affected = true;
                local_affected_del.insert(c_global);
            }
            MPI_Iprobe(MPI_ANY_SOURCE, MPI_TAG_DEL_PROPAGATE, MPI_COMM_WORLD, &flag, &status);
        }
    }

    PRINT_DEBUG(1, "Deletion propagation phase completed.");
}

void relax_edges_distributed(LocalGraphPartition& lp, int rank, int size, bool& any_local_change_made) {
    PRINT_DEBUG(1, "Starting relaxation phase.");

    any_local_change_made = false;

    // Step 1: Identify ghost vertices whose distances are needed
    std::set<idx_t> requested_ghosts;
    for (idx_t i = 0; i < lp.total_vertices; ++i) {
        if (lp.vertices[i].affected && !lp.vertices[i].is_ghost) {
            for (const auto& edge : lp.vertices[i].edges) {
                idx_t neighbor_global = edge.target_global;
                int owner = find_owner_rank(lp, neighbor_global);
                if (owner != rank) {
                    requested_ghosts.insert(neighbor_global);
                }
            }
        }
    }

    // Step 2: Communicate number of requests each rank will send
    int local_request_count = requested_ghosts.size();
    std::vector<int> request_counts(size);
    MPI_Allgather(&local_request_count, 1, MPI_INT, request_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
    int total_requests_expected = 0;
    for (int i = 0; i < size; ++i) {
        total_requests_expected += request_counts[i];
    }
    PRINT_DEBUG(3, "Total requests expected: %d", total_requests_expected);

    // Step 3: Send requests for ghost vertex distances
    std::vector<MPI_Request> send_requests;
    std::map<idx_t, int> request_destinations;
    for (const auto& ghost_global : requested_ghosts) {
        int owner = find_owner_rank(lp, ghost_global);
        PRINT_DEBUG(3, "Sending DIST_REQUEST for ghost %ld to rank %d", static_cast<long>(ghost_global), owner);
        MPI_Request req;
        MPI_Isend(&ghost_global, 1, MPI_INT, owner, MPI_TAG_GHOST_DIST_REQUEST, MPI_COMM_WORLD, &req);
        send_requests.push_back(req);
        request_destinations[ghost_global] = owner;
    }

    // Step 4: Process incoming requests and send replies
    std::vector<MPI_Request> reply_requests;
    int requests_received = 0;
    while (requests_received < total_requests_expected) {
        MPI_Status status;
        int flag;
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_TAG_GHOST_DIST_REQUEST, MPI_COMM_WORLD, &flag, &status);
        if (flag) {
            idx_t requested_global;
            MPI_Recv(&requested_global, 1, MPI_INT, status.MPI_SOURCE, MPI_TAG_GHOST_DIST_REQUEST, MPI_COMM_WORLD, &status);
            idx_t local_idx = find_local_index(lp, requested_global);
            if (local_idx != -1 && !lp.vertices[local_idx].is_ghost) {
                int distance = lp.vertices[local_idx].distance;
                PRINT_DEBUG(3, "Replying to DIST_REQUEST for vertex %ld to rank %d with distance %d", static_cast<long>(requested_global), status.MPI_SOURCE, distance);
                MPI_Request req;
                MPI_Isend(&distance, 1, MPI_INT, status.MPI_SOURCE, MPI_TAG_GHOST_DIST_REQUEST, MPI_COMM_WORLD, &req);
                reply_requests.push_back(req);
            }
            requests_received++;
        }
    }

    // Wait for all send and reply requests to complete
    for (auto& req : send_requests) {
        MPI_Wait(&req, MPI_STATUS_IGNORE);
    }
    for (auto& req : reply_requests) {
        MPI_Wait(&req, MPI_STATUS_IGNORE);
    }

    // Synchronize before receiving replies
    MPI_Barrier(MPI_COMM_WORLD);

    // Step 5: Receive replies for requested ghost vertices
    std::map<idx_t, int> ghost_distances;
    std::vector<MPI_Request> recv_requests;
    std::vector<int*> recv_buffers;
    for (const auto& ghost_global : requested_ghosts) {
        int owner = request_destinations[ghost_global];
        int* distance = new int;
        recv_buffers.push_back(distance);
        MPI_Request req;
        MPI_Irecv(distance, 1, MPI_INT, owner, MPI_TAG_GHOST_DIST_REPLY, MPI_COMM_WORLD, &req);
        recv_requests.push_back(req);
    }

    // Wait for all receives to complete
    for (size_t i = 0; i < recv_requests.size(); ++i) {
        MPI_Status status;
        MPI_Wait(&recv_requests[i], &status);
        idx_t ghost_global = *std::next(requested_ghosts.begin(), i);
        ghost_distances[ghost_global] = *recv_buffers[i];
        PRINT_DEBUG(3, "Received DIST_REPLY for ghost %ld from rank %d: distance=%d", static_cast<long>(ghost_global), status.MPI_SOURCE, *recv_buffers[i]);
    }

    // Clean up receive buffers
    for (auto ptr : recv_buffers) {
        delete ptr;
    }

    // Step 6: Relax edges
    bool local_change_this_phase = false;
    std::vector<bool> currently_affected(lp.total_vertices);
    for (idx_t i = 0; i < lp.total_vertices; ++i) {
        currently_affected[i] = lp.vertices[i].affected;
    }

    #pragma omp parallel for schedule(dynamic) reduction(||:local_change_this_phase)
    for (idx_t i = 0; i < lp.total_vertices; ++i) {
        if (!currently_affected[i]) continue;
        lp.vertices[i].affected = false;
        idx_t v_global = lp.vertices[i].global_id;
        int current_dist_v = lp.vertices[i].distance;
        int thread_id = omp_get_thread_num();

        for (const auto& edge : lp.vertices[i].edges) {
            idx_t neighbor_global = edge.target_global;
            int weight = edge.weight;
            int current_dist_n = INF;

            idx_t neighbor_local_idx = find_local_index(lp, neighbor_global);
            if (neighbor_local_idx != -1) {
                current_dist_n = lp.vertices[neighbor_local_idx].distance;
            } else {
                auto it = ghost_distances.find(neighbor_global);
                if (it != ghost_distances.end()) {
                    current_dist_n = it->second;
                }
            }

            if (current_dist_v != INF) {
                int potential_dist_n = (current_dist_v == INF || weight == INF) ? INF : current_dist_v + weight;
                if (potential_dist_n < current_dist_n) {
                    int owner = find_owner_rank(lp, neighbor_global);
                    if (owner == rank && neighbor_local_idx != -1) {
                        #pragma omp critical(sssp_update)
                        {
                            if (potential_dist_n < lp.vertices[neighbor_local_idx].distance) {
                                PRINT_DEBUG(3, "[Thread %d] Relax %ld->%ld: Dist[%ld] %d -> %d", thread_id, static_cast<long>(v_global), static_cast<long>(neighbor_global), static_cast<long>(neighbor_global), lp.vertices[neighbor_local_idx].distance, potential_dist_n);
                                lp.vertices[neighbor_local_idx].distance = potential_dist_n;
                                lp.vertices[neighbor_local_idx].parent_global = v_global;
                                lp.vertices[neighbor_local_idx].affected = true;
                                lp.vertices[neighbor_local_idx].just_updated = true;
                                local_change_this_phase = true;
                            }
                        }
                    }
                }
            }

            if (current_dist_n != INF) {
                int potential_dist_v = (current_dist_n == INF || weight == INF) ? INF : current_dist_n + weight;
                if (potential_dist_v < current_dist_v && !lp.vertices[i].is_ghost) {
                    #pragma omp critical(sssp_update)
                    {
                        if (potential_dist_v < lp.vertices[i].distance) {
                            PRINT_DEBUG(3, "[Thread %d] Relax %ld->%ld: Dist[%ld] %d -> %d", thread_id, static_cast<long>(neighbor_global), static_cast<long>(v_global), static_cast<long>(v_global), lp.vertices[i].distance, potential_dist_v);
                            lp.vertices[i].distance = potential_dist_v;
                            lp.vertices[i].parent_global = neighbor_global;
                            lp.vertices[i].affected = true;
                            lp.vertices[i].just_updated = true;
                            local_change_this_phase = true;
                        }
                    }
                }
            }
        }
    }

    // Step 7: Send updates for ghost vertices
    std::vector<MPI_Request> update_requests;
    for (idx_t i = 0; i < lp.total_vertices; ++i) {
        if (lp.vertices[i].just_updated && !lp.vertices[i].is_ghost) {
            for (int r = 0; r < size; r++) {
                if (r != rank) {
                    PRINT_DEBUG(3, "Sending GHOST_UPDATE for vertex %ld to rank %d", static_cast<long>(lp.vertices[i].global_id), r);
                    MPI_Request req1, req2, req3;
                    MPI_Isend(&lp.vertices[i].global_id, 1, MPI_INT, r, MPI_TAG_GHOST_UPDATE, MPI_COMM_WORLD, &req1);
                    MPI_Isend(&lp.vertices[i].distance, 1, MPI_INT, r, MPI_TAG_GHOST_UPDATE, MPI_COMM_WORLD, &req2);
                    MPI_Isend(&lp.vertices[i].parent_global, 1, MPI_INT, r, MPI_TAG_GHOST_UPDATE, MPI_COMM_WORLD, &req3);
                    update_requests.push_back(req1);
                    update_requests.push_back(req2);
                    update_requests.push_back(req3);
                }
            }
            lp.vertices[i].just_updated = false;
        }
    }

    // Step 8: Receive updates for ghost vertices
    int update_flag;
    MPI_Status update_status;
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_TAG_GHOST_UPDATE, MPI_COMM_WORLD, &update_flag, &update_status);
    while (update_flag) {
        idx_t global_id;
        int distance;
        idx_t parent_global;
        MPI_Recv(&global_id, 1, MPI_INT, update_status.MPI_SOURCE, MPI_TAG_GHOST_UPDATE, MPI_COMM_WORLD, &update_status);
        MPI_Recv(&distance, 1, MPI_INT, update_status.MPI_SOURCE, MPI_TAG_GHOST_UPDATE, MPI_COMM_WORLD, &update_status);
        MPI_Recv(&parent_global, 1, MPI_INT, update_status.MPI_SOURCE, MPI_TAG_GHOST_UPDATE, MPI_COMM_WORLD, &update_status);
        idx_t local_idx = find_local_index(lp, global_id);
        if (local_idx != -1 && lp.vertices[local_idx].is_ghost) {
            PRINT_DEBUG(3, "Received GHOST_UPDATE for ghost %ld: distance=%d, parent=%ld", static_cast<long>(global_id), distance, static_cast<long>(parent_global));
            lp.vertices[local_idx].distance = distance;
            lp.vertices[local_idx].parent_global = parent_global;
            lp.vertices[local_idx].affected = true;
            any_local_change_made = true; // Trigger further relaxation
        }
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_TAG_GHOST_UPDATE, MPI_COMM_WORLD, &update_flag, &update_status);
    }

    // Wait for all update sends to complete
    for (auto& req : update_requests) {
        MPI_Wait(&req, MPI_STATUS_IGNORE);
    }

    any_local_change_made |= local_change_this_phase;
    PRINT_DEBUG(1, "Relaxation phase finished. Local change detected: %s", any_local_change_made ? "YES" : "NO");
}

void apply_changes_step1(LocalGraphPartition& local_partition, const std::vector<Change>& changes, int rank) {
    PRINT_DEBUG(2, "Applying %zu changes (Step 1).", changes.size());

    // Process deletions sequentially to avoid threading issues
    PRINT_DEBUG(2, "Processing deletions in Step 1.");
    for (size_t i = 0; i < changes.size(); ++i) {
        if (changes[i].type != "DEL") continue;

        PRINT_DEBUG(2, "Processing change %zu: %s %lld %lld", i, 
                    changes[i].type.c_str(), static_cast<long long>(changes[i].u), 
                    static_cast<long long>(changes[i].v));
        idx_t u = changes[i].u;
        idx_t v = changes[i].v;

        if (u < 0 || u >= local_partition.num_global_vertices || 
            v < 0 || v >= local_partition.num_global_vertices) {
            PRINT_DEBUG(2, "Skipping invalid vertices u=%lld or v=%lld", 
                        static_cast<long long>(u), static_cast<long long>(v));
            continue;
        }

        idx_t u_local = local_partition.global_to_local.count(u) ? local_partition.global_to_local[u] : -1;
        idx_t v_local = local_partition.global_to_local.count(v) ? local_partition.global_to_local[v] : -1;

        PRINT_DEBUG(2, "u_local=%lld, v_local=%lld, total_vertices=%lld", 
                    static_cast<long long>(u_local), static_cast<long long>(v_local), 
                    static_cast<long long>(local_partition.total_vertices));

        if (u_local != -1 && u_local < local_partition.total_vertices && 
            local_partition.vertices[u_local].partition_id == rank && !local_partition.vertices[u_local].is_ghost) {
            auto& edges = local_partition.vertices[u_local].edges;
            PRINT_DEBUG(2, "Before deletion from %lld to %lld: %zu edges", 
                        static_cast<long long>(u), static_cast<long long>(v), edges.size());
            for (const auto& e : edges) {
                PRINT_DEBUG(2, "Edge: %lld -> %lld (weight=%d)", 
                            static_cast<long long>(u), static_cast<long long>(e.target_global), e.weight);
            }

            std::vector<Edge> new_edges;
            bool edge_deleted = false;
            for (const auto& e : edges) {
                if (e.target_global != v) {
                    new_edges.push_back(e);
                } else {
                    edge_deleted = true;
                }
            }
            edges = std::move(new_edges);

            if (edge_deleted) {
                PRINT_DEBUG(2, "Deleted edge from %lld to %lld", 
                            static_cast<long long>(u), static_cast<long long>(v));
                local_partition.vertices[u_local].affected_del = true;
                local_partition.vertices[u_local].affected = true;
            } else {
                PRINT_DEBUG(2, "Edge from %lld to %lld does not exist", 
                            static_cast<long long>(u), static_cast<long long>(v));
            }

            PRINT_DEBUG(2, "After deletion from %lld to %lld: %zu edges", 
                        static_cast<long long>(u), static_cast<long long>(v), edges.size());
            for (const auto& e : edges) {
                PRINT_DEBUG(2, "Edge: %lld -> %lld (weight=%d)", 
                            static_cast<long long>(u), static_cast<long long>(e.target_global), e.weight);
            }
        }

        if (v_local != -1 && v_local < local_partition.total_vertices && 
            local_partition.vertices[v_local].partition_id == rank && !local_partition.vertices[v_local].is_ghost) {
            auto& edges = local_partition.vertices[v_local].edges;
            PRINT_DEBUG(2, "Before deletion from %lld to %lld: %zu edges", 
                        static_cast<long long>(v), static_cast<long long>(u), edges.size());
            for (const auto& e : edges) {
                PRINT_DEBUG(2, "Edge: %lld -> %lld (weight=%d)", 
                            static_cast<long long>(v), static_cast<long long>(e.target_global), e.weight);
            }

            std::vector<Edge> new_edges;
            bool edge_deleted = false;
            for (const auto& e : edges) {
                if (e.target_global != u) {
                    new_edges.push_back(e);
                } else {
                    edge_deleted = true;
                }
            }
            edges = std::move(new_edges);

            if (edge_deleted) {
                PRINT_DEBUG(2, "Deleted edge from %lld to %lld", 
                            static_cast<long long>(v), static_cast<long long>(u));
                local_partition.vertices[v_local].affected_del = true;
                local_partition.vertices[v_local].affected = true;
            } else {
                PRINT_DEBUG(2, "Edge from %lld to %lld does not exist", 
                            static_cast<long long>(v), static_cast<long long>(u));
            }

            PRINT_DEBUG(2, "After deletion from %lld to %lld: %zu edges", 
                        static_cast<long long>(v), static_cast<long long>(u), edges.size());
            for (const auto& e : edges) {
                PRINT_DEBUG(2, "Edge: %lld -> %lld (weight=%d)", 
                            static_cast<long long>(v), static_cast<long long>(e.target_global), e.weight);
            }
        }
    }

    // Process insertions in parallel
    PRINT_DEBUG(2, "Processing insertions in Step 1.");
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        #pragma omp for schedule(dynamic)
        for (size_t i = 0; i < changes.size(); ++i) {
            if (changes[i].type != "INS") continue;

            PRINT_DEBUG(2, "[Thread %d] Processing change %zu: %s %lld %lld weight=%d", tid, i, 
                        changes[i].type.c_str(), static_cast<long long>(changes[i].u), 
                        static_cast<long long>(changes[i].v), changes[i].weight);
            idx_t u = changes[i].u;
            idx_t v = changes[i].v;

            if (u < 0 || u >= local_partition.num_global_vertices || 
                v < 0 || v >= local_partition.num_global_vertices) {
                PRINT_DEBUG(2, "[Thread %d] Skipping invalid vertices u=%lld or v=%lld", 
                            tid, static_cast<long long>(u), static_cast<long long>(v));
                continue;
            }

            idx_t u_local = local_partition.global_to_local.count(u) ? local_partition.global_to_local[u] : -1;
            idx_t v_local = local_partition.global_to_local.count(v) ? local_partition.global_to_local[v] : -1;

            if (u_local != -1 && u_local < local_partition.total_vertices && 
                local_partition.vertices[u_local].partition_id == rank && !local_partition.vertices[u_local].is_ghost) {
                #pragma omp critical
                {
                    auto& edges = local_partition.vertices[u_local].edges;
                    Edge e;
                    e.target_global = v;
                    e.weight = changes[i].weight;
                    edges.push_back(e);
                    PRINT_DEBUG(2, "[Thread %d] Inserted edge from %lld to %lld weight=%d", tid, 
                                static_cast<long long>(u), static_cast<long long>(v), e.weight);
                    local_partition.vertices[u_local].affected = true;
                }
            }

            if (v_local != -1 && v_local < local_partition.total_vertices && 
                local_partition.vertices[v_local].partition_id == rank && !local_partition.vertices[v_local].is_ghost) {
                #pragma omp critical
                {
                    auto& edges = local_partition.vertices[v_local].edges;
                    Edge e;
                    e.target_global = u;
                    e.weight = changes[i].weight;
                    edges.push_back(e);
                    PRINT_DEBUG(2, "[Thread %d] Inserted edge from %lld to %lld weight=%d", tid, 
                                static_cast<long long>(v), static_cast<long long>(u), e.weight);
                    local_partition.vertices[v_local].affected = true;
                }
            }
        }
    }

    PRINT_DEBUG(2, "Finished processing changes locally in Step 1.");
}

void update_distances_step2(LocalGraphPartition& lp, int rank, int size) {
    PRINT_DEBUG(1, "Starting Step 2 iterative update process.");

    PRINT_DEBUG(1, "Step 2a: Propagating deletion effects (distributed)...");
    propagate_deletions_distributed(lp, rank, size);
    PRINT_DEBUG(1, "Step 2a: Finished deletion propagation phase.");

    MPI_Barrier(MPI_COMM_WORLD);

    PRINT_DEBUG(1, "Step 2b: Starting iterative edge relaxation...");
    bool global_change_occurred_in_iteration = true;
    int iteration = 0;
    int max_iterations = lp.num_global_vertices * 2; // Ensure enough iterations for distributed graph

    while (global_change_occurred_in_iteration && iteration < max_iterations) {
        iteration++;
        PRINT_DEBUG(1, "--- Relaxation Iteration %d ---", iteration);
        bool local_change_occurred = false;

        relax_edges_distributed(lp, rank, size, local_change_occurred);

        // Synchronize ghost updates across ranks
        MPI_Allreduce(&local_change_occurred, &global_change_occurred_in_iteration, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);

        if (rank == ROOT_RANK) {
            PRINT_DEBUG(1, "Global change in iteration %d: %s", iteration, global_change_occurred_in_iteration ? "YES" : "NO");
        }
    }

    if (iteration >= max_iterations && global_change_occurred_in_iteration) {
        if (rank == ROOT_RANK) {
            std::cerr << "Warning: Reached maximum iterations (" << max_iterations << ") without convergence." << std::endl;
        }
    }
    PRINT_DEBUG(1, "Step 2b: Iterative relaxation finished after %d iterations.", iteration);
}

void gather_and_print_results(const LocalGraphPartition& lp, idx_t num_global_vertices, int rank, int size) {
    PRINT_DEBUG(1, "Gathering final results to Rank 0.");

    std::vector<long> send_buffer(lp.num_local_vertices * 3);
    int send_count = 0;
    for (idx_t i = 0; i < lp.total_vertices; ++i) {
        if (!lp.vertices[i].is_ghost) {
            send_buffer[send_count++] = static_cast<long>(lp.vertices[i].global_id);
            send_buffer[send_count++] = static_cast<long>(lp.vertices[i].distance);
            send_buffer[send_count++] = static_cast<long>(lp.vertices[i].parent_global);
        }
    }

    int items_per_vertex = 3;
    int send_items = lp.num_local_vertices * items_per_vertex;
    std::vector<int> recv_counts(size);
    MPI_Gather(&send_items, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, ROOT_RANK, MPI_COMM_WORLD);

    std::vector<long> recv_buffer;
    std::vector<int> displs;
    int total_items_received = 0;

    if (rank == ROOT_RANK) {
        displs.resize(size);
        displs[0] = 0;
        total_items_received = recv_counts[0];
        for (int i = 1; i < size; ++i) {
            displs[i] = displs[i - 1] + recv_counts[i - 1];
            total_items_received += recv_counts[i];
        }
        recv_buffer.resize(total_items_received);
        PRINT_DEBUG(2, "Root receiving %d total items for results.", total_items_received);
    }

    MPI_Gatherv(send_buffer.data(), send_items, MPI_LONG, recv_buffer.data(), recv_counts.data(), displs.data(), MPI_LONG, ROOT_RANK, MPI_COMM_WORLD);

    if (rank == ROOT_RANK) {
        std::cout << "\n--- Final SSSP Results (gathered at Rank 0) ---\n";
        std::cout << "GlobalID\tDistance\tParentID\n";
        std::cout << "------------------------------------------------\n";
        for (int i = 0; i < total_items_received; i += 3) {
            long gid = recv_buffer[i];
            long dist = recv_buffer[i + 1];
            long parent = recv_buffer[i + 2];
            std::cout << gid << "\t\t" << (dist == INF ? -1 : dist) << "\t\t" << parent << "\n";
        }
        std::cout << "------------------------------------------------\n";
    }
    PRINT_DEBUG(1, "Result gathering complete.");
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if (argc < 4) {
        if (world_rank == ROOT_RANK) {
            std::cerr << "Usage: mpirun -np <num_procs> " << argv[0] << " <graph_edge_file> <changes_file> <source_vertex_id>" << std::endl;
        }
        MPI_Finalize();
        return 1;
    }

    std::string graph_filename = argv[1];
    std::string changes_filename = argv[2];
    idx_t source_global_id = std::atoi(argv[3]);

    PRINT_DEBUG(1, "Application started with %d processes.", world_size);

    std::vector<Change> changes;
    LocalGraphPartition local_partition;
    double start_time, end_time, setup_time, initial_sssp_time, update_time;

    start_time = MPI_Wtime();

    if (world_rank == ROOT_RANK) {
        PRINT_DEBUG(1, "Reading changes from %s", changes_filename.c_str());
        if (!read_changes(changes_filename, changes)) {
            std::cerr << "[Rank " << world_rank << "] Error reading changes file." << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        PRINT_DEBUG(1, "Read %zu change operations.", changes.size());
    }

    distribute_graph(world_rank, world_size, local_partition);
    PRINT_DEBUG(1, "Graph distribution complete. Local partition has %ld local and %ld ghost vertices (total %ld).",
                static_cast<long>(local_partition.num_local_vertices), static_cast<long>(local_partition.num_ghost_vertices),
                static_cast<long>(local_partition.total_vertices));

    // Broadcast number of changes
    int num_changes = changes.size();
    MPI_Bcast(&num_changes, 1, MPI_INT, ROOT_RANK, MPI_COMM_WORLD);
    PRINT_DEBUG(2, "Received broadcast: num_changes = %d", num_changes);

    // Serialize and broadcast changes
    std::vector<char> buffer;
    if (world_rank == ROOT_RANK) {
        for (const auto& ch : changes) {
            ch.serialize(buffer);
        }
    }
    int buffer_size = buffer.size();
    MPI_Bcast(&buffer_size, 1, MPI_INT, ROOT_RANK, MPI_COMM_WORLD);
    buffer.resize(buffer_size);
    MPI_Bcast(buffer.data(), buffer_size, MPI_CHAR, ROOT_RANK, MPI_COMM_WORLD);

    if (world_rank != ROOT_RANK) {
        changes.clear();
        size_t offset = 0;
        for (int i = 0; i < num_changes; ++i) {
            Change ch;
            ch.deserialize(buffer.data(), offset);
            changes.push_back(ch);
        }
    }
    PRINT_DEBUG(2, "Received broadcast of changes data.");

    setup_time = MPI_Wtime() - start_time;
    PRINT_DEBUG(1, "Setup phase completed in %.4f seconds.", setup_time);

    initialize_sssp(local_partition, source_global_id);
    PRINT_DEBUG(1, "Initial SSSP distances set (source = %ld).", static_cast<long>(source_global_id));

    // Compute initial SSSP distances
    MPI_Barrier(MPI_COMM_WORLD);
    initial_sssp_time = MPI_Wtime();
    PRINT_DEBUG(1, "*** Starting Initial SSSP Computation ***");
    update_distances_step2(local_partition, world_rank, world_size);
    PRINT_DEBUG(1, "*** Finished Initial SSSP Computation ***");
    initial_sssp_time = MPI_Wtime() - initial_sssp_time;
    PRINT_DEBUG(1, "Initial SSSP phase completed in %.4f seconds.", initial_sssp_time);

    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();

    PRINT_DEBUG(1, "*** Starting SSSP Update Step 1: Applying Changes ***");
    apply_changes_step1(local_partition, changes, world_rank);
    PRINT_DEBUG(1, "*** Finished SSSP Update Step 1 ***");

    MPI_Barrier(MPI_COMM_WORLD);

    PRINT_DEBUG(1, "*** Starting SSSP Update Step 2: Iterative Update ***");
    update_distances_step2(local_partition, world_rank, world_size);
    PRINT_DEBUG(1, "*** Finished SSSP Update Step 2 ***");

    update_time = MPI_Wtime() - start_time;
    PRINT_DEBUG(1, "SSSP Update phase completed in %.4f seconds.", update_time);

    MPI_Barrier(MPI_COMM_WORLD);

    idx_t final_num_global_vertices = local_partition.num_global_vertices;
    gather_and_print_results(local_partition, final_num_global_vertices, world_rank, world_size);

    PRINT_DEBUG(1, "Cleaning up resources...");
    cleanup_local_partition(local_partition);

    end_time = MPI_Wtime();
    if (world_rank == ROOT_RANK) {
        std::cout << "\n----------------------------------------\n";
        std::cout << "Total execution time: " << (end_time - (start_time - setup_time - initial_sssp_time)) << " seconds\n";
        std::cout << "  Setup time:         " << setup_time << " seconds\n";
        std::cout << "  Initial SSSP time:  " << initial_sssp_time << " seconds\n";
        std::cout << "  Update time:        " << update_time << " seconds\n";
        std::cout << "----------------------------------------\n";
    }

    MPI_Finalize();
    return 0;
}