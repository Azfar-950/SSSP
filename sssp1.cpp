// #include <iostream>
// #include <fstream>
// #include <sstream>
// #include <string>
// #include <vector>
// #include <set>
// #include <map>
// #include <limits>
// #include <ctime>
// #include <unistd.h>
// #include <stdarg.h>
// #include <mpi.h>
// #include <omp.h>
// #include <metis.h>

// // Forward declarations to resolve "not declared in scope" errors
// struct LocalGraphPartition;
// idx_t find_local_index(const LocalGraphPartition& lp, idx_t global_id);
// void propagate_deletions_distributed(LocalGraphPartition& lp, int rank, int size);
// void relax_edges_distributed(LocalGraphPartition& lp, int rank, int size, bool& any_local_change_made);

// #define INF std::numeric_limits<int>::max()
// #define ROOT_RANK 0
// #define MPI_TAG_GHOST_DIST_REQUEST 1001
// #define MPI_TAG_GHOST_DIST_REPLY   1002
// #define MPI_TAG_GHOST_UPDATE       1003
// #define MPI_TAG_DEL_PROPAGATE      1004

// #define DEBUG_LEVEL 2

// void print_debug_impl(int level, int rank, const char *func, int line, const char *fmt, ...) {
//     if (level <= DEBUG_LEVEL) {
//         char hostname[MPI_MAX_PROCESSOR_NAME];
//         int len;
//         MPI_Get_processor_name(hostname, &len);
//         std::cout << "[Rank " << rank << " on " << hostname << " | " << func << ":" << line << "] ";
//         va_list args;
//         va_start(args, fmt);
//         vprintf(fmt, args);
//         va_end(args);
//         std::cout << std::endl;
//     }
// }
// #define PRINT_DEBUG(level, ...) print_debug_impl(level, world_rank, __func__, __LINE__, __VA_ARGS__)

// int world_rank = -1;

// using idx_t = int;  // Assuming idx_t is int, adjust if necessary

// struct Edge {
//     idx_t target_global;
//     int weight;
// };

// struct LocalVertex {
//     idx_t global_id;
//     int partition_id;
//     int distance;
//     idx_t parent_global;
//     std::vector<Edge> edges;
//     bool affected;
//     bool affected_del;
//     bool just_updated;
//     bool is_ghost;
// };

// struct Change {
//     idx_t u;
//     idx_t v;
//     int weight;
//     std::string type;
// };

// struct LocalGraphPartition {
//     idx_t num_local_vertices;
//     idx_t num_ghost_vertices;
//     idx_t total_vertices;
//     std::vector<LocalVertex> vertices;
//     std::vector<idx_t> local_to_global_id;
//     idx_t num_global_vertices;
//     std::vector<idx_t> partition_assignment;
// };

// struct GraphCSR {
//     idx_t nvtxs;
//     idx_t nedges;
//     std::vector<idx_t> xadj;
//     std::vector<idx_t> adjncy;
//     std::vector<idx_t> adjwgt;
// };

// bool read_graph_edges(const std::string& filename, idx_t& num_vertices, idx_t& num_edges, std::vector<idx_t>& u, std::vector<idx_t>& v, std::vector<int>& w) {
//     std::ifstream f(filename);
//     if (!f.is_open()) {
//         std::cerr << "Error opening graph file: " << filename << std::endl;
//         return false;
//     }
//     PRINT_DEBUG(2, "Reading edges from %s", filename.c_str());

//     idx_t u_i, v_i, max_id = 0;
//     int w_i;
//     idx_t edge_count = 0;
//     std::string line;
//     std::vector<idx_t> u_temp, v_temp;
//     std::vector<int> w_temp;

//     while (std::getline(f, line)) {
//         std::istringstream iss(line);
//         if (iss >> u_i >> v_i >> w_i) {
//             edge_count++;
//             u_temp.push_back(u_i);
//             v_temp.push_back(v_i);
//             w_temp.push_back(w_i);
//             if (u_i > max_id) max_id = u_i;
//             if (v_i > max_id) max_id = v_i;
//         }
//     }

//     num_vertices = max_id + 1;
//     num_edges = edge_count;
//     u = std::move(u_temp);
//     v = std::move(v_temp);
//     w = std::move(w_temp);

//     f.close();
//     PRINT_DEBUG(2, "Read %ld edges, max vertex ID %ld", static_cast<long>(edge_count), static_cast<long>(max_id));
//     return true;
// }

// void convert_edges_to_csr(idx_t nvtxs, idx_t nedges, const std::vector<idx_t>& u, const std::vector<idx_t>& v, const std::vector<int>& w, GraphCSR& csr) {
//     PRINT_DEBUG(2, "Converting %ld vertices, %ld edges to CSR", static_cast<long>(nvtxs), static_cast<long>(nedges));
//     csr.nvtxs = nvtxs;
//     csr.nedges = nedges;
//     csr.xadj.resize(nvtxs + 1, 0);
//     csr.adjncy.resize(nedges);
//     csr.adjwgt.resize(nedges);

//     for (idx_t i = 0; i < nedges; ++i) {
//         csr.xadj[u[i]]++;
//     }

//     idx_t sum = 0;
//     for (idx_t i = 0; i < nvtxs; ++i) {
//         idx_t temp = csr.xadj[i];
//         csr.xadj[i] = sum;
//         sum += temp;
//     }
//     csr.xadj[nvtxs] = sum;

//     std::vector<idx_t> pos(nvtxs, 0);
//     for (idx_t i = 0; i < nedges; ++i) {
//         idx_t src = u[i];
//         idx_t edge_idx = csr.xadj[src] + pos[src];
//         csr.adjncy[edge_idx] = v[i];
//         csr.adjwgt[edge_idx] = w[i];
//         pos[src]++;
//     }

//     PRINT_DEBUG(2, "CSR conversion complete.");
// }

// bool read_changes(const std::string& filename, std::vector<Change>& changes) {
//     std::ifstream f(filename);
//     if (!f.is_open()) {
//         std::cerr << "Error opening changes file: " << filename << std::endl;
//         return false;
//     }
//     PRINT_DEBUG(2, "Reading changes from %s", filename.c_str());

//     std::string line, type;
//     idx_t u, v;
//     int w;
//     while (std::getline(f, line)) {
//         std::istringstream iss(line);
//         if (iss >> type) {
//             Change ch;
//             ch.type = type;
//             if (type == "INS" && iss >> u >> v >> w) {
//                 ch.u = u;
//                 ch.v = v;
//                 ch.weight = w;
//             } else if (type == "DEL" && iss >> u >> v) {
//                 ch.u = u;
//                 ch.v = v;
//                 ch.weight = 0;
//             } else {
//                 std::cerr << "Unknown or malformed change: " << line << std::endl;
//                 continue;
//             }
//             changes.push_back(ch);
//         }
//     }

//     f.close();
//     PRINT_DEBUG(2, "Read %zu changes.", changes.size());
//     return true;
// }

// void partition_graph(const GraphCSR& csr, int num_partitions, std::vector<idx_t>& partition_result) {
//     if (num_partitions <= 1) {
//         PRINT_DEBUG(1, "Number of partitions is 1, skipping METIS call.");
//         partition_result.resize(csr.nvtxs, 0);
//         return;
//     }

//     PRINT_DEBUG(2, "Calling METIS_PartGraphKway...");
//     partition_result.resize(csr.nvtxs);
//     idx_t ncon = 1;
//     idx_t objval;
//     idx_t options[METIS_NOPTIONS];
//     METIS_SetDefaultOptions(options);
//     options[METIS_OPTION_DBGLVL] = 0;
//     options[METIS_OPTION_CONTIG] = 1;

//     // Create mutable copies of vector data
//     std::vector<idx_t> xadj_copy = csr.xadj;
//     std::vector<idx_t> adjncy_copy = csr.adjncy;
//     std::vector<idx_t> adjwgt_copy = csr.adjwgt;

//     int ret = METIS_PartGraphKway(&csr.nvtxs, &ncon, xadj_copy.data(), adjncy_copy.data(),
//                                   adjwgt_copy.data(), nullptr, adjwgt_copy.data(),
//                                   &num_partitions, nullptr, nullptr, options, &objval, partition_result.data());

//     if (ret != METIS_OK) {
//         std::cerr << "METIS_PartGraphKway failed with error code: " << ret << std::endl;
//         MPI_Abort(MPI_COMM_WORLD, 1);
//     }
//     PRINT_DEBUG(1, "METIS partitioning successful. Edge cut: %ld", static_cast<long>(objval));
// }

// void distribute_graph(const GraphCSR& csr, const std::vector<idx_t>& partition_assignment, int rank, int size, LocalGraphPartition& local_partition) {
//     PRINT_DEBUG(2, "Starting graph distribution process.");

//     std::vector<idx_t> local_counts(size, 0);
//     std::vector<idx_t> ghost_counts(size, 0);
//     idx_t num_global_vertices = 0;

//     if (rank == ROOT_RANK) {
//         num_global_vertices = csr.nvtxs;
//         std::vector<bool> is_ghost_for_rank(csr.nvtxs * size, false);
//         std::vector<idx_t> ghost_node_indices_count(size, 0);

//         for (idx_t i = 0; i < csr.nvtxs; ++i) {
//             local_counts[partition_assignment[i]]++;
//         }

//         for (idx_t u_global = 0; u_global < csr.nvtxs; ++u_global) {
//             int owner_rank = partition_assignment[u_global];
//             for (idx_t edge_idx = csr.xadj[u_global]; edge_idx < csr.xadj[u_global + 1]; ++edge_idx) {
//                 idx_t v_global = csr.adjncy[edge_idx];
//                 int target_rank = partition_assignment[v_global];
//                 if (owner_rank != target_rank) {
//                     if (!is_ghost_for_rank[v_global * size + owner_rank]) {
//                         is_ghost_for_rank[v_global * size + owner_rank] = true;
//                         ghost_node_indices_count[owner_rank]++;
//                     }
//                     if (!is_ghost_for_rank[u_global * size + target_rank]) {
//                         is_ghost_for_rank[u_global * size + target_rank] = true;
//                         ghost_node_indices_count[target_rank]++;
//                     }
//                 }
//             }
//         }
//         for (int r = 0; r < size; ++r) ghost_counts[r] = ghost_node_indices_count[r];
//         PRINT_DEBUG(2, "Calculated local/ghost counts on root.");
//     }

//     MPI_Bcast(&num_global_vertices, 1, MPI_LONG_LONG_INT, ROOT_RANK, MPI_COMM_WORLD);
//     MPI_Bcast(local_counts.data(), size, MPI_LONG_LONG_INT, ROOT_RANK, MPI_COMM_WORLD);
//     MPI_Bcast(ghost_counts.data(), size, MPI_LONG_LONG_INT, ROOT_RANK, MPI_COMM_WORLD);

//     local_partition.num_local_vertices = local_counts[rank];
//     local_partition.num_ghost_vertices = ghost_counts[rank];
//     local_partition.total_vertices = local_partition.num_local_vertices + local_partition.num_ghost_vertices;
//     local_partition.num_global_vertices = num_global_vertices;

//     if (rank == ROOT_RANK) {
//         MPI_Bcast(const_cast<idx_t*>(partition_assignment.data()), num_global_vertices, MPI_LONG_LONG_INT, ROOT_RANK, MPI_COMM_WORLD);
//     } else {
//         local_partition.partition_assignment.resize(num_global_vertices);
//         MPI_Bcast(local_partition.partition_assignment.data(), num_global_vertices, MPI_LONG_LONG_INT, ROOT_RANK, MPI_COMM_WORLD);
//     }

//     PRINT_DEBUG(2, "Received counts: local=%ld, ghost=%ld, total=%ld", static_cast<long>(local_counts[rank]), static_cast<long>(ghost_counts[rank]), static_cast<long>(local_partition.total_vertices));

//     local_partition.vertices.resize(local_partition.total_vertices);
//     local_partition.local_to_global_id.resize(local_partition.total_vertices);

//     if (rank == ROOT_RANK) {
//         idx_t local_idx = 0;
//         std::set<idx_t> ghost_set;

//         for (idx_t v_global = 0; v_global < csr.nvtxs; v_global++) {
//             if (partition_assignment[v_global] == rank) {
//                 local_partition.local_to_global_id[local_idx] = v_global;
//                 local_partition.vertices[local_idx].global_id = v_global;
//                 local_partition.vertices[local_idx].partition_id = rank;
//                 local_partition.vertices[local_idx].is_ghost = false;
//                 idx_t num_edges = csr.xadj[v_global + 1] - csr.xadj[v_global];
//                 local_partition.vertices[local_idx].edges.resize(num_edges);
//                 idx_t start = csr.xadj[v_global];
//                 for (idx_t e = 0; e < num_edges; e++) {
//                     local_partition.vertices[local_idx].edges[e].target_global = csr.adjncy[start + e];
//                     local_partition.vertices[local_idx].edges[e].weight = csr.adjwgt[start + e];
//                     if (partition_assignment[csr.adjncy[start + e]] != rank) {
//                         ghost_set.insert(csr.adjncy[start + e]);
//                     }
//                 }
//                 local_idx++;
//             }
//         }

//         idx_t ghost_idx = local_partition.num_local_vertices;
//         for (const auto& ghost : ghost_set) {
//             local_partition.local_to_global_id[ghost_idx] = ghost;
//             local_partition.vertices[ghost_idx].global_id = ghost;
//             local_partition.vertices[ghost_idx].partition_id = partition_assignment[ghost];
//             local_partition.vertices[ghost_idx].is_ghost = true;
//             ghost_idx++;
//         }

//         for (int r = 1; r < size; r++) {
//             idx_t num_local_r = local_counts[r];
//             std::vector<idx_t> local_global_ids;
//             for (idx_t i = 0; i < csr.nvtxs; i++) {
//                 if (partition_assignment[i] == r) {
//                     local_global_ids.push_back(i);
//                 }
//             }
//             MPI_Send(local_global_ids.data(), num_local_r, MPI_LONG_LONG_INT, r, 0, MPI_COMM_WORLD);

//             for (idx_t k = 0; k < num_local_r; k++) {
//                 idx_t v_global = local_global_ids[k];
//                 idx_t num_edges = csr.xadj[v_global + 1] - csr.xadj[v_global];
//                 MPI_Send(&num_edges, 1, MPI_LONG_LONG_INT, r, 0, MPI_COMM_WORLD);
//                 idx_t start = csr.xadj[v_global];
//                 MPI_Send(&csr.adjncy[start], num_edges, MPI_LONG_LONG_INT, r, 0, MPI_COMM_WORLD);
//                 MPI_Send(&csr.adjwgt[start], num_edges, MPI_INT, r, 0, MPI_COMM_WORLD);
//             }

//             std::set<idx_t> ghost_set_r;
//             for (idx_t v_global = 0; v_global < csr.nvtxs; v_global++) {
//                 if (partition_assignment[v_global] == r) {
//                     idx_t start = csr.xadj[v_global];
//                     idx_t end = csr.xadj[v_global + 1];
//                     for (idx_t idx = start; idx < end; idx++) {
//                         idx_t target = csr.adjncy[idx];
//                         if (partition_assignment[target] != r) {
//                             ghost_set_r.insert(target);
//                         }
//                     }
//                 }
//             }
//             idx_t num_ghost_r = ghost_set_r.size();
//             MPI_Send(&num_ghost_r, 1, MPI_LONG_LONG_INT, r, 0, MPI_COMM_WORLD);
//             std::vector<idx_t> ghost_global_ids(ghost_set_r.begin(), ghost_set_r.end());
//             MPI_Send(ghost_global_ids.data(), num_ghost_r, MPI_LONG_LONG_INT, r, 0, MPI_COMM_WORLD);
//         }
//     } else {
//         std::vector<idx_t> local_global_ids(local_partition.num_local_vertices);
//         MPI_Recv(local_global_ids.data(), local_partition.num_local_vertices, MPI_LONG_LONG_INT, ROOT_RANK, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

//         for (idx_t k = 0; k < local_partition.num_local_vertices; k++) {
//             idx_t v_global = local_global_ids[k];
//             local_partition.local_to_global_id[k] = v_global;
//             local_partition.vertices[k].global_id = v_global;
//             local_partition.vertices[k].partition_id = rank;
//             local_partition.vertices[k].is_ghost = false;
//             idx_t num_edges;
//             MPI_Recv(&num_edges, 1, MPI_LONG_LONG_INT, ROOT_RANK, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//             local_partition.vertices[k].edges.resize(num_edges);
//             std::vector<idx_t> targets(num_edges);
//             std::vector<int> weights(num_edges);
//             MPI_Recv(targets.data(), num_edges, MPI_LONG_LONG_INT, ROOT_RANK, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//             MPI_Recv(weights.data(), num_edges, MPI_INT, ROOT_RANK, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//             for (idx_t e = 0; e < num_edges; e++) {
//                 local_partition.vertices[k].edges[e].target_global = targets[e];
//                 local_partition.vertices[k].edges[e].weight = weights[e];
//             }
//         }

//         idx_t num_ghost_r;
//         MPI_Recv(&num_ghost_r, 1, MPI_LONG_LONG_INT, ROOT_RANK, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//         std::vector<idx_t> ghost_global_ids(num_ghost_r);
//         MPI_Recv(ghost_global_ids.data(), num_ghost_r, MPI_LONG_LONG_INT, ROOT_RANK, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//         for (idx_t k = 0; k < num_ghost_r; k++) {
//             idx_t ghost_idx = local_partition.num_local_vertices + k;
//             local_partition.local_to_global_id[ghost_idx] = ghost_global_ids[k];
//             local_partition.vertices[ghost_idx].global_id = ghost_global_ids[k];
//             local_partition.vertices[ghost_idx].partition_id = local_partition.partition_assignment[ghost_global_ids[k]];
//             local_partition.vertices[ghost_idx].is_ghost = true;
//         }
//     }

//     PRINT_DEBUG(1, "Graph distribution completed successfully.");
// }

// void cleanup_csr(GraphCSR& csr) {
//     csr.xadj.clear();
//     csr.adjncy.clear();
//     csr.adjwgt.clear();
//     PRINT_DEBUG(2, "Cleaned up global CSR data.");
// }

// void cleanup_local_partition(LocalGraphPartition& lp) {
//     lp.vertices.clear();
//     lp.local_to_global_id.clear();
//     lp.partition_assignment.clear();
//     PRINT_DEBUG(2, "Cleaned up local graph partition data.");
// }

// void initialize_sssp(LocalGraphPartition& lp, idx_t source_global_id) {
//     PRINT_DEBUG(2, "Initializing SSSP distances and parents.");
//     #pragma omp parallel for
//     for (idx_t i = 0; i < lp.total_vertices; ++i) {
//         lp.vertices[i].distance = (lp.vertices[i].global_id == source_global_id) ? 0 : INF;
//         lp.vertices[i].parent_global = -1;
//         lp.vertices[i].affected = (lp.vertices[i].global_id == source_global_id);
//         lp.vertices[i].affected_del = false;
//         lp.vertices[i].just_updated = false;
//     }
//     PRINT_DEBUG(2, "SSSP initialization complete.");
// }

// void apply_changes_step1(const std::vector<Change>& changes, LocalGraphPartition& lp, int rank, int size) {
//     PRINT_DEBUG(1, "Applying %zu changes (Step 1).", changes.size());

//     #pragma omp parallel for schedule(dynamic)
//     for (size_t i = 0; i < changes.size(); ++i) {
//         const Change& ch = changes[i];
//         idx_t u_global = ch.u;
//         idx_t v_global = ch.v;
//         int thread_id = omp_get_thread_num();
//         PRINT_DEBUG(3, "[Thread %d] Processing change %zu: %s %ld %ld", thread_id, i, ch.type.c_str(), static_cast<long>(u_global), static_cast<long>(v_global));

//         idx_t u_local = find_local_index(lp, u_global);
//         idx_t v_local = find_local_index(lp, v_global);

//         if (ch.type == "DEL") {
//             bool changed = false;
//             if (v_local != -1 && lp.vertices[v_local].parent_global == u_global) {
//                 #pragma omp critical(sssp_update)
//                 {
//                     if (lp.vertices[v_local].distance != INF) {
//                         PRINT_DEBUG(2, "[Thread %d] DEL %ld->%ld detected as SSSP edge. Disconnecting %ld.", thread_id, static_cast<long>(u_global), static_cast<long>(v_global), static_cast<long>(v_global));
//                         lp.vertices[v_local].distance = INF;
//                         lp.vertices[v_local].parent_global = -1;
//                         lp.vertices[v_local].affected_del = true;
//                         lp.vertices[v_local].affected = true;
//                         changed = true;
//                     }
//                 }
//             } else if (u_local != -1 && lp.vertices[u_local].parent_global == v_global) {
//                 #pragma omp critical(sssp_update)
//                 {
//                     if (lp.vertices[u_local].distance != INF) {
//                         PRINT_DEBUG(2, "[Thread %d] DEL %ld->%ld detected as SSSP edge. Disconnecting %ld.", thread_id, static_cast<long>(v_global), static_cast<long>(u_global), static_cast<long>(u_global));
//                         lp.vertices[u_local].distance = INF;
//                         lp.vertices[u_local].parent_global = -1;
//                         lp.vertices[u_local].affected_del = true;
//                         lp.vertices[u_local].affected = true;
//                         changed = true;
//                     }
//                 }
//             }
//         } else if (ch.type == "INS") {
//             int weight = ch.weight;
//             int dist_u = (u_local != -1) ? lp.vertices[u_local].distance : INF;
//             int dist_v = (v_local != -1) ? lp.vertices[v_local].distance : INF;

//             if (v_local != -1 && dist_u != INF) {
//                 int potential_dist_v = (dist_u == INF || weight == INF) ? INF : dist_u + weight;
//                 if (potential_dist_v < dist_v) {
//                     #pragma omp critical(sssp_update)
//                     {
//                         if (potential_dist_v < lp.vertices[v_local].distance) {
//                             PRINT_DEBUG(2, "[Thread %d] INS %ld->%ld (w=%d) updates Dist[%ld]: %d -> %d",
//                                         thread_id, static_cast<long>(u_global), static_cast<long>(v_global), weight, static_cast<long>(v_global), lp.vertices[v_local].distance, potential_dist_v);
//                             lp.vertices[v_local].distance = potential_dist_v;
//                             lp.vertices[v_local].parent_global = u_global;
//                             lp.vertices[v_local].affected = true;
//                         }
//                     }
//                 }
//             }

//             if (u_local != -1 && dist_v != INF) {
//                 int potential_dist_u = (dist_v == INF || weight == INF) ? INF : dist_v + weight;
//                 if (potential_dist_u < dist_u) {
//                     #pragma omp critical(sssp_update)
//                     {
//                         if (potential_dist_u < lp.vertices[u_local].distance) {
//                             PRINT_DEBUG(2, "[Thread %d] INS %ld->%ld (w=%d) updates Dist[%ld]: %d -> %d",
//                                         thread_id, static_cast<long>(v_global), static_cast<long>(u_global), weight, static_cast<long>(u_global), lp.vertices[u_local].distance, potential_dist_u);
//                             lp.vertices[u_local].distance = potential_dist_u;
//                             lp.vertices[u_local].parent_global = v_global;
//                             lp.vertices[u_local].affected = true;
//                         }
//                     }
//                 }
//             }
//         }
//     }
//     PRINT_DEBUG(1, "Finished processing changes locally in Step 1.");
// }

// void update_distances_step2(LocalGraphPartition& lp, int rank, int size) {
//     PRINT_DEBUG(1, "Starting Step 2 iterative update process.");

//     PRINT_DEBUG(1, "Step 2a: Propagating deletion effects (distributed)...");
//     propagate_deletions_distributed(lp, rank, size);
//     PRINT_DEBUG(1, "Step 2a: Finished deletion propagation phase.");

//     MPI_Barrier(MPI_COMM_WORLD);

//     PRINT_DEBUG(1, "Step 2b: Starting iterative edge relaxation...");
//     bool global_change_occurred_in_iteration = true;
//     int iteration = 0;

//     while (global_change_occurred_in_iteration) {
//         iteration++;
//         PRINT_DEBUG(1, "--- Relaxation Iteration %d ---", iteration);
//         bool local_change_occurred = false;

//         relax_edges_distributed(lp, rank, size, local_change_occurred);

//         MPI_Allreduce(&local_change_occurred, &global_change_occurred_in_iteration, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);

//         if (rank == ROOT_RANK) {
//             PRINT_DEBUG(1, "Global change in iteration %d: %s", iteration, global_change_occurred_in_iteration ? "YES" : "NO");
//         }
//         if (iteration > lp.total_vertices * 2 && global_change_occurred_in_iteration) {
//             if (rank == ROOT_RANK) std::cerr << "Warning: Potential non-convergence after " << iteration << " iterations." << std::endl;
//             break;
//         }
//     }
//     PRINT_DEBUG(1, "Step 2b: Iterative relaxation finished after %d iterations.", iteration);
// }

// idx_t find_local_index(const LocalGraphPartition& lp, idx_t global_id) {
//     for (idx_t i = 0; i < lp.total_vertices; ++i) {
//         if (lp.vertices[i].global_id == global_id) {
//             return i;
//         }
//     }
//     return -1;
// }

// int find_owner_rank(const LocalGraphPartition& lp, idx_t global_id) {
//     if (global_id < 0 || global_id >= lp.num_global_vertices) {
//         PRINT_DEBUG(1, "Error: Global ID %ld out of bounds (0 to %ld)", static_cast<long>(global_id), static_cast<long>(lp.num_global_vertices - 1));
//         return -1;
//     }
//     return lp.partition_assignment[global_id];
// }

// void propagate_deletions_distributed(LocalGraphPartition& lp, int rank, int size) {
//     PRINT_DEBUG(1, "Starting deletion propagation phase.");

//     std::set<idx_t> local_affected_del;
//     for (idx_t i = 0; i < lp.total_vertices; ++i) {
//         if (lp.vertices[i].affected_del && !lp.vertices[i].is_ghost) {
//             local_affected_del.insert(lp.vertices[i].global_id);
//         }
//     }

//     while (!local_affected_del.empty()) {
//         std::set<idx_t> next_level;
//         for (const auto& v_global : local_affected_del) {
//             idx_t v_local = find_local_index(lp, v_global);
//             if (v_local == -1) continue;

//             for (idx_t j = 0; j < lp.total_vertices; ++j) {
//                 if (lp.vertices[j].parent_global == v_global) {
//                     idx_t c_global = lp.vertices[j].global_id;
//                     int owner = find_owner_rank(lp, c_global);
//                     if (owner == rank) {
//                         lp.vertices[j].distance = INF;
//                         lp.vertices[j].parent_global = -1;
//                         lp.vertices[j].affected_del = true;
//                         lp.vertices[j].affected = true;
//                         next_level.insert(c_global);
//                     } else {
//                         MPI_Send(&c_global, 1, MPI_LONG_LONG_INT, owner, MPI_TAG_DEL_PROPAGATE, MPI_COMM_WORLD);
//                     }
//                 }
//             }
//         }
//         local_affected_del = std::move(next_level);

//         MPI_Status status;
//         int flag;
//         MPI_Iprobe(MPI_ANY_SOURCE, MPI_TAG_DEL_PROPAGATE, MPI_COMM_WORLD, &flag, &status);
//         while (flag) {
//             idx_t c_global;
//             MPI_Recv(&c_global, 1, MPI_LONG_LONG_INT, status.MPI_SOURCE, MPI_TAG_DEL_PROPAGATE, MPI_COMM_WORLD, &status);
//             idx_t c_local = find_local_index(lp, c_global);
//             if (c_local != -1 && !lp.vertices[c_local].is_ghost) {
//                 lp.vertices[c_local].distance = INF;
//                 lp.vertices[c_local].parent_global = -1;
//                 lp.vertices[c_local].affected_del = true;
//                 lp.vertices[c_local].affected = true;
//                 local_affected_del.insert(c_global);
//             }
//             MPI_Iprobe(MPI_ANY_SOURCE, MPI_TAG_DEL_PROPAGATE, MPI_COMM_WORLD, &flag, &status);
//         }
//     }

//     PRINT_DEBUG(1, "Deletion propagation phase completed.");
// }

// void relax_edges_distributed(LocalGraphPartition& lp, int rank, int size, bool& any_local_change_made) {
//     PRINT_DEBUG(1, "Starting relaxation phase.");

//     any_local_change_made = false;

//     std::set<idx_t> requested_ghosts;
//     for (idx_t i = 0; i < lp.total_vertices; ++i) {
//         if (lp.vertices[i].affected && !lp.vertices[i].is_ghost) {
//             for (const auto& edge : lp.vertices[i].edges) {
//                 idx_t neighbor_global = edge.target_global;
//                 int owner = find_owner_rank(lp, neighbor_global);
//                 if (owner != rank) {
//                     requested_ghosts.insert(neighbor_global);
//                 }
//             }
//         }
//     }

//     for (const auto& ghost_global : requested_ghosts) {
//         int owner = find_owner_rank(lp, ghost_global);
//         MPI_Send(&ghost_global, 1, MPI_LONG_LONG_INT, owner, MPI_TAG_GHOST_DIST_REQUEST, MPI_COMM_WORLD);
//     }

//     MPI_Status status;
//     int flag;
//     MPI_Iprobe(MPI_ANY_SOURCE, MPI_TAG_GHOST_DIST_REQUEST, MPI_COMM_WORLD, &flag, &status);
//     while (flag) {
//         idx_t requested_global;
//         MPI_Recv(&requested_global, 1, MPI_LONG_LONG_INT, status.MPI_SOURCE, MPI_TAG_GHOST_DIST_REQUEST, MPI_COMM_WORLD, &status);
//         idx_t local_idx = find_local_index(lp, requested_global);
//         if (local_idx != -1) {
//             int distance = lp.vertices[local_idx].distance;
//             MPI_Send(&distance, 1, MPI_INT, status.MPI_SOURCE, MPI_TAG_GHOST_DIST_REPLY, MPI_COMM_WORLD);
//         }
//         MPI_Iprobe(MPI_ANY_SOURCE, MPI_TAG_GHOST_DIST_REQUEST, MPI_COMM_WORLD, &flag, &status);
//     }

//     std::map<idx_t, int> ghost_distances;
//     for (const auto& ghost_global : requested_ghosts) {
//         int distance;
//         MPI_Recv(&distance, 1, MPI_INT, find_owner_rank(lp, ghost_global), MPI_TAG_GHOST_DIST_REPLY, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//         ghost_distances[ghost_global] = distance;
//     }

//     bool local_change_this_phase = false;
//     std::vector<bool> currently_affected(lp.total_vertices);
//     for (idx_t i = 0; i < lp.total_vertices; ++i) {
//         currently_affected[i] = lp.vertices[i].affected;
//     }

//     #pragma omp parallel for schedule(dynamic) reduction(||:local_change_this_phase)
//     for (idx_t i = 0; i < lp.total_vertices; ++i) {
//         if (!currently_affected[i]) continue;
//         lp.vertices[i].affected = false;
//         idx_t v_global = lp.vertices[i].global_id;
//         int current_dist_v = lp.vertices[i].distance;
//         int thread_id = omp_get_thread_num();

//         for (const auto& edge : lp.vertices[i].edges) {
//             idx_t neighbor_global = edge.target_global;
//             int weight = edge.weight;
//             int current_dist_n = INF;

//             idx_t neighbor_local_idx = find_local_index(lp, neighbor_global);
//             if (neighbor_local_idx != -1) {
//                 current_dist_n = lp.vertices[neighbor_local_idx].distance;
//             } else {
//                 auto it = ghost_distances.find(neighbor_global);
//                 if (it != ghost_distances.end()) {
//                     current_dist_n = it->second;
//                 }
//             }

//             if (current_dist_v != INF) {
//                 int potential_dist_n = (current_dist_v == INF || weight == INF) ? INF : current_dist_v + weight;
//                 if (potential_dist_n < current_dist_n) {
//                     int owner = find_owner_rank(lp, neighbor_global);
//                     if (owner == rank) {
//                         #pragma omp critical(sssp_update)
//                         {
//                             if (potential_dist_n < lp.vertices[neighbor_local_idx].distance) {
//                                 PRINT_DEBUG(3, "[Thread %d] Relax %ld->%ld: Dist[%ld] %d -> %d", thread_id, static_cast<long>(v_global), static_cast<long>(neighbor_global), static_cast<long>(neighbor_global), lp.vertices[neighbor_local_idx].distance, potential_dist_n);
//                                 lp.vertices[neighbor_local_idx].distance = potential_dist_n;
//                                 lp.vertices[neighbor_local_idx].parent_global = v_global;
//                                 lp.vertices[neighbor_local_idx].affected = true;
//                                 local_change_this_phase = true;
//                             }
//                         }
//                     }
//                 }
//             }

//             if (current_dist_n != INF) {
//                 int potential_dist_v = (current_dist_n == INF || weight == INF) ? INF : current_dist_n + weight;
//                 if (potential_dist_v < current_dist_v) {
//                     #pragma omp critical(sssp_update)
//                     {
//                         if (potential_dist_v < lp.vertices[i].distance) {
//                             PRINT_DEBUG(3, "[Thread %d] Relax %ld->%ld: Dist[%ld] %d -> %d", thread_id, static_cast<long>(neighbor_global), static_cast<long>(v_global), static_cast<long>(v_global), lp.vertices[i].distance, potential_dist_v);
//                             lp.vertices[i].distance = potential_dist_v;
//                             lp.vertices[i].parent_global = neighbor_global;
//                             lp.vertices[i].affected = true;
//                             local_change_this_phase = true;
//                             lp.vertices[i].just_updated = true;
//                         }
//                     }
//                 }
//             }
//         }
//     }

//     for (idx_t i = 0; i < lp.total_vertices; ++i) {
//         if (lp.vertices[i].just_updated && !lp.vertices[i].is_ghost) {
//             for (int r = 0; r < size; r++) {
//                 if (r != rank) {
//                     MPI_Send(&lp.vertices[i].global_id, 1, MPI_LONG_LONG_INT, r, MPI_TAG_GHOST_UPDATE, MPI_COMM_WORLD);
//                     MPI_Send(&lp.vertices[i].distance, 1, MPI_INT, r, MPI_TAG_GHOST_UPDATE, MPI_COMM_WORLD);
//                     MPI_Send(&lp.vertices[i].parent_global, 1, MPI_LONG_LONG_INT, r, MPI_TAG_GHOST_UPDATE, MPI_COMM_WORLD);
//                 }
//             }
//             lp.vertices[i].just_updated = false;
//         }
//     }

//     MPI_Iprobe(MPI_ANY_SOURCE, MPI_TAG_GHOST_UPDATE, MPI_COMM_WORLD, &flag, &status);
//     while (flag) {
//         idx_t global_id;
//         int distance;
//         idx_t parent_global;
//         MPI_Recv(&global_id, 1, MPI_LONG_LONG_INT, status.MPI_SOURCE, MPI_TAG_GHOST_UPDATE, MPI_COMM_WORLD, &status);
//         MPI_Recv(&distance, 1, MPI_INT, status.MPI_SOURCE, MPI_TAG_GHOST_UPDATE, MPI_COMM_WORLD, &status);
//         MPI_Recv(&parent_global, 1, MPI_LONG_LONG_INT, status.MPI_SOURCE, MPI_TAG_GHOST_UPDATE, MPI_COMM_WORLD, &status);
//         idx_t local_idx = find_local_index(lp, global_id);
//         if (local_idx != -1 && lp.vertices[local_idx].is_ghost) {
//             lp.vertices[local_idx].distance = distance;
//             lp.vertices[local_idx].parent_global = parent_global;
//             lp.vertices[local_idx].affected = true;
//         }
//         MPI_Iprobe(MPI_ANY_SOURCE, MPI_TAG_GHOST_UPDATE, MPI_COMM_WORLD, &flag, &status);
//     }

//     any_local_change_made = local_change_this_phase;
//     PRINT_DEBUG(1, "Relaxation phase finished. Local change detected: %s", any_local_change_made ? "YES" : "NO");
// }

// void gather_and_print_results(const LocalGraphPartition& lp, idx_t num_global_vertices, int rank, int size) {
//     PRINT_DEBUG(1, "Gathering final results to Rank 0.");

//     std::vector<long> send_buffer(lp.num_local_vertices * 3);
//     int send_count = 0;
//     for (idx_t i = 0; i < lp.total_vertices; ++i) {
//         if (!lp.vertices[i].is_ghost) {
//             send_buffer[send_count++] = static_cast<long>(lp.vertices[i].global_id);
//             send_buffer[send_count++] = static_cast<long>(lp.vertices[i].distance);
//             send_buffer[send_count++] = static_cast<long>(lp.vertices[i].parent_global);
//         }
//     }

//     int items_per_vertex = 3;
//     int send_items = lp.num_local_vertices * items_per_vertex;
//     std::vector<int> recv_counts(size);
//     MPI_Gather(&send_items, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, ROOT_RANK, MPI_COMM_WORLD);

//     std::vector<long> recv_buffer;
//     std::vector<int> displs;
//     int total_items_received = 0;

//     if (rank == ROOT_RANK) {
//         displs.resize(size);
//         displs[0] = 0;
//         total_items_received = recv_counts[0];
//         for (int i = 1; i < size; ++i) {
//             displs[i] = displs[i - 1] + recv_counts[i - 1];
//             total_items_received += recv_counts[i];
//         }
//         recv_buffer.resize(total_items_received);
//         PRINT_DEBUG(2, "Root receiving %d total items for results.", total_items_received);
//     }

//     MPI_Gatherv(send_buffer.data(), send_items, MPI_LONG, recv_buffer.data(), recv_counts.data(), displs.data(), MPI_LONG, ROOT_RANK, MPI_COMM_WORLD);

//     if (rank == ROOT_RANK) {
//         std::cout << "\n--- Final SSSP Results (gathered at Rank 0) ---\n";
//         std::cout << "GlobalID\tDistance\tParentID\n";
//         std::cout << "------------------------------------------------\n";
//         for (int i = 0; i < total_items_received; i += 3) {
//             long gid = recv_buffer[i];
//             long dist = recv_buffer[i + 1];
//             long parent = recv_buffer[i + 2];
//             std::cout << gid << "\t\t" << (dist == INF ? -1 : dist) << "\t\t" << parent << "\n";
//         }
//         std::cout << "------------------------------------------------\n";
//     }
//     PRINT_DEBUG(1, "Result gathering complete.");
// }

// int main(int argc, char *argv[]) {
//     MPI_Init(&argc, &argv);
//     MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
//     int world_size;
//     MPI_Comm_size(MPI_COMM_WORLD, &world_size);

//     if (argc < 4) {
//         if (world_rank == ROOT_RANK) {
//             std::cerr << "Usage: mpirun -np <num_procs> " << argv[0] << " <graph_edge_file> <changes_file> <source_vertex_id>" << std::endl;
//         }
//         MPI_Finalize();
//         return 1;
//     }

//     std::string graph_filename = argv[1];
//     std::string changes_filename = argv[2];
//     idx_t source_global_id = std::atoi(argv[3]);

//     PRINT_DEBUG(1, "Application started with %d processes.", world_size);

//     GraphCSR csr;
//     std::vector<idx_t> partition_assignment;
//     std::vector<Change> changes;
//     LocalGraphPartition local_partition;
//     double start_time, end_time, setup_time, update_time;

//     start_time = MPI_Wtime();

//     if (world_rank == ROOT_RANK) {
//         PRINT_DEBUG(1, "Reading graph from %s", graph_filename.c_str());
//         idx_t num_vertices, num_edges;
//         std::vector<idx_t> u, v;
//         std::vector<int> w;
//         if (!read_graph_edges(graph_filename, num_vertices, num_edges, u, v, w)) {
//             std::cerr << "[Rank " << world_rank << "] Error reading graph file." << std::endl;
//             MPI_Abort(MPI_COMM_WORLD, 1);
//         }
//         PRINT_DEBUG(1, "Graph has %ld vertices and %ld directed edges.", static_cast<long>(num_vertices), static_cast<long>(num_edges));

//         PRINT_DEBUG(1, "Converting graph to CSR format...");
//         convert_edges_to_csr(num_vertices, num_edges, u, v, w, csr);

//         PRINT_DEBUG(1, "Partitioning graph using METIS for %d partitions...", world_size);
//         partition_graph(csr, world_size, partition_assignment);
//         PRINT_DEBUG(1, "METIS partitioning complete.");

//         PRINT_DEBUG(1, "Reading changes from %s", changes_filename.c_str());
//         if (!read_changes(changes_filename, changes)) {
//             std::cerr << "[Rank " << world_rank << "] Error reading changes file." << std::endl;
//             MPI_Abort(MPI_COMM_WORLD, 1);
//         }
//         PRINT_DEBUG(1, "Read %zu change operations.", changes.size());
//     }

//     distribute_graph(csr, partition_assignment, world_rank, world_size, local_partition);
//     PRINT_DEBUG(1, "Graph distribution complete. Local partition has %ld local and %ld ghost vertices (total %ld).",
//                 static_cast<long>(local_partition.num_local_vertices), static_cast<long>(local_partition.num_ghost_vertices), static_cast<long>(local_partition.total_vertices));

//     int num_changes = changes.size();
//     MPI_Bcast(&num_changes, 1, MPI_INT, ROOT_RANK, MPI_COMM_WORLD);
//     PRINT_DEBUG(2, "Received broadcast: num_changes = %d", num_changes);
//     if (world_rank != ROOT_RANK) {
//         changes.resize(num_changes);
//     }
//     MPI_Bcast(changes.data(), num_changes * sizeof(Change), MPI_BYTE, ROOT_RANK, MPI_COMM_WORLD);
//     PRINT_DEBUG(2, "Received broadcast of changes data.");

//     if (world_rank == ROOT_RANK) {
//         cleanup_csr(csr);
//         partition_assignment.clear();
//     }

//     setup_time = MPI_Wtime() - start_time;
//     PRINT_DEBUG(1, "Setup phase completed in %.4f seconds.", setup_time);

//     initialize_sssp(local_partition, source_global_id);
//     PRINT_DEBUG(1, "Initial SSSP distances set (source = %ld).", static_cast<long>(source_global_id));

//     MPI_Barrier(MPI_COMM_WORLD);
//     start_time = MPI_Wtime();

//     PRINT_DEBUG(1, "*** Starting SSSP Update Step 1: Applying Changes ***");
//     apply_changes_step1(changes, local_partition, world_rank, world_size);
//     PRINT_DEBUG(1, "*** Finished SSSP Update Step 1 ***");

//     MPI_Barrier(MPI_COMM_WORLD);

//     PRINT_DEBUG(1, "*** Starting SSSP Update Step 2: Iterative Update ***");
//     update_distances_step2(local_partition, world_rank, world_size);
//     PRINT_DEBUG(1, "*** Finished SSSP Update Step 2 ***");

//     update_time = MPI_Wtime() - start_time;
//     PRINT_DEBUG(1, "SSSP Update phase completed in %.4f seconds.", update_time);

//     MPI_Barrier(MPI_COMM_WORLD);

//     idx_t final_num_global_vertices = local_partition.num_global_vertices;
//     gather_and_print_results(local_partition, final_num_global_vertices, world_rank, world_size);

//     PRINT_DEBUG(1, "Cleaning up resources...");
//     cleanup_local_partition(local_partition);

//     end_time = MPI_Wtime();
//     if (world_rank == ROOT_RANK) {
//         std::cout << "\n----------------------------------------\n";
//         std::cout << "Total execution time: " << (end_time - (start_time - update_time)) << " seconds\n";
//         std::cout << "  Setup time:         " << setup_time << " seconds\n";
//         std::cout << "  Update time:        " << update_time << " seconds\n";
//         std::cout << "----------------------------------------\n";
//     }

//     MPI_Finalize();
//     return 0;
// }
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
#include <metis.h>

// Forward declarations to resolve "not declared in scope" errors
struct LocalGraphPartition;
idx_t find_local_index(const LocalGraphPartition& lp, idx_t global_id);
void propagate_deletions_distributed(LocalGraphPartition& lp, int rank, int size);
void relax_edges_distributed(LocalGraphPartition& lp, int rank, int size, bool& any_local_change_made);

#define INF std::numeric_limits<int>::max()
#define ROOT_RANK 0
#define MPI_TAG_GHOST_DIST_REQUEST 1001
#define MPI_TAG_GHOST_DIST_REPLY   1002
#define MPI_TAG_GHOST_UPDATE       1003
#define MPI_TAG_DEL_PROPAGATE      1004

#define DEBUG_LEVEL 2

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

using idx_t = int;  // Assuming idx_t is int, adjust if necessary

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
};

struct LocalGraphPartition {
    idx_t num_local_vertices;
    idx_t num_ghost_vertices;
    idx_t total_vertices;
    std::vector<LocalVertex> vertices;
    std::vector<idx_t> local_to_global_id;
    idx_t num_global_vertices;
    std::vector<idx_t> partition_assignment;
};

struct GraphCSR {
    idx_t nvtxs;
    idx_t nedges;
    std::vector<idx_t> xadj;
    std::vector<idx_t> adjncy;
    std::vector<idx_t> adjwgt;
};

bool read_graph_edges(const std::string& filename, idx_t& num_vertices, idx_t& num_edges, std::vector<idx_t>& u, std::vector<idx_t>& v, std::vector<int>& w) {
    std::ifstream f(filename);
    if (!f.is_open()) {
        std::cerr << "Error opening graph file: " << filename << std::endl;
        return false;
    }
    PRINT_DEBUG(2, "Reading edges from %s", filename.c_str());

    idx_t u_i, v_i, max_id = 0;
    int w_i;
    idx_t edge_count = 0;
    std::string line;
    std::vector<idx_t> u_temp, v_temp;
    std::vector<int> w_temp;

    while (std::getline(f, line)) {
        std::istringstream iss(line);
        if (iss >> u_i >> v_i >> w_i) {
            edge_count++;
            u_temp.push_back(u_i);
            v_temp.push_back(v_i);
            w_temp.push_back(w_i);
            if (u_i > max_id) max_id = u_i;
            if (v_i > max_id) max_id = v_i;
        }
    }

    num_vertices = max_id + 1;
    num_edges = edge_count;
    u = std::move(u_temp);
    v = std::move(v_temp);
    w = std::move(w_temp);

    f.close();
    PRINT_DEBUG(2, "Read %ld edges, max vertex ID %ld", static_cast<long>(edge_count), static_cast<long>(max_id));
    return true;
}

void convert_edges_to_csr(idx_t nvtxs, idx_t nedges, const std::vector<idx_t>& u, const std::vector<idx_t>& v, const std::vector<int>& w, GraphCSR& csr) {
    PRINT_DEBUG(2, "Converting %ld vertices, %ld edges to CSR", static_cast<long>(nvtxs), static_cast<long>(nedges));
    csr.nvtxs = nvtxs;
    csr.nedges = nedges;
    csr.xadj.resize(nvtxs + 1, 0);
    csr.adjncy.resize(nedges);
    csr.adjwgt.resize(nedges);

    for (idx_t i = 0; i < nedges; ++i) {
        csr.xadj[u[i]]++;
    }

    idx_t sum = 0;
    for (idx_t i = 0; i < nvtxs; ++i) {
        idx_t temp = csr.xadj[i];
        csr.xadj[i] = sum;
        sum += temp;
    }
    csr.xadj[nvtxs] = sum;

    std::vector<idx_t> pos(nvtxs, 0);
    for (idx_t i = 0; i < nedges; ++i) {
        idx_t src = u[i];
        idx_t edge_idx = csr.xadj[src] + pos[src];
        csr.adjncy[edge_idx] = v[i];
        csr.adjwgt[edge_idx] = w[i];
        pos[src]++;
    }

    PRINT_DEBUG(2, "CSR conversion complete.");
}

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
    while (std::getline(f, line)) {
        std::istringstream iss(line);
        if (iss >> type) {
            Change ch;
            ch.type = type;
            if (type == "INS" && iss >> u >> v >> w) {
                ch.u = u;
                ch.v = v;
                ch.weight = w;
            } else if (type == "DEL" && iss >> u >> v) {
                ch.u = u;
                ch.v = v;
                ch.weight = 0;
            } else {
                std::cerr << "Unknown or malformed change: " << line << std::endl;
                continue;
            }
            changes.push_back(ch);
        }
    }

    f.close();
    PRINT_DEBUG(2, "Read %zu changes.", changes.size());
    return true;
}
void partition_graph(const GraphCSR& csr, int num_partitions, std::vector<idx_t>& partition_result) {
    if (num_partitions <= 1) {
        PRINT_DEBUG(1, "Number of partitions is 1, skipping METIS call.");
        partition_result.resize(csr.nvtxs, 0);
        return;
    }

    PRINT_DEBUG(2, "Calling METIS_PartGraphKway...");
    partition_result.resize(csr.nvtxs);
    idx_t ncon = 1;
    idx_t objval;
    idx_t options[METIS_NOPTIONS];
    // Explicitly initialize options to zero before setting
    std::fill_n(options, METIS_NOPTIONS, 0);
    int set_options_ret = METIS_SetDefaultOptions(options);
    if (set_options_ret != METIS_OK) {
        std::cerr << "[Rank " << world_rank << "] METIS_SetDefaultOptions failed with error: " << set_options_ret << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    options[METIS_OPTION_DBGLVL] = 0;
    options[METIS_OPTION_CONTIG] = 1;

    // Create mutable copies of all required data
    idx_t nvtxs_copy = csr.nvtxs;
    idx_t ncon_copy = ncon;
    std::vector<idx_t> xadj_copy = csr.xadj;
    std::vector<idx_t> adjncy_copy = csr.adjncy;
    std::vector<idx_t> adjwgt_copy = csr.adjwgt;
    idx_t num_partitions_copy = num_partitions;

    int ret = METIS_PartGraphKway(&nvtxs_copy, &ncon_copy, xadj_copy.data(), adjncy_copy.data(),
                                  adjwgt_copy.data(), nullptr, adjwgt_copy.data(),
                                  &num_partitions_copy, nullptr, nullptr, options, &objval, partition_result.data());

    if (ret != METIS_OK) {
        std::cerr << "[Rank " << world_rank << "] METIS_PartGraphKway failed with error code: " << ret << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    PRINT_DEBUG(1, "METIS partitioning successful. Edge cut: %ld", static_cast<long>(objval));
}

void distribute_graph(const GraphCSR& csr, const std::vector<idx_t>& partition_assignment, int rank, int size, LocalGraphPartition& local_partition) {
    PRINT_DEBUG(2, "Starting graph distribution process.");

    std::vector<idx_t> local_counts(size, 0);
    std::vector<idx_t> ghost_counts(size, 0);
    idx_t num_global_vertices = 0;

    if (rank == ROOT_RANK) {
        num_global_vertices = csr.nvtxs;
        std::vector<bool> is_ghost_for_rank(csr.nvtxs * size, false);
        std::vector<idx_t> ghost_node_indices_count(size, 0);

        for (idx_t i = 0; i < csr.nvtxs; ++i) {
            local_counts[partition_assignment[i]]++;
        }

        for (idx_t u_global = 0; u_global < csr.nvtxs; ++u_global) {
            int owner_rank = partition_assignment[u_global];
            for (idx_t edge_idx = csr.xadj[u_global]; edge_idx < csr.xadj[u_global + 1]; ++edge_idx) {
                idx_t v_global = csr.adjncy[edge_idx];
                int target_rank = partition_assignment[v_global];
                if (owner_rank != target_rank) {
                    if (!is_ghost_for_rank[v_global * size + owner_rank]) {
                        is_ghost_for_rank[v_global * size + owner_rank] = true;
                        ghost_node_indices_count[owner_rank]++;
                    }
                    if (!is_ghost_for_rank[u_global * size + target_rank]) {
                        is_ghost_for_rank[u_global * size + target_rank] = true;
                        ghost_node_indices_count[target_rank]++;
                    }
                }
            }
        }
        for (int r = 0; r < size; ++r) ghost_counts[r] = ghost_node_indices_count[r];
        PRINT_DEBUG(2, "Calculated local/ghost counts on root.");
    }

    MPI_Bcast(&num_global_vertices, 1, MPI_LONG_LONG_INT, ROOT_RANK, MPI_COMM_WORLD);
    MPI_Bcast(local_counts.data(), size, MPI_LONG_LONG_INT, ROOT_RANK, MPI_COMM_WORLD);
    MPI_Bcast(ghost_counts.data(), size, MPI_LONG_LONG_INT, ROOT_RANK, MPI_COMM_WORLD);

    local_partition.num_local_vertices = local_counts[rank];
    local_partition.num_ghost_vertices = ghost_counts[rank];
    local_partition.total_vertices = local_partition.num_local_vertices + local_partition.num_ghost_vertices;
    local_partition.num_global_vertices = num_global_vertices;

    if (rank == ROOT_RANK) {
        MPI_Bcast(const_cast<idx_t*>(partition_assignment.data()), num_global_vertices, MPI_LONG_LONG_INT, ROOT_RANK, MPI_COMM_WORLD);
    } else {
        local_partition.partition_assignment.resize(num_global_vertices);
        MPI_Bcast(local_partition.partition_assignment.data(), num_global_vertices, MPI_LONG_LONG_INT, ROOT_RANK, MPI_COMM_WORLD);
    }

    PRINT_DEBUG(2, "Received counts: local=%ld, ghost=%ld, total=%ld", static_cast<long>(local_counts[rank]), static_cast<long>(ghost_counts[rank]), static_cast<long>(local_partition.total_vertices));

    local_partition.vertices.resize(local_partition.total_vertices);
    local_partition.local_to_global_id.resize(local_partition.total_vertices);

    if (rank == ROOT_RANK) {
        idx_t local_idx = 0;
        std::set<idx_t> ghost_set;

        for (idx_t v_global = 0; v_global < csr.nvtxs; v_global++) {
            if (partition_assignment[v_global] == rank) {
                local_partition.local_to_global_id[local_idx] = v_global;
                local_partition.vertices[local_idx].global_id = v_global;
                local_partition.vertices[local_idx].partition_id = rank;
                local_partition.vertices[local_idx].is_ghost = false;
                idx_t num_edges = csr.xadj[v_global + 1] - csr.xadj[v_global];
                local_partition.vertices[local_idx].edges.resize(num_edges);
                idx_t start = csr.xadj[v_global];
                for (idx_t e = 0; e < num_edges; e++) {
                    local_partition.vertices[local_idx].edges[e].target_global = csr.adjncy[start + e];
                    local_partition.vertices[local_idx].edges[e].weight = csr.adjwgt[start + e];
                    if (partition_assignment[csr.adjncy[start + e]] != rank) {
                        ghost_set.insert(csr.adjncy[start + e]);
                    }
                }
                local_idx++;
            }
        }

        idx_t ghost_idx = local_partition.num_local_vertices;
        for (const auto& ghost : ghost_set) {
            local_partition.local_to_global_id[ghost_idx] = ghost;
            local_partition.vertices[ghost_idx].global_id = ghost;
            local_partition.vertices[ghost_idx].partition_id = partition_assignment[ghost];
            local_partition.vertices[ghost_idx].is_ghost = true;
            ghost_idx++;
        }

        for (int r = 1; r < size; r++) {
            idx_t num_local_r = local_counts[r];
            std::vector<idx_t> local_global_ids;
            for (idx_t i = 0; i < csr.nvtxs; i++) {
                if (partition_assignment[i] == r) {
                    local_global_ids.push_back(i);
                }
            }
            MPI_Send(local_global_ids.data(), num_local_r, MPI_LONG_LONG_INT, r, 0, MPI_COMM_WORLD);

            for (idx_t k = 0; k < num_local_r; k++) {
                idx_t v_global = local_global_ids[k];
                idx_t num_edges = csr.xadj[v_global + 1] - csr.xadj[v_global];
                MPI_Send(&num_edges, 1, MPI_LONG_LONG_INT, r, 0, MPI_COMM_WORLD);
                idx_t start = csr.xadj[v_global];
                MPI_Send(&csr.adjncy[start], num_edges, MPI_LONG_LONG_INT, r, 0, MPI_COMM_WORLD);
                MPI_Send(&csr.adjwgt[start], num_edges, MPI_INT, r, 0, MPI_COMM_WORLD);
            }

            std::set<idx_t> ghost_set_r;
            for (idx_t v_global = 0; v_global < csr.nvtxs; v_global++) {
                if (partition_assignment[v_global] == r) {
                    idx_t start = csr.xadj[v_global];
                    idx_t end = csr.xadj[v_global + 1];
                    for (idx_t idx = start; idx < end; idx++) {
                        idx_t target = csr.adjncy[idx];
                        if (partition_assignment[target] != r) {
                            ghost_set_r.insert(target);
                        }
                    }
                }
            }
            idx_t num_ghost_r = ghost_set_r.size();
            MPI_Send(&num_ghost_r, 1, MPI_LONG_LONG_INT, r, 0, MPI_COMM_WORLD);
            std::vector<idx_t> ghost_global_ids(ghost_set_r.begin(), ghost_set_r.end());
            MPI_Send(ghost_global_ids.data(), num_ghost_r, MPI_LONG_LONG_INT, r, 0, MPI_COMM_WORLD);
        }
    } else {
        std::vector<idx_t> local_global_ids(local_partition.num_local_vertices);
        MPI_Recv(local_global_ids.data(), local_partition.num_local_vertices, MPI_LONG_LONG_INT, ROOT_RANK, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (idx_t k = 0; k < local_partition.num_local_vertices; k++) {
            idx_t v_global = local_global_ids[k];
            local_partition.local_to_global_id[k] = v_global;
            local_partition.vertices[k].global_id = v_global;
            local_partition.vertices[k].partition_id = rank;
            local_partition.vertices[k].is_ghost = false;
            idx_t num_edges;
            MPI_Recv(&num_edges, 1, MPI_LONG_LONG_INT, ROOT_RANK, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            local_partition.vertices[k].edges.resize(num_edges);
            std::vector<idx_t> targets(num_edges);
            std::vector<int> weights(num_edges);
            MPI_Recv(targets.data(), num_edges, MPI_LONG_LONG_INT, ROOT_RANK, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(weights.data(), num_edges, MPI_INT, ROOT_RANK, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (idx_t e = 0; e < num_edges; e++) {
                local_partition.vertices[k].edges[e].target_global = targets[e];
                local_partition.vertices[k].edges[e].weight = weights[e];
            }
        }

        idx_t num_ghost_r;
        MPI_Recv(&num_ghost_r, 1, MPI_LONG_LONG_INT, ROOT_RANK, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        std::vector<idx_t> ghost_global_ids(num_ghost_r);
        MPI_Recv(ghost_global_ids.data(), num_ghost_r, MPI_LONG_LONG_INT, ROOT_RANK, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (idx_t k = 0; k < num_ghost_r; k++) {
            idx_t ghost_idx = local_partition.num_local_vertices + k;
            local_partition.local_to_global_id[ghost_idx] = ghost_global_ids[k];
            local_partition.vertices[ghost_idx].global_id = ghost_global_ids[k];
            local_partition.vertices[ghost_idx].partition_id = local_partition.partition_assignment[ghost_global_ids[k]];
            local_partition.vertices[ghost_idx].is_ghost = true;
        }
    }

    PRINT_DEBUG(1, "Graph distribution completed successfully.");
}

void cleanup_csr(GraphCSR& csr) {
    csr.xadj.clear();
    csr.adjncy.clear();
    csr.adjwgt.clear();
    PRINT_DEBUG(2, "Cleaned up global CSR data.");
}

void cleanup_local_partition(LocalGraphPartition& lp) {
    lp.vertices.clear();
    lp.local_to_global_id.clear();
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

void apply_changes_step1(const std::vector<Change>& changes, LocalGraphPartition& lp, int rank, int size) {
    PRINT_DEBUG(1, "Applying %zu changes (Step 1).", changes.size());

    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < changes.size(); ++i) {
        const Change& ch = changes[i];
        idx_t u_global = ch.u;
        idx_t v_global = ch.v;
        int thread_id = omp_get_thread_num();
        PRINT_DEBUG(3, "[Thread %d] Processing change %zu: %s %ld %ld", thread_id, i, ch.type.c_str(), static_cast<long>(u_global), static_cast<long>(v_global));

        idx_t u_local = find_local_index(lp, u_global);
        idx_t v_local = find_local_index(lp, v_global);

        if (ch.type == "DEL") {
            bool changed = false;
            if (v_local != -1 && lp.vertices[v_local].parent_global == u_global) {
                #pragma omp critical(sssp_update)
                {
                    if (lp.vertices[v_local].distance != INF) {
                        PRINT_DEBUG(2, "[Thread %d] DEL %ld->%ld detected as SSSP edge. Disconnecting %ld.", thread_id, static_cast<long>(u_global), static_cast<long>(v_global), static_cast<long>(v_global));
                        lp.vertices[v_local].distance = INF;
                        lp.vertices[v_local].parent_global = -1;
                        lp.vertices[v_local].affected_del = true;
                        lp.vertices[v_local].affected = true;
                        changed = true;
                    }
                }
            } else if (u_local != -1 && lp.vertices[u_local].parent_global == v_global) {
                #pragma omp critical(sssp_update)
                {
                    if (lp.vertices[u_local].distance != INF) {
                        PRINT_DEBUG(2, "[Thread %d] DEL %ld->%ld detected as SSSP edge. Disconnecting %ld.", thread_id, static_cast<long>(v_global), static_cast<long>(u_global), static_cast<long>(u_global));
                        lp.vertices[u_local].distance = INF;
                        lp.vertices[u_local].parent_global = -1;
                        lp.vertices[u_local].affected_del = true;
                        lp.vertices[u_local].affected = true;
                        changed = true;
                    }
                }
            }
        } else if (ch.type == "INS") {
            int weight = ch.weight;
            int dist_u = (u_local != -1) ? lp.vertices[u_local].distance : INF;
            int dist_v = (v_local != -1) ? lp.vertices[v_local].distance : INF;

            if (v_local != -1 && dist_u != INF) {
                int potential_dist_v = (dist_u == INF || weight == INF) ? INF : dist_u + weight;
                if (potential_dist_v < dist_v) {
                    #pragma omp critical(sssp_update)
                    {
                        if (potential_dist_v < lp.vertices[v_local].distance) {
                            PRINT_DEBUG(2, "[Thread %d] INS %ld->%ld (w=%d) updates Dist[%ld]: %d -> %d",
                                        thread_id, static_cast<long>(u_global), static_cast<long>(v_global), weight, static_cast<long>(v_global), lp.vertices[v_local].distance, potential_dist_v);
                            lp.vertices[v_local].distance = potential_dist_v;
                            lp.vertices[v_local].parent_global = u_global;
                            lp.vertices[v_local].affected = true;
                        }
                    }
                }
            }

            if (u_local != -1 && dist_v != INF) {
                int potential_dist_u = (dist_v == INF || weight == INF) ? INF : dist_v + weight;
                if (potential_dist_u < dist_u) {
                    #pragma omp critical(sssp_update)
                    {
                        if (potential_dist_u < lp.vertices[u_local].distance) {
                            PRINT_DEBUG(2, "[Thread %d] WEL %ld->%ld (w=%d) updates Dist[%ld]: %d -> %d",
                                        thread_id, static_cast<long>(v_global), static_cast<long>(u_global), weight, static_cast<long>(u_global), lp.vertices[u_local].distance, potential_dist_u);
                            lp.vertices[u_local].distance = potential_dist_u;
                            lp.vertices[u_local].parent_global = v_global;
                            lp.vertices[u_local].affected = true;
                        }
                    }
                }
            }
        }
    }
    PRINT_DEBUG(1, "Finished processing changes locally in Step 1.");
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

    while (global_change_occurred_in_iteration) {
        iteration++;
        PRINT_DEBUG(1, "--- Relaxation Iteration %d ---", iteration);
        bool local_change_occurred = false;

        relax_edges_distributed(lp, rank, size, local_change_occurred);

        MPI_Allreduce(&local_change_occurred, &global_change_occurred_in_iteration, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);

        if (rank == ROOT_RANK) {
            PRINT_DEBUG(1, "Global change in iteration %d: %s", iteration, global_change_occurred_in_iteration ? "YES" : "NO");
        }
        if (iteration > lp.total_vertices * 2 && global_change_occurred_in_iteration) {
            if (rank == ROOT_RANK) std::cerr << "Warning: Potential non-convergence after " << iteration << " iterations." << std::endl;
            break;
        }
    }
    PRINT_DEBUG(1, "Step 2b: Iterative relaxation finished after %d iterations.", iteration);
}

idx_t find_local_index(const LocalGraphPartition& lp, idx_t global_id) {
    for (idx_t i = 0; i < lp.total_vertices; ++i) {
        if (lp.vertices[i].global_id == global_id) {
            return i;
        }
    }
    return -1;
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
                        MPI_Send(&c_global, 1, MPI_LONG_LONG_INT, owner, MPI_TAG_DEL_PROPAGATE, MPI_COMM_WORLD);
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
            MPI_Recv(&c_global, 1, MPI_LONG_LONG_INT, status.MPI_SOURCE, MPI_TAG_DEL_PROPAGATE, MPI_COMM_WORLD, &status);
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

    for (const auto& ghost_global : requested_ghosts) {
        int owner = find_owner_rank(lp, ghost_global);
        MPI_Send(&ghost_global, 1, MPI_LONG_LONG_INT, owner, MPI_TAG_GHOST_DIST_REQUEST, MPI_COMM_WORLD);
    }

    MPI_Status status;
    int flag;
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_TAG_GHOST_DIST_REQUEST, MPI_COMM_WORLD, &flag, &status);
    while (flag) {
        idx_t requested_global;
        MPI_Recv(&requested_global, 1, MPI_LONG_LONG_INT, status.MPI_SOURCE, MPI_TAG_GHOST_DIST_REQUEST, MPI_COMM_WORLD, &status);
        idx_t local_idx = find_local_index(lp, requested_global);
        if (local_idx != -1) {
            int distance = lp.vertices[local_idx].distance;
            MPI_Send(&distance, 1, MPI_INT, status.MPI_SOURCE, MPI_TAG_GHOST_DIST_REPLY, MPI_COMM_WORLD);
        }
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_TAG_GHOST_DIST_REQUEST, MPI_COMM_WORLD, &flag, &status);
    }

    std::map<idx_t, int> ghost_distances;
    for (const auto& ghost_global : requested_ghosts) {
        int distance;
        MPI_Recv(&distance, 1, MPI_INT, find_owner_rank(lp, ghost_global), MPI_TAG_GHOST_DIST_REPLY, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        ghost_distances[ghost_global] = distance;
    }

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
                    if (owner == rank) {
                        #pragma omp critical(sssp_update)
                        {
                            if (potential_dist_n < lp.vertices[neighbor_local_idx].distance) {
                                PRINT_DEBUG(3, "[Thread %d] Relax %ld->%ld: Dist[%ld] %d -> %d", thread_id, static_cast<long>(v_global), static_cast<long>(neighbor_global), static_cast<long>(neighbor_global), lp.vertices[neighbor_local_idx].distance, potential_dist_n);
                                lp.vertices[neighbor_local_idx].distance = potential_dist_n;
                                lp.vertices[neighbor_local_idx].parent_global = v_global;
                                lp.vertices[neighbor_local_idx].affected = true;
                                local_change_this_phase = true;
                            }
                        }
                    }
                }
            }

            if (current_dist_n != INF) {
                int potential_dist_v = (current_dist_n == INF || weight == INF) ? INF : current_dist_n + weight;
                if (potential_dist_v < current_dist_v) {
                    #pragma omp critical(sssp_update)
                    {
                        if (potential_dist_v < lp.vertices[i].distance) {
                            PRINT_DEBUG(3, "[Thread %d] Relax %ld->%ld: Dist[%ld] %d -> %d", thread_id, static_cast<long>(neighbor_global), static_cast<long>(v_global), static_cast<long>(v_global), lp.vertices[i].distance, potential_dist_v);
                            lp.vertices[i].distance = potential_dist_v;
                            lp.vertices[i].parent_global = neighbor_global;
                            lp.vertices[i].affected = true;
                            local_change_this_phase = true;
                            lp.vertices[i].just_updated = true;
                        }
                    }
                }
            }
        }
    }

    for (idx_t i = 0; i < lp.total_vertices; ++i) {
        if (lp.vertices[i].just_updated && !lp.vertices[i].is_ghost) {
            for (int r = 0; r < size; r++) {
                if (r != rank) {
                    MPI_Send(&lp.vertices[i].global_id, 1, MPI_LONG_LONG_INT, r, MPI_TAG_GHOST_UPDATE, MPI_COMM_WORLD);
                    MPI_Send(&lp.vertices[i].distance, 1, MPI_INT, r, MPI_TAG_GHOST_UPDATE, MPI_COMM_WORLD);
                    MPI_Send(&lp.vertices[i].parent_global, 1, MPI_LONG_LONG_INT, r, MPI_TAG_GHOST_UPDATE, MPI_COMM_WORLD);
                }
            }
            lp.vertices[i].just_updated = false;
        }
    }

    MPI_Iprobe(MPI_ANY_SOURCE, MPI_TAG_GHOST_UPDATE, MPI_COMM_WORLD, &flag, &status);
    while (flag) {
        idx_t global_id;
        int distance;
        idx_t parent_global;
        MPI_Recv(&global_id, 1, MPI_LONG_LONG_INT, status.MPI_SOURCE, MPI_TAG_GHOST_UPDATE, MPI_COMM_WORLD, &status);
        MPI_Recv(&distance, 1, MPI_INT, status.MPI_SOURCE, MPI_TAG_GHOST_UPDATE, MPI_COMM_WORLD, &status);
        MPI_Recv(&parent_global, 1, MPI_LONG_LONG_INT, status.MPI_SOURCE, MPI_TAG_GHOST_UPDATE, MPI_COMM_WORLD, &status);
        idx_t local_idx = find_local_index(lp, global_id);
        if (local_idx != -1 && lp.vertices[local_idx].is_ghost) {
            lp.vertices[local_idx].distance = distance;
            lp.vertices[local_idx].parent_global = parent_global;
            lp.vertices[local_idx].affected = true;
        }
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_TAG_GHOST_UPDATE, MPI_COMM_WORLD, &flag, &status);
    }

    any_local_change_made = local_change_this_phase;
    PRINT_DEBUG(1, "Relaxation phase finished. Local change detected: %s", any_local_change_made ? "YES" : "NO");
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

    GraphCSR csr;
    std::vector<idx_t> partition_assignment;
    std::vector<Change> changes;
    LocalGraphPartition local_partition;
    double start_time, end_time, setup_time, update_time;

    start_time = MPI_Wtime();

    if (world_rank == ROOT_RANK) {
        PRINT_DEBUG(1, "Reading graph from %s", graph_filename.c_str());
        idx_t num_vertices, num_edges;
        std::vector<idx_t> u, v;
        std::vector<int> w;
        if (!read_graph_edges(graph_filename, num_vertices, num_edges, u, v, w)) {
            std::cerr << "[Rank " << world_rank << "] Error reading graph file." << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        PRINT_DEBUG(1, "Graph has %ld vertices and %ld directed edges.", static_cast<long>(num_vertices), static_cast<long>(num_edges));

        PRINT_DEBUG(1, "Converting graph to CSR format...");
        convert_edges_to_csr(num_vertices, num_edges, u, v, w, csr);

        PRINT_DEBUG(1, "Partitioning graph using METIS for %d partitions...", world_size);
        partition_graph(csr, world_size, partition_assignment);
        PRINT_DEBUG(1, "METIS partitioning complete.");

        PRINT_DEBUG(1, "Reading changes from %s", changes_filename.c_str());
        if (!read_changes(changes_filename, changes)) {
            std::cerr << "[Rank " << world_rank << "] Error reading changes file." << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        PRINT_DEBUG(1, "Read %zu change operations.", changes.size());
    }

    distribute_graph(csr, partition_assignment, world_rank, world_size, local_partition);
    PRINT_DEBUG(1, "Graph distribution complete. Local partition has %ld local and %ld ghost vertices (total %ld).",
                static_cast<long>(local_partition.num_local_vertices), static_cast<long>(local_partition.num_ghost_vertices), static_cast<long>(local_partition.total_vertices));

    int num_changes = changes.size();
    MPI_Bcast(&num_changes, 1, MPI_INT, ROOT_RANK, MPI_COMM_WORLD);
    PRINT_DEBUG(2, "Received broadcast: num_changes = %d", num_changes);
    if (world_rank != ROOT_RANK) {
        changes.resize(num_changes);
    }
    MPI_Bcast(changes.data(), num_changes * sizeof(Change), MPI_BYTE, ROOT_RANK, MPI_COMM_WORLD);
    PRINT_DEBUG(2, "Received broadcast of changes data.");

    if (world_rank == ROOT_RANK) {
        cleanup_csr(csr);
        partition_assignment.clear();
    }

    setup_time = MPI_Wtime() - start_time;
    PRINT_DEBUG(1, "Setup phase completed in %.4f seconds.", setup_time);

    initialize_sssp(local_partition, source_global_id);
    PRINT_DEBUG(1, "Initial SSSP distances set (source = %ld).", static_cast<long>(source_global_id));

    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();

    PRINT_DEBUG(1, "*** Starting SSSP Update Step 1: Applying Changes ***");
    apply_changes_step1(changes, local_partition, world_rank, world_size);
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
        std::cout << "Total execution time: " << (end_time - (start_time - update_time)) << " seconds\n";
        std::cout << "  Setup time:         " << setup_time << " seconds\n";
        std::cout << "  Update time:        " << update_time << " seconds\n";
        std::cout << "----------------------------------------\n";
    }

    MPI_Finalize();
    return 0;
}