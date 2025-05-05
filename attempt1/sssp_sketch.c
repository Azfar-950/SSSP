// dynamic_sssp_mpi_omp.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <omp.h>
#include <metis.h>
#include <limits.h> // For INT_MAX as infinity
#include <stdbool.h> // For bool type

// --- Data Structures ---
typedef struct {
    int target;
    int weight;
} Edge;

typedef struct {
    int id; // Global ID
    // Local data
    int parent;
    int distance;
    bool affected;
    bool affected_del;
    // Edges originating from this node (only store local edges)
    Edge *edges;
    int num_edges;
    // Ghost node info (if needed for communication)
    bool is_ghost_for_rank[/*MAX_PROCS*/]; // Example
} VertexData;

typedef struct {
    int u;
    int v;
    int weight; // Only for INS
    char type[4]; // "INS" or "DEL"
} Change;

// Function Prototypes (Simplified)
void read_graph(const char* filename, int rank, int size, /*...*/);
void partition_and_distribute_graph(int num_global_vertices, /*...*/);
void apply_changes_step1(Change* changes, int num_changes, /*...*/);
void update_distances_step2(/*...*/);
void print_debug(int rank, const char* message);
void gather_results(/*...*/);


// --- Main ---
int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 4) {
        if (rank == 0) {
            fprintf(stderr, "Usage: mpirun -np <num_procs> %s <graph_file> <changes_file> <source_vertex>\n", argv[0]);
        }
        MPI_Finalize();
        return 1;
    }

    const char* graph_filename = argv[1];
    const char* changes_filename = argv[2];
    int source_vertex = atoi(argv[3]);

    // --- Variables ---
    int num_global_vertices = 0; // Master reads this
    int num_local_vertices = 0;
    VertexData *local_vertices = NULL; // Array for local + potentially ghost nodes
    Change *changes = NULL; // Master reads this
    int num_changes = 0;

    // --- Master (Rank 0) ---
    if (rank == 0) {
        print_debug(rank, "Reading graph and changes...");
        // read_full_graph(...) to get num_global_vertices
        // read_changes_file(...)

        print_debug(rank, "Partitioning graph using METIS...");
        // Prepare graph data in METIS format (adjacency structure: xadj, adjncy, adjwgt)
        // Call METIS_PartGraphKway(...) to get partition assignments (part array)

        print_debug(rank, "Distributing graph partitions...");
        // partition_and_distribute_graph(num_global_vertices, part, ...);
        // This involves complex MPI sends to slaves

        print_debug(rank, "Broadcasting changes...");
        // Broadcast num_changes
        // Broadcast changes array
    } else {
        // --- Slaves (Rank > 0) ---
        print_debug(rank, "Waiting for partition data...");
        // Receive local graph structure (num_local_vertices, edges, ghost info) from master

        print_debug(rank, "Waiting for changes...");
        // Receive num_changes
        // Receive changes array
    }

    // --- All Ranks: Initialization ---
    print_debug(rank, "Initializing local SSSP data...");
    // Allocate local_vertices if not done already
    // Initialize distances (0 for source if local, infinity otherwise), parents (-1), affected flags (false)
    #pragma omp parallel for
    for (int i = 0; i < num_local_vertices; ++i) {
        local_vertices[i].distance = (local_vertices[i].id == source_vertex) ? 0 : INT_MAX;
        local_vertices[i].parent = -1;
        local_vertices[i].affected = false;
        local_vertices[i].affected_del = false;
        // Potentially mark source vertex as affected if we start update immediately
    }


    // --- SSSP Update ---
    MPI_Barrier(MPI_COMM_WORLD); // Ensure setup is complete
    double start_time = MPI_Wtime();

    print_debug(rank, "Applying changes (Step 1)...");
    apply_changes_step1(changes, num_changes, local_vertices, num_local_vertices, rank, size /*, ghost info...*/);

    print_debug(rank, "Updating distances iteratively (Step 2)...");
    update_distances_step2(local_vertices, num_local_vertices, rank, size /*, ghost info...*/);

    MPI_Barrier(MPI_COMM_WORLD);
    double end_time = MPI_Wtime();

    if (rank == 0) {
         printf("----------------------------------------\n");
         printf("Update process took %f seconds.\n", end_time - start_time);
         printf("----------------------------------------\n");
    }

    // --- Gather and Print Results (Master) ---
    print_debug(rank, "Gathering results...");
    gather_results(local_vertices, num_local_vertices, num_global_vertices, rank, size);

    // --- Cleanup ---
    print_debug(rank, "Cleaning up...");
    // Free allocated memory (local_vertices, edges, changes if master, etc.)

    MPI_Finalize();
    return 0;
}

// --- Helper Functions Implementation (Highly Simplified) ---

void print_debug(int rank, const char* message) {
    // Get hostname (requires unistd.h, etc.) - basic version:
    char hostname[MPI_MAX_PROCESSOR_NAME];
    int len;
    MPI_Get_processor_name(hostname, &len);
    printf("[Rank %d on %s]: %s\n", rank, hostname, message);
    fflush(stdout); // Ensure message appears immediately
}

void apply_changes_step1(Change* changes, int num_changes, VertexData* local_verts, int n_local, int rank, int size /*,...*/) {
    print_debug(rank, "Starting Step 1: Identifying affected nodes from changes.");
    // Use OpenMP to parallelize processing of changes relevant to this rank
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < num_changes; ++i) {
        Change current_change = changes[i];
        int u_global = current_change.u;
        int v_global = current_change.v;
        int thread_id = omp_get_thread_num(); // For debug

        // Determine if u or v are local or ghost nodes for this rank
        // Find local indices for u and v if they exist locally
        int u_local_idx = -1; /* find_local_index(u_global); */
        int v_local_idx = -1; /* find_local_index(v_global); */

        // Simplified logic - needs proper mapping and ghost handling
        if (strcmp(current_change.type, "DEL") == 0) {
             // --- Handle Deletion (Adapting Algo 2, lines 4-9) ---
             // Check if edge (u,v) or (v,u) is in the *current* SSSP tree locally
             // Requires checking parent pointers
             int affected_node_local_idx = -1;
             if (u_local_idx != -1 && local_verts[u_local_idx].parent == v_global) {
                 affected_node_local_idx = u_local_idx;
             } else if (v_local_idx != -1 && local_verts[v_local_idx].parent == u_global) {
                 affected_node_local_idx = v_local_idx;
             }

             if (affected_node_local_idx != -1) {
                 #pragma omp critical // Protect concurrent updates if needed, though paper avoids locks
                 {
                     if (local_verts[affected_node_local_idx].distance != INT_MAX) { // Avoid redundant work
                        printf("[Rank %d, Thr %d] DEL affects local node %d (global %d)\n", rank, thread_id, affected_node_local_idx, local_verts[affected_node_local_idx].id);
                        local_verts[affected_node_local_idx].distance = INT_MAX;
                        local_verts[affected_node_local_idx].parent = -1; // Disconnect
                        local_verts[affected_node_local_idx].affected_del = true;
                        local_verts[affected_node_local_idx].affected = true;
                        // If affected_node is a ghost node for others, need MPI Send
                     }
                 }
             }
        } else { // INS
            // --- Handle Insertion (Adapting Algo 2, lines 10-19) ---
            int weight = current_change.weight;
            // Need potentially current distances of u AND v, may require MPI comms if one is external
            int dist_u = (u_local_idx != -1) ? local_verts[u_local_idx].distance : INT_MAX; // Fetch if external via MPI
            int dist_v = (v_local_idx != -1) ? local_verts[v_local_idx].distance : INT_MAX; // Fetch if external via MPI

            // Check u -> v
            if (v_local_idx != -1 && dist_u != INT_MAX && dist_v > dist_u + weight) {
                 #pragma omp critical
                 {
                    // Check again inside critical section to avoid race
                    if (local_verts[v_local_idx].distance > dist_u + weight) {
                         printf("[Rank %d, Thr %d] INS %d->%d improves dist for %d\n", rank, thread_id, u_global, v_global, v_global);
                         local_verts[v_local_idx].distance = dist_u + weight;
                         local_verts[v_local_idx].parent = u_global;
                         local_verts[v_local_idx].affected = true;
                         // If v is a ghost node, need MPI Send update
                    }
                 }
            }
            // Check v -> u (symmetric check if needed)
            if (u_local_idx != -1 && dist_v != INT_MAX && dist_u > dist_v + weight) {
                #pragma omp critical
                {
                   if (local_verts[u_local_idx].distance > dist_v + weight) {
                        printf("[Rank %d, Thr %d] INS %d->%d improves dist for %d\n", rank, thread_id, v_global, u_global, u_global);
                        local_verts[u_local_idx].distance = dist_v + weight;
                        local_verts[u_local_idx].parent = v_global;
                        local_verts[u_local_idx].affected = true;
                        // If u is a ghost node, need MPI Send update
                   }
                }
            }
        }
    }
    print_debug(rank, "Finished Step 1 processing.");
     // MPI Barrier might be needed here to ensure all Step 1 effects (esp. MPI messages) are processed before Step 2
    MPI_Barrier(MPI_COMM_WORLD);
}


void update_distances_step2(VertexData* local_verts, int n_local, int rank, int size /*,...*/) {
    print_debug(rank, "Starting Step 2: Iterative updates.");
    bool global_change_occurred = true;

    while (global_change_occurred) {
        bool local_change_occurred = false;

        // --- Part 1: Propagate Deletions (Algo 3, lines 2-8) ---
        // This needs careful distributed implementation. A breadth-first or depth-first traversal
        // starting from affected_del=true nodes. If traversal hits a ghost node boundary,
        // send MPI message to the owner rank to continue traversal. Use flags to avoid cycles.
        // Set distances to INT_MAX and mark affected=true.
        // This part can be complex and might need iterations itself.
        // Simplified: Assume propagation happens somehow and sets affected=true.
        print_debug(rank, "Step 2a: Propagating deletions (simplified)...");


        // --- Part 2: Iterative Relaxation (Algo 3, lines 9-20) ---
        print_debug(rank, "Step 2b: Relaxing edges for affected nodes...");
        local_change_occurred = false; // Reset for this relaxation phase

        // Create a copy of the 'affected' status before the parallel loop,
        // as 'affected' might be set true *within* the loop for the *next* iteration.
        bool *currently_affected = malloc(n_local * sizeof(bool));
        for(int i=0; i<n_local; ++i) currently_affected[i] = local_verts[i].affected;

        #pragma omp parallel for schedule(dynamic) reduction(||:local_change_occurred)
        for (int i = 0; i < n_local; ++i) {
             int thread_id = omp_get_thread_num();

             if (!currently_affected[i]) continue; // Process only nodes marked affected at start of iteration

             local_verts[i].affected = false; // Reset flag for current node

             int v_global = local_verts[i].id;
             int current_dist_v = local_verts[i].distance;

             // Iterate over neighbors 'n' of vertex 'v' (local index i)
             for (int j = 0; j < local_verts[i].num_edges; ++j) {
                 int neighbor_global = local_verts[i].edges[j].target;
                 int weight = local_verts[i].edges[j].weight;

                 // Determine if neighbor is local or external (ghost)
                 int neighbor_local_idx = -1; /* find_local_index(neighbor_global); */
                 int current_dist_n = INT_MAX;
                 int neighbor_owner_rank = -1; /* find_owner_rank(neighbor_global); */

                 if (neighbor_owner_rank == rank) { // Neighbor is local
                     neighbor_local_idx = /* find_local_index(neighbor_global); */; // Should be quick lookup
                     current_dist_n = local_verts[neighbor_local_idx].distance;
                 } else { // Neighbor is external (ghost node)
                     // ** MPI Communication needed here! **
                     // Send request for distance of neighbor_global to neighbor_owner_rank
                     // Receive the distance current_dist_n
                     // This is a major simplification - needs non-blocking comms or careful structuring
                     printf("[Rank %d, Thr %d] Needs distance for external neighbor %d from Rank %d\n", rank, thread_id, neighbor_global, neighbor_owner_rank);
                     // Placeholder: current_dist_n = get_distance_via_mpi(neighbor_global, neighbor_owner_rank);
                     current_dist_n = INT_MAX; // Simulate not having it easily
                 }


                 // --- Relaxation Checks (Algo 3, lines 13-20) ---
                 // Check if path through v improves n
                 if (current_dist_v != INT_MAX && current_dist_n > current_dist_v + weight) {
                     if (neighbor_local_idx != -1) { // Can only update if n is local
                         #pragma omp critical
                         {
                             // Re-check distance inside critical section
                             if (local_verts[neighbor_local_idx].distance > current_dist_v + weight) {
                                 printf("[Rank %d, Thr %d] Relax %d->%d: Dist[%d] %d -> %d\n", rank, thread_id, v_global, neighbor_global, neighbor_global, local_verts[neighbor_local_idx].distance, current_dist_v + weight);
                                 local_verts[neighbor_local_idx].distance = current_dist_v + weight;
                                 local_verts[neighbor_local_idx].parent = v_global;
                                 local_verts[neighbor_local_idx].affected = true; // Mark for next iteration
                                 local_change_occurred = true;
                             }
                         }
                     } else {
                         // If n is external, we can't update it directly.
                         // The owner rank (neighbor_owner_rank) needs to perform this check based on *its* local distance for v (which might be a ghost for it).
                         // This highlights the complexity of distributed relaxation.
                     }
                 }

                 // Check if path through n improves v (symmetric check)
                 if (current_dist_n != INT_MAX && current_dist_v > current_dist_n + weight) {
                      #pragma omp critical
                      {
                          if (local_verts[i].distance > current_dist_n + weight) {
                               printf("[Rank %d, Thr %d] Relax %d->%d: Dist[%d] %d -> %d\n", rank, thread_id, neighbor_global, v_global, v_global, local_verts[i].distance, current_dist_n + weight);
                               local_verts[i].distance = current_dist_n + weight;
                               local_verts[i].parent = neighbor_global;
                               local_verts[i].affected = true; // Mark for next iteration
                               local_change_occurred = true;
                          }
                      }
                 }
             } // End neighbor loop
        } // End parallel for loop over local vertices
        free(currently_affected);

        // --- Global Synchronization ---
        print_debug(rank, "Finished relaxation iteration. Checking for global convergence...");
        // Use Allreduce to see if *any* process made a change
        MPI_Allreduce(&local_change_occurred, &global_change_occurred, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);
        if (rank == 0) printf("Global change status this iteration: %s\n", global_change_occurred ? "YES" : "NO");

    } // End while loop

    print_debug(rank, "Finished Step 2 iterative updates.");
}


void gather_results(VertexData* local_verts, int n_local, int n_global, int rank, int size) {
    // Simplified: Each process prints its local results
    MPI_Barrier(MPI_COMM_WORLD); // Wait for all ranks to finish computation

    for (int r = 0; r < size; ++r) {
        if (rank == r) {
            printf("\n--- Results from Rank %d ---\n", rank);
            printf("GlobalID\tDistance\tParentID\n");
            for (int i = 0; i < n_local; ++i) {
                // Only print non-ghost nodes if structure distinguishes them
                printf("%d\t\t%d\t\t%d\n",
                       local_verts[i].id,
                       (local_verts[i].distance == INT_MAX) ? -1 : local_verts[i].distance, // Print -1 for infinity
                       local_verts[i].parent);
            }
             printf("---------------------------\n");
             fflush(stdout);
        }
         MPI_Barrier(MPI_COMM_WORLD); // Synchronize printing
    }

    // A proper gather would involve:
    // Slaves sending their local_verts data (ID, dist, parent) to Rank 0.
    // Rank 0 receiving data using MPI_Gatherv (since n_local can vary)
    // Rank 0 assembling the final global distance/parent arrays.
}

// ... Implementations for read_graph, partition_and_distribute_graph etc. are complex ...
// These would involve file I/O, METIS API calls, and significant MPI communication logic.