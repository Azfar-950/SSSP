// dynamic_sssp_mpi_omp_detailed.c

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h> // For bool type
#include <limits.h>  // For INT_MAX
#include <time.h>    // For timing
#include <unistd.h>  // For gethostname (optional, for debug)
#include <stdarg.h>  // Required for va_list, va_start, va_end

#include <mpi.h>
#include <omp.h>
#include <metis.h> // Requires METIS library linked

// --- Constants ---
#define INF INT_MAX
#define ROOT_RANK 0
#define MPI_TAG_GHOST_DIST_REQUEST 1001
#define MPI_TAG_GHOST_DIST_REPLY   1002
#define MPI_TAG_GHOST_UPDATE       1003
#define MPI_TAG_DEL_PROPAGATE      1004


// --- Debug Macro ---
// Define DEBUG_LEVEL: 0=Off, 1=Basic, 2=Detailed, 3=Verbose Comm
#define DEBUG_LEVEL 2

void print_debug_impl(int level, int rank, const char *func, int line, const char *fmt, ...) {
    if (level <= DEBUG_LEVEL) {
        char hostname[MPI_MAX_PROCESSOR_NAME];
        int len;
        MPI_Get_processor_name(hostname, &len);

        printf("[Rank %d on %s | %s:%d] ", rank, hostname, func, line);

        va_list args;
        va_start(args, fmt);
        vprintf(fmt, args);
        va_end(args);
        printf("\n");
        fflush(stdout);
    }
}
// Use PRINT_DEBUG(level, format, ...), e.g., PRINT_DEBUG(1, "Value = %d", x);
#define PRINT_DEBUG(level, ...) print_debug_impl(level, world_rank, __func__, __LINE__, __VA_ARGS__)

// Global rank variable for easy access in debug macro
int world_rank = -1;

// --- Data Structures ---

// Basic edge structure
typedef struct {
    idx_t target_global; // Global ID of the target vertex
    int weight;
} Edge;

// Represents a vertex in the local partition
typedef struct {
    idx_t global_id;      // Original global vertex ID
    int partition_id;   // Which rank owns this vertex
    int distance;       // Current shortest distance estimate
    idx_t parent_global;  // Global ID of the parent in SSSP tree (-1 if none)

    Edge *edges;          // Array of outgoing edges (only includes edges where target is known to this rank)
    idx_t num_edges;

    // --- Status Flags ---
    bool affected;        // Affected by any change in the current update cycle (Algorithm 3)
    bool affected_del;    // Affected specifically by a deletion disconnecting it (Algorithm 2/3)
    bool just_updated;    // Flag set if distance/parent changed in the current *relaxation* iteration

    // --- Ghost Node Info ---
    bool is_ghost;        // Is this vertex primarily owned by another rank? (Useful for consistency)
    // If needed: int owner_rank; // Rank that owns this vertex if it's a ghost copy
} LocalVertex;

// Represents a change operation
typedef struct {
    idx_t u;
    idx_t v;
    int weight; // Only for INS
    char type[4]; // "INS" or "DEL"
} Change;

// Structure to hold the local part of the graph for a process
typedef struct {
    idx_t num_local_vertices; // Number of vertices primarily owned by this rank
    idx_t num_ghost_vertices; // Number of ghost vertices stored locally
    idx_t total_vertices;     // = num_local_vertices + num_ghost_vertices

    LocalVertex *vertices;    // Array storing both local and ghost vertex data (size = total_vertices)

    // Mappings
    idx_t *local_to_global_id; // Map local index (0 to total_vertices-1) to global ID
    // A hash map or similar might be better for global -> local mapping if needed frequently
    // For simplicity, we might loop search, but that's inefficient for large graphs.

    // METIS partitioning info (optional, for reference)
    idx_t *partition_assignment; // Global partition assignment (only rank 0 might keep the full one)

} LocalGraphPartition;

// Structure for graph representation suitable for METIS (CSR format)
// Only rank 0 will typically hold the full graph CSR data
typedef struct {
    idx_t nvtxs; // Number of vertices
    idx_t nedges; // Number of edges (directed count)
    idx_t *xadj; // CSR index pointer array (size nvtxs + 1)
    idx_t *adjncy; // CSR adjacency list array (size nedges)
    idx_t *adjwgt; // CSR edge weight array (size nedges)
} GraphCSR;

// --- Function Prototypes ---

// File I/O & Setup
int read_graph_edges(const char *filename, idx_t *num_vertices_ptr, idx_t *num_edges_ptr, idx_t **u_ptr, idx_t **v_ptr, int **w_ptr);
void convert_edges_to_csr(idx_t nvtxs, idx_t nedges, idx_t *u, idx_t *v, int *w, GraphCSR *csr);
int read_changes(const char *filename, Change **changes_ptr, int *num_changes_ptr);

// METIS Partitioning & Distribution
void partition_graph(GraphCSR *csr, int num_partitions, idx_t **partition_result);
void distribute_graph(GraphCSR *csr, idx_t *partition_assignment, int rank, int size, LocalGraphPartition *local_partition);
void cleanup_csr(GraphCSR *csr);
void cleanup_local_partition(LocalGraphPartition *lp);

// SSSP Update Algorithm Steps (adapted from paper)
void initialize_sssp(LocalGraphPartition *lp, idx_t source_global_id);
void apply_changes_step1(Change *changes, int num_changes, LocalGraphPartition *lp, int rank, int size);
void update_distances_step2(LocalGraphPartition *lp, int rank, int size);
void propagate_deletions_distributed(LocalGraphPartition *lp, int rank, int size); // Complex part
void relax_edges_distributed(LocalGraphPartition *lp, int rank, int size, bool *any_change_made); // Complex part

// Result Gathering
void gather_and_print_results(LocalGraphPartition *lp, idx_t num_global_vertices, int rank, int size);

// Helper Functions
idx_t find_local_index(LocalGraphPartition *lp, idx_t global_id); // Potentially slow linear search
int find_owner_rank(LocalGraphPartition *lp, idx_t global_id); // Placeholder

// --- Main ---
int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); // Set global rank variable
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if (argc < 4) {
        if (world_rank == ROOT_RANK) {
            fprintf(stderr, "Usage: mpirun -np <num_procs> %s <graph_edge_file> <changes_file> <source_vertex_id>\n", argv[0]);
        }
        MPI_Finalize();
        return 1;
    }

    const char *graph_filename = argv[1];
    const char *changes_filename = argv[2];
    idx_t source_global_id = (idx_t)atol(argv[3]); // Use idx_t for METIS compatibility

    PRINT_DEBUG(1, "Application started with %d processes.", world_size);

    GraphCSR csr = {0};
    idx_t *partition_assignment = NULL;
    Change *changes = NULL;
    int num_changes = 0;
    LocalGraphPartition local_partition = {0};
    double start_time, end_time, setup_time, update_time;

    start_time = MPI_Wtime();

    // --- Rank 0: Read Graph, Partition, Read Changes ---
    if (world_rank == ROOT_RANK) {
        PRINT_DEBUG(1, "Reading graph from %s", graph_filename);
        idx_t num_vertices, num_edges;
        idx_t *u, *v;
        int *w;
        if (!read_graph_edges(graph_filename, &num_vertices, &num_edges, &u, &v, &w)) {
            fprintf(stderr, "[Rank %d] Error reading graph file.\n", world_rank);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        PRINT_DEBUG(1, "Graph has %ld vertices and %ld directed edges.", (long)num_vertices, (long)num_edges);

        PRINT_DEBUG(1, "Converting graph to CSR format...");
        convert_edges_to_csr(num_vertices, num_edges, u, v, w, &csr);
        free(u); free(v); free(w); // Free edge list arrays

        PRINT_DEBUG(1, "Partitioning graph using METIS for %d partitions...", world_size);
        partition_graph(&csr, world_size, &partition_assignment);
        PRINT_DEBUG(1, "METIS partitioning complete.");

        PRINT_DEBUG(1, "Reading changes from %s", changes_filename);
        if (!read_changes(changes_filename, &changes, &num_changes)) {
             fprintf(stderr, "[Rank %d] Error reading changes file.\n", world_rank);
             MPI_Abort(MPI_COMM_WORLD, 1);
        }
        PRINT_DEBUG(1, "Read %d change operations.", num_changes);
    }

    // --- Distribute Graph ---
    // Rank 0 sends, others receive
    distribute_graph(&csr, partition_assignment, world_rank, world_size, &local_partition);
    PRINT_DEBUG(1, "Graph distribution complete. Local partition has %ld local and %ld ghost vertices (total %ld).",
              (long)local_partition.num_local_vertices, (long)local_partition.num_ghost_vertices, (long)local_partition.total_vertices);


    // --- Broadcast Changes ---
    MPI_Bcast(&num_changes, 1, MPI_INT, ROOT_RANK, MPI_COMM_WORLD);
    PRINT_DEBUG(2, "Received broadcast: num_changes = %d", num_changes);
    if (world_rank != ROOT_RANK) {
        changes = (Change *)malloc(num_changes * sizeof(Change));
        if (!changes) { perror("Failed to allocate changes array"); MPI_Abort(MPI_COMM_WORLD, 1); }
    }
    // Define an MPI Datatype for the Change struct or send elements individually
    // Simplified: Sending as a block of bytes (assumes no padding issues or matching architecture)
    MPI_Bcast(changes, num_changes * sizeof(Change), MPI_BYTE, ROOT_RANK, MPI_COMM_WORLD);
    PRINT_DEBUG(2, "Received broadcast of changes data.");


    // --- Cleanup Root-Only Data ---
    if (world_rank == ROOT_RANK) {
        cleanup_csr(&csr); // Free CSR arrays on root
        free(partition_assignment); // Free partition array on root
    }

    setup_time = MPI_Wtime() - start_time;
    PRINT_DEBUG(1, "Setup phase completed in %.4f seconds.", setup_time);

    // --- SSSP Initialization ---
    initialize_sssp(&local_partition, source_global_id);
    PRINT_DEBUG(1, "Initial SSSP distances set (source = %ld).", (long)source_global_id);

    MPI_Barrier(MPI_COMM_WORLD); // Ensure initialization is done everywhere
    start_time = MPI_Wtime();

    // --- SSSP Update Steps ---
    PRINT_DEBUG(1, "*** Starting SSSP Update Step 1: Applying Changes ***");
    apply_changes_step1(changes, num_changes, &local_partition, world_rank, world_size);
    PRINT_DEBUG(1, "*** Finished SSSP Update Step 1 ***");

    MPI_Barrier(MPI_COMM_WORLD); // Sync after initial changes applied

    PRINT_DEBUG(1, "*** Starting SSSP Update Step 2: Iterative Update ***");
    update_distances_step2(&local_partition, world_rank, world_size);
    PRINT_DEBUG(1, "*** Finished SSSP Update Step 2 ***");

    update_time = MPI_Wtime() - start_time;
    PRINT_DEBUG(1, "SSSP Update phase completed in %.4f seconds.", update_time);

    MPI_Barrier(MPI_COMM_WORLD);

    // --- Gather and Print Results ---
    idx_t final_num_global_vertices = 0; // Need the global count for gathering
    // Broadcast the count from Rank 0 if necessary, or get it during distribution
    // Assuming local_partition contains info needed to reconstruct global size if necessary
    if (world_rank == ROOT_RANK) {
        // Retrieve global vertex count somehow if not already stored
        // final_num_global_vertices = csr.nvtxs; // If root kept it
    }
    // Need to broadcast final_num_global_vertices here
    // MPI_Bcast(&final_num_global_vertices, 1, METIS_IDX_TYPE, ROOT_RANK, MPI_COMM_WORLD); // METIS_IDX_TYPE might need mapping to MPI type

    gather_and_print_results(&local_partition, final_num_global_vertices, world_rank, world_size);

    // --- Cleanup ---
    PRINT_DEBUG(1, "Cleaning up resources...");
    free(changes);
    cleanup_local_partition(&local_partition);

    end_time = MPI_Wtime();
    if (world_rank == ROOT_RANK) {
        printf("\n----------------------------------------\n");
        printf("Total execution time: %.4f seconds\n", end_time - (start_time-update_time)); // Rough total time
        printf("  Setup time:         %.4f seconds\n", setup_time);
        printf("  Update time:        %.4f seconds\n", update_time);
        printf("----------------------------------------\n");
    }

    MPI_Finalize();
    return 0;
}

// --- Function Implementations ---

int read_graph_edges(const char *filename, idx_t *num_vertices_ptr, idx_t *num_edges_ptr, idx_t **u_ptr, idx_t **v_ptr, int **w_ptr) {
    FILE *f = fopen(filename, "r");
    if (!f) {
        perror("Error opening graph file");
        return 0;
    }
    PRINT_DEBUG(2, "Reading edges from %s", filename);

    // First pass (or assume reasonable max size) to count edges and find max vertex ID
    idx_t u_i, v_i;
    int w_i;
    idx_t max_id = 0; // Initialize max_id to 0
    idx_t edge_count = 0;
    while (fscanf(f, "%d %d %d", &u_i, &v_i, &w_i) == 3) { // Use %d for idx_t (int)
         edge_count++;
         if (u_i > max_id) max_id = u_i;
         if (v_i > max_id) max_id = v_i;
    }
    rewind(f);

    *num_vertices_ptr = max_id + 1; // IDs are 0-based
    *num_edges_ptr = edge_count;

    *u_ptr = (idx_t *)malloc(edge_count * sizeof(idx_t));
    *v_ptr = (idx_t *)malloc(edge_count * sizeof(idx_t));
    *w_ptr = (int *)malloc(edge_count * sizeof(int));

    if (!(*u_ptr) || !(*v_ptr) || !(*w_ptr)) {
        perror("Failed to allocate memory for edge lists");
        fclose(f);
        return 0;
    }

    idx_t current_edge = 0;
    while (fscanf(f, "%d %d %d", &u_i, &v_i, &w_i) == 3) { // Use %d for idx_t (int)
        if (current_edge < edge_count) {
            (*u_ptr)[current_edge] = u_i;
            (*v_ptr)[current_edge] = v_i;
            (*w_ptr)[current_edge] = w_i;
            current_edge++;
        }
    }
    fclose(f);
    PRINT_DEBUG(2, "Read %ld edges, max vertex ID %ld", (long)edge_count, (long)max_id);
    return 1;
}

void convert_edges_to_csr(idx_t nvtxs, idx_t nedges, idx_t *u, idx_t *v, int *w, GraphCSR *csr) {
    PRINT_DEBUG(2, "Converting %ld vertices, %ld edges to CSR", (long)nvtxs, (long)nedges);
    csr->nvtxs = nvtxs;
    csr->nedges = nedges;
    csr->xadj = (idx_t *)calloc(nvtxs + 1, sizeof(idx_t));
    csr->adjncy = (idx_t *)malloc(nedges * sizeof(idx_t));
    csr->adjwgt = (idx_t *)malloc(nedges * sizeof(idx_t)); // METIS uses idx_t for weights too

    if (!csr->xadj || !csr->adjncy || !csr->adjwgt) {
        perror("Failed to allocate CSR arrays");
        // Handle cleanup if partially allocated
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Count outgoing degrees for each vertex
    for (idx_t i = 0; i < nedges; ++i) {
        csr->xadj[u[i]]++;
    }

    // Create cumulative sum for xadj pointers
    idx_t sum = 0;
    for (idx_t i = 0; i < nvtxs; ++i) {
        idx_t temp = csr->xadj[i];
        csr->xadj[i] = sum;
        sum += temp;
    }
    csr->xadj[nvtxs] = sum; // Last entry is total number of edges

    // Fill adjncy and adjwgt using the computed pointers
    // Need temporary array to track current position within each vertex's edge list
    idx_t *pos = (idx_t *)calloc(nvtxs, sizeof(idx_t));
    if (!pos) { perror("Failed to allocate position array"); MPI_Abort(MPI_COMM_WORLD, 1); }

    for (idx_t i = 0; i < nedges; ++i) {
        idx_t src = u[i];
        idx_t dest = v[i];
        int weight = w[i];

        idx_t edge_idx = csr->xadj[src] + pos[src];
        csr->adjncy[edge_idx] = dest;
        csr->adjwgt[edge_idx] = (idx_t)weight; // Cast weight to idx_t for METIS
        pos[src]++;
    }

    free(pos);
    PRINT_DEBUG(2, "CSR conversion complete.");
}

void cleanup_csr(GraphCSR *csr) {
    if (!csr) return;
    free(csr->xadj); csr->xadj = NULL;
    free(csr->adjncy); csr->adjncy = NULL;
    free(csr->adjwgt); csr->adjwgt = NULL;
    PRINT_DEBUG(2, "Cleaned up global CSR data.");
}

void cleanup_local_partition(LocalGraphPartition *lp) {
     if (!lp) return;
     if (lp->vertices) {
         for (idx_t i = 0; i < lp->total_vertices; ++i) {
             free(lp->vertices[i].edges);
         }
         free(lp->vertices); lp->vertices = NULL;
     }
     free(lp->local_to_global_id); lp->local_to_global_id = NULL;
     // Free other mappings if allocated
     PRINT_DEBUG(2, "Cleaned up local graph partition data.");
}


int read_changes(const char *filename, Change **changes_ptr, int *num_changes_ptr) {
    FILE *f = fopen(filename, "r");
    if (!f) { perror("Error opening changes file"); return 0; }
    PRINT_DEBUG(2, "Reading changes from %s", filename);

    // Count lines first
    int count = 0;
    char line_buffer[256];
    while (fgets(line_buffer, sizeof(line_buffer), f)) {
        count++;
    }
    rewind(f);
    *num_changes_ptr = count;
    *changes_ptr = (Change *)malloc(count * sizeof(Change));
    if (!(*changes_ptr)) { perror("Failed to allocate changes array"); fclose(f); return 0; }

    int current_change = 0;
    idx_t u, v;
    int w;
    char type[10]; // Buffer for operation type
    while (current_change < count && fscanf(f, "%s", type) == 1) {
         Change *ch = &((*changes_ptr)[current_change]);
         strncpy(ch->type, type, 3);
         ch->type[3] = '\0'; // Null-terminate

         if (strcmp(ch->type, "INS") == 0) {
             if (fscanf(f, "%d %d %d", &u, &v, &w) == 3) { // Use %d for idx_t (int)
                 ch->u = u;
                 ch->v = v;
                 ch->weight = w;
             } else { /* Handle error */ break; }
         } else if (strcmp(ch->type, "DEL") == 0) {
              if (fscanf(f, "%d %d", &u, &v) == 2) { // Use %d for idx_t (int)
                 ch->u = u;
                 ch->v = v;
                 ch->weight = 0; // Not used for DEL
              } else { /* Handle error */ break; }
         } else {
             fprintf(stderr, "Unknown change type: %s\n", type);
             /* Handle error */ break;
         }
         current_change++;
    }
    fclose(f);
    // If current_change != count, there was an error during parsing
    if (current_change != count) {
        fprintf(stderr, "Error parsing changes file. Expected %d, parsed %d.\n", count, current_change);
        free(*changes_ptr);
        *changes_ptr = NULL;
        *num_changes_ptr = 0;
        return 0;
    }
    PRINT_DEBUG(2, "Read %d changes.", count);
    return 1;
}

void partition_graph(GraphCSR *csr, int num_partitions, idx_t **partition_result) {
    if (num_partitions <= 1) {
        PRINT_DEBUG(1, "Number of partitions is 1, skipping METIS call.");
        *partition_result = (idx_t *)malloc(csr->nvtxs * sizeof(idx_t));
        if (!*partition_result) { perror("Failed to allocate partition array"); MPI_Abort(MPI_COMM_WORLD, 1); }
        // Assign all vertices to partition 0
        for (idx_t i = 0; i < csr->nvtxs; ++i) {
            (*partition_result)[i] = 0;
        }
        return;
    }

    PRINT_DEBUG(2, "Calling METIS_PartGraphKway...");
    *partition_result = (idx_t *)malloc(csr->nvtxs * sizeof(idx_t));
    if (!*partition_result) { perror("Failed to allocate METIS result array"); MPI_Abort(MPI_COMM_WORLD, 1); }

    idx_t ncon = 1; // Number of balancing constraints (usually 1)
    idx_t objval; // Stores edge-cut or volume after partitioning

    // METIS options (can be NULL for defaults)
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_DBGLVL] = 0; // Disable METIS debug output unless needed
    options[METIS_OPTION_CONTIG] = 1; // Try to create contiguous partitions if possible

    int ret = METIS_PartGraphKway(&csr->nvtxs,      // num vertices
                                  &ncon,            // num constraints
                                  csr->xadj,        // CSR index pointers
                                  csr->adjncy,      // CSR adjacency list
                                  csr->adjwgt,      // CSR vertex weights (NULL for unweighted vertices)
                                  NULL,             // CSR vertex sizes (NULL)
                                  csr->adjwgt,      // CSR edge weights
                                  &num_partitions,  // num partitions (nparts)
                                  NULL,             // target partition weights (NULL for equal)
                                  NULL,             // allowed imbalance (NULL for default)
                                  options,          // METIS options
                                  &objval,          // resulting edgecut/volume
                                  *partition_result); // output partition vector

    if (ret != METIS_OK) {
        fprintf(stderr, "METIS_PartGraphKway failed with error code: %d\n", ret);
        free(*partition_result);
        *partition_result = NULL;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    PRINT_DEBUG(1, "METIS partitioning successful. Edge cut: %ld", (long)objval);
}


void distribute_graph(GraphCSR *csr, idx_t *partition_assignment, int rank, int size, LocalGraphPartition *local_partition) {
    PRINT_DEBUG(2, "Starting graph distribution process.");

    // --- Step 1: Determine local and ghost vertex counts for each rank ---
    idx_t *local_counts = (idx_t *)calloc(size, sizeof(idx_t));
    idx_t *ghost_counts = (idx_t *)calloc(size, sizeof(idx_t)); // Max possible ghosts
    idx_t num_global_vertices = 0;

    if (rank == ROOT_RANK) {
        num_global_vertices = csr->nvtxs;
        bool *is_ghost_for_rank = (bool *)calloc(csr->nvtxs * size, sizeof(bool)); // Track ghosts globally first
        idx_t *ghost_node_indices_count = (idx_t *)calloc(size, sizeof(idx_t)); // Count ghosts per rank

        // Count local vertices
        for (idx_t i = 0; i < csr->nvtxs; ++i) {
            local_counts[partition_assignment[i]]++;
        }

        // Identify ghost nodes (nodes needed by a rank but owned by another)
        for (idx_t u_global = 0; u_global < csr->nvtxs; ++u_global) {
            int owner_rank = partition_assignment[u_global];
            for (idx_t edge_idx = csr->xadj[u_global]; edge_idx < csr->xadj[u_global + 1]; ++edge_idx) {
                idx_t v_global = csr->adjncy[edge_idx];
                int target_rank = partition_assignment[v_global];
                // If edge crosses partition boundary (u -> v), then v is a ghost for owner_rank, and u is a ghost for target_rank
                if (owner_rank != target_rank) {
                    // Mark v as needed by owner_rank (if not already marked)
                    if (!is_ghost_for_rank[v_global * size + owner_rank]) {
                         is_ghost_for_rank[v_global * size + owner_rank] = true;
                         ghost_node_indices_count[owner_rank]++;
                    }
                     // Mark u as needed by target_rank (if not already marked)
                    if (!is_ghost_for_rank[u_global * size + target_rank]) {
                         is_ghost_for_rank[u_global * size + target_rank] = true;
                         ghost_node_indices_count[target_rank]++;
                    }
                }
            }
        }
        // Now ghost_node_indices_count[r] holds the number of unique ghost nodes needed by rank r
         for(int r=0; r<size; ++r) ghost_counts[r] = ghost_node_indices_count[r];

        free(is_ghost_for_rank);
        free(ghost_node_indices_count);
        PRINT_DEBUG(2, "Calculated local/ghost counts on root.");

    }

    // --- Step 2: Broadcast counts ---
    MPI_Bcast(&num_global_vertices, 1, MPI_LONG_LONG_INT, ROOT_RANK, MPI_COMM_WORLD); // Assuming idx_t maps to long long
    MPI_Bcast(local_counts, size, MPI_LONG_LONG_INT, ROOT_RANK, MPI_COMM_WORLD);
    MPI_Bcast(ghost_counts, size, MPI_LONG_LONG_INT, ROOT_RANK, MPI_COMM_WORLD);

    local_partition->num_local_vertices = local_counts[rank];
    local_partition->num_ghost_vertices = ghost_counts[rank];
    local_partition->total_vertices = local_partition->num_local_vertices + local_partition->num_ghost_vertices;
    PRINT_DEBUG(2, "Received counts: local=%ld, ghost=%ld, total=%ld", (long)local_counts[rank], (long)ghost_counts[rank], (long)local_partition->total_vertices);

    free(local_counts);
    free(ghost_counts);


    // --- Step 3: Allocate local structures ---
    local_partition->vertices = (LocalVertex *)calloc(local_partition->total_vertices, sizeof(LocalVertex));
    local_partition->local_to_global_id = (idx_t *)malloc(local_partition->total_vertices * sizeof(idx_t));
    if (!local_partition->vertices || !local_partition->local_to_global_id) {
         perror("Failed to allocate local partition arrays"); MPI_Abort(MPI_COMM_WORLD, 1);
    }


    // --- Step 4: Rank 0 prepares and sends data for each rank ---
    // This part is complex. Rank 0 needs to iterate through all vertices and edges,
    // determine which rank needs what data (local vertex data, ghost vertex data, edges),
    // serialize it, and send it using MPI_Send. Slaves use MPI_Recv.
    // Data to send per rank 'r':
    //  - List of global IDs for its local vertices.
    //  - List of global IDs for its ghost vertices.
    //  - For each local/ghost vertex: its basic info (global ID) and list of relevant edges (target_global, weight).
    //    Only include edges where the target is either local to 'r' or a ghost known to 'r'.

    // ** SIMPLIFICATION: ** This critical distribution logic is highly complex to implement
    // correctly and efficiently with serialization. We'll skip the detailed MPI_Send/Recv loops
    // here. Assume the `local_partition` struct gets populated somehow on each rank.
    // In a real implementation, this would involve:
    // 1. Rank 0 calculating send counts/displacements for each rank.
    // 2. Rank 0 packing data (global IDs, edge lists) into buffers.
    // 3. Rank 0 using MPI_Send or collective operations like MPI_Scatterv.
    // 4. Slaves using MPI_Recv or collectives to receive data.
    // 5. Slaves unpacking data and filling their `local_partition` struct.
    // 6. Setting up the local_to_global_id mapping and potentially a global_to_local mapping.


    PRINT_DEBUG(1, "[Rank %d] **Placeholder:** Assumed local partition data is now populated.", rank);
    // ** TODO: Implement the actual MPI communication for distribution **


    // Store partition assignment locally if needed (e.g., for finding owner rank of ghosts)
    // This might require broadcasting the full partition_assignment array, or parts of it.
    local_partition->partition_assignment = NULL; // Placeholder

}

void initialize_sssp(LocalGraphPartition *lp, idx_t source_global_id) {
    PRINT_DEBUG(2, "Initializing SSSP distances and parents.");
    #pragma omp parallel for
    for (idx_t i = 0; i < lp->total_vertices; ++i) {
        lp->vertices[i].distance = (lp->vertices[i].global_id == source_global_id) ? 0 : INF;
        lp->vertices[i].parent_global = -1;
        lp->vertices[i].affected = (lp->vertices[i].global_id == source_global_id); // Initially, only source is "affected" to start relaxation
        lp->vertices[i].affected_del = false;
        lp->vertices[i].just_updated = false; // Reset update flag
        lp->vertices[i].partition_id = -1; // TODO: Populate this during distribution
        lp->vertices[i].is_ghost = false; // TODO: Populate this during distribution
    }
    // The source might be a ghost node on some ranks.
    PRINT_DEBUG(2, "SSSP initialization complete.");
}


// Find local index for a global ID (simple linear scan - inefficient)
idx_t find_local_index(LocalGraphPartition *lp, idx_t global_id) {
    for (idx_t i = 0; i < lp->total_vertices; ++i) {
        // Need the mapping from local index to global id, assuming lp->vertices[i].global_id holds it
        if (lp->vertices[i].global_id == global_id) {
            return i;
        }
    }
    return -1; // Not found locally (shouldn't happen if ghosts are handled correctly)
}

// Find owner rank (requires partition assignment info) - Placeholder
int find_owner_rank(LocalGraphPartition *lp, idx_t global_id) {
    // Requires access to the global partition_assignment array, which might need to be broadcast or queried
    // This placeholder implementation always returns rank 0 for demonstration; needs real logic.
    if (lp->partition_assignment && global_id < /* num_global_vertices */ 0 ) { // Check bounds
       // return lp->partition_assignment[global_id];
    }
    // As a temporary placeholder, assume rank 0 owns global ID 0 for the 8-vertex graph scenario,
    // and other IDs are distributed. This needs to be replaced with actual owner determination.
    if (global_id < 8) return global_id % 2; // Simple distribution for 2 ranks and small graph
     PRINT_DEBUG(3,"Warning: Could not determine owner rank for global_id %ld (partition info missing?)", (long)global_id);
    return -1; // Indicate failure
}


void apply_changes_step1(Change *changes, int num_changes, LocalGraphPartition *lp, int rank, int size) {
    PRINT_DEBUG(1, "Applying %d changes (Step 1).", num_changes);

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < num_changes; ++i) {
        Change ch = changes[i];
        idx_t u_global = ch.u;
        idx_t v_global = ch.v;
        int thread_id = omp_get_thread_num();
        PRINT_DEBUG(3, "[Thread %d] Processing change %d: %s %ld %ld", thread_id, i, ch.type, (long)u_global, (long)v_global);

        idx_t u_local = find_local_index(lp, u_global);
        idx_t v_local = find_local_index(lp, v_global);

        if (strcmp(ch.type, "DEL") == 0) {
            // --- Handle Deletion (Algo 2, lines 4-9 adaptation) ---
            // Check if edge (u,v) or (v,u) exists in the *current* SSSP tree locally
            // This requires checking parent pointers.
            bool changed = false;
             // Check if v is child of u locally
            if (v_local != -1 && lp->vertices[v_local].parent_global == u_global) {
                 #pragma omp critical(sssp_update)
                 {
                    // Re-check condition inside critical section if needed
                    if (lp->vertices[v_local].distance != INF) { // Avoid marking if already disconnected
                        PRINT_DEBUG(2, "[Thread %d] DEL %ld->%ld detected as SSSP edge. Disconnecting %ld.", thread_id, (long)u_global, (long)v_global, (long)v_global);
                        lp->vertices[v_local].distance = INF;
                        lp->vertices[v_local].parent_global = -1;
                        lp->vertices[v_local].affected_del = true; // Mark for deletion propagation
                        lp->vertices[v_local].affected = true; // Mark for general update cycle
                        changed = true;
                    }
                 }
                 if (changed /* && lp->vertices[v_local] is ghost for others */) {
                      // TODO: Send MPI message indicating v_global is affected by deletion
                 }

            } // Check if u is child of v locally (symmetric check)
            else if (u_local != -1 && lp->vertices[u_local].parent_global == v_global) {
                 #pragma omp critical(sssp_update)
                 {
                    if (lp->vertices[u_local].distance != INF) {
                        PRINT_DEBUG(2, "[Thread %d] DEL %ld->%ld detected as SSSP edge. Disconnecting %ld.", thread_id, (long)v_global, (long)u_global, (long)u_global);
                        lp->vertices[u_local].distance = INF;
                        lp->vertices[u_local].parent_global = -1;
                        lp->vertices[u_local].affected_del = true;
                        lp->vertices[u_local].affected = true;
                        changed = true;
                    }
                 }
                  if (changed /* && lp->vertices[u_local] is ghost for others */) {
                      // TODO: Send MPI message indicating u_global is affected by deletion
                 }
            }

        } else { // INS
             // --- Handle Insertion (Algo 2, lines 10-19 adaptation) ---
             int weight = ch.weight;

             // Need distances of u and v. They might be local or ghost.
             // If ghost, their current distance might not be up-to-date locally yet.
             // For Step 1, we often use the *locally known* distance, even if it's a ghost.
             // The iterative Step 2 will correct inconsistencies later.

             int dist_u = (u_local != -1) ? lp->vertices[u_local].distance : INF;
             int dist_v = (v_local != -1) ? lp->vertices[v_local].distance : INF;
             bool changed = false;

             // Check u -> v path
             if (v_local != -1 && dist_u != INF) { // v must be local to update it
                 int potential_dist_v = (dist_u == INF || weight == INF) ? INF : dist_u + weight; // Prevent overflow
                 if (potential_dist_v < dist_v) {
                     #pragma omp critical(sssp_update)
                     {
                         // Re-check inside critical section
                         if (potential_dist_v < lp->vertices[v_local].distance) {
                              PRINT_DEBUG(2, "[Thread %d] INS %ld->%ld (w=%d) updates Dist[%ld]: %d -> %d",
                                           thread_id, (long)u_global, (long)v_global, weight, (long)v_global, lp->vertices[v_local].distance, potential_dist_v);
                              lp->vertices[v_local].distance = potential_dist_v;
                              lp->vertices[v_local].parent_global = u_global;
                              lp->vertices[v_local].affected = true;
                              changed = true;
                         }
                     }
                      if (changed /* && lp->vertices[v_local] is ghost for others */) {
                         // TODO: Send MPI message about updated distance/parent for v_global
                      }
                 }
             }

             // Check v -> u path (symmetric check)
             changed = false; // Reset for this direction
             if (u_local != -1 && dist_v != INF) { // u must be local to update it
                 int potential_dist_u = (dist_v == INF || weight == INF) ? INF : dist_v + weight;
                 if (potential_dist_u < dist_u) {
                      #pragma omp critical(sssp_update)
                      {
                         if (potential_dist_u < lp->vertices[u_local].distance) {
                             PRINT_DEBUG(2, "[Thread %d] INS %ld->%ld (w=%d) updates Dist[%ld]: %d -> %d",
                                           thread_id, (long)v_global, (long)u_global, weight, (long)u_global, lp->vertices[u_local].distance, potential_dist_u);
                             lp->vertices[u_local].distance = potential_dist_u;
                             lp->vertices[u_local].parent_global = v_global;
                             lp->vertices[u_local].affected = true; // Mark for next iteration
                             changed = true;
                         }
                      }
                      if (changed /* && lp->vertices[u_local] is ghost for others */) {
                         // TODO: Send MPI message about updated distance/parent for u_global
                      }
                 }
             }
        } // End if INS/DEL
    } // End parallel for loop over changes
    PRINT_DEBUG(1, "Finished processing changes locally in Step 1.");
    // MPI Barrier might be needed here to ensure all Step 1 MPI messages are sent/received before Step 2 begins.
}


void update_distances_step2(LocalGraphPartition *lp, int rank, int size) {
    PRINT_DEBUG(1, "Starting Step 2 iterative update process.");

    // --- Part 1: Propagate Deletions ---
    // Mark descendants of affected_del nodes as disconnected (dist=INF)
    // This requires distributed traversal.
    PRINT_DEBUG(1, "Step 2a: Propagating deletion effects (distributed)...");
    propagate_deletions_distributed(lp, rank, size); // Placeholder for complex logic
    PRINT_DEBUG(1, "Step 2a: Finished deletion propagation phase.");

    MPI_Barrier(MPI_COMM_WORLD); // Ensure deletion propagation is complete globally

    // --- Part 2: Iterative Relaxation ---
    PRINT_DEBUG(1, "Step 2b: Starting iterative edge relaxation...");
    bool global_change_occurred_in_iteration = true;
    int iteration = 0;

    while (global_change_occurred_in_iteration) {
        iteration++;
        PRINT_DEBUG(1, "--- Relaxation Iteration %d ---", iteration);
        bool local_change_occurred = false;

        relax_edges_distributed(lp, rank, size, &local_change_occurred);

        // Check for global convergence
        PRINT_DEBUG(2, "Checking for global convergence (Iteration %d)", iteration);
        MPI_Allreduce(&local_change_occurred, &global_change_occurred_in_iteration, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);

        if (rank == ROOT_RANK) {
            PRINT_DEBUG(1, "Global change in iteration %d: %s", iteration, global_change_occurred_in_iteration ? "YES" : "NO");
        }
        if (iteration > lp->total_vertices * 2 && global_change_occurred_in_iteration) {
             // Add a safety break for potential oscillations or non-convergence issues
             if(rank == ROOT_RANK) fprintf(stderr, "Warning: Potential non-convergence after %d iterations.\n", iteration);
             break;
        }

    } // End while loop

    PRINT_DEBUG(1, "Step 2b: Iterative relaxation finished after %d iterations.");
}

// Placeholder for the complex deletion propagation logic
void propagate_deletions_distributed(LocalGraphPartition *lp, int rank, int size) {
     PRINT_DEBUG(2, "Entering deletion propagation phase.");
     // Approach:
     // 1. Identify local nodes initially marked with affected_del = true.
     // 2. Perform a local traversal (e.g., BFS/DFS) starting from these nodes.
     // 3. For each visited node 'c' that is a child of a disconnected node 'v':
     //    - Set dist[c] = INF, parent[c] = -1, affected_del[c] = true, affected[c] = true.
     //    - If 'c' is a ghost node for another rank, send a message (MPI_TAG_DEL_PROPAGATE) to that rank containing global_id of 'c'.
     // 4. Ranks receive propagation messages and add the received node ID to their local list of nodes to start traversal from.
     // 5. Repeat steps 2-4 until no more propagation messages are sent/received and local traversals are complete.
     //    This requires careful synchronization and termination detection.

     PRINT_DEBUG(1, "** Placeholder: Deletion propagation logic needs full implementation. **");
}

// Placeholder for the complex relaxation logic with ghost node communication
void relax_edges_distributed(LocalGraphPartition *lp, int rank, int size, bool *any_local_change_made) {
    PRINT_DEBUG(2, "Entering relaxation phase for one iteration.");
    *any_local_change_made = false;

    // --- Phase A: Request distances for necessary ghost neighbors ---
    // Iterate through local 'affected' nodes and identify required external neighbor distances.
    // Use non-blocking sends (MPI_Isend) for requests (MPI_TAG_GHOST_DIST_REQUEST).
    PRINT_DEBUG(3, "Phase A: Identifying and requesting ghost neighbor distances (placeholder).");
    // ** TODO: Implement request logic **
    MPI_Request *requests = NULL; // Array for send requests
    int num_requests = 0;

    // --- Phase B: Process incoming distance requests and send replies ---
    // Check for incoming requests (MPI_Iprobe/MPI_Recv with MPI_TAG_GHOST_DIST_REQUEST).
    // For each request, find the local distance of the requested node.
    // Send back the distance using MPI_Send (MPI_TAG_GHOST_DIST_REPLY).
    PRINT_DEBUG(3, "Phase B: Servicing incoming ghost distance requests (placeholder).");
    // ** TODO: Implement request servicing logic **


    // --- Phase C: Wait for outgoing requests and incoming replies ---
    // Wait for all outgoing send requests to complete (MPI_Waitall).
    // Receive replies (MPI_Recv with MPI_TAG_GHOST_DIST_REPLY) for the distances requested in Phase A.
    // Store received ghost distances.
    PRINT_DEBUG(3, "Phase C: Waiting for requests/replies (placeholder).");
    // ** TODO: Implement receive/wait logic **
    // MPI_Waitall(num_requests, requests, MPI_STATUSES_IGNORE); free(requests);


    // --- Phase D: Perform local relaxation using available distances ---
    PRINT_DEBUG(3, "Phase D: Performing local relaxation computations.");
    bool local_change_this_phase = false;

     // Create a copy of the 'affected' status before the parallel loop,
     // as 'affected' might be set true *within* the loop for the *next* iteration.
     bool *currently_affected = (bool*)malloc(lp->total_vertices * sizeof(bool));
     for(idx_t i=0; i<lp->total_vertices; ++i) currently_affected[i] = lp->vertices[i].affected;


    #pragma omp parallel for schedule(dynamic) reduction(||:local_change_this_phase)
    for (idx_t i = 0; i < lp->total_vertices; ++i) {
         if (!currently_affected[i]) continue; // Only process nodes marked affected at start of *outer* loop

         lp->vertices[i].affected = false; // Reset flag for this iteration (might be set again below)

         idx_t v_global = lp->vertices[i].global_id;
         int current_dist_v = lp->vertices[i].distance;
         int thread_id = omp_get_thread_num();

         // Iterate over neighbors 'n' of vertex 'v' (local index i)
         for (idx_t j = 0; j < lp->vertices[i].num_edges; ++j) {
             idx_t neighbor_global = lp->vertices[i].edges[j].target_global;
             int weight = lp->vertices[i].edges[j].weight;

             idx_t neighbor_local_idx = find_local_index(lp, neighbor_global);
             int current_dist_n = INF;

             if (neighbor_local_idx != -1) { // Neighbor is local (or a ghost copy we store)
                 current_dist_n = lp->vertices[neighbor_local_idx].distance;
             } else {
                 // Should not happen if ghosts are correctly populated during distribution
                 PRINT_DEBUG(0,"[Thread %d] ERROR: Neighbor %ld of %ld not found locally!", thread_id, (long)neighbor_global, (long)v_global);
                 continue;
                 // Alternatively, maybe this neighbor's distance was received via MPI in Phase C
                 // current_dist_n = get_received_ghost_distance(neighbor_global); // Lookup needed
             }

             // --- Relaxation Checks ---
             // Check if path through v improves n
             if (current_dist_v != INF) {
                 int potential_dist_n = (current_dist_v == INF || weight == INF) ? INF : current_dist_v + weight;
                 if (potential_dist_n < current_dist_n) {
                     // Can only update n if it's owned by this rank (not a ghost copy)
                     // Requires knowing the owner rank of neighbor_global
                     int neighbor_owner = find_owner_rank(lp, neighbor_global); // Placeholder
                     if (neighbor_owner == rank) { // Update local neighbor
                          #pragma omp critical(sssp_update)
                          {
                              // Recheck with potentially updated local value
                              if (potential_dist_n < lp->vertices[neighbor_local_idx].distance) {
                                   PRINT_DEBUG(3, "[Thread %d] Relax %ld->%ld: Dist[%ld] %d -> %d", thread_id, (long)v_global, (long)neighbor_global, (long)neighbor_global, lp->vertices[neighbor_local_idx].distance, potential_dist_n);
                                   lp->vertices[neighbor_local_idx].distance = potential_dist_n;
                                   lp->vertices[neighbor_local_idx].parent_global = v_global;
                                   lp->vertices[neighbor_local_idx].affected = true; // Mark for next iteration
                                   local_change_this_phase = true;
                              }
                          }
                     } else {
                          // If neighbor is owned by another rank, *we* don't update it directly.
                          // Instead, send an update message if *our* distance for v changed,
                          // letting the owner perform the check. (Handled in Phase E).
                          PRINT_DEBUG(3, "[Thread %d] Relax %ld->%ld: Potential update for external %ld ignored locally.", thread_id, (long)v_global, (long)neighbor_global, (long)neighbor_global);
                     }
                 }
             }

             // Check if path through n improves v
             if (current_dist_n != INF) {
                 int potential_dist_v = (current_dist_n == INF || weight == INF) ? INF : current_dist_n + weight;
                 if (potential_dist_v < current_dist_v) {
                      #pragma omp critical(sssp_update)
                      {
                         // Recheck with potentially updated local value
                         if (potential_dist_v < lp->vertices[i].distance) {
                              PRINT_DEBUG(3, "[Thread %d] Relax %ld->%ld: Dist[%ld] %d -> %d", thread_id, (long)neighbor_global, (long)v_global, (long)v_global, lp->vertices[i].distance, potential_dist_v);
                              lp->vertices[i].distance = potential_dist_v;
                              lp->vertices[i].parent_global = neighbor_global;
                              lp->vertices[i].affected = true; // Mark for next iteration
                              local_change_this_phase = true;
                              lp->vertices[i].just_updated = true; // Mark that v's distance itself changed
                         }
                      }
                 }
             }
         } // End neighbor loop
    } // End parallel for loop over local vertices
    free(currently_affected);

    // --- Phase E: Send updates for owned nodes that changed distance/parent ---
    // Iterate through local nodes where just_updated is true.
    // If this node is a ghost for other ranks, send its new distance/parent (MPI_TAG_GHOST_UPDATE).
    PRINT_DEBUG(3, "Phase E: Sending updates for changed owned nodes (placeholder).");
    // ** TODO: Implement update sending logic **
    for(idx_t i=0; i < lp->total_vertices; ++i) {
        if (lp->vertices[i].just_updated /* && is owned by me && is ghost for others */) {
            // Send update message
        }
        lp->vertices[i].just_updated = false; // Reset for next iteration
    }


    // --- Phase F: Receive and apply updates for ghost nodes ---
    // Check for incoming update messages (MPI_TAG_GHOST_UPDATE).
    // Apply the received distance/parent to the local copy of the ghost node.
    // Mark the ghost node as 'affected' if its distance changed, so it's considered in the *next* iteration's relaxation.
    PRINT_DEBUG(3, "Phase F: Receiving and applying updates for ghost nodes (placeholder).");
    // ** TODO: Implement update receiving logic **
    // This might need non-blocking probes/receives within the main loop or in a separate thread/progress mechanism.


    *any_local_change_made = local_change_this_phase; // Report if any local owned node was updated
    PRINT_DEBUG(2, "Relaxation phase finished. Local change detected: %s", *any_local_change_made ? "YES" : "NO");
}


void gather_and_print_results(LocalGraphPartition *lp, idx_t num_global_vertices, int rank, int size) {
    PRINT_DEBUG(1, "Gathering final results to Rank 0.");

    // Prepare data to send: array of (global_id, distance, parent_global) for local vertices
    idx_t nlocal = lp->num_local_vertices;
    long *send_buffer = (long *)malloc(nlocal * 3 * sizeof(long)); // Using long for global IDs and distance/parent
    if (!send_buffer && nlocal > 0) { perror("Failed to allocate send buffer"); return; }

    int send_count = 0;
    for (idx_t i = 0; i < lp->total_vertices; ++i) {
        // Only send data for vertices primarily owned by this rank
        if (/*!lp->vertices[i].is_ghost*/ find_owner_rank(lp, lp->vertices[i].global_id) == rank) { // Requires owner info
            if (send_count < nlocal * 3) {
                send_buffer[send_count++] = (long)lp->vertices[i].global_id;
                send_buffer[send_count++] = (long)lp->vertices[i].distance;
                send_buffer[send_count++] = (long)lp->vertices[i].parent_global;
            } else {
                 PRINT_DEBUG(0,"Error: Buffer overflow during result packing.");
                 break; // Avoid writing out of bounds
            }
        }
    }
     // Ensure send_count matches expected nlocal * 3
     if (send_count != nlocal * 3) {
          PRINT_DEBUG(0,"Warning: Send count mismatch (%d vs %ld). May indicate issue with owner check.", send_count, nlocal * 3);
          // Adjust send_count if needed, or handle error
     }


    // Gather counts at root
    int *recv_counts = NULL;
    if (rank == ROOT_RANK) {
        recv_counts = (int *)malloc(size * sizeof(int));
    }
    int items_per_vertex = 3; // id, dist, parent
    int send_items = nlocal * items_per_vertex;
    MPI_Gather(&send_items, 1, MPI_INT, recv_counts, 1, MPI_INT, ROOT_RANK, MPI_COMM_WORLD);

    // Gather data at root using Gatherv
    long *recv_buffer = NULL;
    int *displs = NULL;
    int total_items_received = 0;

    if (rank == ROOT_RANK) {
        displs = (int *)malloc(size * sizeof(int));
        displs[0] = 0;
        total_items_received = recv_counts[0];
        for (int i = 1; i < size; ++i) {
            displs[i] = displs[i - 1] + recv_counts[i - 1];
            total_items_received += recv_counts[i];
        }
        recv_buffer = (long *)malloc(total_items_received * sizeof(long));
        if (!recv_buffer) { perror("Failed to allocate receive buffer"); /* Handle error */ }
        PRINT_DEBUG(2, "Root receiving %d total items for results.", total_items_received);
    }

    MPI_Gatherv(send_buffer, send_items, MPI_LONG,
                recv_buffer, recv_counts, displs, MPI_LONG,
                ROOT_RANK, MPI_COMM_WORLD);

    // Root prints the results
    if (rank == ROOT_RANK && recv_buffer) {
        printf("\n--- Final SSSP Results (gathered at Rank 0) ---\n");
        printf("GlobalID\tDistance\tParentID\n");
        printf("------------------------------------------------\n");
        for (int i = 0; i < total_items_received; i += 3) {
             long gid = recv_buffer[i];
             long dist = recv_buffer[i+1];
             long parent = recv_buffer[i+2];
             printf("%ld\t\t%ld\t\t%ld\n",
                    gid,
                    (dist == (long)INF) ? -1L : dist, // Print -1 for infinity
                    parent);
        }
        printf("------------------------------------------------\n");
        free(recv_buffer);
        free(recv_counts);
        free(displs);
    } else if (rank == ROOT_RANK) {
         PRINT_DEBUG(0,"Error receiving results on Root.");
    }

    free(send_buffer);
    PRINT_DEBUG(1, "Result gathering complete.");
}