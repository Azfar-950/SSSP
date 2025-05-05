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


// CSR representation on all ranks
int *g_xadj = NULL, *g_adjncy = NULL, *g_adjwgt = NULL;
int g_nvtxs = 0;
int *partition_assignment = NULL;
LocalVertex *local_verts = NULL;
int local_count, ghost_count, total_count;
int *local_map = NULL; // maps local index -> global id
// --- Function Prototypes ---

// File I/O & Setup
int read_graph_edges(const char *filename, idx_t *num_vertices_ptr, idx_t *num_edges_ptr, idx_t **u_ptr, idx_t **v_ptr, int **w_ptr);
void convert_edges_to_csr(idx_t nvtxs, idx_t nedges, idx_t *u, idx_t *v, int *w, GraphCSR *csr);
int read_changes(const char *filename, Change **changes_ptr, int *num_changes_ptr);

// METIS Partitioning & Distribution
void partition_graph(GraphCSR *csr, int num_partitions, idx_t **partition_result);
void distribute_graph();
void cleanup_csr(GraphCSR *csr);
void cleanup_local_partition(LocalGraphPartition *lp);

// SSSP Update Algorithm Steps (adapted from paper)
void initialize_sssp(int source_global);
void apply_changes_step1(Change *changes, int num_changes, LocalGraphPartition *lp, int rank, int size);
void update_distances_step2(LocalGraphPartition *lp, int rank, int size);

// Result Gathering
void gather_and_print_results(LocalGraphPartition *lp, idx_t num_global_vertices, int rank, int size);

// Helper Functions
idx_t find_local_index(LocalGraphPartition *lp, idx_t global_id); // Potentially slow linear search
//int find_owner_rank(LocalGraphPartition *lp, idx_t global_id); // Placeholder
int find_owner_rank(int gid) {
    if (gid<0 || gid>=g_nvtxs) return -1;
    return partition_assignment[gid];
}
int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if (argc < 4) {
        if (world_rank == ROOT_RANK)
            fprintf(stderr, "Usage: %s <graph> <changes> <source>\n", argv[0]);
        MPI_Finalize(); return 1;
    }
    char *graph_file = argv[1];
    char *changes_file = argv[2];
    int source = atoi(argv[3]);

    // Rank 0: read CSR and partition
    if (world_rank == ROOT_RANK) {
        FILE *f = fopen(graph_file, "r");
        int u,v,w; int edge_cnt=0, max_id=0;
        while(fscanf(f, "%d %d %d", &u,&v,&w)==3) { edge_cnt++; max_id=max(max_id,u); max_id=max(max_id,v);} 
        rewind(f);
        g_nvtxs = max_id+1;
        // read edges into arrays
        int *eu = malloc(edge_cnt*sizeof(int));
        int *ev = malloc(edge_cnt*sizeof(int));
        int *ew = malloc(edge_cnt*sizeof(int));
        for(int i=0;i<edge_cnt;i++) fscanf(f, "%d %d %d", &eu[i],&ev[i],&ew[i]);
        fclose(f);
        // build CSR
        g_xadj = calloc(g_nvtxs+1,sizeof(int));
        for(int i=0;i<edge_cnt;i++) g_xadj[eu[i]+1]++;
        for(int i=1;i<=g_nvtxs;i++) g_xadj[i]+=g_xadj[i-1];
        g_adjncy = malloc(edge_cnt*sizeof(int));
        g_adjwgt = malloc(edge_cnt*sizeof(int));
        int *pos = calloc(g_nvtxs,sizeof(int));
        for(int i=0;i<edge_cnt;i++){
            int s=eu[i], idx=g_xadj[s]+pos[s]++;
            g_adjncy[idx]=ev[i]; g_adjwgt[idx]=ew[i];
        }
        free(eu); free(ev); free(ew); free(pos);
        // partition (round-robin for simplicity)
        partition_assignment = malloc(g_nvtxs*sizeof(int));
        for(int i=0;i<g_nvtxs;i++) partition_assignment[i] = i % world_size;
    }
    // Broadcast CSR and partition
    MPI_Bcast(&g_nvtxs,1,MPI_INT,ROOT_RANK,MPI_COMM_WORLD);
    if (world_rank!=ROOT_RANK) {
        g_xadj = malloc((g_nvtxs+1)*sizeof(int));
    }
    MPI_Bcast(g_xadj, g_nvtxs+1, MPI_INT, ROOT_RANK, MPI_COMM_WORLD);
    int nedges = g_xadj[g_nvtxs];
    if (world_rank!=ROOT_RANK) {
        g_adjncy = malloc(nedges*sizeof(int));
        g_adjwgt = malloc(nedges*sizeof(int));
        partition_assignment = malloc(g_nvtxs*sizeof(int));
    }
    MPI_Bcast(g_adjncy, nedges, MPI_INT, ROOT_RANK, MPI_COMM_WORLD);
    MPI_Bcast(g_adjwgt, nedges, MPI_INT, ROOT_RANK, MPI_COMM_WORLD);
    MPI_Bcast(partition_assignment, g_nvtxs, MPI_INT, ROOT_RANK, MPI_COMM_WORLD);

    distribute_graph();
    initialize_sssp(source);

    // Read changes
    FILE *cf = fopen(changes_file,"r"); char op[4]; int cu,cv,cw;
    while(fscanf(cf,"%s",op)==1){ if(strcmp(op,"INS")==0){ fscanf(cf,"%d %d %d",&cu,&cv,&cw);} else fscanf(cf,"%d %d",&cu,&cv);} fclose(cf);

    // Apply changes (simple rebuild distances/invoke dynamic update)
    bool local_changed;
    propagate_deletions_distributed();
    do {
        relax_edges_distributed(&local_changed);
        bool global_changed;
        MPI_Allreduce(&local_changed, &global_changed, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);
        local_changed = global_changed;
    } while(local_changed);

    gather_and_print_results();

    // Cleanup
    free(g_xadj); free(g_adjncy); free(g_adjwgt); free(partition_assignment);
    free(local_verts); free(local_map);
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


void distribute_graph() {
    // Build local and ghost lists
    local_count = ghost_count = 0;
    bool *is_ghost = calloc(g_nvtxs, sizeof(bool));
    for(int gid=0; gid<g_nvtxs; gid++){
        int owner = partition_assignment[gid];
        if (owner == world_rank) local_count++;
        else {
            // check if any neighbor owned locally
            for(int ei=g_xadj[gid]; ei<g_xadj[gid+1]; ei++){
                if (partition_assignment[g_adjncy[ei]] == world_rank) { is_ghost[gid]=true; break; }
            }
            if(is_ghost[gid]) ghost_count++;
        }
    }
    total_count = local_count + ghost_count;
    local_verts = calloc(total_count, sizeof(LocalVertex));
    local_map = malloc(total_count*sizeof(int));
    int li=0;
    // add locals
    for(int gid=0; gid<g_nvtxs; gid++){
        if(partition_assignment[gid]==world_rank){
            local_map[li]=gid;
            li++;
        }
    }
    // add ghosts
    for(int gid=0; gid<g_nvtxs; gid++){
        if(is_ghost[gid]){
            local_map[li]=gid;
            li++;
        }
    }
    free(is_ghost);
    // populate local vertices
    #pragma omp parallel for
    for(int i=0;i<total_count;i++){
        int gid = local_map[i];
        LocalVertex *v = &local_verts[i];
        v->global_id = gid;
        int start = g_xadj[gid], end = g_xadj[gid+1];
        v->num_edges = end-start;
        v->edges = malloc(v->num_edges * sizeof(Edge));
        for(int e=0;e<v->num_edges;e++){
            v->edges[e].target_global = g_adjncy[start+e];
            v->edges[e].weight = g_adjwgt[start+e];
        }
    }
}

void initialize_sssp(int source_global) {
    #pragma omp parallel for
    for(int i=0;i<total_count;i++){
        local_verts[i].distance = (local_verts[i].global_id==source_global?0:INF);
        local_verts[i].parent_global = -1;
        local_verts[i].affected = (local_verts[i].distance==0);
        local_verts[i].just_updated = false;
    }
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
void relax_edges_distributed(bool *local_changed) {
    *local_changed = false;
    bool *affected = malloc(total_count*sizeof(bool));
    for(int i=0;i<total_count;i++) affected[i] = local_verts[i].affected;

    #pragma omp parallel for reduction(||:*local_changed)
    for(int i=0;i<total_count;i++){
        if(!affected[i]) continue;
        LocalVertex *v = &local_verts[i];
        v->affected = false;
        int dv = v->distance;
        if(dv==INF) continue;
        for(int e=0;e<v->num_edges;e++){
            int ng = v->edges[e].target_global;
            int w = v->edges[e].weight;
            int owner = find_owner_rank(ng);
            // locate neighbor local
            for(int j=0;j<total_count;j++){
                if(local_verts[j].global_id==ng && owner==world_rank){
                    int nd = local_verts[j].distance;
                    int pd = dv + w;
                    if(pd < nd){
                        #pragma omp critical
                        {
                            local_verts[j].distance = pd;
                            local_verts[j].parent_global = v->global_id;
                            local_verts[j].affected = true;
                            local_verts[j].just_updated = true;
                            *local_changed = true;
                        }
                    }
                }
            }
        }
    }
    free(affected);
    // share updated distances across ranks (ghosts)
    int *all_dists = malloc(g_nvtxs*sizeof(int));
    int *my_dists = malloc(g_nvtxs*sizeof(int));
    for(int i=0;i<g_nvtxs;i++) my_dists[i] = INF;
    for(int i=0;i<total_count;i++) my_dists[local_verts[i].global_id] = local_verts[i].distance;
    MPI_Allreduce(my_dists, all_dists, g_nvtxs, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    for(int i=0;i<total_count;i++){
        local_verts[i].distance = all_dists[local_verts[i].global_id];
    }
    free(all_dists); free(my_dists);
}
void gather_and_print_results() {
    int *sendbuf = malloc(total_count*2*sizeof(int));
    int sc=0;
    for(int i=0;i<total_count;i++){
        int gid = local_verts[i].global_id;
        if(find_owner_rank(gid)==world_rank){
            sendbuf[sc++]=gid;
            sendbuf[sc++]=local_verts[i].distance;
        }
    }
    int *recvcounts = NULL, *displs = NULL;
    if(world_rank==ROOT_RANK){ recvcounts=malloc(world_size*sizeof(int)); displs=malloc(world_size*sizeof(int)); }
    MPI_Gather(&sc,1,MPI_INT,recvcounts,1,MPI_INT,ROOT_RANK,MPI_COMM_WORLD);
    int total=0;
    if(world_rank==ROOT_RANK){ displs[0]=0; total=recvcounts[0];
        for(int i=1;i<world_size;i++){ displs[i]=displs[i-1]+recvcounts[i-1]; total+=recvcounts[i]; }
    }
    int *recvbuf = NULL;
    if(world_rank==ROOT_RANK) recvbuf=malloc(total*sizeof(int));
    MPI_Gatherv(sendbuf, sc, MPI_INT, recvbuf, recvcounts, displs, MPI_INT, ROOT_RANK, MPI_COMM_WORLD);
    if(world_rank==ROOT_RANK){
        printf("--- Final SSSP Results ---\nGlobalID  Distance\n");
        for(int i=0;i<total;i+=2) printf("%8d  %8d\n", recvbuf[i], recvbuf[i+1]);
    }
    free(sendbuf);
    if(world_rank==ROOT_RANK){ free(recvcounts); free(displs); free(recvbuf); }
}
