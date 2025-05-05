// generate_graph.c
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define NUM_NODES 100
#define EDGE_FACTOR 5 // Average number of edges per node
#define MAX_WEIGHT 100
#define NUM_CHANGES 20 // Number of updates

int main() {
    srand(time(NULL));

    // --- Generate Graph File (graph.txt) ---
    FILE *graph_file = fopen("graph.txt", "w");
    if (!graph_file) {
        perror("Error opening graph.txt");
        return 1;
    }

    int num_edges = 0;
    // Header (optional, depends on how you read it): num_nodes num_edges
    // We'll write edges first, then update the header if needed, or just omit it.
    // fprintf(graph_file, "%d %d\n", NUM_NODES, 0); // Placeholder count

    printf("Generating graph with %d nodes...\n", NUM_NODES);
    for (int i = 0; i < NUM_NODES; ++i) {
        // Add some edges for each node
        int edges_for_node = rand() % (2 * EDGE_FACTOR) + 1;
        for (int j = 0; j < edges_for_node; ++j) {
            int neighbor = rand() % NUM_NODES;
            if (neighbor != i) { // Avoid self-loops for simplicity
                int weight = (rand() % MAX_WEIGHT) + 1;
                fprintf(graph_file, "%d %d %d\n", i, neighbor, weight);
                // For undirected, add the reverse edge? Assume directed for now based on paper's SSSP logic.
                // If undirected, you might add: fprintf(graph_file, "%d %d %d\n", neighbor, i, weight);
                num_edges++;
            }
        }
    }
    fclose(graph_file);
    printf("Generated graph.txt with %d edges.\n", num_edges);


    // --- Generate Changes File (changes.txt) ---
    FILE *changes_file = fopen("changes.txt", "w");
     if (!changes_file) {
        perror("Error opening changes.txt");
        return 1;
    }
    printf("Generating changes file with %d operations...\n", NUM_CHANGES);
    for (int i = 0; i < NUM_CHANGES; ++i) {
        int op_type = rand() % 2; // 0 for DEL, 1 for INS
        int u = rand() % NUM_NODES;
        int v = rand() % NUM_NODES;

        if (u == v) continue; // Skip self-loops in changes

        if (op_type == 0) { // Deletion
             // Note: Deleting a non-existent edge is handled by the algorithm check (line 5, Algo 2)
            fprintf(changes_file, "DEL %d %d\n", u, v);
        } else { // Insertion
            int weight = (rand() % MAX_WEIGHT) + 1;
            fprintf(changes_file, "INS %d %d %d\n", u, v, weight);
        }
    }
    fclose(changes_file);
    printf("Generated changes.txt.\n");

    return 0;
}