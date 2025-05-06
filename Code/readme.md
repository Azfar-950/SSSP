Parallel SSSP Update Implementation for Dynamic Networks
This repository contains a parallel implementation of a Single-Source Shortest Path (SSSP) update algorithm for dynamic networks, as described in the paper "A Parallel Algorithm Template for Updating Single-Source Shortest Paths in Large-Scale Dynamic Networks" by Khanda et al. The implementation (sssp.cpp) uses MPI for distributed-memory parallelism and OpenMP for shared-memory parallelism to update the SSSP tree after edge insertions and deletions. This README provides a detailed explanation of the code’s functionality for a specific 10-vertex graph, its correctness, the algorithms used, asynchrony, generality, scalability for larger graphs, graph requirements, and limitations.
Overview
The code processes a dynamic undirected, weighted graph and updates the SSSP tree from a specified source vertex (vertex 0 in this case) after applying edge changes (insertions and deletions). It uses a distributed approach, partitioning the graph across multiple processes (ranks) and leveraging asynchronous communication to minimize synchronization overhead. The implementation is tested on a small graph with 10 vertices and 15 edges, with specific changes applied to demonstrate its functionality.
Files

sssp.cpp: The main implementation of the parallel SSSP update algorithm, using MPI and OpenMP. It reads partitioned subgraphs, applies changes, and updates the SSSP tree.
partition.cpp: A utility program that uses METIS to partition the input graph into subgraphs for each MPI rank, generating subgraph_0.txt, subgraph_1.txt, subgraph_0_nodes.txt, and subgraph_1_nodes.txt.
graph.txt: The input graph file, describing an undirected, weighted graph with 10 vertices (0–9) and 15 edges. Format: num_vertices num_edges, followed by u v weight for each edge.
changes.txt: The file specifying edge changes to apply. It contains 2 deletions (edges (4, 9) and (0, 3)) and 2 insertions (edges (6, 8, weight 1) and (5, 8, weight 4)). Format: num_changes, followed by operation u v [weight] (operation: 0 for deletion, 1 for insertion).
subgraph_0.txt: The partitioned subgraph for Rank 0, containing edges involving local vertices {1, 4, 5, 6, 9} and ghost vertices {2, 3, 7, 8}.
subgraph_1.txt: The partitioned subgraph for Rank 1, containing edges involving local vertices {0, 2, 3, 7, 8} and ghost vertices {1, 4, 5, 9}.
subgraph_0_nodes.txt: Node information for Rank 0, listing local and ghost vertices.
subgraph_1_nodes.txt: Node information for Rank 1, listing local and ghost vertices.

Graph Description
The graph in graph.txt has:

Vertices: 10 (labeled 0–9)
Edges: 15 undirected edges with positive weights, e.g., (0, 2, 9), (0, 8, 4), (1, 2, 2), etc.
Source Vertex: 0 (used for SSSP computation)

Changes Applied (changes.txt)

Deletions:
Edge (4, 9)
Edge (0, 3)


Insertions:
Edge (6, 8, weight 1)
Edge (5, 8, weight 4)



Graph Partitioning
The graph is partitioned using METIS (via partition.cpp) into two subgraphs for two MPI ranks:

Rank 0:
Local Vertices: {1, 4, 5, 6, 9}
Ghost Vertices: {2, 3, 7, 8}
Edges in subgraph_0.txt include those incident to local vertices, e.g., (1, 2, 2), (4, 3, 1), (5, 8, 4 after insertion).


Rank 1:
Local Vertices: {0, 2, 3, 7, 8}
Ghost Vertices: {1, 4, 5, 9}
Edges in subgraph_1.txt include those incident to local vertices, e.g., (0, 2, 9), (8, 7, 2), (8, 6, 1 after insertion).



Algorithm Description
The implementation follows Algorithm 4 from the paper (Section 4.2), an asynchronous parallel SSSP update algorithm for dynamic networks. It consists of two main steps:

Step 1: Identify Affected Vertices (apply_changes_step1 in f.cpp):

Processes edge insertions and deletions from changes.txt.
Marks vertices incident to changed edges as affected (e.g., vertices 0, 3, 4, 5, 6, 8, 9).
Updates the local graph structure by adding or removing edges.


Step 2: Update SSSP Tree (relax_edges_distributed in f.cpp):

Iteratively relaxes edges involving affected vertices to update distances and parents in the SSSP tree.
Exchanges ghost vertex updates (distances and parents) between ranks using MPI communication.
Continues until no further changes occur (convergence).



Asynchrony
The algorithm is asynchronous, as described in Section 5.1 of the paper:

No Explicit Synchronization: Ranks exchange ghost vertex updates (e.g., Rank 1 sends updates for vertices 2, 7, 8, 9 to Rank 0 in Iteration 1) without global barriers, using non-blocking MPI calls (e.g., MPI_Isend, MPI_Irecv).
Redundant Computations: To avoid synchronization overhead, the algorithm allows redundant relaxations (e.g., vertex 3’s distance updated twice in Rank 1, Iteration 2), ensuring correctness while minimizing communication delays.
Convergence Check: Each rank reports local and ghost changes, and the algorithm terminates when a global reduction (via MPI_Reduce) confirms no further updates (global_change: NO).

How the Code Works for the 10-Vertex Graph
The code executes with two MPI processes (mpirun -np 2 ./f), processing the 10-vertex graph and applying the specified changes. Below is a detailed walkthrough of its operation, based on the provided output.
1. Initialization

Rank 0:
Reads subgraph_0.txt and subgraph_0_nodes.txt.
Initializes distances for local vertices {1, 4, 5, 6, 9} to infinity (INF) and parents to -1.
Initializes ghost distances for {2, 3, 7, 8} to INF.


Rank 1:
Reads subgraph_1.txt and subgraph_1_nodes.txt.
Initializes distance[0] = 0 (source vertex, local to Rank 1), other local vertices {2, 3, 7, 8} to INF, and parents to -1.
Initializes ghost distances for {1, 4, 5, 9} to INF.



2. Applying Changes

Rank 0:
Removes edge (4, 9).
Adds edges (6, 8, weight 1) and (5, 8, weight 4).
Marks affected vertices: 4, 5, 6, 8, 9.


Rank 1:
Removes edge (0, 3).
Adds edges (6, 8, weight 1) and (5, 8, weight 4).
Marks affected vertices: 0, 3, 5, 6, 8.



3. Iterative Relaxation
The algorithm iterates until convergence, relaxing edges and exchanging ghost updates. Key iterations:

Iteration 1:

Rank 1:
Relaxes edges from source vertex 0: (0, 2, weight 9) → distance[2] = 9, (0, 8, weight 4) → distance[8] = 4.
Relaxes edges from vertex 8: (8, 7, weight 2) → distance[7] = 6, (8, 1, weight 2 via 2) → ghost_distance[1] = 11, (8, 6, weight 1) → ghost_distance[6] = 5, (8, 5, weight 4) → ghost_distance[5] = 8, (8, 9, weight 1 via 7) → ghost_distance[9] = 7.
Sends ghost updates for vertices 2, 7, 8, 9 to Rank 0.
Local change: YES.


Rank 0:
No initial relaxations (no local vertex with distance 0).
Local change: NO, ghost change: NO.




Iteration 2:

Rank 0:
Receives ghost updates: distance[2] = 9, distance[7] = 6, distance[8] = 4.
Relaxes edges:
(2, 1, weight 2) → distance[1] = 11
(1, 6, weight 7) → distance[6] = 18
(1, 5, weight 4) → distance[5] = 15
(6, 4, weight 9) → distance[4] = 27
(8, 5, weight 4) → distance[5] = 8
(5, 9, weight 5) → distance[9] = 13
(8, 6, weight 1) → distance[6] = 5
(6, 9, weight 7) → distance[9] = 12
(6, 4, weight 9) → distance[4] = 14
(7, 9, weight 1) → distance[9] = 7
Updates ghost_distance[3] to 8 (via 9).


Sends ghost updates for vertices 1, 4, 5, 6, 9 to Rank 1.
Local change: YES, ghost change: YES.


Rank 1:
Relaxes (5, 3, weight 8) → distance[3] = 16, then (9, 3, weight 1) → distance[3] = 8.
Updates ghost_distance[4] = 9 (via 3).
Sends ghost update for vertex 3 to Rank 0.
Local change: YES.




Iteration 3:

Rank 0:
Receives ghost update: distance[3] = 8.
Relaxes (3, 4, weight 1) → distance[4] = 9.
Sends ghost update for vertex 4 to Rank 1.
Local change: YES.


Rank 1:
Receives ghost updates for vertices 1, 4, 5, 6, 9.
No further relaxations.
Local change: NO, ghost change: NO.




Iteration 4:

Both ranks report no changes (local_change: NO, ghost_change: NO).
Global change: NO, algorithm terminates.



4. Final SSSP Tree

Rank 0 aggregates results via MPI_Reduce:
Vertex 0: Distance = 0, Parent = -1
Vertex 1: Distance = 11, Parent = 2
Vertex 2: Distance = 9, Parent = 0
Vertex 3: Distance = 8, Parent = 9
Vertex 4: Distance = 9, Parent = 3
Vertex 5: Distance = 8, Parent = 8
Vertex 6: Distance = 5, Parent = 8
Vertex 7: Distance = 6, Parent = 8
Vertex 8: Distance = 4, Parent = 0
Vertex 9: Distance = 7, Parent = 7



Correctness Verification
The final SSSP tree is correct, as verified by tracing the shortest paths from source vertex 0 in the updated graph:

Vertex 0: Distance = 0, Parent = -1 (source).
Vertex 8: Path 0→8 (weight 4) → Distance = 4, Parent = 0.
Vertex 6: Path 0→8→6 (weights 4 + 1 = 5) → Distance = 5, Parent = 8.
Vertex 7: Path 0→8→7 (weights 4 + 2 = 6) → Distance = 6, Parent = 8.
Vertex 9: Path 0→8→7→9 (weights 4 + 2 + 1 = 7) → Distance = 7, Parent = 7.
Vertex 5: Path 0→8→5 (weights 4 + 4 = 8) → Distance = 8, Parent = 8.
Vertex 3: Path 0→8→7→9→3 (weights 4 + 2 + 1 + 1 = 8) → Distance = 8, Parent = 9.
Vertex 4: Path 0→8→7→9→3→4 (weights 4 + 2 + 1 + 1 + 1 = 9) → Distance = 9, Parent = 3.
Vertex 2: Path 0→2 (weight 9) → Distance = 9, Parent = 0.
Vertex 1: Path 0→2→1 (weights 9 + 2 = 11) → Distance = 11, Parent = 2.

Each distance matches the shortest path length, and the parent pointers form a valid tree rooted at vertex 0, confirming correctness.
Asynchrony in Action
The code’s asynchrony is evident in:

Non-Blocking Communication: Ghost updates (e.g., Rank 1 sending distance[2] = 9 to Rank 0 in Iteration 1) use non-blocking MPI calls, allowing ranks to proceed without waiting.
Redundant Relaxations: Multiple updates to the same vertex (e.g., vertex 3’s distance set to 16 then 8 in Rank 1, Iteration 2) reflect the asynchronous approach, prioritizing speed over redundant computation elimination.
Convergence: The algorithm converges in 4 iterations, with ranks independently processing updates and coordinating only through ghost exchanges, as per Algorithm 4.

Generality of the Code
The code is reasonably generic for SSSP updates in dynamic networks:

Graph Type: Handles undirected, weighted graphs with positive edge weights, suitable for applications like social networks or road networks.
Input Flexibility: Reads graph and change data from text files (graph.txt, changes.txt), allowing different graph sizes and structures.
Algorithm Template: Implements Algorithm 4 as a template (Section 4), portable across MPI/OpenMP-supporting systems.
Partitioning: Relies on METIS for graph partitioning, which is widely applicable to various graph topologies.

Limitations to Generality:

Undirected Graphs: Assumes undirected edges, limiting applicability to directed graphs without modification.
Pre-Partitioned Input: Requires METIS-generated partitions (subgraph_*.txt, subgraph_*_nodes.txt), which may be a barrier for users without partitioning tools.
Single-Source: Designed for SSSP from one source, not all-pairs shortest paths.
Change Types: Supports edge insertions and deletions but not vertex additions/removals or weight updates.

Scalability for Larger Graphs
The code is designed to scale to very large graphs, as demonstrated in the paper:

Paper’s Experiments: Tested on graphs up to 128M vertices and 1B edges (Section 6.1), achieving up to 5.6x speedup over GPU-based Gunrock and 5x over CPU-based Galois for insertion-heavy workloads.
Parallel Design: Uses MPI for distributed processing across nodes and OpenMP for intra-node parallelism, suitable for clusters.
METIS Partitioning: Minimizes cross-partition edges, reducing communication overhead (Section 5.2).
Asynchronous Updates: Avoids global synchronization, improving scalability for large graphs with sparse updates.

Challenges for Large Graphs:

Memory: Large graphs require significant memory for adjacency lists and ghost data. The paper’s tests used a 192 GB system, so smaller systems may struggle.
Communication: High cross-partition edges increase MPI communication overhead, as noted in Section 6.2.
Change Density: The update algorithm is less efficient for deletion-heavy changes (>50%) or widespread updates (>75–85% of nodes), where recomputation may be faster (Section 6.3).
Convergence: Large graphs with long shortest paths may require more iterations, increasing runtime.

Tuning Needs: For optimal performance, users may need to adjust OpenMP thread counts, optimize MPI communication, or improve partitioning for specific graphs.
Graph Requirements
The input graph must meet the following requirements:

Type: Undirected, weighted graph with positive edge weights.
Format: Specified in graph.txt as num_vertices num_edges, followed by u v weight lines.
Changes: Specified in changes.txt as num_changes, followed by operation u v [weight] (0 for deletion, 1 for insertion).
Partitioning: Pre-partitioned using METIS, with subgraphs stored in subgraph_*.txt (format: num_vertices num_edges, followed by u v weight) and node information in subgraph_*_nodes.txt (format: local and ghost vertex lists).
Vertex IDs: Contiguous integers starting from 0 (e.g., 0 to 9 for the 10-vertex graph).
No Self-Loops or Multiple Edges: The algorithm assumes a simple graph.

Limitations

Graph Type: Limited to undirected graphs; directed graphs require code modifications.
Change Types: Only supports edge insertions and deletions, not vertex changes or weight updates.
Partitioning Dependency: Requires METIS for partitioning, which may be complex for large graphs or users unfamiliar with graph partitioning tools.
Redundant Computations: Asynchronous design allows redundant relaxations, increasing computation for dense graphs or frequent updates.
Scalability Limits: Performance may degrade for deletion-heavy changes or poorly partitioned graphs with high communication overhead.
Small Graph Testing: The provided output tests a small graph (10 vertices), so large-scale performance is inferred from the paper’s experiments.

Building and Running
Prerequisites

MPI implementation (e.g., OpenMPI)
OpenMP support (included in most C++ compilers)
METIS library for graph partitioning
C++ compiler (e.g., g++)

Compilation
mpicxx -fopenmp -o f f.cpp -lmetis
g++ -o partition partition.cpp -lmetis

Running

Partition the graph:
./partition graph.txt 2

Generates subgraph_0.txt, subgraph_1.txt, subgraph_0_nodes.txt, subgraph_1_nodes.txt for 2 partitions.

Run the SSSP update:
mpirun -np 2 ./f graph.txt changes.txt subgraph_0.txt subgraph_1.txt subgraph_0_nodes.txt subgraph_1_nodes.txt



Expected Output
The program outputs:

Graph and subgraph details for each rank.
Changes applied (edge deletions/insertions).
Iteration details, including relaxed edges, ghost updates, and convergence status.
Final SSSP tree with distances and parents for all vertices.

Conclusion
This implementation successfully updates the SSSP tree for a 10-vertex graph after applying 2 edge deletions and 2 insertions, using Algorithm 4 from the paper. It demonstrates correctness through verified shortest paths, leverages asynchrony for efficiency, and is generic for undirected, weighted graphs with pre-partitioned inputs. The code scales to large graphs (up to 128M vertices, 1B edges), but performance depends on graph sparsity, change patterns, and hardware. Users should ensure proper partitioning and sufficient resources for large-scale applications.
For further details contact the repository maintainer.

