import numpy as np

# Parameters
n = 100  # Number of nodes
avg_degree = 2  # Average degree
p = avg_degree / (n - 1)  # Edge probability for Erdos-Renyi
num_changes = 30  # Number of edge changes
output_graph = "graph.metis"
output_changes = "changes.txt"

# Generate random graph
np.random.seed(42)
edges = []
for i in range(n):
    for j in range(i + 1, n):
        if np.random.random() < p:
            weight = np.random.randint(1, 100)  # Random weight between 1 and 99
            edges.append((i + 1, j + 1, weight))  # 1-based indexing for METIS

# Write graph to METIS format
with open(output_graph, "w") as f:
    # Header: number of vertices, number of edges
    f.write(f"{n} {len(edges)}\n")
    # Adjacency list
    adj_list = [[] for _ in range(n)]
    for u, v, w in edges:
        adj_list[u - 1].append(f"{v} {w}")
        adj_list[v - 1].append(f"{u} {w}")
    for adj in adj_list:
        f.write(f"{' '.join(adj)}\n")

# Generate changes file
changes = []
existing_edges = set((min(u, v), max(u, v)) for u, v, _ in edges)
non_edges = [(i, j) for i in range(1, n + 1) for j in range(i + 1, n + 1) if (i, j) not in existing_edges]
np.random.shuffle(non_edges)

# 50% insertions
for _ in range(num_changes // 2):
    if non_edges:
        u, v = non_edges.pop()
        w = np.random.randint(1, 100)
        changes.append(f"I {u} {v} {w}")

# 50% deletions
delete_edges = list(existing_edges)
np.random.shuffle(delete_edges)
for i in range(num_changes // 2):
    if i < len(delete_edges):
        u, v = delete_edges[i]
        w = next(w for x, y, w in edges if (x == u and y == v) or (x == v and y == u))
        changes.append(f"D {u} {v} {w}")

# Write changes file
with open(output_changes, "w") as f:
    for change in changes:
        f.write(f"{change}\n")

print(f"Generated graph: {output_graph}")
print(f"Generated changes: {output_changes}")