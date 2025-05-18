import networkx as nx
import matplotlib.pyplot as plt
import os

# Try importing pygraphviz for graphviz_layout
try:
    from networkx.drawing.nx_agraph import graphviz_layout
    use_graphviz = True
except ImportError:
    use_graphviz = False
    print("Pygraphviz not installed. Falling back to spring_layout.")

# Define input and output paths
input_file = "sssp_output.txt"
output_dir = os.path.join("gemsec_facebook_dataset", "facebook_clean_data")
output_image = os.path.join(output_dir, "sssp_graph.png")

# Parameters for large graph
max_distance = 50  # Only include nodes with distance <= max_distance (set to None for all nodes)
label_nodes = False  # Set to True to label nodes (may clutter large graphs)

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Check if input file exists
if not os.path.exists(input_file):
    raise FileNotFoundError(f"Input file not found: {input_file}")

# Create a directed graph
G = nx.DiGraph()

# Read SSSP output
vertices = []
edges = []
with open(input_file, 'r') as f:
    # Skip header
    next(f)
    for line in f:
        try:
            vertex, distance, parent = map(int, line.strip().split())
            if max_distance is None or distance <= max_distance:
                vertices.append((vertex, distance))
                if parent != -1 and parent != vertex:  # Exclude source's self-loop or invalid parent
                    edges.append((parent, vertex))
        except ValueError:
            print(f"Skipping invalid line: {line.strip()}")

# Add nodes with labels (vertex:distance)
for vertex, distance in vertices:
    G.add_node(vertex, label=f"{vertex}:{distance}")

# Add edges based on parent relationships
G.add_edges_from(edges)

# Set up the plot
plt.figure(figsize=(15, 12))

# Choose layout
if use_graphviz:
    pos = graphviz_layout(G, prog="dot")  # Hierarchical layout
else:
    pos = nx.spring_layout(G, seed=42, k=1.0/len(G)**0.5)  # Adjust k for node separation

# Draw nodes
nx.draw_networkx_nodes(G, pos, node_color="lightblue", node_size=50, alpha=0.7)

# Draw edges
nx.draw_networkx_edges(G, pos, edge_color="gray", arrows=True, arrowsize=10, width=0.5)

# Draw labels (optional)
if label_nodes:
    labels = nx.get_node_attributes(G, 'label')
    nx.draw_networkx_labels(G, pos, labels, font_size=6)

# Set title and layout
plt.title(f"SSSP Shortest Path Tree (Max Distance: {max_distance if max_distance else 'All'})")
plt.tight_layout()

# Save and close
plt.savefig(output_image, format="png", dpi=300, bbox_inches="tight")
plt.close()
print(f"Graph saved as {output_image}")