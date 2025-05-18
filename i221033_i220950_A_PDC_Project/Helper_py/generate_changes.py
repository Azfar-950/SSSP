import random
import os

# Define input and output file paths
input_file = 'graph.txt'
output_file = 'changes.txt'

# Parameters
total_changes = 1000  # Total number of changes (deletions + insertions)
deletion_percentage = 30  # Percentage of changes that are deletions (0 to 100)

# Validate parameters
if total_changes < 0:
    raise ValueError("Total changes must be non-negative")
if not 0 <= deletion_percentage <= 100:
    raise ValueError("Deletion percentage must be between 0 and 100")

# Calculate number of deletions and insertions
num_deletions = int(total_changes * (deletion_percentage / 100))
num_insertions = total_changes - num_deletions

# Ensure input file exists
if not os.path.exists(input_file):
    raise FileNotFoundError(f"Input file not found at: {input_file}")

# Read edges and nodes from input file
edges = set()
nodes = set()
with open(input_file, 'r') as f:
    for line in f:
        try:
            u, v, w = map(int, line.strip().split())
            edges.add((min(u, v), max(u, v)))  # Store as sorted tuple to avoid duplicates
            nodes.add(u)
            nodes.add(v)
        except ValueError:
            print(f"Skipping invalid line: {line.strip()}")

# Validate number of deletions
if num_deletions > len(edges):
    print(f"Warning: Requested {num_deletions} deletions, but only {len(edges)} edges exist. Adjusting to {len(edges)}.")
    num_deletions = len(edges)
    num_insertions = total_changes - num_deletions

# Select edges to delete
deletions = random.sample(list(edges), num_deletions)

# Generate new edges for insertion
insertions = []
nodes_list = list(nodes)
existing_edges = set(edges)  # Copy for checking
attempts = 0
max_attempts = num_insertions * 100  # Prevent infinite loops
while len(insertions) < num_insertions and attempts < max_attempts:
    u = random.choice(nodes_list)
    v = random.choice(nodes_list)
    if u != v:  # Avoid self-loops
        edge = (min(u, v), max(u, v))
        if edge not in existing_edges:
            weight = random.randint(1, 10)  # Random weight between 1 and 10
            insertions.append((u, v, weight))
            existing_edges.add(edge)  # Prevent duplicates
    attempts += 1

if len(insertions) < num_insertions:
    print(f"Warning: Could only generate {len(insertions)} insertions due to limited new edge possibilities.")

# Write changes to output file
with open(output_file, 'w') as f:
    # Write deletions
    for u, v in deletions:
        f.write(f"DEL {u} {v}\n")
    # Write insertions
    for u, v, w in insertions:
        f.write(f"INS {u} {v} {w}\n")

print(f"Generated {len(deletions)} deletions and {len(insertions)} insertions in {output_file}")