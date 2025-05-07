import pandas as pd
import random

# Read the CSV file
df = pd.read_csv('politician_edges.csv')

# Create a list to store edges with weights
weighted_edges = []

# Generate random weights for each edge
for index, row in df.iterrows():
    node1 = row['node_1']
    node2 = row['node_2']
    weight = random.randint(1, 10)  # Assign random weight between 1 and 10
    weighted_edges.append((node1, node2, weight))

# Convert to DataFrame for easier handling
weighted_df = pd.DataFrame(weighted_edges, columns=['node1', 'node2', 'weight'])

# Save to a text file with space-separated values
weighted_df.to_csv('weighted_politician_edges.txt', sep=' ', index=False, header=False)

print("Weighted graph has been saved to 'weighted_politician_edges.txt'")