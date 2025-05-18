import csv
import random
import os

# Define input and output file paths
input_path = os.path.join('gemsec_facebook_dataset', 'facebook_clean_data', 'tvshow_edges.csv')
output_file = 'graph.txt'

# Ensure input file exists
if not os.path.exists(input_path):
    raise FileNotFoundError(f"Input file not found at: {input_path}")

# Read the CSV file and write to a text file with random weights
with open(input_path, 'r') as csv_file, open(output_file, 'w') as txt_file:
    # Initialize CSV reader
    reader = csv.reader(csv_file)
    # Skip header row
    next(reader)
    # Process each edge
    for row in reader:
        if len(row) >= 2:  # Ensure row has at least two columns
            node_1, node_2 = row[0], row[1]
            # Generate random weight between 1 and 10
            weight = random.randint(1, 10)
            # Write to text file in the format: node_1 node_2 weight
            txt_file.write(f"{node_1} {node_2} {weight}\n")