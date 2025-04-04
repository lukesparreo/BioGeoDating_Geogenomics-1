import numpy as np

# Load original NEXUS file
with open("modified_sequences.nex", "r") as f:
    lines = f.readlines()

# Locate the matrix section
start_index = next(i for i, line in enumerate(lines) if "matrix" in line.lower()) + 1
end_index = next(i for i, line in enumerate(lines) if ";" in line and i > start_index)

# Parse sequence data
matrix_lines = [line.strip() for line in lines[start_index:end_index] if line.strip() and not line.strip().startswith("[")]
sequences = []
sample_names = []

for line in matrix_lines:
    parts = line.split()
    if len(parts) >= 2:
        sample_names.append(parts[0])
        sequences.append(parts[1])

# Convert sequences to numpy array
seq_array = np.array([list(seq) for seq in sequences])
num_samples, num_sites = seq_array.shape

# Replace all "?" in fully missing columns with a random invariant base
nucleotides = ["A", "C", "G", "T"]
for col in range(num_sites):
    if np.all(seq_array[:, col] == "?"):
        seq_array[:, col] = np.random.choice(nucleotides)

# Rebuild the NEXUS file with filled-in sequences
output_lines = lines[:start_index]
for name, seq in zip(sample_names, seq_array):
    output_lines.append(f"{name}    {''.join(seq)}\n")
output_lines.append("  ;\nEND;\n")

# Save to a new NEXUS file
with open("modified_sequences_filled.nex", "w") as f:
    f.writelines(output_lines)