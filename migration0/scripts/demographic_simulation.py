# Phylogeny simulation without gene flow in Python
# We will create a three-taxon demographic scenario with AB and C splitting 1 Ma and A and B splitting 3 Ma in order to generate sequences data
# Phylogeny created from selecting taxa in the simulated dataset (generates for n=60, select 1 randomly from each population)

# General setup:
# conda activate msprime-env
# python
import msprime
import matplotlib.pyplot as plt
import demesdraw

# Define generation time (in years per generation)
generation_time = 1

# Define the demographic model
demography = msprime.Demography()
demography.add_population(name="A", initial_size=50)
demography.add_population(name="B", initial_size=50)
demography.add_population(name="C", initial_size=100)
demography.add_population(name="AB", initial_size=100)
demography.add_population(name="ABC", initial_size=200)
demography.add_population_split(time=1000000, derived=["A", "B"], ancestral="AB")
demography.add_population_split(time=3000000, derived=["AB", "C"], ancestral="ABC")

# Simulate tree sequence
ts = msprime.sim_ancestry(samples={"A": 10, "B": 10, "C": 10}, 
                          demography=demography, 
                          sequence_length=10000,  #Adjust to a higher length if computationally feasbile in later steps
                          recombination_rate=4.8e-8)

# Define the mutation rate and HKY model parameters
mutation_rate = 1e-8
kappa = 2.0  # Transition/transversion ratio for HKY model
equilibrium_frequencies = [0.25, 0.25, 0.25, 0.25]  #A, C, G, T frequencies

# Apply HKY mutations to the tree sequence
ts = msprime.sim_mutations(
    ts,
    rate=mutation_rate, keep=True,
    model=msprime.HKY(kappa=kappa, equilibrium_frequencies=equilibrium_frequencies)
)

# Function to get MRCA for each population
def get_population_mrca(ts, pop_name):
    pop_id = next((i.id for i in ts.populations() if i.metadata["name"] == pop_name), None)
    if pop_id is None:
        raise ValueError(f"Population {pop_name} not found in tree sequence.")
    sample_nodes = [n for n in ts.samples() if ts.node(n).population == pop_id]
    if len(sample_nodes) < 2:
        raise ValueError(f"Not enough samples for population {pop_name} to compute MRCA.")
    return ts.first().mrca(*sample_nodes)

# Get MRCA nodes for A, B, and C
pop_mrcas = {}
for pop in ["A", "B", "C"]:
    try:
        pop_mrcas[pop] = get_population_mrca(ts, pop)
    except ValueError as e:
        print(e)

# Compute MRCA for AB and ABC
tree = ts.first()
try:
    pop_mrcas["AB"] = tree.mrca(pop_mrcas["A"], pop_mrcas["B"])
    pop_mrcas["ABC"] = tree.mrca(pop_mrcas["AB"], pop_mrcas["C"])
except KeyError as e:
    print(f"Error computing MRCA: {e}")

print("Population MRCAs:", pop_mrcas)

# Get population-level Newick tree
pop_newick = tree.as_newick(root=pop_mrcas["ABC"])
print("Population-Level Newick Tree:", pop_newick)

# Simplify the tree sequence to retain only species-level structure
simplified_ts = ts.simplify(samples=[pop_mrcas["A"], pop_mrcas["B"], pop_mrcas["C"]], keep_input_roots=True)

# Get the first tree from the simplified tree sequence
collapsed_tree = simplified_ts.first()

# Convert to Newick format
collapsed_newick = collapsed_tree.as_newick()
print("Species-Level Newick Tree:", collapsed_newick)

# Save the Newick tree to file
with open("collapsed_newick.tre", "w") as f:
    f.write(collapsed_newick)
    
# Write the full Nexus format sequence data
with open("simulated_sequences.nex", "w") as nexus_file:
    ts.write_nexus(nexus_file)

# Simplify the tree sequence
simplified_ts_single = ts.simplify(samples=[pop_mrcas["A"], pop_mrcas["B"], pop_mrcas["C"]], keep_input_roots=True)

# Write the reduced dataset to another Nexus file
with open("simulated_sequences_single.nex", "w") as nexus_file:
    simplified_ts_single.write_nexus(nexus_file)
    
# Draw demographic model
fig, ax = plt.subplots(figsize=(6, 4))
demesdraw.tubes(demography.to_demes(), ax=ax)
plt.show()

print("Simulated sequence data has been written to:")
print("- Full dataset: 'simulated_sequences.nex'")
print("- Single-individual dataset: 'simulated_sequences_single.nex'")

# Reformat Nexus file correct format for later analysis
input_file = "simulated_sequences_single.nex"
output_file = "modified_sequences.nex"

with open(input_file, "r") as infile:
    lines = infile.readlines()

# Extract the matrix content and remove anything after it
start_writing = False
modified_lines = ["BEGIN DATA;\n", "DIMENSIONS  NTAX=3 NCHAR=10000;\n", "FORMAT DATATYPE=DNA GAP=- MISSING=?;\n", "MATRIX\n"]

for line in lines:
    if "MATRIX" in line:
        start_writing = True
    if start_writing:
        modified_lines.append(line)
    if ";" in line and start_writing:  # Stop at the end of MATRIX block
        break

# Write to a new file
with open(output_file, "w") as outfile:
    outfile.writelines(modified_lines)

print(f"Modified NEXUS file saved as {output_file}")

# NOTE: right now you have to edit the newick tree by hand using a text editor like nano; simply move the decimal place so that the tree is scaled down by a factor of a million (eg: 3.0 instead of 3,000,000)
