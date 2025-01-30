#Some notes

A. Time-Dependent Birth-Death Models
You can modify speciation (λ) and extinction (μ) rates over geological time.
This can reflect tectonic shifts, sea level changes, or climate events.
Tools:
##R: TreeSim, TESS
##Python: DendroPy

import numpy as np
import dendropy

# Function to simulate a phylogenetic tree with divergence events due to barriers
def simulate_barrier_phylogeny(num_species=10, divergence_time=5.0, barrier_effect=0.5):
    """
    Simulates a phylogenetic tree with species divergence influenced by geographic barriers.
    
    Parameters:
        num_species (int): Number of species to simulate.
        divergence_time (float): Time when the barrier appears, leading to divergence.
        barrier_effect (float): Probability of divergence when encountering a barrier.
    
    Returns:
        DendroPy Tree object.
    """
    tree = dendropy.Tree()  # Create an empty tree
    root = tree.seed_node  # Root node
    
    # Create initial species lineages
    for i in range(num_species):
        new_node = root.new_child(taxon=dendropy.Taxon(label=f"Species_{i+1}"))
    
    # Introduce a barrier-driven divergence event
    for node in tree.leaf_nodes():
        if np.random.rand() < barrier_effect:
            # Create a new lineage due to barrier-driven divergence
            new_child = node.new_child(taxon=dendropy.Taxon(label=f"{node.taxon.label}_Diverged"))
            new_child.edge.length = divergence_time  # Set divergence time

    return tree

# Simulate the tree with 10 species and a barrier effect at time 5.0
simulated_tree = simulate_barrier_phylogeny(num_species=10, divergence_time=5.0, barrier_effect=0.7)

# Convert the tree to Newick format for visualization
tree_newick = simulated_tree.as_string(schema="newick")
print(tree_newick)
