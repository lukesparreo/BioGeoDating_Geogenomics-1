# Notes on possible steps in the project

Time-Dependent Birth-Death Models:

You can modify speciation (λ) and extinction (μ) rates over geological time.
This can reflect tectonic shifts, sea level changes, or climate events.

## Function to simulate a phylogenetic tree with divergence events due to barriers with DendroPy

    import numpy as np

    import dendropy


    def simulate_barrier_phylogeny(num_species=10, divergence_time=5.0, barrier_effect=0.5)
    
Simulates a phylogenetic tree with species divergence influenced by geographic barriers.
    
Parameters:
        num_species (int): Number of species to simulate.
        divergence_time (float): Time when the barrier appears, leading to divergence.
        barrier_effect (float): Probability of divergence when encountering a barrier.
    
Returns:
        DendroPy Tree object.

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

### Simulate the tree with 10 species and a barrier effect at time 5.0
    simulated_tree = simulate_barrier_phylogeny(num_species=10, divergence_time=5.0, barrier_effect=0.7)

### Convert the tree to Newick format for visualization
    tree_newick = simulated_tree.as_string(schema="newick")
    print(tree_newick)

## Steps to Simulate Divergence Due to a Geographic Barrier with msprime (preferred)

1. Define a Population Split Model:

Before the barrier, populations exchange genes (migration).
When the barrier appears, migration drops to zero.
Over time, populations accumulate differences due to genetic drift.

2. Use msprime to Simulate the Process

### msprime Simulation of Divergence Due to a Geographic Barrier

    import msprime
    import demesdraw
    import matplotlib.pyplot as plt

### Define simulation parameters

    N_ancestral = 10000  # Ancestral population size
    N_population = 5000   # Population sizes after split
    T_divergence = 1000   # Time when the barrier appears
    migration_rate = 0.01 # Migration rate before barrier

### Define demographic events

    demography = msprime.Demography()
    demography.add_population(name="ancestral", initial_size=N_ancestral)
    demography.add_population(name="pop1", initial_size=N_population)
    demography.add_population(name="pop2", initial_size=N_population)

### Define the split event when the barrier appears
    demography.add_population_split(time=T_divergence, derived=["pop1", "pop2"], ancestral="ancestral")

### Migration before the barrier
    demography.set_migration_rate(source="pop1", dest="pop2", rate=migration_rate)
    demography.set_migration_rate(source="pop2", dest="pop1", rate=migration_rate)

### Stop migration at time of barrier
    demography.add_migration_rate_change(time=T_divergence, rate=0)

### Simulate tree sequence
    ts = msprime.sim_ancestry(samples={"pop1": 10, "pop2": 10}, demography=demography, sequence_length=1e6, recombination_rate=1e-8)

### Draw the demographic model
    fig, ax = plt.subplots(figsize=(6,4))
    demesdraw.tubes(demography.to_demes(), ax=ax)
    plt.show()

### Explanation of the Model

Ancestral Population ("ancestral"): Starts with N_ancestral = 10,000.
Barrier Forms at T_divergence = 1000:
pop1 and pop2 split from a shared ancestral population.
Migration is initially 0.01, then drops to 0 when the barrier appears.
Two Descendant Populations Evolve Independently:
Without gene flow, they accumulate genetic differences over time.

### Visualizing Divergence Times

To extract divergence times from the simulated coalescent tree:

    import numpy as np

### Get divergence times from the first tree in the sequence
    tree = ts.first()
    divergence_times = [tree.time(root) for root in tree.roots]

#### Print average divergence time
    print(f"Average divergence time: {np.mean(divergence_times)} generations")

### Next Steps

Simulate genetic diversity with mutations:

    ts = msprime.sim_mutations(ts, rate=1e-8, model="jc69")

Convert to Newick format for visualization in R (ape):

    tree = ts.first().as_newick()
    print(tree)  # Save this to a .nwk file for R visualization


### This is great but what if we want to simulate over millions of years? msprime is optimized for simulating genetic divergence over millions of years while efficiently handling large-scale coalescent events

You can extend the time scale by adjusting parameters like:

Divergence time (T_divergence) to reflect millions of years.
Mutation rate (μ) based on realistic substitution rates.
Recombination rate (r) to model genomic evolution.
Effective population sizes (N) reflecting long-term trends.

### Simulating Species Divergence Over 5 Millions of Years due to a geograohic barrier in msprime

    import msprime
    import demesdraw
    import matplotlib.pyplot as plt

### Define simulation parameters

    N_ancestral = 50_000  # Effective population size before divergence
    N_population = 20_000  # Population sizes after split
    T_divergence = 5_000_000  # Divergence time (5 million years)
    generation_time = 1  # Assume 1 year per generation (adjust if needed)
    migration_rate = 0.005  # Migration rate before barrier

### Convert divergence time into generations
    T_generations = T_divergence / generation_time

### Define demographic model
    demography = msprime.Demography()
    demography.add_population(name="ancestral", initial_size=N_ancestral)
    demography.add_population(name="pop1", initial_size=N_population)
    demography.add_population(name="pop2", initial_size=N_population)

### Define the split event at 5 million years
    demography.add_population_split(time=T_generations, derived=["pop1", "pop2"], ancestral="ancestral")

### Migration before the barrier
    demography.set_migration_rate(source="pop1", dest="pop2", rate=migration_rate)
    demography.set_migration_rate(source="pop2", dest="pop1", rate=migration_rate)

### Stop migration at the barrier formation
    demography.add_migration_rate_change(time=T_generations, rate=0)

### Simulate the tree sequence
    ts = msprime.sim_ancestry(samples={"pop1": 20, "pop2": 20}, demography=demography, sequence_length=1e7, recombination_rate=1e-9)

### Visualize the demographic history
    fig, ax = plt.subplots(figsize=(6, 4))
    demesdraw.tubes(demography.to_demes(), ax=ax)
    plt.show()

### Extract and display divergence times
    tree = ts.first()
    divergence_times = [tree.time(root) for root in tree.roots]
    print(f"Average divergence time: {sum(divergence_times) / len(divergence_times)} generations")

### Adjusting for Realism

Use Realistic Mutation Rates

Mammalian genome: ~1 × 10⁻⁸ mutations per site per generation.

Angiosperms: ~1 × 10⁻⁸ — ~1 × 10⁻9 (Check the literature and Brown et al., 2016. AJB https://doi.org/10.3732/ajb.1500117)

Example:
    ts = msprime.sim_mutations(ts, rate=1e-8, model="jc69")

This will introduce substitutions over evolutionary time.

Use Large Chromosome Lengths for Genome-Wide Simulations

Example: Simulate 100 Mbp genome:

    sequence_length = 1e8  # 100 million bases

#### Change Generation Time

If species have longer generation times (e.g., 10 years per generation), adjust:

    T_generations = T_divergence / 10

Example:

Bacteria: ~1-hour generation time → 10⁷ generations per million years.
Mammals: ~10-year generation time → 10⁵ generations per million years.

#### Extracting Divergence Events

You can extract divergence times and convert them into a phylogenetic tree:

    tree_newick = ts.first().as_newick()
    print(tree_newick)  # Save this to "tree.nwk" for visualization in R


## Simulating 10 million years in msprime is computationally feasible, but it depends on three key factors:

Population Size (N): Larger populations take longer to simulate.
Sample Size (n): The number of individuals sampled affects memory and runtime.
Sequence Length (L): Whole-genome simulations increase computational cost.

msprime is optimized for large-scale simulations, and since it's based on the coalescent model, it is far more efficient than forward-time simulators (like SLiM or fwdpy11). However, some adjustments may be needed.

Factors impacting runtime and memory usage:

Population Size (N)	10,000 (fast)	1,000,000 (slow)

Number of Samples (n)	10 - 100 samples	1,000+ samples

Divergence Time (T)	1M - 10M years (manageable)	100M+ years (longer)

Sequence Length (L)	1 Mbp (fast)	Whole genome (slow)

Recombination Rate (r)	1e-9 (efficient)	1e-8 (higher cost)

A 10-million-year simulation is feasible if:

N < 100,000
n < 500
L < 1e9 (1 billion bases)
Moderate mutation rate (1e-8)

## Example script for simulating 10 Million Years Efficiently

    import msprime
    import demesdraw
    import matplotlib.pyplot as plt

### Parameters

    N_ancestral = 50_000  # Effective ancestral population size
    N_population = 20_000  # Population size after split
    T_divergence = 10_000_000  # 10 million years
    generation_time = 10  # 10 years per generation
    migration_rate = 0.001  # Pre-barrier migration rate

### Convert time from years to generations
    T_generations = T_divergence / generation_time

### Define demographic history
    demography = msprime.Demography()
    demography.add_population(name="ancestral", initial_size=N_ancestral)
    demography.add_population(name="pop1", initial_size=N_population)
    demography.add_population(name="pop2", initial_size=N_population)

### Add divergence event
    demography.add_population_split(time=T_generations, derived=["pop1", "pop2"], ancestral="ancestral")

### Pre-barrier migration, then stop at divergence
    demography.set_migration_rate(source="pop1", dest="pop2", rate=migration_rate)
    demography.set_migration_rate(source="pop2", dest="pop1", rate=migration_rate)
    demography.add_migration_rate_change(time=T_generations, rate=0)

### Simulate tree sequence (with realistic sequence length)
    ts = msprime.sim_ancestry(
    samples={"pop1": 20, "pop2": 20}, 
    demography=demography, 
    sequence_length=1e7,  # 10 million bases
    recombination_rate=1e-9  # Realistic recombination rate
)

### Plot demographic model
    fig, ax = plt.subplots(figsize=(6, 4))
    demesdraw.tubes(demography.to_demes(), ax=ax)
    plt.show()

### Print divergence times
    tree = ts.first()
    divergence_times = [tree.time(root) for root in tree.roots]
    print(f"Average divergence time: {sum(divergence_times) / len(divergence_times)} generations")

##Optimizations for Faster Computation

1. Use Coalescent Advantage (Lower Sample Size)

msprime uses coalescent theory, so it doesn't simulate all individuals, only the ancestors of your sampled genomes.
Reducing samples={"pop1": X, "pop2": Y} to X=10-50 greatly reduces runtime.

2. Limit Genome Length

If you don't need a whole-genome simulation, set sequence_length = 1e6 (1 Mbp) instead of 1e9 (1 Gbp).

3. Use Lower Recombination Rate

Instead of the human recombination rate (1e-8), use 1e-9 to reduce computational overhead.

4. Use Fewer Population Events

Too many migration rate changes and population splits slow down the simulation.
