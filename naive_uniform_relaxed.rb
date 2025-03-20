### REVBAYES CODE FOR GEO-INFORMED RELAXED MODEL ###

range_fn = "simulated_range.nex"
mol_fn = "modified_sequences.nex"
tree_fn = "collapsed_newick.tre"
out_fn = "output_naive_uniform_relaxed_1/simulationoutput" #MODIFY EACH RUN!
geo_fn = "/Users/lukesparreo/simulated_data/simulated"
times_fn = geo_fn + ".times.naiverelaxed.txt" #MODIFY EACH RUN!
dist_fn = geo_fn + ".distances.txt"

# Analysis helper variables
n_gen = 1000000

# Read in molecular alignment
dat_mol = readDiscreteCharacterData(mol_fn)

# Read in species ranges
dat_range_01 = readDiscreteCharacterData(range_fn)

# Compute the number of ranges when ranges may only be one or two areas in size
n_areas <- dat_range_01.nchar()
max_areas <- 3
n_states <- 0
for (k in 0:max_areas) n_states += choose(n_areas, k)
#format
dat_range_n = formatDiscreteCharacterData(dat_range_01, "DEC", n_states)

# Record list of ranges to a file
state_desc = dat_range_n.getStateDescriptions()
state_desc_str = "state,range\n"
for (i in 1:state_desc.size())
{
    state_desc_str += (i-1) + "," + state_desc[i] + "\n"
}
write(state_desc_str, file=out_fn+".state_labels.txt")

# Read the minimum and maximum ages of the barrier events
time_bounds <- readDataDelimitedFile(file=times_fn, delimiter=" ")
n_epochs <- time_bounds.size()

# Read in connectivity matrices
for (i in 1:n_epochs) {
  epoch_fn[i] = geo_fn + ".connectivity." + i + ".txt"
  connectivity[i] <- readDataDelimitedFile(file=epoch_fn[i], delimiter=" ")
}

# Read the distances
distances <- readDataDelimitedFile(file=dist_fn, delimiter=" ")

# Read the tree file
tree_init = readTrees(tree_fn)[1]

# And record information about tree file
taxa = tree_init.taxa()
n_taxa = taxa.size()
n_branches = 2 * n_taxa - 2

# TREE MODEL
root_age ~ dnUniform(3, 4)
moves = VectorMoves()
moves.append( mvScale(root_age, weight=2.5) )

# Assign birth/death priors
birth ~ dnExp(10)
moves.append( mvScale(birth, weight=1.5) )
death ~ dnExp(10)
moves.append( mvScale(death, weight=1.5) )

# Initiate tree
tree ~ dnBDP(lambda=birth, mu=death, rho=rho, rootAge=root_age, taxa=taxa)

# Topology and branch lengths
moves.append( mvNNI(tree, weight=n_branches/2) )
moves.append( mvFNPR(tree, weight=n_branches/8) )
moves.append( mvNodeTimeSlideUniform(tree, weight=n_branches/2) )
moves.append( mvSubtreeScale(tree, weight=n_branches/8) )
moves.append( mvTreeScale(tree, root_age, weight=n_branches/8) )

# Provide starting tree for biogeographic model 
tree.setValue(tree_init)
root_age.setValue(tree_init.rootAge())

# Base rate for molecular clock
rate_mol ~ dnLoguniform(1E-6, 1E0)
rate_mol.setValue(1E-2)
moves.append( mvScale(rate_mol, lambda=0.2, weight=2.5) )
moves.append( mvScale(rate_mol, lambda=1.0, weight=1.5) )

# Assign log-normal relaxed clock rate
branch_sd <- 1.0
branch_mean <- 0.0 - 0.5 * branch_sd^2
for (i in 1:n_branches) {
    branch_rate_multiplier[i] ~ dnLognormal(mean=branch_mean, sd=branch_sd)
    branch_rates[i] := rate_mol * branch_rate_multiplier[i]
    moves.append( mvScale(branch_rate_multiplier[i], weight=2) )
}
moves.append( mvVectorScale(branch_rate_multiplier, weight=2.5) )

# Create HKY rate matrix
kappa ~ dnGamma(2,2)
moves.append( mvScale(kappa, weight=2) )

bf ~ dnDirichlet([1,1,1,1])
moves.append( mvSimplexElementScale(bf, alpha=10, weight=1.5) )

alpha ~ dnUniform(0,50)
moves.append( mvScale(alpha, weight=2) )

site_rates := fnDiscretizeGamma(alpha, alpha, 4)

m_mol ~ dnPhyloCTMC(Q=Q_mol,
                    tree=tree,
                    branchRates=branch_rates,
                    siteRates=site_rates,
                    type="DNA",
                    nSites=dat_mol.nchar())

m_mol.clamp(dat_mol)

# Adjust epoch time moves
for (i in 1:n_epochs) {
  time_max[i] <- time_bounds[i][1]
  time_min[i] <- time_bounds[i][2]
  if (i != n_epochs) {
    epoch_times[i] ~ dnUniform(time_min[i], time_max[i])
    epoch_width = time_bounds[i][1] - time_bounds[i][2]
    moves.append( mvSlide(epoch_times[i], delta=epoch_width/4) )
  } else {
    epoch_times[i] <- 0.0
  }
}

# Adjust dispersal rate scaling move
moves.append( mvScale(distance_scale, weight=2) )

# Adjust cladogenetic event probability move
moves.append( mvSlide(p_sympatry, delta=0.05, weight=1.5) )

# Run MCMC
mymodel = model(m_bg, ingroup_older_island)
mymcmc = mcmc(mymodel, moves, monitors)
mymcmc.run(n_gen)

# Results in /simulated_data/output
###

##Summarizing output

out_str = "output_naive_uniform_relaxed_1/simulationoutput" #MODIFY EACH RUN!
out_state_fn = out_str + ".states.log"
out_tree_fn = out_str + ".tre"
out_mcc_fn = out_str + ".mcc.tre"

tree_trace = readTreeTrace(file=out_tree_fn, treetype="clock")
tree_trace.setBurnin(0.25)
n_burn = tree_trace.getBurnin()
mcc_tree = mccTree(tree_trace, file=out_mcc_fn)

state_trace = readAncestralStateTrace(file=out_state_fn)

tree_trace = readAncestralStateTreeTrace(file=out_tree_fn, treetype="clock")

anc_tree = ancestralStateTree(tree=mcc_tree,
                               ancestral_state_trace_vector=state_trace,
                               tree_trace=tree_trace,
                               include_start_states=true,
                               file=out_str+".ase.tre",
                               burnin=n_burn,
                               site=1)
