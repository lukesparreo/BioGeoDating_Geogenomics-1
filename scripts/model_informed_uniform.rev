# RevBayes code for geology informed (correct) model with uniform distribution

range_fn = "simulated_range.nex"
mol_fn = "modified_sequences_filled.nex"
tree_fn = "collapsed_newick.tre"
out_fn = "output_informed_uniform_migration2" #Modify each run if needed
geo_fn = "simulated"
times_fn = geo_fn + ".times.informed.txt"
dist_fn = geo_fn + ".distances.txt"

# Read in molecular alignment
dat_mol = readDiscreteCharacterData(mol_fn)

# Read in species ranges
dat_range_01 = readDiscreteCharacterData(range_fn)

# Compute the number of ranges when ranges may only be one or two areas in size
n_areas <- dat_range_01.nchar()
max_areas <- 3
n_states <- 0
for (k in 0:max_areas) n_states += choose(n_areas, k)
# Format
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

# Get the converted state descriptions
state_desc = dat_range_n.getStateDescriptions()

# Write the state descriptions to file
state_desc_str = "state,range\n"
for (i in 1:state_desc.size())
{
    state_desc_str += (i-1) + "," + state_desc[i] + "\n"
}
write(state_desc_str, file=out_fn+".state_labels.txt")

# TREE MODEL
# Define the root age

root_age ~ dnUniform(0, 10.5)

moves = VectorMoves()
moves.append( mvScale(root_age, weight=5) )

# Assign the proportion of sampled taxa (changed from non-uniform sampling scheme in Landis to complete sampling here
rho <- 3/3

# Assign birth/death priors
birth ~ dnExp(10)
moves.append( mvScale(birth, weight=2) )
death ~ dnExp(10)
moves.append( mvScale(death, weight=2) )

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

# MOLECULAR MODEL

# Base rate for molcular clock
rate_mol ~ dnLoguniform(1E-6, 1E0)
rate_mol.setValue(1E-2)
moves.append( mvScale(rate_mol, lambda=0.2, weight=4) )
moves.append( mvScale(rate_mol, lambda=1.0, weight=2) )

# Assign log-normal relaxed clock rate
branch_sd <- 1.0
branch_mean <- 0.0 - 0.5 * branch_sd^2
for (i in 1:n_branches) {
    branch_rate_multiplier[i] ~ dnLognormal(mean=branch_mean, sd=branch_sd)
    branch_rates[i] := rate_mol * branch_rate_multiplier[i]
    moves.append( mvScale(branch_rate_multiplier[i]) )
}
moves.append( mvVectorScale(branch_rate_multiplier, weight=3) )

# Create HKY rate matrix like we did for Newick
kappa ~ dnGamma(2,2)
moves.append( mvScale(kappa) )

bf ~ dnDirichlet([1,1,1,1])
moves.append( mvSimplexElementScale(bf, alpha=10, weight=2) )

Q_mol := fnHKY(kappa, bf)

alpha ~ dnUniform(0,50)
moves.append( mvScale(alpha) )

site_rates := fnDiscretizeGamma(alpha, alpha, 4)

m_mol ~ dnPhyloCTMC(Q=Q_mol,
                    tree=tree,
                    branchRates=branch_rates,
                    siteRates=site_rates,
                    type="DNA",
                    nSites=dat_mol.nchar())

m_mol.clamp(dat_mol)

# BIOGEOGRAPHIC MODEL

rate_bg ~ dnLoguniform(1E-4,1E2)
rate_bg.setValue(1E-2)
moves.append( mvScale(rate_bg, lambda=0.2, weight=4) )
moves.append( mvScale(rate_bg, lambda=1.0, weight=2) )

# Fix dispersal rate
dispersal_rate <- 0.2
distance_scale ~ dnUnif(0,20)
distance_scale.setValue(0.001)
moves.append( mvScale(distance_scale, weight=3) )

# Then, the dispersal rate matrix
for (i in 1:n_epochs) {
  for (j in 1:n_areas) {
    for (k in 1:n_areas) {
     dr[i][j][k] <- 0.0
     if (connectivity[i][j][k] > 0) {
       dr[i][j][k] := dispersal_rate * exp(-distance_scale * distances[j][k])
     }
    }
  }
}
            
# Extirpation rate
log_sd <- 0.5
log_mean <- ln(1) - 0.5*log_sd^2
extirpation_rate ~ dnLognormal(mean=log_mean, sd=log_sd)
moves.append( mvScale(extirpation_rate, weight=2) )

for (i in 1:n_epochs) {
  for (j in 1:n_areas) {
    for (k in 1:n_areas) {
      er[i][j][k] <- 0.0
    }
    er[i][j][j] := extirpation_rate
  }
}

# Build DEC rate matrices
for (i in 1:n_epochs) {
  Q_DEC[i] := fnDECRateMatrix(dispersalRates=dr[i],
                          extirpationRates=er[i],
                          maxRangeSize=max_areas)
}
            
# Build the epoch times
for (i in 1:n_epochs) {
  time_max[i] <- time_bounds[i][1]
  time_min[i] <- time_bounds[i][2]
  if (i != n_epochs) {
    epoch_times[i] ~ dnUniform(time_min[i], time_max[i])
    epoch_width = time_bounds[i][1] - time_bounds[i][2]
    moves.append( mvSlide(epoch_times[i], delta=epoch_width/2) )
  } else {
    epoch_times[i] <- 0.0
  }
}
                           
# Combine the epoch rate matrices and times
Q_DEC_epoch := fnEpoch(Q=Q_DEC, times=epoch_times, rates=rep(1, n_epochs))

# Build cladogenetic transition probabilities
clado_event_types <- [ "s", "a" ]
p_sympatry ~ dnUniform(0,1)
p_allopatry := abs(1.0 - p_sympatry)
moves.append( mvSlide(p_sympatry, delta=0.1, weight=2) )
clado_event_probs := simplex(p_sympatry, p_allopatry)
# Note: P_DEC is defined, but you can't view it from print(). Instead check print(type(P_DEC)) and str(P_DEC)
P_DEC := fnDECCladoProbs(eventProbs=clado_event_probs,
                         eventTypes=clado_event_types,
                         numCharacters=n_areas,
                         maxRangeSize=max_areas)
                       
# Root frequencies
rf_DEC_raw            <- rep(0, n_states)
rf_DEC_raw[n_areas+1] <- 1  # "Mainland" (original river) is the only possible starting state
rf_DEC                <- simplex(rf_DEC_raw)
    
# The phylogenetic CTMC with cladogenetic events
m_bg ~ dnPhyloCTMCClado(tree=tree,
                           Q=Q_DEC_epoch,
                           cladoProbs=P_DEC,
                           branchRates=rate_bg,
                           rootFrequencies=rf_DEC,
                           type="NaturalNumbers",
                           nSites=1)     

# Attach the range data
m_bg.clamp(dat_range_n)

# Monitor the age of the ingroup
ingroup_clade <- clade("n0",
                       "n1",
                       "n2")

# Set ingroup age
ingroup_age := tmrca(tree, ingroup_clade)

for (i in 1:n_epochs) {
    ingroup_older_island[i] := ifelse(ingroup_age > epoch_times[i], 1, 0)
}

monitors = VectorMonitors()
monitors.append( mnScreen(printgen=100, ingroup_age) )
monitors.append( mnModel(file=out_fn+".model.log", printgen=100) )
monitors.append( mnFile(tree, filename=out_fn+".tre", printgen=100) )
monitors.append( mnJointConditionalAncestralState(tree=tree,
                                                  ctmc=m_bg,
                                                  type="NaturalNumbers",
                                                  withTips=true,
                                                  withStartStates=true,
                                                  filename=out_fn+".states.log",
                                                  printgen=100) )
monitors.append( mnStochasticCharacterMap(ctmc=m_bg,
                                          filename=out_fn+".stoch.log",
                                          printgen=100) )

# Analysis generations
n_gen = 3000000
    
# Create model
mymodel = model(m_bg, ingroup_older_island)

# Run
mymcmc = mcmc(mymodel, moves, monitors)
mymcmc.run(n_gen)

# Results in output

###

# Summarizing output
# Go to folder of output once run is completed
out_str = "output_informed_uniform_migration2" #May need to modify depending on output filename
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

