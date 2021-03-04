#INSTALLATION
install.packages("devtools", repos = "http://cloud.r-project.org", clean =
                   T)
devtools::install_github(
  "hoffmann-lab/PhyloEpiGenomics",
  upgrade = "never",
  force = T,
  clean = T
)

#DATA ORGANIZATION AND PREPROCESSING

#load library and example data
library(PhyloEpiGenomics)
library(ape)
data(PhyloEpiGenomics_example_data)

#original data
head(nucl_aln)
head(meth_fraction_aln)

#data converted to required numeric states
nucl_states_aln = nucl_states_aln = sapply(nucl_aln, function(x)
  as.numeric(factor(x, levels = c("A", "C", "G", "T"))))
head(nucl_states_aln)

discretization = list(c(-0.01, 0.2), c(0.2, 0.4), c(0.4, 0.6), c(0.6, 0.8), c(0.8, 1))
meth_states_aln = discretize(meth_fraction_aln, discretization)
head(meth_states_aln)

#MAXIMUM LIKELIHOOD

#model parameters are estimated from data
my_JC69_model = make_evolutionary_model(nucl_states_aln, model = "JC69")
#evolutionary models consist of a transition rate matrix Q and an equilibrium
#frequency vector pi
my_JC69_model

# model parameters can also specified manually
my_HKY85_model = make_evolutionary_model(model = "HKY85",
                                         pi = rep(0.25, 4),
                                         kappa = 2)

# make a 5 states methylation data model according to the discretization used
# above
my_noJump_model = make_evolutionary_model(meth_states_aln,
                                          model = "noJump",
                                          nstates = length(discretization))
my_noJump_model

#you can, of course, also define your own models for nucleotide, amino acid,
#codon, epigenomic or other data!


#tree reconstruction functions expect a list of tree topologies that are to be
#examined/compared
unrooted_tree_topologies = all_unrooted_tree_topologies(colnames(nucl_states_aln))
#one of the wrong topologies
plot(unrooted_tree_topologies[[2]],
     type = "unrooted",
     lab4ut = "axial")
#correct topology
plot(unrooted_tree_topologies[[3]],
     type = "unrooted",
     lab4ut = "axial")

#tree reconstruction via maximum likelihood, nucleotide example
ml_best_nucl_tree = find_optimal_tree(nucl_states_aln,
                                      trees = unrooted_tree_topologies,
                                      Q = my_HKY85_model$Q,
                                      pi = my_HKY85_model$pi)
plot(
  ml_best_nucl_tree,
  type = "unrooted",
  lab4ut = "axial",
  rotate.tree = 270
)

#tree reconstruction via maximum likelihood using multiple CPU cores,
#methylation example
library(parallel)
#using up to as many CPUs as examined trees reduces the runtime
cluster = makeCluster(max(3, detectCores()))
ml_best_meth_tree = find_optimal_tree(
  cluster = cluster,
  meth_states_aln,
  trees = unrooted_tree_topologies,
  Q = my_noJump_model$Q,
  pi = my_noJump_model$pi
)
plot(
  ml_best_meth_tree,
  type = "unrooted",
  lab4ut = "axial",
  rotate.tree = 270
)

#tree reconstruction via maximum likelihood using the molecular clock constraint

#using a clock implies rooted trees
rooted_tree_topologies = all_rooted_tree_topologies(colnames(nucl_states_aln))
plot(rooted_tree_topologies[[15]], lab4ut = "axial")#correct topology
plot(rooted_tree_topologies[[11]], lab4ut = "axial")#one of the wrong topologies

stopCluster(cluster)
cluster = makeCluster(max(15, detectCores()))#15 tree topologies to examine

ml_best_nucl_tree_clock = find_optimal_tree(
  cluster = cluster,
  nucl_states_aln,
  trees = rooted_tree_topologies,
  Q = my_HKY85_model$Q,
  pi = my_HKY85_model$pi,
  clock = T
)
plot(ml_best_nucl_tree_clock)

#use the same rate class "1" for all 5 edges of the unrooted tree
branches_to_evolutionary_rate_classes = rep(1, nrow(ml_best_nucl_tree$edge))
nrow(ml_best_nucl_tree$edge)
branches_to_evolutionary_rate_classes

#reconstructs a tree on methylation data expanding/contracting the form of
#nucleotide tree
ml_meth_tree_based_on_nucl = maximize_tree_log_likelihood_extended(
  meth_states_aln,
  tree = ml_best_nucl_tree,
  Q = my_noJump_model$Q,
  my_noJump_model$pi,
  branches_to_evolutionary_rate_classes = branches_to_evolutionary_rate_classes
)

#resulting meth tree (left bottom) has the same proportions as the nucl tree
#(right top)
kronoviz(
  list(ml_best_nucl_tree, ml_meth_tree_based_on_nucl$tree),
  type = "unrooted",
  lab4ut = "axial",
  rotate.tree = 270
)
add.scale.bar()

#now do the same as above but let the human branch vary freely:
#which edge represents human?
h = which(ml_best_nucl_tree$edge[, 2] == which(ml_best_nucl_tree$tip.label ==
                                                 "human"))
#use rate class for human that differs from the rest
branches_to_evolutionary_rate_classes[h] = 2
h
branches_to_evolutionary_rate_classes

ml_meth_tree_based_on_nucl_2 = maximize_tree_log_likelihood_extended(
  meth_states_aln,
  tree = ml_best_nucl_tree,
  Q = my_noJump_model$Q,
  my_noJump_model$pi,
  branches_to_evolutionary_rate_classes = branches_to_evolutionary_rate_classes
)

#while the other 4 branches maintain their proportions, the human branches gets
#longer in the second scenario (down)
kronoviz(
  list(
    ml_meth_tree_based_on_nucl$tree,
    ml_meth_tree_based_on_nucl_2$tree
  ),
  type = "unrooted",
  lab4ut = "axial",
  rotate.tree = 270
)
add.scale.bar()

#applying likelihood-ratio tests

#the scenario with the free human branch has one degree of freedom more than the
#scenario that maintained the proportions of the nucl tree
pchisq(
  2 * (
    ml_meth_tree_based_on_nucl_2$lnL - ml_meth_tree_based_on_nucl$lnL
  ),
  1,
  lower.tail = F
)
#-> the scenario s significantly more likely than the scenario that maintained
#the proportions of the nucl tree

#the scenario that let all branches vary freely has four degrees of freedom more
#than the scenario that maintained the proportions of the nucl tree
pchisq(2 * (ml_best_meth_tree$lnL - ml_meth_tree_based_on_nucl$lnL),
       4,
       lower.tail = F)
#-> the scenario is also significantly more likely than the scenario that
#maintained the proprtions of the nucl tree

#DISTANCE BASED

#get the distance matrix using maximum likelihood
distance_matrix_nucl = maximize_distance_log_likelihoods(nucl_states_aln, Q =
                                                           my_HKY85_model$Q, pi = my_HKY85_model$pi)

#reconstruct the tree using neighbour joining from the ape package
nj_best_nucl_tree = nj(distance_matrix_nucl)
plot(
  nj_best_nucl_tree,
  type = "unrooted",
  lab4ut = "axial",
  rotate.tree = 90
)

#reconstruct the tree using fitch-margoliash
fm_best_nucl_tree = best_fitch_margoliash(distance_matrix_nucl, unrooted_tree_topologies)

#top: ml, middle: nj, bottom: fm
kronoviz(
  list(ml_best_nucl_tree,
       nj_best_nucl_tree,
       fm_best_nucl_tree),
  type = "unrooted",
  lab4ut = "axial",
  rotate.tree = 270
)

#PARSIMONY

#parsimony requires rooted trees as input
pars_best_meth_tree = maximum_parsimony(meth_states_aln, rooted_trees = rooted_tree_topologies)
plot(pars_best_meth_tree)

#but all rooted trees that represent the same unrooted tree result in the same
#costs
pars_best_meth_trees = all_parsimony_rooted_trees(meth_states_aln, rooted_trees = rooted_tree_topologies)
sapply(pars_best_meth_trees, function(tree)
  tree$cost_sum)
#these solution are therefore equivalent
layout(matrix(c(1:5, 0), nrow = 2))
sapply(pars_best_meth_trees[11:15], plot)

#SIMULATION 

#for the simulation a rooted tree is needed
sim_tree_nucl=ml_best_nucl_tree_clock
#the example takes the tree calculated before but multiplies the time with 7 to
#result in more substitutions
sim_tree_nucl$edge.length=sim_tree_nucl$edge.length*7
#simulate the evolution of 5000 nucleotides
sim_nucl = simulate_evolution(
  nstates = 5000,
  tree = sim_tree_nucl,
  Q = my_HKY85_model$Q,
  pi = my_HKY85_model$pi
)
#the simulation contains also the states at the ancestral nodes
head(sim_nucl$states_alignment)
layout(matrix(1))
ml_best_nucl_tree_clock$node.label=setNames(5:7,5:7)
plot(ml_best_nucl_tree_clock,show.node.label = T)
#convert to ACGT
sim_nucl_aln=apply(sim_nucl$states_alignment,c(1,2),function(x) c("A","C","G","T")[x])
head(sim_nucl_aln)

#simulate the evolution of 3000 methylation fractions, if the discretization is
#passed instead of states fractions are returned directly - assuming a uniform
#within state distribution
sim_meth = simulate_evolution(
  nstates = 3000,
  tree = ml_best_meth_tree,
  Q = my_noJump_model$Q,
  pi = my_noJump_model$pi,
  discretization = discretization
)
head(sim_meth$states_alignment)
head(sim_meth$frequency_alignment)



