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
#examined/compared/optimized; in this case unrooted topologies are used
unrooted_tree_topologies = all_unrooted_tree_topologies(colnames(nucl_states_aln))
#one of the wrong topologies
plot(unrooted_tree_topologies[[2]],
     type = "unrooted",
     lab4ut = "axial")
#correct topology
plot(unrooted_tree_topologies[[3]],
     type = "unrooted",
     lab4ut = "axial",
     rotate.tree = 270)

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

#15 tree topologies to examine
cluster2 = makeCluster(max(length(rooted_tree_topologies), detectCores()))

ml_best_nucl_tree_clock = find_optimal_tree(
  cluster = cluster2,
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
#the proportions of the nucleotide tree

#the scenario that let all branches vary freely has four degrees of freedom more
#than the scenario that maintained the proportions of the nucl tree
pchisq(2 * (ml_best_meth_tree$lnL - ml_meth_tree_based_on_nucl$lnL),
       4,
       lower.tail = F)
#-> the scenario is also significantly more likely than the scenario that
#maintained the proportions of the nucleotide tree

#DISTANCE BASED

#get the distance matrix using maximum likelihood
distance_matrix_nucl = maximize_distance_log_likelihoods(nucl_states_aln, Q =
                                                           my_HKY85_model$Q, pi = my_HKY85_model$pi)
distance_matrix_nucl

#reconstruct the tree using neighbour joining from the ape package
require(ape)
nj_best_nucl_tree = nj(distance_matrix_nucl)

#reconstruct the tree using fitch-margoliash
fm_best_nucl_tree = best_fitch_margoliash(distance_matrix_nucl, unrooted_tree_topologies)

kronoviz(
  list(ml_best_nucl_tree,
       nj_best_nucl_tree,
       fm_best_nucl_tree),
  type = "unrooted",
  lab4ut = "axial",
  rotate.tree = 270
)

#PARSIMONY

#parsimony requires rooted trees as input but all rooted trees that represent
#the same unrooted tree result in the same costs, i.e. the method cannot
#determine a good outgroup
pars_best_meth_trees = all_parsimony(meth_states_aln, rooted_trees = rooted_tree_topologies)
sapply(pars_best_meth_trees, function(tree)
  tree$cost_sum)
#these solution are therefore equivalent
layout(matrix(c(1:5, 0), nrow = 2))
sapply(pars_best_meth_trees[11:15], plot)

#Parsimony does not assign branch lengths, you can use, e.g., maximum likelihood
#for that purpose
pars_best_meth_tree = maximize_tree_log_likelihood(
  meth_states_aln,
  tree = pars_best_meth_trees[[15]],
  Q = my_noJump_model$Q,
  pi = my_noJump_model$pi
)

layout(matrix(1))
plot(pars_best_meth_tree)

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

ml_best_nucl_tree_clock$node.label=setNames(5:7,5:7)
plot(ml_best_nucl_tree_clock,show.node.label = T)

#convert to ACGT
sim_nucl_aln=apply(sim_nucl$states_alignment,c(1,2),function(x) c("A","C","G","T")[x])
head(sim_nucl_aln)

#simulate the evolution of 3000 methylation fractions, if the discretization is
#passed instead of states also fractions are returned - assuming a uniform
#within state distribution; in the example below there is also some noise added
#to simulate, e.g., sampling errors due to low coverage

#get a tree for the simulation from real data

sim_meth_noise = simulate_evolution(
  nstates = 3000,
  tree = pars_best_meth_tree$tree,
  Q = my_noJump_model$Q,
  pi = my_noJump_model$pi,
  discretization = discretization,
  noise_sd = 0.1
)

head(sim_meth_noise$states_alignment)
head(sim_meth_noise$frequency_alignment)
head(sim_meth_noise$frequency_alignment_noise)

#disturbs the alignment to simulate, e.g., occasional alignment problems

sim_meth_noise_disturbed = exchange_states(
  sim_meth_noise$frequency_alignment_noise[,1:4],
  exchange_dist = 2,
  exchange_prob = 0.2
)
head(sim_meth_noise_disturbed)

#GRAPHICS FOR README

dir.create("PhyloEpiGenomics/read.me_plots")
jpeg("PhyloEpiGenomics/readme_plots/wrong_unrooted_topology.jpg")
plot(unrooted_tree_topologies[[2]],
     type = "unrooted",
     lab4ut = "axial")
dev.off()
jpeg("PhyloEpiGenomics/readme_plots/correct_unrooted_topology.jpg")
plot(unrooted_tree_topologies[[3]],
     type = "unrooted",
     lab4ut = "axial",
     rotate.tree = 270)
dev.off()
jpeg("PhyloEpiGenomics/readme_plots/ml_best_nucl_tree.jpg")
plot(
  ml_best_nucl_tree,
  type = "unrooted",
  lab4ut = "axial",
  rotate.tree = 270
)
dev.off()
jpeg("PhyloEpiGenomics/readme_plots/ml_best_meth_tree.jpg")
plot(
  ml_best_meth_tree,
  type = "unrooted",
  lab4ut = "axial",
  rotate.tree = 270
)
dev.off()
jpeg("PhyloEpiGenomics/readme_plots/wrong_rooted_topology.jpg")
plot(rooted_tree_topologies[[15]], lab4ut = "axial")#correct topology
dev.off()
jpeg("PhyloEpiGenomics/readme_plots/correct_rooted_topology.jpg")
plot(rooted_tree_topologies[[11]], lab4ut = "axial")#one of the wrong topologies
dev.off()
jpeg("PhyloEpiGenomics/readme_plots/ml_best_nucl_tree_clock.jpg")
plot(ml_best_nucl_tree_clock)
dev.off()
jpeg("PhyloEpiGenomics/readme_plots/ml_meth_tree_based_on_nucl.jpg")
kronoviz(
  list(ml_best_nucl_tree, ml_meth_tree_based_on_nucl$tree),
  type = "unrooted",
  lab4ut = "axial",
  rotate.tree = 270,
)
add.scale.bar()
dev.off()
jpeg("PhyloEpiGenomics/readme_plots/ml_meth_tree_based_on_nucl_2.jpg")
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
dev.off()
jpeg("PhyloEpiGenomics/readme_plots/distance_based.jpg")
kronoviz(
  list(ml_best_nucl_tree,
       nj_best_nucl_tree,
       fm_best_nucl_tree),
  type = "unrooted",
  lab4ut = "axial",
  rotate.tree = 270,
  cex=1.2
)
dev.off()
jpeg("PhyloEpiGenomics/readme_plots/parsimony.jpg")
layout(matrix(c(1:5, 0), nrow = 2))
sapply(pars_best_meth_trees[11:15], function(x) plot(x,cex=1.2))
dev.off()
jpeg("PhyloEpiGenomics/readme_plots/parsimony_2.jpg")
layout(matrix(1))
plot(pars_best_meth_tree$tree)
dev.off()
jpeg("PhyloEpiGenomics/readme_plots/inner_nodes_labeled.jpg")
plot(ml_best_nucl_tree_clock,show.node.label = T)
dev.off()