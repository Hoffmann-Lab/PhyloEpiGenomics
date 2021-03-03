#installation
install.packages("devtools", repos="http://cloud.r-project.org", clean=T)
devtools::install_github("hoffmann-lab/PhyloEpiGenomics", upgrade="never", force=T, clean=T)

data(PhyloEpiGenomics_example_data)

head(nucl_aln)
head(meth_fraction_aln)

nucl_states_aln=nucl_states_aln=sapply(nucl_aln,function(x) as.numeric(factor(x,levels=c("A","C","G","T"))))
head(nucl_states_aln)

discretization=list(c(-0.01,0.2),c(0.2,0.4),c(0.4,0.6),c(0.6,0.8),c(0.8,1))
meth_states_aln=discretize(meth_fraction_aln,discretization)
head(meth_states_aln)

#model parameters are estimated from data
my_JC69_model=make_evolutionary_model(nucl_states_aln,model="JC69")
my_JC69_model#evolutionary models consist of a transition rate matrix Q and an equilibrium frequency vector pi

# model parameters can also specified manually
my_HKY85_model=make_evolutionary_model(model="HKY85",pi=rep(0.25,4),kappa = 2)

# make a 5 states methylation data model according to the discretization used above
my_noJump_model=make_evolutionary_model(meth_states_aln,model="noJump",nstates=length(discretization))
my_noJump_model

#Tree reconstruction functions expect a list of tree topologies to examine/compare
unrooted_tree_topologies=all_unrooted_tree_topologies(colnames(nucl_states_aln))
plot(unrooted_tree_topologies[[3]],type="unrooted",lab4ut="axial")
ml_best_nucl_tree=find_optimal_tree(nucl_states_aln,trees=unrooted_tree_topologies,Q=my_HKY85_model$Q,pi=my_HKY85_model$pi)
plot(ml_best_nucl_tree,type="unrooted",lab4ut="axial") 

ml_best_meth_tree=find_optimal_tree(meth_states_aln,trees=unrooted_tree_topologies,Q=my_noJump_model$Q,pi=my_noJump_model$pi)
plot(ml_best_meth_tree,type="unrooted",lab4ut="axial") 