#' @author Arne Sahm \email{arne.sahm@@leibniz-fli.de}
#' @import expm
#' @import parallel
#' @import ape
#' @import data.table
NULL

#' Tree reconstruction via maximum parsimony
#'
#' `maximum_parsimony` determines the single, best fitting among the given tree topologies by maximum parsimony for given data. The respective determined costs wil be attached to the returning phylo object.\cr\cr
#' `all_parsimony` performs parsimony for all tree topologies given. The respective determined costs wil be attached to the phylo objects.\cr\cr
#' `small_parsimony performs` parsimony for one tree topology. The respective determined costs wil be attached to the phylo object.
#' @param states_table Matrix or data frame containing numericals representing an alignment of interval or ordinal scaled states with sites as rows and species/strains as columns.
#' @param rooted_trees List of rooted tree topologies of class phylo (library ape) to be analyzed.
#' @param rooted_tree Rooted tree topology of class phylo (library ape) to be analyzed.
#' @details 
#' These functions work on any interval scaled data such as methylation fractions or ordinal scaled data such as descritized methylation data. They do not work, however, on nominal scaled data, such as nucleotides.
#' @examples
#' test
#' @export
maximum_parsimony=function(states_table,rooted_trees){
  rooted_trees=all_parsimony(states_table,rooted_trees)
  rooted_trees[[which.min(sapply(simplify=F,rooted_trees,function(tree) tree$cost_sum))]]
}

#' @rdname  maximum_parsimony
#' @export
all_parsimony=function(states_table,rooted_trees) sapply(simplify=F,rooted_trees,function(tree) small_parsimony(states_table,tree))

#' @rdname  maximum_parsimony
#' @export
small_parsimony=function(states_table,rooted_tree){
  tree.root=rooted_tree$edge[1,1]
  interval_table=apply(states_table,c(1,2),function(x) list(x,x))
  colnames(interval_table)=1:ncol(interval_table)
  rooted_tree$interval_table=small_parsimony_bottom_up(interval_table,rooted_tree,tree.root)
  rooted_tree$states_table=small_parsimony_top_down(rooted_tree$interval_table,rooted_tree,tree.root)
  rooted_tree$cost_table=small_parsimony_costs(rooted_tree$states_table,rooted_tree,tree.root)
  rooted_tree$cost_sum=sum(rooted_tree$cost_table)
  return(rooted_tree)
}

small_parsimony_costs=function(states_table,rooted_tree,node,parent_node=NULL){
  if (is.null(parent_node)) current_node_costs=rep(0,nrow(states_table))
  else current_node_costs=apply(states_table[,as.character(c(node,parent_node))],1,function(states) abs(states[1]-states[2]))
  child_node_edge_indices=which(rooted_tree$edge[,1]==node)
  child_nodes=sapply(child_node_edge_indices,function(child_node_edge_index) rooted_tree$edge[child_node_edge_index,2])
  new_costs_table=do.call(cbind,sapply(simplify = F,child_nodes,function(child_node) small_parsimony_costs(states_table,rooted_tree,child_node,node)))
  new_costs_table=cbind(current_node_costs,new_costs_table)
  colnames(new_costs_table)[1]=node
  return(new_costs_table)
}

small_parsimony_top_down=function(interval_table,rooted_tree,node,parent_node_states=NULL){
  if (is.null(parent_node_states)) current_node_states=sapply(interval_table[,as.character(node)],function(interval) (interval[[1]]+interval[[2]])/2)
  else current_node_states=sapply(1:nrow(interval_table),function(i) {
    current_interval=interval_table[i,as.character(node)][[1]]
    if (parent_node_states[i]>=current_interval[[1]]){
      if (parent_node_states[i]<=current_interval[[2]]) return(parent_node_states[i])
      else return(current_interval[[2]])
    } else return(current_interval[[1]])
  })
  child_node_edge_indices=which(rooted_tree$edge[,1]==node)
  child_nodes=sapply(child_node_edge_indices,function(child_node_edge_index) rooted_tree$edge[child_node_edge_index,2])
  new_states_table=do.call(cbind,sapply(simplify = F,child_nodes,function(child_node) small_parsimony_top_down(interval_table,rooted_tree,child_node,current_node_states)))
  new_states_table=cbind(current_node_states,new_states_table)
  colnames(new_states_table)[1]=node
  return(new_states_table)
}

small_parsimony_bottom_up=function(interval_table,rooted_tree,node){
  if (node<=length(rooted_tree$tip.label)) return(matrix(interval_table[,as.character(node)],dimnames=list(NULL,node)))
  else{
    child_node_edge_indices=which(rooted_tree$edge[,1]==node)
    child_nodes=sapply(child_node_edge_indices,function(child_node_edge_index) rooted_tree$edge[child_node_edge_index,2])
    new_interval_table=do.call(cbind,sapply(simplify = F,child_nodes,function(child_node) small_parsimony_bottom_up(interval_table,rooted_tree,child_node)))
    current_node_intervals=apply(new_interval_table[,as.character(child_nodes)],1,function(intervals) {
      new_interval=intersect_intervals(intervals[[1]],intervals[[2]])
      if(is.null(new_interval)) new_interval=list(min(intervals[[1]][[2]],intervals[[2]][[2]]),max(intervals[[1]][[1]],intervals[[2]][[1]]))
      return(new_interval)
    })
    new_interval_table=cbind(current_node_intervals,new_interval_table)
    colnames(new_interval_table)[1]=node
    return(new_interval_table)
  }
}

intersect_intervals=function(interval_a,interval_b){
  if ((interval_a[[2]]>=interval_b[[1]]) && (interval_a[[2]]<=interval_b[[2]])) return(list(max(interval_b[[1]],interval_a[[1]]),interval_a[[2]]))
  else if((interval_a[[2]]>=interval_b[[2]]) && (interval_a[[1]]<=interval_b[[2]])) return(list(max(interval_b[[1]],interval_a[[1]]),interval_b[[2]]))
  else return(NULL)
}

#' Tree reconstruction via maximum parsimony
#'
#' Determines the best fitting tree topology by maximum parsimony for given data
#' @param states_table Matrix or data frame containing numericals representing an alignment of states with sites as rows and species/strains as columns
#' @param rooted_trees List of rooted tree topologies to be analyzed
#' @export
#' @examples
#' test
maximize_clock_tree_log_likelihood=function(states_alignment,tree,Q,pi,initial_depth,is_states_alignment_compact=F){
  initial_branch_fractions=initialize_clock_tree_branch_fractions(tree$edge[1,1],tree)
  newTree=tree
  if (!is_states_alignment_compact) states_alignment=make_states_alignment_compact(states_alignment[,tree$tip.label])
  optimal_branch_lengths=optim(c(initial_branch_fractions,depth=initial_depth),method="L-BFGS-B",lower=c(rep(0.001,length(initial_branch_fractions)),0.1),upper=c(rep(0.999,length(initial_branch_fractions)),10000), function(params) {
    tree=make_tree_from_branch_fractions(tree,params[1:(length(params)-1)],params[length(params)])
    newTree<<-get_tree_log_likelihood(states_alignment,tree,Q,pi,is_states_alignment_compact=T)
    return(-newTree$lnL)
  })
  return(list(tree=newTree,lnL=-optimal_branch_lengths$value,convergence=(optimal_branch_lengths$convergence==0)))
}

make_tree_from_branch_fractions=function(tree,branch_fractions,depth){
  tree$edge.length=branch_fractions_to_lengths(tree$edge[1,1],tree,branch_fractions,depth)
  return(tree)
}

branch_fractions_to_lengths=function(node,tree,branch_fractions,depth){
  child_node_edge_indices=which(tree$edge[,1]==node)
  ret=c()
  sapply(simplify = F,child_node_edge_indices,function(child_node_edge_index){
    child_node=tree$edge[child_node_edge_index,2]
    if (child_node<=length(tree$tip.label)) ret[as.character(child_node_edge_index)]<<-depth
    else{
      ret[as.character(child_node_edge_index)]<<-depth*branch_fractions[as.character(child_node_edge_index)]
      ret<<-c(ret,branch_fractions_to_lengths(child_node,tree,branch_fractions,depth-ret[as.character(child_node_edge_index)]))
    }
  })
  return(ret)
}

initialize_clock_tree_branch_fractions=function(node,tree){
  new_parameter_list=c()
  child_node_edge_indices=which(tree$edge[,1]==node)
  sapply(simplify = F,child_node_edge_indices,function(child_node_edge_index){
    child_node=tree$edge[child_node_edge_index,2]
    if (child_node>length(tree$tip.label)){
      old_parameter_list=initialize_clock_tree_branch_fractions(child_node,tree)
      if (length(old_parameter_list)==0) h=1
      else h=min(old_parameter_list)
      new_parameter_list[as.character(child_node_edge_index)]<<-1/(1/h+1)
      new_parameter_list<<-c(new_parameter_list,old_parameter_list)
    }
  })
  return(new_parameter_list)
}

#' @export
exchange_states=function(alignment,exchange_prob,exchange_dist,gene_info=NULL,overwrite=F){
  if(is.null(gene_info)) gene_info=rep(1,nrow(alignment))
  na.omit(do.call(rbind,sapply(simplify=F,unique(gene_info), function(gene) {
    gene_states_alignment=alignment[which(gene_info==gene),]
    if(nrow(gene_states_alignment)<2) return(NA)
    exchange_table_1=matrix(data=sample(c("exchange","no_exchange"),nrow(gene_states_alignment)*ncol(gene_states_alignment),prob=c(exchange_prob,1-exchange_prob),replace=T),nrow=nrow(gene_states_alignment),ncol=ncol(gene_states_alignment),dimnames=list(1:nrow(gene_states_alignment),colnames(gene_states_alignment)))
    exchange_table=matrix(data=F,nrow=nrow(gene_states_alignment),ncol=ncol(gene_states_alignment),dimnames = list(rownames(gene_states_alignment),colnames(gene_states_alignment)))
    exchanged_table=matrix(data=NA,nrow=nrow(gene_states_alignment),ncol=ncol(gene_states_alignment),dimnames = list(rownames(gene_states_alignment),colnames(gene_states_alignment)))
    sapply(1:nrow(gene_states_alignment), function(row) sapply(1:ncol(gene_states_alignment), function(col){
      if (is.na(exchanged_table[row,col]) && exchange_table_1[row,col]=="exchange"){
        possible_exchange_rows=unique(c(row,intersect(max(row-exchange_dist,1):(min(row+exchange_dist,nrow(gene_states_alignment))),which(!exchange_table[,col]))))
        exchange_row=sample(possible_exchange_rows,1)
        exchanged_table[row,col]<<-gene_states_alignment[exchange_row,col]
        if (!overwrite) {
          exchanged_table[exchange_row,col]<<-gene_states_alignment[row,col]
          exchange_table[row,col]<<-T
          exchange_table[exchange_row,col]<<-T
        }
      } else exchanged_table[row,col]<<-gene_states_alignment[row,col]
    }))
    exchanged_table
  })))
}

#' Generation of a states alignment
#'
#' This function simulates evolution given an evolutionary model and a phylogenetic tree resulting in the generation of a states table (alignment). It is possible to add a noise to the final result that could represent, e.g., sequencing errors.
#' @param nstates Number of sites to be simulated (= number of rows of resulting states table).
#' @param tree Phylogenetic tree of class phylo used for simulation of evolution (branch lengths required!).
#' @param Q Numerical n x n matrix representing the substittion rate matrix of the evolutionary model to be used, with n being the number of states in the model.
#' @param pi Numerical vector of size n representing the equilibrium frequency, with n being the number of states in the model.
#' @export
#' @examples
#' test
simulate_evolution=function(nstates,tree,Q,pi,noise_sd=0,discretization=NA,small_num=0.0001){
  tree$edge.P=sapply(simplify = F,tree$edge.length,function(t) expm(t*Q))
  root_states_seq=sample(1:length(pi),nstates,replace=T,prob=pi)
  tree.root=tree$edge[1,1]
  tree$states_alignment=simulate_evolution_node(tree.root,tree,1:length(pi),root_states_seq)
  tree$states_alignment=tree$states_alignment[,c(tree$tip.label,setdiff(colnames(tree$states_alignment),tree$tip.label))]
    if (is.na(discretization)){
      tree$states_alignment_noise=tree$states_alignment+round(rnorm(nrow(tree$states_alignment)*ncol(tree$states_alignment),0,noise_sd))
      tree$states_alignment_noise=ifelse(tree$states_alignment_noise<1,1,tree$states_alignment_noise)
      tree$states_alignment_noise=ifelse(tree$states_alignment_noise>length(pi),length(pi),tree$states_alignment_noise)
    } else {
      old_discretization=discretization
      sapply(1:length(discretization), function (x) discretization[[x]]<<-ifelse(discretization[[x]]<0,0,discretization[[x]]))
      sapply(1:length(discretization), function (x) discretization[[x]]<<-ifelse(discretization[[x]]>1,1,discretization[[x]]))
      sapply(2:length(discretization), function (x) discretization[[x]][1]<<-ifelse(discretization[[x]][1]==discretization[[x-1]][2],discretization[[x]][1]+small_num,discretization[[x]][1]))
      tree$frequency_alignment=apply(tree$states_alignment,c(1,2),function(x) runif(1,discretization[[x]][1],discretization[[x]][2]))
      tree$frequency_alignment_noise=tree$frequency_alignment+rnorm(nrow(tree$frequency_alignment)*ncol(tree$frequency_alignment),0,noise_sd)
      tree$frequency_alignment_noise=ifelse(tree$frequency_alignment_noise<0,0,tree$frequency_alignment_noise)
      tree$frequency_alignment_noise=ifelse(tree$frequency_alignment_noise>1,1,tree$frequency_alignment_noise)
      tree$states_alignment_noise=apply(tree$frequency_alignment_noise,c(1,2),function(x) which(sapply(old_discretization, function(y) x>y[1] && x<=y[2])))
    }
  return(tree)
}

simulate_evolution_node=function(node,tree,states,states_seq){
  if (node<=length(tree$tip.label)) return(matrix(states_seq,dimnames=list(NULL,tree$tip.label[node])))
  else{
    child_node_edge_indices=which(tree$edge[,1]==node)
    ret=cbind(states_seq,do.call(cbind,sapply(simplify = F,child_node_edge_indices,function(child_node_edge_index){
      child_node=tree$edge[child_node_edge_index,2]
      child_node_states_seq=sapply(states_seq,function(state) sample(states,1,replace=T,prob=tree$edge.P[[child_node_edge_index]][state,]))
      simulate_evolution_node(child_node,tree,states,child_node_states_seq)
    })))
    colnames(ret)[1]=node
    return(ret)
  }
}

#' @export
make_states_alignment_compact=function(states_alignment){
  setDTthreads(1)
  states_alignment=data.table(states_alignment)
  as.matrix(states_alignment[, .(COUNT = .N), by = names(states_alignment)])
}

#' @export
get_tree_log_likelihood=function(states_alignment,tree,Q,pi,is_states_alignment_compact=F){
  if (!is_states_alignment_compact) tree$states_alignment=make_states_alignment_compact(states_alignment[,tree$tip.label])
  else tree$states_alignment=states_alignment
  tree$edge.P=sapply(simplify = F,tree$edge.length,function(t) expm(t*Q))
  tree$states_alignment=cbind(tree$states_alignment,lnL=apply(tree$states_alignment,1,function (homolog_states){
    tree$tip.states=homolog_states[1:(length(homolog_states)-1)]
    tree.root=tree$edge[1,1]
    root_likelihoods=get_node_likelihoods(tree.root,tree,1:length(pi))
    log(sum(sapply(1:length(pi), function(root_state) pi[root_state]*root_likelihoods[root_state])))
  }))
  tree$states_alignment=cbind(tree$states_alignment,`lnL*COUNT`=tree$states_alignment[,"lnL"]*tree$states_alignment[,"COUNT"])
  tree$lnL=sum(tree$states_alignment[,"lnL*COUNT"])
  return(tree)
}

get_node_likelihoods=function(node,tree,states){
  if (node<=length(tree$tip.label)) return(sapply(states,function(state) as.numeric(state==tree$tip.states[node])))
  else {
    child_node_infos=sapply(simplify = F,which(tree$edge[,1]==node),function(child_node_edge_index) {
      child_node=tree$edge[child_node_edge_index,2]
      child_likelihoods=get_node_likelihoods(child_node,tree,states)
      list(child_node_edge_index=child_node_edge_index,child_node=child_node,child_likelihoods=child_likelihoods)
    })
    return(sapply(states,function(state_i){
      prod(sapply(child_node_infos,function(child_node_info) {
        sum(sapply(states, function(state_j) tree$edge.P[[child_node_info$child_node_edge_index]][state_i,state_j]*child_node_info$child_likelihoods[state_j]))
      }))
    }))
  }
}

#' @export
maximize_tree_log_likelihood=function(states_alignment,tree,Q,pi,is_states_alignment_compact=F){
  newTree=tree
  if (!is_states_alignment_compact) states_alignment=make_states_alignment_compact(states_alignment[,tree$tip.label])
  optimal_branch_lengths=optim(control=list(factr=1e10),tree$edge.length,method="L-BFGS-B",lower=rep(0.1,nrow(tree$edge)),upper=rep(100000,nrow(tree$edge)), function(branch_lengths) {
    tree$edge.length=branch_lengths
    newTree<<-get_tree_log_likelihood(states_alignment,tree,Q,pi,is_states_alignment_compact=T)
    return(-newTree$lnL)
  })
  return(list(tree=newTree,lnL=-optimal_branch_lengths$value,convergence=(optimal_branch_lengths$convergence==0)))
}

#' @export
maximize_tree_log_likelihood_extended=function(states_alignment,tree,Q,pi,branches_to_evolutionary_rate_classes=NULL,is_states_alignment_compact=F){
  newTree=tree
  if (is.null(branches_to_evolutionary_rate_classes)) branches_to_evolutionary_rate_classes=1:length(tree$edge.length)
  else {
    branches_to_evolutionary_rate_classes=as.integer(branches_to_evolutionary_rate_classes)
    if (!(sum(is.na(branches_to_evolutionary_rate_classes)==0) && length(branches_to_evolutionary_rate_classes)==length(tree$edge.length))) stop("branches_to_evolutionary_rate_classes must be an integer and has to contain as many elemements as tree has branches")
  }
  if (!is_states_alignment_compact) states_alignment=make_states_alignment_compact(states_alignment[,tree$tip.label])
  evolutionary_rate_classes=unique(branches_to_evolutionary_rate_classes)
  evolutionary_rates=rep(1,length(evolutionary_rate_classes))
  ret=optim(control=list(factr=1e9),evolutionary_rates,method="L-BFGS-B",lower=tree$edge.length/100000,upper=tree$edge.length/0.1, function(evolutionary_rates) {
    tree$edge.length=sapply(1:length(tree$edge.length),function(i) tree$edge.length[i]*evolutionary_rates[branches_to_evolutionary_rate_classes[i]])
    newTree<<-get_tree_log_likelihood(states_alignment,tree,Q,pi,is_states_alignment_compact=T)
    return(-newTree$lnL)
  })
  return(list(tree=newTree,lnL=-ret$value,convergence=(ret$convergence==0),message=ret$message,evolutionary_rates=ret$par))
}

#' @export
find_each_optimal_tree=function(states_alignment,trees,Q,pi,cluster=NULL,clock=F,initial_depth=1){
  if (is.null(cluster)) sapply(simplify = F, trees, function(tree)
    if (clock) maximize_clock_tree_log_likelihood(states_alignment,tree,Q,pi,initial_depth)
    else maximize_tree_log_likelihood(states_alignment,tree,Q,pi)
  )
  else {
    clusterExport(cluster,c("get_tree_log_likelihood","get_node_likelihoods","maximize_tree_log_likelihood","maximize_clock_tree_log_likelihood","initialize_clock_tree_branch_fractions","branch_fractions_to_lengths","make_tree_from_branch_fractions","Q","pi","states_alignment"),envir=environment())
    if (clock) parSapply(cluster,simplify = F, trees, function(tree) maximize_clock_tree_log_likelihood(states_alignment,tree,Q,pi,initial_depth))
    else parSapply(cluster,simplify = F, trees, function(tree) maximize_tree_log_likelihood(states_alignment,tree,Q,pi))
  }
}

#' @export
find_optimal_tree=function(states_alignment,trees,Q,pi,cluster=NULL,clock=F,initial_depth=1){
  best_tree_results=find_each_optimal_tree(states_alignment,trees,Q,pi,cluster,clock,initial_depth)
  best_tree_results[[which.max(sapply(best_tree_results, function(best_tree_result) best_tree_result$lnL))]]$tree
}

#' @export
all_unrooted_tree_topologies=function(leaves,i=NULL,baseTree=NULL,branch_length=1){
  if (is.null(i) || (i<=3) || is.null(baseTree)) {
    baseTree=stree(3,"star",leaves[1:3])
    i=3
  }
  if (i<length(leaves)) {
    newBaseTrees=sapply(simplify = F, 1:nrow(baseTree$edge),function(edge_i) {
      newBaseTree=baseTree
      newBaseTree$tip.label=leaves[1:(i+1)]
      newBaseTree$Nnode=newBaseTree$Nnode+1
      newBaseTree$edge=ifelse(baseTree$edge>i,baseTree$edge+1,baseTree$edge)
      if (edge_i==1) pre=c() else pre=newBaseTree$edge[1:(edge_i-1),]
      if (edge_i==nrow(baseTree$edge)) post=c() else post=newBaseTree$edge[(edge_i+1):nrow(baseTree$edge),]
      newBaseTree$edge=rbind(pre,
                             c(newBaseTree$edge[edge_i,1],newBaseTree$Nnode+i+1),
                             c(newBaseTree$Nnode+i+1,newBaseTree$edge[edge_i,2]),
                             c(newBaseTree$Nnode+i+1,i+1),
                             post)
      newBaseTree$edge.length=rep(branch_length,nrow(newBaseTree$edge))
      return(newBaseTree)
    })
    unlist(sapply(simplify = F,newBaseTrees,function(newBaseTree) all_unrooted_tree_topologies(leaves,i+1,newBaseTree,branch_length)),recursive = F)
  } else list(baseTree)
}

#' @export
all_rooted_tree_topologies=function(leaves,branch_length=1){
  trees=all_unrooted_tree_topologies(c(leaves,"root_indicator"),branch_length=branch_length)
  sapply(simplify = F, trees, function(tree){
    root_indicator_index=which(tree$edge[,2]==length(tree$tip.label))
    tree=root(tree,node=tree$edge[root_indicator_index,1])
    root_indicator_index=which(tree$edge[,2]==length(tree$tip.label))
    tree$tip.label=tree$tip.label[1:(length(tree$tip.label)-1)]
    tree$edge.length=tree$edge.length[1:(length(tree$edge.length)-1)]
    tree$edge=tree$edge[-root_indicator_index,]
    tree$edge=ifelse(tree$edge>length(tree$tip.label),tree$edge-1,tree$edge)
    return(tree)
  })
}

#' @export
maximize_distance_log_likelihoods=function(states_table,Q,pi){
  species_pairs=combn(colnames(states_table),2)
  colnames(species_pairs)=apply(species_pairs,2,function(x) paste0(x,collapse = "_"))
  distance_matrix=matrix(0,nrow=ncol(states_table),ncol=ncol(states_table),dimnames=list(colnames(states_table),colnames(states_table)))
  h=apply(species_pairs,2,function(species_pair)
    distance_matrix[species_pair[1],species_pair[2]]<<-distance_matrix[species_pair[2],species_pair[1]]<<-optim(1,function(t)
      get_distance_log_likelihood(states_table[,species_pair[1]],states_table[,species_pair[2]],Q,t,pi),method="L-BFGS-B",lower=0.001,upper=100000)$par
  )
  #distance_matrix=distance_matrix*2
  return(distance_matrix)
}

#' @export
get_distance_log_likelihood=function(statesA,statesB,Q,t,pi){
  P=expm(Q*t)
  -sum(log(sapply(1:length(statesA),function(pos) pi[statesA[pos]]*P[statesA[pos],statesB[pos]] )))
}

#' @export
tree_to_distances=function(tree){
  node_chains=sapply(simplify = F, tree$tip.label,function(species){
    cur_edge_idx=which(tree$edge[,2]==which(tree$tip.label==species))
    ret=c(cur_edge_idx)
    while(tree$edge[cur_edge_idx,1]!=tree$edge[1,1]) {cur_edge_idx=which(tree$edge[,2]==tree$edge[cur_edge_idx,1]);ret=c(cur_edge_idx,ret)}
    ret
  })
  species_pairs=combn(tree$tip.label,2)
  colnames(species_pairs)=apply(species_pairs,2,function(x) paste0(x,collapse = "_"))
  distances=apply(species_pairs,2,function(pair) {
    edge_overlap=intersect(node_chains[[pair[1]]],node_chains[[pair[2]]])
    sum(tree$edge.length[setdiff(node_chains[[pair[1]]],edge_overlap)])+sum(tree$edge.length[setdiff(node_chains[[pair[2]]],edge_overlap)])
  })
  return(distances)
}

#' @export
tree_to_distance_matrix=function(tree){
  node_chains=sapply(simplify = F, tree$tip.label,function(species){
    cur_edge_idx=which(tree$edge[,2]==which(tree$tip.label==species))
    ret=c(cur_edge_idx)
    while(tree$edge[cur_edge_idx,1]!=tree$edge[1,1]) {cur_edge_idx=which(tree$edge[,2]==tree$edge[cur_edge_idx,1]);ret=c(cur_edge_idx,ret)}
    ret
  })
  distance_matrix=matrix(0,nrow=length(tree$tip.label),ncol=length(tree$tip.label),dimnames=list(tree$tip.label,tree$tip.label))
  species_pairs=combn(tree$tip.label,2)
  colnames(species_pairs)=apply(species_pairs,2,function(x) paste0(x,collapse = "_"))
  apply(species_pairs,2,function(pair) {
    edge_overlap=intersect(node_chains[[pair[1]]],node_chains[[pair[2]]])
    distance_matrix[pair[1],pair[2]]<<-distance_matrix[pair[2],pair[1]]<<-sum(tree$edge.length[setdiff(node_chains[[pair[1]]],edge_overlap)])+sum(tree$edge.length[setdiff(node_chains[[pair[2]]],edge_overlap)])
  })
  return(distance_matrix)
}

#' @export
best_fitch_margoliash=function(distance_matrix,trees){
  trees=all_fitch_margoliash(distance_matrix,trees)
  trees[[which.min(sapply(simplify=F,trees,function(tree) tree$deviation))]]$tree
}

#' @export
all_fitch_margoliash=function(distance_matrix,trees) sapply(simplify=F,trees,function(tree) fitch_margoliash(tree,distance_matrix))

#' @export
fitch_margoliash=function(tree,distance_matrix){
  if ("edge.length" %in% names(tree))initial_lengths=tree$edge.length
  else initial_lengths=rep(100,nrow(tree$edge))
  newTree=tree
  optimal_branch_lengths=optim(initial_lengths,method="L-BFGS-B",lower=0.1,upper=100000,function(branch_lengths){
    newTree$edge.length<<-branch_lengths
    deviation=tree_to_distance_matrix(newTree)-distance_matrix
    return(sum(deviation^2))
  })
  return(list(tree=newTree,deviation=optimal_branch_lengths$value^0.5,convergence=(optimal_branch_lengths$convergence==0)))
}

#' @export
discretize=function(frequency_alignment,discretization) apply(frequency_alignment,c(1,2), function(x) which(sapply(discretization, function(y) x>y[1] && x<=y[2])))

#' @export
make_evolutionary_model=function(states_table=NULL,nstates=NULL,model="JC69",pi=NULL,kappa=NULL){
  model=tolower(model)
  if(model!="jc69" && model!="k80" && model!="f81" && model!="hky85" && model!="cooc" && model!="nojump") stop("\"",paste0(model, "\" is no valid model."))
  if (is.null(pi)) {
    states_distr_abs=sapply(colnames(states_table), function(species) as.numeric(table(states_table[,species])))
    states_distr_rel=apply(states_distr_abs,2,function(x) x/sum(x))
    pi=apply(states_distr_rel,1,function(x) sum(x)/ncol(states_distr_rel))
  }
  if ((model=="k80" || model=="hky85") && is.null(kappa)) kappa=estimate_kappa(states_table)
  if (model=="nojump"){
    Q=sapply(1:nstates,function(x) sapply(1:nstates, function(y)  if (abs(x-y)==1) pi[x] else 0))
  } else if (model=="cooc"){
    species_pairs=combn(colnames(states_table),2)
    colnames(species_pairs)=apply(species_pairs,2,function(x) paste0(x,collapse = "_"))
    cooccurrence_matrices=sapply(simplify=F, 1:ncol(species_pairs),function(x) matrix(0,nrow=nstates,ncol=nstates))
    names(cooccurrence_matrices)=colnames(species_pairs)
    h=apply(states_table,1, function(line){
      sapply(colnames(species_pairs), function(species_pair)
        cooccurrence_matrices[[species_pair]][line[species_pairs[1,species_pair]],line[species_pairs[2,species_pair]]]<<-cooccurrence_matrices[[species_pair]][line[species_pairs[1,species_pair]],line[species_pairs[2,species_pair]]]+1 )
    })
    
    cooccurrence_matrices_balanced=cooccurrence_matrices
    h=sapply(1:length(cooccurrence_matrices_balanced), function(i){ 
      sapply(1:(nstates-1), function(a) {
        sapply((a+1):nstates, function(b) cooccurrence_matrices_balanced[[i]][a,b]<<-cooccurrence_matrices_balanced[[i]][b,a]<<-cooccurrence_matrices_balanced[[i]][a,b]+cooccurrence_matrices_balanced[[i]][b,a] )
        cooccurrence_matrices_balanced[[i]][a,a]<<-0
      })
      cooccurrence_matrices_balanced[[i]][nstates,nstates]<<-0
    })  
    
    cooccurrence_matrices_balanced_frac=sapply(simplify = F, cooccurrence_matrices_balanced, function(x) x/sum(x)*2)
    cooccurrence_matrices_balanced_frac_all=Reduce("+",cooccurrence_matrices_balanced_frac)/length(cooccurrence_matrices_balanced_frac)
    
    Q=sapply(1:nstates,function(x) sapply(1:nstates, function(y)  pi[x]*cooccurrence_matrices_balanced_frac_all[x,y]))
  } else if (model=="jc69") Q=matrix(c(0,1,1,1,1,0,1,1,1,1,0,1,1,1,1,0),ncol=4,dimnames=list(c("A","C","G","T"),c("A","C","G","T")))
  else if (model=="k80") Q=matrix(c(0,1,kappa,1,1,0,1,kappa,kappa,1,0,1,1,kappa,1,0),ncol=4,dimnames=list(c("A","C","G","T"),c("A","C","G","T")))
  else if (model=="f81") Q=matrix(c(0,pi[2],pi[3],pi[4],pi[1],0,pi[3],pi[4],pi[1],pi[2],0,pi[4],pi[1],pi[2],pi[3],0),ncol=4,dimnames=list(c("A","C","G","T"),c("A","C","G","T")))
  else if (model=="hky85") Q=matrix(c(0,pi[2],pi[3]*kappa,pi[4],pi[1],0,pi[3],pi[4]*kappa,pi[1]*kappa,pi[2],0,pi[4],pi[1],pi[2]*kappa,pi[3],0),ncol=4,dimnames=list(c("A","C","G","T"),c("A","C","G","T")))
  Q=Q/sum(rowSums(Q)*pi)/100#PEM calibration
  sapply(1:nrow(Q), function (x) Q[x,x]<<- -sum(Q[x,]))
  return(list(Q=Q,pi=pi))
}

estimate_kappa=function(states_table){
  species_pairs=combn(colnames(states_table),2)
  transitions=0
  transversions=0
  apply(species_pairs,2,function(species_pair) apply(states_table[,species_pair],1,function(x) if (x[1]!=x[2]){
    x=sort(x)
    if ((x[1]==1 && x[2]==3) || (x[1]==2 && x[2]==4) ) transitions<<-transitions+1
    else transversions<<-transversions+1
  }))
  2*transitions/transversions
}