require('igraph')
require('Matrix')

# The algorithm for generating cascades works like follows: 
# 1.    For each node n and each factor k we generate S[n,k] and A[n,k]
# 2.    We randomly choose the initial nodes for each cascade
# 3.    for each cascade to generate, we proceed like follows: 
# 3.1   randomly pick a node n among those active nodes and remove it from the list of active
# 3.2   try to activate all his neighbors
# 3.3   for each neighbor m
# 3.3.1 generate a delay D using S[m,k]*A[n,k]. Add the delay to the activation time of n to obtain the activation time of m
# 3.3.2 if the activation time is greater than Tmax, then we skip it
# 3.3.3 otherwise we add it to the list of active nodes
# 4     the algorithm terminates when the list of active nodes is empty
generate_cascades <- function(g,numNodes,numWords,numCascades,numFactors,Tmax,Tmin,S, A, Phi){
  cascades = Matrix(0,
                    nrow = numCascades,
                    ncol = numNodes,
                    sparse = TRUE)
  
  
  content = Matrix(0,
                    nrow = numCascades,
                    ncol = numWords,
                    sparse = TRUE)
#  initial_nodes = sample(numNodes,numCascades, replace = TRUE,prob = in.d.g)
  totnodes = 1:numNodes
  
  
  assignments = rep(0,numCascades)
  
  for (c in 1:numCascades) {
    triggers = rep(0, numNodes)
    K = sample(numFactors, 1)
    assignments[c] = K
    
    # Generate content by sampling from a Poisson
    content[c,] = rpois(numWords, Phi[,K])
    
    #Generate the cascade     
#    currnode = initial_nodes[c]
    currnode = sample(numNodes,1,prob = A[,K]/sum(A[,K]))
    
    
    timestamp = runif(1, 0, Tmin)
    cascades[c, currnode] = timestamp
    canExpand = TRUE
    while (canExpand) {
      cond = cascades[c,] > 0 & triggers == 0
      currnodes = totnodes[cond]
      if (length(currnodes) > 0) {
        if (length(currnodes) > 1) {
          p = A[currnodes, K] / sum(A[currnodes, K])
          currnode = sample(currnodes, size = 1, prob = p)
        } else {
          currnode = currnodes[1]
        }
        triggers[currnode] = 1
        nb = as.numeric(neighbors(g, currnode, mode = "in"))
        for (nextnode in nb) {
          if (cascades[c, nextnode] == 0) {
            rate = S[nextnode, K] * A[currnode, K]
            
            timestamp = cascades[c, currnode] + rexp(1, rate)
            if (timestamp < Tmax) {
              cascades[c, nextnode] <- timestamp
            }
          }
        }
      } else {
        canExpand = FALSE
      }
    }
  }
  output<-list(assignments,cascades,content)
  return(output)
}


# Generate the expected frequency vectors for each word and each topic
# words are randomly assigned to a topic, with a bias for non-assignment.
# For each topic, we sample the expected frequency of the related words
# from Gamma(tshape, tscale)
# The remaining expected frequencies are sampled from Gamma(blshape, blscale)
generateFrequencyVecs <- function(numFactors,numWords, bias, blshape, blscale, tshape, tscale){
  Phi = matrix(rgamma(numWords * numFactors, shape = blshape, scale = blscale), ncol = numFactors)
  
  probs = c(bias,rep((1-bias)/numFactors, numFactors))
  assignments <- sample(0:numFactors,numWords,replace=TRUE, prob=probs)
  
  for (K in 1:numFactors){
    active_nodes = assignments==K
    Phi[active_nodes,K] = rgamma(sum(active_nodes), shape = tshape, scale = tscale)
  }
  return(Phi)
}

# Generate the influence and susceptibility matrices. 
# For each topic, we choose a given number of nodes and assign high authoritativeness in that topic. 
# All nodes with outgoing edges to these nodes are given high susceptibility in that topic. 
# The remaining nodes are given default (low) susceptibility and authoritativeness
generateInfluenceMatrices <- function(g,numNodes,numFactors){
  # Generate the matrix with default (low) values
  S = matrix(runif(numNodes * numFactors, min = 0.001, max = 0.01), ncol = numFactors)
  
  A = matrix(runif(numNodes * numFactors, min = 0.001, max = 0.01), ncol = numFactors)
  
  for (k in 1:numFactors){
    #For each factor we choose a number of nodes whose influence is boosted
    # The nodes are split into batches 
    numNodesToChoose = ceiling(numNodes/(2*numFactors))
    if (numNodesToChoose == 0){
      numNodesToChoose = 1
    }
    nodesToChoose = sample(numNodes,numNodesToChoose)
    
    # boost the influence degree for those nodes
    A[nodesToChoose,k] = runif(numNodesToChoose, min = 0.1, max = 1)
    
    # Collect all incoming nodes
    neighborsOfNodesToChoose = c() 
    for (i in 1:length(nodesToChoose)) {
      neighborsOfNodesToChoose = c(neighborsOfNodesToChoose,as.numeric(neighbors(g,nodesToChoose[i],mode="in")))
    }
    neighborsOfNodesToChoose = unique(neighborsOfNodesToChoose)
  
    # and boost their susceptibility as well
    if (length(neighborsOfNodesToChoose) > 0){
      S[neighborsOfNodesToChoose,k] = runif(length(neighborsOfNodesToChoose), min = 0.1, max = 1)
    }
    
  }
  
  output<-list(S,A)
  return(output)
}


# Alternate method for generating influence and susceptibility matrices
# The strength of each link is determined by considering both the out-degree (k-out) 
# of the source and the in-degree (k-in) of the destination: 
#
#  weight(u, v) propto lambda * (1 - k-out(u)/k-outmax) * k-in(v)/k-inmax + (1 − lambda) * rand(0.1, 1)
#  
#  # and κin are the maximum out-degree and in-degree respectively, and lambda
# is used to introduce a random effect. The idea is the following: nodes susceptibility 
# of a node within a community is inversely proportional to its outdegree, whereas
#authoritativenes is proportional to its indegree
generateInfluenceMatricesFromCommunity <- function(g, numNodes, numFactors,node_communities,lambda){
   in.d.g = degree(g, mode = "in")
   in.d.g = in.d.g / max(in.d.g)
   
   out.d.g = degree(g, mode = "out")
   out.d.g = out.d.g / max(out.d.g)
  
   maxval = min(out.d.g[out.d.g>0],in.d.g[in.d.g>0])*10^-1
   minval = maxval*10^-1
   # Generate the matrix with default (low) values
   S = matrix(runif(numNodes * numFactors, min = minval, max = maxval), ncol = numFactors)
   A = matrix(runif(numNodes * numFactors, min = minval, max = maxval), ncol = numFactors)
   
   sz = dim(node_communities)
  for (i in 1:sz[1]){
    K = node_communities[i,2]
    node = node_communities[i,1]
    
    #toss the coin based on lambda
    p = rbinom(1,1,prob = lambda)
    if (p==1){
      S[node,K] = max(0.01,rnorm(1,mean = 1-out.d.g,sd = 0.05))
      A[node,K] = in.d.g
      
    } else {
      S[node,K] = runif(1,min= 0.1,max=1)
      A[node,K] = runif(1,min= 0.1,max=1)
    }
  }
  output<-list(S,A)
  return(output)
}





# Remove from the graph all the nodes not appearing in any cascade, and remap the node names. 
clean_graph <- function(g,cascades,S,A){
  # Check the nodes who did not become active in any cascade. These nodes should be removed from the network
  nodes_not_in_cascades = apply(t(cascades), 1, function(row) all(row ==0 ))
  
  # a vector for node renaming
  names = rep(0,numNodes)
  
  # Now we associate a name to each node, by skipping those nodes to remove
  k = 1
  for (i in 1:numNodes){
    if (!nodes_not_in_cascades[i]){
      names[i] = k
      k = k + 1
    }
  }
  
  # the new graph now has new names
  gm <- set.vertex.attribute(g, "name", value=names)
  
  # and the nodes to remove are deleted from the graph
  ids = 1:numNodes
  gm <-  delete.vertices(gm,V(gm)[ids[nodes_not_in_cascades]])
  
  # and next from the cascade matrix
  newcascades = cascades[,!nodes_not_in_cascades]
  
  Sn = S[!nodes_not_in_cascades,]
  An = S[!nodes_not_in_cascades,]
  
  output<-list(newcascades,gm,Sn,An)
  return(output)
}

