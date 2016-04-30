library(igraph)
library(Matrix)
library(gplots)

#### This file generates cascades from a network

#Some utility functions
source(
  "/Users/manco/Dropbox/shared ICAR/Influence_PP/exp/datasets/synth/cascade_generator/generate_cascades.R"
)
source(
  "/Users/manco/Dropbox/shared ICAR/Influence_PP/exp/datasets/synth/cascade_generator/output_format.R"
)
source(
  '/Users/manco/Dropbox/shared ICAR/Influence_PP/exp/datasets/synth/cascade_generator/process_data.R'
)

#set the path and choose the file to build cascades upon
setwd("/Users/manco/Dropbox/shared ICAR/Influence_PP/exp/datasets/synth/random")
edgeListFile = "random.txt"

use_communities = TRUE

prefix = "data/random"
dat = read.table(file = edgeListFile, header = FALSE, sep = "\t") # choose an adjacency matrix from a .csv file
communityFile = ""

#Set generation parameters
#numCascades = 4096
#numFactors = 8
Tmax = 100
Tmin = 5

# Specify the number of datasets to generate
cmin = 9
cmax = 9#15
kmin = 1
kmax = 1#7

# When generating content, we assume the following:
# bias is the percentage of words not associated with any topic
# blshape/blscale represent the shape and scale parameters upon which to sample the baseline mean for the poisson distribution
# tshape/tscale represent the shape and scale parameters upon which to sample the topic-based mean for the poisson distribution
# Assumption: baseline poisson distributions are random, while topic-based poisson distributions are biased
bias = 0.5
blshape = 1
blscale = 1
tshape = 3
tscale = 3

numWords = 1024

gamma_params = c(bias, blshape, blscale, tshape, tscale)

set.seed(0.124532)


for (u in cmin:cmax) {
  numCascades = 2 ^ u
  
  
  # generate the graph
  # NOTICE: the node names are those denoted in the edgeListFile
  g = graph.data.frame(dat, directed = TRUE)
  
  numNodes = vcount(g)
  #OPTIONAL: plot the graph to verify
  #igraph.options(vertex.size=0.2, vertex.label=NA, edge.arrow.size=0.1)
  #plot(g, layout=layout.kamada.kawai)
  
  
  ### when a community detection is enabled the
  ### number of factors is given by the number of communities
  if (use_communities) {
    if (communityFile == "") {
      ## if the community file is not given, generate communities using a
      ## community detection algorithm
      fc <- fastgreedy.community(simplify(as.undirected(g)))
      memb <- as.matrix(membership(fc))
      node_communities <- cbind(1:numNodes, memb)
    } else{
      # Otherwise read the file into a matrix
      node_communities = read.table(communityFile, sep = "\t", header = F)
    }
    prefix = paste(prefix,".communities",sep="")
    numFactors = max(node_communities[, 2])
    process_data(
      g,
      numFactors,
      numCascades,
      numWords,
      Tmax,
      Tim,
      gamma_params,
      prefix,
      node_communities
    )
    
  } else {
    ### number of factors is iterated
    for (v in kmin:kmax) {
      numFactors = 2 ^ v
      
      process_data(g,
                   numFactors,
                   numCascades,
                   numWords,
                   Tmax,
                   Tim,
                   gamma_params,
                   prefix,
                   NULL)
      
    }
  }
}
