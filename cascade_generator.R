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

#set the path and choose the file to build cascades upon
setwd("/Users/manco/Dropbox/shared ICAR/Influence_PP/exp/datasets/synth/random")
edgeListFile = "random.txt"

use_communities = FALSE

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

gamma_params = c(bias,blshape,blscale,tshape,tscale)

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
      node_communities <- cbind(1:numNodes, memb[, 2])
    } else{
      # Otherwise read the file into a matrix
      node_communities = read.table(communityFile, sep = "\t", header = F)
    }
    numFactors = max(node_communities[, 2])
    process_data(g, numFactors, numCascades, numWords, Tmax,Tim,gamma_params, prefix,node_communities)
    
  } else {
    ### number of factors is iterated
    for (v in kmin:kmax) {
      numFactors = 2 ^ v
      
      process_data(g, numFactors, numCascades, numWords,  Tmax,Tim,gamma_params,prefix,NULL)
      
    }
    
    
    
    
    process_data <-
      function(g,
               numFactors,
               numCascades,
               numWords,
               prefix,
               Tmax,
               Tim,
               gamma_params,
               node_communities) {
        # Generate the data iteratively
        
        
        filename = paste(prefix,
                         ".nn",
                         numNodes,
                         ".nc",
                         numCascades,
                         ".k",
                         numFactors,
                         sep = "")
        
        
        # Generate the influence matrices to be used for the cascade generation
        if (node_communities == NULL) {
          out <- generateInfluenceMatrices(g, numNodes, numFactors)
        }
        else {
          out <-
            generateInfluenceMatricesFromCommunity(g, numNodes, numFactors, node_communities)
        }
        
        
        Phi <-
          generateFrequencyVecs(numFactors,
                                numWords,
                                gamma_params[1],
                                gamma_params[2],
                                gamma_params[3],
                                gamma_params[4],
                                gamma_params[5])
        
        S = out[[1]]
        A = out[[2]]
        
        # Cascades are generated as a matrix numCascades x numNodes. Each cell contains the activation time. cells with value 0 are not active
        out = generate_cascades(g,
                                numNodes,
                                numWords,
                                numCascades,
                                numFactors,
                                Tmax,
                                Tmin,
                                S,
                                A,
                                Phi)
        
        assignments = out[[1]]
        cascades = out[[2]]
        content = out[[3]]
        # Clean up the graph by removing those nodes not included in any cascade
        out = clean_graph(g, cascades, S, A)
        
        cascades = out[[1]]
        g = out[[2]]
        S = out[[3]]
        A = out[[4]]
        
        
        ## We are done: let's generate the output files
        
        # Cascades (Relational format)
        cascades_relational = paste(filename, ".cascades", sep = "")
        generate_relational_format(cascades,
                                   cascades_relational,
                                   header_names = c("User", "Item", "TimeStamp"))
        
        # Cascades (Netrate format)
        cascade_netrate = paste(filename, ".cascades.netrate", sep = "")
        network_netrate = paste(filename, ".network.netrate", sep = "")
        generate_netrate_format(cascades, g, cascade_netrate, network_netrate)
        
        
        # Content (Relational format)
        content_relational = paste(filename, ".content", sep = "")
        generate_relational_format(content,
                                   content_relational,
                                   header_names = c("Word", "Item", "Frequency"))
        
        # The graph
        
        network_file = paste(filename, ".network", sep = "")
        #write.graph(g,network_file,format="edgelist")
        edges = get.edgelist(g)
        class(edges) <- "numeric"
        write.table(
          edges,
          file = network_file,
          sep = "\t",
          col.names = c("Src", "Dst"),
          row.names = FALSE,
          quote = FALSE
        )
        
        # the graph of two hops (for testing the accuracy)
        network_hops = paste(filename, ".network.2_hops", sep = "")
        generate_two_hops(g, network_hops, 1)
        
        # The matrices
        a_file = network_file = paste(filename, ".A", sep = "")
        s_file = network_file = paste(filename, ".S", sep = "")
        p_file = network_file = paste(filename, ".P", sep = "")
        write.table(
          A,
          file = a_file,
          sep = "\t",
          col.names = FALSE,
          row.names = FALSE
        )
        write.table(
          S,
          file = s_file,
          sep = "\t",
          col.names = FALSE,
          row.names = FALSE
        )
        write.table(
          Phi,
          file = p_file,
          sep = "\t",
          col.names = FALSE,
          row.names = FALSE
        )
        
        # The assignments
        ass_file = network_file = paste(filename, ".clusters", sep = "")
        write.table(
          assignments,
          file = ass_file,
          sep = "\t",
          col.names = FALSE,
          row.names = FALSE
        )
      }
  }
}
