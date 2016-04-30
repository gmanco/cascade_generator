process_data <-
  function(g,
           numFactors,
           numCascades,
           numWords,
           Tmax,
           Tim,
           gamma_params,
           prefix,
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
    if (is.null(node_communities) ){
      out <- generateInfluenceMatrices(g, numNodes, numFactors)
    }
    else {
      out <-
        generateInfluenceMatricesFromCommunity(g, numNodes, numFactors, node_communities)
    }
    
    S = out[[1]]
    A = out[[2]]
    
    
    Phi <-
      generateFrequencyVecs(numFactors,
                            numWords,
                            bias = gamma_params[1],
                            blshape = gamma_params[2],
                            blscale = gamma_params[3],
                            tshape = gamma_params[4],
                            tscale = gamma_params[5])
    
    
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