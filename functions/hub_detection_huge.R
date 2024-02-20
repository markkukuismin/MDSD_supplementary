
hub_detection_huge = function(huge_est = NULL, plot = TRUE){
  
  if(is.null(huge_est)) stop("Error: huge_est is missing")
  
  p = ncol(huge_est$data)
  
  node_names = colnames(huge_est$data)
  
  if(is.null(node_names)) node_names = 1:p
  
  lambda = huge_est$lambda
  
  nlambda = length(lambda)
  
  Degree_lambda = matrix(0, nrow = p, ncol = nlambda)
  
  colnames(Degree_lambda) = lambda
  rownames(Degree_lambda) = node_names
  
  for(j in 1:nlambda){
    
    d_temp = colSums(as.matrix(huge_est$path[[j]]))
    
    names(d_temp) = node_names
    
    Degree_lambda[node_names, j] = d_temp
    
  }
  
  BC = matrix(0, nrow = p, ncol = nlambda)
  
  for(j in 1:nlambda){
    
    G = huge_est$path[[j]]
    
    G = as.matrix(G)
    
    G = igraph::graph_from_adjacency_matrix(G, mode = "undirected", diag = FALSE)
    
    BC[ , j] = igraph::betweenness(G, directed = FALSE)
    
  }
  
  bcss = rep(0, p)
  
  for(j in 1:p){
    
    bcss[j] = mean((BC[rep(j, p - 1), ] - BC[-j, ])^2)
    
  }
  
  degrees = c(t(Degree_lambda))
  
  BC = c(t(BC))
  
  lambdas = rep(lambda, times = p)
  
  all_nodes = rep(node_names, each = ncol(Degree_lambda))
  
  MDSD = rep(0, p)
  
  for(j in 1:p){
    
    MDSD[j] = mean((Degree_lambda[rep(j, p - 1), ] - Degree_lambda[-j, ])^2)
    
  }
  
  names(MDSD) = rownames(Degree_lambda)
  
  names(bcss) = rownames(Degree_lambda)
  
  hub_nodes_MDSD = node_names[MDSD > 3*mean(MDSD)]
  
  hub_nodes_BC = node_names[bcss > 3*mean(bcss)]
  
  raw_data = data.frame(node = all_nodes,
                        Degree = degrees,
                        BC = BC,
                        lambda = lambdas
                        )
  
  sol_path_data = data.frame(node = node_names,
                             MDSD = MDSD,
                             bcss = bcss
                             )
  
  p1 = ggplot(raw_data, 
              aes(x = lambda, y = Degree, group = node)
  ) +
    geom_line() +
    ylab("Node degree")
  
  p2 = ggplot(raw_data, 
              aes(x = lambda, y = BC, group = node)
  ) +
    geom_line() +
    ylab("Betweenness centrality")
  
  p3 = ggplot(sol_path_data, aes(x = node, y = MDSD)) +
    geom_segment(
      aes(x = node, xend = node, y = 0, yend = MDSD)
    ) +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    ) +
    ylab("Mean Degree Squared Difference of degree path")
  
  p4 = ggplot(sol_path_data, aes(x = node, y = bcss)) +
    geom_segment(
      aes(x = node, xend = node, y = 0, yend = bcss)
      ) +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
      ) +
    ylab("Mean Degree Squared Difference of betweenness centrality path")
  
  if(!plot) p1 = NULL
  if(!plot) p2 = NULL
  if(!plot) p3 = NULL
  if(!plot) p4 = NULL
  
  results = list(bcss= bcss,
                 lambda = lambda,
                 hub_nodes_MDSD = hub_nodes_MDSD,
                 hub_nodes_BC = hub_nodes_BC,
                 raw_data = raw_data,
                 sol_path_data = sol_path_data,
                 degree_plot = p1,
                 BC_plot = p2,
                 MDSD_plot = p3,
                 BC_AUC_plot = p4
                 )
  
  return(results)
  
}