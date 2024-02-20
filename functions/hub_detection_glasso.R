
hub_detection_glasso = function(glasso_path = NULL, node_names = NULL){
  
  p = ncol(glasso_path$wi)
  
  if(is.null(node_names)) node_names = 1:p
  
  lambda = glasso_path$rholist
  
  nlambda = length(lambda)
  
  Degree_lambda = matrix(0, nrow = p, ncol = nlambda)
  
  colnames(Degree_lambda) = lambda
  rownames(Degree_lambda) = node_names
  
  for(j in 1:nlambda){
    
    A_temp = glasso_path$wi[,, j]
    
    A_temp[A_temp != 0] = 1
    
    diag(A_temp) = 0
    
    d_temp = colSums(A_temp)
    
    names(d_temp) = node_names
    
    Degree_lambda[node_names, j] = d_temp
    
  }
  
  rm(A_temp)
  
  MDSD = rep(0, p)
  
  for(j in 1:p){
    
    MDSD[j] = mean((Degree_lambda[rep(j, p - 1), ] - Degree_lambda[-j, ])^2)
    
  }
  
  names(MDSD) = rownames(Degree_lambda)
  
  hub_nodes_MDSD = node_names[MDSD > 3*mean(MDSD)]
  
  return(list(hub_nodes_MDSD = hub_nodes_MDSD,
              MDSD = MDSD,
              Degree = Degree_lambda))
  
}