
hub_detection_hglasso = function(hglasso_path = NULL, node_names = NULL, gamma = NULL){
  
  if(is.null(gamma)) gamma = 3
  
  p = ncol(hglasso_path$wi)
  
  if(is.null(node_names)) node_names = 1:p
  
  lambda = hglasso_path$rholist
  
  nlambda = length(lambda)
  
  Degree_lambda = matrix(0, nrow = p, ncol = nlambda)
  
  colnames(Degree_lambda) = lambda
  rownames(Degree_lambda) = node_names
  
  for(j in 1:nlambda){
    
    A_temp = hglasso_path$wi[,, j]
    
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
  
  d_skewness = apply(Degree_lambda, 2, e1071::skewness)
  
  d_skewness[is.na(d_skewness)] = 0
  
  burn = which(d_skewness <= 1)
  
  if(!rlang::is_empty(burn)){
    
    MDSD_burn = rep(0, p)
    
    for(j in 1:p){
      
      MDSD_burn[j] = mean((Degree_lambda[rep(j, p - 1), -burn] - Degree_lambda[-j, -burn])^2)
      
    }
    
  }else{
    
    MDSD_burn = MDSD
    
  }
  
  names(MDSD) = rownames(Degree_lambda)
  names(MDSD_burn) = rownames(Degree_lambda)
  
  hub_nodes_MDSD = node_names[MDSD > gamma*mean(MDSD)]
  hub_nodes_MDSD_burn = node_names[MDSD_burn > gamma*mean(MDSD_burn)]
  
  degrees = c(t(Degree_lambda))
  
  return(list(hub_nodes_MDSD = hub_nodes_MDSD,
              hub_nodes_MDSD_burn = hub_nodes_MDSD_burn,
              MDSD = MDSD,
              MDSD_burn = MDSD_burn,
              Degree = degrees))
  
}
