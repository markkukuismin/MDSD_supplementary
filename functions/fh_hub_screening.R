
fh_hub_screening = function(X, rho = NULL, delta = NULL, alpha = 0.01){

  if(is.null(delta)) stop("delta is missing with no default")
  
  n = nrow(X)
  
  p = ncol(X)
  
  node_names = colnames(X)
  
  if(is.null(node_names)) node_names = 1:p
  
  if(is.null(rho)){
    
    cn = 2^(-2/(n - 4))
    
    rho = sqrt(1 - cn)
    
  }
  
  S = abs(cor(X))
  
  S = S - diag(1, p)
  
  S[S < rho] = 0
  
  rho_delta = apply(S, 2, function(x) sort(x, decreasing = TRUE)[delta])
  
  Lambda = (p - 1)*zipfR::Rbeta(x = 1 - rho_delta^2, a = (n - 2)/2, b = 1/2, lower = TRUE, log = FALSE)

  p_values = rep(0, p)
  
  for(i in 1:p){
    
    p_values[i] = 1 - ppois(delta - 1, lambda = Lambda[i])
    
  }
  
  hub_nodes_fh = node_names[p_values < alpha]
  
  return(list(hub_nodes_fh = hub_nodes_fh, p_values = p_values))
  
}