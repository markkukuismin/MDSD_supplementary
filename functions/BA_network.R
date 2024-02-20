
BA_network = function(p = 500){
  
  hubcol = c(1, 2)
  
  while(length(hubcol) < 3){
    
    G = sample_pa(n = p, directed = FALSE, power = 1.5) 
    
    d = degree(G)
    
    hubcol = which(d >= (p - 1)*0.05)
    
  }
  
  A = igraph::as_adjacency_matrix(G, type = "both")
  
  A = as.matrix(A)
  
  a = rbinom(p^2, 1, 0.5)
  
  u = a*runif(p^2, min = -0.75, max = -0.25) + (1 - a)*runif(p^2, min = 0.25, max = 0.75)
  
  E = matrix(u, p, p)
  
  E = (E + t(E))/2
  
  E = A*E
  
  rho = min(eigen(E, only.values = TRUE)$values)
  
  Theta = E + (0.1 - rho)*diag(p)
  
  Sigma = solve(Theta)
  
  return(list(Theta = Theta, Sigma = Sigma, A = A, hubcol = hubcol))
  
}