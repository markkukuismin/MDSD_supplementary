
# Generates a graph with inter- and intramodular hub nodes

inter_network = function(p, v = 0.3, u = 0.1){
  
  p = p - 1
  
  if(p > 40) g = ceiling(p/20)
  if(p <= 40) g = 2
  
  ds = c(rep((p - p %% g)/g, g), p %% g)
  
  ds = ds[ds != 0]
  
  if(ds[length(ds)] < ds[1]){
    
    ds[length(ds) - 1] = ds[1] + ds[length(ds)]
    
    ds = head(ds, length(ds) - 1)
    
  }
  
  lst = list()
  
  i = 1
  
  for(d in ds){
    
    A = matrix(0, nrow = d, ncol = d)
    
    A[1, ] = A[, 1] = 1
    
    lst[[i]] = A
    
    i = i + 1
    
  }
  
  A = Matrix::bdiag(lst)
  
  diag(A) = 0 
  
  A = rbind(0, A)
  
  A = cbind(0, A)
  
  ind = which(colSums(as.matrix(A)) > 1)
  
  A[ind, 1] = A[1, ind] = 1
  
  G = igraph::graph_from_adjacency_matrix(A, mode = "undirected")
  
  hubs = which(colSums(as.matrix(A)) > 1)
  
  Theta = A*v
  
  min_evalue = min(eigen(Theta, only.values = TRUE)$values)
  
  diag(Theta) = abs(min_evalue) + 0.1 + u
  Sigma = solve(Theta)
  Sigma = cov2cor(Sigma)
  Theta = solve(Theta)
  
  ls = list(A = A, 
            G = G, 
            hubs = hubs, 
            Sigma = Sigma, 
            Theta = Theta
            )
  
  return(ls)
  
}