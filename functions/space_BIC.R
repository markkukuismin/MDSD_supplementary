
space_BIC = function(X = NULL, space_res = NULL){
  
  if(is.null(X)) stop("X is missing with no default.")
  if(is.null(space_res)) stop("space_res is missing with no default.")
  
  p = ncol(X)
  n = nrow(X)
  
  A = space_res$ParCor
  
  A[A != 0] = 1
  
  diag(A) = 0
  
  d = colSums(A)
  
  RSS = BIC = rep(0, p)
  
  for(i in 1:p){
    
    sqrt_sig = sqrt(space_res$sig.fit[-i]/space_res$sig.fit[i])
    
    D = sweep(X[, -i], 2, sqrt_sig, "*")
    
    D = sweep(D, 2, space_res$ParCor[i, -i] , "*")
    
    D = rowSums(D)
    
    RSS[i] = sum((X[, i] - D)^2)
    
    BIC[i] = n*log(RSS[i]) + log(n)*d[i]
    
  }
  
  BIC = sum(BIC)
  
  return(BIC)
  
}