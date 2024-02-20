
# This can be described as a lossy screening method
cor_screening = function(X, nlambda = NULL, lambda.min.ratio = NULL){
  
  if(is.null(nlambda)) stop("nlambda is missing with no default")
  
  X = scale(X)
  
  S = cor(X)
  
  p = ncol(S)
  
  if(is.null(lambda.min.ratio)) lambda.min.ratio = 0.05
  
  cor_path = list(wi = array(0, dim = c(p, p, nlambda)), rholist = rep(0, nlambda))
  
  lambda.max = max(max(S - diag(p)), -min(S - diag(p)))
  lambda.min = lambda.min.ratio * lambda.max
  lambda = exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
  S = abs(S)
  for(i in 1:nlambda){
    cor_path$wi[,,i][S > lambda[i]] = 1
  }
  
  cor_path$rholist = lambda
  
  return(cor_path)
  
}
