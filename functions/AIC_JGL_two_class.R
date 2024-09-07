
AIC_JGL_two_class = function(data, jgl_res){
  
  n1 = dim(data[[1]])[1]
  n2 = dim(data[[2]])[1]
  
  S1 = cov(data[[1]])
  S2 = cov(data[[2]])
  
  Theta1 = jgl_res$theta[[1]]
  Theta2 = jgl_res$theta[[2]]
  
  E1 = sum(Theta1 != 0)
  E2 = sum(Theta1 != 0)
  
  a1 = n1*sum(diag(S1%*%Theta1)) - n1*determinant(Theta1, logarithm = TRUE)$modulus + 2*E1
  a2 = n2*sum(diag(S2%*%Theta2)) - n2*determinant(Theta2, logarithm = TRUE)$modulus + 2*E2
  
  c(a1 + a2)
  
}