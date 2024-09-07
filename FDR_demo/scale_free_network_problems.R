
# 05.09.2024

library(igraph)
library(huge)
library(tidyverse)
library(pracma)
library(hglasso)
library(glmnet)
library(doParallel)
library(hglasso)

source("functions/BA_network.R")

p = 500 # The number of nodes
n = 510 # The number of samples

nlambda = 50

set.seed(1)
    
Data = BA_network(p = p)
    
true_hubs = Data$hubcol
    
true_nmb_hubs = length(true_hubs)
    
Sigma = Data$Sigma
    
A_true = Data$A
    
A_true = as.matrix(A_true)
  
# glasso analysis
  
X = MASS::mvrnorm(n, rep(0, p), Sigma = Sigma)
  
X = scale(X)
  
res = huge(X, 
           method = "glasso",
           nlambda = nlambda,
           scr = TRUE,
           verbose = FALSE)

true_hubs

lambda = res$lambda

node_names = 1:p

Degree_lambda = matrix(0, nrow = p, ncol = nlambda)

colnames(Degree_lambda) = lambda
rownames(Degree_lambda) = node_names

for(j in 1:nlambda){
  
  A_temp = res$path[[j]]
  
  A_temp[A_temp != 0] = 1
  
  diag(A_temp) = 0
  
  d_temp = colSums(A_temp)
  
  names(d_temp) = node_names
  
  Degree_lambda[node_names, j] = d_temp
  
}

plot(lambda, 
     Degree_lambda[1, ], 
     pch = " ",
     ylim = c(min(Degree_lambda), max(Degree_lambda)))

for(i in 1:p){
  
  if(i %in% true_hubs){
    lines(lambda, Degree_lambda[i, ], lwd = 2, col = "red")
  }else{
    lines(lambda, Degree_lambda[i, ], lty = 2)
  }
  
}
