
# 09.02.2024

library(igraph)
library(huge)
library(tidyverse)
library(hglasso)

source("functions/inter_network.R")
source("functions/BA_network.R")

p = 500 # The number of nodes
n = 100 # The number of samples

set.seed(1)

G_type = "star" # hglasso, hglassotwo, inter, scale-free, or star

if(G_type == "hglasso"){
  
  if(p < 1500) true_nmb_hubs = 10
  if(p >= 1500) true_nmb_hubs = 30
  
  network = HubNetwork(p, 0.99, true_nmb_hubs, 0.1)
  
  true_hubs = sort(network$hubcol)
  
  Sigma = solve(network$Theta)
  
  A_true = as.matrix(network$Theta) 
  
}

if(G_type == "hglassotwo"){
  
  if(p < 1500) true_nmb_hubs = 10
  if(p >= 1500) true_nmb_hubs = 30
  
  network1 = HubNetwork(p/2, 0.99, true_nmb_hubs/2, 0.1)
  network2 = HubNetwork(p/2, 0.99, true_nmb_hubs/2, 0.1)
  
  Theta = Matrix::bdiag(network1$Theta, network2$Theta)
  
  hubcol = c(network1$hubcol, network2$hubcol + p/2)
  
  network = list(Theta = Theta, hubcol = hubcol)
  
  rm(hubcol)
  rm(Theta)
  
  true_hubs = sort(network$hubcol)
  
  Sigma = solve(network$Theta)
  
  A_true = as.matrix(network$Theta) 
  
}

if(G_type == "star"){
  
  Data = huge.generator(n = n, d = p, graph = "hub")
  
  true_nmb_hubs = ifelse(p > 40, ceiling(p/20), 2) 
  
  Sigma = Data$sigma
  
  A_true = Data$theta
  
  A_true = as.matrix(A_true)
  
  true_hubs = colSums(A_true)
  
  true_hubs = order(-true_hubs)[1:true_nmb_hubs]
  
}

if(G_type == "inter"){
  
  Data = inter_network(p = p)
  
  true_hubs = Data$hubs
  
  true_nmb_hubs = length(true_hubs)
  
  Sigma = Data$Sigma
  
  A_true = Data$A
  
  A_true = as.matrix(A_true)
  
}

if(G_type == "scale-free"){
  
  Data = BA_network(p = p)
  
  true_hubs = Data$hubcol
  
  true_nmb_hubs = length(true_hubs)
  
  Sigma = Data$Sigma
  
  A_true = Data$A
  
  A_true = as.matrix(A_true)
  
}

A_true[A_true != 0] = 1

diag(A_true) = 0

G = igraph::graph_from_adjacency_matrix(A_true, mode = "undirected", diag = F)

node_names = as.character(1:p)

vertex_attr(G) = list(name = node_names)

node_shape = rep("circle", p)

node_color = rep("gray", p)

node_color[true_hubs] = "black"

node_size = rep(2, p)

node_size[true_hubs] = 3

plot(G, vertex.shape = node_shape, vertex.color = node_color,
     vertex.size = node_size, vertex.label = NA)

fpath = paste0("manuscript_figures/Sup_Fig_", G_type, ".png")

png(fpath,
    units = "in",
    res = 300,
    width = 6,
    height = 6)

plot(G, vertex.shape = node_shape, vertex.color = node_color,
     vertex.size = node_size, vertex.label = NA)

dev.off()
