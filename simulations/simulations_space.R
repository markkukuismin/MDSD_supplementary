
# 08.02.2024

library(igraph)
library(huge)
library(tidyverse)
library(hglasso)
library(space)

source("functions/inter_network.R")
source("functions/BA_network.R")
source("functions/space_BIC.R")

source("functions/hub_detection_hglasso.R")
source("functions/evaluation_metrics.R") # For hubs

# Only (n, p) %in% (100, 500) due to computational requirements
# of space

p = 500 # The number of nodes
n = 100 # The number of samples

set.seed(1)

# Graphical model construction,

G_type = "hglassotwo" # hglasso, hglassotwo, inter, scale-free, or star

if(G_type == "hglasso"){
  
  true_nmb_hubs = 10
  
  network = HubNetwork(p, 0.99, true_nmb_hubs, 0.1)
  
  true_hubs = sort(network$hubcol)
  
  Sigma = solve(network$Theta)
  
  A_true = as.matrix(network$Theta) 
  
}

if(G_type == "hglassotwo"){
  
  true_nmb_hubs = 10
  
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

G_true = igraph::graph.adjacency(A_true, mode = "undirected", diag = F)

node_names = as.character(1:p)

vertex_attr(G_true) = list(name = node_names)

d_true = igraph::degree(G_true)

# Allocate vectors, tables and matrices to collect simulation
# results,

M = 100 # The number of simulation rounds

temp = evaluation_metrics(1, 1:2, 1:2)

temp = names(temp)

Results_bic = Results_detection = matrix(0, M, length(temp))

temp[temp == "Precision"] = "Pre"

colnames(Results_bic) = 
  colnames(Results_detection) = temp

##

# Tuning parameters for space,

nlambda = 10
lam1 = seq(0.1, 2, length.out = nlambda)
alpha = 1
l1 = 1/sqrt(n)*qnorm(1-alpha/(2*p^2))
iter = 3
weight = 2

Results_space = matrix(0, nlambda, length(temp))

temp[temp == "Precision"] = "Pre"

colnames(Results_space) = temp

rm(temp)

##

space_path = list(wi = array(0, c(p, p, length(lam1))),
                  rholist = lam1)

BIC_sp = rep(0, nlambda)

# Simulation loop,

for(i in 1:M){
  
  # Simulate data,
  
  X = MASS::mvrnorm(n, rep(0, p), Sigma = Sigma)
  
  X = scale(X)
  
  # Compute space solution path,
  
  for(j in 1:nlambda){
    
    ddpcr::quiet(space_est <- space.joint(X, lam1=l1*n*lam1[j], lam2=0, weight=weight, iter=iter))
    
    space_path$wi[,,j] = space_est$ParCor
    
    A_space = space_est$ParCor
    
    A_space[A_space != 0] = 1
    
    diag(A_space) = 0
    
    space_hubs = colSums(A_space)
    
    space_hubs = order(space_hubs, decreasing = TRUE)[1:true_nmb_hubs]
    
    space_res = evaluation_metrics(est_hubs = space_hubs, 
                             true_hubs = true_hubs,
                             node_names = node_names)
    
    Results_space[j, ] = unlist(space_res)
    
    BIC_sp[j] = space_BIC(X = X, space_res = space_est)
    
  }
  
  ind = which.min(BIC_sp)
  
  ind = tail(ind, 1)
  
  detection_res = hub_detection_hglasso(hglasso_path = space_path,
                                        node_names = node_names)
  
  space_hubs_MDSD = evaluation_metrics(est_hubs = detection_res$hub_nodes_MDSD, 
                                         true_hubs = true_hubs,
                                         node_names = node_names)
  
  Results_bic[i, ] = Results_space[ind, ]
  Results_detection[i, ] = unlist(space_hubs_MDSD)
  
  cat("\r", i)
  
}

# Data cleaning before saving,

Results_bic = as.data.frame(Results_bic)
Results_detection = as.data.frame(Results_detection)

Results_bic$Criterion = "BIC"
Results_detection$Criterion = "MDSD"

Results_hub = rbind(Results_bic, Results_detection)

Results_temp = Results_hub %>%
  tidyr::pivot_longer(!Criterion, names_to = "Metric", values_to = "Value")

Results_temp %>%
  group_by(Criterion, Metric) %>%
  summarize(mean = mean(Value, na.rm = TRUE)) %>%
  dplyr::arrange(Metric) 

Results_temp = Results_temp %>%
  dplyr::filter(Metric %in% c("FDR", "TPR", "Pre", "FPR", "MCC"))

# Visualize results,

box_plot = ggplot(Results_temp, aes(x = Metric, y = Value)) +
  geom_boxplot(aes(fill = Criterion), position=position_dodge(.9)) +
  facet_wrap(~ Criterion) +
  theme(legend.position = "none")

box_plot

# Save results,

write.table(Results_hub, 
            paste0(
              "simulations/space_results/space_hub_results_", G_type, "_n=", n, "_p=", p, ".txt"
            ), 
            row.names = FALSE
)
