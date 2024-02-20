
# 30.01.2024

library(igraph)
library(huge)
library(tidyverse)
library(hglasso)

source("functions/inter_network.R")
source("functions/BA_network.R")

source("functions/fh_hub_screening.R")
source("functions/evaluation_metrics.R") # For hubs

# (n, p) %in% (100, 500) and (500, 1500)

p = 1500 # The number of nodes
n = 500 # The number of samples

set.seed(1)

G_type = "hglassotwo" # hglasso, hglassotwo, inter, scale-free, or star

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

G_true = igraph::graph_from_adjacency_matrix(A_true, mode = "undirected", diag = F)

node_names = as.character(1:p)

vertex_attr(G_true) = list(name = node_names)

d_true = igraph::degree(G_true)

M = 100 # The number of simulation rounds

Results_fh = Results_detection = matrix(0, M, 9)

temp = evaluation_metrics(1, 1:2, 1:2)

temp = names(temp)

temp[temp == "Precision"] = "Pre"

colnames(Results_fh) = temp

rm(temp)

# Threshold parameter for Firouzi & Hero hub screening,

delta = 20

for(i in 1:M){

  X = MASS::mvrnorm(n, rep(0, p), Sigma = Sigma)
  
  X = scale(X)
  
  fh_res = fh_hub_screening(X, delta = delta)
  
  fh_hubs = evaluation_metrics(est_hubs = fh_res$hub_nodes_fh, 
                                        true_hubs = true_hubs,
                                        node_names = node_names
                               )
  
  Results_fh[i, ] = unlist(fh_hubs)
  
  cat("\r", i)
  
}


Results_fh = as.data.frame(Results_fh)

Results_fh$Criterion = "FH"

Results_temp = Results_fh %>%
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

write.table(Results_fh, 
            paste0(
              "simulations/Firouzi_Hero_results/FH_hub_results_", G_type, "_n=", n, "_p=", p, ".txt"
            ), 
            row.names = FALSE
)
