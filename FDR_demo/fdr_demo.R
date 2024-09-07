
# 05.09.2024

library(igraph)
library(huge)
library(tidyverse)
library(pracma)
library(hglasso)
library(glmnet)
library(doParallel)
library(hglasso)
library(gridExtra)
library(patchwork)

source("functions/evaluation_metrics.R") # For hubs
source("functions/inter_network.R")
source("functions/BA_network.R")

# How the hubs detected with MDSD change when gamma
# changes?

gamma = seq(1.5, 4, length.out = 10)
gamma = c(gamma, 3)
gamma = sort(gamma)
burn_thr = 0.5

# n in (100, 300, 510)

p = 500 # The number of nodes
n = 100 # The number of samples

nlambda = 50

MDSD_data = data.frame()

set.seed(1)

G_models = c("hglasso", "hglassotwo", "Inter", 
             "Scale-free", "Star")
  
for(G_type in G_models){
  
  if(G_type == "hglasso"){
    
    true_nmb_hubs = 30
    
    network = HubNetwork(p, 0.99, true_nmb_hubs, 0.1)
    
    true_hubs = sort(network$hubcol)
    
    Sigma = solve(network$Theta)
    
    A_true = as.matrix(network$Theta) 
    
  }
  
  if(G_type == "hglassotwo"){
    
    true_nmb_hubs = 30
    
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
  
  if(G_type == "Star"){
    
    Data = huge.generator(n = n, d = p, graph = "hub")
    
    true_nmb_hubs = ifelse(p > 40, ceiling(p/20), 2) 
    
    Sigma = Data$sigma
    
    A_true = Data$theta
    
    A_true = as.matrix(A_true)
    
    true_hubs = colSums(A_true)
    
    true_hubs = order(-true_hubs)[1:true_nmb_hubs]
    
  }
  
  if(G_type == "Inter"){
    
    Data = inter_network(p = p)
    
    true_hubs = Data$hubs
    
    true_nmb_hubs = length(true_hubs)
    
    Sigma = Data$Sigma
    
    A_true = Data$A
    
    A_true = as.matrix(A_true)
    
  }
  
  if(G_type == "Scale-free"){
    
    Data = BA_network(p = p)
    
    true_hubs = Data$hubcol
    
    true_nmb_hubs = length(true_hubs)
    
    Sigma = Data$Sigma
    
    A_true = Data$A
    
    A_true = as.matrix(A_true)
    
  }
  
  # glasso analysis
  
  X = MASS::mvrnorm(n, rep(0, p), Sigma = Sigma)
  
  X = scale(X)
  
  res = huge(X, 
             method = "glasso",
             nlambda = nlambda,
             scr = TRUE,
             verbose = FALSE)
  
  # Compute MDSD,
  
  node_names = 1:p
  
  glasso_path = list(wi = array(0, c(p, p, nlambda)),
                     rholist = res$lambda)
  
  lambda = glasso_path$rholist
  
  Degree_lambda = matrix(0, nrow = p, ncol = nlambda)
  
  colnames(Degree_lambda) = lambda
  rownames(Degree_lambda) = node_names
  
  for(j in 1:nlambda){
    
    glasso_path$wi[,,j] = res$path[[j]]
    
    A_temp = glasso_path$wi[,, j]
    
    A_temp[A_temp != 0] = 1
    
    diag(A_temp) = 0
    
    d_temp = colSums(A_temp)
    
    names(d_temp) = node_names
    
    Degree_lambda[node_names, j] = d_temp
    
  }
  
  rm(A_temp)
  
  MDSD = rep(0, p)
  
  for(j in 1:p){
    
    MDSD[j] = mean((Degree_lambda[rep(j, p - 1), ] - Degree_lambda[-j, ])^2)
    
  }
  
  d_skewness = apply(Degree_lambda, 2, e1071::skewness)
  
  d_skewness[is.na(d_skewness)] = 0
  
  burn = which(d_skewness <= burn_thr)
  
  if(!rlang::is_empty(burn)){
    
    MDSD_burn = rep(0, p)
    
    for(j in 1:p){
      
      MDSD_burn[j] = mean((Degree_lambda[rep(j, p - 1), -burn] - Degree_lambda[-j, -burn])^2)
      
    }
    
  }else{
    
    MDSD_burn = MDSD
    
  }
  
  
  names(MDSD) = rownames(Degree_lambda)
  names(MDSD_burn) = rownames(Degree_lambda)
  
  hub_nodes_MDSD = list()
  hub_nodes_MDSD_burn = list()
  
  FDR = FDR_burn = list()
  
  for(i in 1:length(gamma)){
    
    hub_nodes_MDSD[[i]] = node_names[MDSD > gamma[i]*mean(MDSD)]
    
    hub_nodes_MDSD_burn[[i]] = node_names[MDSD_burn > gamma[i]*mean(MDSD_burn)]
    
    FDR[[i]] = evaluation_metrics(hub_nodes_MDSD[[i]],
                                  true_hubs,
                                  node_names = node_names)$FDR
    
    FDR_burn[[i]] = evaluation_metrics(hub_nodes_MDSD_burn[[i]],
                                       true_hubs,
                                       node_names = node_names)$FDR
    
  }
  
  FDR = unlist(FDR)
  
  ##
  
  FDR_burn = unlist(FDR_burn)
  
  l = length(gamma)
  
  temp_MDSD_data = data.frame(gamma = gamma,
                              FDR = FDR,
                              G = rep(G_type, l),
                              Burn = rep("No", l))
  
  temp_MDSD_burn_data = data.frame(gamma = gamma,
                                   FDR = FDR_burn,
                                   G = rep(G_type, l),
                                   Burn = rep("Yes", l))
  
  MDSD_data = dplyr::bind_rows(MDSD_data, 
                               temp_MDSD_data,
                               temp_MDSD_burn_data)
  
  cat("\n", G_type)
  
}

MDSD_data$G = as.factor(MDSD_data$G)
MDSD_data$Burn = as.factor(MDSD_data$Burn)

write.table(MDSD_data, 
            "FDR_demo/FDR_results.txt",
            row.names = FALSE)

MDSD_data = MDSD_data %>%
  mutate(G = case_when(G == "hglasso" ~ "Hub network", 
                       G == "hglassotwo" ~ "Hub network, two components",
                       TRUE ~ G))

MDSD_data = MDSD_data %>%
  rename(Graph = G)

p1 = ggplot(data = subset(MDSD_data, Burn %in% "No"), 
            aes(x = gamma, 
                y = FDR, 
                color = Graph, 
                group = Graph)) +
  geom_point(aes(color = Graph, 
                 shape = Graph)) +
  geom_line() +
  ylim(c(0, 1)) +
  xlab(expression(gamma)) +
  theme(
    axis.title.x = element_text(size = 15)
  ) + 
  labs(title = "(A)")

p2 = ggplot(data = subset(MDSD_data, Burn %in% "Yes"),
       aes(x = gamma, 
           y = FDR, 
           color = Graph, 
           group = Graph)) +
  geom_point(aes(color = Graph, 
                 shape = Graph)) +
  geom_line() +
  ylim(c(0, 1)) +
  xlab(expression(gamma)) +
  theme(
    axis.title.x = element_text(size = 15)
    ) + 
  labs(title = "(B)")

final_plot = p1 + 
  p2 & theme(legend.position = "bottom",
                             legend.title = element_blank())

final_plot = final_plot + patchwork::plot_layout(guides = "collect")

final_plot

png(paste0("fdr_demo/figures/fdr_vs_gamma_", n, ".png"),
    units = "in",
    res = 300,
    width = 9,
    height = 4)

final_plot

dev.off()

tiff(paste0("fdr_demo/figures/fdr_vs_gamma_", n, ".tif"),
     units = "in",
     res = 300,
     width = 9,
     height = 4)

final_plot

dev.off()

pdf(paste0("fdr_demo/figures/fdr_vs_gamma_", n, ".pdf"),
    width = 9,
    height = 4)

final_plot

dev.off()
