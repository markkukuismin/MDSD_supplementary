
# 01.03.2024

library(hglasso)
library(huge)
library(ggplot2)
library(dplyr)
library(gridExtra)

source("functions/hub_detection_huge.R")

p = 500
n = 200

set.seed(1)

true_nmb_hubs = 10

Data = huge.generator(n = n, d = p, 
                      g = true_nmb_hubs, graph = "hub")

Sigma = Data$sigma

A_true = Data$theta

A_true = as.matrix(A_true)

true_hubs = colSums(A_true)

true_hubs = order(-true_hubs)[1:true_nmb_hubs]

X = MASS::mvrnorm(n = n, mu = rep(0, p), Sigma)

nlambda = 50

res = huge(X, nlambda = nlambda, method = "glasso")

hub_res = hub_detection_huge(huge_est = res, plot = FALSE)

degree_data = hub_res$raw_data

node_h = rep("No hub", p)
node_h[true_hubs] = "Hub"
node_h = as.factor(node_h)

degree_data$node_type = rep(node_h, each = nlambda)

p1 = ggplot(degree_data, 
            aes(x = lambda, y = Degree, group = node, color = node_type)
) +
  geom_line(aes(lty = node_type)) +
  theme(legend.position = "bottom", 
        plot.title = element_text(face = "bold")
  ) + 
  ylab("Node degree") +
  xlab("Tuning parameter value") +
  labs(title = "(A)")

p1

sol_path_data = hub_res$sol_path_data

sol_path_data$node_type = node_h

p2 = ggplot(sol_path_data, aes(x = node, y = MDSD)) +
  geom_segment(
    aes(x = node, xend = node, y = 0, yend = MDSD, 
        colour = node_type, lty = node_type)
  ) +
  theme(axis.ticks.x=element_blank(),
        legend.position = "bottom",
        plot.title = element_text(face = "bold")
  ) +
  ylab("MDSD") +
  xlab("Node label") +
  labs(title = "(B)")

p2

###

# Add vertical line to plot p1 indicating models where the 
# degree distribution is left skewed,

Degree_lambda = matrix(0, nrow = p, ncol = nlambda)

node_names = 1:p

colnames(Degree_lambda) = res$lambda
rownames(Degree_lambda) = node_names

for(j in 1:nlambda){
  
  d_temp = colSums(as.matrix(res$path[[j]]))
  
  names(d_temp) = node_names
  
  Degree_lambda[node_names, j] = d_temp
  
}

d_skewness = apply(Degree_lambda, 2, e1071::skewness)

d_skewness[is.na(d_skewness)] = 0

round(d_skewness, 3)

v = res$lambda[which(d_skewness < 0)[1]]

p1 = p1 + geom_vline(xintercept = v, linetype = "dashed")

###

MDSD_p = mean(sol_path_data$MDSD)

p2 = p2 + geom_hline(yintercept = 3*MDSD_p)

grid.arrange(p1, p2, ncol = 2)
