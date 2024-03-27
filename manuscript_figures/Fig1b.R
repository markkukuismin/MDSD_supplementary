
# 21.03.2024

library(hglasso)
library(huge)
library(ggplot2)
library(dplyr)
library(gridExtra)

source("functions/hub_detection_huge.R")

p = 500
n = 50

set.seed(1)

true_nmb_hubs = 5

Data = huge.generator(n = n, d = p, graph = "hub", g = true_nmb_hubs)

Sigma = Data$sigma

A_true = Data$theta

A_true = as.matrix(A_true)

true_hubs = colSums(A_true)

true_hubs = order(-true_hubs)[1:true_nmb_hubs]

X = MASS::mvrnorm(n = n, mu = rep(0, p), Sigma)

nlambda = 50

res = huge(X, nlambda = nlambda, method = "glasso")

burn_thr = 0.5

hub_res = hub_detection_huge(huge_est = res, burn_thr = burn_thr, plot = FALSE)

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

MDSD_p = mean(sol_path_data$MDSD)

p2 = p2 + geom_hline(yintercept = 3*MDSD_p)

p2

## When degree distribution skewness is taken into account,

node_names = as.character(1:p)

Degree_lambda = matrix(0, nrow = p, ncol = nlambda)

colnames(Degree_lambda) = res$lambda
rownames(Degree_lambda) = node_names

for(j in 1:nlambda){
  
  d_temp = colSums(as.matrix(res$path[[j]]))
  
  names(d_temp) = node_names
  
  Degree_lambda[node_names, j] = d_temp
  
}

d_skewness = apply(Degree_lambda, 2, e1071::skewness)

d_skewness[is.na(d_skewness)] = 0

burn = which(d_skewness < burn_thr)

thr = res$lambda[burn[2]]

p3 = ggplot(degree_data, 
            aes(x = lambda, y = Degree, group = node, color = node_type)
) +
  geom_line(aes(lty = node_type)) +
  theme(legend.position = "bottom", 
        plot.title = element_text(face = "bold")
  ) + 
  ylab("Node degree") +
  xlab("Tuning parameter value") +
  labs(title = "(C)") +
  geom_vline(xintercept = thr, linetype = "dashed")

p3

p4 = ggplot(sol_path_data, aes(x = node, y = MDSD_burn)) +
  geom_segment(
    aes(x = node, xend = node, y = 0, yend = MDSD_burn, 
        colour = node_type, lty = node_type)
  ) +
  theme(axis.ticks.x=element_blank(),
        legend.position = "bottom",
        plot.title = element_text(face = "bold")
  ) +
  ylab("MDSD") +
  xlab("Node label") +
  labs(title = "(D)")

p4

MDSD_burn_p = mean(sol_path_data$MDSD_burn)

p4 = p4 + geom_hline(yintercept = 3*MDSD_burn_p)

p4

##

grid.arrange(p1, p2, p3, p4, ncol = 2)

tiff("manuscript_figures/Fig1b.tif",
     units = "in",
     res = 300,
     width = 7,
     height = 6)

grid.arrange(p1, p2, p3, p4, ncol = 2)

dev.off()

pdf("manuscript_figures/Fig1b.pdf",
    width = 7,
    height = 6)

grid.arrange(p1, p2, p3, p4, ncol = 2)

dev.off()

png("manuscript_figures/Fig1b.png",
    units = "in",
    res = 300,
    width = 7,
    height = 6)

grid.arrange(p1, p2, p3, p4, ncol = 2)

dev.off()
