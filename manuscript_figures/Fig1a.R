
# 24.01.2024

library(hglasso)
library(huge)
library(ggplot2)
library(dplyr)
library(gridExtra)

source("functions/hub_detection_huge.R")
source("functions/BA_network.R")

p = 500
n = 110

set.seed(1)

true_nmb_hubs = 5

Data = huge.generator(n = n, d = p, g = true_nmb_hubs, graph = "hub")

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

##

ric_res = huge.select(res, criterion = "ric") 

ric_A = ric_res$refit

ric_d = colSums(ric_A)

summary(ric_d)

ric_df = data.frame(node = 1:p, degree = ric_d, node_type = node_h)

p2 = ggplot(ric_df, aes(x = node, y = degree)) +
  geom_segment(
    aes(x = node, xend = node, y = 0, yend = degree, 
        colour = node_type, lty = node_type)
  ) +
  theme(axis.ticks.x=element_blank(),
        legend.position = "bottom",
        plot.title = element_text(face = "bold")
  ) +
  ylab("Degree") +
  xlab("Node label") +
  labs(title = "(B)")

p2

##

sol_path_data = hub_res$sol_path_data

sol_path_data$node_type = node_h

p3 = ggplot(sol_path_data, aes(x = node, y = MDSD)) +
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
  labs(title = "(C)")

p3

MDSD_p = mean(sol_path_data$MDSD)

p3 = p3 + geom_hline(yintercept = 3*MDSD_p)

p3

grid.arrange(p1, p2, p3, ncol = 3)

tiff("manuscript_figures/Fig1.tif",
     units = "in",
     res = 300,
     width = 9,
     height = 4)

grid.arrange(p1, p2, p3, ncol = 3)

dev.off()

pdf("manuscript_figures/Fig1.pdf",
    width = 9,
    height = 4)

grid.arrange(p1, p2, p3, ncol = 3)

dev.off()

png("manuscript_figures/Fig1.png",
    units = "in",
    res = 300,
    width = 9,
    height = 4)

grid.arrange(p1, p2, p3, ncol = 3)

dev.off()
