
# 15.8.2024

# Is there a workaround for the non-monotonic solution
# path?

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

Data = huge.generator(n = n, 
                      d = p, 
                      g = true_nmb_hubs, 
                      graph = "hub")

Sigma = Data$sigma

A_true = Data$theta

A_true = as.matrix(A_true)

true_hubs = colSums(A_true)
true_hubs = order(-true_hubs)[1:true_nmb_hubs]

X = MASS::mvrnorm(n = n, mu = rep(0, p), Sigma)

nlambda = 50

res = huge(X, 
           nlambda = nlambda, 
           method = "glasso",
           scr = TRUE)

MDSD_res = hub_detection_huge(res, gamma = 3, plot = F)

MDSD = MDSD_res$sol_path_data

##

node_names = colnames(X)

if(is.null(node_names)) node_names = 1:p

lambda = res$lambda

nlambda = length(lambda)

Degree_lambda = matrix(0, nrow = p, ncol = nlambda)

colnames(Degree_lambda) = lambda
rownames(Degree_lambda) = node_names

for(j in 1:nlambda){
  
  d_temp = colSums(as.matrix(res$path[[j]]))
  
  names(d_temp) = node_names
  
  Degree_lambda[node_names, j] = d_temp
  
}

degrees = c(t(Degree_lambda))

lambdas = rep(lambda, times = p)

all_nodes = rep(node_names, each = ncol(Degree_lambda))

cumMDSD = matrix(0, p, nlambda)

colnames(cumMDSD) = lambda

for(j in 1:p){
  
  d = (Degree_lambda[rep(j, p - 1), ] - Degree_lambda[-j, ])^2
  
  cumMDSD[j, ] = cumsum(colSums(d))/((p-1)*(1:nlambda))
  
}

##

node_h = rep("No hub", p)
node_h[true_hubs] = "Hub"
node_h = as.factor(node_h)

D = data.frame(Degree = as.vector(t(Degree_lambda)),
               cumMDSD = as.vector(t(cumMDSD)),
               node = rep(node_names, each = nlambda),
               node_type = rep(node_h, each = nlambda),
               lambda = rep(lambda, times = p))

D$node = as.factor(D$node)

p1 = ggplot(D, 
            aes(x = lambda, y = cumMDSD, 
                group = node, color = node_type)
) +
  geom_line(aes(lty = node_type)) +
  theme(legend.position = "bottom", 
        plot.title = element_text(face = "bold")
  ) + 
  ylab("Cumulative MDSD") +
  xlab("Tuning parameter value")

##

p2 = ggplot(D, 
            aes(x = Degree, y = cumMDSD, 
                group = node, color = node_type)
) +
  geom_line(aes(lty = node_type)) +
  theme(legend.position = "bottom", 
        plot.title = element_text(face = "bold")
  ) + 
  ylab("Cumulative MDSD") +
  xlab("Degree") + 
  labs(title = "(A)")

###

# The solution path is not monotonic,

x = Degree_lambda[101, ]

all(diff(x) >= 0)

plot(lambda, x, type = "l")

##

# How does the burning phase change the curve?

d_skewness = apply(Degree_lambda, 2, e1071::skewness)

d_skewness[is.na(d_skewness)] = 0

burn_thr = 0.5

burn = which(d_skewness <= burn_thr)

b = length(burn)

cumMDSD_burn = matrix(0, p, nlambda - b)

colnames(cumMDSD_burn) = lambda[-burn]

lambda_burn = lambda[-burn]
nlambdab = length(lambda_burn)

for(j in 1:p){
  
  d = (Degree_lambda[rep(j, p - 1), -burn] - Degree_lambda[-j, -burn])^2
  
  cumMDSD_burn[j, ] = cumsum(colSums(d))/((p-1)*(1:(nlambda - b)))
  
}

##

node_h = rep("No hub", p)
node_h[true_hubs] = "Hub"
node_h = as.factor(node_h)

D_burn = data.frame(Degree = as.vector(t(Degree_lambda[,-burn])),
                    cumMDSD = as.vector(t(cumMDSD_burn)),
                    node = rep(node_names, each = nlambdab),
                    node_type = rep(node_h, each = nlambdab),
                    lambda = rep(lambda_burn, times = p))

D_burn$node = as.factor(D_burn$node)

p1_burn = ggplot(D_burn, 
            aes(x = lambda, y = cumMDSD, 
                group = node, color = node_type)
) +
  geom_line(aes(lty = node_type)) +
  theme(legend.position = "bottom", 
        plot.title = element_text(face = "bold")
  ) + 
  ylab("Cumulative MDSD") +
  xlab("Tuning parameter value") +
  labs(title = "(A)")

grid.arrange(p1, p1_burn, ncol = 2)

##

p2_burn = ggplot(D_burn, 
            aes(x = Degree, y = cumMDSD, 
                group = node, color = node_type)
) +
  geom_line(aes(lty = node_type)) +
  theme(legend.position = "bottom", 
        plot.title = element_text(face = "bold")
  ) + 
  ylab("Cumulative MDSD") +
  xlab("Degree") +
  labs(title = "(B)")

grid.arrange(p2, p2_burn, ncol = 2)

tiff("burn-in_demo/Fig1c.tif",
     units = "in",
     res = 400,
     width = 9,
     height = 4)

grid.arrange(p2, p2_burn, ncol = 2)

dev.off()

pdf("burn-in_demo/Fig1c.pdf",
    width = 9,
    height = 4)

grid.arrange(p2, p2_burn, ncol = 2)

dev.off()

png("burn-in_demo/Fig1c.png",
    units = "in",
    res = 400,
    width = 9,
    height = 4)

grid.arrange(p2, p2_burn, ncol = 2)

dev.off()
