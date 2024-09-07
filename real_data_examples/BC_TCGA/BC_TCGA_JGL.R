
# 19.08.2024

# Examine Breast cancer tissues vs normal tissues

# Data sets are publicly available here,

# https://data.mendeley.com/datasets/v3cc2p38hb/1

library(tidyverse)
library(DESeq2)
library(genefilter)

library(JGL)

source("functions/hub_detection_hglasso.R")
source("functions/AIC_JGL_two_class.R")

normal_data = read.csv("real_data_examples/BC_TCGA/BC-TCGA-Normal.txt", 
                       header = TRUE, sep = "\t", na.strings = "null")

tumor_data = read.csv("real_data_examples/BC_TCGA/BC-TCGA-Tumor.txt", 
                      header = TRUE, sep = "\t", na.strings = "null")

all.equal(tumor_data[, 1], normal_data[, 1])

dim(tumor_data)
dim(normal_data)

gene_names = tumor_data[, 1]

tumor_data = tumor_data[, -1]
normal_data = normal_data[, -1]

all(apply(tumor_data, 2, is.numeric))
all(apply(normal_data, 2, is.numeric))

# Remove genes with over 30% missing values,

tumor_mis = apply(tumor_data, 1, function(x) mean(is.na(x)))

summary(100*tumor_mis)

normal_mis = apply(normal_data, 1, function(x) mean(is.na(x)))

summary(100*normal_mis)

# No need to delete anything

# Impute missing values. Because the proportion of 
# missing values is less than 0.1% overall we use simple
# mean imputation

set.seed(1)

100*mean(is.na(normal_data))
100*mean(is.na(tumor_data))

f = function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x)

tumor_data = apply(tumor_data, 1, f)

normal_data = apply(normal_data, 1, f)

dim(tumor_data)
dim(normal_data)

##

qqnorm(tumor_data[1, ])
qqnorm(normal_data[1, ])

y = seq(-5, 5, length.out = 1000)

par(mfrow = c(1, 2))
hist(tumor_data[11, ], probability = TRUE, main = "tumor")
lines(y, dnorm(y), lwd = 2)

hist(normal_data[11, ], probability = TRUE, main = "normal")
lines(y, dnorm(y), lwd = 2)

###

# Select gene probes with a large variance to reduce the number
# of genes. In particular, select 1000 genes with the highest
# variance,

temp_data = rbind(normal_data, tumor_data)

dim(temp_data)

colnames(temp_data) = gene_names

top_vars = apply(temp_data, 2, function(x) var(x, na.rm = TRUE))
top_vars = sort(top_vars, decreasing = TRUE)[1:1000]
top_vars = names(top_vars)

colnames(normal_data) = gene_names
colnames(tumor_data) = gene_names

normal_data = normal_data[ , colnames(normal_data) %in%
                             top_vars]

tumor_data = tumor_data[ , colnames(tumor_data) %in%
                           top_vars]

dim(normal_data)
dim(tumor_data)

all.equal(colnames(normal_data), colnames(tumor_data))

# The fused network analysis,

data = list(normal_data, tumor_data)

p = ncol(normal_data)

s_temp = cor(normal_data)

s_temp = s_temp[upper.tri(s_temp)]

round(quantile(abs(s_temp)), 3)

nlambda = 10

lambda1 = 0.9 # the graphical lasso penalty
lambda2s = seq(0.1, 0.5, length.out = nlambda) # fused lasso penalty

fglassopath = list(normal_theta = array(0, c(p, p, nlambda)),
                   tumor_theta = array(0, c(p, p, nlambda)))

AIC_values = rep(0, nlambda)
                     
k = 1

start_time = proc.time()

for(lambda2 in lambda2s){
  
  fgl_res = JGL(data, 
                penalty = "fused",
                lambda1 = lambda1, 
                lambda2 = lambda2,
                weights = "equal",
                return.whole.theta = TRUE)
  
  AIC_values[k] = AIC_JGL_two_class(data, fgl_res)
  
  fglassopath$normal_theta[,,k] = fgl_res$theta[[1]]
  fglassopath$tumor_theta[,,k] = fgl_res$theta[[2]]
  
  cat("\r", k)
  
  k = k + 1
  
}

stop_time = proc.time()

#save.image(file = "real_data_examples/BC_TCGA/JGLexample.RData")
#load(file = "real_data_examples/BC_TCGA/JGLexample.RData")    

(stop_time - start_time)/(60*60)

plot(lambda2s, AIC_values)

ind = which.min(AIC_values)

ind

# Normal data hubs,

node_names = colnames(normal_data)

fglassopath_normal = list(wi = fglassopath$normal_theta,
                          rholist = 1:nlambda)

detection_res_normal = hub_detection_hglasso(hglasso_path = fglassopath_normal,
                                             node_names = node_names,
                                             gamma = 3)

detection_res_normal$hub_nodes_MDSD

detection_res_normal$hub_nodes_MDSD_burn

# Tumor data hubs,

fglassopath_tumor = list(wi = fglassopath$tumor_theta,
                         rholist = 1:nlambda)

detection_res_tumor = hub_detection_hglasso(hglasso_path = fglassopath_tumor,
                                            node_names = node_names,
                                            gamma = 3)

detection_res_tumor$hub_nodes_MDSD

detection_res_tumor$hub_nodes_MDSD_burn

# Solution path degree data,

degrees_normal = c(t(detection_res_normal$Degree))
degrees_tumor = c(t(detection_res_tumor$Degree))

Degree_data_normal = detection_res_normal$Degree
Degree_data_tumor = detection_res_tumor$Degree

all_nodes = rep(node_names, each = ncol(Degree_data_normal))

top_nodes_normal = detection_res_normal$hub_nodes_MDSD
top_nodes_tumor = detection_res_tumor$hub_nodes_MDSD

node_col_normal = node_col_tumor = 
  rep("#CCCCCC", length(all_nodes))

node_col_normal[all_nodes %in% top_nodes_normal] = "#F8766D"
node_col_tumor[all_nodes %in% top_nodes_tumor] = "#F8766D"

node_col_normal = as.factor(node_col_normal)
node_col_tumor = as.factor(node_col_tumor)

node_col_normal = relevel(node_col_normal, ref = "#F8766D")
node_col_tumor = relevel(node_col_tumor, ref = "#F8766D")

lty_normal = ifelse(node_col_normal == "#F8766D", "solid", "dashed")
lty_tumor = ifelse(node_col_tumor == "#F8766D", "solid", "dashed")

lty_normal = as.factor(lty_normal)
lty_tumor = as.factor(lty_tumor)

lty_normal = relevel(lty_normal, ref = "solid")
lty_tumor = relevel(lty_tumor, ref = "solid")

Degree_data_normal = data.frame(node = all_nodes,
                                Degree = degrees_normal,
                                lambda = lambda2s,
                                node_col = node_col_normal,
                                lty = lty_normal
)

Degree_data_tumor = data.frame(node = all_nodes,
                               Degree = degrees_tumor,
                               lambda = lambda2s,
                               node_col = node_col_tumor,
                               lty = lty_tumor
)

# Normal data: Solution path degree plot,

p_degree_normal = ggplot(Degree_data_normal, 
                         aes(x = lambda, y = Degree, group = node)
) +
  geom_line(aes(colour = node_col, linetype = lty)) +
  scale_colour_identity() +
  ylab("Node degree") +
  xlab("Tuning parameter value")

p_degree_normal = p_degree_normal + 
  labs(title = "(A)") + 
  theme_classic() +
  theme(legend.position = "none")

# Normal data: MDSD 

p = ncol(normal_data)

node_col_normal = node_col_tumor = rep("#CCCCCC", p)

names(node_col_normal) = names(node_col_tumor) = 
  node_names

node_col_normal[top_nodes_normal] = "#F8766D"
node_col_tumor[top_nodes_tumor] = "#F8766D"

node_col_normal = as.factor(node_col_normal)
node_col_tumor = as.factor(node_col_tumor)

node_col_normal = relevel(node_col_normal, ref = "#F8766D")
node_col_tumor = relevel(node_col_tumor, ref = "#F8766D")

lty_normal = ifelse(node_col_normal == "#F8766D", "solid", "dashed")
lty_tumor = ifelse(node_col_tumor == "#F8766D", "solid", "dashed")

lty_normal = as.factor(lty_normal)
lty_tumor = as.factor(lty_tumor)

lty_normal = relevel(lty_normal, ref = "solid")
lty_tumor = relevel(lty_tumor, ref = "solid")

MDSD_data_normal = data.frame(node = node_names,
                              MDSD = detection_res_normal$MDSD,
                              node_col = node_col_normal,
                              lty = lty_normal
)

MDSD_data_tumor = data.frame(node = node_names,
                              MDSD = detection_res_tumor$MDSD,
                              node_col = node_col_tumor,
                              lty = lty_tumor
)

p_MDSD_normal = ggplot(MDSD_data_normal, 
                       aes(x = node, y = MDSD)) +
  geom_segment(
    aes(x = node, xend = node, y = 0, yend = MDSD, 
        colour = node_col, linetype = lty)
  ) +
  scale_x_discrete(breaks = top_nodes_normal,
                   labels = top_nodes_normal
  ) +
  scale_colour_identity() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 7)
  ) +
  ylab("MDSD") +
  xlab("Gene name")

p_MDSD_normal = p_MDSD_normal + 
  labs(title = "(B)") + 
  theme_classic() +
  theme(legend.position = "none")

p_MDSD_normal = p_MDSD_normal + 
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank()
  )

gridExtra::grid.arrange(p_degree_normal, 
                        p_MDSD_normal, ncol = 2)

# tumor data: Solution path degree plot,

p_degree_tumor = ggplot(Degree_data_tumor, 
                         aes(x = lambda, y = Degree, group = node)
) +
  geom_line(aes(colour = node_col, linetype = lty)) +
  scale_colour_identity() +
  ylab("Node degree") +
  xlab("Tuning parameter value")

p_degree_tumor = p_degree_tumor + 
  labs(title = "(C)") + 
  theme_classic() +
  theme(legend.position = "none")

# tumor data: MDSD 

p_MDSD_tumor = ggplot(MDSD_data_tumor, 
                       aes(x = node, y = MDSD)) +
  geom_segment(
    aes(x = node, xend = node, y = 0, yend = MDSD, 
        colour = node_col, linetype = lty)
  ) +
  scale_x_discrete(breaks = top_nodes_normal,
                   labels = top_nodes_normal
  ) +
  scale_colour_identity() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 7)
  ) +
  ylab("MDSD") +
  xlab("Gene name")

p_MDSD_tumor = p_MDSD_tumor + 
  labs(title = "(D)") + 
  theme_classic() +
  theme(legend.position = "none")

p_MDSD_tumor = p_MDSD_tumor + 
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank()
  )

gridExtra::grid.arrange(p_degree_tumor, 
                        p_MDSD_tumor, ncol = 2)

gridExtra::grid.arrange(p_degree_normal, 
                        p_MDSD_normal, 
                        p_degree_tumor, 
                        p_MDSD_tumor, 
                        ncol = 2)

png("real_data_examples/BC_TCGA/Degree_MDSD.png",
    units = "in",
    res = 300,
    width = 7,
    height = 4)

gridExtra::grid.arrange(p_degree_normal, 
                        p_MDSD_normal, 
                        p_degree_tumor, 
                        p_MDSD_tumor, 
                        ncol = 2)

dev.off()

tiff("real_data_examples/BC_TCGA/Degree_MDSD.tif",
    units = "in",
    res = 300,
    width = 7,
    height = 4)

gridExtra::grid.arrange(p_degree_normal, 
                        p_MDSD_normal, 
                        p_degree_tumor, 
                        p_MDSD_tumor, 
                        ncol = 2)

dev.off()


#####

# ROC analysis,

library(pROC)

hubs_normal = detection_res_normal$hub_nodes_MDSD
hubs_tumor = detection_res_tumor$hub_nodes_MDSD

normal_tumor_data = rbind(normal_data, 
                          tumor_data)

tumor = rep(c(0, 1), c(nrow(normal_data),
                       nrow(tumor_data)))

#tumor = as.factor(tumor)

normal_tumor_data = cbind(tumor, normal_tumor_data)

normal_tumor_data = as.data.frame(normal_tumor_data)

str(normal_tumor_data[, 1:5])

## All tumor specific MDSD hubs in the same model,

set.seed(1)

colnames(normal_tumor_data) = 
  str_replace_all(colnames(normal_tumor_data), '-', '')

N = nrow(normal_tumor_data)

train_ind = sample(1:N, floor(0.7*N))

train = normal_tumor_data[train_ind, ]
test = normal_tumor_data[-train_ind, ]

hubs_tumor_specific = setdiff(hubs_tumor, hubs_normal)

d = as.formula(paste0("tumor ~ ", 
                      paste(hubs_tumor_specific, collapse=" + ")))

glm_tumor = glm(d, family = binomial(logit), 
              data = train)

summary(glm_tumor)

# train accuracy

roc_train = pROC::roc(glm_tumor$y, glm_tumor$fitted.values, ci = TRUE)
roc_train

# test accuracy

tumor_pred = predict(glm_tumor, 
                     newdata = test,
                     type = "response")

tumor_pred = ifelse(tumor_pred > 0.5, 1, 0)

roc_test = pROC::roc(tumor_pred, test$tumor, ci = TRUE)
roc_test

test_AUC_t = roc_test$auc

roc_tumor = ggroc(roc_test, size = 1) +
  geom_segment(aes(x = 1, 
                   xend = 0, 
                   y = 0, 
                   yend = 1),
               linetype = 2,
               linewidth = 1) +
  geom_text(aes(0, 0.1, label = round(test_AUC_t, 3)),
            hjust = 1,
            size = 7) 
  

roc_tumor = roc_tumor + labs(title = "(A)")

roc_tumor

## All normal specific MDSD hubs in the same model

hubs_normal_specific = setdiff(hubs_normal, hubs_tumor)

hubs_normal_specific = 
  str_replace_all(hubs_normal_specific, '-', '')

d = as.formula(paste0("tumor ~ ", 
                      paste(hubs_normal_specific, collapse=" + ")))

glm_normal = glm(d, family = binomial(logit), 
                 data = train)

summary(glm_normal)

# train accuracy

roc_train = pROC::roc(glm_normal$y, 
                      glm_normal$fitted.values, ci = TRUE)
roc_train

# test accuracy

tumor_pred = predict(glm_normal, 
                     newdata = test,
                     type = "response")

tumor_pred = ifelse(tumor_pred > 0.5, 1, 0)

roc_test = pROC::roc(tumor_pred, test$tumor, ci = TRUE)
roc_test

test_AUC_n = roc_test$auc

roc_normal = ggroc(roc_test, size = 1) +
  geom_segment(aes(x = 1, 
                   xend = 0, 
                   y = 0, 
                   yend = 1),
               linetype = 2,
               linewidth = 1) +
  geom_text(aes(0, 0.1, label = round(test_AUC_n, 3)),
            hjust = 1,
            size = 7) 


roc_normal = roc_normal + labs(title = "(B)")

roc_normal

## Hubs from both classes in the same model,

hubs_normal_tumor = union(hubs_normal, hubs_tumor)

hubs_normal_tumor = 
  str_replace_all(hubs_normal_tumor, '-', '')

d = as.formula(paste0("tumor ~ ", 
                      paste(hubs_normal_tumor, collapse=" + ")))

glm_tumor_normal = glm(d, family = binomial(logit), 
                       data = train)

summary(glm_tumor_normal)

# train accuracy

roc_train = pROC::roc(glm_tumor_normal$y, 
                      glm_tumor$fitted.values, ci = TRUE)
roc_train

# test accuracy

tumor_pred = predict(glm_tumor_normal, 
                     newdata = test,
                     type = "response")

tumor_pred = ifelse(tumor_pred > 0.5, 1, 0)

roc_test = pROC::roc(tumor_pred, test$tumor, ci = TRUE)
roc_test

test_AUC_tn = roc_test$auc

roc_tumor_normal = ggroc(roc_test, size = 1) +
  geom_segment(aes(x = 1, 
                   xend = 0, 
                   y = 0, 
                   yend = 1),
               linetype = 2,
               linewidth = 1) +
  geom_text(aes(0, 0.1, label = round(test_AUC_tn, 3)),
            hjust = 1,
            size = 7) 


roc_tumor_normal = roc_tumor_normal + labs(title = "(C)")

roc_tumor_normal

## All genes in one model,

glm_all = glm(tumor ~ ., family = binomial(logit), 
              data = train)

# train accuracy

roc_train = pROC::roc(glm_all$y, 
                      glm_all$fitted.values, ci = TRUE)
roc_train

# test accuracy

tumor_pred = predict(glm_all, 
                     newdata = test,
                     type = "response")

tumor_pred = ifelse(tumor_pred > 0.5, 1, 0)

roc_test = pROC::roc(tumor_pred, test$tumor, ci = TRUE)
roc_test

test_AUC_all = roc_test$auc

roc_all = ggroc(roc_test, size = 1) +
  geom_segment(aes(x = 1, 
                   xend = 0, 
                   y = 0, 
                   yend = 1),
               linetype = 2,
               linewidth = 1) +
  geom_text(aes(0, 0.1, label = round(test_AUC_all, 3)),
            hjust = 1,
            size = 7) 


roc_all = roc_all + labs(title = "(D)")

roc_all

##

roc_tumor = roc_tumor + 
  theme(axis.text.x = element_text(angle=45, size=7))

roc_normal = roc_normal + 
  theme(axis.text.x = element_text(angle=45, size=7))

roc_tumor_normal = roc_tumor_normal + 
  theme(axis.text.x = element_text(angle=45, size=7))

roc_all = roc_all + 
  theme(axis.text.x = element_text(angle=45, size=7))

gridExtra::grid.arrange(roc_tumor, 
                        roc_normal,
                        roc_tumor_normal,
                        roc_all,
                        ncol = 4)

png("real_data_examples/BC_TCGA/tumor_normal_ROC.png",
    units = "in",
    res = 300,
    width = 7,
    height = 4)

gridExtra::grid.arrange(roc_tumor, 
                        roc_normal,
                        roc_tumor_normal,
                        roc_all,
                        ncol = 4)

dev.off()

tiff("real_data_examples/BC_TCGA/tumor_normal_ROC.tif",
     units = "in",
     res = 300,
     width = 7,
     height = 4)

gridExtra::grid.arrange(roc_tumor, 
                        roc_normal,
                        roc_tumor_normal,
                        roc_all,
                        ncol = 4)

dev.off()

# What about if we look each hub gene one by one?

all_formulas = lapply(colnames(normal_tumor_data[, -1]),
                      function(gname){
                        formula(paste0('tumor ~ ', gname) ,env = globalenv())
                      } 
)

tumor_formulas = lapply(colnames(normal_tumor_data[, hubs_tumor_specific]),
                       function(gname){
                         formula(paste0('tumor ~ ', gname) ,env = globalenv())
                       } 
)

normal_formulas = lapply(colnames(normal_tumor_data[, hubs_normal_specific]),
                         function(gname){
                           formula(paste0('tumor ~ ', gname) ,env = globalenv())
                        } 
)

# common hub genes,

hubs_common = intersect(hubs_normal, hubs_tumor)

common_formulas = lapply(colnames(normal_tumor_data[, hubs_common]),
                         function(gname){
                           formula(paste0('tumor ~ ', gname) ,env = globalenv())
                         } 
)

AUC_all = rep(0, ncol(tumor_data[, -1]))

for(i in 1:length(all_formulas)){
  
  roc_temp = suppressMessages(roc(all_formulas[[i]], data = normal_tumor_data))
  
  AUC_all[i] = c(auc(roc_temp))
  
}

AUC_tumor = rep(0, length(tumor_formulas))

for(i in 1:length(tumor_formulas)){
  
  roc_temp = suppressMessages(roc(tumor_formulas[[i]], data = normal_tumor_data))
  
  AUC_tumor[i] = c(auc(roc_temp))
  
}

AUC_normal = rep(0, length(normal_formulas))

for(i in 1:length(normal_formulas)){
  
  roc_temp = suppressMessages(roc(normal_formulas[[i]], data = normal_tumor_data))
  
  AUC_normal[i] = c(auc(roc_temp))
  
}

AUC_common = rep(0, length(common_formulas))

for(i in 1:length(common_formulas)){
  
  roc_temp = suppressMessages(roc(common_formulas[[i]], data = normal_tumor_data))
  
  AUC_common[i] = c(auc(roc_temp))
  
}

mean(AUC_tumor)
mean(AUC_normal)
mean(AUC_common)
mean(AUC_all)

AUC_g = c(AUC_all, AUC_tumor, AUC_normal, AUC_common)

b = c(length(AUC_all), 
      length(AUC_tumor), 
      length(AUC_normal),
      length(AUC_common))

geneset = rep(c("All", 
                "Tumor-specific", 
                "Normal-specific", 
                "Common-hubs"), b)

data_box = data.frame(geneset = geneset, AUC = AUC_g)

# Blot violin plots (= the kernel probability density)
# Median and quartiles are illustrated using box-plots,

p_boxplot = ggplot(data = data_box,
                   mapping = aes(x = reorder(geneset, AUC, FUN = median, decreasing = TRUE), 
                                 y = AUC, group = geneset)) +
  geom_violin() +
  geom_boxplot(width = 0.2) +
  #ggtitle("AUC of MDSD and hglasso genes") + 
  ylab("AUC") + 
  xlab("")

p_boxplot

png("real_data_examples/BC_TCGA/AUC_boxplot.png",
    units = "in",
    res = 300,
    width = 7,
    height = 4)

p_boxplot

dev.off()

tiff("real_data_examples/BC_TCGA/AUC_boxplot.tif",
     units = "in",
     res = 300,
     width = 7,
     height = 4)

p_boxplot

dev.off()

##

t.test(AUC_tumor, AUC_all)
wilcox.test(AUC_tumor, AUC_all, alternative = "two.sided")
