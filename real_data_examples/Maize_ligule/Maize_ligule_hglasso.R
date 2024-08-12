
# 31.01.2024

# Examine Maize and Ligule Development 
# (https://doi.org/10.1105/tpc.114.132688)

# The data set is publicly available here,

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE61333

library(tidyverse)
library(DESeq2)
library(genefilter)

library(hglasso)

source("functions/hub_detection_hglasso.R")

data = readr::read_delim("real_data_examples/Maize_ligule/GSE61333_ligule_count.txt", delim = "\t")

dim(data)

names(data)

names(data)[1] = "GeneId"

names(data)

# Data modifications: normalize expressions,

DE_input_data = as.matrix(data[, -1])

row.names(DE_input_data) = data$GeneId

meta_data = data.frame( Sample = names(data[-1])) %>%
  mutate(
    Type = gsub("-.*","", Sample) %>% gsub("[.].*","", .)
  )

meta_data$Type = as.factor(meta_data$Type)

dds = DESeq2::DESeqDataSetFromMatrix(round(DE_input_data),
                                     meta_data,
                                     design = ~Type)

dds = DESeq(dds)

# Apply a variance stabilizing transformation (VST) 
# to the count data

vsd = varianceStabilizingTransformation(dds)

wpn_vsd = getVarianceStabilizedData(dds)
rv_wpn = rowVars(wpn_vsd)

summary(rv_wpn)

dim(wpn_vsd)

# There is just a handful of treatments...

# Select gene probes with a large variance to reduce the number
# of genes. In particular, select 1000 genes with the highest
# variance,

top_vars = apply(wpn_vsd, 1, var)
top_vars = sort(top_vars, decreasing = TRUE)[1:1000]
top_vars = names(top_vars)
expr_normalized = wpn_vsd[rownames(wpn_vsd) %in% top_vars, ]

dim(expr_normalized)

head(expr_normalized)

# The hub network analysis,

expr_normalized = t(expr_normalized)

dim(expr_normalized)

p = ncol(expr_normalized)
n = nrow(expr_normalized)

S = cor(expr_normalized)

S_abs = abs(S[upper.tri(S)])

min(S_abs)
max(S_abs)

boxplot(S_abs)

quantile(S_abs, probs = c(0.1, 0.5, 0.9))

lambda1 = 0.9 # Off-diagonal elements
lambda2s = seq(0.7, 0.9, length.out = 10) # Sparsity of the hub nodes
lambda3 = 20 # Selection of hub nodes

BICcriterion = rep(0, length(lambda2s))

hglassopath = list(wi = array(0, c(p, p, length(lambda2s))),
                   rholist = 1:length(lambda2s))

hglassohubs = list()

k = 1

for(lambda2 in lambda2s){
  
  hglasso_est = hglasso(S, 
                        lambda1 = lambda1,
                        lambda2 = lambda2,
                        lambda3 = lambda3
  )
  
  BICcriterion[k] = hglassoBIC(hglasso_est, S)$BIC
  
  hglassopath$wi[,,k] = hglasso_est$Theta
  
  hglassohubs[[k]] = hglasso_est$hubind
  
  cat("\r", k)
  
  k = k + 1
  
}

#save.image(file = "real_data_examples/Maize_ligule/example.RData")
#load(file = "real_data_examples/Maize_ligule/example.RData")    

plot(lambda2s, BICcriterion)

ind = which.min(BICcriterion)

lambda2s[ind]

hglassohubs[[ind]]

hglassoTheta = hglassopath$wi[ , , ind]

node_names = colnames(expr_normalized)

detection_res = hub_detection_hglasso(hglasso_path = hglassopath,
                                      node_names = node_names,
                                      gamma = 3)

detection_res$hub_nodes_MDSD

degrees = c(t(detection_res$Degree))

Degree_data = detection_res$Degree

all_nodes = rep(node_names, each = ncol(Degree_data))

#top_nodes = sort(detection_res$MDSD, decreasing = TRUE)[1:10]
#top_nodes = names(top_nodes)

top_nodes = detection_res$hub_nodes_MDSD

node_col = rep("#CCCCCC", length(all_nodes))

node_col[all_nodes %in% top_nodes] = "#F8766D"

node_col = as.factor(node_col)

node_col = relevel(node_col, ref = "#F8766D")

lty = ifelse(node_col == "#F8766D", "solid", "dashed")

lty = as.factor(lty)

lty = relevel(lty, ref = "solid")

Degree_data = data.frame(node = all_nodes,
                         Degree = degrees,
                         lambda = lambda2s,
                         node_col = node_col,
                         lty = lty
                         )

# Plot results,

p_degree = ggplot(Degree_data, 
                  aes(x = lambda, y = Degree, group = node)
                  ) +
  geom_line(aes(colour = node_col, linetype = lty)) +
  scale_colour_identity() +
  ylab("Node degree") +
  xlab("Tuning parameter value")

p_degree = p_degree + 
  labs(title = "(A)") + 
  theme_classic() +
  theme(legend.position = "none")

p_degree

#

node_names_trim = stringr::str_remove(node_names, "GRMZM")

top_nodes = stringr::str_remove(top_nodes, "GRMZM")

node_col = rep("#CCCCCC", p)

names(node_col) = node_names_trim

node_col[top_nodes] = "#F8766D"

node_col = as.factor(node_col)

node_col = relevel(node_col, ref = "#F8766D")

lty = ifelse(node_col == "#F8766D", "solid", "dashed")

lty = as.factor(lty)

lty = relevel(lty, ref = "solid")

MDSD_data = data.frame(node = node_names_trim,
                       MDSD = detection_res$MDSD,
                       node_col = node_col,
                       lty = lty
)

p_MDSD = ggplot(MDSD_data, aes(x = node, y = MDSD)) +
  geom_segment(
    aes(x = node, xend = node, y = 0, yend = MDSD, 
        colour = node_col, linetype = lty)
  ) +
  scale_x_discrete(breaks = top_nodes,
                   labels = top_nodes
                   ) +
  scale_colour_identity() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 7)
  ) +
  #ylab("Mean Degree Squared Difference of degree path") +
  ylab("MDSD") +
  xlab("Gene name")

p_MDSD = p_MDSD + 
  labs(title = "(B)") + 
  theme_classic() +
  theme(legend.position = "none")

p_MDSD

p_MDSD = p_MDSD + 
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank()
        )

p_MDSD

gridExtra::grid.arrange(p_degree, p_MDSD, ncol = 2)

tiff("manuscript_figures/Fig4.tif",
     units = "in",
     res = 300,
     width = 7,
     height = 4)

gridExtra::grid.arrange(p_degree, p_MDSD, ncol = 2)

dev.off()

pdf("manuscript_figures/Fig4.pdf",
    width = 7,
    height = 4)

gridExtra::grid.arrange(p_degree, p_MDSD, ncol = 2)

dev.off()

png("manuscript_figures/Fig4.png",
    units = "in",
    res = 300,
    width = 7,
    height = 4)

gridExtra::grid.arrange(p_degree, p_MDSD, ncol = 2)

dev.off()

#

hubs_hglasso = node_names[hglassohubs[[ind]]]
hubs_MDSD = detection_res$hub_nodes_MDSD

length(hubs_hglasso)
length(hubs_MDSD)

length(lubridate::intersect(hubs_hglasso, hubs_MDSD))

lubridate::setdiff(hubs_MDSD, hubs_hglasso)

lubridate::setdiff(hubs_hglasso, hubs_MDSD)

#

cat(
  paste0(stringr::str_remove(hubs_MDSD, "GRMZM"), ",")
)

cat(
  paste0("\\textbf{",
    stringr::str_remove(
      lubridate::setdiff(hubs_hglasso, hubs_MDSD), "GRMZM"
      ), "},"
    )
)

###

expr_normalized_hubs = expr_normalized[, hubs_hglasso]

sample_names = rownames(expr_normalized_hubs)

mutants = stringr::str_detect(sample_names, "lg")

mutants = 1*mutants

mutants_df = data.frame(mutant = mutants)

dim(mutants_df)
dim(expr_normalized_hubs)

mutants_df = dplyr::bind_cols(mutants_df, expr_normalized_hubs)

# ROC analysis,

library(pROC)

pROC::roc(mutants_df$mutant, 
          mutants_df$GRMZM2G005818, ci = TRUE)

# This is the same as,

glm_test = glm(mutant ~ GRMZM2G005818, family = binomial(logit), 
               data = mutants_df)

pROC::roc(glm_test$y, glm_test$fitted.values, ci = TRUE)

# Add interactions,

glm_test = glm(mutant ~ GRMZM2G005818*GRMZM2G007953, family = binomial(logit), 
              data = mutants_df)

summary(glm_test)

roc_test = pROC::roc(glm_test$y, glm_test$fitted.values, ci = TRUE)

roc_test

ggroc(roc_test) +
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), linetype="dashed")


## All variables in the same model (high-dimensional data),

glm_all = glm(mutant ~ ., family = binomial(logit), 
              data = mutants_df)

summary(glm_all)

pROC::roc(glm_all$y, glm_all$fitted.values, ci = TRUE)

#

df_temp = mutants_df[, c("mutant", hubs_MDSD)]

glm_MDSD = glm(mutant ~ ., family = binomial(logit), 
               data = df_temp)

pROC::roc(glm_MDSD$y, glm_MDSD$fitted.values, ci = TRUE)

#

df_temp = mutants_df[, c("mutant", hglasso_only)]

glm_hglasso = glm(mutant ~ ., family = binomial(logit), 
                  data = df_temp)

pROC::roc(glm_hglasso$y, glm_hglasso$fitted.values, ci = TRUE)

# What about if we look each hub gene one by one?

roc_test = pROC::roc(mutant ~ GRMZM2G005818, data = mutants_df)

auc(roc_test)

ggroc(roc_test, size = 1) +
  xlab("Specificity") +
  ylab("Sensitivity") +
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), linetype="dashed")

# or,

roc_test = pROC::roc(mutants_df$mutant, 
                     mutants_df$GRMZM2G005818)

auc(roc_test)

ggroc(roc_test, size = 1) +
  xlab("Specificity") +
  ylab("Sensitivity") +
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), linetype="dashed")

###

# GRMZM2G030123 found only with hglasso,

roc_test = pROC::roc(mutant ~ GRMZM2G005818 + GRMZM2G030123, 
                     data = mutants_df)

ggroc(roc_test)

auc(roc_test$GRMZM2G005818)
auc(roc_test$GRMZM2G030123)

roc_test = pROC::roc(mutant ~ GRMZM2G005818, 
                     data = mutants_df)

auc(roc_test)

roc_test = pROC::roc(mutant ~ GRMZM2G030123, 
                     data = mutants_df)

auc(roc_test)

ggroc(roc_test, size = 1) +
  xlab("Specificity") +
  ylab("Sensitivity") +
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), linetype="dashed")


all_formulas = lapply(colnames(mutants_df[, -1]),
                       function(gname){
                         formula(paste0('mutant ~ ', gname) ,env = globalenv())
                       } 
)

MDSD_formulas = lapply(colnames(mutants_df[, hubs_MDSD]),
       function(gname){
         formula(paste0('mutant ~ ', gname) ,env = globalenv())
       } 
)

hglasso_only = lubridate::setdiff(hubs_hglasso, hubs_MDSD)

hglasso_formulas = lapply(colnames(mutants_df[, hglasso_only]),
       function(gname){
         formula(paste0('mutant ~ ', gname), env = globalenv())
       } 
)

AUC_all = rep(0, ncol(mutants_df[, -1]))

for(i in 1:length(all_formulas)){
  
  roc_temp = suppressMessages(roc(all_formulas[[i]], data = mutants_df))
  
  AUC_all[i] = c(auc(roc_temp))
  
}

AUC_MDSD = rep(0, length(MDSD_formulas))

for(i in 1:length(MDSD_formulas)){
  
  roc_temp = suppressMessages(roc(MDSD_formulas[[i]], data = mutants_df))

  AUC_MDSD[i] = c(auc(roc_temp))
    
}

AUC_hglasso = rep(0, length(hglasso_formulas))

for(i in 1:length(hglasso_formulas)){
  
  roc_temp = suppressMessages(roc(hglasso_formulas[[i]], data = mutants_df))
  
  AUC_hglasso[i] = c(auc(roc_temp))
  
}

AUC_MDSD
AUC_hglasso

AUC_g = c(AUC_all, AUC_hglasso, AUC_MDSD)

b = c(length(AUC_all), length(AUC_hglasso), length(AUC_MDSD))

geneset = rep(c("All", "hglasso", "MDSD"), b)

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

png("manuscript_figures/Fig5.png",
    units = "in",
    res = 300,
    width = 7,
    height = 4)

p_boxplot

dev.off()
