
# 09.02.2024

library(tidyverse)
library(gridExtra)
library(ggpubr)

files_space = dir("simulations/space_results/")

files_space = files_space[str_detect(files_space, ".txt")]

G_type = c("hglasso_n=100_p=500", 
           "hglassotwo_n=100_p=500",
           "inter_n=100_p=500",
           "scale-free_n=100_p=500",
           "star_n=100_p=500")

gtype = c("hglasso", "hglasso_two", "inter", "scale-free", "star")

n = rep("100", times = 5)

p = rep("500", times = 5)

Results = tibble()

Results$Gtype = " "

Results$p = 0

Results$n = 0

for(i in 1:length(G_type)){
  
  temp_res_space = read.table(paste0("simulations/space_results/", files_space[i]), header = TRUE)

  temp_res_space = temp_res_space %>%
    mutate(Criterion = if_else(Criterion == "BIC", "space (BIC)", "space (MDSD)"))
  
  temp_res_space$Gtype = gtype[i]

  temp_res_space$p = p[i]

  temp_res_space$n = n[i]

  temp_res_space$p = factor(temp_res_space$p, levels = c("500", "1500"))

  temp_res_space$n = factor(temp_res_space$n, levels = c("100", "1000"))

  Results = rbind(Results, temp_res_space)
  
}

Results = Results %>%
  mutate(BM = TPR - FPR)

Results = Results %>%
  select(!c(FP, FN, TN, TP))

Results_temp = Results %>%
  tidyr::pivot_longer(cols = c(FDR, TPR, Pre, FPR, MCC, BM), 
                      names_to = "Metric", values_to = "Value")


Results_temp = Results_temp %>%
  rename(Method = Criterion)

# Suplementary Fig 4(?),

i = 1

p_box = list()

gtype = unique(gtype)

pdf(file = "manuscript_figures/Sup_Fig4.pdf", 
    onefile = TRUE, paper = "a4")

for(g in gtype){
  
  temp_res = Results_temp %>% dplyr::filter(Gtype == g)
  
  p_box[[i]] = ggplot(temp_res, aes(x = Metric, y = Value)) +
    geom_boxplot(aes(fill = Method), position = position_dodge(0.8)) +
    facet_wrap(~ Method) +
    facet_grid(vars(p)) +
    labs(title = g)
  
  if(g == "hglasso_two") p_box[[i]] = p_box[[i]] + labs(title = "Hub network, two components")
  if(g == "hglasso") p_box[[i]] = p_box[[i]] + labs(title = "Hub network")
  if(g == "star") p_box[[i]] = p_box[[i]] + labs(title = "Star network")
  if(g == "scale-free") p_box[[i]] = p_box[[i]] + labs(title = "Scale-free")
  if(g == "inter") p_box[[i]] = p_box[[i]] + labs(title = "Inter-hub network")
  
  p_box[[i]] = p_box[[i]] + 
    stat_compare_means(aes(group = Method, 
                           label = ifelse(
                             p < 0.01,
                             "p < 0.01",
                             sprintf("p = %5.2f", as.numeric(..p.format..)))), 
                       label.y = 1.2,
                       method = "wilcox.test")
  
  print(p_box[[i]])
  
  i = i + 1
  
}

dev.off()

# Estimated values,

# MCC

D = Results_temp %>% 
  filter(Metric == "MCC" & p == 500) %>% 
  group_by(Method, Gtype) %>%
  summarise(MCC_avg = mean(Value), sd = sd(Value))

D %>%
  group_by(Gtype) %>%
  arrange(Gtype, desc(MCC_avg))

# Power,

D = Results_temp %>% 
  filter(Metric == "TPR" & p == 500) %>% 
  group_by(Method, Gtype) %>%
  summarise(TPR_avg = mean(Value), sd = sd(Value))

D %>%
  group_by(Gtype) %>%
  arrange(Gtype, desc(TPR_avg))

# Bookmaker Informedness,

D = Results_temp %>% 
  filter(Metric == "BM" & p == 500) %>% 
  group_by(Method, Gtype) %>%
  summarise(BM_avg = mean(Value), sd = sd(Value))

D %>%
  group_by(Gtype) %>%
  arrange(Gtype, desc(BM_avg))
