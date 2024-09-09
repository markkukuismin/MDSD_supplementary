
# 31.01.2024

library(tidyverse)
library(gridExtra)
library(ggpubr)

files_hglasso = dir("simulations/hglasso_results/")
files_cor = dir("simulations/cor_screening_results/")
files_fh = dir("simulations/Firouzi_Hero_results/")

files_hglasso = files_hglasso[str_detect(files_hglasso, ".txt")]
files_cor = files_cor[str_detect(files_cor, ".txt")]
files_fh = files_fh[str_detect(files_fh, ".txt")]

G_type = c("hglasso_n=100_p=500", "hglasso_n=500_p=1500", 
           "hglassotwo_n=100_p=500", "hglassotwo_n=500_p=1500",
           "inter_n=100_p=500", "inter_n=500_p=1500",
           "scale-free_n=100_p=500", "scale-free_n=500_p=1500",
           "star_n=100_p=500", "star_n=500_p=1500")

gtype = rep(c("hglasso", "hglasso_two", "inter", "scale-free", "star"), each = 2)

n = rep(c("100", "1000"), times = 5)

p = rep(c("500", "1500"), times = 5)

Results = tibble()

Results$Gtype = " "

Results$p = 0

Results$n = 0

for(i in 1:length(G_type)){

  temp_res_hglasso = read.table(paste0("simulations/hglasso_results/", files_hglasso[i]), header = TRUE)
  temp_res_cor = read.table(paste0("simulations/cor_screening_results/", files_cor[i]), header = TRUE)
  temp_res_fh = read.table(paste0("simulations/Firouzi_Hero_results/", files_fh[i]), header = TRUE)
  
  temp_res_hglasso = temp_res_hglasso %>%
    mutate(Criterion = if_else(Criterion == "BIC", "hglasso (BIC)", "hglasso (MDSD)"))
  
  temp_res_cor$Criterion = "cor (MDSD)"
  
  temp_res_fh$Criterion = "cor (FH)"
  
  temp_res_cor_fh = rbind(temp_res_cor, temp_res_fh)
  
  temp_res_hglasso$Gtype = gtype[i]
  temp_res_cor_fh$Gtype = gtype[i]
  
  temp_res_hglasso$p = p[i]
  temp_res_cor_fh$p = p[i]
  
  temp_res_hglasso$n = n[i]
  temp_res_cor_fh$n = n[i]
 
  temp_res_hglasso$p = factor(temp_res_hglasso$p, levels = c("500", "1500"))
  temp_res_cor_fh$p = factor(temp_res_cor_fh$p, levels = c("500", "1500"))
  
  temp_res_hglasso$n = factor(temp_res_hglasso$n, levels = c("100", "1000"))
  temp_res_cor_fh$n = factor(temp_res_cor_fh$n, levels = c("100", "1000"))
   
  Results = rbind(Results, temp_res_hglasso, temp_res_cor_fh)
    
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

# Manuscript Fig 2,

Results_temp_f2 = Results_temp %>%
  filter(Method %in% c("hglasso (BIC)", "hglasso (MDSD)"))

i = 1

p_box = p_box2 = list()

gtype = unique(gtype)

pdf(file = "manuscript_figures/Sup_Fig2.pdf", 
    onefile = TRUE, paper = "a4")

for(g in gtype){
  
  temp_res = Results_temp_f2 %>% dplyr::filter(Gtype == g)
  
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
  
  p_box2[[i]] = p_box[[i]]
  
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

##

dummy = temp_res

dummy$Value = runif(nrow(dummy))

p_box2[[6]] = ggplot(dummy, aes(Value, fill = Method)) +
  geom_bar() +
  theme_void() +
  theme(legend.position = c(0.5, 0.5),
        legend.direction = "horizontal",
        legend.title = element_text(face = "bold")) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5))

p_box2[[6]]

for(i in 1:length(gtype)){
  p_box2[[i]] = p_box2[[i]] + 
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  if(i <= 5){
    p_box2[[i]] = p_box2[[i]] + 
      theme(legend.position = "none")
  }
  
  if(i %in% c(2, 3, 5)){
    p_box2[[i]] = p_box2[[i]] + 
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
  }
  
  if(i %in% c(1, 2, 4)){
    p_box2[[i]] = p_box2[[i]] + 
      theme(strip.background = element_blank(),
            strip.text.y = element_blank())
  }
  
}

grid.arrange(p_box2[[1]], p_box2[[2]], p_box2[[3]], 
             p_box2[[4]], p_box2[[5]], p_box2[[6]], ncol = 3)

tiff("manuscript_figures/Fig2.tif",
     units = "in",
     res = 300,
     width = 8.5,
     height = 5.5)

grid.arrange(p_box2[[1]], p_box2[[2]], p_box2[[3]], 
             p_box2[[4]], p_box2[[5]], p_box2[[6]], ncol = 3)

dev.off()

pdf("manuscript_figures/Fig2.pdf",
    width = 10,
    height = 6)

grid.arrange(p_box2[[1]], p_box2[[2]], p_box2[[3]], 
             p_box2[[4]], p_box2[[5]], p_box2[[6]], ncol = 3)

dev.off()

png("manuscript_figures/Fig2.png",
    units = "in",
    res = 300,
    width = 8.5,
    height = 5.5)

grid.arrange(p_box2[[1]], p_box2[[2]], p_box2[[3]], 
             p_box2[[4]], p_box2[[5]], p_box2[[6]], ncol = 3)

dev.off()

# Manuscript Fig 3,

Results_temp_f3 = Results_temp %>%
  filter(Method %in% c("cor (FH)", "cor (MDSD)"))

i = 1

p_box = p_box2 = list()

gtype = unique(gtype)

pdf(file = "manuscript_figures/Sup_Fig3.pdf", 
    onefile = TRUE, paper = "a4")

for(g in gtype){
  
  temp_res = Results_temp_f3 %>% dplyr::filter(Gtype == g)
 
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
  
  p_box2[[i]] = p_box[[i]]
  
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

##

dummy = temp_res

dummy$Value = runif(nrow(dummy))

p_box2[[6]] = ggplot(dummy, aes(Value, fill = Method)) +
  geom_bar() +
  theme_void() +
  theme(legend.position = c(0.5, 0.5),
        legend.direction = "horizontal",
        legend.title = element_text(face = "bold")) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5))

p_box2[[6]]

for(i in 1:length(gtype)){
  p_box2[[i]] = p_box2[[i]] + 
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  if(i <= 5){
    p_box2[[i]] = p_box2[[i]] + 
      theme(legend.position = "none")
  }
  
  if(i %in% c(2, 3, 5)){
    p_box2[[i]] = p_box2[[i]] + 
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
  }
  
  if(i %in% c(1, 2, 4)){
    p_box2[[i]] = p_box2[[i]] + 
      theme(strip.background = element_blank(),
            strip.text.y = element_blank())
  }
  
}

grid.arrange(p_box2[[1]], p_box2[[2]], p_box2[[3]], 
             p_box2[[4]], p_box2[[5]], p_box2[[6]], ncol = 3)

tiff("manuscript_figures/Fig3.tif",
     units = "in",
     res = 300,
     width = 8.5,
     height = 5.5)

grid.arrange(p_box2[[1]], p_box2[[2]], p_box2[[3]], 
             p_box2[[4]], p_box2[[5]], p_box2[[6]], ncol = 3)

dev.off()

pdf("manuscript_figures/Fig3.pdf",
    width = 10,
    height = 6)

grid.arrange(p_box2[[1]], p_box2[[2]], p_box2[[3]], 
             p_box2[[4]], p_box2[[5]], p_box2[[6]], ncol = 3)

dev.off()

png("manuscript_figures/Fig3.png",
    units = "in",
    res = 300,
    width = 8.5,
    height = 5.5)

grid.arrange(p_box2[[1]], p_box2[[2]], p_box2[[3]], 
             p_box2[[4]], p_box2[[5]], p_box2[[6]], ncol = 3)

dev.off()


##

# Estimated values,

# MCC

D = Results_temp %>% 
  filter(Metric == "MCC" & p == 500) %>% 
  group_by(Method, Gtype) %>%
  summarise(MCC_avg = mean(Value), sd = sd(Value))

D %>%
  group_by(Gtype) %>%
  arrange(Gtype, desc(MCC_avg))

D = Results_temp %>% 
  filter(Metric == "MCC" & p == 1500) %>% 
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

D = Results_temp %>% 
  filter(Metric == "TPR" & p == 1500) %>% 
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

D = Results_temp %>% 
  filter(Metric == "BM" & p == 1500) %>% 
  group_by(Method, Gtype) %>%
  summarise(BM_avg = mean(Value), sd = sd(Value))

D %>%
  group_by(Gtype) %>%
  arrange(Gtype, desc(BM_avg))
