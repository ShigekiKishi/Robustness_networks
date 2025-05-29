

#  ------------------------------------------------------------------
#* Author: Shigeki KISHI
#* Affiliation: RCAIT, NARO
#* Address: 1-31-1, Kan-non-dai, Tsukuba, Ibaraki, Japan
#* Last date: May 29, 2025
#* 

## libraries
library(tidyverse)
library(magrittr)
library(cowplot)
library(patchwork)
library(bipartite)
library(scales) # to use parcent 
library(lme4)
library(lmerTest)
library(furrr)
plan(multisession)

# -------------------------------------------------------------------------

# setting working directory
#setwd("C:/Users/kishis201/Desktop/R訪花昆虫データ/jbpress")

# file reading, preparing datasets

lf <- list.files(path = "./networkdata/", full.names = T)
lf
lf[1] %>% read_csv() %>% column_to_rownames(colnames(.)[1])

net_names <- lf %>% 
  str_remove("./networkdata/") %>% 
  str_remove(".csv") %>% 
  str_sub(1, str_length(.)-1) %>% 
  unique()
net_names
# [1] "Ashu84"   "Ashu85"   "Ashu86"   "Ashu87"   "Kibune84" "Kibune85" "Kibune86" "Kibune87" "KyotoU85"
# [10] "KyotoU86" "KyotoU87"

## original networks
datA <- lf %>% str_subset("A.csv") %>% map(\(x){
  read_csv(x) %>% 
    column_to_rownames(colnames(.)[1]) %>% 
    as.matrix()
})
## female visitor networks
datF <- lf %>% str_subset("F.csv") %>% map(\(x){
  read_csv(x) %>% 
    column_to_rownames(colnames(.)[1]) %>% 
    as.matrix()
})
## male visitor networks
datM <- lf %>% str_subset("M.csv") %>% map(\(x){
  read_csv(x) %>% 
    column_to_rownames(colnames(.)[1]) %>% 
    as.matrix()
})


#########################

## Robustness

########################

# Removing insect visitors ------------------------------------------------------------------

## female visitors
se_df_F <- 1:11 %>% map(\(x){
  flower_sp <- nrow(datF[[x]])
  se_df <- datF[[x]] %>% 
    second.extinct(method = "random", nrep = 100) %>% 
    asplit(1) %>% map(\(x)x %>% t() %>% as.data.frame()) %>% list_rbind() %>% 
    mutate(net = net_names[x]) %>% # 
    mutate(trophic = "upper") %>% 
    mutate(cm_lower = flower_sp - cumsum(ext.lower)) %>% 
    mutate(cm_lower_ratio = cm_lower/flower_sp) %>% 
    mutate(no = no/nrow(.))
  return(se_df)
}) %>% list_rbind()

p_F <- se_df_F %>% 
  ggplot(aes(no, cm_lower_ratio, group = net))+
  theme_cowplot()+
  geom_line()
p_F

## male visitors
se_df_M <- 1:11 %>% map(\(x){
  flower_sp <- nrow(datM[[x]])
  se_df <- datM[[x]] %>% 
    second.extinct(method = "random", nrep = 100) %>% 
    asplit(1) %>% map(\(x)x %>% t() %>% as.data.frame()) %>% list_rbind() %>% 
    mutate(net = net_names[x]) %>% # 
    mutate(trophic = "upper") %>% 
    mutate(cm_lower = flower_sp - cumsum(ext.lower)) %>% 
    mutate(cm_lower_ratio = cm_lower/flower_sp) %>% 
    mutate(no = no/nrow(.))
  return(se_df)
}) %>% list_rbind()

se_df_M %>% head()

p_M <- se_df_M %>% 
  ggplot(aes(no, cm_lower_ratio, group = net))+
  theme_cowplot()+
  geom_line()
p_M

p_F + p_M

#
se_df_F %>% 
  mutate(sex = "F") %>% 
  rbind(se_df_M %>% mutate(sex = "M")) %>% 
  ggplot(aes(no, cm_lower_ratio, fill = net, color = sex))+
  theme_cowplot()+
  geom_line()

# original networks
se_df_A <- 1:11 %>% map(\(x){
  flower_sp <- nrow(datA[[x]])
  se_df <- datA[[x]] %>% 
    second.extinct(method = "random", nrep = 100) %>% 
    asplit(1) %>% map(\(x)x %>% t() %>% as.data.frame()) %>% list_rbind() %>% 
    mutate(net = net_names[x]) %>% # 芦生とか
    mutate(trophic = "upper") %>% 
    mutate(cm_lower = flower_sp - cumsum(ext.lower)) %>% 
    mutate(cm_lower_ratio = cm_lower/flower_sp) %>% 
    mutate(no = no/nrow(.))
  return(se_df)
}) %>% list_rbind()

se_df_F %>% 
  mutate(sex = "F") %>% 
  rbind(se_df_M %>% mutate(sex = "M")) %>% 
  rbind(se_df_A %>% mutate(sex = "A")) %>% 
  ggplot(aes(no, cm_lower_ratio, fill = net, color = sex))+
  theme_cowplot()+
  geom_line()

##

# calc. area ---------------------------------------------------------------


se_df_F %>% 
  filter(net == "Ashu87") %>% 
  select(no, cm_lower_ratio) %>% 
  rbind(data.frame(no = 0, cm_lower_ratio = 1)) %>% 
  arrange(no) %>% 
  mutate(diagonal_diff = cm_lower_ratio - (1-no)) %>% 
  mutate(no_diff = no - lag(no, 1)) %>% 
  mutate(area = no_diff*diagonal_diff) %>% 
  drop_na() %>% 
  .$area %>% sum()
## [1] 0.2267133

lorenz_area <- se_df_F %>% 
  mutate(sex = "F") %>% 
  rbind(se_df_M %>% mutate(sex = "M")) %>% 
  rbind(se_df_A %>% mutate(sex = "A")) %>% 
  group_by(sex, net) %>% 
  nest() %>% 
  mutate(area = data %>% map(\(x){
    res <- x %>% 
      select(no, cm_lower_ratio) %>% 
      rbind(data.frame(no = 0, cm_lower_ratio = 1)) %>% 
      arrange(no) %>% 
      mutate(diagonal_diff = cm_lower_ratio - (1-no)) %>% 
      mutate(no_diff = no - lag(no, 1)) %>% 
      mutate(area = no_diff*diagonal_diff) %>% 
      drop_na() %>% 
      .$area %>% sum()
    return(res)
  })) %>% 
  select(net, sex, area) %>% 
  unnest(c(net, sex, area))

lorenz_area %>% head()


## statistical test
lm_res <- lorenz_area %>% lmer(data = ., area ~ sex+(1|net))
anova(lm_res)
# Type III Analysis of Variance Table with Satterthwaite's method
#        Sum Sq   Mean Sq NumDF DenDF F value   Pr(>F)   
# sex 0.0055449 0.0027725     2    20  6.3727 0.007224 **

## bonferroni correction
lorenz_area$sex %>% unique()

lorenz_area %>% 
  filter(sex != "A") %>% 
  lmer(data = ., area ~ sex+(1|net)) %>% 
  anova()
## no difference between F and M

lorenz_area %>% 
  filter(sex != "M") %>% 
  lmer(data = ., area ~ sex+(1|net)) %>% 
  anova()
# Sum Sq   Mean Sq NumDF DenDF F value   Pr(>F)   
# sex 0.0018945 0.0018945     1    10  14.214 0.003661 **
# dif. btw F and A

lorenz_area %>% 
  filter(sex != "F") %>% 
  lmer(data = ., area ~ sex+(1|net)) %>% 
  anova()
# Sum Sq   Mean Sq NumDF DenDF F value   Pr(>F)   
# sex 0.0054889 0.0054889     1    10  12.718 0.005128 **
## dif. btw A and M


# Removing flowers ------------------------------------------------------------------

## female networks

se_df_F2 <- 1:11 %>% map(\(x){
  col_num <- ncol(datF[[x]])
  se_df <- datF[[x]] %>% 
    second.extinct(participant = "lower", method = "random", nrep = 100) %>% 
    asplit(1) %>% map(\(x)x %>% t() %>% as.data.frame()) %>% list_rbind() %>% 
    mutate(net = net_names[x]) %>% # 芦生とか
    # mutate(trophic = "lower") %>% # .$ext.lower %>% sum()
    mutate(cm_higher = col_num - cumsum(ext.higher)) %>% 
    mutate(cm_higher_ratio = cm_higher/col_num) %>% 
    mutate(no = no/nrow(.))
  return(se_df)
}) %>% list_rbind()

se_df_F %>% head()

## male networks

se_df_M2 <- 1:11 %>% map(\(x){
  col_num <- ncol(datM[[x]])
  se_df <- datM[[x]] %>% 
    second.extinct(participant = "lower", method = "random", nrep = 100) %>% 
    asplit(1) %>% map(\(x)x %>% t() %>% as.data.frame()) %>% list_rbind() %>% 
    mutate(net = net_names[x]) %>% # 芦生とか
    # mutate(trophic = "lower") %>% # .$ext.lower %>% sum()
    mutate(cm_higher = col_num - cumsum(ext.higher)) %>% 
    mutate(cm_higher_ratio = cm_higher/col_num) %>% 
    mutate(no = no/nrow(.))
  return(se_df)
}) %>% list_rbind()

## original networks

se_df_A2 <- 1:11 %>% map(\(x){
  col_num <- ncol(datA[[x]])
  se_df <- datA[[x]] %>% 
    second.extinct(participant = "lower", method = "random", nrep = 100) %>% 
    asplit(1) %>% map(\(x)x %>% t() %>% as.data.frame()) %>% list_rbind() %>% 
    mutate(net = net_names[x]) %>% # 芦生とか
    # mutate(trophic = "lower") %>% # .$ext.lower %>% sum()
    mutate(cm_higher = col_num - cumsum(ext.higher)) %>% 
    mutate(cm_higher_ratio = cm_higher/col_num) %>% 
    mutate(no = no/nrow(.))
  return(se_df)
}) %>% list_rbind()


se_df_F2 %>% 
  mutate(sex = "F") %>% 
  rbind(se_df_M2 %>% mutate(sex = "M")) %>% 
  rbind(se_df_A2 %>% mutate(sex = "A")) %>% 
  ggplot(aes(no, cm_higher_ratio, fill = net, color = sex))+
  theme_cowplot()+
  geom_line()


lorenz_area_flower <- se_df_F2 %>% 
  mutate(sex = "F") %>% 
  rbind(se_df_M2 %>% mutate(sex = "M")) %>% 
  rbind(se_df_A2 %>% mutate(sex = "A")) %>% 
  group_by(sex, net) %>% 
  nest() %>% 
  mutate(area = data %>% map(\(x){
    res <- x %>% 
      select(no, cm_higher_ratio) %>% 
      rbind(data.frame(no = 0, cm_higher_ratio = 1)) %>% 
      arrange(no) %>% 
      mutate(diagonal_diff = cm_higher_ratio - (1-no)) %>% 
      mutate(no_diff = no - lag(no, 1)) %>% 
      mutate(area = no_diff*diagonal_diff) %>% 
      drop_na() %>% 
      .$area %>% sum()
    return(res)
  })) %>% 
  select(net, sex, area) %>% 
  unnest(c(net, sex, area))

lorenz_area_flower %>% head()


# RESULTS ----------------------------------------------------------------------

#* when removing flowers, the robustness was not different between 
#* male and female sub-networks and original networks
#
# -------------------------------------------------------------------------

## removing insect visitors

## box plot

p <- lorenz_area %>% 
  ggplot(aes(sex, area))+
  theme_cowplot()+
  geom_jitter(height = 0, width = 0.2, alpha = 0.8, size = 3)+
  geom_boxplot(fill = "transparent", coef = Inf, outliers = F)+
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    legend.position = "none"
  )+
  labs(y = "Bowl-shaped area")+
  scale_x_discrete(labels = c("A" = "All", "F" = "Female", "M" = "Male"))+
  #scale_color_manual(values = c("A" = "black", "F" = "tomato", "M" = "darkgreen"))+
  annotate("segment", x = 1, xend = 2, y = 0.12)+
  annotate("text", label = "**", x = 1.5, y = 0.11, size = 10)+
  annotate("segment", x = 1.1, xend = 2.9, y = 0.25)+
  annotate("text", label = "**", x = 2, y = 0.26, size = 10)
p  

p %>% ggsave(filename = "lorenz_areaplot.jpeg", units = "cm", 
             height = 10, width = 15, dpi = 300, bg = "white")

##


p <- se_df_F %>% 
  mutate(sex = "F") %>% 
  rbind(se_df_M %>% mutate(sex = "M")) %>% 
  rbind(se_df_A %>% mutate(sex = "A")) %>% 
  filter(net == "Ashu86") %>% 
  ggplot(aes(no, cm_lower_ratio))+
  theme_cowplot()+
  geom_line(color = "black")+
  geom_abline(slope = -1, intercept = 1, linetype = "dashed", color = "black")+
  labs(x = "Ratio of removed visitors", y = "Ratio of survived flowers")+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 10),
    strip.background = element_rect(fill = "transparent", color = "transparent"),
    strip.text = element_text(size = 10, hjust = 0),
    panel.spacing = unit(5, units = "mm"),
    #plot.margin = unit(c(5,10,5,5), units = "mm")
  )+
  coord_fixed(ratio = 1)+
  geom_ribbon(aes(x = no, ymax = cm_lower_ratio, ymin = 1-no, alpha = 0.5))+
  facet_wrap(sex~., nrow = 1, labeller = as_labeller(c("A" = "Original", "F" = "Female", "M"= "Male")))
  #scale_color_manual(values = c("A" = "black", "F" = "tomato", "M" = "green3"))+
  #scale_fill_manual(values = c("A" = "black", "F" = "tomato", "M" = "green3"))
p
p %>% ggsave(filename = "ashu86_plotarea_wrap.jpeg", 
             height = 7, width = 15, units = "cm", dpi = 300, bg = "white")



# -------------------------------------------------------------------------


## removing flowers

## test
lm_res <- lorenz_area_flower %>% lmer(data = ., area ~ sex+(1|net))
anova(lm_res)
# Type III Analysis of Variance Table with Satterthwaite's method
#        Sum Sq    Mean Sq NumDF DenDF F value    Pr(>F)    
# sex 0.0017117 0.00085587     2    20  10.009 0.0009721 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## bonferroni 補正
lorenz_area_flower$sex %>% unique()

lorenz_area_flower %>% 
  filter(sex != "A") %>% 
  lmer(data = ., area ~ sex+(1|net)) %>% 
  anova()
# Type III Analysis of Variance Table with Satterthwaite's method
#       Sum Sq  Mean Sq NumDF DenDF F value  Pr(>F)  
# sex 0.001074 0.001074     1    10   8.676 0.01465 *
## FとMには違いあり

lorenz_area_flower %>% 
  filter(sex != "M") %>% 
  lmer(data = ., area ~ sex+(1|net)) %>% 
  anova()
# Type III Analysis of Variance Table with Satterthwaite's method
#         Sum Sq    Mean Sq NumDF DenDF F value Pr(>F)
# sex 3.0075e-05 3.0075e-05     1    10  0.6618 0.4349
## 全体とFで有意差なし

lorenz_area_flower %>% 
  filter(sex != "F") %>% 
  lmer(data = ., area ~ sex+(1|net)) %>% 
  anova()
# Type III Analysis of Variance Table with Satterthwaite's method
#        Sum Sq   Mean Sq NumDF DenDF F value   Pr(>F)   
# sex 0.0014635 0.0014635     1    10  16.766 0.002163 **
## 全体とMで有意差あり

## box plot
p <- lorenz_area_flower %>% 
  ggplot(aes(sex, area))+
  theme_cowplot()+
  geom_jitter(height = 0, width = 0.2, alpha = 0.8, size = 3)+
  geom_boxplot(fill = "transparent", coef = Inf, outliers = F)+
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    legend.position = "none"
  )+
  labs(y = "Bowl-shaped area")+
  scale_x_discrete(labels = c("A" = "All", "F" = "Female", "M" = "Male"))+
  #scale_color_manual(values = c("A" = "black", "F" = "tomato", "M" = "darkgreen"))+
  annotate("segment", x = 2, xend = 2.8, y = 0.016)+
  annotate("text", label = "**", x = 2.45, y = 0.011, size = 10)+
  annotate("segment", x = 1.1, xend = 2.9, y = 0.1)+
  annotate("text", label = "**", x = 2, y = 0.103, size = 10)
p  

p %>% ggsave(filename = "lorenz_area_flower_plot.jpeg", units = "cm", 
             height = 10, width = 15, dpi = 300, bg = "white")


# -------------------------------------------------------------------------


## box plot when removing insect visitors
p_visitor <- lorenz_area %>% 
  ggplot(aes(sex, area))+
  theme_bw()+
  geom_jitter(height = 0, width = 0.2, size = 0.7)+
  geom_boxplot(fill = "transparent", coef = Inf, outliers = F)+
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    legend.position = "none",
    title = element_text(size = 10)
  )+
  labs(y = "Bowl-shaped area", title = "(a) Removing Visitors")+
  scale_x_discrete(labels = c("A" = "Original", "F" = "Female", "M" = "Male"))+
  #scale_color_manual(values = c("A" = "black", "F" = "tomato", "M" = "darkgreen"))+
  annotate("segment", x = 1, xend = 2, y = 0.12)+
  annotate("text", label = "*", x = 1.5, y = 0.11, size = 5)+
  annotate("segment", x = 1.1, xend = 2.9, y = 0.25)+
  annotate("text", label = "**", x = 2, y = 0.26, size = 5)
  #ylim(c(0, 0.3))
p_visitor  


## box plot when removing flowers

p_flower <- lorenz_area_flower %>% 
  ggplot(aes(sex, area))+
  theme_bw()+
  geom_jitter(height = 0, width = 0.2, size = 0.7)+
  geom_boxplot(fill = "transparent", coef = Inf, outliers = F)+
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    legend.position = "none",
    title = element_text(size = 10)
  )+
  labs(y = "Bowl-shaped area", title = "(b) Removing Flowers")+
  scale_x_discrete(labels = c("A" = "Original", "F" = "Female", "M" = "Male"))+
  #scale_color_manual(values = c("A" = "black", "F" = "tomato", "M" = "darkgreen"))+
  annotate("segment", x = 2, xend = 2.8, y = 0.016)+
  annotate("text", label = "**", x = 2.45, y = 0.005, size = 5)+
  annotate("segment", x = 1.1, xend = 2.9, y = 0.1)+
  annotate("text", label = "***", x = 2, y = 0.103, size = 5)
  #ylim(c(0, 0.3))
p_flower  

p_2 <- p_visitor + p_flower + plot_layout(ncol = 1, axes = "collect")
p_2

p_2 %>% 
  ggsave(filename = "ashu86_areacomparison_wrap.jpeg", 
             height = 12, width = 12, units = "cm", dpi = 300, bg = "white")



# -------------------------------------------------------------------------


## model fitting

lorenz_dat <- lorenz_area %>% 
  mutate(remove = "visitor") %>% 
  rbind(lorenz_area_flower %>% 
          mutate(remove = "flower"))

lmer_res <- lmer(data = lorenz_dat, area ~ sex + remove + (1|net))
lmer_res %>% anova()

# Type III Analysis of Variance Table with Satterthwaite's method
#          Sum Sq  Mean Sq NumDF DenDF  F value  Pr(>F)    
# sex    0.006345 0.003172     2    52   3.5505 0.03586 *  
# remove 0.310251 0.310251     1    52 347.2393 < 2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

lmer_res %>% summary()

# Fixed effects:
#   Estimate Std. Error        df t value Pr(>|t|)    
# (Intercept)    0.066156   0.008248 43.294702   8.021 4.22e-10 ***
#   sexF          -0.010449   0.009013 52.000000  -1.159   0.2516    
# sexM          -0.023952   0.009013 52.000000  -2.658   0.0104 *  
#   removevisitor  0.137124   0.007359 52.000000  18.634  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr) sexF   sexM  
# sexF        -0.546              
# sexM        -0.546  0.500       
# removevistr -0.446  0.000  0.000

TukeyHSD(lmer_res)

##
lorenz_dat %>% head()
lorenz_dat %>% 
  filter(sex != "A") %>% 
  lmer(data = ., area ~ sex + remove + (1|net)) %>% 
  anova()
# Type III Analysis of Variance Table with Satterthwaite's method
#          Sum Sq  Mean Sq NumDF DenDF  F value    Pr(>F)    
# sex    0.002006 0.002006     1    31   1.9367    0.1739    
# remove 0.191299 0.191299     1    31 184.7358 1.332e-14 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
## 性では差はない。つまり全体だとオスとメスでは違いがない。

lorenz_dat %>% 
  filter(sex != "F") %>% 
  lmer(data = ., area ~ sex + remove + (1|net)) %>% 
  anova()
# Type III Analysis of Variance Table with Satterthwaite's method
#         Sum Sq Mean Sq NumDF DenDF  F value    Pr(>F)    
# sex    0.00631 0.00631     1    31   5.9933   0.02022 *  
# remove 0.21555 0.21555     1    31 204.7212 3.352e-15 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
## 全体とオスでは違いがある。

lorenz_dat %>% 
  filter(sex != "M") %>% 
  lmer(data = ., area ~ sex + remove + (1|net)) %>% 
  anova()

# Type III Analysis of Variance Table with Satterthwaite's method
#          Sum Sq  Mean Sq NumDF DenDF F value Pr(>F)    
# sex    0.001201 0.001201     1    31   1.518 0.2272    
# remove 0.214104 0.214104     1    31 270.623 <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# -------------------------------------------------------------------------

lmer_res <- lorenz_dat %>% 
  ungroup() %>% 
  select(area, everything()) %>% 
  mutate(across(net:remove, \(x)as.factor(x))) %>% 
  lmer(data = ., area ~ remove + sex + (1|net))

lmer_res %>% anova()
ls_means(lmer_res, pairwise = T)

# Estimate Std. Error df  t value      lower      upper Pr(>|t|)    
# removeflower - removevisitor -0.1371243  0.0073587 52 -18.6344 -0.1518906 -0.1223580  < 2e-16 ***
#   sexA - sexF                   0.0104490  0.0090125 52   1.1594 -0.0076359  0.0285340  0.25159    
# sexA - sexM                   0.0239516  0.0090125 52   2.6576  0.0058667  0.0420365  0.01043 *  
#   sexF - sexM                   0.0135026  0.0090125 52   1.4982 -0.0045824  0.0315875  0.14013    

lmer_res %>% difflsmeans()


# -------------------------------------------------------------------------

lmer_visitorres <- lorenz_area %>% 
  ungroup() %>% 
  mutate(across(net:sex, \(x)as.factor(x))) %>% 
  lmer(data = ., area ~ sex + (1|net))

lmer_visitorres %>% anova()
# Type III Analysis of Variance Table with Satterthwaite's method
#        Sum Sq   Mean Sq NumDF DenDF F value   Pr(>F)   
# sex 0.0055449 0.0027725     2    20  6.3727 0.007224 **

ls_means(lmer_visitorres, pairwise = T)
# Least Squares Means table:
#   
#   Estimate  Std. Error df t value       lower       upper Pr(>|t|)   
# sexA - sexF  1.8560e-02  8.8938e-03 20  2.0868  7.4936e-06  3.7112e-02 0.049916 * 
#   sexA - sexM  3.1591e-02  8.8938e-03 20  3.5520  1.3039e-02  5.0143e-02 0.001999 **
#   sexF - sexM  1.3031e-02  8.8938e-03 20  1.4652 -5.5211e-03  3.1583e-02 0.158416  



lmer_flowerres <- lorenz_area_flower %>% 
  ungroup() %>% 
  mutate(across(net:sex, \(x)as.factor(x))) %>% 
  lmer(data = ., area ~ sex + (1|net))

lmer_flowerres %>% anova()
# Type III Analysis of Variance Table with Satterthwaite's method
#        Sum Sq    Mean Sq NumDF DenDF F value    Pr(>F)    
# sex 0.0017117 0.00085587     2    20  10.009 0.0009721 ***

ls_means(lmer_flowerres, pairwise = T)
# Estimate Std. Error df t value      lower      upper  Pr(>|t|)    
# sexA - sexF  0.0023384  0.0039430 20  0.5931 -0.0058865  0.0105633 0.5597899    
# sexA - sexM  0.0163125  0.0039430 20  4.1371  0.0080875  0.0245374 0.0005106 ***
#   sexF - sexM  0.0139740  0.0039430 20  3.5440  0.0057491  0.0221990 0.0020364 ** 



# -------------------------------------------------------------------------


# removing flowers

## connectance

connec_res <- c(datA, datF, datM) %>% map_vec(\(x)networklevel(x, index = "connectance"))
connec_df <- connec_res %>% as_tibble() %>% 
  set_colnames("connectance") %>% 
  mutate(network = rep(c("A", "F", "M"), each = 11))

connec_df %>% 
  ggplot(data=., aes(x = network, y = connectance))+
  geom_boxplot()

connec_df %>% 
  kruskal.test(data =., connectance ~ network)

# Kruskal-Wallis rank sum test
# 
# data:  connectance by network
# Kruskal-Wallis chi-squared = 4.5639, df = 2, p-value = 0.1021

connec_df %>% 
  group_by(network) %>% 
  summarise(connec_mean = mean(connectance),
            connec_sd = sd(connectance))

# # A tibble: 3 × 3
# network connec_mean connec_sd
# <chr>         <dbl>     <dbl>
#   1 A            0.0362   0.00852
# 2 F            0.0424   0.0117 
# 3 M            0.0469   0.0157 

connec_df
connec_df %>% 
  pull(connectance) %>% matrix(ncol = 3) %>% ## byrow = F
  friedman.test(y = .)

# Friedman rank sum test
# 
# data:  .
# Friedman chi-squared = 17.636, df = 2, p-value = 0.000148

connec_df
connec_df2 <- connec_df %>% 
  mutate(num = rep(1:11, 3)) %>% 
  pivot_wider(id_cols = num, names_from = network, values_from = connectance) %>% 
  select(-num) %>% 
  set_colnames(c("All", "Female", "Male"))

## all - female
connec_df2 %$% 
  wilcox.test(All, Female, paired = T)

# Wilcoxon signed rank exact test
# 
# data:  All and Female
# V = 0, p-value = 0.0009766
# alternative hypothesis: true location shift is not equal to 0

## all - male
connec_df2 %$% 
  wilcox.test(All, Male, paired = T)

# Wilcoxon signed rank exact test
# 
# data:  All and Male
# V = 0, p-value = 0.0009766
# alternative hypothesis: true location shift is not equal to 0

## female - male
connec_df2 %$% 
  wilcox.test(Female, Male, paired = T)

# Wilcoxon signed rank exact test
# 
# data:  Female and Male
# V = 13, p-value = 0.08301
# alternative hypothesis: true location shift is not equal to 0

connec_df %>% head()
connec_df2 %>% head()

lorenz_area %>% ungroup()
area_connec_df <- lorenz_area %>% 
  ungroup() %>% 
  mutate(sex = sex %>% as.factor() %>% fct_relevel(c("A", "F", "M"))) %>% 
  arrange(sex) %>% 
  cbind(connectance = connec_df$connectance)

## plot
p <- area_connec_df %>% 
  ggplot(aes(connectance, area, shape = sex, fill = sex))+
  theme(
    legend.position = "none",
    axis.text = element_text(size = 20, color = "black"),
    axis.title = element_text(size = 20, color = "black")
  )+
  geom_point(size = 5)+
  labs(x = "Connectance", y = "Area")+
  scale_shape_manual(values = c("A" = 22, "F" = 21, "M" = 21))+
  scale_fill_manual(values = c("A" = "gray", "F" = "white", "M" = "black"))
p  

## save
ggsave(p, filename = "area_connectance.png", height = 6, width = 6, bg = "white")


### connectance, area
area_connec_df %>% as_tibble() %>% head()

lmer_res <- area_connec_df %>% 
  as_tibble() %>% 
  lmer(data = ., area ~ connectance + sex +(1|net))


# -------------------------------------------------------------------------
# areas are gathered into a single dataframe


whole_df <- lorenz_area %>% 
  ungroup() %>% 
  rbind(lorenz_area_flower %>% ungroup()) %>% 
  mutate(remove = rep(c("visitor", "flower"), each = 33)) %>% 
  mutate(sex = sex %>% as.factor() %>% fct_relevel(c("A", "F", "M"))) %>% 
  arrange(sex, remove) %>% 
  mutate(connec = rep(connec_df$connectance, 2)) %>% 
  mutate(across(where(is.character), \(x)as.factor(x)))

whole_df %>%  
  filter(remove == "flower") %$% 
  cor.test(area, connec)

# Pearson's product-moment correlation
# 
# data:  area and connec
# t = -0.74406, df = 31, p-value = 0.4624
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.4550747  0.2208940
# sample estimates:
#        cor 
# -0.1324601 

whole_df %>%  
  filter(remove == "visitor") %$% 
  cor.test(area, connec)

# Pearson's product-moment correlation
# 
# data:  area and connec
# t = -1.3286, df = 31, p-value = 0.1937
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.5329497  0.1208284
# sample estimates:
#        cor 
# -0.2321088 

## plot
whole_df %>% head()

p <- whole_df %>% 
  mutate(remove = remove %>% fct_relevel(c("visitor", "flower"))) %>% 
  ggplot(aes(connec, area, shape = sex, fill = sex))+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = 20, color = "black"),
    axis.title = element_text(size = 20, color = "black"),
    strip.background = element_rect(fill = "transparent", color = "transparent"),
    strip.text = element_text(size = 24, color = "black", hjust = 0)
  )+
  geom_point(size = 5)+
  labs(x = "Connectance", y = "Area")+
  scale_shape_manual(values = c("A" = 22, "F" = 21, "M" = 21))+
  scale_fill_manual(values = c("A" = "gray", "F" = "white", "M" = "black"))+
  facet_wrap(remove~., ncol = 2, scale = "free",
             labeller = as_labeller(c("visitor" = "(a) Visitor removal", "flower" = "(b) Flower removal")))
p  

## save
ggsave(p, filename = "area_connectance.png", height = 6, width = 12, bg = "white")




## visitor removal

lmer_res <- whole_df %>% 
  #filter(remove == "visitor") %>% 
  lmer(data =., area ~ sex + (1|connec) + remove + (1|net))
anova(lmer_res)
# Type III Analysis of Variance Table with Satterthwaite's method
#           Sum Sq    Mean Sq NumDF  DenDF F value   Pr(>F)   
# sex    0.0056863 0.00284314     2 21.090  6.7797 0.005325 **
# connec 0.0001750 0.00017501     1 28.199  0.4173 0.523494  

# -------------------------------------------------------------------------

##
lorenz_area %>% ungroup() %>% 
  group_by(sex) %>% 
  summarise(area_mean = mean(area),
            area_median = median(area),
            area_sd = sd(area))

##
lorenz_area_flower %>% ungroup() %>% 
  group_by(sex) %>% 
  summarise(area_mean = mean(area),
            area_median = median(area),
            area_sd = sd(area))

connec_df2 %>% head()

netlevel_A <- datA %>% future_map(networklevel)
netlevel_F <- datF %>% future_map(networklevel)
netlevel_M <- datM %>% future_map(networklevel)


# -------------------------------------------------------------------------



## removing insects box plot

p_visitor <- lorenz_area %>% 
  ggplot(aes(sex, area))+
  theme_bw()+
  geom_jitter(height = 0, width = 0.2, size = 0.7)+
  geom_boxplot(fill = "transparent", coef = Inf, outliers = F)+
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    legend.position = "none",
    title = element_text(size = 10)
  )+
  labs(y = "Bowl-shaped area", title = "(a) Visitor Removal")+
  scale_x_discrete(labels = c("A" = "Original", "F" = "Female", "M" = "Male"))+
  #scale_color_manual(values = c("A" = "black", "F" = "tomato", "M" = "darkgreen"))+
  annotate("segment", x = 1, xend = 2, y = 0.12)+
  annotate("text", label = "*", x = 1.5, y = 0.11, size = 5)+
  annotate("segment", x = 1.1, xend = 2.9, y = 0.25)+
  annotate("text", label = "**", x = 2, y = 0.26, size = 5)
#ylim(c(0, 0.3))
p_visitor  


## removing flowers

p_flower <- lorenz_area_flower %>% 
  ggplot(aes(sex, area))+
  theme_bw()+
  geom_jitter(height = 0, width = 0.2, size = 0.7)+
  geom_boxplot(fill = "transparent", coef = Inf, outliers = F)+
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    legend.position = "none",
    title = element_text(size = 10)
  )+
  labs(y = "Bowl-shaped area", title = "(b) Flower Removal")+
  scale_x_discrete(labels = c("A" = "Original", "F" = "Female", "M" = "Male"))+
  #scale_color_manual(values = c("A" = "black", "F" = "tomato", "M" = "darkgreen"))+
  annotate("segment", x = 2, xend = 2.8, y = 0.016)+
  annotate("text", label = "**", x = 2.45, y = 0.005, size = 5)+
  annotate("segment", x = 1.1, xend = 2.9, y = 0.1)+
  annotate("text", label = "***", x = 2, y = 0.103, size = 5)
#ylim(c(0, 0.3))
p_flower  

p_2 <- p_visitor + p_flower + plot_layout(ncol = 1, axes = "collect")
p_2

p_2 %>% 
  ggsave(filename = "areacomparison_patchwork.jpeg", 
         height = 12, width = 12, units = "cm", dpi = 300, bg = "white")




