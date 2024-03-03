library(tidyverse)
library(ggridges)
# data from the propensity_score.R file
sampledData <- readRDS('sampledData.Rds')
sampledData<- sampledData %>% 
  mutate(Treatment = factor(ifelse(ppi == 1, 'Intervention', 'Comparator'), levels = c('Intervention', 'Comparator'))) 

##########################################
#                                        #
#     Stratification plots               #
#                                        #   
##########################################

#----------------------------------------
#       Scatterplot plot by strata
#----------------------------------------

stra_plot <- sampledData |> 
  ggplot(aes(x= Treatment,y = ps_m, color = Treatment, alpha = ps_m))+
  geom_point(size = 1,position = position_jitter(w = 0.25, h = 0)) +
  geom_hline(yintercept = c(0.6, .8, 1), color = "grey30", linetype = "dashed") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = 'none') +
  labs(y= 'Propensity Score')

ggsave(plot = stra_plot, "plot/stra_plot.png", width = 5, height = 3)

#----------------------------------------
#       Density plot by strata
#----------------------------------------

sampledData <- sampledData |> 
  mutate(strata_ps = factor(case_when(
    ps_m < 0.2 ~ "<0.2",
    ps_m >= 0.2 & ps_m < 0.4 ~ "0.2-0.4",
    ps_m >= 0.4 & ps_m < 0.6 ~ "0.4-0.6",
    ps_m >= 0.6 & ps_m < 0.8 ~ "0.6-0.8",
    ps_m >= 0.8 ~ ">0.8"
), levels = c("<0.2","0.2-0.4","0.4-0.6","0.6-0.8",">0.8")))


# ggplot(sampledData, aes(x = ps_m, y = strata_ps)) +
#   geom_density_ridges(aes(fill = paste(strata_ps, ppi))) +
#  # scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +
#   theme(legend.position = "none")

plot_strata<-sampledData |> 
  ggplot(aes(x = ps_m, y = strata_ps, fill = Treatment)) +
  geom_density_ridges(alpha = 0.4) +
  scale_fill_manual(values =c("#404080", "#69b3a2")) +
  theme_bw() +
  #guides(fill = guide_legend(title = "Treatment"), alpha = FALSE) +
  labs(x = " Propensity score value",
       y = "Propensity Score strata")

ggsave(plot = plot_strata, "plot/stra_plot2.svg", width = 8, height = 4)


##########################################
#                                        #
#          Weighting plots               #
#                                        #   
##########################################

#----------------------------------------
#       Density plot by strata
#----------------------------------------
sampledData$ATEwg <- ifelse(sampledData$ppi== 1, 1/sampledData$ps_m, 1/(1- sampledData$ps_m)) 

## No weight
noWeights <- sampledData |> 
  ggplot(aes(x = ps_m, fill = Treatment)) +
  geom_density(alpha = 0.4) +
   scale_fill_manual(values =c("#404080", "#69b3a2")) +
  theme_bw() +
  labs(title= "Pre-weighting",
      x = "Propensity score value",
       y = "Propensity Score strata")
ggsave(plot = noWeights, "plot/noWeights.svg", width = 5, height = 3)  

## Weight
Weights <-sampledData |> 
  mutate(Treatment = factor(ifelse(ppi == 1, 'Intervention', 'Comparator'), levels = c('Intervention', 'Comparator'))) |> 
  ggplot(aes(x = ps_m, fill = Treatment, weight = ATEwg)) +
  geom_density(alpha = 0.4) +
   scale_fill_manual(values =c("#404080", "#69b3a2")) +
  theme_bw() +
  labs(title= "Post-weighting",
  x = "Propensity score value",
       y = "Propensity Score strata")

ggsave(plot = Weights, "plot/Weights.svg", width = 5, height = 3)

#----------------------------------------
#       Scatterplot plot by strata
#----------------------------------------
#Treatment same value as ppi
wt_plot <- sampledData |> 
  mutate(ATEwg = ifelse(ppi== 1, 1/ps_m, 1/(1-ps_m))) |> 
  ggplot(aes(x= Treatment,y = ps_m, color = Treatment, alpha = ps_m, weight = ATEwg))+
  geom_point(size = 1,position = position_jitter(w = 0.25, h = 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = 'none') +
  labs(title= "Post-weighting",
  y= 'Propensity Score')
ggsave(plot = wt_plot, "plot/wt_plot.png", width = 5, height = 3)


##########################################
#                                        #
#          Matching  plots               #
#                                        #   
##########################################


library(MatchIt)
nnmatch <- matchit(ppi~ps_m, data = sampledData, distance="glm",
                   link="linear.logit", method= "nearest", ratio = 1 , caliper =0.2, replace = FALSE,estimand = "ATT") 
summary(nnmatch, standardize=TRUE) 
ps_match_data <- match.data(nnmatch) 


# Matrix of pairs
matched_pairs <- nnmatch$match.matrix
# To a data frame
pairs <- data.frame(a =rownames(matched_pairs), b = matched_pairs[,1]) |> 
  mutate(pair = row_number()) |> 
  pivot_longer(1:2,
               values_to = 'position',
               names_to = 'id') |> 
  select(position, pair) |> 
  drop_na() |> 
  group_by(pair) |> 
  mutate(n = n()) |> 
  filter(n == 2)

df_pairs <-sampledData[pairs$position,]
df_pairs$pair <- pairs$pair
df_pairs <- df_pairs[c('Treatment', 'ps_m', 'pair')]

match <- df_pairs |>  
  slice(1:60) |> 
  ggplot(aes(x= Treatment,y = ps_m, color = Treatment))+
  geom_point(shape = 16,size = 2) +
  scale_color_manual(values =c("#404080", "#69b3a2")) +
  #geom_point(size = 1,position = position_jitter(w = 0.25, h = 0)) +
  geom_line(aes(group = pair), color= 'grey45',linetype = "dashed") +
  theme_bw() +
  theme(legend.position = 'none')+
  labs(y = "Propensity Score strata")
ggsave(plot = match, "plot/match.svg", width = 6, height = 6)
