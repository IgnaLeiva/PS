#Library
library(tidyverse)
library(survival)
library(survey)
library(patchwork)
setwd("/Users/igna/Documents/R projects/PS")
# Data
dat <- readRDS("datwithps.Rds")
vars <- c("age","gender", "eth5", "imd_person", "calendarperiod", 
          "prior_gast_cancer", "recent_gerd", "recent_peptic", "bmi", "n_consult")

# Andjusted effect
andjusted_effect <- coxph(Surv(time, died) ~ ppi, data = dat)
broom::tidy(andjusted_effect, exponentiate=TRUE, conf.int=TRUE) 

# Adjustment set
formula0 <- as.formula(paste0("Surv(time, died) ~ ppi +", paste(vars, collapse = "+")))
         
adj_set <-coxph(formula0, data = dat)
broom::tidy(adj_set, exponentiate=TRUE, conf.int=TRUE) 

# Add propensity score

formula <- as.formula(paste0("Surv(time, died) ~ ppi +", paste(vars, collapse = "+"), "+ps_m"))
         
adj_setps <-coxph(formula, data = dat)
broom::tidy(adj_setps, exponentiate=TRUE, conf.int=TRUE) 

# weight

dat$ATEw <- ifelse(dat$ppi== 1, 1/dat$ps_m, 1/(1- dat$ps_m))

weight <- coxph(Surv(time, died) ~ ppi, weights = ATEw, data = dat)
broom::tidy(weight, exponentiate = TRUE, conf.int = TRUE)  

# Design packege for std error
cox_design_w <-svydesign(id=~1, weights=~ATEw, data=dat)
#rpbc<-as.svrepdesign(ATTwg_design)

weight2 <- svycoxph(Surv(time, died)~ ppi, design=cox_design_w)
broom::tidy(weight2, exponentiate=TRUE, conf.int=TRUE)


# matching
library(MatchIt)
nnmatch <- matchit(ppi~ps_m, data = dat, distance="glm",
                   link="linear.logit", method= "nearest", ratio = 1 , caliper =0.1) 
summary(nnmatch, standardize=TRUE) 
ps_match_data <- match.data(nnmatch) 

matched_model <- coxph(Surv(time, died) ~ ppi, data = ps_match_data)
broom::tidy(matched_model, exponentiate=TRUE, conf.int=TRUE)

# HR and confident intervals

models <- c("andjusted_effect", "adj_set", "adj_set_ps", "weight", "matched")
and <-broom::tidy(andjusted_effect, exponentiate=TRUE, conf.int=TRUE)[c(2,6:7)]
addj <- broom::tidy(adj_set, exponentiate=TRUE, conf.int=TRUE)[1,c(2,6:7)]
addj_ps <- broom::tidy(adj_setps, exponentiate=TRUE, conf.int=TRUE)[1,c(2,6:7)]
weights <- broom::tidy(weight, exponentiate = TRUE, conf.int = TRUE)[c(2,7:8)]
match <- broom::tidy(matched_model, exponentiate=TRUE, conf.int=TRUE)[1,c(2,6:7)]

## HR and confident intervals
data_hr<- data.frame(rbind(and, addj, add_ps, weights, match ), models)
# Modify some results

dta_hr2 <- data.frame(
Model = factor(c("Undajusted", "Adjusment set", "PS inclusion", "Weights", "Matched"), levels = c("Undajusted", "Adjusment set", "PS inclusion", "Weights", "Matched")),
HR = c(1.68, 1.30, 1.1, 1.05, 1.34),
conf.low = c(1.61, 1.26, 0.99, 1.02, 1.28),
conf.high =c(1.76, 1.40, 1.18, 1.12, 1.42))

data_plot <- dta_hr2 %>%
mutate(estimate_lab = case_when(
    Model == 'Model' ~ 'OR (95% CI)',
    TRUE ~ paste0(sprintf( '%.2f',(round(HR,2))), ' (',sprintf( '%.2f',(round(conf.low,2))), '-', sprintf( '%.2f',(round(conf.high,2))), ')')))


left <- data_plot |>
ggplot(aes(y = fct_rev(Model))) +
  geom_text(
    aes(x = 0, label = fct_rev(Model)),
    hjust = 0,
    fontface = ifelse(data_plot$Model == "Model", "bold", "plain"))+
  theme_void() +
  coord_cartesian(xlim = c(0, 1))+
annotate("text", x = 0, y = 5.5, label = "Model", hjust = 0, fontface = "bold")




# Central

central <- data_plot |>
  ggplot(aes(y = fct_rev(Model),x = HR), na.rm = TRUE) +
  geom_vline(xintercept = 1, linetype = 'dashed',linewidth = 0.3 ) +
  geom_point(shape = 19,size = 2, 
             position = position_dodge(.5), 
             show.legend = T,
             color = '#471164FF') +
  geom_errorbar(aes(xmin= conf.low, xmax = conf.high),
                width =.25,
                size = 0.8,
                color = '#471164FF',
                position = position_dodge(.5), 
                show.legend = T) +
  labs(x= '\n\n  Hazard Ratio (95% Confidence Interval)')+
  
  # Zooming the axis
  coord_cartesian(ylim = c(1,5), xlim=c(0.4,1.8), clip="off") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y =  element_blank(),
        #axis.text.x =  element_blank(), 
        axis.title.x = element_text(size = 9),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(size = 1, color="grey50") ,
        strip.text = element_text(colour = 'white', face = 'bold', size = 12.5),
        legend.position = 'right',
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 11, face = 'bold'))+
  annotate("text", x = 0.4, y = 0, label = "Increased risk", size = 3, color = 'grey40',fontface = 'italic') +
  annotate("text", x = 1.6, y = 0, label = "Decreased risk",size = 3,  color = 'grey40',fontface = 'italic')

#Right


right <- data_plot |>
  ggplot(aes(y = fct_rev(Model))) +
  geom_label(aes(y = fct_rev(Model), x = 0, label = estimate_lab),
             size = 3.5, fill = 'white',
             label.size = NA,
             #fontface = ifelse(dta$estimate_lab == "OR (95% CI)", "bold", "plain"),
             #hjust = 0,
             position = position_dodge(.5),
             show.legend = F,
             inherit.aes = F) +

  annotate("text", x = 0, y = 5.5, label = "HR (95% CI)", fontface = "bold")+

  theme_void()

# Leyaut one
layout <- c(
  area(t = 0, l = 0, b = 30, r = 2), # left plot, starts at the top of the page (0) and goes 30 units down and 3 units to the right
  area(t = 1, l = 3, b = 30, r = 3.5), # middle plot starts a little lower (t=1) because there's no title. starts 1 unit right of the left plot (l=4, whereas left plot is r=3), goes to the bottom of the page (30 units), and 6 units further over from the left plot (r=9 whereas left plot is r=3)
  area(t = 0, l = 4, b = 30, r = 5) # right most plot starts at top of page, begins where middle plot ends (l=9, and middle plot is r=9), goes to bottom of page (b=30), and extends two units wide (r=11)
)
# final plot arrangement
left + central + right + plot_layout(design = layout)

ggsave("plot/plot.svg", width = 8, height = 4)





