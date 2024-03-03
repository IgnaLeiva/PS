library(tidyverse)
library(patchwork)
detectable <- data.frame(
outcome = rep("Virally non-suppressed", 3),
Model = factor(c("Unadjusted", "Matched", "Weights"), levels = c("Unadjusted", "Matched", "Weights")),
RR = c(0.54, 0.78, 0.61),
conf.low = c(0.33, 0.35, 0.34),
conf.high =c(0.77, 0.98, 0.89))

adherent <- data.frame(
outcome = rep("Adherent", 3),
Model = factor(c("Unadjusted", "Matched", "Weights"), levels = c("Unadjusted", "Matched", "Weights")),
RR = c(1.81, 1.30, 1.20),
conf.low = c(1.44, 1.20, 1.10),
conf.high =c(2.30, 1.95, 1.67))

data.outcomes <- rbind(detectable, adherent)

data_plot <- data.outcomes  %>%
mutate(estimate_lab = case_when(
    Model == 'Model' ~ 'RR (95% CI)',
    TRUE ~ paste0(sprintf( '%.2f',(round(RR,2))), ' (',sprintf( '%.2f',(round(conf.low,2))), '-', sprintf( '%.2f',(round(conf.high,2))), ')')))


left <- data_plot |>
ggplot(aes(y = fct_rev(Model))) +
  geom_text(
    aes(x = 0, label = fct_rev(Model)),
    hjust = 0,
    fontface = ifelse(data_plot$Model == "Model", "bold", "plain"))+
  theme_void() +
  coord_cartesian(xlim = c(0, 1))+
annotate("text", x = 0, y = 6.5, label = "Model", hjust = 0, fontface = "bold")




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
  annotate("text", x = 0.4, y = 0, label = "Decreased risk", size = 3, color = 'grey40',fontface = 'italic') +
  annotate("text", x = 1.6, y = 0, label = "Increased risk",size = 3,  color = 'grey40',fontface = 'italic')

#Right


right <- data_plot |>
  ggplot(aes(y = fct_rev(Model))) +
  geom_label(aes(y = fct_rev(Model), x = 0, label = estimate_lab),
             size = 3.5, fill = 'white',
             label.size = NA,
             fontface = ifelse(data_plot$estimate_lab  == "1.05 (1.02-1.12)", "bold", "plain"),
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

#ggsave("plot/plot.svg", width = 8, height = 4)

