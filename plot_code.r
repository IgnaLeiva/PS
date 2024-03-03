library(tidyverse)
library(viridis)
library(hrbrthemes)
library(cobalt)


##################################
#.        Love plot              #
##################################


model <- glm(treat ~ re74 + race + married + I(re74^2) + re74:race, data = lalonde, family = "binomial")
eps <- predict(model, type = "response")

n.treated <- sum(lalonde$treat == 1)
n.control <- sum(lalonde$treat == 0)
weights <- ifelse(lalonde$treat == 1, 1/eps, 1/(1 - eps))


lalonde$black <- lalonde$race == "black"
lalonde$hispan <- lalonde$race == "hispan"
lalonde$white <- lalonde$race == "white"

## Draw love plot
love.plot = function(cov, treat,  ## cov is the matrix of covariates and treat is a vector of treatment assignment
                     weights = rep(1, length(treat)),
                     plot = F) 
{
  
  ## mean with normalized weights \sum w_i x_i / (\sum w_i)
  treat.means <- colSums(cov[treat == 1,] * weights[treat == 1])/sum(weights[treat == 1])
  treat.var <- colSums(t(t(cov[treat == 1,]) - treat.means)^2 *
                         weights[treat == 1])/sum(weights[treat == 1])
  
  control.means <- colSums(cov[treat == 0,] * weights[treat == 0])/sum(weights[treat == 0])
  control.var <- colSums(t(t(cov[treat == 0,]) - control.means)^2 *
                           weights[treat == 0])/sum(weights[treat == 0])
  
  ## the standardized mean differences for every covariate
  smd <- (treat.means - control.means)/sqrt((treat.var + control.var)/2)
  names(smd) <- colnames(cov)
  
  if (plot == T) {
    plot.data <- data.frame(smd = smd, covariates = names(smd))
    range <- max(abs(smd))
    ggplot(plot.data) + geom_point(aes(x = as.numeric(smd), y = covariates)) +
      geom_vline(xintercept = 0) + xlim(-range, range) +
      labs(x = 'Standardized Difference in Means')
  }
  return(smd)
}

raw.smd <- love.plot(lalonde[, c(2:3, 5:9, 10:12)], lalonde$treat)
weighted.smd <- love.plot(lalonde[, c(2:3, 5:9, 10:12)], lalonde$treat, weights = weights)


plot.data <- data.frame(smd = c(raw.smd, weighted.smd), 
                        covariates = c(names(raw.smd), names(weighted.smd)),
                        category = c(rep("Original", length(raw.smd)), rep("IPW", length(weighted.smd))))
range <- max(abs(plot.data$smd))

plot_A <- plot.data %>%
ggplot() + geom_point(aes(x = as.numeric(smd), y = covariates, color = category)) +
  geom_vline(xintercept = c(-0.1, 0,  0.1),
             linetype = c("dashed", "solid","dashed")) + 
  coord_cartesian(xlim = c(-0.5, 0.5)) +
 # xlim(-range, range) +
  labs(x = 'Standardized Difference in Means') +
  theme_ipsum() +
   theme(legend.position = c(0.9, 0.9))

ggsave(plot = plot_A, "plot/plot_A.svg", width = 4, height = 4)


##################################
#.        Density              #
##################################

currentDataset <- read_csv("https://raw.githubusercontent.com/gckc123/ExampleData/main/smoking_psyc_distress.csv")

formula <- as.formula(smoker ~ sex + age + indigeneity + high_school + partnered + remoteness + language + risky_alcohol)
data <- glm(formula , family = "binomial", currentDataset)

currentDataset$ps <- predict(data, type = "response")


plot_B <- currentDataset %>% 
    mutate(Smoker = ifelse(smoker == 1, "smoker", "non-smoker")) %>%
  ggplot(aes(x = ps, group = Smoker, fill = Smoker)) +
         geom_density(alpha = .6) +
  scale_fill_manual(values =c("#404080", "#69b3a2")) +
  # scale_fill_viridis(discrete = TRUE, alpha=0.6, option="D") +
   theme_ipsum() +
   theme(legend.position = c(0.8, 0.6),
         aspect.ratio = 18/22,
         legend.title = element_text(size = 10)) +
  labs(fill = "Smoking\nStatus")

ggsave(plot = plot_B, "plot/plot_B.svg", width = 4, height = 4)

##################################
#.        Boxplot              #
##################################

plot_c <-
currentDataset %>%
  mutate(Smoker = ifelse(smoker == 1, "smoker", "non-smoker")) %>%
ggplot(aes(x = Smoker, y = ps, fill = Smoker)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6, option="D") +
    theme_ipsum() +
  coord_flip() +
    theme(
      legend.position="none",
      plot.title = element_text(size=11),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      aspect.ratio = 18/22
    )

    ggsave(plot = plot_c, "plot/plot_C.svg", width = 4, height = 4)


pagedown::chrome_print("template.html")
