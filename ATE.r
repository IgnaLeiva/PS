
library(tidyverse)
library(twang)

dat_hiv <- readRDS("hiv.Rds")

dat_hiv <- dat_hiv %>% 
  mutate(sex = (as.numeric(sex) -1),
         employed = as.numeric(employed)-1,
         insurance = as.numeric(insurance)-1,
         education.level = as.numeric(education.level)-1,
         VL.dectable = factor(VL.dectable))


#class
sapply(dat_hiv, class)


table(dat_hiv$STR, dat_hiv$adherent)
table(dat_hiv$STR, dat_hiv$adherent2)

# Formula for the propensity score model
# Outcome: adherent
# Exposure: STR

names(dat_hiv)
covs <- c(
    "sex", 
    "age", 
    "adherent2",
    #"VL.dectable", 
    "employed", 
    "insurance", 
    "comorbidities", 
    "migrant", 
    "education.level", 
    "timeOnART", 
    "housing", 
    "self.stigma", 
    "addiction", 
    "mental.health",
    "comuna", 
    #"total_pill", # this was used to create SRT
    "burden_pills"
)


formula <- as.formula(paste0("STR ~", paste(covs, collapse = "+")))

# Logsitic regresion
regression <-glm(formula, family = "binomial", data = dat_hiv)
dat_hiv$PS <- regression %>% predict(dat_hiv,  type = "response") 

boxplot(PS ~ STR, data= dat_hiv, ylab="Propensity score", xlab="PPI prescription") 


# Propensity score twang package
set.seed(123)
ps.gbm = ps(formula,
              data = dat_hiv,
              n.trees=5000,
              interaction.depth=2,
              shrinkage=0.01,
              estimand = "ATE",
              stop.method=c("es.mean","ks.max"),
              n.minobsinnode = 10,
              n.keep = 1,
              n.grid = 25,
              ks.exact = NULL,
              verbose=FALSE)


dat_hiv$ps_m <- ps.gbm$ps$es.mean.ATE     


ggplot(dat_hiv, aes(x=ps_m, color=factor(STR),fill = factor(STR))) +  
  geom_density(position="identity",bins = 50, alpha=0.5)+    
  labs(x = "Propensity score", y="Number of individuals" )


boxplot(ps_m ~ STR, data= dat_hiv, ylab="Propensity score", xlab="PPI prescription") 

x <- glm(STR ~ 1, family = 'binomial', data = dat_hiv)
p <- predict(x, type = "response", newdata = dat_hiv) %>% unique() 
dat_hiv$sw <- ifelse(dat_hiv$STR == 1, p/dat_hiv$ps_m, (1-p)/(1-dat_hiv$ps_m))

dta.survey <- svydesign(ids=~1, # 1 indicating no cluster
                          weights=~sw, 
                          data=dat_hiv) 

detectable <- svyglm(VL.dectable ~STR, 
                  # family = stats::quasibinomial(), #stats::binomial()
                  family = stats::binomial(), #stats::binomial()
                  design=dta.survey)

broom::tidy(detectable, conf.int = T, exponentiate = TRUE) 

adherent <- svyglm(adherent2 ~STR, 
                  # family = stats::quasibinomial(), #stats::binomial()
                  family = stats::binomial(), #stats::binomial()
                  design=dta.survey)
broom::tidy(adherent, conf.int = T, exponentiate = TRUE) 




# Adeherence
#unadjusted

adherent_un <- glm(VL.dectable ~ STR, family = 'binomial', data = dat_hiv)
broom::tidy(adherent_un, conf.int = T, exponentiate = TRUE) 

library(MatchIt)
nnmatch <- matchit(STR~ps_m, data = dat_hiv, distance="glm",
                   link="linear.logit", method= "nearest", ratio = 1 , caliper =0.1) 
summary(nnmatch, standardize=TRUE) 
ps_match_data <- match.data(nnmatch) 
tab_match <- CreateTableOne(vars, strata="ppi", ps_match_data, test=FALSE)  
print(tab_match, smd = TRUE) 

matched_model <- glm(VL.dectable ~ STR, family = 'binomial', data = ps_match_data)
broom::tidy(matched_model, exponentiate=TRUE, conf.int=TRUE) 
